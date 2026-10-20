"""Research contracts: retain identities, context and the evidence behind reviews."""

import csv
import hashlib
import json
import sqlite3

import pytest

from promoter_discovery import build_database as importer
from promoter_discovery.database import open_database
from promoter_discovery.db_api import Database
from promoter_discovery.provenance import preserve_snapshot


def test_source_identifiers_resolve_without_collapsing_loci(database_root, tmp_path):
    with Database(importer.build_database(database_root, tmp_path / "db.sqlite"), read_only=True) as db:
        for identifier in ("genea", "B0001", "b0002"):
            assert db.get_candidate(identifier)["candidate_id"] == "genea"
            assert [row["candidate_id"] for row in db.search_candidates(query=identifier)] == ["genea"]
        assert {row["source_gene_id"] for row in db.get_candidate("genea")["evidence"]} == {"b0001", "b0002"}
        assert db.candidate_network("b0001")["candidate_id"] == "genea"
        assert db.candidate_paths("b0001") == db.candidate_paths("genea")
        assert db.get_candidate("unknown") is None


def test_ambiguous_identifier_returns_all_search_hits_but_never_picks_one(database_root, tmp_path):
    path = importer.build_database(database_root, tmp_path / "db.sqlite")
    with sqlite3.connect(path) as connection:
        run = connection.execute("SELECT run_id FROM analysis_runs").fetchone()[0]
        connection.execute("INSERT INTO candidates(run_id,candidate_id,gene_id,antibiotic_class) VALUES (?,'geneb','geneb','beta_lactam')", (run,))
        connection.execute(
            "INSERT INTO gene_identifiers SELECT run_id,'geneb',namespace,identifier,source_asset_id,source_row,source_release,mapping_basis "
            "FROM gene_identifiers WHERE identifier='b0001'"
        )
    with Database(path) as db:
        assert {row["candidate_id"] for row in db.search_candidates(query="b0001")} == {"genea", "geneb"}
        with pytest.raises(ValueError, match="Ambiguous identifier"):
            db.get_candidate("b0001")
        assert db.get_candidate("genea")["candidate_id"] == "genea"


def test_cross_run_de_and_replicate_contrasts_are_rejected(tmp_path):
    c = open_database(tmp_path / "db.sqlite")
    try:
        c.execute("INSERT INTO analysis_runs(run_id,created_at,status) VALUES ('a','now','complete'),('b','now','complete')")
        c.execute("INSERT INTO genes(gene_id) VALUES ('g')")
        c.execute("INSERT INTO contrasts(contrast_id,run_id,name,numerator,denominator,comparison_type) VALUES ('contrast','a','drug/control','drug','water','drug_control')")
        with pytest.raises(sqlite3.IntegrityError):
            c.execute("INSERT INTO de_results(run_id,contrast_id,gene_id) VALUES ('b','contrast','g')")
        with pytest.raises(sqlite3.IntegrityError):
            c.execute("INSERT INTO replicate_metrics(run_id,contrast_id) VALUES ('b','contrast')")
        for value in (-0.01, 1.01):
            with pytest.raises(sqlite3.IntegrityError):
                c.execute("INSERT INTO de_results(run_id,contrast_id,gene_id,padj) VALUES ('a','contrast','g',?)", (value,))
    finally:
        c.close()


def test_sample_context_preserves_unknown_values_and_records_its_source(database_root, tmp_path):
    with Database(importer.build_database(database_root, tmp_path / "db.sqlite"), read_only=True) as db:
        row = db.sample_qc()[0]
        assert row["replicate"] == "1"
        assert row["antibiotic_role"] == "discovery"
        assert row["medium"] == "test medium"
        assert row["source_reference"] == "synthetic fixture"
        assert row["concentration"] is None and row["exposure_minutes"] is None
        assert db.reference_coverage()["candidate_genes_with_promoters"] == 1
        assert db.runs()[0]["code_sha256"]
        assert json.loads(db.runs()[0]["code_manifest_json"])["promoter_discovery/database.py"]


@pytest.mark.parametrize("field,value", [("study_id", "wrong-study"), ("source_reference", "")])
def test_invalid_context_preserves_previous_database(database_root, tmp_path, field, value):
    path = importer.build_database(database_root, tmp_path / "db.sqlite")
    original = path.read_bytes()
    context_path = database_root / "config/study_context.json"
    context = json.loads(context_path.read_text())
    context[field] = value
    context_path.write_text(json.dumps(context))
    with pytest.raises(ValueError):
        importer.build_database(database_root, path)
    assert path.read_bytes() == original


def test_qc_metadata_mismatch_is_not_silently_merged(database_root, tmp_path):
    metadata = database_root / "results/differential_expression/model/joint_metadata.csv"
    metadata.write_text("sample_id,condition,source_file\nS1,water,S1.txt\n")
    with pytest.raises(ValueError, match="conflicts with QC"):
        importer.build_database(database_root, tmp_path / "db.sqlite")


@pytest.mark.parametrize("dose,unit", [("unknown", "mg/L"), (float("nan"), "mg/L"), (2, None)])
def test_invalid_doses_are_rejected_instead_of_becoming_evidence(database_root, tmp_path, dose, unit):
    path = database_root / "config/study_context.json"
    context = json.loads(path.read_text())
    context["conditions"]["benchmark_ceftazidime"].update(concentration=dose, concentration_unit=unit)
    path.write_text(json.dumps(context))
    with pytest.raises(ValueError):
        importer.build_database(database_root, tmp_path / "db.sqlite")


def test_legacy_sample_index_and_configured_replicate_labels(database_root, tmp_path):
    metadata = database_root / "results/differential_expression/model/joint_metadata.csv"
    metadata.write_text(",condition,source_file\nS1,benchmark_ceftazidime,S1.txt\n")
    config_path = database_root / "config/benchmark.json"
    config = json.loads(config_path.read_text())
    config["datasets"][0]["sample_name_regex"] = r"S(?P<replicate>\d+)"
    config_path.write_text(json.dumps(config))
    with Database(importer.build_database(database_root, tmp_path / "db.sqlite"), read_only=True) as db:
        assert db.sample_qc()[0]["replicate"] == "1"


def test_code_content_change_invalidates_approval_and_retains_evidence(database_root, tmp_path, monkeypatch):
    path = importer.build_database(database_root, tmp_path / "db.sqlite")
    with Database(path) as db:
        old_run = db.runs()[0]["run_id"]
        db.record_candidate_review(old_run, "genea", "approved", "curator", "reviewed old evidence")
    provenance = importer.source_provenance()
    monkeypatch.setattr(importer, "source_provenance", lambda: {**provenance, "code_sha256": "f" * 64})
    importer.build_database(database_root, path)
    with Database(path, read_only=True) as db:
        assert db.runs()[0]["run_id"] != old_run
        assert db.get_candidate("genea")["review_status"] is None
        assert db.review_history("candidate", "genea")[0]["run_id"] == old_run
    snapshots = list((path.parent / "database_snapshots" / path.stem).glob("*.sqlite"))
    assert len(snapshots) == 1
    with Database(snapshots[0], read_only=True) as old:
        assert old.runs()[0]["run_id"] == old_run
        assert old.get_candidate("genea")["review_status"] == "approved"
    importer.build_database(database_root, path)
    assert len(list(snapshots[0].parent.glob("*.sqlite"))) == 1


def test_snapshot_uses_sqlite_backup_including_committed_wal(database_root, tmp_path):
    path = importer.build_database(database_root, tmp_path / "db.sqlite")
    with Database(path) as db:
        run = db.runs()[0]["run_id"]
        db.record_candidate_review(run, "genea", "reviewed", "A", "in WAL")
        saved = preserve_snapshot(path, "different", 12)
        with sqlite3.connect(saved) as backup:
            assert backup.execute("SELECT COUNT(*) FROM candidate_reviews").fetchone()[0] == 1


def test_snapshot_failure_aborts_replacement(database_root, tmp_path, monkeypatch):
    path = importer.build_database(database_root, tmp_path / "db.sqlite")
    original = path.read_bytes()
    def fail(*args):
        raise OSError("archive unavailable")
    monkeypatch.setattr(importer, "preserve_snapshot", fail)
    with pytest.raises(OSError, match="archive unavailable"):
        importer.build_database(database_root, path)
    assert path.read_bytes() == original


def test_summary_lists_do_not_create_fictitious_tus_or_genes(database_root, tmp_path):
    reference = database_root / "data/references/regulondb/TUSet.tsv"
    with reference.open(newline="") as handle:
        tus = list(csv.DictReader(handle, delimiter="\t"))
    tus.append({**tus[0], "1)id": "T3", "2)name": "alternative TU"})
    with reference.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=tus[0], delimiter="\t")
        writer.writeheader()
        writer.writerows(tus)
    lock_path = database_root / "config/regulondb.lock.json"
    lock = json.loads(lock_path.read_text())
    lock["assets"]["TUSet.tsv"]["sha256"] = hashlib.sha256(reference.read_bytes()).hexdigest()
    lock_path.write_text(json.dumps(lock))
    summary_path = database_root / "results/promoter_review/promoter_review.csv"
    with summary_path.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    rows[0]["tu_ids"] = "T1, T3"
    rows[0]["candidate_genes"] = "genea, genea"
    with summary_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0])
        writer.writeheader()
        writer.writerows(rows)
    with Database(importer.build_database(database_root, tmp_path / "db.sqlite"), read_only=True) as db:
        assert db.connection.execute("SELECT COUNT(*) FROM transcription_units").fetchone()[0] == 3
        assert not db.connection.execute("SELECT 1 FROM transcription_units WHERE tu_id LIKE '%,%'").fetchone()
        assert not db.connection.execute("SELECT 1 FROM genes WHERE gene_id LIKE '%,%'").fetchone()


def test_unknown_reference_tu_is_rejected(database_root, tmp_path):
    path = database_root / "results/promoter_review/candidate_promoter_mapping.csv"
    path.write_text(path.read_text().replace("T1", "unknown-tu"))
    with pytest.raises(ValueError, match="reference TU membership"):
        importer.build_database(database_root, tmp_path / "db.sqlite")
