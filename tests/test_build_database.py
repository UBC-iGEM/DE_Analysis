import sqlite3
import pickle

import pytest

from promoter_discovery.build_database import build_database
from promoter_discovery.db_api import Database


def test_build_database_links_all_fixture_evidence(database_root, tmp_path):
    database_path = tmp_path / "project.sqlite"
    assert build_database(database_root, database_path) == database_path
    connection = sqlite3.connect(database_path)
    counts = {
        table: connection.execute(f"SELECT COUNT(*) FROM {table}").fetchone()[0]
        for table in ("assets", "contrasts", "de_results", "candidates", "candidate_evidence",
                      "regulatory_edges", "regulatory_interactions", "transcription_units", "tu_genes", "candidate_tus",
                      "candidate_promoters", "candidate_panels")
    }
    assert counts == {
        "assets": 20, "contrasts": 3, "de_results": 6, "candidates": 1,
        "candidate_evidence": 3, "regulatory_edges": 1, "regulatory_interactions": 2,
        "transcription_units": 2, "tu_genes": 3, "candidate_tus": 1, "candidate_promoters": 1,
        "candidate_panels": 1,
    }
    assert connection.execute("PRAGMA foreign_key_check").fetchall() == []
    assert connection.execute(
        "SELECT COUNT(*) FROM candidate_evidence e JOIN de_results d "
        "ON d.run_id=e.run_id AND d.contrast_id=e.contrast_id AND d.gene_id=e.source_gene_id"
    ).fetchone()[0] == 3
    assert connection.execute(
        "SELECT regulator_id, candidate_id, tu_id, promoter_id FROM regulator_candidate_paths"
    ).fetchall() == [("crp", "genea", "T1", "P1")]
    assert connection.execute(
        "SELECT COUNT(*) FROM assets WHERE LENGTH(sha256) = 64"
    ).fetchone()[0] == counts["assets"]
    connection.close()


def test_full_reference_paths_keep_promoter_specific_effects(database_root, tmp_path):
    database_path = build_database(database_root, tmp_path / "project.sqlite")
    database = Database(database_path)
    paths = database.reference_promoter_paths("CRP")
    assert {(row["interaction_id"], row["effect"], row["promoter_id"], row["gene_id"])
            for row in paths} == {
        ("I1", "activator", "P1", "genea"),
        ("I1", "activator", "P1", "geneb"),
        ("I2", "repressor", "P2", "genec"),
    }
    assert database.reference_promoter_paths("missing") == []
    interactions = database.reference_interactions("R1", target_kind="promoter")
    assert [(row["interaction_id"], row["effect"]) for row in interactions] == [
        ("I1", "activator"), ("I2", "repressor")
    ]
    with pytest.raises(ValueError, match="target_kind"):
        database.reference_interactions("CRP", target_kind="unknown")
    database.close()


def test_reference_hash_mismatch_preserves_previous_database(database_root, tmp_path):
    database_path = build_database(database_root, tmp_path / "project.sqlite")
    stable = database_path.read_bytes()
    with (database_root / "data/references/regulondb/RISet.tsv").open("a") as handle:
        handle.write("\n")
    with pytest.raises(ValueError, match="not verified"):
        build_database(database_root, database_path)
    assert database_path.read_bytes() == stable


def test_rebuild_is_idempotent_and_failure_preserves_previous_database(database_root, tmp_path):
    database_path = tmp_path / "project.sqlite"
    build_database(database_root, database_path)
    with sqlite3.connect(database_path) as connection:
        run_id = connection.execute("SELECT run_id FROM analysis_runs").fetchone()[0]
    build_database(database_root, database_path)
    with sqlite3.connect(database_path) as connection:
        assert connection.execute("SELECT run_id FROM analysis_runs").fetchone()[0] == run_id
        assert connection.execute("SELECT COUNT(*) FROM candidate_evidence").fetchone()[0] == 3
    stable = database_path.read_bytes()
    (database_root / "results/regulatory_network/regulatory_network.pkl").write_bytes(b"invalid graph")
    with pytest.raises(pickle.UnpicklingError):
        build_database(database_root, database_path)
    assert database_path.read_bytes() == stable


def test_missing_graph_is_an_error_not_an_empty_network(database_root, tmp_path):
    (database_root / "results/regulatory_network/regulatory_network.pkl").unlink()
    with pytest.raises(FileNotFoundError, match="regulatory_network.pkl"):
        build_database(database_root, tmp_path / "project.sqlite")
