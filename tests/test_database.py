import sqlite3

import pytest

from promoter_discovery.database import (
    SCHEMA_VERSION,
    close_database,
    open_database,
    schema_version,
    transaction,
)


def test_open_database_initializes_schema_and_enables_foreign_keys(tmp_path):
    connection = open_database(tmp_path / "project.sqlite")

    assert schema_version(connection) == SCHEMA_VERSION
    assert connection.execute("PRAGMA foreign_keys").fetchone()[0] == 1
    tables = {
        row["name"]
        for row in connection.execute(
            "SELECT name FROM sqlite_master WHERE type = 'table'"
        )
    }
    assert {"database_metadata", "analysis_runs", "assets", "promoters"} <= tables
    close_database(connection)


def test_foreign_keys_reject_orphan_records(tmp_path):
    connection = open_database(tmp_path / "project.sqlite")

    with pytest.raises(sqlite3.IntegrityError):
        connection.execute(
            "INSERT INTO samples(sample_id, study_id, condition) VALUES (?, ?, ?)",
            ("sample-1", "missing-study", "control"),
        )

    close_database(connection)


def test_transaction_rolls_back_on_error(tmp_path):
    connection = open_database(tmp_path / "project.sqlite")

    with pytest.raises(RuntimeError):
        with transaction(connection):
            connection.execute(
                "INSERT INTO studies(study_id, accession) VALUES (?, ?)",
                ("study-1", "GSE-test"),
            )
            raise RuntimeError("stop")

    assert connection.execute("SELECT COUNT(*) FROM studies").fetchone()[0] == 0
    close_database(connection)


def test_de_results_require_run_contrast_and_gene(tmp_path):
    connection = open_database(tmp_path / "project.sqlite")
    connection.execute(
        "INSERT INTO analysis_runs(run_id, created_at, status) VALUES (?, ?, ?)",
        ("run-1", "2026-01-01T00:00:00Z", "complete"),
    )
    connection.execute(
        "INSERT INTO contrasts(contrast_id, run_id, name, numerator, denominator, comparison_type) "
        "VALUES (?, ?, ?, ?, ?, ?)",
        ("contrast-1", "run-1", "drug_vs_water", "drug", "water", "drug_control"),
    )
    connection.execute("INSERT INTO genes(gene_id) VALUES (?)", ("b0001",))
    connection.execute(
        "INSERT INTO de_results(run_id, contrast_id, gene_id, padj, regulation) VALUES (?, ?, ?, ?, ?)",
        ("run-1", "contrast-1", "b0001", 0.01, "upregulated"),
    )
    with pytest.raises(sqlite3.IntegrityError):
        connection.execute(
            "INSERT INTO de_results(run_id, contrast_id, gene_id) VALUES (?, ?, ?)",
            ("run-1", "contrast-1", "missing-gene"),
        )
    assert connection.execute("SELECT COUNT(*) FROM de_results").fetchone()[0] == 1
    close_database(connection)


def test_regulatory_relationships_can_traverse_regulator_to_promoter(tmp_path):
    connection = open_database(tmp_path / "project.sqlite")
    connection.execute("INSERT INTO regulators(regulator_id, name) VALUES (?, ?)", ("crp", "CRP"))
    connection.execute("INSERT INTO genes(gene_id, canonical_name) VALUES (?, ?)", ("b0001", "geneA"))
    connection.execute("INSERT INTO transcription_units(tu_id, name) VALUES (?, ?)", ("tu-1", "tuA"))
    connection.execute("INSERT INTO promoters(promoter_id, tu_id, name) VALUES (?, ?, ?)", ("p-1", "tu-1", "pA"))
    connection.execute("INSERT INTO tu_genes(tu_id, gene_id) VALUES (?, ?)", ("tu-1", "b0001"))
    connection.execute("INSERT INTO tu_promoters(tu_id, promoter_id) VALUES (?, ?)", ("tu-1", "p-1"))
    connection.execute("INSERT INTO regulatory_edges(edge_id, regulator_id, target_gene_id, edge_type) VALUES (?, ?, ?, ?)", ("e-1", "crp", "b0001", "transcriptional"))
    row = connection.execute(
        "SELECT r.name, p.promoter_id FROM regulatory_edges e "
        "JOIN regulators r ON r.regulator_id = e.regulator_id "
        "JOIN tu_genes tg ON tg.gene_id = e.target_gene_id "
        "JOIN tu_promoters tp ON tp.tu_id = tg.tu_id "
        "JOIN promoters p ON p.promoter_id = tp.promoter_id"
    ).fetchone()
    assert tuple(row) == ("CRP", "p-1")
    close_database(connection)


def test_panel_rejects_duplicate_promoter_and_transcription_group(tmp_path):
    connection = open_database(tmp_path / "project.sqlite")
    connection.execute("INSERT INTO analysis_runs(run_id, created_at, status) VALUES ('run-1', 'now', 'complete')")
    connection.execute("INSERT INTO genes(gene_id) VALUES ('b0001'), ('b0002')")
    connection.execute("INSERT INTO candidates(candidate_id, run_id, gene_id, antibiotic_class) VALUES ('c1', 'run-1', 'b0001', 'beta_lactam'), ('c2', 'run-1', 'b0002', 'beta_lactam')")
    connection.execute("INSERT INTO candidate_panels(run_id, panel_id, candidate_id, promoter_id, antibiotic_class, transcription_group) VALUES ('run-1', 'panel-1', 'c1', 'p1', 'beta_lactam', 'tu1')")
    with pytest.raises(sqlite3.IntegrityError):
        connection.execute("INSERT INTO candidate_panels(run_id, panel_id, candidate_id, promoter_id, antibiotic_class, transcription_group) VALUES ('run-1', 'panel-1', 'c2', 'p1', 'beta_lactam', 'tu2')")
    with pytest.raises(sqlite3.IntegrityError):
        connection.execute("INSERT INTO candidate_panels(run_id, panel_id, candidate_id, promoter_id, antibiotic_class, transcription_group) VALUES ('run-1', 'panel-1', 'c2', 'p2', 'beta_lactam', 'tu1')")
    close_database(connection)


def test_qc_and_fragment_reviews_require_valid_source_entities(tmp_path):
    connection = open_database(tmp_path / "project.sqlite")
    connection.execute("INSERT INTO studies(study_id, accession) VALUES ('study-1', 'GSE-test')")
    connection.execute("INSERT INTO samples(sample_id, study_id, condition) VALUES ('sample-1', 'study-1', 'control')")
    connection.execute("INSERT INTO analysis_runs(run_id, created_at, status) VALUES ('run-1', 'now', 'complete')")
    connection.execute("INSERT INTO sample_qc(run_id, study_id, sample_id, library_size) VALUES ('run-1', 'study-1', 'sample-1', 1000)")
    connection.execute("INSERT INTO promoters(promoter_id) VALUES ('p-1')")
    connection.execute("INSERT INTO fragment_reviews(fragment_id, run_id, promoter_id, sequence_match) VALUES ('f-1', 'run-1', 'p-1', 1)")
    assert connection.execute("SELECT COUNT(*) FROM sample_qc").fetchone()[0] == 1
    assert connection.execute("SELECT COUNT(*) FROM fragment_reviews WHERE sequence_match = 1").fetchone()[0] == 1
    close_database(connection)


def test_candidate_evidence_preserves_two_loci_and_requires_matching_de_rows(tmp_path):
    connection = open_database(tmp_path / "project.sqlite")
    connection.execute("INSERT INTO analysis_runs(run_id, created_at, status) VALUES ('run-1', 'now', 'complete')")
    connection.execute("INSERT INTO genes(gene_id) VALUES ('insi2'), ('b4284'), ('b4708')")
    connection.execute("INSERT INTO contrasts(contrast_id, run_id, name, numerator, denominator, comparison_type) VALUES ('drug', 'run-1', 'drug_vs_water', 'drug', 'water', 'drug_control')")
    connection.execute("INSERT INTO candidates(run_id, candidate_id, gene_id, antibiotic_class) VALUES ('run-1', 'insi2', 'insi2', 'aminoglycoside')")
    for locus in ('b4284', 'b4708'):
        connection.execute("INSERT INTO de_results(run_id, contrast_id, gene_id) VALUES ('run-1', 'drug', ?)", (locus,))
        connection.execute("INSERT INTO candidate_evidence(run_id, candidate_id, contrast_id, source_gene_id) VALUES ('run-1', 'insi2', 'drug', ?)", (locus,))
    assert connection.execute("SELECT COUNT(*) FROM candidate_evidence").fetchone()[0] == 2
    with pytest.raises(sqlite3.IntegrityError):
        connection.execute("INSERT INTO candidate_evidence(run_id, candidate_id, contrast_id, source_gene_id) VALUES ('run-1', 'insi2', 'drug', 'insi2')")
    close_database(connection)


def test_old_generated_schema_is_not_silently_relabelled(tmp_path):
    path = tmp_path / "old.sqlite"
    connection = sqlite3.connect(path)
    connection.execute("CREATE TABLE database_metadata(key TEXT PRIMARY KEY, value TEXT NOT NULL)")
    connection.execute("INSERT INTO database_metadata VALUES ('schema_version', '7')")
    connection.commit()
    connection.close()
    with pytest.raises(ValueError, match="rebuild the generated database"):
        open_database(path)
