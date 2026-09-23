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
