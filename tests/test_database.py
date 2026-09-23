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
