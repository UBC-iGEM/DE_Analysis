from promoter_discovery.build_database import build_database
from promoter_discovery.db_api import Database
from promoter_discovery.curation import archive_path

import json
import sqlite3
from types import SimpleNamespace

import pytest


def test_candidate_review_is_versioned_and_visible_in_summary(database_root, tmp_path):
    database_path = tmp_path / "project.sqlite"
    build_database(database_root, database_path)
    database = Database(database_path)
    run_id = database.connection.execute("SELECT run_id FROM analysis_runs").fetchone()[0]
    edit_id = database.record_candidate_review(run_id, "genea", "approved", "reviewer", "ready for construct design", priority=1, construct_ready=True)
    assert edit_id.startswith("edit-")
    row = database.get_candidate("genea", run_id)
    assert row["review_status"] == "approved"
    assert row["construct_ready"] == 1
    assert database.connection.execute("SELECT COUNT(*) FROM curation_edits").fetchone()[0] == 1
    database.close()


def test_review_order_uses_timestamps_and_audits_all_fields(database_root, tmp_path, monkeypatch):
    path = build_database(database_root, tmp_path / "project.sqlite")
    database = Database(path)
    run_id = database.runs()[0]["run_id"]
    ids = iter(("zzzz", "aaaa"))
    monkeypatch.setattr("promoter_discovery.db_api.uuid4", lambda: SimpleNamespace(hex=next(ids)))
    database.record_candidate_review(run_id, "genea", "reviewed", "A", "first", priority=5, notes="first note")
    database.record_candidate_review(run_id, "genea", "approved", "B", "second", priority=1, notes="second note")
    assert database.get_candidate("genea")["review_status"] == "approved"
    assert database.get_candidate("genea")["priority"] == 1
    history = database.review_history("candidate", "genea")
    assert [row["review"]["review_status"] for row in history] == ["reviewed", "approved"]
    assert json.loads(history[-1]["old_value"])["notes"] == "first note"
    assert json.loads(history[-1]["new_value"])["priority"] == 1
    database.close()


def test_same_run_rebuild_and_database_deletion_preserve_both_review_types(database_root, tmp_path):
    path = build_database(database_root, tmp_path / "project.sqlite")
    database = Database(path)
    run_id = database.runs()[0]["run_id"]
    database.record_candidate_review(run_id, "genea", "approved", "A", "ready")
    database.record_promoter_review(run_id, "P1", "ready_for_order", "A", "boundaries checked", True)
    database.close()
    assert archive_path(path).exists()
    build_database(database_root, path)
    database = Database(path)
    assert database.get_candidate("genea")["review_status"] == "approved"
    assert database.get_promoter("P1")["review"]["construct_ready"] == 1
    assert len(database.review_history("promoter", "P1")) == 1
    database.close()
    path.unlink()
    build_database(database_root, path)
    database = Database(path)
    assert database.get_candidate("genea")["review_status"] == "approved"
    database.close()


def test_changed_run_keeps_history_without_reapplying_approval(database_root, tmp_path):
    path = build_database(database_root, tmp_path / "project.sqlite")
    database = Database(path)
    old_run = database.runs()[0]["run_id"]
    database.record_candidate_review(old_run, "genea", "approved", "A", "old evidence")
    database.record_promoter_review(old_run, "P1", "ready_for_order", "A", "old annotation")
    database.close()
    source = database_root / "results/differential_expression/model/benchmark_provenance.json"
    source.write_text('{"new_input": true}')
    build_database(database_root, path)
    database = Database(path)
    assert database.runs()[0]["run_id"] != old_run
    assert database.get_candidate("genea")["review_status"] is None
    assert database.get_promoter("P1")["review"] is None
    assert database.review_history("candidate", "genea")[0]["run_id"] == old_run
    assert len(database.review_history("promoter", "P1")) == 1
    database.close()


def test_schema_nine_reviews_are_archived_before_upgrade(database_root, tmp_path):
    path = build_database(database_root, tmp_path / "project.sqlite")
    with sqlite3.connect(path) as connection:
        run_id = connection.execute("SELECT run_id FROM analysis_runs").fetchone()[0]
        connection.execute("UPDATE database_metadata SET value='9' WHERE key='schema_version'")
        connection.execute("INSERT INTO curation_edits(edit_id, run_id, entity_type, entity_id, field_name, new_value, editor, edited_at, reason) VALUES ('legacy', ?, 'candidate', 'genea', 'review_status', 'reviewed', 'A', '2026-01-01', 'legacy review')", (run_id,))
        connection.execute("INSERT INTO candidate_reviews(edit_id, run_id, candidate_id, review_status) VALUES ('legacy', ?, 'genea', 'reviewed')", (run_id,))
    build_database(database_root, path)
    database = Database(path)
    assert database.get_candidate("genea")["review_status"] == "reviewed"
    assert database.review_history("candidate", "genea")[0]["edit_id"] == "legacy"
    database.close()


@pytest.mark.parametrize("changes", [{"review_status": "invalid"}, {"candidate_id": "unknown"},
                                   {"editor": " "}, {"reason": ""}, {"priority": -1},
                                   {"construct_ready": "yes"}, {"run_id": "unknown"}])
def test_invalid_reviews_do_not_write_history(database_root, tmp_path, changes):
    path = build_database(database_root, tmp_path / "project.sqlite")
    database = Database(path)
    args = dict(run_id=database.runs()[0]["run_id"], candidate_id="genea", review_status="approved",
                editor="A", reason="valid")
    args.update(changes)
    with pytest.raises(ValueError):
        database.record_candidate_review(**args)
    assert not database.review_history("candidate", "genea")
    assert database.connection.execute("SELECT COUNT(*) FROM curation_edits").fetchone()[0] == 0
    database.close()
