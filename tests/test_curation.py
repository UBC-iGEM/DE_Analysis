from promoter_discovery.build_database import build_database
from promoter_discovery.db_api import Database


def test_candidate_review_is_versioned_and_visible_in_summary(tmp_path):
    database_path = tmp_path / "project.sqlite"
    build_database(".", database_path)
    database = Database(database_path)
    run_id = database.connection.execute("SELECT run_id FROM analysis_runs").fetchone()[0]
    edit_id = database.record_candidate_review(run_id, "gfcc", "approved", "reviewer", "ready for construct design", priority=1, construct_ready=True)
    assert edit_id.startswith("edit-")
    row = database.get_candidate("gfcc", run_id)
    assert row["review_status"] == "approved"
    assert row["construct_ready"] == 1
    assert database.connection.execute("SELECT COUNT(*) FROM curation_edits").fetchone()[0] == 1
    database.close()
