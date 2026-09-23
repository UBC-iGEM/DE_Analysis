from promoter_discovery.build_database import build_database
from promoter_discovery.db_api import Database


def test_database_api_searches_candidates_and_panel(tmp_path):
    database_path = tmp_path / "project.sqlite"
    build_database(".", database_path)
    database = Database(database_path)
    candidates = database.search_candidates(antibiotic_class="beta_lactam")
    assert candidates
    assert all(row["antibiotic_class"] == "beta_lactam" for row in candidates)
    panel = database.panel()
    assert len(panel) == 12
    assert database.get_candidate(candidates[0]["candidate_id"]) is not None
    database.close()
