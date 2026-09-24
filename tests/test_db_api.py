from promoter_discovery.build_database import build_database
from promoter_discovery.db_api import Database


def test_database_api_searches_candidates_and_panel(database_root, tmp_path):
    database_path = tmp_path / "project.sqlite"
    build_database(database_root, database_path)
    database = Database(database_path)
    candidates = database.search_candidates(antibiotic_class="beta_lactam")
    assert len(candidates) == 1
    assert all(row["antibiotic_class"] == "beta_lactam" for row in candidates)
    assert candidates[0]["promoter_count"] == 1
    assert len(database.search_candidates(padj_max=0.05, effect_min=2)) == 1
    assert not database.search_candidates(padj_max=0.005, effect_min=4)
    panel = database.panel()
    assert len(panel) == 1
    candidate = database.get_candidate("genea")
    assert len(candidate["evidence"]) == 3
    assert [row["promoter_id"] for row in candidate["promoters"]] == ["P1"]
    assert [row["tu_id"] for row in candidate["transcription_units"]] == ["T1"]
    assert database.search_regulator("crp")[0]["candidate_id"] == "genea"
    database.close()
