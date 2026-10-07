import pytest

from promoter_discovery.build_database import build_database
from promoter_discovery.db_api import Database
from promoter_discovery.explorer import query_database


def test_explorer_search_detail_and_graph_routes(extended_database_root, tmp_path):
    database = Database(build_database(extended_database_root, tmp_path / "project.sqlite"))
    assert len(query_database(database, "/api/candidates", {"query": ["GeneA"], "class": ["beta_lactam"]})) == 1
    assert query_database(database, "/api/candidates", {"query": ["missing"]}) == []
    assert query_database(database, "/api/candidate/genea", {})["candidate_id"] == "genea"
    graph = query_database(database, "/api/network/genea", {})
    assert {edge["id"] for edge in graph["edges"] if edge["kind"] == "regulation"} == {"I1", "I3", "I4"}
    assert query_database(database, "/api/coverage", {})["regulatory_interactions"] == 7
    assert query_database(database, "/api/promoter/P1", {})["fragments"][0]["sequence_match"] == 1
    with pytest.raises(LookupError):
        query_database(database, "/api/candidate/missing", {})
    with pytest.raises(LookupError):
        query_database(database, "/api/unknown", {})
    with pytest.raises(ValueError):
        query_database(database, "/api/candidates", {"class": ["invalid"]})
    with pytest.raises(ValueError):
        query_database(database, "/api/candidates", {"padj_max": ["2"]})
    database.close()
