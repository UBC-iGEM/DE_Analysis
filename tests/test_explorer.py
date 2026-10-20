"""Verify local explorer routes, validation, and static assets."""

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


@pytest.mark.parametrize("field,label,invalid", [
    ("padj_max", "Maximum adjusted p-value", "abc"),
    ("padj_max", "Maximum adjusted p-value", "1.1"),
    ("padj_max", "Maximum adjusted p-value", "-0.1"),
    ("padj_max", "Maximum adjusted p-value", "NaN"),
    ("effect_min", "Minimum absolute log2 change", "-1"),
    ("effect_min", "Minimum absolute log2 change", "Infinity"),
    ("effect_min", "Minimum absolute log2 change", "1e999"),
])
def test_invalid_filter_names_field_before_querying(field, label, invalid):
    class QuerySpy:
        def search_candidates(self, **kwargs):
            pytest.fail("Invalid filters must not reach the database")

    with pytest.raises(ValueError, match=label):
        query_database(QuerySpy(), "/api/candidates", {field: [invalid]})


@pytest.mark.parametrize("value,expected", [("0", 0), ("1", 1), ("1e-3", .001), (" ", None)])
def test_numeric_filter_boundaries_and_optional_values(value, expected):
    class QuerySpy:
        def search_candidates(self, **kwargs):
            return kwargs

    result = query_database(QuerySpy(), "/api/candidates", {"padj_max": [value]})
    assert result["padj_max"] == expected
