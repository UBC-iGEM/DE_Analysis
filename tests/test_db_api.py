import sqlite3

import pytest

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
    assert candidate["replicates"][0]["replicate_count"] == 3
    assert candidate["replicates"][0]["coefficient_of_variation"] == 0.1
    assert {row["assessment_type"] for row in candidate["assessments"]} == {"stability", "cross_reactivity"}
    database.close()


def test_read_only_queries_leave_database_unchanged_and_reject_writes(database_root, tmp_path):
    path = build_database(database_root, tmp_path / "project.sqlite")
    before = path.read_bytes()
    with Database(path, read_only=True) as database:
        assert database.get_candidate("genea")
        with pytest.raises(sqlite3.OperationalError, match="readonly"):
            database.connection.execute("DELETE FROM candidates")
        for value in (float("nan"), float("inf"), -1):
            with pytest.raises(ValueError, match="effect_min"):
                database.search_candidates(effect_min=value)
    assert path.read_bytes() == before


def test_targets_are_separate_from_observed_class_proxies(database_root, tmp_path):
    with Database(build_database(database_root, tmp_path / "project.sqlite")) as database:
        targets = database._rows("SELECT * FROM antibiotics WHERE proxy_role='target'")
        assert {row["target_name"] for row in targets} == {"amoxicillin", "cephalexin", "gentamicin", "tobramycin"}
        assert all(row["observed_compound"] is None for row in targets)
        proxies = database._rows("SELECT * FROM antibiotics WHERE proxy_role='proxy'")
        assert all(row["target_name"] == row["observed_compound"] for row in proxies)
