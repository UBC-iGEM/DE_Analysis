"""Verify command-line searches, JSON output, and review commands."""

from promoter_discovery.build_database import build_database
from promoter_discovery.db_cli import main

import json
import pytest
from promoter_discovery.db_api import Database


def test_cli_json_candidate_search(database_root, tmp_path, capsys):
    database_path = tmp_path / "project.sqlite"
    build_database(database_root, database_path)
    main(["--database", str(database_path), "--json", "candidates", "--class", "beta_lactam"])
    output = capsys.readouterr().out
    assert '"antibiotic_class": "beta_lactam"' in output


def test_cli_reference_paths(database_root, tmp_path, capsys):
    database_path = tmp_path / "project.sqlite"
    build_database(database_root, database_path)
    main(["--database", str(database_path), "--json", "promoter-paths", "CRP"])
    output = capsys.readouterr().out
    assert '"interaction_id": "I1"' in output
    assert '"promoter_id": "P2"' in output


def test_cli_review_history_and_graph(extended_database_root, tmp_path, capsys):
    path = build_database(extended_database_root, tmp_path / "project.sqlite")
    database = Database(path)
    run_id = database.runs()[0]["run_id"]
    database.close()
    main(["--database", str(path), "--json", "review-candidate", "genea", "approved",
          "--run-id", run_id, "--editor", "A", "--reason", "checked", "--priority", "1"])
    assert json.loads(capsys.readouterr().out)[0]["edit_id"].startswith("edit-")
    main(["--database", str(path), "--json", "history", "candidate", "genea"])
    assert json.loads(capsys.readouterr().out)[0]["review"]["priority"] == 1
    main(["--database", str(path), "--json", "network", "genea"])
    assert json.loads(capsys.readouterr().out)[0]["nodes"]


def test_cli_bad_filters_and_missing_database_are_clear_errors(database_root, tmp_path, capsys):
    path = build_database(database_root, tmp_path / "project.sqlite")
    with pytest.raises(SystemExit) as error:
        main(["--database", str(path), "candidates", "--padj-max", "2"])
    assert error.value.code == 2
    assert "padj_max must be between 0 and 1" in capsys.readouterr().err
    missing = tmp_path / "missing.sqlite"
    with pytest.raises(SystemExit):
        main(["--database", str(missing), "candidates"])
    assert "Database not found" in capsys.readouterr().err
    assert not missing.exists()
