from promoter_discovery.build_database import build_database
from promoter_discovery.db_cli import main


def test_cli_json_candidate_search(database_root, tmp_path, capsys):
    database_path = tmp_path / "project.sqlite"
    build_database(database_root, database_path)
    main(["--database", str(database_path), "--json", "candidates", "--class", "beta_lactam"])
    output = capsys.readouterr().out
    assert '"antibiotic_class": "beta_lactam"' in output
