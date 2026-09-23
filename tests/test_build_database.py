from pathlib import Path

from promoter_discovery.build_database import build_database
from promoter_discovery.database import open_database


def test_build_database_loads_current_artifacts(tmp_path):
    root = Path(__file__).parents[1]
    database_path = tmp_path / "results" / "project.sqlite"
    output = build_database(root, database_path)

    assert output == database_path
    connection = open_database(output)
    assert connection.execute("SELECT status FROM analysis_runs").fetchone()[0] == "complete"
    assert connection.execute("SELECT COUNT(*) FROM contrasts").fetchone()[0] == 5
    assert connection.execute("SELECT COUNT(*) FROM de_results").fetchone()[0] > 1000
    assert connection.execute("SELECT COUNT(*) FROM candidates").fetchone()[0] > 0
    assert connection.execute("SELECT COUNT(*) FROM fragment_reviews").fetchone()[0] > 0
    connection.close()
