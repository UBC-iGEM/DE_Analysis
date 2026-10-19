"""Release-facing contracts: offline reuse and an ordered, fail-fast pipeline."""

import subprocess
import sys

import pytest

from promoter_discovery.demo import build_demo
from promoter_discovery.db_api import Database
from promoter_discovery import workflow


def test_demo_runs_without_scientific_dependencies_and_keeps_evidence(tmp_path):
    path = tmp_path / 'demo.sqlite'
    subprocess.run([sys.executable, '-S', '-m', 'promoter_discovery.demo', '--database', str(path)], check=True)
    with Database(path, read_only=True) as database:
        assert len(database.search_candidates()) == 2
        assert len(database.panel()) == 2
        graph = database.candidate_network('demo-gene-1')
        assert {node['kind'] for node in graph['nodes']} == {'actor', 'promoter', 'tu', 'gene'}
        assert database.runs()[0]['reference_release'] == 'synthetic-demo'
    original = path.read_bytes()
    with pytest.raises(FileExistsError):
        build_demo(path)
    assert path.read_bytes() == original


def test_workflow_orders_stages_builds_database_once_and_stops_on_failure(monkeypatch):
    stages = []
    def execute(command, **kwargs):
        assert kwargs['check'] and kwargs['cwd'] == workflow.ROOT
        stages.append(command[2].split('.')[-1])
    monkeypatch.setattr(workflow.subprocess, 'run', execute)
    workflow.run_workflow()
    assert stages == ['setup_data', 'run_benchmark', 'build_network', 'score_candidates',
                      'promoter_selection', 'assess_candidates', 'visualize_network']
    stages.clear()
    workflow.run_workflow(reuse_fit=True)
    assert 'run_benchmark' not in stages
    def failed(command, **kwargs):
        raise subprocess.CalledProcessError(1, command)
    monkeypatch.setattr(workflow.subprocess, 'run', failed)
    with pytest.raises(subprocess.CalledProcessError):
        workflow.run_workflow()
