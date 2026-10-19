"""Create a tiny, fictional example database using only the Python standard library."""

from __future__ import annotations

import argparse
import json
import os
import tempfile
from pathlib import Path

from .database import close_database, open_database, transaction


DEFAULT_DEMO = Path("results/demo.sqlite")


def build_demo(path: str | Path = DEFAULT_DEMO) -> Path:
    """Create a separate demo; never overwrite scientific results or reviews."""
    path = Path(path)
    if path.exists():
        raise FileExistsError(f"Demo already exists: {path}. Open it or choose --database <new-path>.")
    path.parent.mkdir(parents=True, exist_ok=True)
    handle, name = tempfile.mkstemp(prefix="prism-demo-", suffix=".sqlite", dir=path.parent)
    os.close(handle)
    temporary = Path(name)
    connection = open_database(temporary)
    try:
        with transaction(connection):
            connection.execute(
                "INSERT INTO analysis_runs(run_id, created_at, status, reference_release, summary_json) VALUES (?, ?, 'complete', ?, ?)",
                ("demo-v1", "2026-01-01T00:00:00+00:00", "synthetic-demo",
                 json.dumps({"description": "Fictional software example; no biological claims", "reference_coverage": {
                     "regulatory_interactions": 2, "promoters": 2, "transcription_units": 2,
                     "interactions_without_promoter": 0}})),
            )
            connection.execute("INSERT INTO regulatory_actors VALUES ('DEMO-R1', 'Example regulator', 'TF', 'synthetic-demo')")
            for index, antibiotic_class in enumerate(("beta_lactam", "aminoglycoside"), 1):
                gene, promoter, tu, contrast = (f"demo-{kind}-{index}" for kind in ("gene", "promoter", "tu", "contrast"))
                connection.execute("INSERT INTO genes(gene_id, canonical_name) VALUES (?, ?)", (gene, f"Example gene {index}"))
                connection.execute("INSERT INTO transcription_units(tu_id, name, source_release) VALUES (?, ?, 'synthetic-demo')", (tu, f"Example TU {index}"))
                connection.execute("INSERT INTO promoters(promoter_id, name, sequence, annotation_status, source_release) VALUES (?, ?, 'ACGTACGT', 'fictional_example', 'synthetic-demo')", (promoter, f"Example promoter {index}"))
                connection.execute("INSERT INTO tu_genes(tu_id, gene_id) VALUES (?, ?)", (tu, gene))
                connection.execute("INSERT INTO tu_promoters(tu_id, promoter_id) VALUES (?, ?)", (tu, promoter))
                connection.execute("INSERT INTO promoter_genes(promoter_id, gene_id) VALUES (?, ?)", (promoter, gene))
                connection.execute("INSERT INTO contrasts(contrast_id, run_id, name, numerator, denominator, comparison_type) VALUES (?, 'demo-v1', ?, 'fictional exposure', 'fictional control', 'drug_control')", (contrast, f"Example comparison {index}"))
                connection.execute("INSERT INTO de_results(run_id, contrast_id, gene_id, log2_fold_change, padj) VALUES ('demo-v1', ?, ?, 2.5, 0.01)", (contrast, gene))
                connection.execute("INSERT INTO candidates(run_id, candidate_id, gene_id, promoter_id, antibiotic_class, support_tier, ranking_score) VALUES ('demo-v1', ?, ?, ?, ?, 'fictional_example', 2.5)", (gene, gene, promoter, antibiotic_class))
                connection.execute("INSERT INTO candidate_evidence(run_id, candidate_id, contrast_id, source_gene_id, effect, padj, evidence_role) VALUES ('demo-v1', ?, ?, ?, 2.5, 0.01, 'discovery')", (gene, contrast, gene))
                connection.execute("INSERT INTO candidate_tus VALUES ('demo-v1', ?, ?)", (gene, tu))
                connection.execute("INSERT INTO candidate_promoters VALUES ('demo-v1', ?, ?)", (gene, promoter))
                connection.execute("INSERT INTO candidate_panels(run_id, panel_id, candidate_id, promoter_id, antibiotic_class, selection_order, selection_rationale, transcription_group) VALUES ('demo-v1', 'experimental_panel', ?, ?, ?, ?, 'Fictional demonstration only', ?)", (gene, promoter, antibiotic_class, index, tu))
                connection.execute("INSERT INTO regulatory_interactions(interaction_id, run_id, actor_id, interaction_type, target_kind, promoter_id, effect, interaction_evidence, source_release, source_json) VALUES (?, 'demo-v1', 'DEMO-R1', 'TF-promoter', 'promoter', ?, ?, 'Fictional example', 'synthetic-demo', '{}')", (f"demo-interaction-{index}", promoter, "activator" if index == 1 else "repressor"))
            if connection.execute("PRAGMA foreign_key_check").fetchone():
                raise ValueError("Demo has broken relationships")
        connection.execute("PRAGMA journal_mode = DELETE")
        close_database(connection)
        connection = None
        # Hard-link publication fails if another process created the destination.
        os.link(temporary, path)
    finally:
        if connection is not None:
            connection.close()
        temporary.unlink(missing_ok=True)
        Path(f"{temporary}-wal").unlink(missing_ok=True)
        Path(f"{temporary}-shm").unlink(missing_ok=True)
    return path


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--database", type=Path, default=DEFAULT_DEMO)
    args = parser.parse_args(argv)
    try:
        path = build_demo(args.database)
    except (OSError, ValueError) as error:
        parser.error(str(error))
    print(f"Fictional demo: {path}")
    print(f"python -m promoter_discovery.explorer --database {path}")


if __name__ == "__main__":
    main()
