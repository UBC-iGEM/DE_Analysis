"""Descriptive comparison of panel selection rules on the same observed data.

This is an internal design comparison, not an independent biological benchmark.
Both selectors use the same mapped, sequence-matching pool and class budget.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import pandas as pd


CLASSES = ("beta_lactam", "aminoglycoside")


def effect_only_panel(options: pd.DataFrame, per_class: int) -> pd.DataFrame:
    """Rank by shrunken induction; break ties by annotation confidence and IDs.

    Choose distinct genes across classes, but impose no promoter/TU diversity
    criterion. Missing effects and sequence mismatches are ineligible.
    """
    chosen, genes = [], set()
    for _ in range(per_class):
        for cls in CLASSES:
            available = options[options.antibiotic_class.eq(cls) & options.reference_sequence_matches.eq(True)
                                & options.shrunk_effect.notna() & ~options.gene.isin(genes)]
            if available.empty:
                continue
            row = available.sort_values(
                ["shrunk_effect", "confidence_rank", "gene", "promoter_id"],
                ascending=[False, False, True, True],
            ).iloc[0].to_dict()
            chosen.append(row)
            genes.add(row["gene"])
    return pd.DataFrame(chosen, columns=options.columns)


def panel_metrics(panel: pd.DataFrame) -> dict:
    def median(column):
        values = panel[column].dropna()
        return float(values.median()) if len(values) else None
    return {
        "selected_slots": len(panel), "unique_genes": int(panel.gene.nunique()),
        "unique_promoters": int(panel.promoter_id.nunique()),
        "unique_tu_groups": int(panel.transcription_group.nunique()),
        "class_counts": class_counts(panel),
        "median_shrunk_log2_effect": median("shrunk_effect"),
        "median_screen_fraction": median("screen_fraction"),
        "median_max_counterpart_abs_log2_effect": median("challenge_max_abs_log2FC"),
        "missing_counterpart_effects": int(panel.challenge_max_abs_log2FC.isna().sum()),
    }


def class_counts(panel: pd.DataFrame) -> dict:
    """Count selected rows for each supported antibiotic class."""
    counts = {}
    for antibiotic_class in CLASSES:
        counts[antibiotic_class] = int(panel.antibiotic_class.eq(antibiotic_class).sum())
    return counts


def compare_selection(options: pd.DataFrame, panel: pd.DataFrame, per_class: int) -> tuple[dict, pd.DataFrame]:
    if type(per_class) is not int or per_class < 1:
        raise ValueError("per_class must be a positive integer")
    if not panel.antibiotic_class.isin(CLASSES).all() or any(panel.antibiotic_class.eq(cls).sum() > per_class for cls in CLASSES):
        raise ValueError("Existing panel exceeds the requested class budget; use its original per-class budget")
    keys = ["gene", "promoter_id", "antibiotic_class"]
    eligible = options[options.reference_sequence_matches.eq(True)].set_index(keys)
    if not eligible.index.is_unique:
        raise ValueError("Panel options contain duplicate gene/promoter/class records")
    for row in panel.to_dict("records"):
        key = tuple(row[name] for name in keys)
        if key not in eligible.index:
            raise ValueError(f"Selected panel entry is outside the current eligible pool: {key}")
        for field in ("shrunk_effect", "screen_fraction", "challenge_max_abs_log2FC"):
            current, stored = eligible.loc[key, field], row[field]
            if pd.isna(current) and pd.isna(stored):
                continue
            if pd.isna(current) or pd.isna(stored) or not math.isclose(float(current), float(stored), rel_tol=1e-9, abs_tol=1e-12):
                raise ValueError(f"Panel metrics disagree with current source artifacts: {key}, {field}")
    baseline = effect_only_panel(options, per_class)
    report = {
        "comparison_type": "descriptive_same_data_design_comparison",
        "requested_per_class": per_class,
        "eligible_gene_promoter_class_options": int(options.reference_sequence_matches.eq(True).sum()),
        "baseline_rule": "Descending shrunken induction; confidence then gene/promoter IDs break ties; unique genes across classes; no promoter/TU diversity constraint",
        "shared_eligibility": "Mapped candidate-promoter pairs with matching reference sequences; same per-class budget",
        "prism": panel_metrics(panel), "effect_only": panel_metrics(baseline),
        "limitations": [
            "Uses the discovery data and annotation pool; no held-out validation or reporter measurements.",
            "Diversity is a PRISM selection objective, so this comparison demonstrates that constraint rather than independent predictive superiority.",
            "Counterpart effects summarize observed drugs, not specificity against all stresses or the exact wet-lab targets.",
            "Unmapped genes are excluded from both selectors; this comparison cannot remove annotation bias.",
        ],
    }
    return report, baseline


def evaluate(root: Path, per_class: int = 6) -> dict:
    from .assess_candidates import panel_options, transcription_groups
    from .promoter_selection import read_product
    from .provenance import source_provenance

    assessment = root / "results/candidate_assessment"
    files = {
        "stability": assessment / "candidate_stability.csv",
        "candidates": root / "results/promoter_candidates/annotated_candidates.csv",
        "mapped": root / "results/promoter_review/candidate_promoter_mapping.csv",
        "fragments": assessment / "fragment_review.csv",
        "reactivity": assessment / "cross_reactivity.csv",
        "replicates": assessment / "candidate_replicates.csv",
        "panel": assessment / "experimental_panel.csv",
        "tus": root / "data/references/regulondb/TUSet.tsv",
    }
    frames = {}
    for key, path in files.items():
        if key != "tus":
            frames[key] = pd.read_csv(path)
    options = panel_options(*(frames[key] for key in ("stability", "candidates", "mapped", "fragments", "reactivity", "replicates")),
                            transcription_groups(read_product(files["tus"])))
    report, baseline = compare_selection(options, frames["panel"], per_class)
    input_hashes = {}
    for path in files.values():
        relative_path = path.relative_to(root).as_posix()
        input_hashes[relative_path] = hashlib.sha256(path.read_bytes()).hexdigest()
    report["input_sha256"] = input_hashes
    report["code_sha256"] = source_provenance()["code_sha256"]
    baseline.to_csv(assessment / "selection_baseline.csv", index=False)
    (assessment / "selection_comparison.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    return report


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--per-class", type=int, default=6)
    args = parser.parse_args(argv)
    print(json.dumps(evaluate(args.root.resolve(), args.per_class), indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
