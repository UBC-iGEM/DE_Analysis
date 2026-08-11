"""Score typed regulatory edges and annotate candidate evidence.

Only edges whose ``edge_type`` is ``activates``, ``represses``, or ``dual`` are
counted as regulation.  ``co-imodulon`` edges are retained for context but can
never leak into regulator specificity or candidate regulator lists.
"""

from __future__ import annotations

import argparse
import json
import pickle
from pathlib import Path
from typing import Any

import networkx as nx
import pandas as pd

try:
    from .build_network import REGULATORY_EDGE_TYPES
except ImportError:  # pragma: no cover
    from build_network import REGULATORY_EDGE_TYPES  # type: ignore


CLASS_KEYS = ("beta_lactam", "aminoglycoside")
TIER_ORDER = {"conflicted": 0, "limited": 1, "supported": 2, "corroborated": 3}


def load_graph(path: str | Path) -> nx.DiGraph:
    with Path(path).open("rb") as handle:
        return pickle.load(handle)


def _target_classes(data: dict[str, Any]) -> set[str]:
    classes = data.get("significant_classes")
    if isinstance(classes, str):
        try:
            classes = json.loads(classes)
        except json.JSONDecodeError:
            classes = [classes]
    if isinstance(classes, (list, tuple, set)):
        return set(classes) & set(CLASS_KEYS)
    group = data.get("group")
    return {group} if group in CLASS_KEYS else set(CLASS_KEYS) if group == "cross" else set()


def compute_tf_specificity(graph: nx.DiGraph) -> pd.DataFrame:
    """Count only typed regulatory edges for each regulator."""

    rows: list[dict[str, Any]] = []
    for regulator, data in graph.nodes(data=True):
        if data.get("node_type") != "regulator":
            continue
        targets: list[str] = []
        class_counts = {class_key: 0 for class_key in CLASS_KEYS}
        for target in graph.successors(regulator):
            edge = graph.edges[regulator, target]
            if edge.get("edge_type") not in REGULATORY_EDGE_TYPES:
                continue
            target_classes = _target_classes(graph.nodes[target])
            if not target_classes:
                continue
            targets.append(target)
            for class_key in target_classes:
                class_counts[class_key] += 1
        if not targets:
            continue
        total = len(targets)
        beta_count = class_counts["beta_lactam"]
        amino_count = class_counts["aminoglycoside"]
        rows.append({
            "regulator": regulator,
            "regulator_type": data.get("regulator_type", "unknown"),
            "n_beta_lactam_targets": beta_count,
            "n_aminoglycoside_targets": amino_count,
            "n_total_targets": total,
            "specificity_beta_lactam": round(beta_count / total, 3),
            "specificity_aminoglycoside": round(amino_count / total, 3),
            "is_cross_reactive": beta_count > 0 and amino_count > 0,
            "dominant_class": "beta_lactam" if beta_count > amino_count else "aminoglycoside" if amino_count > beta_count else "mixed",
            "target_genes": ", ".join(sorted(set(targets))),
        })
    columns = ["regulator", "regulator_type", "n_beta_lactam_targets", "n_aminoglycoside_targets",
               "n_total_targets", "specificity_beta_lactam", "specificity_aminoglycoside",
               "is_cross_reactive", "dominant_class", "target_genes"]
    return pd.DataFrame(rows, columns=columns).sort_values("n_total_targets", ascending=False, ignore_index=True)


def annotate_candidates(graph: nx.DiGraph, tf_scores: pd.DataFrame) -> pd.DataFrame:
    """Attach typed regulatory context without treating co-membership as regulation."""

    score_map = {row["regulator"]: row for row in tf_scores.to_dict("records")}
    rows: list[dict[str, Any]] = []
    for node, data in graph.nodes(data=True):
        if data.get("node_type") != "candidate":
            continue
        regulators = [
            predecessor for predecessor in graph.predecessors(node)
            if graph.nodes[predecessor].get("node_type") == "regulator"
            and graph.edges[predecessor, node].get("edge_type") in REGULATORY_EDGE_TYPES
        ]
        relevant_classes = _target_classes(data)
        specs: list[tuple[str, float, str]] = []
        cross_flags: list[bool] = []
        for regulator in regulators:
            score = score_map.get(regulator)
            if score is None:
                continue
            values = [score.get(f"specificity_{class_key}") for class_key in relevant_classes]
            values = [float(value) for value in values if value is not None]
            if values:
                specs.append((regulator, max(values), score.get("regulator_type", "unknown")))
            cross_flags.append(bool(score.get("is_cross_reactive")))
        priority = {"TF": 2, "sigma": 2, "effector": 0, "unknown": 1}
        best = max(specs, key=lambda item: (item[1], priority.get(item[2], 0))) if specs else ("", None, "")
        clarity = None if best[1] is None else round(float(best[1]) * (0.6 if any(cross_flags) else 1.0), 3)
        row = {
            "gene": node,
            "group": data.get("group", ""),
            "evidence_tier": data.get("evidence_tier", ""),
            "evidence_tier_by_class": json.dumps(data.get("evidence_tier_by_class", {}), sort_keys=True, default=str),
            "n_datasets_observed": data.get("n_datasets_observed", 0),
            "n_datasets_significant": data.get("n_datasets_significant", 0),
            "direction_consistent": data.get("direction_consistent"),
            "candidate_direction_policy": data.get("candidate_direction_policy", ""),
            "datasets_significant_by_class": json.dumps(data.get("datasets_significant_by_class", {}), sort_keys=True, default=str),
            "support_fraction_by_class": json.dumps(data.get("support_fraction_by_class", {}), sort_keys=True, default=str),
            "dataset_evidence_json": data.get("dataset_evidence_json", "[]"),
            "source_observations_json": data.get("source_observations_json", "[]"),
            "caveats": "; ".join(data.get("caveats", [])) if isinstance(data.get("caveats"), list) else data.get("caveats", ""),
            "evidence_quality_flags": "; ".join(data.get("evidence_quality_flags", [])) if isinstance(data.get("evidence_quality_flags"), list) else data.get("evidence_quality_flags", ""),
            "tobramycin_only": data.get("tobramycin_only", False),
            "max_abs_log2_fold_change": data.get("max_abs_log2_fold_change"),
            "regulators": ", ".join(sorted(regulators)),
            "best_regulator": best[0],
            "best_regulator_type": best[2],
            "best_regulator_specificity": best[1],
            "has_cross_reactive_regulator": any(cross_flags),
            "regulatory_clarity": clarity,
            "biosensor_flag": "?" if clarity is None else "recommended" if clarity >= 0.75 else "review" if clarity >= 0.5 else "weak",
            "gene_expression_json": json.dumps(data.get("gene_expression", {}), sort_keys=True, default=str),
            "gene_expression_basal": json.dumps(data.get("gene_expression_basal", {}), sort_keys=True, default=str),
            "gene_expression_induced": json.dumps(data.get("gene_expression_induced", {}), sort_keys=True, default=str),
            "gene_expression_delta": json.dumps(data.get("gene_expression_delta", {}), sort_keys=True, default=str),
            "imodulon_activity_json": json.dumps(data.get("imodulon_activity", {}), sort_keys=True, default=str),
            "imodulon_activity_basal": json.dumps(data.get("imodulon_activity_basal", {}), sort_keys=True, default=str),
            "imodulon_activity_induced": json.dumps(data.get("imodulon_activity_induced", {}), sort_keys=True, default=str),
            "imodulon_activity_delta": json.dumps(data.get("imodulon_activity_delta", {}), sort_keys=True, default=str),
            "metabolic_burden_proxy": data.get("metabolic_burden_proxy"),
            "translation_burden_proxy": data.get("translation_burden_proxy"),
            "burden_proxy_system_categories": ", ".join(data.get("burden_proxy_system_categories", [])),
            "burden_proxy_imodulons": ", ".join(data.get("burden_proxy_imodulons", [])),
            "burden_proxy_basis": data.get("burden_proxy_basis", ""),
            "note": "no regulator in RegulonDB" if not regulators else "",
        }
        rows.append(row)
    frame = pd.DataFrame(rows)
    if frame.empty:
        return frame
    frame["_tier_order"] = frame["evidence_tier"].map(TIER_ORDER).fillna(99)
    return frame.sort_values(["_tier_order", "regulatory_clarity", "max_abs_log2_fold_change"], ascending=[True, False, False], na_position="last").drop(columns="_tier_order").reset_index(drop=True)


def print_summary(tf_scores: pd.DataFrame, candidate_scores: pd.DataFrame, out_dir: str | Path) -> None:
    lines = ["REGULATORY NETWORK SCORING SUMMARY", f"regulators scored: {len(tf_scores)}", f"candidates annotated: {len(candidate_scores)}"]
    if not candidate_scores.empty:
        lines.append("evidence tiers: " + candidate_scores["evidence_tier"].value_counts().to_string())
        lines.append("Tobramycin-only candidates are limited evidence and require follow-up.")
    (Path(out_dir) / "scoring_summary.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--graph", default=root / "network_analysis" / "output" / "regulatory_network.pkl")
    parser.add_argument("--out", default=root / "network_analysis" / "output")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv)
    graph = load_graph(args.graph)
    tf_scores = compute_tf_specificity(graph)
    candidate_scores = annotate_candidates(graph, tf_scores)
    output = Path(args.out)
    output.mkdir(parents=True, exist_ok=True)
    tf_scores.to_csv(output / "tf_specificity_scores.csv", index=False)
    candidate_scores.to_csv(output / "annotated_candidates.csv", index=False)
    print_summary(tf_scores, candidate_scores, output)


if __name__ == "__main__":
    main()
