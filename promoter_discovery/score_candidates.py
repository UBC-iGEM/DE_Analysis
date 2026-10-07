"""Score typed regulatory edges and annotate candidate evidence.

Only curated regulatory edges count as regulation. Candidate coverage is
descriptive, not antibiotic specificity. Primary ranking uses measured DE
support; PRECISE annotations do not alter it.
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
    from .dataset_registry import CLASS_REGISTRY
    from .regulon_enrichment import compute_regulon_enrichment
except ImportError:  # pragma: no cover
    from build_network import REGULATORY_EDGE_TYPES  # type: ignore
    from dataset_registry import CLASS_REGISTRY  # type: ignore
    from regulon_enrichment import compute_regulon_enrichment  # type: ignore


CLASS_KEYS = tuple(CLASS_REGISTRY)
TIER_ORDER = {"multiple_comparisons_supported": 0, "single_comparison_supported": 1, "within_class_direction_difference": 2}


def load_graph(path: str | Path) -> nx.DiGraph:
    with Path(path).open("rb") as handle:
        return pickle.load(handle)


def _is_regulator(data: dict[str, Any]) -> bool:
    return data.get("node_type") == "regulator" or bool(data.get("is_regulator"))


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


def compute_regulon_coverage(graph: nx.DiGraph) -> pd.DataFrame:
    """Describe candidate coverage of full regulons; this is not specificity."""

    rows: list[dict[str, Any]] = []
    for regulator, data in graph.nodes(data=True):
        if not _is_regulator(data):
            continue
        targets: list[str] = []
        class_counts = {class_key: 0 for class_key in CLASS_KEYS}
        for target in graph.successors(regulator):
            edge = graph.edges[regulator, target]
            if edge.get("edge_type") not in REGULATORY_EDGE_TYPES:
                continue
            target_classes = _target_classes(graph.nodes[target])
            targets.append(target)
            for class_key in target_classes:
                class_counts[class_key] += 1
        if not targets:
            continue
        total = len(targets)
        row = {
            "regulator": regulator,
            "regulator_type": data.get("regulator_type", "unknown"),
            "n_total_targets": total,
            "target_genes": ", ".join(sorted(set(targets))),
        }
        for class_key in CLASS_KEYS:
            row[f"n_{class_key}_targets"] = class_counts[class_key]
            row[f"candidate_coverage_{class_key}"] = round(class_counts[class_key] / total, 3)
        nonzero_classes = [class_key for class_key in CLASS_KEYS if class_counts[class_key] > 0]
        max_count = max(class_counts.values())
        dominant = [class_key for class_key in CLASS_KEYS if class_counts[class_key] == max_count]
        row["has_candidates_in_multiple_classes"] = len(nonzero_classes) > 1
        row["dominant_class"] = dominant[0] if len(dominant) == 1 else "mixed"
        rows.append(row)
    count_columns = [f"n_{class_key}_targets" for class_key in CLASS_KEYS]
    coverage_columns = [f"candidate_coverage_{class_key}" for class_key in CLASS_KEYS]
    columns = ["regulator", "regulator_type", *count_columns, "n_total_targets", *coverage_columns,
               "has_candidates_in_multiple_classes", "dominant_class", "target_genes"]
    return pd.DataFrame(rows, columns=columns).sort_values("n_total_targets", ascending=False, ignore_index=True)


def annotate_candidates(graph: nx.DiGraph, enrichment: pd.DataFrame | None = None) -> pd.DataFrame:
    """Attach typed regulatory context without treating co-membership as regulation."""

    rows: list[dict[str, Any]] = []
    for node, data in graph.nodes(data=True):
        if data.get("node_type") != "candidate":
            continue
        regulators = [
            predecessor for predecessor in graph.predecessors(node)
            if _is_regulator(graph.nodes[predecessor])
            and graph.edges[predecessor, node].get("edge_type") in REGULATORY_EDGE_TYPES
        ]
        profile = json.loads(data.get("dataset_evidence_json", "[]"))
        detected = [entry for entry in profile if entry.get("significant")]
        response_classes = sorted({entry["antibiotic_class"] for entry in detected})
        profile_status = "responses_in_multiple_classes" if len(response_classes) > 1 else "response_detected_in_one_class" if response_classes else "response_profile_unavailable"
        regulatory_evidence = [] if enrichment is None or enrichment.empty else enrichment[enrichment.regulator.isin(regulators)].to_dict("records")
        # A non-significant challenge is not evidence of equivalence or absence.
        row = {
            "gene": node,
            "regulondb_release": graph.graph.get("regulondb_release", "unknown"),
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
            "primary_max_abs_log2_fold_change": data.get("primary_max_abs_log2_fold_change", data.get("max_abs_log2_fold_change")),
            "regulators": ", ".join(sorted(regulators)),
            "response_profile_json": json.dumps(profile, sort_keys=True),
            "response_classes": ", ".join(response_classes),
            "response_profile_status": profile_status,
            "regulon_evidence_json": json.dumps(regulatory_evidence, sort_keys=True),
            "cross_class_direction_difference": data.get("cross_class_direction_difference", False),
            "independent_validation": False,
            "biosensor_flag": "requires_promoter_mapping_and_wetlab_validation",
            "primary_ranking_basis": "DE support tier, then absolute fold change; no calibrated performance prediction",
            "external_annotation_role": "PRECISE expression/activity and burden weights are supporting context only",
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
        return pd.DataFrame(columns=["gene", "evidence_tier", "regulators", "response_profile_json", "biosensor_flag"])
    frame["_tier_order"] = frame["evidence_tier"].map(TIER_ORDER).fillna(99)
    return frame.sort_values(["_tier_order", "primary_max_abs_log2_fold_change", "gene"], ascending=[True, False, True], na_position="last").drop(columns="_tier_order").reset_index(drop=True)


def print_summary(tf_scores: pd.DataFrame, candidate_scores: pd.DataFrame, out_dir: str | Path) -> None:
    lines = ["REGULATORY NETWORK SCORING SUMMARY", f"regulators scored: {len(tf_scores)}", f"candidates annotated: {len(candidate_scores)}"]
    if not candidate_scores.empty:
        lines.append("evidence tiers: " + candidate_scores["evidence_tier"].value_counts().to_string())
        lines.append("Support labels describe comparisons, not independent validation. Non-significance does not establish specificity.")
    (Path(out_dir) / "scoring_summary.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))


def class_shortlist(candidate_scores, graph, class_key):
    keep = [class_key in graph.nodes[gene].get("significant_classes", []) for gene in candidate_scores.get("gene", [])]
    shortlist = candidate_scores.loc[keep].copy()
    shortlist["shortlist_class"] = class_key
    shortlist["class_support_tier"] = shortlist.gene.map(lambda gene: graph.nodes[gene].get("evidence_tier_by_class", {}).get(class_key, "unknown"))
    def class_effect(gene):
        profile = json.loads(graph.nodes[gene].get("dataset_evidence_json", "[]"))
        values = [abs(value) for entry in profile if entry["antibiotic_class"] == class_key
                  and entry.get("qualifying") and entry.get("dataset_role") != "cross_reactivity_challenge"
                  for value in entry.get("log2FoldChange", [])]
        return max(values, default=float("nan"))
    shortlist["class_max_abs_log2_fold_change"] = shortlist.gene.map(class_effect)
    shortlist["_class_tier"] = shortlist.class_support_tier.map(TIER_ORDER).fillna(99)
    return shortlist.sort_values(["_class_tier", "class_max_abs_log2_fold_change", "gene"], ascending=[True, False, True]).drop(columns="_class_tier")


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--graph", default=root / "results" / "regulatory_network" / "regulatory_network.pkl")
    parser.add_argument("--out", default=root / "results", help="Results root; writes regulatory_network/ and promoter_candidates/ stages")
    parser.add_argument("--pairwise", type=Path, default=root / "results/differential_expression/model/drug_pairwise_contrasts.csv",
                        help="Direct drug contrasts from the joint benchmark fit; optional for non-benchmark runs")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv)
    graph = load_graph(args.graph)
    if graph.graph.get("network_scope") != "full_reference":
        raise ValueError("Legacy candidate-only graph cannot support full-regulon statistics; rerun build_network")
    tf_scores = compute_regulon_coverage(graph)
    enrichment = compute_regulon_enrichment(graph)
    candidate_scores = annotate_candidates(graph, enrichment)
    output = Path(args.out) / "regulatory_network"
    output.mkdir(parents=True, exist_ok=True)
    candidates_output = Path(args.out) / "promoter_candidates"
    candidates_output.mkdir(parents=True, exist_ok=True)
    config = graph.graph.get("dataset_config", {})
    if config.get("benchmark"):
        if not args.pairwise.exists():
            raise FileNotFoundError("Benchmark direct drug contrasts missing; run promoter_discovery.run_benchmark first")
        try:
            from .setup_data import sha256_file
        except ImportError:  # pragma: no cover
            from setup_data import sha256_file
        if graph.graph.get("benchmark_provenance", {}).get("pairwise_sha256") != sha256_file(args.pairwise):
            raise ValueError("Direct drug contrasts do not match the graph's joint fit; rebuild the benchmark and network")
        pairs = pd.read_csv(args.pairwise)
        required = {"gene_id", "numerator_dataset", "denominator_dataset", "log2FoldChange", "lfcSE", "padj", "lfc_ci_low", "lfc_ci_high"}
        if required - set(pairs):
            raise ValueError("Direct drug contrasts are missing required fields")
        observations = pd.DataFrame(graph.graph["source_observations"])
        lookup = observations.drop_duplicates("source_gene_id").set_index("source_gene_id").canonical_gene
        pairs["gene"] = pairs.gene_id.map(lookup)
        pairs = pairs[pairs.gene.isin(candidate_scores.get("gene", []))].copy()
        pairs.to_csv(output / "candidate_drug_pairwise_contrasts.csv", index=False)
        pair_lookup = {gene: group.to_dict("records") for gene, group in pairs.groupby("gene")}
        candidate_scores["direct_drug_contrasts_json"] = candidate_scores.gene.map(lambda gene: json.dumps(pair_lookup.get(gene, []), sort_keys=True))
    enrichment.to_csv(output / "regulon_enrichment.csv", index=False)
    tf_scores.to_csv(output / "regulon_candidate_coverage.csv", index=False)
    stale = output / "tf_specificity_scores.csv"
    if stale.exists():
        stale.unlink()
    candidate_scores.to_csv(candidates_output / "annotated_candidates.csv", index=False)
    for class_key in ("beta_lactam", "aminoglycoside"):
        class_shortlist(candidate_scores, graph, class_key).to_csv(candidates_output / f"{class_key}_candidates.csv", index=False)
    provenance = {"regulondb_release": graph.graph.get("regulondb_release"),
                  "reference_provenance": graph.graph.get("reference_provenance", []),
                  "benchmark_provenance": graph.graph.get("benchmark_provenance"),
                  "enrichment_background": graph.graph.get("enrichment_background"),
                  "enrichment_test": "one-sided Fisher exact; BH across all regulators x contrasts x up/down directions",
                  "n_enrichment_tests": len(enrichment), "n_candidates": len(candidate_scores),
                  "interpretation": "Exploratory gene-level enrichment; operon correlation and annotation coverage limit interpretation. No calibrated biosensor recommendation."}
    (output / "scoring_provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    print_summary(tf_scores, candidate_scores, output)


if __name__ == "__main__":
    main()
