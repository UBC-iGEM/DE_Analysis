"""Annotate a regulatory graph with PRECISE/iModulon evidence.

The module keeps two evidence families separate:

* candidate-gene expression from ``IcaData.log_tpm``, ``IcaData.X``, or a
  companion PRECISE expression matrix;
* primary iModulon activity from ``IcaData.A``.

Both retain basal, induced, delta, per-background counts, units, and
availability.  Metabolism and Translation are reported as burden proxies, not
as measurements of cellular burden.
"""

from __future__ import annotations

import argparse
import gzip
import html
import json
import os
import pickle
from dataclasses import dataclass
from io import StringIO
from itertools import combinations
from pathlib import Path
from typing import Any

import networkx as nx
import numpy as np
import pandas as pd

from pymodulon.core import IcaData
from pymodulon.enrichment import compute_trn_enrichment

try:
    from .dataset_registry import CLASS_REGISTRY, DATASET_CAVEATS
    from .build_network import REGULATORY_EDGE_TYPES
except ImportError:  # pragma: no cover
    from dataset_registry import CLASS_REGISTRY, DATASET_CAVEATS  # type: ignore
    from build_network import REGULATORY_EDGE_TYPES  # type: ignore


CONDITION_KEYWORDS: dict[str, list[str]] = {
    "beta_lactam": ["cef", "mero", "imipenem", "ampicillin", "amoxicillin", "penicillin"],
    "aminoglycoside": ["kan", "tobramycin", "gentamicin", "amikacin", "streptomycin"],
}
CONTROL_TOKEN = "ctrl"
BURDEN_RISK_SYSTEM_CATEGORIES = {"Metabolism", "Translation"}


@dataclass
class GeneLookup:
    by_name: dict[str, str]
    by_id: dict[str, str]


def _read_json_table(value: Any) -> pd.DataFrame | None:
    if isinstance(value, pd.DataFrame):
        return value
    if value is None:
        return None
    return pd.read_json(StringIO(value))


def load_precise_model(path: str | Path) -> IcaData:
    """Load PRECISE JSON.GZ with the pymodulon 0.2.1/pandas compatibility fix."""

    with gzip.open(path, "rt", encoding="utf-8") as handle:
        serial = json.load(handle)
    cutoff_optimized = serial.pop("_cutoff_optimized", False)
    dagostino_cutoff = serial.pop("_dagostino_cutoff", None)
    for key in ("M", "A", "X", "log_tpm", "gene_table", "sample_table", "imodulon_table", "trn"):
        if key in serial and serial[key] is not None:
            serial[key] = _read_json_table(serial[key])
    original_astype = pd.Index.astype

    def patched_astype(index: pd.Index, dtype: Any, **kwargs: Any) -> pd.Index:
        try:
            return original_astype(index, dtype, **kwargs)
        except (TypeError, ValueError):
            return index

    pd.Index.astype = patched_astype  # type: ignore[assignment]
    try:
        model = IcaData(**serial)
    finally:
        pd.Index.astype = original_astype  # type: ignore[assignment]
    model._cutoff_optimized = cutoff_optimized
    model._dagostino_cutoff = dagostino_cutoff
    return model


def load_gene_mapping(path: str | Path | None) -> pd.DataFrame | None:
    if not path:
        return None
    frame = pd.read_csv(path, sep=None, engine="python")
    frame.columns = [str(c).strip().lower() for c in frame.columns]
    return frame


def _canonical(value: Any) -> str:
    if value is None or pd.isna(value):
        return ""
    return str(value).strip().lower()


def _split_aliases(value: Any) -> list[str]:
    text = _canonical(value)
    return [part.strip() for part in text.replace("|", ";").replace(",", ";").split(";") if part.strip()]


def build_gene_lookup(ica: IcaData, mapping: pd.DataFrame | None = None) -> GeneLookup:
    by_name: dict[str, str] = {}
    by_id: dict[str, str] = {}
    if mapping is not None and not mapping.empty:
        gene_cols = [c for c in ("canonical_gene", "gene", "gene_name", "symbol") if c in mapping.columns]
        id_cols = [c for c in ("canonical_locus_tag", "locus_tag", "b_number", "gene_id") if c in mapping.columns]
        for _, row in mapping.iterrows():
            gene = _canonical(row[gene_cols[0]]) if gene_cols else ""
            gene_id = _canonical(row[id_cols[0]]) if id_cols else ""
            if gene and gene_id:
                by_name.setdefault(gene, gene_id)
                by_id.setdefault(gene_id, gene_id)
    table = ica.gene_table.copy()
    table.index = table.index.astype(str).str.lower()
    for gene_id, row in table.iterrows():
        by_id.setdefault(gene_id, gene_id)
        if "gene_name" in table.columns:
            name = _canonical(row.get("gene_name"))
            if name:
                by_name.setdefault(name, gene_id)
        for alias in _split_aliases(row.get("synonyms")):
            by_name.setdefault(alias, gene_id)
    return GeneLookup(by_name=by_name, by_id=by_id)


def resolve_gene_id(name: str, lookup: GeneLookup, ica: IcaData) -> str | None:
    query = _canonical(name)
    if query in lookup.by_id:
        return lookup.by_id[query]
    if query in lookup.by_name:
        return lookup.by_name[query]
    index = {str(value).lower(): str(value) for value in ica.gene_table.index}
    return index.get(query)


def build_membership_tables(ica: IcaData) -> tuple[pd.DataFrame, dict[str, set[str]]]:
    thresholds = pd.Series(ica.thresholds, dtype=float).reindex(ica.M.columns).fillna(0.0)
    membership = ica.M.abs().gt(thresholds, axis=1)
    return membership, {imod: set(membership.index[membership[imod]]) for imod in membership.columns}


def build_condition_groups(ica: IcaData) -> dict[str, dict[str, Any]]:
    """Pair ``{media}_{drug}`` samples to the matching ``{media}_ctrl``."""

    result = {key: {"treated": [], "control": [], "pairs": [], "backgrounds": []} for key in CONDITION_KEYWORDS}
    sample_table = ica.sample_table
    if sample_table is None or "condition" not in sample_table.columns:
        return result
    conditions = sample_table["condition"].astype(str)
    for class_key, keywords in CONDITION_KEYWORDS.items():
        by_background: dict[str, dict[str, list[str]]] = {}
        for sample_id, condition in conditions.items():
            condition_lower = condition.lower()
            tokens = condition_lower.split("_")
            drug = next((token for token in tokens if any(keyword in token for keyword in keywords)), None)
            if not drug:
                continue
            background = condition_lower.replace(f"_{drug}", "").strip("_") or tokens[0]
            control_condition = f"{background}_{CONTROL_TOKEN}"
            controls = conditions.index[conditions.str.lower() == control_condition].tolist()
            if not controls:
                continue
            slot = by_background.setdefault(background, {"treated": [], "control": []})
            slot["treated"].append(str(sample_id))
            slot["control"].extend(map(str, controls))
        for background, slot in sorted(by_background.items()):
            treated = sorted(set(slot["treated"]))
            control = sorted(set(slot["control"]))
            result[class_key]["treated"].extend(treated)
            result[class_key]["control"].extend(control)
            result[class_key]["pairs"].extend((treated_id, control_id) for treated_id in treated for control_id in control)
            result[class_key]["backgrounds"].append({"background": background, "treated": treated, "control": control})
    return result


def _paired_summary(frame: pd.DataFrame, groups: dict[str, Any], units: str, normalization: str) -> dict[str, Any]:
    backgrounds: list[dict[str, Any]] = []
    for background in groups.get("backgrounds", []):
        treated = [sample for sample in background["treated"] if sample in frame.columns]
        control = [sample for sample in background["control"] if sample in frame.columns]
        if not treated or not control:
            backgrounds.append(
                {"background": background["background"], "basal": None, "induced": None, "delta": None,
                 "n_treated": len(treated), "n_control": len(control), "available": False,
                 "units": units, "normalization": normalization}
            )
            continue
        basal = float(frame[control].median(axis=1).iloc[0])
        induced = float(frame[treated].median(axis=1).iloc[0])
        backgrounds.append(
            {"background": background["background"], "basal": basal, "induced": induced,
             "delta": induced - basal, "n_treated": len(treated), "n_control": len(control),
             "available": True, "units": units, "normalization": normalization}
        )
    available = [item for item in backgrounds if item["available"]]
    if not available:
        return {"basal": None, "induced": None, "delta": None, "n_treated": 0, "n_control": 0,
                "available": False, "per_background": backgrounds, "units": units, "normalization": normalization}
    basal = float(np.median([item["basal"] for item in available]))
    induced = float(np.median([item["induced"] for item in available]))
    return {"basal": basal, "induced": induced, "delta": induced - basal,
            "n_treated": int(sum(item["n_treated"] for item in available)),
            "n_control": int(sum(item["n_control"] for item in available)),
            "available": True, "per_background": backgrounds, "units": units, "normalization": normalization}


def compute_dynamic_range(ica: IcaData, condition_groups: dict[str, dict[str, Any]]) -> dict[str, dict[str, dict[str, Any]]]:
    """Compute basal/induced/delta activity from ``ica.A``."""

    result: dict[str, dict[str, dict[str, Any]]] = {}
    activity = ica.A
    for imodulon in activity.index:
        result[str(imodulon)] = {}
        for class_key, groups in condition_groups.items():
            one_row = activity.loc[[imodulon]]
            result[str(imodulon)][class_key] = _paired_summary(one_row, groups, "iModulon activity", "PRECISE A")
    return result


def _prepare_expression_frame(frame: pd.DataFrame) -> pd.DataFrame:
    frame = frame.copy()
    frame.index = frame.index.astype(str).str.strip().str.lower()
    frame.columns = [str(column).strip() for column in frame.columns]
    return frame.apply(pd.to_numeric, errors="coerce")


def load_expression_matrix(
    path: str | Path,
    gene_id_column: str | None = None,
    units: str = "unknown",
    normalization: str = "provided",
) -> tuple[pd.DataFrame, dict[str, str]]:
    """Load the companion expression matrix (rows genes, columns PRECISE samples)."""

    frame = pd.read_csv(path, sep=None, engine="python")
    frame.columns = [str(column).strip() for column in frame.columns]
    id_column = gene_id_column or next((c for c in frame.columns if c.lower() in {"gene", "gene_id", "locus_tag", "b_number"}), frame.columns[0])
    frame = frame.set_index(id_column)
    frame = _prepare_expression_frame(frame)
    metadata = {"units": units, "normalization": normalization, "gene_id_column": id_column, "source": str(path)}
    return frame, metadata


def load_embedded_expression(ica: IcaData) -> tuple[pd.DataFrame | None, dict[str, str]]:
    """Prefer gene-level expression embedded in an IcaData object.

    ``A`` is intentionally excluded: it contains iModulon activities rather
    than gene-level expression.  ``X`` is accepted only as a fallback because
    its centered-expression semantics vary by dataset.
    """

    for attribute, units, normalization in (
        ("log_tpm", "log2(TPM)", "embedded log2(TPM)"),
        ("X", "unknown", "embedded centered expression; verify source semantics"),
    ):
        value = getattr(ica, attribute, None)
        if value is None:
            continue
        frame = value if isinstance(value, pd.DataFrame) else pd.DataFrame(value)
        if frame.empty:
            continue
        return _prepare_expression_frame(frame), {
            "units": units,
            "normalization": normalization,
            "source": f"IcaData.{attribute}",
            "gene_id_column": "index",
        }
    return None, {"units": "unknown", "normalization": "unavailable", "source": "unavailable", "gene_id_column": "index"}


def validate_expression_alignment(expression: pd.DataFrame, ica: IcaData) -> None:
    """Reject an expression matrix that cannot support the model's evidence."""

    model_samples = {str(value).strip() for value in ica.A.columns}
    expression_samples = {str(value).strip() for value in expression.columns}
    if not model_samples.intersection(expression_samples):
        raise ValueError("Expression matrix has no sample IDs overlapping IcaData.A")
    model_genes = {str(value).strip().lower() for value in ica.M.index}
    expression_genes = {str(value).strip().lower() for value in expression.index}
    if not model_genes.intersection(expression_genes):
        raise ValueError("Expression matrix has no gene IDs overlapping IcaData.M")


def compute_gene_expression_evidence(
    expression: pd.DataFrame,
    gene_id: str | None,
    condition_groups: dict[str, dict[str, Any]],
    units: str = "unknown",
    normalization: str = "provided",
) -> dict[str, dict[str, Any]]:
    """Compute matched-background expression evidence for one candidate gene."""

    result: dict[str, dict[str, Any]] = {}
    key = _canonical(gene_id)
    for class_key, groups in condition_groups.items():
        if not key or key not in expression.index:
            result[class_key] = {"basal": None, "induced": None, "delta": None, "n_treated": 0,
                                 "n_control": 0, "available": False, "per_background": [],
                                 "units": units, "normalization": normalization}
            continue
        result[class_key] = _paired_summary(expression.loc[[key]], groups, units, normalization)
    return result


def host_burden_lookup(imodulon_table: pd.DataFrame) -> dict[str, dict[str, str]]:
    result: dict[str, dict[str, str]] = {}
    for imodulon, row in imodulon_table.iterrows():
        system = str(row.get("system_category", ""))
        functional = str(row.get("functional_category", ""))
        result[str(imodulon)] = {"system_category": system, "functional_category": functional}
    return result


def build_burden_proxy_lookup(ica: IcaData, membership: pd.DataFrame) -> dict[str, dict[str, Any]]:
    """Aggregate every weighted membership into transparent burden proxies."""

    categories = host_burden_lookup(ica.imodulon_table)
    result: dict[str, dict[str, Any]] = {}
    for gene_id in membership.index:
        hits = [str(imod) for imod, included in membership.loc[gene_id].items() if bool(included)]
        weights = {imod: float(ica.M.at[gene_id, imod]) for imod in hits}
        metabolic = sum(abs(weights[imod]) for imod in hits if categories.get(imod, {}).get("system_category") == "Metabolism")
        translation = sum(abs(weights[imod]) for imod in hits if categories.get(imod, {}).get("system_category") == "Translation")
        result[str(gene_id)] = {
            "metabolic_burden_proxy": float(metabolic),
            "translation_burden_proxy": float(translation),
            "burden_proxy_system_categories": sorted({categories.get(imod, {}).get("system_category", "") for imod in hits if categories.get(imod, {}).get("system_category", "")}),
            "burden_proxy_imodulons": hits,
            "burden_proxy_basis": "sum(abs(M gene membership weights)) by PRECISE system_category; lower is preferable; heuristic",
        }
    return result


def summarize_imodulon_enrichment(ica: IcaData, members_by_imodulon: dict[str, set[str]], fdr: float, max_regs: int, method: str) -> tuple[pd.DataFrame, dict[str, dict[str, Any]]]:
    rows: list[pd.DataFrame] = []
    summary: dict[str, dict[str, Any]] = {}
    all_genes = set(map(str, ica.M.index))
    for imodulon, gene_set in members_by_imodulon.items():
        try:
            enrichment = compute_trn_enrichment(gene_set=gene_set, all_genes=all_genes, trn=ica.trn, max_regs=max_regs, fdr=fdr, method=method)
        except Exception:
            enrichment = pd.DataFrame()
        if enrichment.empty:
            summary[imodulon] = {"top_regulator": "", "qvalue": None, "precision": None, "recall": None, "f1score": None, "n_significant": 0}
            continue
        enrichment = enrichment.reset_index(names="regulon")
        enrichment.insert(0, "imodulon", imodulon)
        rows.append(enrichment)
        top = enrichment.sort_values(["qvalue", "f1score", "precision"], ascending=[True, False, False]).iloc[0]
        summary[imodulon] = {"top_regulator": top["regulon"], "qvalue": float(top["qvalue"]),
                             "precision": float(top["precision"]), "recall": float(top["recall"]),
                             "f1score": float(top["f1score"]), "n_significant": int(len(enrichment))}
    return (pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()), summary


def _format_weight_map(values: dict[str, float]) -> str:
    return "; ".join(f"{key}:{value:.4f}" for key, value in values.items())


def _format_list(values: list[str]) -> str:
    return ", ".join(values)


def _percentile(values: list[float], value: float | None) -> float | None:
    if value is None or not values:
        return None
    return round(float((np.asarray(values) <= value).mean() * 100), 2)


def annotate_graph(
    graph: nx.DiGraph,
    ica: IcaData,
    lookup: GeneLookup,
    membership: pd.DataFrame,
    imodulon_summary: dict[str, dict[str, Any]],
    condition_groups: dict[str, dict[str, Any]] | None = None,
    expression: pd.DataFrame | None = None,
    expression_units: str = "unknown",
    expression_normalization: str = "provided",
    add_discovered: bool = False,
) -> tuple[nx.DiGraph, pd.DataFrame, pd.DataFrame | None, int]:
    """Attach iModulon, expression, activity, burden, and co-membership evidence."""

    condition_groups = condition_groups or {key: {"backgrounds": []} for key in CLASS_REGISTRY}
    activity = compute_dynamic_range(ica, condition_groups)
    burden = build_burden_proxy_lookup(ica, membership)
    candidates = [node for node, data in graph.nodes(data=True) if data.get("node_type") == "candidate"]
    candidate_gene_ids: dict[str, str] = {}
    discovered_rows: list[dict[str, Any]] = []
    for node in candidates:
        data = graph.nodes[node]
        gene_id = resolve_gene_id(str(node), lookup, ica)
        if gene_id:
            candidate_gene_ids[str(node)] = gene_id
        data.update({"precise_gene_id": gene_id or "", "imodulon_mapping_status": "mapped" if gene_id else "not_found"})
        data.update({"imodulons": [], "imodulon_weights": {}, "imodulon_primary": "", "imodulon_primary_weight": None,
                     "imodulon_primary_abs_weight": None, "imodulon_primary_regulon": "", "imodulon_count": 0})
        data["gene_expression"] = {key: {"basal": None, "induced": None, "delta": None, "available": False, "per_background": [], "units": expression_units, "normalization": expression_normalization} for key in CLASS_REGISTRY}
        data["imodulon_activity"] = {key: activity.get("", {}).get(key) for key in CLASS_REGISTRY}
        data["metabolic_burden_proxy"] = None
        data["translation_burden_proxy"] = None
        data["burden_proxy_system_categories"] = []
        data["burden_proxy_imodulons"] = []
        data["burden_proxy_basis"] = ""
        if not gene_id or gene_id not in membership.index:
            continue
        hits = [str(imod) for imod, included in membership.loc[gene_id].items() if bool(included)]
        weights = {imod: float(ica.M.at[gene_id, imod]) for imod in hits}
        hits.sort(key=lambda item: abs(weights[item]), reverse=True)
        data["imodulons"] = hits
        data["imodulon_weights"] = weights
        data["imodulon_count"] = len(hits)
        if hits:
            primary = hits[0]
            data["imodulon_primary"] = primary
            data["imodulon_primary_weight"] = weights[primary]
            data["imodulon_primary_abs_weight"] = abs(weights[primary])
            primary_summary = imodulon_summary.get(primary, {})
            data["imodulon_primary_regulon"] = primary_summary.get("top_regulator", "")
            proxy = burden.get(gene_id, {})
            data.update(proxy)
            data["imodulon_activity"] = activity.get(primary, {key: None for key in CLASS_REGISTRY})
        if expression is not None:
            expression_id = gene_id if gene_id in expression.index else str(node).lower() if str(node).lower() in expression.index else gene_id
            data["gene_expression"] = compute_gene_expression_evidence(expression, expression_id, condition_groups, expression_units, expression_normalization)
        data["gene_expression_basal"] = {key: value.get("basal") for key, value in data["gene_expression"].items()}
        data["gene_expression_induced"] = {key: value.get("induced") for key, value in data["gene_expression"].items()}
        data["gene_expression_delta"] = {key: value.get("delta") for key, value in data["gene_expression"].items()}
        data["gene_expression_available"] = {key: bool(value.get("available")) for key, value in data["gene_expression"].items()}
        data["imodulon_activity_basal"] = {key: value.get("basal") if value else None for key, value in data["imodulon_activity"].items()}
        data["imodulon_activity_induced"] = {key: value.get("induced") if value else None for key, value in data["imodulon_activity"].items()}
        data["imodulon_activity_delta"] = {key: value.get("delta") if value else None for key, value in data["imodulon_activity"].items()}
        data["imodulon_activity_available"] = {key: bool(value and value.get("available")) for key, value in data["imodulon_activity"].items()}
        if add_discovered:
            for imodulon in hits:
                for discovered_id in membership.index[membership[imodulon]].tolist():
                    if discovered_id == gene_id or discovered_id in candidate_gene_ids.values():
                        continue
                    row = ica.gene_table.loc[discovered_id]
                    discovered_rows.append({"gene_id": discovered_id, "gene_name": row.get("gene_name", ""),
                                            "imodulon": imodulon, "weight": float(ica.M.at[discovered_id, imodulon]),
                                            "source_candidate": node})
    # Add within-class expression percentiles without hiding the raw values.
    for class_key in CLASS_REGISTRY:
        basal_values = [graph.nodes[node]["gene_expression"][class_key]["basal"] for node in candidates if graph.nodes[node]["gene_expression"][class_key].get("available")]
        induced_values = [graph.nodes[node]["gene_expression"][class_key]["induced"] for node in candidates if graph.nodes[node]["gene_expression"][class_key].get("available")]
        for node in candidates:
            evidence = graph.nodes[node]["gene_expression"][class_key]
            evidence["basal_percentile"] = _percentile(basal_values, evidence.get("basal"))
            evidence["induced_percentile"] = _percentile(induced_values, evidence.get("induced"))
    # Co-iModulon edges are intentionally separate from regulatory edges.
    by_imodulon: dict[str, list[str]] = {}
    for node in candidates:
        for imodulon in graph.nodes[node].get("imodulons", []):
            by_imodulon.setdefault(imodulon, []).append(node)
    co_edge_count = 0
    for imodulon, nodes in by_imodulon.items():
        for source, target in combinations(sorted(set(nodes)), 2):
            if graph.has_edge(source, target) and graph.edges[source, target].get("edge_type") in REGULATORY_EDGE_TYPES:
                continue
            if graph.has_edge(source, target):
                graph.edges[source, target].setdefault("imodulons", []).append(imodulon)
                graph.edges[source, target]["n_shared_imodulons"] = len(set(graph.edges[source, target]["imodulons"]))
            else:
                graph.add_edge(source, target, edge_type="co-imodulon", interaction_type="co-membership", source="PRECISE", imodulons=[imodulon], n_shared_imodulons=1)
                co_edge_count += 1
    rows = []
    for node in candidates:
        data = graph.nodes[node]
        rows.append({"gene": node, **data, "imodulons": _format_list(data.get("imodulons", [])),
                     "imodulon_weights": _format_weight_map(data.get("imodulon_weights", {})),
                     "gene_expression_json": json.dumps(data.get("gene_expression", {}), sort_keys=True, default=str),
                     "imodulon_activity_json": json.dumps(data.get("imodulon_activity", {}), sort_keys=True, default=str),
                     "gene_expression_basal": json.dumps(data.get("gene_expression_basal", {}), sort_keys=True, default=str),
                     "gene_expression_induced": json.dumps(data.get("gene_expression_induced", {}), sort_keys=True, default=str),
                     "gene_expression_delta": json.dumps(data.get("gene_expression_delta", {}), sort_keys=True, default=str),
                     "gene_expression_available": json.dumps(data.get("gene_expression_available", {}), sort_keys=True, default=str),
                     "imodulon_activity_basal": json.dumps(data.get("imodulon_activity_basal", {}), sort_keys=True, default=str),
                     "imodulon_activity_induced": json.dumps(data.get("imodulon_activity_induced", {}), sort_keys=True, default=str),
                     "imodulon_activity_delta": json.dumps(data.get("imodulon_activity_delta", {}), sort_keys=True, default=str),
                     "imodulon_activity_available": json.dumps(data.get("imodulon_activity_available", {}), sort_keys=True, default=str),
                     "burden_proxy_system_categories": _format_list(data.get("burden_proxy_system_categories", [])),
                     "burden_proxy_imodulons": _format_list(data.get("burden_proxy_imodulons", []))})
    candidate_df = pd.DataFrame(rows)
    discovered_df = pd.DataFrame(discovered_rows) if discovered_rows else None
    return graph, candidate_df, discovered_df, co_edge_count


def save_outputs(graph: nx.DiGraph, candidate_df: pd.DataFrame, enrichment_df: pd.DataFrame, out_dir: str | Path, graph_path: str | Path, discovered_df: pd.DataFrame | None, summary_lines: list[str]) -> None:
    output = Path(out_dir)
    output.mkdir(parents=True, exist_ok=True)
    with Path(graph_path).open("wb") as handle:
        pickle.dump(graph, handle)
    node_df = pd.DataFrame([{"gene": node, **data} for node, data in graph.nodes(data=True)])
    for frame in (node_df, candidate_df):
        for column in frame.columns:
            frame[column] = frame[column].map(lambda value: json.dumps(value, sort_keys=True, default=str) if isinstance(value, (dict, list, tuple, set)) else value)
    node_df.to_csv(output / "imodulon_node_table.csv", index=False)
    candidate_df.to_csv(output / "imodulon_candidate_annotations.csv", index=False)
    enrichment_df.to_csv(output / "imodulon_regulon_enrichment.csv", index=False)
    if discovered_df is not None:
        discovered_df.to_csv(output / "imodulon_discovered_candidates.csv", index=False)
    (output / "imodulon_analysis_summary.txt").write_text("\n".join(summary_lines).rstrip() + "\n", encoding="utf-8")


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--graph", default=root / "network_analysis" / "output" / "regulatory_network.pkl")
    parser.add_argument("--precise", required=True, help="PRECISE/iModulon model JSON.GZ")
    parser.add_argument("--expression", default=None, help="Optional same-release companion expression matrix; embedded log_tpm/X is preferred when available")
    parser.add_argument("--mapping", default=None, help="Same-release gene/locus mapping")
    parser.add_argument("--out", default=root / "network_analysis" / "output")
    parser.add_argument("--out-graph", default=None)
    parser.add_argument("--expression-units", default=None)
    parser.add_argument("--expression-normalization", default=None)
    parser.add_argument("--fdr", type=float, default=0.01)
    parser.add_argument("--max-regs", type=int, default=1)
    parser.add_argument("--method", choices=["and", "or", "both"], default="both")
    parser.add_argument("--add-discovered", action="store_true")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv)
    for path, label in ((args.graph, "--graph"), (args.precise, "--precise")):
        if not Path(path).exists():
            raise FileNotFoundError(f"{label} file not found: {path}")
    if args.expression and not Path(args.expression).exists():
        raise FileNotFoundError(f"--expression file not found: {args.expression}")
    graph = pickle.loads(Path(args.graph).read_bytes())
    ica = load_precise_model(args.precise)
    if args.expression:
        expression, metadata = load_expression_matrix(
            args.expression,
            units=args.expression_units or "unknown",
            normalization=args.expression_normalization or "provided",
        )
    else:
        expression, metadata = load_embedded_expression(ica)
    if expression is None:
        raise ValueError(
            "No gene-level expression is available: provide --expression or use an IcaData model containing log_tpm/X"
        )
    validate_expression_alignment(expression, ica)
    mapping = load_gene_mapping(args.mapping)
    lookup = build_gene_lookup(ica, mapping)
    membership, members_by_imodulon = build_membership_tables(ica)
    enrichment, summary = summarize_imodulon_enrichment(ica, members_by_imodulon, args.fdr, args.max_regs, args.method)
    groups = build_condition_groups(ica)
    graph, candidates, discovered, co_edges = annotate_graph(
        graph, ica, lookup, membership, summary, groups, expression,
        args.expression_units or metadata["units"], args.expression_normalization or metadata["normalization"], args.add_discovered,
    )
    graph_path = args.out_graph or str(Path(args.out) / "regulatory_network_imodulon.pkl")
    save_outputs(graph, candidates, enrichment, args.out, graph_path, discovered,
                 [f"candidate nodes: {len(candidates)}", f"co-iModulon edges: {co_edges}",
                  f"Expression source: {metadata['source']} ({metadata['units']}; {metadata['normalization']})",
                  "Expression and activity values are evidence proxies; burden proxies are heuristic and lower is preferable."])


if __name__ == "__main__":
    main()
