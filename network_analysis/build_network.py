"""Build a provenance-preserving regulatory network from configured DE data.

The input contract is the repository's ``data/<dataset>/standardized/de_results.csv``
schema.  Every source row is retained as evidence; only candidate seeding is
direction-selectable.  RegulonDB edges and co-iModulon edges have distinct
``edge_type`` values so downstream scoring cannot treat co-membership as
regulation.
"""

from __future__ import annotations

import argparse
import json
import os
import pickle
import re
from collections import defaultdict
from pathlib import Path
from typing import Any, Iterable

import networkx as nx
import pandas as pd

try:  # Support both ``python -m network_analysis.build_network`` and direct scripts.
    from .dataset_registry import (
        CLASS_REGISTRY,
        DATASET_CAVEATS,
        REQUIRED_DE_COLUMNS,
        configured_datasets_by_class,
        dataset_caveats,
        load_dataset_config,
    )
except ImportError:  # pragma: no cover - exercised by direct CLI use
    from dataset_registry import (  # type: ignore
        CLASS_REGISTRY,
        DATASET_CAVEATS,
        REQUIRED_DE_COLUMNS,
        configured_datasets_by_class,
        dataset_caveats,
        load_dataset_config,
    )


EFFECT_MAP = {"+": "activates", "-": "represses", "-+": "dual", "+-": "dual"}
CONFIDENCE_RANK = {"C": 3, "S": 2, "W": 1, "?": 0, "": 0}
SIGMA_NAME_MAP = {
    "sigma19": "feci",
    "sigma24": "rpoe",
    "sigma28": "flia",
    "sigma32": "rpoh",
    "sigma38": "rpos",
    "sigma54": "rpon",
    "sigma70": "rpod",
}
REGULATORY_EDGE_TYPES = {"activates", "represses", "dual"}
_LOCUS_RE = re.compile(r"^b\d{4}$", re.IGNORECASE)


def _is_header_or_comment(parts: list[str]) -> bool:
    if not parts:
        return True
    first = parts[0].strip()
    return not first or first.startswith("#") or (
        len(first) >= 2 and first[0].isdigit() and ")" in first[:4]
    )


def _parse_regulator_gene_file(path: str | Path) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    with Path(path).open(encoding="utf-8") as handle:
        for line in handle:
            parts = line.rstrip("\n\r").split("\t")
            if _is_header_or_comment(parts) or len(parts) < 7:
                continue
            regulator, regulator_gene, target, effect, confidence = (
                parts[1].strip(), parts[2].strip(), parts[4].strip(), parts[5].strip(), parts[6].strip()
            )
            if regulator and target:
                rows.append(
                    {
                        "regulator": regulator.lower(),
                        "target": target.lower(),
                        "effect": effect,
                        "confidence": confidence or "?",
                        "regulator_type": "TF" if regulator_gene else "effector",
                        "source": "regulator-gene",
                    }
                )
    return rows


def _parse_sigma_gene_file(path: str | Path) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    with Path(path).open(encoding="utf-8") as handle:
        for line in handle:
            parts = line.rstrip("\n\r").split("\t")
            if _is_header_or_comment(parts) or len(parts) < 5:
                continue
            sigma, target, effect, confidence = parts[0].strip(), parts[1].strip(), parts[2].strip(), parts[4].strip()
            if sigma and target:
                rows.append(
                    {
                        "regulator": SIGMA_NAME_MAP.get(sigma.lower(), sigma.lower()),
                        "target": target.lower(),
                        "effect": effect,
                        "confidence": confidence or "?",
                        "regulator_type": "sigma",
                        "source": "sigma-gene",
                    }
                )
    return rows


def load_regulondb(
    regulator_gene_path: str | Path,
    sigma_gene_path: str | Path | None = None,
    min_confidence: str = "W",
) -> pd.DataFrame:
    """Load same-release RegulonDB flat files and collapse duplicate evidence."""

    rows = _parse_regulator_gene_file(regulator_gene_path)
    if sigma_gene_path and Path(sigma_gene_path).exists():
        rows.extend(_parse_sigma_gene_file(sigma_gene_path))
    elif sigma_gene_path:
        print(f"[warn] sigma-gene file not found: {sigma_gene_path}")
    if not rows:
        raise ValueError("No regulatory interactions parsed; check RegulonDB paths and release format")
    frame = pd.DataFrame(rows)
    frame["edge_type"] = frame["effect"].map(EFFECT_MAP).fillna("unknown")
    frame["conf_rank"] = frame["confidence"].map(CONFIDENCE_RANK).fillna(0)
    minimum = CONFIDENCE_RANK.get(min_confidence, 0)
    frame = frame[frame["conf_rank"] >= minimum]
    frame = frame[frame["edge_type"].isin(REGULATORY_EDGE_TYPES)].copy()
    frame = frame.sort_values("conf_rank", ascending=False).drop_duplicates(["regulator", "target"])
    return frame.reset_index(drop=True)


def _clean(value: Any) -> str:
    if value is None or pd.isna(value):
        return ""
    return str(value).strip().lower()


def _first_value(row: pd.Series, columns: Iterable[str]) -> str:
    for column in columns:
        if column in row and _clean(row[column]):
            return _clean(row[column])
    return ""


def load_identity_mapping(path: str | Path | None) -> dict[str, tuple[str, str]]:
    """Load a same-release gene/locus mapping as ``token -> (gene, locus)``."""

    if not path:
        return {}
    mapping_path = Path(path)
    if not mapping_path.exists():
        raise FileNotFoundError(
            f"Identity mapping not found: {mapping_path}. "
            "Ceftazidime b-number identifiers and probe aliases require a same-release mapping."
        )
    frame = pd.read_csv(mapping_path, sep=None, engine="python")
    frame.columns = [str(column).strip().lower() for column in frame.columns]
    result: dict[str, tuple[str, str]] = {}
    gene_columns = ("canonical_gene", "gene", "gene_name", "symbol", "gene_symbol", "name")
    locus_columns = ("canonical_locus_tag", "locus_tag", "b_number", "gene_id", "locus")
    alias_columns = ("alias", "aliases", "probe_id", "gene_id", "locus_tag", "b_number", "gene")
    for _, row in frame.iterrows():
        gene = _first_value(row, gene_columns)
        locus = _first_value(row, locus_columns)
        aliases: set[str] = set()
        for column in alias_columns:
            if column not in row or not _clean(row[column]):
                continue
            aliases.update(_clean(part) for part in re.split(r"[|;,]", str(row[column])) if _clean(part))
        aliases.update({gene, locus})
        aliases.discard("")
        if not gene and locus:
            gene = locus
        for alias in aliases:
            result.setdefault(alias, (gene, locus))
    return result


def _identity(row: pd.Series, mapping: dict[str, tuple[str, str]]) -> tuple[str, str, str, str]:
    source_gene = _clean(row.get("gene"))
    source_gene_id = _first_value(row, ("gene_id", "gene", "probe_id"))
    source_locus = _first_value(row, ("locus_tag", "b_number", "locus"))
    tokens = [source_locus, source_gene_id, source_gene]
    for token in tokens:
        if token and token in mapping:
            gene, locus = mapping[token]
            return gene or token, locus or (token if _LOCUS_RE.match(token) else ""), "same_release_mapping", source_locus
    fallback = source_gene or source_gene_id or source_locus
    return fallback, source_locus or (source_gene_id if _LOCUS_RE.match(source_gene_id) else ""), "source_gene", source_locus


def _normalize_regulation(row: pd.Series, fc_threshold: float) -> str:
    raw = _clean(row.get("regulation"))
    if raw in {"upregulated", "up", "induced"}:
        return "upregulated"
    if raw in {"downregulated", "down", "repressed"}:
        return "downregulated"
    if raw in {"not_regulated", "not regulated", "ns", "none"}:
        fc = pd.to_numeric(row.get("log2FoldChange"), errors="coerce")
        padj = pd.to_numeric(row.get("padj"), errors="coerce")
        if pd.notna(fc) and pd.notna(padj) and padj <= 0.05 and abs(fc) >= fc_threshold:
            return "upregulated" if fc > 0 else "downregulated"
        return "not_regulated"
    fc = pd.to_numeric(row.get("log2FoldChange"), errors="coerce")
    return "upregulated" if pd.notna(fc) and fc > 0 else "downregulated" if pd.notna(fc) else "unknown"


def load_de_results(
    config_path: str | Path,
    data_dir: str | Path | None = None,
    mapping_path: str | Path | None = None,
    candidate_direction: str = "upregulated",
    top_n: int | None = None,
    min_fc: float | None = None,
) -> pd.DataFrame:
    """Load all configured standardized DE files and derive candidate seeds.

    ``candidate_direction`` controls only seeding.  Non-qualifying observations
    remain in the returned frame and are attached to candidate nodes as evidence.
    """

    if candidate_direction not in {"upregulated", "either-direction"}:
        raise ValueError("candidate_direction must be 'upregulated' or 'either-direction'")
    config = load_dataset_config(config_path)
    root = Path(data_dir) if data_dir else Path(config["_path"]).parents[1]
    mapping = load_identity_mapping(mapping_path)
    fc_threshold = float(min_fc if min_fc is not None else config["thresholds"]["log2_fold_change"])
    padj_threshold = float(config["thresholds"]["padj"])
    frames: list[pd.DataFrame] = []
    for dataset in config["datasets"]:
        name = dataset["name"]
        path = root / "data" / name / "standardized" / "de_results.csv"
        if not path.exists():
            raise FileNotFoundError(
                f"Standardized DE file missing for {name}: {path}. "
                "Run the upstream DE analysis or setup documented in network_analysis/README.md."
            )
        frame = pd.read_csv(path)
        frame.columns = [str(column).strip() for column in frame.columns]
        missing = REQUIRED_DE_COLUMNS - set(frame.columns)
        if missing:
            raise ValueError(f"{path} missing required DE columns: {sorted(missing)}")
        frame = frame.copy()
        frame["source_dataset"] = name
        frame["antibiotic_class"] = dataset["antibiotic_class"]
        frame["source_gene"] = frame["gene"].astype(str).str.strip()
        frame["source_gene_id"] = frame.apply(lambda row: _first_value(row, ("gene_id", "gene", "probe_id")), axis=1)
        identities = frame.apply(lambda row: _identity(row, mapping), axis=1, result_type="expand")
        identities.columns = ["canonical_gene", "canonical_locus_tag", "identity_source", "source_locus_tag"]
        frame = pd.concat([frame, identities], axis=1)
        frame["canonical_gene"] = frame["canonical_gene"].replace("", pd.NA).fillna(frame["source_gene"].map(_clean))
        frame["canonical_gene"] = frame["canonical_gene"].map(_clean)
        frame["canonical_locus_tag"] = frame["canonical_locus_tag"].map(_clean)
        frame["source_locus_tag"] = frame["source_locus_tag"].map(_clean)
        frame["log2FoldChange"] = pd.to_numeric(frame["log2FoldChange"], errors="coerce")
        frame["padj"] = pd.to_numeric(frame["padj"], errors="coerce")
        frame["regulation"] = frame.apply(lambda row: _normalize_regulation(row, fc_threshold), axis=1)
        frame["is_significant"] = frame["padj"].le(padj_threshold) & frame["log2FoldChange"].abs().ge(fc_threshold)
        frame["is_qualifying_direction"] = frame["is_significant"] & (
            frame["regulation"].eq("upregulated")
            if candidate_direction == "upregulated"
            else frame["regulation"].isin(["upregulated", "downregulated"])
        )
        if "signal_strength" in frame.columns:
            frame["signal_strength"] = pd.to_numeric(frame["signal_strength"], errors="coerce").fillna(
                frame["log2FoldChange"].abs()
            )
        else:
            frame["signal_strength"] = frame["log2FoldChange"].abs()
        if "padj_source" in frame.columns:
            frame["padj_source"] = frame["padj_source"].fillna("unknown")
        else:
            frame["padj_source"] = "unknown"
        frame["source_row_id"] = [f"{name}:{index}" for index in frame.index]
        frames.append(frame)
    observations = pd.concat(frames, ignore_index=True)
    observations["candidate_seed"] = observations.groupby("canonical_gene")["is_qualifying_direction"].transform("any")
    if top_n and top_n > 0:
        candidates = observations[observations["candidate_seed"]].drop_duplicates("canonical_gene")
        keep = set(
            candidates.sort_values("signal_strength", ascending=False)
            .groupby("antibiotic_class", sort=False)
            .head(top_n)["canonical_gene"]
        )
        observations["candidate_seed"] &= observations["canonical_gene"].isin(keep)
    observations.attrs["candidate_direction_policy"] = candidate_direction
    return observations


def _json(value: Any) -> str:
    return json.dumps(value, sort_keys=True, default=str)


def aggregate_candidate_evidence(observations: pd.DataFrame, config: dict[str, Any]) -> pd.DataFrame:
    """Collapse observations by canonical gene without discarding provenance."""

    if observations.empty:
        return pd.DataFrame()
    configured = configured_datasets_by_class(config)
    rows: list[dict[str, Any]] = []
    for gene, group in observations[observations["candidate_seed"]].groupby("canonical_gene", sort=True):
        significant = group[group["is_significant"]]
        qualifying = group[group["is_qualifying_direction"]]
        up_count = int((significant["regulation"] == "upregulated").sum())
        down_count = int((significant["regulation"] == "downregulated").sum())
        direction_conflict = up_count > 0 and down_count > 0
        sig_classes = sorted(qualifying["antibiotic_class"].unique().tolist())
        evidence_by_dataset: list[dict[str, Any]] = []
        for dataset, dataset_group in group.groupby("source_dataset", sort=True):
            evidence_by_dataset.append(
                {
                    "source_dataset": dataset,
                    "antibiotic_class": str(dataset_group["antibiotic_class"].iloc[0]),
                    "observed": True,
                    "significant": bool(dataset_group["is_significant"].any()),
                    "qualifying": bool(dataset_group["is_qualifying_direction"].any()),
                    "regulations": sorted(dataset_group["regulation"].dropna().unique().tolist()),
                    "log2FoldChange": [float(x) for x in dataset_group["log2FoldChange"].dropna().tolist()],
                    "padj": [float(x) for x in dataset_group["padj"].dropna().tolist()],
                    "source_row_ids": dataset_group["source_row_id"].tolist(),
                }
            )
        tiers: dict[str, str] = {}
        support_fraction: dict[str, float] = {}
        observed_by_class: dict[str, list[str]] = {}
        significant_by_class: dict[str, list[str]] = {}
        for class_key, datasets in configured.items():
            class_group = group[group["antibiotic_class"] == class_key]
            class_qualifying = qualifying[qualifying["antibiotic_class"] == class_key]
            class_significant = significant[significant["antibiotic_class"] == class_key]
            observed_by_class[class_key] = sorted(class_group["source_dataset"].unique().tolist())
            significant_by_class[class_key] = sorted(class_qualifying["source_dataset"].unique().tolist())
            support_fraction[class_key] = round(len(significant_by_class[class_key]) / len(datasets), 3) if datasets else 0.0
            if not significant_by_class[class_key]:
                continue
            if direction_conflict:
                tiers[class_key] = "conflicted"
            elif set(significant_by_class[class_key]) == {"tobramycin"}:
                tiers[class_key] = "limited"
            elif len(significant_by_class[class_key]) == len(datasets):
                tiers[class_key] = "corroborated"
            else:
                tiers[class_key] = "supported"
        if direction_conflict:
            tier = "conflicted"
        elif any(value == "corroborated" for value in tiers.values()):
            tier = "corroborated"
        elif qualifying["source_dataset"].nunique() == 1 and qualifying["source_dataset"].iloc[0] == "tobramycin":
            tier = "limited"
        else:
            tier = "supported"
        caveats = sorted({c for dataset in group["source_dataset"].unique() for c in dataset_caveats(dataset)})
        flags: list[str] = []
        if direction_conflict:
            flags.append("direction_conflict")
        if set(qualifying["source_dataset"]) == {"tobramycin"}:
            flags.append("tobramycin_only")
        if (group["identity_source"] != "same_release_mapping").any():
            flags.append("identity_mapping_not_verified")
        rows.append(
            {
                "canonical_gene": gene,
                "canonical_locus_tag": next((x for x in group["canonical_locus_tag"] if x), ""),
                "group": "cross" if len(sig_classes) > 1 else (sig_classes[0] if sig_classes else "unknown"),
                "node_type": "candidate",
                "n_datasets_observed": int(group["source_dataset"].nunique()),
                "n_datasets_significant": int(qualifying["source_dataset"].nunique()),
                "n_significant_source_rows": int(len(significant)),
                "n_source_rows": int(len(group)),
                "upregulated_dataset_count": int(qualifying[qualifying["regulation"] == "upregulated"]["source_dataset"].nunique()),
                "downregulated_dataset_count": int(qualifying[qualifying["regulation"] == "downregulated"]["source_dataset"].nunique()),
                "datasets_observed_by_class": observed_by_class,
                "datasets_significant_by_class": significant_by_class,
                "configured_datasets_by_class": configured,
                "support_fraction_by_class": support_fraction,
                "evidence_tier_by_class": tiers,
                "evidence_tier": tier,
                "significant_classes": sig_classes,
                "direction_consistent": not direction_conflict,
                "candidate_direction_policy": observations.attrs.get("candidate_direction_policy", "upregulated"),
                "dataset_evidence_json": _json(evidence_by_dataset),
                "source_observations_json": _json(
                    group[
                        [
                            "source_dataset", "antibiotic_class", "source_gene", "source_gene_id",
                            "source_locus_tag", "canonical_gene", "canonical_locus_tag", "identity_source",
                            "log2FoldChange", "padj", "signal_strength", "regulation", "padj_source",
                        ]
                    ].to_dict("records")
                ),
                "source_row_ids": group["source_row_id"].tolist(),
                "direction_conflict": direction_conflict,
                "collapse_method": "canonical_gene; preserve all source rows",
                "signal_strength": float(group["signal_strength"].max()),
                "max_abs_log2_fold_change": float(group["log2FoldChange"].abs().max()),
                "caveats": caveats,
                "evidence_quality_flags": flags,
                "tobramycin_only": "tobramycin_only" in flags,
                "regulation_by_dataset": {
                    dataset: sorted(dataset_group["regulation"].unique().tolist())
                    for dataset, dataset_group in group.groupby("source_dataset", sort=True)
                },
            }
        )
    result = pd.DataFrame(rows)
    if not result.empty:
        result["candidate_direction_policy"] = observations.attrs.get(
            "candidate_direction_policy", "upregulated"
        )
    return result


def build_graph(
    observations: pd.DataFrame,
    regulondb: pd.DataFrame,
    config: dict[str, Any],
    candidate_direction: str = "upregulated",
) -> nx.DiGraph:
    """Build a graph with typed regulatory edges and candidate evidence."""

    observations = observations.copy()
    observations.attrs["candidate_direction_policy"] = candidate_direction
    candidates = aggregate_candidate_evidence(observations, config)
    graph = nx.DiGraph()
    if candidates.empty:
        return graph
    for row in candidates.to_dict("records"):
        attrs = dict(row)
        attrs["candidate_direction_policy"] = candidate_direction
        graph.add_node(row["canonical_gene"], **attrs)
    candidate_names = set(candidates["canonical_gene"])
    for edge in regulondb.to_dict("records"):
        target = str(edge["target"]).lower()
        if target not in candidate_names:
            continue
        regulator = str(edge["regulator"]).lower()
        graph.add_node(
            regulator,
            node_type="regulator",
            group="tf",
            regulator_type=edge.get("regulator_type", "TF"),
            label=regulator,
        )
        graph.add_edge(
            regulator,
            target,
            edge_type=edge["edge_type"],
            interaction_type="regulatory",
            confidence=edge.get("confidence", "?"),
            source=edge.get("source", "RegulonDB"),
        )
    return graph


def save_graph(graph: nx.DiGraph, out_dir: str | Path, observations: pd.DataFrame | None = None) -> tuple[Path, Path]:
    output = Path(out_dir)
    output.mkdir(parents=True, exist_ok=True)
    pickle_path = output / "regulatory_network.pkl"
    node_path = output / "node_table.csv"
    with pickle_path.open("wb") as handle:
        pickle.dump(graph, handle)
    records: list[dict[str, Any]] = []
    for node, attrs in graph.nodes(data=True):
        record = {"node": node, **attrs}
        for key, value in list(record.items()):
            if isinstance(value, (dict, list, tuple, set)):
                record[key] = _json(value)
        records.append(record)
    pd.DataFrame(records).to_csv(node_path, index=False)
    if observations is not None:
        observations.to_csv(output / "source_observations.csv", index=False)
    return pickle_path, node_path


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    root = _repo_root()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default=root / "config" / "datasets.json")
    parser.add_argument("--data-dir", default=root)
    parser.add_argument("--mapping", default=None, help="Same-release gene/locus mapping TSV or CSV")
    parser.add_argument("--regulator", default=root / "data" / "network_regulator_gene.tsv")
    parser.add_argument("--sigma", default=root / "data" / "network_sigma_gene.tsv")
    parser.add_argument("--min-confidence", choices=sorted(CONFIDENCE_RANK), default="W")
    parser.add_argument("--candidate-direction", choices=["upregulated", "either-direction"], default="upregulated")
    parser.add_argument("--top-n", type=int, default=None)
    parser.add_argument("--min-fc", type=float, default=None)
    parser.add_argument("--out", default=root / "network_analysis" / "output")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv)
    config = load_dataset_config(args.config)
    observations = load_de_results(
        args.config,
        args.data_dir,
        args.mapping,
        args.candidate_direction,
        args.top_n,
        args.min_fc,
    )
    if not Path(args.regulator).exists():
        raise FileNotFoundError(
            f"RegulonDB regulator file missing: {args.regulator}. "
            "Run network_analysis/setup_data.py with a >=14.5.0 release asset."
        )
    regulondb = load_regulondb(args.regulator, args.sigma, args.min_confidence)
    graph = build_graph(observations, regulondb, config, args.candidate_direction)
    pickle_path, node_path = save_graph(graph, args.out, observations)
    print(f"[network] {graph.number_of_nodes():,} nodes, {graph.number_of_edges():,} edges")
    print(f"[save] {pickle_path}")
    print(f"[save] {node_path}")


if __name__ == "__main__":
    main()
