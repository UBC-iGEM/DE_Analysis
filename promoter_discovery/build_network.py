"""Build a provenance-preserving regulatory network from configured DE data.

The input contract is the repository's ``results/differential_expression/contrasts/<comparison>/de_results.csv``
schema.  Every source row is retained as evidence; only candidate seeding is
direction-selectable.  RegulonDB edges and co-iModulon edges have distinct
``edge_type`` values so downstream scoring cannot treat co-membership as
regulation.
"""

from __future__ import annotations

import argparse
import json
import pickle
import re
from pathlib import Path
from typing import Any, Iterable

import networkx as nx
import pandas as pd

try:  # Support both ``python -m promoter_discovery.build_network`` and direct scripts.
    from .layout import contrast_directory
    from .dataset_registry import (
        REQUIRED_DE_COLUMNS,
        configured_datasets_by_class,
        dataset_caveats,
        load_dataset_config,
    )
except ImportError:  # pragma: no cover - exercised by direct CLI use
    from layout import contrast_directory
    from dataset_registry import (  # type: ignore
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
    # A gene-level DiGraph cannot represent individual promoter contexts.
    # Preserve every evidence record; opposing signs become an explicit dual
    # relationship rather than whichever record happened to sort first.
    collapsed = []
    for _, group in frame.groupby(["regulator", "target"], sort=True):
        row = group.sort_values("conf_rank", ascending=False).iloc[0].to_dict()
        signs = set(group["edge_type"])
        row["edge_type"] = next(iter(signs)) if len(signs) == 1 else "dual"
        row["interaction_evidence"] = group.to_dict("records")
        collapsed.append(row)
    frame = pd.DataFrame(collapsed, columns=list(frame.columns) + ["interaction_evidence"])
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


def load_de_results(
    config_path: str | Path,
    data_dir: str | Path | None = None,
    mapping_path: str | Path | None = None,
    candidate_direction: str = "upregulated",
    top_n: int | None = None,
    min_fc: float | None = None,
    provenance_path: str | Path | None = None,
) -> pd.DataFrame:
    """Load all configured standardized DE files and derive candidate seeds.

    ``candidate_direction`` controls only seeding.  Non-qualifying observations
    remain in the returned frame and are attached to candidate nodes as evidence.
    """

    if candidate_direction not in {"upregulated", "either-direction"}:
        raise ValueError("candidate_direction must be 'upregulated' or 'either-direction'")
    config = load_dataset_config(config_path)
    root = Path(data_dir) if data_dir else Path(config["_path"]).parents[1]
    provenance = None
    if config.get("benchmark"):
        try:
            from .setup_data import sha256_file
        except ImportError:  # pragma: no cover
            from setup_data import sha256_file
        provenance_file = Path(provenance_path) if provenance_path else root / "results/differential_expression/model/benchmark_provenance.json"
        provenance = json.loads(provenance_file.read_text())
        if provenance.get("config_sha256") != sha256_file(config_path):
            raise ValueError("Benchmark config differs from the joint fit; rerun promoter_discovery.run_benchmark")
    mapping = load_identity_mapping(mapping_path)
    fc_threshold = float(min_fc if min_fc is not None else config["thresholds"]["log2_fold_change"])
    padj_threshold = float(config["thresholds"]["padj"])
    frames: list[pd.DataFrame] = []
    for dataset in config["datasets"]:
        name = dataset["name"]
        path = contrast_directory(root, dataset) / "de_results.csv"
        if not path.exists():
            raise FileNotFoundError(
                f"Standardized DE file missing for {name}: {path}. "
                "Run the upstream DE analysis or setup documented in README.md."
            )
        if provenance and provenance.get("de_exports", {}).get(name) != sha256_file(path):
            raise ValueError(f"Benchmark DE export is not verified by joint-fit provenance: {name}; rerun the benchmark")
        frame = pd.read_csv(path)
        if dataset.get("de_method") == "pydeseq2" and ("de_method" not in frame or not frame["de_method"].astype(str).str.startswith("PyDESeq2 ").all()):
            raise ValueError(f"{name} is not a count-model benchmark export; run promoter_discovery.run_benchmark")
        frame.columns = [str(column).strip() for column in frame.columns]
        missing = REQUIRED_DE_COLUMNS - set(frame.columns)
        if missing:
            raise ValueError(f"{path} missing required DE columns: {sorted(missing)}")
        frame = frame.copy()
        frame["source_dataset"] = name
        frame["source_study"] = dataset.get("source_study", "unknown")
        frame["dataset_role"] = dataset.get("role", "supplemental")
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
        # Derive direction from the current numeric estimate, not stale labels.
        frame["regulation"] = frame["log2FoldChange"].map(
            lambda fc: "upregulated" if pd.notna(fc) and fc > 0 else "downregulated" if pd.notna(fc) and fc < 0 else "not_regulated"
        )
        finite_test = frame["padj"].between(0, 1) & frame["log2FoldChange"].map(lambda x: pd.notna(x) and abs(x) != float("inf"))
        if "eligible_for_network" in frame:
            eligibility = frame["eligible_for_network"].astype(str).str.lower().isin(["true", "1"])
            finite_test &= eligibility
        frame["eligible_for_network"] = finite_test
        frame["is_significant"] = (
            finite_test & frame["padj"].lt(padj_threshold)
            & frame["log2FoldChange"].abs().gt(fc_threshold)
        )
        frame.loc[~frame["is_significant"], "regulation"] = "not_regulated"
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
    # Challenge drugs assess responses of discovery candidates; their hits alone
    # do not seed the promoter shortlist.
    seed_evidence = observations["is_qualifying_direction"] & observations["dataset_role"].ne("cross_reactivity_challenge")
    observations["candidate_seed"] = seed_evidence.groupby(observations["canonical_gene"]).transform("any")
    if top_n and top_n > 0:
        candidates = (
            observations[seed_evidence]
            .groupby(["antibiotic_class", "canonical_gene"], as_index=False)["signal_strength"]
            .max()
        )
        keep = set(
            candidates.sort_values(
                ["antibiotic_class", "signal_strength", "canonical_gene"],
                ascending=[True, False, True],
            )
            .groupby("antibiotic_class", sort=False)
            .head(top_n)["canonical_gene"]
        )
        observations["candidate_seed"] &= observations["canonical_gene"].isin(keep)
    observations.attrs["candidate_direction_policy"] = candidate_direction
    if provenance:
        observations.attrs["benchmark_provenance"] = provenance
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
        def opposed(frame):
            return {"upregulated", "downregulated"}.issubset(set(frame["regulation"]))
        within_class_conflicts = {key for key, subset in significant.groupby("antibiotic_class") if opposed(subset)}
        direction_conflict = bool(within_class_conflicts)
        cross_class_direction_difference = opposed(significant) and not direction_conflict
        sig_classes = sorted(qualifying["antibiotic_class"].unique().tolist())
        evidence_by_dataset: list[dict[str, Any]] = []
        for dataset, dataset_group in group.groupby("source_dataset", sort=True):
            evidence_by_dataset.append(
                {
                    "source_dataset": dataset,
                    "source_study": str(dataset_group["source_study"].iloc[0]),
                    "dataset_role": str(dataset_group["dataset_role"].iloc[0]),
                    "antibiotic_class": str(dataset_group["antibiotic_class"].iloc[0]),
                    "observed": True,
                    "significant": bool(dataset_group["is_significant"].any()),
                    "qualifying": bool(dataset_group["is_qualifying_direction"].any()),
                    "regulations": sorted(dataset_group["regulation"].dropna().unique().tolist()),
                    "log2FoldChange": [float(x) for x in dataset_group["log2FoldChange"].dropna().tolist()],
                    "padj": [float(x) for x in dataset_group["padj"].dropna().tolist()],
                    "source_row_ids": dataset_group["source_row_id"].tolist(),
                    "eligible": bool(dataset_group["eligible_for_network"].any()),
                    "response_status": "significant_response" if dataset_group["is_significant"].any() else "no_detected_response" if dataset_group["eligible_for_network"].any() else "unavailable",
                    **{column: dataset_group[column].dropna().tolist() for column in ("lfcSE", "lfc_ci_low", "lfc_ci_high", "benchmark_basal", "benchmark_induced", "benchmark_expression_units") if column in dataset_group},
                }
            )
        tiers: dict[str, str] = {}
        support_fraction: dict[str, float] = {}
        observed_by_class: dict[str, list[str]] = {}
        significant_by_class: dict[str, list[str]] = {}
        for class_key, datasets in configured.items():
            class_group = group[group["antibiotic_class"] == class_key]
            class_qualifying = qualifying[qualifying["antibiotic_class"] == class_key]
            observed_by_class[class_key] = sorted(class_group["source_dataset"].unique().tolist())
            significant_by_class[class_key] = sorted(class_qualifying["source_dataset"].unique().tolist())
            support_fraction[class_key] = round(len(significant_by_class[class_key]) / len(datasets), 3) if datasets else 0.0
            if not significant_by_class[class_key]:
                continue
            if class_key in within_class_conflicts:
                tiers[class_key] = "within_class_direction_difference"
            elif len(significant_by_class[class_key]) > 1:
                tiers[class_key] = "multiple_comparisons_supported"
            else:
                tiers[class_key] = "single_comparison_supported"
        if direction_conflict:
            tier = "within_class_direction_difference"
        elif "multiple_comparisons_supported" in tiers.values():
            tier = "multiple_comparisons_supported"
        else:
            tier = "single_comparison_supported"
        caveats = sorted({c for dataset in group["source_dataset"].unique() for c in dataset_caveats(dataset)})
        caveats += [d["caveat"] for d in config["datasets"] if d["name"] in set(group.source_dataset) and d.get("caveat")]
        flags: list[str] = []
        if direction_conflict:
            flags.append("within_class_direction_difference")
        if cross_class_direction_difference:
            flags.append("cross_class_direction_difference")
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
                "cross_class_direction_difference": cross_class_direction_difference,
                "independent_validation": False,
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
                "primary_max_abs_log2_fold_change": float(qualifying[qualifying["dataset_role"].ne("cross_reactivity_challenge")]["log2FoldChange"].abs().max()),
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
    """Keep the full reference network and overlay candidate evidence.

    Statistical observations are stored on the graph independently of candidate
    selection, so top-n display filters cannot change enrichment denominators.
    """

    observations = observations.copy()
    observations.attrs["candidate_direction_policy"] = candidate_direction
    candidates = aggregate_candidate_evidence(observations, config)
    graph = nx.DiGraph()
    graph.graph["source_observations"] = observations.to_dict("records")
    if "benchmark_provenance" in observations.attrs:
        graph.graph["benchmark_provenance"] = observations.attrs["benchmark_provenance"]
    graph.graph["dataset_config"] = config
    graph.graph["network_scope"] = "full_reference"
    for row in candidates.to_dict("records"):
        attrs = dict(row)
        attrs["candidate_direction_policy"] = candidate_direction
        graph.add_node(row["canonical_gene"], **attrs)
    for edge in regulondb.to_dict("records"):
        target = str(edge["target"]).lower()
        if target not in graph:
            graph.add_node(target, node_type="reference_gene", group="reference")
        regulator = str(edge["regulator"]).lower()
        if regulator not in graph:
            graph.add_node(regulator, node_type="regulator", group="tf")
        # A responsive gene can also regulate other candidates (or itself).
        # Preserve its candidate identity and class evidence while marking
        # the additional regulatory role for downstream scoring.
        graph.nodes[regulator].update(
            is_regulator=True,
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
            interaction_evidence=edge.get("interaction_evidence", []),
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
    parser.add_argument("--config", default=root / "config" / "benchmark.json")
    parser.add_argument("--data-dir", default=root)
    parser.add_argument("--mapping", default=root / "data/references/regulondb/network_gene_mapping.tsv", help="Same-release gene/locus mapping TSV or CSV")
    parser.add_argument("--regulator", default=root / "data/references/regulondb/network_regulator_gene.tsv")
    parser.add_argument("--sigma", default=root / "data/references/regulondb/network_sigma_gene.tsv")
    parser.add_argument("--min-confidence", choices=sorted(CONFIDENCE_RANK), default="W")
    parser.add_argument("--reference-lock", type=Path, default=root / "config/regulondb.lock.json")
    parser.add_argument("--benchmark-provenance", type=Path, default=None)
    parser.add_argument("--candidate-direction", choices=["upregulated", "either-direction"], default="upregulated")
    parser.add_argument("--top-n", type=int, default=None)
    parser.add_argument("--min-fc", type=float, default=None)
    parser.add_argument("--out", default=root / "results" / "regulatory_network")
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
        args.benchmark_provenance,
    )
    if not Path(args.regulator).exists():
        raise FileNotFoundError(
            f"RegulonDB regulator file missing: {args.regulator}. "
            "Run promoter_discovery/setup_data.py with a >=14.5.0 release asset."
        )
    regulondb = load_regulondb(args.regulator, args.sigma, args.min_confidence)
    graph = build_graph(observations, regulondb, config, args.candidate_direction)
    try:
        from .setup_data import sha256_file
    except ImportError:  # pragma: no cover
        from setup_data import sha256_file
    lock = json.loads(args.reference_lock.read_text())
    releases = set()
    verified = []
    for path in (args.regulator, args.sigma, args.mapping):
        path = Path(path)
        record = next((record for record in lock["assets"].values()
                       if (_repo_root() / record["path"]).resolve() == path.resolve()), None)
        if not record or record["sha256"] != sha256_file(path):
            raise ValueError(f"Reference is not verified by the same-release lock: {path}")
        releases.add(record.get("regulondb_release") or record.get("source_release") or record.get("dataset_release"))
        verified.append({"path": str(path), "sha256": record["sha256"]})
    if len(releases) != 1 or None in releases:
        raise ValueError("Network and gene mapping must share one verified RegulonDB release")
    graph.graph["regulondb_release"] = releases.pop()
    graph.graph["reference_provenance"] = verified
    graph.graph["enrichment_background"] = "All measured genes with valid DE tests and verified same-release gene mapping; not candidate subset"
    pickle_path, node_path = save_graph(graph, args.out, observations)
    print(f"[network] {graph.number_of_nodes():,} nodes, {graph.number_of_edges():,} edges")
    print(f"[save] {pickle_path}")
    print(f"[save] {node_path}")


if __name__ == "__main__":
    main()
