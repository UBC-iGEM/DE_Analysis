"""
i_modulon_analysis.py
─────────────────────
Annotates the RegulonDB candidate graph with PRECISE-1K iModulon metadata,
computes iModulon -> regulon enrichment, and optionally reports additional
genes that co-occur with the current candidates inside the same iModulons.

Usage
-----
    python i_modulon_analysis.py \
        --graph output/regulatory_network.pkl \
        --precise ../data/precise1k.json.gz \
        --mapping ../caz_kan_DE/gene_mapping.csv \
        --out output/

Outputs
-------
    output/regulatory_network_imodulon.pkl
    output/imodulon_candidate_annotations.csv
    output/imodulon_regulon_enrichment.csv
    output/imodulon_node_table.csv
    output/imodulon_analysis_summary.txt
    output/imodulon_discovered_candidates.csv   (only with --add-discovered)
"""

from __future__ import annotations

import argparse
import gzip
import json
import os
import pickle
import sys
from dataclasses import dataclass
from io import StringIO
from itertools import combinations

import networkx as nx
import pandas as pd

from pymodulon.core import IcaData
from pymodulon.enrichment import compute_trn_enrichment


MATRIX_COLUMNS = ("M", "A", "gene_table", "sample_table", "imodulon_table", "trn")

# Keyword groups for matching PRECISE-1K sample_table['condition'] strings to
# our two DE-derived reporter classes. Confirmed against the real
# precise1k.json.gz sample_table (1,035 samples): the only antibiotic
# conditions present at all live under project == "abx_media" (30 samples,
# 2 reps x 3 media backgrounds x 5 drugs), and there is NO kanamycin or
# other aminoglycoside condition anywhere in the compendium. "cef" there is
# a cephalosporin (beta-lactam class) but not confirmed to be literally
# ceftazidime — treat this as a beta-lactam-class proxy, not an exact match.
CONDITION_KEYWORDS: dict[str, list[str]] = {
    "caz": ["cef", "mero", "imipenem", "ampicillin", "amoxicillin", "penicillin"],
    "kan": ["kan", "tobramycin", "gentamicin", "amikacin", "streptomycin"],
}

# Each PRECISE-1K abx_media condition is named "{media}_{drug}", with a
# matched "{media}_ctrl" for the same media background — e.g. "camhb_cef"
# pairs with "camhb_ctrl". This lets us compute a real basal-vs-induced
# delta per media background rather than comparing against a global mean.
CONTROL_TOKEN = "ctrl"

# PRECISE-1K's imodulon_table.system_category values (confirmed from the
# real table: Metabolism, Stress Responses, Genetic Alterations,
# Single Gene, Translation, ALE Effects, Unknown). Metabolism/Translation
# iModulons are broad housekeeping/growth-linked programs — a reporter gene
# whose primary iModulon falls here risks confounding the fluorescence
# readout with a general growth-rate effect rather than a specific
# antibiotic stress response. "Stress Responses" is the desirable category
# (this is where LexA/DNA Damage and RpoE·RpoH/Temperature Shock live).
BURDEN_RISK_SYSTEM_CATEGORIES = {"Metabolism", "Translation"}
BURDEN_OK_SYSTEM_CATEGORIES = {"Stress Responses"}


@dataclass
class GeneLookup:
    by_name: dict[str, str]
    by_id: dict[str, str]


def _read_json_table(value):
    if isinstance(value, pd.DataFrame):
        return value
    if value is None:
        return None
    return pd.read_json(StringIO(value))


def load_precise_model(path: str) -> IcaData:
    """
    Load a PRECISE-1K IcaData object from its JSON.GZ export.

    Two compatibility fixes are applied here, both confirmed against the
    actual precise1k.json.gz file and the installed pymodulon==0.2.1 /
    pandas==3.x versions pinned in this project:

    1. Private-cutoff kwargs.  The raw JSON includes "_cutoff_optimized" and
       "_dagostino_cutoff" keys, which are NOT accepted by
       IcaData.__init__() (pymodulon's own load_json_model pops them first
       and assigns them as private attributes afterward). Passing them
       straight through raises:
           TypeError: IcaData.__init__() got an unexpected keyword
           argument '_cutoff_optimized'

    2. M.columns int-cast bug (pymodulon 0.2.1 x pandas 3.x).  Inside
       IcaData.__init__, pymodulon does:
           try:
               M.columns = M.columns.astype(int)
           except TypeError:
               pass
       ...to silently skip the cast when iModulons have curated string
       names (e.g. PRECISE-1K's "Sugar Diacid", "Translation") rather than
       bare integer indices. Under pandas >= 2.x — including the pandas 3.x
       pinned here — that same failed cast raises ValueError, not
       TypeError, so pymodulon's except clause doesn't catch it and the
       whole load crashes. This is an upstream pymodulon/pandas version
       incompatibility, not a problem with the data file.

       We work around it with a *scoped* monkeypatch of pd.Index.astype
       that also swallows ValueError, active only for the duration of the
       IcaData(**serial) call below, and restored immediately afterward
       (in a finally block, so it's restored even if construction fails).
       This must never be left patched globally.
    """
    with gzip.open(path, "rt", encoding="utf-8") as fh:
        serial = json.load(fh)

    # Fix 1: pop private-cutoff keys before constructing IcaData; restore
    # them as attributes afterward, matching pymodulon's own load_json_model.
    cutoff_optimized = serial.pop("_cutoff_optimized", False)
    dagostino_cutoff = serial.pop("_dagostino_cutoff", None)

    for key in MATRIX_COLUMNS:
        if key in serial and serial[key] is not None:
            serial[key] = _read_json_table(serial[key])

    # Fix 2: scoped astype patch, restored immediately after use.
    _orig_index_astype = pd.Index.astype

    def _patched_index_astype(self, dtype, **kw):
        try:
            return _orig_index_astype(self, dtype, **kw)
        except (TypeError, ValueError):
            return self

    pd.Index.astype = _patched_index_astype
    try:
        ica = IcaData(**serial)
    finally:
        pd.Index.astype = _orig_index_astype

    ica._cutoff_optimized = cutoff_optimized
    ica._dagostino_cutoff = dagostino_cutoff
    return ica


def load_gene_mapping(path: str) -> pd.DataFrame:
    df = pd.read_csv(path)
    df.columns = [c.strip().lower() for c in df.columns]
    return df


def _canonical(name: str | None) -> str:
    return str(name).strip().lower() if name is not None else ""


def _split_aliases(value) -> list[str]:
    if pd.isna(value):
        return []
    text = str(value).strip()
    if not text:
        return []
    parts = []
    for chunk in text.replace("|", ";").replace(",", ";").split(";"):
        alias = chunk.strip().lower()
        if alias:
            parts.append(alias)
    return parts


def build_gene_lookup(ica: IcaData, mapping: pd.DataFrame | None = None) -> GeneLookup:
    by_name: dict[str, str] = {}
    by_id: dict[str, str] = {}

    if mapping is not None and not mapping.empty:
        cols = set(mapping.columns)
        gene_name_cols = [c for c in ("gene", "gene_name", "symbol") if c in cols]
        gene_id_cols = [c for c in ("gene_id", "locus_tag", "b_number") if c in cols]
        if gene_name_cols and gene_id_cols:
            gene_name_col = gene_name_cols[0]
            gene_id_col = gene_id_cols[0]
            for _, row in mapping.iterrows():
                gene_name = _canonical(row[gene_name_col])
                gene_id = _canonical(row[gene_id_col])
                if gene_name and gene_id:
                    by_name.setdefault(gene_name, gene_id)
                    by_id.setdefault(gene_id, gene_id)

    gene_table = ica.gene_table.copy()
    gene_table.index = gene_table.index.astype(str).str.lower()
    if "gene_name" in gene_table.columns:
        for gene_id, row in gene_table.iterrows():
            gene_name = _canonical(row.get("gene_name"))
            if gene_name:
                by_name.setdefault(gene_name, gene_id)
                by_id.setdefault(gene_id, gene_id)
            for alias in _split_aliases(row.get("synonyms")):
                by_name.setdefault(alias, gene_id)

    trn = ica.trn.copy()
    if not trn.empty:
        if "gene_name" in trn.columns and "gene_id" in trn.columns:
            for _, row in trn.iterrows():
                gene_name = _canonical(row.get("gene_name"))
                gene_id = _canonical(row.get("gene_id"))
                if gene_name and gene_id:
                    by_name.setdefault(gene_name, gene_id)
                    by_id.setdefault(gene_id, gene_id)

    return GeneLookup(by_name=by_name, by_id=by_id)


def resolve_gene_id(name: str, lookup: GeneLookup, ica: IcaData) -> str | None:
    query = _canonical(name)
    if not query:
        return None

    gene_index = {idx.lower(): idx for idx in ica.gene_table.index.astype(str)}
    if query in gene_index:
        return gene_index[query]

    if query in lookup.by_name:
        return lookup.by_name[query]

    if query in lookup.by_id:
        return lookup.by_id[query]

    if query.startswith("b") and query[1:].isdigit():
        return gene_index.get(query)

    return None


def build_membership_tables(ica: IcaData) -> tuple[pd.DataFrame, dict[str, set[str]]]:
    thresholds = pd.Series(ica.thresholds, dtype=float).reindex(ica.M.columns).fillna(0.0)
    membership = ica.M.abs().gt(thresholds, axis=1)
    members_by_imodulon = {
        imod: set(membership.index[membership[imod]])
        for imod in membership.columns
    }
    return membership, members_by_imodulon


def build_condition_groups(ica: IcaData) -> dict[str, dict[str, list[str]]]:
    """
    Identify treated-vs-matched-control sample groups per DE class, using
    PRECISE-1K's real "{media}_{drug}" / "{media}_ctrl" naming convention.

    Returns:
        {
          "caz": {"treated": [sample_ids...], "control": [sample_ids...],
                  "pairs": [(treated_id, control_id), ...]},
          "kan": {"treated": [], "control": [], "pairs": []},   # expected
                  empty — see CONDITION_KEYWORDS note above
        }

    Samples are paired by shared media prefix (the part of the condition
    string before the drug/ctrl token) so that dynamic range is computed as
    a same-background delta, not a comparison against an unrelated global
    control.
    """
    st = ica.sample_table
    if "condition" not in st.columns:
        return {cls: {"treated": [], "control": [], "pairs": []} for cls in CONDITION_KEYWORDS}

    groups: dict[str, dict[str, list]] = {
        cls: {"treated": [], "control": [], "pairs": []} for cls in CONDITION_KEYWORDS
    }

    conditions = st["condition"].astype(str)

    for cls, keywords in CONDITION_KEYWORDS.items():
        for sample_id, cond in conditions.items():
            cond_l = cond.lower()
            tokens = cond_l.split("_")
            drug_token = next((t for t in tokens if any(kw in t for kw in keywords)), None)
            if drug_token is None:
                continue
            media = cond_l.replace(f"_{drug_token}", "").strip("_") or cond_l.split("_")[0]
            control_cond = f"{media}_{CONTROL_TOKEN}"
            control_samples = conditions.index[conditions.str.lower() == control_cond].tolist()

            groups[cls]["treated"].append(sample_id)
            groups[cls]["control"].extend(control_samples)
            for ctrl_id in control_samples:
                groups[cls]["pairs"].append((sample_id, ctrl_id))

        groups[cls]["control"] = sorted(set(groups[cls]["control"]))
        n_t, n_c = len(groups[cls]["treated"]), len(groups[cls]["control"])
        print(f"[dynamic-range] {cls}: {n_t} treated samples, {n_c} matched-control samples "
              f"{'(no PRECISE-1K data for this class)' if n_t == 0 else ''}")

    return groups


def compute_dynamic_range(ica: IcaData,
                          condition_groups: dict[str, dict[str, list]]) -> dict[str, dict[str, dict]]:
    """
    For every iModulon, compute basal / induced / delta activity per DE
    class, using the matched-control pairs from build_condition_groups.

    Returns:
        {imodulon_name: {
            "caz": {"basal": float|None, "induced": float|None,
                    "delta": float|None, "n_treated": int, "n_control": int},
            "kan": {...},   # will be all-None if no PRECISE-1K coverage
        }}

    basal  = median A-matrix activity across matched control samples
    induced = median A-matrix activity across treated samples
    delta  = induced - basal  (the actual dynamic-range proxy: how far the
             iModulon moves from its own paired baseline, not just its raw
             magnitude in one condition)
    """
    A = ica.A
    result: dict[str, dict[str, dict]] = {}

    for imodulon in A.index:
        result[imodulon] = {}
        for cls, grp in condition_groups.items():
            treated = [s for s in grp["treated"] if s in A.columns]
            control = [s for s in grp["control"] if s in A.columns]
            if not treated or not control:
                result[imodulon][cls] = {
                    "basal": None, "induced": None, "delta": None,
                    "n_treated": len(treated), "n_control": len(control),
                }
                continue
            basal = float(A.loc[imodulon, control].median())
            induced = float(A.loc[imodulon, treated].median())
            result[imodulon][cls] = {
                "basal": round(basal, 4),
                "induced": round(induced, 4),
                "delta": round(induced - basal, 4),
                "n_treated": len(treated),
                "n_control": len(control),
            }

    return result


def host_burden_lookup(imodulon_table: pd.DataFrame) -> dict[str, dict[str, str]]:
    """
    Build a per-iModulon host-burden risk lookup from
    imodulon_table.system_category (PRECISE-1K's own curated taxonomy).

    Returns {imodulon_name: {"system_category": str, "functional_category": str,
                              "burden_flag": "risk" | "ok" | "unknown"}}

    "risk"    -> system_category in {Metabolism, Translation}: broad
                 housekeeping/growth-linked programs. A reporter here risks
                 confounding fluorescence with a general growth-rate effect.
    "ok"      -> system_category == "Stress Responses": a specific stress
                 program, not a general growth/metabolic one — this is
                 where LexA/DNA-Damage and RpoE·RpoH/Temperature-Shock live.
    "unknown" -> anything else (Genetic Alterations, Single Gene, ALE
                 Effects, Unknown) — not enough signal either way.
    """
    lookup: dict[str, dict[str, str]] = {}
    has_sys = "system_category" in imodulon_table.columns
    has_func = "functional_category" in imodulon_table.columns

    for imodulon, row in imodulon_table.iterrows():
        sys_cat = str(row.get("system_category", "")) if has_sys else ""
        func_cat = str(row.get("functional_category", "")) if has_func else ""
        if sys_cat in BURDEN_RISK_SYSTEM_CATEGORIES:
            flag = "risk"
        elif sys_cat in BURDEN_OK_SYSTEM_CATEGORIES:
            flag = "ok"
        else:
            flag = "unknown"
        lookup[imodulon] = {
            "system_category": sys_cat,
            "functional_category": func_cat,
            "burden_flag": flag,
        }
    return lookup


def summarize_imodulon_enrichment(ica: IcaData,
                                  members_by_imodulon: dict[str, set[str]],
                                  fdr: float,
                                  max_regs: int,
                                  method: str) -> tuple[pd.DataFrame, dict[str, dict]]:
    all_genes = set(ica.M.index.astype(str))
    rows = []
    summary: dict[str, dict] = {}

    for imodulon, gene_set in members_by_imodulon.items():
        enrich = compute_trn_enrichment(
            gene_set=gene_set,
            all_genes=all_genes,
            trn=ica.trn,
            max_regs=max_regs,
            fdr=fdr,
            method=method,
        )

        if enrich.empty:
            summary[imodulon] = {
                "top_regulator": "",
                "qvalue": None,
                "precision": None,
                "recall": None,
                "f1score": None,
                "n_significant": 0,
            }
            continue

        enrich = enrich.reset_index(names="regulon")
        enrich.insert(0, "imodulon", imodulon)
        rows.append(enrich)

        top = enrich.sort_values(
            ["qvalue", "f1score", "precision"], ascending=[True, False, False]
        ).iloc[0]
        summary[imodulon] = {
            "top_regulator": top["regulon"],
            "qvalue": float(top["qvalue"]),
            "precision": float(top["precision"]),
            "recall": float(top["recall"]),
            "f1score": float(top["f1score"]),
            "n_significant": int(len(enrich)),
        }

    enrich_df = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()
    return enrich_df, summary


def _format_weight_map(weight_map: dict[str, float]) -> str:
    if not weight_map:
        return ""
    parts = [f"{name}:{weight:.4f}" for name, weight in weight_map.items()]
    return "; ".join(parts)


def _format_modulon_list(modulons: list[str]) -> str:
    return ", ".join(modulons) if modulons else ""


def annotate_graph(G: nx.DiGraph,
                   ica: IcaData,
                   lookup: GeneLookup,
                   membership: pd.DataFrame,
                   imodulon_summary: dict[str, dict],
                   dynamic_range: dict[str, dict[str, dict]] | None = None,
                   burden_lookup: dict[str, dict[str, str]] | None = None,
                   add_discovered: bool = False) -> tuple[nx.DiGraph, pd.DataFrame, pd.DataFrame | None, int]:
    dynamic_range = dynamic_range or {}
    burden_lookup = burden_lookup or {}
    candidate_nodes = [node for node, data in G.nodes(data=True) if data.get("group") != "tf"]
    graph_candidate_gene_ids: dict[str, str] = {}
    discovered_rows = []

    for node in candidate_nodes:
        gene_id = resolve_gene_id(node, lookup, ica)
        if gene_id:
            graph_candidate_gene_ids[gene_id] = node
        node_data = G.nodes[node]
        node_data["precise_gene_id"] = gene_id or ""
        node_data["imodulon_mapping_status"] = "mapped" if gene_id else "not_found"
        node_data["imodulons"] = []
        node_data["imodulon_weights"] = {}
        node_data["imodulon_primary"] = ""
        node_data["imodulon_primary_weight"] = None
        node_data["imodulon_primary_abs_weight"] = None
        node_data["imodulon_primary_regulon"] = ""
        node_data["imodulon_primary_regulon_qvalue"] = None
        node_data["imodulon_primary_regulon_precision"] = None
        node_data["imodulon_primary_regulon_recall"] = None
        node_data["imodulon_primary_regulon_f1score"] = None
        node_data["imodulon_count"] = 0

        # Dynamic range / host-burden defaults (filled in below once the
        # primary iModulon is known; left as None/"" if no iModulon match
        # or no PRECISE-1K coverage for this gene's DE class).
        node_data["dr_class"] = ""
        node_data["dr_basal"] = None
        node_data["dr_induced"] = None
        node_data["dr_delta"] = None
        node_data["dr_n_treated"] = 0
        node_data["dr_n_control"] = 0
        node_data["dr_data_available"] = False
        node_data["burden_system_category"] = ""
        node_data["burden_functional_category"] = ""
        node_data["burden_flag"] = "unknown"

        if not gene_id or gene_id not in membership.index:
            continue

        row = membership.loc[gene_id]
        hits = [imod for imod, flag in row.items() if bool(flag)]
        weight_map = {imod: float(ica.M.at[gene_id, imod]) for imod in hits}
        hits = sorted(hits, key=lambda imod: abs(weight_map[imod]), reverse=True)

        node_data["imodulons"] = hits
        node_data["imodulon_weights"] = weight_map
        node_data["imodulon_count"] = len(hits)

        if hits:
            primary = hits[0]
            primary_weight = weight_map[primary]
            node_data["imodulon_primary"] = primary
            node_data["imodulon_primary_weight"] = float(primary_weight)
            node_data["imodulon_primary_abs_weight"] = float(abs(primary_weight))

            summary = imodulon_summary.get(primary, {})
            node_data["imodulon_primary_regulon"] = summary.get("top_regulator", "")
            node_data["imodulon_primary_regulon_qvalue"] = summary.get("qvalue")
            node_data["imodulon_primary_regulon_precision"] = summary.get("precision")
            node_data["imodulon_primary_regulon_recall"] = summary.get("recall")
            node_data["imodulon_primary_regulon_f1score"] = summary.get("f1score")

            # ── Dynamic range (basal / induced / delta) ─────────────────────
            # Matched to the node's DE class (group == 'caz' or 'kan'). Cross-
            # reactive candidates (group == 'cross') check both classes and
            # keep whichever has actual PRECISE-1K coverage.
            group = node_data.get("group", "")
            dr_candidates = ["caz", "kan"] if group == "cross" else [group]
            dr_row = None
            dr_cls_used = ""
            for cls in dr_candidates:
                candidate_dr = dynamic_range.get(primary, {}).get(cls)
                if candidate_dr and candidate_dr.get("delta") is not None:
                    dr_row, dr_cls_used = candidate_dr, cls
                    break
                if candidate_dr and dr_row is None:
                    dr_row, dr_cls_used = candidate_dr, cls  # keep as fallback (no-data) record

            if dr_row is not None:
                node_data["dr_class"] = dr_cls_used
                node_data["dr_basal"] = dr_row.get("basal")
                node_data["dr_induced"] = dr_row.get("induced")
                node_data["dr_delta"] = dr_row.get("delta")
                node_data["dr_n_treated"] = dr_row.get("n_treated", 0)
                node_data["dr_n_control"] = dr_row.get("n_control", 0)
                node_data["dr_data_available"] = dr_row.get("delta") is not None

            # ── Host burden ──────────────────────────────────────────────────
            burden = burden_lookup.get(primary, {})
            node_data["burden_system_category"] = burden.get("system_category", "")
            node_data["burden_functional_category"] = burden.get("functional_category", "")
            node_data["burden_flag"] = burden.get("burden_flag", "unknown")

        if add_discovered:
            for imod in hits:
                source_candidates = sorted(
                    {
                        graph_candidate_gene_ids[cand_id]
                        for cand_id in membership.index[membership[imod]].tolist()
                        if cand_id in graph_candidate_gene_ids
                    }
                )
                for discovered_id in membership.index[membership[imod]].tolist():
                    if discovered_id == gene_id or discovered_id in graph_candidate_gene_ids:
                        continue
                    gene_row = ica.gene_table.loc[discovered_id]
                    discovered_rows.append(
                        {
                            "gene_id": discovered_id,
                            "gene_name": gene_row.get("gene_name", ""),
                            "imodulon": imod,
                            "weight": float(ica.M.at[discovered_id, imod]),
                            "source_candidates": ", ".join(source_candidates),
                        }
                    )

    imodulon_to_candidates: dict[str, list[str]] = {}
    for node in candidate_nodes:
        node_data = G.nodes[node]
        for imod in node_data.get("imodulons", []):
            imodulon_to_candidates.setdefault(imod, []).append(node)

    co_edge_count = 0
    for imodulon, nodes in imodulon_to_candidates.items():
        if len(nodes) < 2:
            continue
        for src, tgt in combinations(sorted(set(nodes)), 2):
            if G.has_edge(src, tgt) and G.edges[src, tgt].get("edge_type") != "co-imodulon":
                continue
            if G.has_edge(src, tgt):
                existing = G.edges[src, tgt]
                imod_list = existing.setdefault("imodulons", [])
                if imodulon not in imod_list:
                    imod_list.append(imodulon)
                existing["n_shared_imodulons"] = len(imod_list)
            else:
                G.add_edge(
                    src,
                    tgt,
                    edge_type="co-imodulon",
                    source="precise1k",
                    imodulons=[imodulon],
                    n_shared_imodulons=1,
                )
                co_edge_count += 1

    candidate_rows = []
    for node in candidate_nodes:
        data = G.nodes[node]
        candidate_rows.append(
            {
                "gene": node,
                "group": data.get("group", ""),
                "precise_gene_id": data.get("precise_gene_id", ""),
                "imodulon_mapping_status": data.get("imodulon_mapping_status", ""),
                "imodulons": _format_modulon_list(data.get("imodulons", [])),
                "imodulon_weights": _format_weight_map(data.get("imodulon_weights", {})),
                "imodulon_count": data.get("imodulon_count", 0),
                "imodulon_primary": data.get("imodulon_primary", ""),
                "imodulon_primary_weight": data.get("imodulon_primary_weight"),
                "imodulon_primary_abs_weight": data.get("imodulon_primary_abs_weight"),
                "imodulon_primary_regulon": data.get("imodulon_primary_regulon", ""),
                "imodulon_primary_regulon_qvalue": data.get("imodulon_primary_regulon_qvalue"),
                "imodulon_primary_regulon_precision": data.get("imodulon_primary_regulon_precision"),
                "imodulon_primary_regulon_recall": data.get("imodulon_primary_regulon_recall"),
                "imodulon_primary_regulon_f1score": data.get("imodulon_primary_regulon_f1score"),
                "dr_class": data.get("dr_class", ""),
                "dr_basal": data.get("dr_basal"),
                "dr_induced": data.get("dr_induced"),
                "dr_delta": data.get("dr_delta"),
                "dr_n_treated": data.get("dr_n_treated", 0),
                "dr_n_control": data.get("dr_n_control", 0),
                "dr_data_available": data.get("dr_data_available", False),
                "burden_system_category": data.get("burden_system_category", ""),
                "burden_functional_category": data.get("burden_functional_category", ""),
                "burden_flag": data.get("burden_flag", "unknown"),
            }
        )

    candidate_df = pd.DataFrame(candidate_rows)
    discovered_df = pd.DataFrame(discovered_rows) if add_discovered and discovered_rows else None
    return G, candidate_df, discovered_df, co_edge_count


def save_outputs(G: nx.DiGraph,
                 candidate_df: pd.DataFrame,
                 enrichment_df: pd.DataFrame,
                 out_dir: str,
                 graph_path: str,
                 discovered_df: pd.DataFrame | None,
                 summary_lines: list[str]) -> None:
    os.makedirs(out_dir, exist_ok=True)

    with open(graph_path, "wb") as fh:
        pickle.dump(G, fh)
    print(f"[save] graph → {graph_path}")

    node_df = pd.DataFrame([{"gene": node, **data} for node, data in G.nodes(data=True)])
    node_path = os.path.join(out_dir, "imodulon_node_table.csv")
    node_df.to_csv(node_path, index=False)
    print(f"[save] node table ({len(node_df)} rows) → {node_path}")

    cand_path = os.path.join(out_dir, "imodulon_candidate_annotations.csv")
    candidate_df.to_csv(cand_path, index=False)
    print(f"[save] candidate annotations ({len(candidate_df)} rows) → {cand_path}")

    enrich_path = os.path.join(out_dir, "imodulon_regulon_enrichment.csv")
    enrichment_df.to_csv(enrich_path, index=False)
    print(f"[save] enrichment table ({len(enrichment_df)} rows) → {enrich_path}")

    if discovered_df is not None:
        discovered_path = os.path.join(out_dir, "imodulon_discovered_candidates.csv")
        discovered_df.to_csv(discovered_path, index=False)
        print(f"[save] discovered candidates ({len(discovered_df)} rows) → {discovered_path}")

    summary_path = os.path.join(out_dir, "imodulon_analysis_summary.txt")
    with open(summary_path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(summary_lines).rstrip() + "\n")
    print(f"[save] summary → {summary_path}")


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--graph", default="output/regulatory_network.pkl",
                   help="Pickled graph from build_network.py")
    p.add_argument("--precise", default="../data/precise1k.json.gz",
                   help="PRECISE-1K IcaData JSON.gz")
    p.add_argument("--mapping", default="../caz_kan_DE/gene_mapping.csv",
                   help="Gene mapping CSV from the DE pipeline")
    p.add_argument("--out", default="output/",
                   help="Output directory")
    p.add_argument("--out-graph", default="output/regulatory_network_imodulon.pkl",
                   help="Annotated graph output path")
    p.add_argument("--fdr", type=float, default=0.01,
                   help="FDR threshold for iModulon regulon enrichment")
    p.add_argument("--max-regs", type=int, default=1,
                   help="Maximum regulators to include in enrichment tests")
    p.add_argument("--method", default="both", choices=["and", "or", "both"],
                   help="How to combine complex regulons when max-regs > 1")
    p.add_argument("--add-discovered", action="store_true",
                   help="Write a report of extra genes found in the same iModulons")
    return p.parse_args()


def main():
    args = parse_args()

    for path, label in [
        (args.graph, "--graph"),
        (args.precise, "--precise"),
        (args.mapping, "--mapping"),
    ]:
        if not os.path.exists(path):
            sys.exit(f"[error] {label} file not found: {path}")

    G = pickle.load(open(args.graph, "rb"))
    ica = load_precise_model(args.precise)
    mapping = load_gene_mapping(args.mapping)
    lookup = build_gene_lookup(ica, mapping)
    membership, members_by_imodulon = build_membership_tables(ica)
    enrichment_df, imodulon_summary = summarize_imodulon_enrichment(
        ica,
        members_by_imodulon,
        fdr=args.fdr,
        max_regs=args.max_regs,
        method=args.method,
    )

    print("\n[dynamic-range] matching PRECISE-1K conditions to DE classes …")
    condition_groups = build_condition_groups(ica)
    dynamic_range = compute_dynamic_range(ica, condition_groups)
    burden_lookup = host_burden_lookup(ica.imodulon_table)

    G, candidate_df, discovered_df, co_edge_count = annotate_graph(
        G,
        ica,
        lookup,
        membership,
        imodulon_summary,
        dynamic_range=dynamic_range,
        burden_lookup=burden_lookup,
        add_discovered=args.add_discovered,
    )

    mapped_candidates = candidate_df["precise_gene_id"].astype(str).str.len().gt(0).sum()
    enriched_imodulons = sum(1 for v in imodulon_summary.values() if v.get("qvalue") is not None)
    dr_covered = int(candidate_df["dr_data_available"].sum()) if "dr_data_available" in candidate_df else 0
    burden_risk = int((candidate_df["burden_flag"] == "risk").sum()) if "burden_flag" in candidate_df else 0
    summary_lines = [
        "REGULATORY NETWORK iModulon SUMMARY",
        f"candidate nodes annotated: {len(candidate_df)}",
        f"candidate nodes mapped to PRECISE: {mapped_candidates}",
        f"iModulons with significant regulon enrichment: {enriched_imodulons}",
        f"co-iModulon edges added: {co_edge_count}",
        f"candidates with dynamic-range data available: {dr_covered}/{len(candidate_df)}",
        f"candidates flagged host-burden risk (Metabolism/Translation iModulon): {burden_risk}",
        f"output graph: {args.out_graph}",
    ]

    save_outputs(
        G=G,
        candidate_df=candidate_df,
        enrichment_df=enrichment_df,
        out_dir=args.out,
        graph_path=args.out_graph,
        discovered_df=discovered_df,
        summary_lines=summary_lines,
    )


if __name__ == "__main__":
    main()

