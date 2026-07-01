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
    with gzip.open(path, "rt", encoding="utf-8") as fh:
        serial = json.load(fh)

    for key in MATRIX_COLUMNS:
        if key in serial and serial[key] is not None:
            serial[key] = _read_json_table(serial[key])

    return IcaData(**serial)


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
                   add_discovered: bool = False) -> tuple[nx.DiGraph, pd.DataFrame, pd.DataFrame | None, int]:
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

    G, candidate_df, discovered_df, co_edge_count = annotate_graph(
        G,
        ica,
        lookup,
        membership,
        imodulon_summary,
        add_discovered=args.add_discovered,
    )

    mapped_candidates = candidate_df["precise_gene_id"].astype(str).str.len().gt(0).sum()
    enriched_imodulons = sum(1 for v in imodulon_summary.values() if v.get("qvalue") is not None)
    summary_lines = [
        "REGULATORY NETWORK iModulon SUMMARY",
        f"candidate nodes annotated: {len(candidate_df)}",
        f"candidate nodes mapped to PRECISE: {mapped_candidates}",
        f"iModulons with significant regulon enrichment: {enriched_imodulons}",
        f"co-iModulon edges added: {co_edge_count}",
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
