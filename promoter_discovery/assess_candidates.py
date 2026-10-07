"""Assess expression quality, candidate stability, regulatory context, and a review panel."""
from __future__ import annotations

import argparse
from collections import Counter
import itertools
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd

from .promoter_selection import read_product
from .regulon_enrichment import bh_adjust


CLASSES = ("beta_lactam", "aminoglycoside")


def sample_quality(raw, vst, metadata, out):
    samples = metadata.index
    matrix = vst[samples].T.to_numpy()
    matrix -= matrix.mean(axis=0)
    u, singular, _ = np.linalg.svd(matrix, full_matrices=False)
    variance = singular**2 / np.sum(singular**2)
    qc = metadata.copy()
    qc["library_counts"] = raw[samples].sum()
    qc["genes_detected"] = raw[samples].gt(0).sum()
    qc[["PC1", "PC2"]] = u[:, :2] * singular[:2]
    qc.to_csv(out / "sample_qc.csv", index_label="sample_id")
    fig, ax = plt.subplots(figsize=(10, 6))
    for condition, group in qc.groupby("condition", observed=True):
        ax.scatter(group.PC1, group.PC2, label=condition.replace("benchmark_", ""))
        for i, (sample, row) in enumerate(group.iterrows()):
            dx, dy = [(-8, 10), (8, 0), (-8, -14)][i % 3]
            ax.annotate(sample, (row.PC1, row.PC2), fontsize=8, xytext=(dx, dy), textcoords="offset points", ha="right" if dx < 0 else "left")
    ax.set(xlabel=f"PC1 ({variance[0]:.1%})", ylabel=f"PC2 ({variance[1]:.1%})", title="Sample PCA — variance-stabilized counts")
    ax.legend(fontsize=8)
    ax.margins(x=0.12, y=0.1)
    fig.tight_layout()
    fig.savefig(out / "sample_pca.png", dpi=160)
    plt.close(fig)
    distances = np.linalg.norm(matrix[:, None] - matrix[None, :], axis=2)
    fig, ax = plt.subplots(figsize=(9, 8))
    image = ax.imshow(distances, cmap="viridis")
    ax.set(xticks=range(len(samples)), yticks=range(len(samples)), xticklabels=samples, yticklabels=samples, title="Sample distances — variance-stabilized counts")
    plt.setp(ax.get_xticklabels(), rotation=90)
    fig.colorbar(image, ax=ax, label="Euclidean distance")
    fig.tight_layout()
    fig.savefig(out / "sample_distances.png", dpi=160)
    plt.close(fig)
    return qc


def replicate_summary(observations, normalized, metadata):
    records = []
    for row in observations.itertuples(index=False):
        for condition in ("water", row.source_dataset):
            values = normalized.loc[row.source_gene_id, metadata.index[metadata.condition.eq(condition)]]
            records.append({"gene": row.canonical_gene, "gene_id": row.source_gene_id,
                            "source_dataset": row.source_dataset, "condition": condition,
                            "mean": values.mean(), "sd": values.std(ddof=1),
                            "cv": values.std(ddof=1) / values.mean() if values.mean() > 0 else np.nan,
                            "minimum": values.min(), "maximum": values.max(),
                            "sample_values_json": json.dumps(values.to_dict())})
    return pd.DataFrame(records)


def candidate_stability(observations, raw, config):
    """Vary post-fit screening, keeping the fitted model and p-values fixed."""
    alpha = config["thresholds"]["padj"]
    settings = config["benchmark"]
    count_limits = sorted({settings["min_count"], 2 * settings["min_count"], 5 * settings["min_count"]})
    effect_limits = [config["thresholds"]["log2_fold_change"] + delta for delta in (-0.5, 0, 0.5)]
    rows = []
    discovery = observations[observations.dataset_role.ne("cross_reactivity_challenge")].copy()
    for cls, group in discovery.groupby("antibiotic_class"):
        qualified = group.eligible_for_network & group.padj.lt(alpha) & group.log2FoldChange.gt(config["thresholds"]["log2_fold_change"])
        candidates = set(group.loc[qualified, "canonical_gene"])
        selected = Counter()
        for count, effect in itertools.product(count_limits, effect_limits):
            passes = raw.ge(count).sum(axis=1).ge(settings["min_samples"])
            keep = group.eligible_for_network & group.padj.lt(alpha) & group.log2FoldChange.gt(effect)
            keep &= group.source_gene_id.map(passes).fillna(False)
            selected.update(set(group.loc[keep, "canonical_gene"]))
        for gene in sorted(candidates):
            hits = group[qualified & group.canonical_gene.eq(gene)]
            rows.append({"gene": gene, "antibiotic_class": cls, "n_supporting_comparisons": len(hits),
                         "raw_effect": hits.log2FoldChange.max(), "shrunk_effect": hits.shrunk_log2FoldChange.max(),
                         "n_screens_selected": selected[gene], "n_screens": len(count_limits) * len(effect_limits),
                         "screen_fraction": selected[gene] / (len(count_limits) * len(effect_limits))})
    result = pd.DataFrame(rows)
    for effect, rank in [("raw_effect", "raw_rank"), ("shrunk_effect", "shrunk_rank")]:
        ordered = result.sort_values(["antibiotic_class", "n_supporting_comparisons", effect, "gene"], ascending=[True, False, False, True])
        result.loc[ordered.index, rank] = ordered.groupby("antibiotic_class").cumcount() + 1
    result["rank_change"] = result.raw_rank - result.shrunk_rank
    return result.sort_values(["antibiotic_class", "shrunk_rank"]).reset_index(drop=True)


def cross_reactivity(observations, pairs, config, off_target_limit):
    """Orient direct contrasts and adjust over all genes and ten drug pairs."""
    pairs = pairs.copy()
    valid = pairs.pvalue.notna()
    pairs["family_padj"] = np.nan
    pairs.loc[valid, "family_padj"] = bh_adjust(pairs.loc[valid, "pvalue"])
    pair_lookup = pairs.set_index(["gene_id", "numerator_dataset", "denominator_dataset"])
    obs = observations.set_index(["source_gene_id", "source_dataset"])
    alpha = config["thresholds"]["padj"]
    cutoff = config["thresholds"]["log2_fold_change"]
    seeds = observations[observations.dataset_role.ne("cross_reactivity_challenge") & observations.eligible_for_network
                         & observations.padj.lt(alpha) & observations.log2FoldChange.gt(cutoff)]
    records = []
    for seed in seeds.itertuples(index=False):
        for other in config["datasets"]:
            if other["antibiotic_class"] == seed.antibiotic_class:
                continue
            key = (seed.source_gene_id, seed.source_dataset, other["name"])
            reverse = key not in pair_lookup.index
            pair = pair_lookup.loc[(key[0], key[2], key[1]) if reverse else key]
            sign = -1 if reverse else 1
            effect = sign * pair.log2FoldChange
            low, high = (-pair.lfc_ci_high, -pair.lfc_ci_low) if reverse else (pair.lfc_ci_low, pair.lfc_ci_high)
            status = "inconclusive"
            if pair.eligible_for_network and pair.family_padj < alpha:
                status = "discovery_response_stronger" if low > 0 else "counterpart_response_stronger" if high < 0 else status
            challenge = obs.loc[(seed.source_gene_id, other["name"])]
            bounded = challenge.eligible_for_network and challenge.lfc_ci_low > -off_target_limit and challenge.lfc_ci_high < off_target_limit
            records.append({"gene": seed.canonical_gene, "antibiotic_class": seed.antibiotic_class,
                            "discovery_dataset": seed.source_dataset, "counterpart_dataset": other["name"],
                            "counterpart_role": other["role"], "direct_log2_difference": effect,
                            "direct_ci_low": low, "direct_ci_high": high,
                            "per_contrast_padj": pair.padj, "family_padj": pair.family_padj,
                            "direct_status": status, "counterpart_log2FC": challenge.log2FoldChange,
                            "counterpart_padj": challenge.padj, "counterpart_ci_low": challenge.lfc_ci_low,
                            "counterpart_ci_high": challenge.lfc_ci_high,
                            "counterpart_response_detected": bool(challenge.eligible_for_network and challenge.padj < alpha and abs(challenge.log2FoldChange) > cutoff),
                            "off_target_log2_limit": off_target_limit, "counterpart_interval_within_bound": bool(bounded)})
    return pd.DataFrame(records)


def transcription_groups(tus):
    """Connect overlapping annotated TUs; these groups are descriptive units."""
    graph = nx.Graph()
    for genes in tus.tuGenes:
        names = [gene.strip().lower() for gene in genes.split(";") if gene.strip()]
        graph.add_nodes_from(names)
        graph.add_edges_from(zip(names, names[1:]))
    return {gene: ";".join(sorted(component)) for component in nx.connected_components(graph) for gene in component}


def operon_support(enrichment, groups):
    result = enrichment.copy()
    context = []
    for row in enrichment.itertuples(index=False):
        genes = [gene.strip() for gene in str(row.overlap_genes).split(",") if gene.strip() and gene.strip() != "nan"]
        counts = Counter(groups.get(gene, "unannotated:" + gene) for gene in genes)
        dominant, count = counts.most_common(1)[0] if counts else ("", 0)
        context.append({"n_transcription_groups": len(counts), "n_genes_without_tu": sum(gene not in groups for gene in genes),
                        "dominant_transcription_group": dominant, "dominant_group_gene_fraction": count / len(genes) if genes else np.nan})
    return pd.concat([result.reset_index(drop=True), pd.DataFrame(context)], axis=1)


def select_panel(options, per_class):
    """Greedy, balanced selection with unique promoters and TU groups."""
    chosen, promoters, groups, regulators = [], set(), set(), Counter()
    for _ in range(per_class):
        for cls in CLASSES:
            available = options[options.antibiotic_class.eq(cls) & options.reference_sequence_matches
                                & ~options.promoter_id.isin(promoters) & ~options.transcription_group.isin(groups)].copy()
            if available.empty:
                continue
            available["regulator_reuse"] = available.regulators.fillna("").map(lambda names: sum(regulators[name.strip()] for name in names.split(",") if name.strip()) if names else np.nan)
            available = available.sort_values(["n_supporting_comparisons", "screen_fraction", "confidence_rank", "regulator_reuse",
                                                "direct_stronger_fraction", "challenge_max_abs_log2FC", "shrunk_effect", "gene", "promoter_id"],
                                               ascending=[False, False, False, True, False, True, False, True, True])
            row = available.iloc[0].to_dict()
            row["selection_order"] = len(chosen) + 1
            row["selection_basis"] = "class support; screen stability; annotation confidence; regulator diversity; direct comparisons; challenge effect; shrunk effect"
            row["construct_status"] = "proposed_panel; requires_operator_and_fragment_review_and_wetlab_validation"
            chosen.append(row)
            promoters.add(row["promoter_id"])
            groups.add(row["transcription_group"])
            regulators.update(name.strip() for name in str(row["regulators"]).split(",") if name.strip() and name.strip() != "nan")
    return pd.DataFrame(chosen, columns=[*options.columns, "regulator_reuse", "selection_order", "selection_basis", "construct_status"])


def oriented_sequence(genome, start, end, strand):
    sequence = genome[start - 1:end].upper()
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1] if strand == "reverse" else sequence


def fragment_review(review, promoters, interactions, genes, genome, upstream, downstream):
    """Propose TSS windows extended to include recorded promoter-linked sites."""
    annotations = promoters.set_index("id")
    sites = interactions[interactions.promoterID.isin(review.promoter_id)].copy()
    for column in ("tfrsLeft", "tfrsRight"):
        sites[column] = pd.to_numeric(sites[column], errors="coerce")
    sites = sites[sites.tfrsLeft.gt(0) & sites.tfrsRight.ge(sites.tfrsLeft)].copy()
    genes = genes.copy()
    genes["leftEndPos"] = pd.to_numeric(genes.leftEndPos, errors="coerce")
    genes["rightEndPos"] = pd.to_numeric(genes.rightEndPos, errors="coerce")
    records = []
    for row in review.itertuples(index=False):
        promoter = annotations.loc[row.promoter_id]
        if not promoter.posTSS or promoter.strand not in ("forward", "reverse"):
            records.append({"promoter_id": row.promoter_id, "promoter_name": row.promoter_name,
                            "candidate_genes": row.candidate_genes, "reference_sequence_matches": False,
                            "construct_status": "missing_tss_or_strand; manual_annotation_review_required"})
            continue
        tss, strand = int(promoter.posTSS), promoter.strand
        left, right = (tss - upstream, tss + downstream) if strand == "forward" else (tss - downstream, tss + upstream)
        linked = sites[sites.promoterID.eq(row.promoter_id)]
        outside = linked.tfrsLeft.lt(left) | linked.tfrsRight.gt(right)
        start = max(1, int(min(left, linked.tfrsLeft.min())) if not linked.empty else left)
        end = min(len(genome), int(max(right, linked.tfrsRight.max())) if not linked.empty else right)
        sequence = oriented_sequence(genome, start, end, strand)
        annotated = promoter.sequence
        anchor = next((i for i, base in enumerate(annotated) if base.isupper()), None)
        match = False
        if anchor is not None:
            a, b = (tss - anchor, tss + len(annotated) - anchor - 1) if strand == "forward" else (tss - len(annotated) + anchor + 1, tss + anchor)
            match = oriented_sequence(genome, a, b, strand) == annotated.upper()
        overlapping = genes[genes.leftEndPos.le(end) & genes.rightEndPos.ge(start)]
        records.append({"promoter_id": row.promoter_id, "promoter_name": row.promoter_name,
                        "candidate_genes": row.candidate_genes, "strand": strand, "tss": tss,
                        "sigma_factor": row.sigma_factor, "promoter_confidence": row.promoter_confidence,
                        "fragment_start": start, "fragment_end": end, "fragment_length": len(sequence),
                        "fragment_sequence": sequence, "genome_accession": "U00096.3",
                        "reference_sequence_matches": match, "n_curated_sites": linked.tfrsID.nunique(),
                        "n_sites_outside_initial_window": linked.loc[outside, "tfrsID"].nunique(),
                        "site_regulators": ", ".join(sorted(set(linked.regulatorName))),
                        "overlapping_genes": ", ".join(sorted(set(overlapping.geneName))),
                        "minus10_position": promoter.boxMinus10pos, "minus35_position": promoter.boxMinus35pos,
                        "construct_status": "proposed_fragment; manual_boundary_and_operator_review_required"})
    return pd.DataFrame(records), sites


def panel_options(stability, candidates, mapped, fragments, reactivity, replicates, groups):
    summary = reactivity.groupby(["gene", "antibiotic_class"]).agg(
        direct_stronger_fraction=("direct_status", lambda s: s.eq("discovery_response_stronger").mean()),
        challenge_max_abs_log2FC=("counterpart_log2FC", lambda s: s.abs().max()))
    variability = replicates[replicates.condition.ne("water")].groupby("gene").cv.max().rename("max_treated_cv")
    result = stability.merge(summary, on=["gene", "antibiotic_class"], how="left").merge(variability, on="gene", how="left")
    result = result.merge(candidates[["gene", "regulators"]], on="gene")
    result = result.merge(mapped[mapped.mapping_status.eq("curated_promoter_mapped")][["gene", "promoter_id"]].drop_duplicates(), on="gene")
    result = result.merge(fragments, on="promoter_id")
    result["transcription_group"] = result.gene.map(groups).fillna(result.gene)
    result["confidence_rank"] = result.promoter_confidence.map({"C": 3, "S": 2, "W": 1}).fillna(0)
    return result


def panel_figures(panel, observations, normalized, metadata, sites, out):
    if panel.empty:
        return
    names = panel.gene.tolist()
    effects = observations[observations.canonical_gene.isin(names)].pivot(index="canonical_gene", columns="treatment", values="log2FoldChange").reindex(names)
    limit = np.nanmax(np.abs(effects.to_numpy()))
    fig, ax = plt.subplots(figsize=(9, max(4, len(names) * 0.35)))
    image = ax.imshow(effects, cmap="RdBu_r", vmin=-limit, vmax=limit, aspect="auto")
    ax.set(xticks=range(len(effects.columns)), xticklabels=effects.columns, yticks=range(len(names)), yticklabels=names, title="Proposed panel — original drug/control effect sizes")
    for i in range(len(names)):
        for j in range(len(effects.columns)):
            ax.text(j, i, f"{effects.iloc[i, j]:.1f}", ha="center", va="center", fontsize=8)
    fig.colorbar(image, ax=ax, label="log₂ fold change")
    fig.tight_layout()
    fig.savefig(out / "panel_response_heatmap.png", dpi=160)
    plt.close(fig)
    conditions = metadata.condition.unique().tolist()
    fig, axes = plt.subplots(len(names), 1, figsize=(10, max(5, len(names) * 1.4)), squeeze=False, sharex=True)
    ids = observations.drop_duplicates("canonical_gene").set_index("canonical_gene").source_gene_id
    for gene, ax in zip(names, axes[:, 0]):
        for i, condition in enumerate(conditions):
            values = normalized.loc[ids[gene], metadata.index[metadata.condition.eq(condition)]]
            ax.scatter(i + np.linspace(-0.12, 0.12, len(values)), np.log2(values + 1), s=18)
        ax.set_ylabel(gene, rotation=0, ha="right", fontsize=9)
    axes[0, 0].set_title("Individual replicates — log₂(normalized counts + 1)")
    axes[-1, 0].set(xticks=range(len(conditions)), xticklabels=[str(c).replace("benchmark_", "") for c in conditions])
    fig.tight_layout()
    fig.savefig(out / "panel_replicates.png", dpi=160)
    plt.close(fig)
    fig, axes = plt.subplots(len(panel), 1, figsize=(11, max(5, len(panel) * 1.5)), squeeze=False)
    for row, ax in zip(panel.itertuples(index=False), axes[:, 0]):
        sign = -1 if row.strand == "reverse" else 1
        endpoints = sorted((sign * (row.fragment_start - row.tss), sign * (row.fragment_end - row.tss)))
        ax.plot(endpoints, [0, 0], color="grey")
        ax.axvline(0, color="black", linestyle=":")
        for column, color, level in [("minus10_position", "tab:green", 0.1), ("minus35_position", "tab:blue", 0.2)]:
            value = getattr(row, column)
            if isinstance(value, str) and value:
                positions = sorted(sign * (int(p) - row.tss) for p in value.split("-"))
                ax.plot(positions, [level, level], color=color, linewidth=5)
        linked = sites[sites.promoterID.eq(row.promoter_id)].drop_duplicates(["tfrsID", "regulatorName"])
        for i, (regulator, group) in enumerate(linked.groupby("regulatorName")):
            level = -0.2 - (i % 3) * 0.3
            centers = []
            for site in group.itertuples(index=False):
                positions = sorted((sign * (site.tfrsLeft - row.tss), sign * (site.tfrsRight - row.tss)))
                ax.plot(positions, [level, level], color="tab:orange", linewidth=5)
                centers.append(np.mean(positions))
            ax.text(np.mean(centers), level - 0.18, regulator, ha="center", fontsize=7)
        ax.set(ylim=(-1.1, 0.4), yticks=[], title=f"{row.promoter_name} ({row.gene}); {row.strand}; {int(row.fragment_length)} bp")
    axes[0, 0].text(0.01, 1.5, "TSS: dotted line | −10: green | −35: blue | curated regulatory sites: orange", transform=axes[0, 0].transAxes, fontsize=9)
    axes[-1, 0].set_xlabel("Position relative to TSS, in transcription direction (bp)")
    fig.tight_layout()
    fig.savefig(out / "panel_promoter_context.png", dpi=160)
    plt.close(fig)


def main(argv=None):
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=root)
    parser.add_argument("--out", type=Path, default=None)
    parser.add_argument("--panel-per-class", type=int, default=6)
    parser.add_argument("--off-target-limit", type=float, default=1.0, help="Descriptive log2 effect bound, default twofold; not a specificity test")
    parser.add_argument("--upstream", type=int, default=200)
    parser.add_argument("--downstream", type=int, default=30)
    parser.add_argument("--skip-database", action="store_true", help="do not rebuild the generated SQLite database")
    args = parser.parse_args(argv)
    root = args.root
    out = args.out or root / "results/candidate_assessment"
    out.mkdir(parents=True, exist_ok=True)
    model = root / "results/differential_expression/model"
    network = root / "results/regulatory_network"
    promoter_dir = root / "results/promoter_review"
    references = root / "data/references/regulondb"
    config = json.loads((root / "config/benchmark.json").read_text())
    raw = pd.read_csv(model / "joint_raw_counts.csv", index_col=0)
    normalized = pd.read_csv(model / "joint_normalized_counts.csv", index_col=0)
    vst = pd.read_csv(model / "joint_vst_counts.csv", index_col=0)
    metadata = pd.read_csv(model / "joint_metadata.csv", index_col=0)
    observations = pd.read_csv(network / "source_observations.csv")
    candidates = pd.read_csv(root / "results/promoter_candidates/annotated_candidates.csv")
    mapped = pd.read_csv(promoter_dir / "candidate_promoter_mapping.csv")
    review = pd.read_csv(promoter_dir / "promoter_review.csv")
    tus = read_product(references / "TUSet.tsv")
    promoters = read_product(references / "PromoterSet.tsv")
    interactions = read_product(references / "RISet.tsv")
    genes = read_product(references / "GeneProductAllIdentifiersSet.tsv")
    genome_path = root / "data/references/ecoli_mg1655/U00096.3.fasta"
    genome = "".join(line.strip() for line in genome_path.read_text().splitlines() if not line.startswith(">"))
    sample_quality(raw, vst, metadata, out)
    replicates = replicate_summary(observations[observations.canonical_gene.isin(candidates.gene)], normalized, metadata)
    stability = candidate_stability(observations, raw, config)
    reactivity = cross_reactivity(observations, pd.read_csv(model / "drug_pairwise_contrasts.csv"), config, args.off_target_limit)
    groups = transcription_groups(tus)
    operons = operon_support(pd.read_csv(network / "regulon_enrichment.csv"), groups)
    fragments, sites = fragment_review(review, promoters, interactions, genes, genome, args.upstream, args.downstream)
    options = panel_options(stability, candidates, mapped, fragments, reactivity, replicates, groups)
    panel = select_panel(options, args.panel_per_class)
    for name, frame in [("candidate_replicates", replicates), ("candidate_stability", stability), ("cross_reactivity", reactivity),
                        ("operon_support", operons), ("fragment_review", fragments), ("regulatory_sites", sites), ("experimental_panel", panel)]:
        frame.to_csv(out / f"{name}.csv", index=False)
    with (out / "proposed_fragments.fasta").open("w") as handle:
        for row in fragments.itertuples(index=False):
            if pd.notna(row.fragment_sequence):
                handle.write(f">{row.promoter_id} {row.promoter_name} U00096.3:{int(row.fragment_start)}-{int(row.fragment_end)} {row.strand}; proposed_not_ready_to_order\n{row.fragment_sequence}\n")
    panel_figures(panel, observations, normalized, metadata, sites, out)
    summary = {"n_candidates": len(candidates), "n_panel_promoters": len(panel), "panel_per_class_requested": args.panel_per_class,
               "n_fragment_reviews": len(fragments), "n_reference_sequence_matches": int(fragments.reference_sequence_matches.sum()),
               "n_promoters_missing_coordinates": int(fragments.fragment_start.isna().sum()),
               "off_target_log2_limit": args.off_target_limit, "initial_window": {"upstream": args.upstream, "downstream": args.downstream},
               "screening": "Nine post-fit screens: min count x1/x2/x5; effect cutoff -0.5/0/+0.5; fixed p-values and fitted gene universe; no model refits",
               "direct_comparisons": "BH across all finite tests in all ten drug pairs before candidate selection",
               "operon_support": "Descriptive connected groups of overlapping curated TUs; enrichment p-values are unchanged",
               "panel_selection": "Greedy review proposal; equal class quotas; no repeated promoter or connected TU group; not an optimized or validated biosensor panel",
               "limitations": "Class proxies, three replicates per condition, shared controls. Expression is not reporter fluorescence; CI bounds are descriptive, not equivalence tests. Fragments require manual review."}
    (out / "analysis_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    if not args.skip_database:
        from .build_database import build_database

        database_path = build_database(root)
        print(f"[database] {database_path}")
    print(f"[assessment] {len(candidates)} candidates; {len(fragments)} fragment reviews; {len(panel)} proposed panel promoters")
    print(f"[save] {out}")


if __name__ == "__main__":
    main()
