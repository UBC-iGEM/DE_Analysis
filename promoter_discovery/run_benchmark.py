"""Fit one count-aware model for GSE220559 and export matched drug contrasts.

Shared water controls occur once in the model. Pairwise drug comparisons are
also exported: significance in one control contrast alone is not specificity.
"""
from __future__ import annotations

import argparse
import importlib.metadata
import itertools
import json
from pathlib import Path

import numpy as np
import pandas as pd

from promoter_discovery.input_data import load_tar_gene_tables
from .dataset_registry import load_dataset_config
from .layout import contrast_directory


def assemble_counts(config):
    """Validate raw integer counts and merge identical shared samples once."""
    samples, records = {}, {}
    for dataset in config["datasets"]:
        expression, metadata = load_tar_gene_tables(dataset)
        if expression.gene_id.duplicated().any():
            raise ValueError("Duplicate gene IDs in raw count table")
        expression = expression.set_index("gene_id")
        for row in metadata.to_dict("records"):
            sample = row["sample_id"]
            values = pd.to_numeric(expression[sample], errors="raise")
            if not np.isfinite(values).all() or (values < 0).any() or (values != np.floor(values)).any():
                raise ValueError("DESeq2 requires finite nonnegative raw integer counts; do not round normalized counts")
            values = values.astype("int64").sort_index()
            condition = "water" if row["condition"] == "control" else dataset["name"]
            if sample in samples:
                if records[sample]["condition"] != condition or not samples[sample].equals(values):
                    raise ValueError(f"Inconsistent shared sample: {sample}")
            else:
                samples[sample] = values
                records[sample] = {"sample_id": sample, "condition": condition, "source_file": row["source_file"]}
    counts = pd.DataFrame(samples)
    if counts.isna().any().any():
        raise ValueError("Raw count tables do not share the same gene universe")
    metadata = pd.DataFrame(records.values()).set_index("sample_id").loc[counts.columns]
    if (metadata.condition.value_counts() < 3).any():
        raise ValueError("Benchmark requires at least three samples in every condition")
    return counts, metadata


def contrast_results(dds, numerator, denominator, config, shrink=False):
    from pydeseq2.ds import DeseqStats
    stats = DeseqStats(dds, contrast=["condition", numerator, denominator],
                       alpha=config["thresholds"]["padj"], independent_filter=False,
                       n_cpus=1, quiet=True)
    stats.summary()
    result = stats.results_df.copy()
    result.index.name = "gene_id"
    result = result.reset_index()
    result["gene"] = result["gene_id"]
    # Cook's-filtered tests remain unavailable, never silently become p=1.
    result["eligible_for_network"] = result.pvalue.notna() & result.padj.notna()
    result["de_method"] = "PyDESeq2 " + importlib.metadata.version("pydeseq2")
    result["padj_source"] = "PyDESeq2_BH; no independent filtering; Cook's filtering enabled"
    fc, alpha = config["thresholds"]["log2_fold_change"], config["thresholds"]["padj"]
    sig = result.padj.lt(alpha) & result.log2FoldChange.abs().gt(fc)
    result["regulation"] = np.where(sig, np.where(result.log2FoldChange > 0, "upregulated", "downregulated"), "not_regulated")
    result["lfc_ci_low"] = result.log2FoldChange - 1.96 * result.lfcSE
    result["lfc_ci_high"] = result.log2FoldChange + 1.96 * result.lfcSE
    result["signal_strength"] = result.log2FoldChange.abs() * -np.log10(result.padj.clip(lower=np.finfo(float).tiny))
    if shrink:
        coefficient = f"condition[T.{numerator}]"
        stats.lfc_shrink(coeff=coefficient)
        result["shrunk_log2FoldChange"] = stats.results_df.log2FoldChange.to_numpy()
        result["shrunk_lfcSE"] = stats.results_df.lfcSE.to_numpy()
    return result


def run(config_path, root, out):
    from pydeseq2.dds import DeseqDataSet
    from .setup_data import sha256_file
    config = load_dataset_config(config_path)
    input_config = {**config, "datasets": [{**d, "archive": str(root / d["archive"])} for d in config["datasets"]]}
    counts, metadata = assemble_counts(input_config)
    metadata["condition"] = pd.Categorical(metadata.condition, categories=["water", *[d["name"] for d in config["datasets"]]])
    settings = config["benchmark"]
    retained = counts.ge(settings["min_count"]).sum(axis=1).ge(settings["min_samples"])
    if retained.sum() < 2:
        raise ValueError("Too few genes pass the common count filter")
    dds = DeseqDataSet(counts=counts.loc[retained].T, metadata=metadata,
                       design="~condition", refit_cooks=True, n_cpus=1, quiet=True)
    dds.deseq2()
    normalized = pd.DataFrame(dds.layers["normed_counts"].T, index=dds.var_names, columns=dds.obs_names)
    out.mkdir(parents=True, exist_ok=True)
    metadata.to_csv(out / "joint_metadata.csv")
    pd.DataFrame({"gene_id": counts.index, "passes_common_count_filter": retained}).to_csv(out / "gene_filter.csv", index=False)
    normalized.to_csv(out / "joint_normalized_counts.csv", index_label="gene_id")
    counts.to_csv(out / "joint_raw_counts.csv", index_label="gene_id")
    dds.vst(use_design=True)
    pd.DataFrame(dds.layers["vst_counts"].T, index=dds.var_names, columns=dds.obs_names).to_csv(out / "joint_vst_counts.csv", index_label="gene_id")
    for dataset in config["datasets"]:
        result = contrast_results(dds, dataset["name"], "water", config, shrink=True)
        ctrl = metadata.index[metadata.condition.eq("water")]
        drug = metadata.index[metadata.condition.eq(dataset["name"])]
        result["benchmark_basal"] = result.gene_id.map(normalized[ctrl].mean(axis=1))
        result["benchmark_induced"] = result.gene_id.map(normalized[drug].mean(axis=1))
        result["benchmark_expression_units"] = "joint DESeq2 normalized counts; relative abundance, not promoter fluorescence"
        result["antibiotic_class"] = dataset["antibiotic_class"]
        result["treatment"] = dataset["treatment"]
        result["comparison"] = dataset["comparison"]
        target = contrast_directory(root, dataset)
        target.mkdir(parents=True, exist_ok=True)
        result.to_csv(target / "de_results.csv", index=False)
        counts[list(ctrl) + list(drug)].to_csv(target / "expression.csv", index_label="gene_id")
        metadata.loc[list(ctrl) + list(drug)].assign(condition=["control"] * len(ctrl) + ["treated"] * len(drug)).to_csv(target / "metadata.csv")
    pairwise = []
    for first, second in itertools.combinations(config["datasets"], 2):
        result = contrast_results(dds, first["name"], second["name"], config)
        result["numerator_dataset"], result["denominator_dataset"] = first["name"], second["name"]
        pairwise.append(result)
    pd.concat(pairwise, ignore_index=True).to_csv(out / "drug_pairwise_contrasts.csv", index=False)
    provenance = {"config_sha256": sha256_file(config_path), "archives": {d["archive"]: sha256_file(root / d["archive"]) for d in config["datasets"]},
                  "de_exports": {d["name"]: sha256_file(contrast_directory(root, d) / "de_results.csv") for d in config["datasets"]},
                  "pairwise_sha256": sha256_file(out / "drug_pairwise_contrasts.csv"),
                  "model": "~condition; one joint fit; unique shared controls", "pydeseq2_version": importlib.metadata.version("pydeseq2"),
                  "filter": settings, "n_samples": len(metadata), "n_genes_retained": int(retained.sum()),
                  "multiple_testing": "BH across eligible genes separately per contrast; no independent filtering",
                  "limitations": "Class proxies and shared controls; within-study evidence, not independent validation. Wald 95% intervals are approximate."}
    (out / "benchmark_provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    print(f"[benchmark] {len(metadata)} unique samples; {retained.sum()} retained genes; {len(pairwise)} direct drug contrasts")


def main(argv=None):
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=root / "config/benchmark.json")
    parser.add_argument("--root", type=Path, default=root)
    parser.add_argument("--out", type=Path, default=root / "results/differential_expression/model")
    args = parser.parse_args(argv)
    run(args.config, args.root, args.out)


if __name__ == "__main__":
    main()
