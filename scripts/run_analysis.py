import argparse
import gzip
import io
import json
import math
import os
import re
import tarfile
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[0]
OUTPUT_ROOT = REPO_ROOT / "outputs"
MPL_CACHE_DIR = OUTPUT_ROOT / "_cache" / "matplotlib"
MPL_CACHE_DIR.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(MPL_CACHE_DIR))


LEADING_COLUMNS = [
    "antibiotic_class",
    "treatment",
    "comparison",
    "regulation",
    "gene",
    "gene_id",
    "baseMean",
    "log2FoldChange",
    "FoldChange",
    "pvalue",
    "padj",
    "signal_strength",
    "padj_source",
]


def repo_path(path):
    return REPO_ROOT / path


def benjamini_hochberg(pvalues):
    pvalues = pd.Series(pd.to_numeric(pvalues, errors="coerce")).fillna(1.0)
    ranked = pvalues.rank(method="first").astype(int)
    adjusted = pvalues * len(pvalues) / ranked
    adjusted = adjusted.sort_values(ascending=False).cummin().sort_index()
    return adjusted.clip(upper=1.0)


def signal_strength(log2fc, padj):
    if pd.isna(log2fc):
        return np.nan
    if pd.isna(padj):
        padj = 1.0
    safe_padj = max(float(padj), np.finfo(float).tiny)
    return abs(float(log2fc)) * -np.log10(safe_padj)


def classify(log2fc, padj, log2fc_threshold, padj_threshold):
    if pd.isna(log2fc) or pd.isna(padj):
        return "not_regulated"
    if padj < padj_threshold and log2fc > log2fc_threshold:
        return "upregulated"
    if padj < padj_threshold and log2fc < -log2fc_threshold:
        return "downregulated"
    return "not_regulated"


def readable(df):
    rounded = df.copy()
    numeric = rounded.select_dtypes(include=[np.number]).columns
    rounded[numeric] = rounded[numeric].round(2)
    return rounded


def order_columns(df):
    return df[
        [column for column in LEADING_COLUMNS if column in df.columns]
        + [column for column in df.columns if column not in LEADING_COLUMNS]
    ]


def finalize_results(results, dataset, thresholds):
    results = results.copy()
    results["antibiotic_class"] = dataset["antibiotic_class"]
    results["treatment"] = dataset["treatment"]
    results["comparison"] = results.get(
        "comparison",
        f"{dataset['treatment']}_vs_control",
    )
    results["log2FoldChange"] = pd.to_numeric(
        results["log2FoldChange"], errors="coerce")
    results["pvalue"] = pd.to_numeric(results.get("pvalue"), errors="coerce")
    results["padj"] = pd.to_numeric(results.get("padj"), errors="coerce")

    fill_padj = benjamini_hochberg(results["pvalue"])
    results["padj_source"] = np.where(
        results["padj"].notna(), "input", "BH_from_pvalue")
    results["padj"] = results["padj"].fillna(fill_padj).fillna(1.0)
    results.loc[results["pvalue"].isna(), "padj_source"] = (
        "no_valid_pvalue_set_to_1"
    )

    results["signal_strength"] = [
        signal_strength(fc, adj)
        for fc, adj in zip(results["log2FoldChange"], results["padj"])
    ]
    results["regulation"] = [
        classify(
            fc,
            adj,
            thresholds["log2_fold_change"],
            thresholds["padj"],
        )
        for fc, adj in zip(results["log2FoldChange"], results["padj"])
    ]
    results = results.dropna(subset=["log2FoldChange"])
    results = results.sort_values(
        ["signal_strength", "log2FoldChange"], ascending=[False, False])
    return order_columns(results)


def expression_de(expression, metadata, dataset):
    control_samples = metadata.loc[
        metadata["condition"] == "control", "sample_id"
    ].tolist()
    treated_samples = metadata.loc[
        metadata["condition"] == "treated", "sample_id"
    ].tolist()
    if len(control_samples) < 2 or len(treated_samples) < 2:
        raise ValueError(
            f"{dataset['name']} needs at least two control and two treated samples."
        )

    expr = expression.set_index("gene_id")
    control = expr[control_samples].apply(pd.to_numeric, errors="coerce")
    treated = expr[treated_samples].apply(pd.to_numeric, errors="coerce")

    if dataset.get("value_scale", "log2") == "linear":
        control_for_fc = np.log2(control + 1)
        treated_for_fc = np.log2(treated + 1)
    else:
        control_for_fc = control
        treated_for_fc = treated

    log2fc = treated_for_fc.mean(axis=1) - control_for_fc.mean(axis=1)
    base_mean = pd.concat([control, treated], axis=1).mean(axis=1)
    pvalues = []
    for gene_id in expr.index:
        treated_values = treated_for_fc.loc[gene_id].dropna()
        control_values = control_for_fc.loc[gene_id].dropna()
        if len(treated_values) < 2 or len(control_values) < 2:
            pvalues.append(np.nan)
            continue
        _, pvalue = stats.ttest_ind(
            treated_values, control_values, equal_var=False, nan_policy="omit")
        pvalues.append(pvalue if math.isfinite(pvalue) else np.nan)

    results = pd.DataFrame({
        "gene_id": expr.index.astype(str),
        "gene": expression.set_index("gene_id").get(
            "gene", pd.Series(expr.index.astype(str), index=expr.index)
        ),
        "baseMean": base_mean,
        "log2FoldChange": log2fc,
        "FoldChange": np.power(2, log2fc),
        "pvalue": pvalues,
    }).reset_index(drop=True)

    annotation_columns = dataset.get("annotation_columns", [])
    sample_columns = set(control_samples + treated_samples)
    annotation_columns += [
        column for column in expression.columns
        if column not in sample_columns
        and column not in ["gene_id", "gene"]
        and column not in annotation_columns
    ]
    for column in annotation_columns:
        if column in expression.columns:
            results[column] = expression.set_index("gene_id")[column].reindex(
                results["gene_id"]
            ).to_numpy()
    return results


def load_tar_processed_text(dataset):
    archive = repo_path(dataset["archive"])
    if not archive.exists():
        raise FileNotFoundError(f"Missing {archive}")

    frames = []
    metadata_rows = []
    annotation = None
    sample_regex = re.compile(dataset["sample_name_regex"])

    with tarfile.open(archive) as tar:
        members = [
            member for member in tar.getmembers()
            if member.isfile()
            and dataset["processed_member_pattern"] in member.name
        ]
        for member in sorted(members, key=lambda item: item.name):
            with tar.extractfile(member) as raw:
                handle = gzip.open(raw, "rt", errors="replace")
                df = pd.read_csv(handle, sep="\t")

            value_column = df.columns[-1]
            sample_match = sample_regex.search(value_column)
            if not sample_match:
                sample_match = sample_regex.search(member.name)
            if not sample_match:
                continue

            sample_group = sample_match.group("sample_group")
            replicate = sample_match.groupdict().get("replicate", "")
            if sample_group in dataset["control_groups"]:
                condition = "control"
            elif sample_group in dataset["treated_groups"]:
                condition = "treated"
            else:
                continue

            sample_id = value_column
            gene_id_column = dataset["gene_id_column"]
            gene_column = dataset.get("gene_column", gene_id_column)
            keep_columns = list(dict.fromkeys(
                [gene_id_column, gene_column] + dataset.get("annotation_columns", [])
            ))
            if annotation is None:
                annotation = df[
                    [column for column in keep_columns if column in df.columns]
                ].rename(columns={gene_id_column: "gene_id"})
                if gene_column == gene_id_column:
                    annotation["gene"] = annotation["gene_id"]
                elif gene_column in annotation.columns:
                    annotation = annotation.rename(columns={gene_column: "gene"})
            frames.append(
                df[[gene_id_column, value_column]]
                .rename(columns={gene_id_column: "gene_id", value_column: sample_id})
                .set_index("gene_id")
            )
            metadata_rows.append({
                "sample_id": sample_id,
                "dataset": dataset["name"],
                "antibiotic_class": dataset["antibiotic_class"],
                "treatment": dataset["treatment"],
                "sample_group": sample_group,
                "condition": condition,
                "replicate": replicate,
                "source_file": member.name,
            })

    if not frames:
        raise ValueError(f"No configured processed samples found in {archive}")

    expression = pd.concat(frames, axis=1).reset_index()
    if annotation is not None:
        expression = annotation.drop_duplicates("gene_id").merge(
            expression, on="gene_id", how="right")
    metadata = pd.DataFrame(metadata_rows)
    expression = expression[["gene_id"] + [
        column for column in expression.columns if column != "gene_id"
    ]]
    return expression, metadata


def parse_series_matrix(path):
    sample_accessions = []
    sample_titles = []
    table_lines = []
    in_table = False
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", errors="replace") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith("!Sample_geo_accession"):
                sample_accessions = [
                    item.strip('"') for item in line.split("\t")[1:]
                ]
            elif line.startswith("!Sample_title"):
                sample_titles = [item.strip('"') for item in line.split("\t")[1:]]
            elif line == "!series_matrix_table_begin":
                in_table = True
            elif line == "!series_matrix_table_end":
                break
            elif in_table:
                table_lines.append(line)

    if not table_lines:
        raise ValueError(f"No series matrix table found in {path}")

    expression = pd.read_csv(io.StringIO("\n".join(table_lines)), sep="\t")
    expression.columns = [
        column.strip('"') for column in expression.columns
    ]
    if "ID_REF" not in expression.columns:
        raise ValueError(f"Series matrix {path} does not contain ID_REF.")
    expression = expression.rename(columns={"ID_REF": "gene_id"})
    return expression, sample_accessions, sample_titles


def load_platform_annotation(path):
    if not path.exists():
        return pd.DataFrame()

    table_lines = []
    in_table = False
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", errors="replace") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line == "!platform_table_begin":
                in_table = True
            elif line == "!platform_table_end":
                break
            elif in_table:
                table_lines.append(line)

    if not table_lines:
        return pd.DataFrame()

    annotation = pd.read_csv(io.StringIO("\n".join(table_lines)), sep="\t")
    annotation = annotation.rename(columns={
        "ID": "gene_id",
        "Gene symbol": "gene",
        "Gene title": "gene_title",
        "Platform_ORF": "locus_tag",
    })
    keep = [
        column for column in ["gene_id", "gene", "locus_tag", "gene_title"]
        if column in annotation.columns
    ]
    annotation = annotation[keep].copy()
    for column in ["gene", "locus_tag", "gene_title"]:
        if column in annotation.columns:
            annotation[column] = (
                annotation[column]
                .astype(str)
                .str.split("///")
                .str[0]
                .replace({"nan": ""})
            )
    return annotation.drop_duplicates("gene_id")


def load_series_matrix(dataset):
    path = repo_path(dataset["series_matrix"])
    raw_archive = repo_path(dataset.get("raw_archive", ""))
    if not path.exists():
        message = (
            f"Missing {path}. The raw gentamicin archive contains Affymetrix CEL "
            "files, which need external Affymetrix/CDF normalization before this "
            "Python pipeline can do gene-level differential expression. Download "
            "the GEO series matrix or another normalized expression matrix to the "
            f"configured path, or change {REPO_ROOT / 'config' / 'datasets.json'}."
        )
        if raw_archive.exists():
            message += f" Raw archive found at {raw_archive}."
        raise FileNotFoundError(message)

    expression, accessions, titles = parse_series_matrix(path)
    sample_labels = dataset.get("sample_labels", {})
    rename_map = {}
    metadata_rows = []
    for accession, title in zip(accessions, titles or accessions):
        if accession in dataset["control_samples"] or accession in dataset["treated_samples"]:
            label = sample_labels.get(accession, title or accession)
            rename_map[accession] = accession
            metadata_rows.append({
                "sample_id": accession,
                "dataset": dataset["name"],
                "antibiotic_class": dataset["antibiotic_class"],
                "treatment": dataset["treatment"],
                "sample_group": label,
                "condition": (
                    "control" if accession in dataset["control_samples"]
                    else "treated"
                ),
                "replicate": re.sub(r".*?([0-9]+)$", r"\1", label),
                "source_file": path.name,
            })

    selected_columns = ["gene_id"] + list(rename_map.values())
    missing = [column for column in selected_columns if column not in expression.columns]
    if missing:
        raise ValueError(f"Missing sample columns in {path}: {missing}")
    expression = expression[selected_columns]
    annotation_path = dataset.get("platform_annotation")
    if annotation_path:
        annotation = load_platform_annotation(repo_path(annotation_path))
        if not annotation.empty:
            expression = annotation.merge(expression, on="gene_id", how="right")
            expression["gene"] = expression["gene"].replace("", np.nan).fillna(
                expression["gene_id"]
            )
    return expression, pd.DataFrame(metadata_rows)


def load_excel_de_results(dataset):
    workbook = repo_path(dataset["workbook"])
    if not workbook.exists():
        raise FileNotFoundError(f"Missing {workbook}")

    df = pd.read_excel(workbook, sheet_name=dataset["de_sheet"])
    gene_id_column = dataset.get("gene_id_column", "gene_id")
    gene_column = dataset.get("gene_column", "gene")
    results = pd.DataFrame({
        "comparison": dataset.get("comparison", f"{dataset['treatment']}_vs_control"),
        "gene": df.get(gene_column, df.get("Name", df.get(gene_id_column, df.index))),
        "gene_id": df.get(gene_id_column, df.get("index", df.index)),
        "baseMean": pd.to_numeric(df.get("baseMean"), errors="coerce"),
        "log2FoldChange": pd.to_numeric(df.get("log2FoldChange"), errors="coerce"),
        "FoldChange": pd.to_numeric(df.get("FoldChange"), errors="coerce"),
        "pvalue": pd.to_numeric(df.get("pvalue"), errors="coerce"),
        "padj": pd.to_numeric(df.get("padj"), errors="coerce"),
    })

    counts = pd.DataFrame()
    metadata = pd.DataFrame()
    if dataset.get("counts_sheet"):
        counts = pd.read_excel(workbook, sheet_name=dataset["counts_sheet"])
        counts = counts.rename(columns={counts.columns[0]: "gene_id"})
        metadata_rows = []
        for sample_id in counts.columns:
            if sample_id == "gene_id":
                continue
            treatment = dataset["treatment"] if dataset["treatment"].lower() in sample_id.lower() else "control"
            metadata_rows.append({
                "sample_id": sample_id,
                "dataset": dataset["name"],
                "antibiotic_class": dataset["antibiotic_class"],
                "treatment": dataset["treatment"],
                "sample_group": treatment,
                "condition": "treated" if treatment == dataset["treatment"] else "control",
                "replicate": "",
                "source_file": workbook.name,
            })
        metadata = pd.DataFrame(metadata_rows)
    return results, counts, metadata


def write_standardized(dataset, expression=None, metadata=None, results=None):
    out_dir = REPO_ROOT / "data" / dataset["name"] / "standardized"
    out_dir.mkdir(parents=True, exist_ok=True)
    sheets = {}
    if expression is not None and not expression.empty:
        expression.to_csv(out_dir / "expression.csv", index=False)
        sheets["expression"] = expression
    if metadata is not None and not metadata.empty:
        metadata.to_csv(out_dir / "metadata.csv", index=False)
        sheets["metadata"] = metadata
    if results is not None and not results.empty:
        results.to_csv(out_dir / "de_results.csv", index=False)
        sheets["de_results"] = results
    if sheets:
        with pd.ExcelWriter(out_dir / "standardized_inputs.xlsx", engine="openpyxl") as writer:
            for name, df in sheets.items():
                df.to_excel(writer, sheet_name=name[:31], index=False)


def write_outputs(dataset, results, make_plots=True):
    final_dir = OUTPUT_ROOT / dataset["name"] / "final"
    plot_dir = OUTPUT_ROOT / dataset["name"] / "plots"
    final_dir.mkdir(parents=True, exist_ok=True)
    plot_dir.mkdir(parents=True, exist_ok=True)

    summary = readable(results)
    summary.to_csv(final_dir / "promoter_summary.csv", index=False)
    split_tables = {}
    for regulation in ["upregulated", "not_regulated", "downregulated"]:
        split = summary[summary["regulation"] == regulation].copy()
        split_tables[regulation] = split
        split.to_csv(final_dir / f"{regulation}_promoters.csv", index=False)

    with pd.ExcelWriter(final_dir / "promoter_summary.xlsx", engine="openpyxl") as writer:
        summary.to_excel(writer, sheet_name="all_promoters", index=False)
        for sheet_name, split in split_tables.items():
            split.to_excel(writer, sheet_name=sheet_name, index=False)

    if make_plots:
        write_volcano(dataset, results, plot_dir)


def write_volcano(dataset, results, plot_dir):
    df = results.dropna(subset=["log2FoldChange", "padj"]).copy()
    df["neg_log10_padj"] = -np.log10(
        pd.to_numeric(df["padj"], errors="coerce")
        .fillna(1.0)
        .clip(lower=np.finfo(float).tiny)
    )
    color_map = {
        "upregulated": "#d62728",
        "downregulated": "#1f77b4",
        "not_regulated": "#808080",
    }
    colors = df["regulation"].map(color_map).fillna("#808080")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(9, 7))
    ax.scatter(df["log2FoldChange"], df["neg_log10_padj"], c=colors, alpha=0.55, s=12)
    ax.axhline(-np.log10(0.05), color="black", linestyle="--", linewidth=0.8)
    ax.axvline(2, color="black", linestyle="--", linewidth=0.8)
    ax.axvline(-2, color="black", linestyle="--", linewidth=0.8)
    ax.set_xlabel("log2 Fold Change")
    ax.set_ylabel("-log10 adjusted p-value")
    ax.set_title(f"{dataset['treatment']} vs control")
    fig.tight_layout()
    fig.savefig(plot_dir / f"volcano_{dataset['name']}.png", dpi=300)
    plt.close(fig)

    try:
        import plotly.express as px
        hover_cols = [
            column for column in ["gene", "gene_id", "regulation", "padj", "signal_strength"]
            if column in df.columns
        ]
        fig_html = px.scatter(
            df,
            x="log2FoldChange",
            y="neg_log10_padj",
            color="regulation",
            color_discrete_map=color_map,
            hover_data=hover_cols,
            title=f"{dataset['treatment']} vs control",
        )
        fig_html.write_html(plot_dir / f"volcano_{dataset['name']}.html")
    except Exception as exc:
        print(f"Could not write interactive plot for {dataset['name']}: {exc}")


def analyze_dataset(dataset, thresholds, make_plots=True):
    print(f"\n=== {dataset['name']} ===")
    input_type = dataset["input_type"]
    if input_type == "tar_processed_text":
        expression, metadata = load_tar_processed_text(dataset)
        raw_results = expression_de(expression, metadata, dataset)
        results = finalize_results(raw_results, dataset, thresholds)
        write_standardized(dataset, expression=expression, metadata=metadata, results=results)
    elif input_type == "series_matrix":
        expression, metadata = load_series_matrix(dataset)
        raw_results = expression_de(expression, metadata, dataset)
        results = finalize_results(raw_results, dataset, thresholds)
        write_standardized(dataset, expression=expression, metadata=metadata, results=results)
    elif input_type == "excel_de_results":
        raw_results, counts, metadata = load_excel_de_results(dataset)
        results = finalize_results(raw_results, dataset, thresholds)
        write_standardized(dataset, expression=counts, metadata=metadata, results=results)
    else:
        raise ValueError(f"Unknown input_type: {input_type}")

    write_outputs(dataset, results, make_plots=make_plots)
    print(f"Rows: {len(results)}")
    print(f"Upregulated: {(results['regulation'] == 'upregulated').sum()}")
    print(f"Downregulated: {(results['regulation'] == 'downregulated').sum()}")
    print(f"Output: {OUTPUT_ROOT / dataset['name'] / 'final'}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--config",
        default="config/datasets.json",
        help="JSON dataset config to run.",
    )
    parser.add_argument(
        "--dataset",
        action="append",
        help="Run only this dataset name. Can be repeated.",
    )
    parser.add_argument(
        "--no-plots",
        action="store_true",
        help="Skip volcano plot generation.",
    )
    parser.add_argument(
        "--strict",
        action="store_true",
        help="Stop on the first dataset error.",
    )
    args = parser.parse_args()

    with open(repo_path(args.config)) as handle:
        config = json.load(handle)

    wanted = set(args.dataset or [])
    errors = []
    for dataset in config["datasets"]:
        if wanted and dataset["name"] not in wanted:
            continue
        try:
            analyze_dataset(
                dataset,
                config["thresholds"],
                make_plots=not args.no_plots,
            )
        except Exception as exc:
            message = f"{dataset['name']}: {exc}"
            errors.append(message)
            print(f"SKIPPED {message}")
            if args.strict:
                raise

    if errors:
        print("\nCompleted with skipped datasets:")
        for error in errors:
            print(f"- {error}")


if __name__ == "__main__":
    main()
