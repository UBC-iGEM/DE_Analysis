from pathlib import Path
import re

import pandas as pd


SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[0]
DATA_ROOT = REPO_ROOT / "data"


STANDARD_COLUMNS = [
    "test_group",
    "comparison",
    "gene",
    "gene_id",
    "baseMean",
    "log2FoldChange",
    "FoldChange",
    "pvalue",
    "padj",
]


def write_standardized_workbook(output_dir, sheets):
    output_dir.mkdir(parents=True, exist_ok=True)
    for name, df in sheets.items():
        if not df.empty:
            df.to_csv(output_dir / f"{name}.csv", index=False)

    workbook_path = output_dir / "standardized_inputs.xlsx"
    with pd.ExcelWriter(workbook_path, engine="openpyxl") as writer:
        for name, df in sheets.items():
            if not df.empty:
                df.to_excel(writer, sheet_name=name[:31], index=False)
    print(f"Wrote {workbook_path}")


def standard_de_results(df, test_group, comparison):
    result = pd.DataFrame(index=df.index)
    result["test_group"] = test_group
    result["comparison"] = comparison
    result["gene"] = df.get("gene", df.get("Name", df.get("locus_tag", df.index)))
    result["gene_id"] = df.get("locus_tag", df.get("index", df.index))

    for column in ["baseMean", "log2FoldChange", "FoldChange", "pvalue", "padj"]:
        result[column] = pd.to_numeric(df.get(column), errors="coerce")

    extra_columns = [
        column for column in df.columns
        if column not in result.columns and column not in ["gene", "locus_tag"]
    ]
    for column in extra_columns:
        result[column] = df[column]

    ordered = STANDARD_COLUMNS + [
        column for column in result.columns if column not in STANDARD_COLUMNS
    ]
    return result[ordered]


def sample_from_caz_kan_path(path):
    match = re.search(r"_(CAZ|KAN|Wu)_([0-9]+)", path.stem)
    if not match:
        return None
    treatment_code, replicate = match.groups()
    treatment = {
        "CAZ": "ceftazidime",
        "KAN": "kanamycin",
        "Wu": "control",
    }[treatment_code]
    sample_id = f"{treatment}_{replicate}"
    return sample_id, treatment, int(replicate)


def standardize_caz_kan():
    raw_dir = DATA_ROOT / "caz_kan" / "GSE220559_RAW"
    output_dir = DATA_ROOT / "caz_kan" / "standardized"
    count_frames = []
    metadata_rows = []
    seen_samples = set()

    count_files = sorted(raw_dir.glob("*.txt"), key=lambda path: (" 2" in path.name, path.name))
    for path in count_files:
        parsed = sample_from_caz_kan_path(path)
        if parsed is None:
            continue

        sample_id, treatment, replicate = parsed
        if sample_id in seen_samples:
            print(f"Skipping duplicate CAZ/KAN sample {sample_id}: {path.name}")
            continue
        seen_samples.add(sample_id)

        df = pd.read_csv(path, sep=r"\s+")
        count_column = next(
            column for column in df.columns if column.endswith("_readcount")
        )
        sample_counts = df[["Gene_id", count_column]].rename(
            columns={"Gene_id": "gene_id", count_column: sample_id}
        )
        count_frames.append(sample_counts.set_index("gene_id"))
        metadata_rows.append({
            "sample_id": sample_id,
            "test_group": "caz_kan",
            "treatment": treatment,
            "condition": "control" if treatment == "control" else "treated",
            "replicate": replicate,
            "source_file": path.name,
        })

    if not count_frames:
        print(f"No CAZ/KAN raw count files found in {raw_dir}")
        return

    counts = pd.concat(count_frames, axis=1).reset_index()
    metadata = pd.DataFrame(metadata_rows).sort_values(
        ["treatment", "replicate"]
    )
    counts = counts[["gene_id"] + metadata["sample_id"].tolist()]
    write_standardized_workbook(
        output_dir,
        {
            "counts": counts,
            "metadata": metadata,
            "de_results": pd.DataFrame(),
        },
    )


def standardize_tobramycin():
    input_path = DATA_ROOT / "tobramycin" / "GSE224240_analysis.xlsx"
    output_dir = DATA_ROOT / "tobramycin" / "standardized"
    if not input_path.exists():
        print(f"No tobramycin workbook found at {input_path}")
        return

    counts = pd.read_excel(input_path, sheet_name="counts")
    counts = counts.rename(columns={counts.columns[0]: "gene_id"})
    metadata_rows = []
    for column in counts.columns:
        if column == "gene_id":
            continue
        strain = "NCM3416" if column.startswith("NCM") else "MG1655"
        treatment = "tobramycin" if "tob" in column.lower() else "control"
        metadata_rows.append({
            "sample_id": column,
            "test_group": "tobramycin",
            "strain": strain,
            "treatment": treatment,
            "condition": "control" if treatment == "control" else "treated",
            "replicate": len([
                row for row in metadata_rows
                if row["strain"] == strain and row["treatment"] == treatment
            ]) + 1,
            "source_file": input_path.name,
        })

    wt_results = pd.read_excel(input_path, sheet_name="MH vs tob")
    de_results = standard_de_results(
        wt_results, "tobramycin", "MG1655_tobramycin_vs_control"
    )
    write_standardized_workbook(
        output_dir,
        {
            "counts": counts,
            "metadata": pd.DataFrame(metadata_rows),
            "de_results": de_results,
        },
    )


def standardize_aminoglycoside():
    input_path = DATA_ROOT / "aminoglycoside" / "aminoglycoside_candidates.csv"
    output_dir = DATA_ROOT / "aminoglycoside" / "standardized"
    if not input_path.exists():
        print(f"No aminoglycoside candidate table found at {input_path}")
        return

    df = pd.read_csv(input_path)
    de_results = standard_de_results(
        df, "aminoglycoside", "imported_candidate_results"
    )
    metadata = pd.DataFrame([{
        "sample_id": "imported_candidate_results",
        "test_group": "aminoglycoside",
        "treatment": "aminoglycoside",
        "condition": "imported_de_results_only",
        "replicate": "",
        "source_file": input_path.name,
    }])
    write_standardized_workbook(
        output_dir,
        {
            "counts": pd.DataFrame(),
            "metadata": metadata,
            "de_results": de_results,
        },
    )


def main():
    standardize_aminoglycoside()
    standardize_caz_kan()
    standardize_tobramycin()


if __name__ == "__main__":
    main()
