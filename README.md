# Promoter Selection Differential Expression Analysis

## Goal
Identify E. coli promoters activated by aminoglycoside and beta-lactam
antibiotics for use as biosensor components.

## Directory Structure

```text
data/
  aminoglycoside/
    aminoglycoside_candidates.csv   # tracked fallback input for tobramycin
    standardized/                   # regenerated standardized import tables
  tobramycin/
    GSE224240_analysis.xlsx         # optional processed tobramycin input
    standardized/
  caz_kan/
    GSE220559_RAW/                  # optional extracted raw count files
    GSE220559_RAW.tar               # optional raw count archive
    standardized/

scripts/
  standardize_inputs.py
  aminoglycoside/
    filter.py
    filter_tr_untr.py
    mapping.py
    compare.py
    summarize.py
  caz_kan/
    mapping.py
    deseq2_caz.py
    deseq2_kan.py
    summarize.py
    interactive_plot.py

outputs/
  tobramycin/                       # regenerated locally, git-ignored
    final/
    plots/
  gentamicin/                       # generated only when gentamicin inputs exist
    final/
  aminoglycoside_shared/
    final/
  aminoglycoside/
    intermediate/
  caz_kan/                          # regenerated locally, git-ignored
    final/
      ceftazidime/
      kanamycin/
      shared/
    intermediate/
    plots/
    scratch/
```

`outputs/` is intentionally ignored by git. Delete it any time you want to run
the analysis from scratch.

## Standardize Inputs

Run this first after adding or extracting dataset files:

```bash
python3 scripts/standardize_inputs.py
```

This creates the same import layout for each test group:

```text
data/aminoglycoside/standardized/
data/caz_kan/standardized/
data/tobramycin/standardized/
```

Each standardized folder contains CSV sheets, plus one Excel workbook:

```text
counts.csv                  # sample count matrix when raw counts are available
metadata.csv                # sample names, treatment labels, controls, replicates
de_results.csv              # imported DE result table when already available
standardized_inputs.xlsx    # the same sheets in one workbook
```

`caz_kan` is treated as one test group with ceftazidime, kanamycin, and water
control samples. `tobramycin` is treated as a separate test group. The older
`aminoglycoside` input currently contains imported candidate-level DE results,
not raw replicate counts, so it standardizes to metadata plus `de_results.csv`.

## Dependencies

```bash
python3 -m pip install pandas numpy pydeseq2 matplotlib plotly openpyxl
```

`matplotlib` and `plotly` are only needed for plots. Use `--no-plots` when you
only want CSV summaries.

## Run CAZ/KAN Analysis

From the repository root:

```bash
python3 scripts/standardize_inputs.py
python3 scripts/caz_kan/mapping.py
python3 scripts/caz_kan/deseq2_caz.py
python3 scripts/caz_kan/deseq2_kan.py
python3 scripts/caz_kan/summarize.py
```

Main readable output:

```text
outputs/caz_kan/final/ceftazidime/promoter_summary.csv
outputs/caz_kan/final/kanamycin/promoter_summary.csv
```

Additional regenerated outputs include DESeq result CSVs, per-category promoter
lists, and optional volcano PNG/HTML plots:

```text
outputs/caz_kan/intermediate/       # DESeq CSVs and mapping files
outputs/caz_kan/final/ceftazidime/
outputs/caz_kan/final/kanamycin/
outputs/caz_kan/final/shared/
outputs/caz_kan/plots/
outputs/caz_kan/scratch/            # extracted raw archive
```

The CAZ/KAN scripts prefer `data/caz_kan/standardized/counts.csv` and
`data/caz_kan/standardized/metadata.csv` when they exist. If those files are
missing, they fall back to extracted files in `data/caz_kan/GSE220559_RAW/`, and
then to the compressed archive in `data/caz_kan/GSE220559_RAW.tar`. If
`data/caz_kan/ecoli_k12.gtf.gz` is absent, `mapping.py` writes an ID-only mapping
so the pipeline can still run; add the GTF file for named genes.

## Run Aminoglycoside Summary

From the repository root:

```bash
python3 scripts/standardize_inputs.py
python3 scripts/aminoglycoside/summarize.py
```

Main readable output:

```text
outputs/tobramycin/final/promoter_summary.csv
```

The tobramycin volcano plot is saved to:

```text
outputs/aminoglycoside/plots/volcano_tobramycin.png
outputs/tobramycin/plots/volcano_tobramycin.png
```

The summary uses `data/tobramycin/GSE224240_analysis.xlsx` when present, and
falls back to `data/aminoglycoside/aminoglycoside_candidates.csv` otherwise. To
regenerate aminoglycoside intermediates from raw data, place these optional
inputs in `data/aminoglycoside/`:

```text
GSE224240_analysis.xlsx
GSE228373_RAW/
ecoli_annotation.gtf
```

Then run:

```bash
python3 scripts/aminoglycoside/filter.py
python3 scripts/aminoglycoside/filter_tr_untr.py
python3 scripts/aminoglycoside/mapping.py
python3 scripts/aminoglycoside/compare.py
python3 scripts/aminoglycoside/summarize.py
```

## Output Conventions

Promoter summaries classify rows as:

- `upregulated`: `log2FoldChange > 2` and, where available, `padj < 0.05`
- `downregulated`: `log2FoldChange < -2` and, where available, `padj < 0.05`
- `not_regulated`: all other rows

Readable summary CSVs round numeric values to two decimal places and sort by
signal strength.
