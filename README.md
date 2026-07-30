# Promoter Selection Differential Expression Analysis

## Goal

Identify E. coli promoters that respond to antibiotic exposure and rank them for
biosensor design.

## Current Datasets

```text
data/
  amoxicillin/
    GSE108190_antibiotics_resistant_mutants.txt.gz
    standardized/
  gentamicin/
    GSE44211_RAW.tar
    GSE44211_series_matrix.txt.gz
    GPL3154.annot.gz
    standardized/
  tobramycin/
    GSE224240_analysis.xlsx
    standardized/

config/
  datasets.json

scripts/
  run_analysis.py

outputs/
  amoxicillin/
  gentamicin/
  tobramycin/
```

`outputs/` is ignored by git and can be deleted whenever you want to regenerate
the analysis from scratch.

## Run Everything

Install dependencies:

```bash
python3 -m pip install pandas numpy scipy matplotlib plotly openpyxl
```

From the repository root:

```bash
python3 scripts/run_analysis.py
```

Run one dataset only:

```bash
python3 scripts/run_analysis.py --dataset gentamicin
```

Skip volcano plots:

```bash
python3 scripts/run_analysis.py --no-plots
```

## Change Antibiotic Classes Or Comparisons

Edit:

```text
config/datasets.json
```

The key parameters are:

```text
name                    # dataset/output folder name
antibiotic_class        # e.g. aminoglycoside, beta_lactam
treatment               # antibiotic name
input_type              # fpkm_matrix, series_matrix, excel_de_results
count_matrix            # read-count CSV for read_counts_csv datasets
expression_matrix       # FPKM matrix for fpkm_matrix datasets
control_groups/samples  # controls
treated_groups/samples  # antibiotic-treated samples
value_scale             # log2 or linear
```

This lets you add or swap datasets without writing a new script for every
antibiotic.

## Outputs

Each dataset gets the same final structure:

```text
outputs/<dataset>/final/
  promoter_summary.csv
  promoter_summary.xlsx
  upregulated_promoters.csv
  not_regulated_promoters.csv
  downregulated_promoters.csv

outputs/<dataset>/plots/
  volcano_<dataset>.png
  volcano_<dataset>.html
```

Each Excel workbook has four sheets:

```text
all_promoters
upregulated
not_regulated
downregulated
```

All promoter summaries are sorted by highest `signal_strength` first. Numeric
outputs are rounded to two decimal places.

## Output Columns

Important columns:

```text
gene
gene_id
log2FoldChange
pvalue
padj
signal_strength
regulation
```

`regulation` is assigned with:

```text
upregulated    log2FoldChange > 2 and padj < 0.05
downregulated  log2FoldChange < -2 and padj < 0.05
not_regulated  everything else
```

`signal_strength` is:

```text
abs(log2FoldChange) * -log10(padj)
```

This ranks promoters by both effect size and statistical confidence.

## Notes

Amoxicillin and gentamicin are microarray datasets, so the pipeline uses
replicate expression values with Welch t-tests and Benjamini-Hochberg adjusted
p-values. Tobramycin already includes a processed DE result sheet, so the
pipeline standardizes and summarizes that existing result.
