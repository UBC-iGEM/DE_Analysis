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
  ceftazidime/
    GSE220559_RAW.tar
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
  ceftazidime/
  gentamicin/
  tobramycin/
```

`outputs/` is ignored by git and can be deleted whenever you want to regenerate
the analysis from scratch.

## Run Everything

Install dependencies:

```bash
python3 -m pip install -r requirements.txt
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
                        # tar_gene_tables, expression_matrix, read_counts_csv
count_matrix            # read-count CSV for read_counts_csv datasets
expression_matrix       # FPKM matrix for fpkm_matrix datasets
archive                 # tar archive for tar_gene_tables datasets
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

The current analysis includes four antibiotic datasets:

```text
amoxicillin   beta-lactam       GSE108190
ceftazidime   beta-lactam       GSE220559
gentamicin    aminoglycoside    GSE44211
tobramycin    aminoglycoside    GSE224240
```

Amoxicillin uses a processed FPKM matrix from fluoxetine-induced
amoxicillin-resistant mutants, so it should be interpreted as a resistance-state
comparison rather than a direct acute amoxicillin exposure. Ceftazidime uses raw
per-sample read-count tables from a tar archive. Gentamicin uses replicate
microarray expression values. Tobramycin already includes a processed
differential expression result sheet, so the pipeline standardizes and summarizes
that existing result.

For datasets without an input adjusted p-value, the pipeline uses Welch t-tests
across replicate expression/count values and Benjamini-Hochberg adjusted
p-values.
