import numpy as np
import pandas as pd
from pathlib import Path
import argparse
import os

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[1]
DATA_DIR = REPO_ROOT / 'data' / 'aminoglycoside'
TOBRAMYCIN_DIR = REPO_ROOT / 'data' / 'tobramycin'
OUTPUT_ROOT = REPO_ROOT / 'outputs'
AMINOGLYCOSIDE_OUTPUT_DIR = OUTPUT_ROOT / 'aminoglycoside'
TOBRAMYCIN_OUTPUT_DIR = OUTPUT_ROOT / 'tobramycin'
AMINOGLYCOSIDE_FINAL_DIR = AMINOGLYCOSIDE_OUTPUT_DIR / 'final'
TOBRAMYCIN_FINAL_DIR = TOBRAMYCIN_OUTPUT_DIR / 'final'
TOBRAMYCIN_PLOT_DIR = TOBRAMYCIN_OUTPUT_DIR / 'plots'
MPL_CACHE_DIR = OUTPUT_ROOT / '_cache' / 'matplotlib'
AMINOGLYCOSIDE_FINAL_DIR.mkdir(parents=True, exist_ok=True)
TOBRAMYCIN_FINAL_DIR.mkdir(parents=True, exist_ok=True)
TOBRAMYCIN_PLOT_DIR.mkdir(parents=True, exist_ok=True)
MPL_CACHE_DIR.mkdir(parents=True, exist_ok=True)
os.environ.setdefault('MPLCONFIGDIR', str(MPL_CACHE_DIR))

parser = argparse.ArgumentParser()
parser.add_argument(
    '--no-plots',
    action='store_true',
    help='Write promoter CSVs without generating volcano PNGs.',
)
args = parser.parse_args()


def readable_output(df):
    rounded = df.copy()
    number_columns = rounded.select_dtypes(include=[np.number]).columns
    rounded[number_columns] = rounded[number_columns].round(2)
    return rounded


def classify(log2fc, padj=None, require_padj=True):
    significant = True if not require_padj else padj < 0.05
    if significant and log2fc > 2:
        return 'upregulated'
    if significant and log2fc < -2:
        return 'downregulated'
    return 'not_regulated'


def signal_strength(log2fc, padj=None):
    if padj is None or pd.isna(padj):
        return abs(log2fc)
    safe_padj = max(float(padj), np.finfo(float).tiny)
    return abs(log2fc) * -np.log10(safe_padj)


def normalize_imported_aminoglycoside():
    standardized_path = DATA_DIR / 'standardized' / 'de_results.csv'
    candidates_path = DATA_DIR / 'aminoglycoside_candidates.csv'

    if standardized_path.exists():
        df = pd.read_csv(standardized_path)
    elif candidates_path.exists():
        df = pd.read_csv(candidates_path)
    else:
        print('No imported aminoglycoside input found; skipping aminoglycoside summary')
        return pd.DataFrame()

    gene_id = df.get('gene_id', df.get('locus_tag', df.get('index', pd.Series(df.index, index=df.index))))
    gene = df.get('gene', df.get('Name', gene_id))
    log2fc = pd.to_numeric(df['log2FoldChange'], errors='coerce')
    padj = pd.to_numeric(df.get('padj'), errors='coerce').fillna(1.0)

    summary = pd.DataFrame({
        'antibiotic_class': 'aminoglycoside',
        'treatment': 'imported_aminoglycoside',
        'regulation': [
            classify(fc, adj) for fc, adj in zip(log2fc, padj)
        ],
        'shared_group': '',
        'gene': gene,
        'gene_id': gene_id,
        'log2FoldChange': log2fc,
        'padj': padj,
        'padj_source': 'input_or_set_to_1_if_missing',
        'signal_strength': [
            signal_strength(fc, adj) for fc, adj in zip(log2fc, padj)
        ],
    })
    return summary.dropna(subset=['log2FoldChange'])


def normalize_tobramycin():
    analysis_path = TOBRAMYCIN_DIR / 'GSE224240_analysis.xlsx'
    candidates_path = DATA_DIR / 'aminoglycoside_candidates.csv'

    if analysis_path.exists():
        df = pd.read_excel(analysis_path)
    elif candidates_path.exists():
        df = pd.read_csv(candidates_path)
    else:
        print('No tobramycin input found; skipping tobramycin summary')
        return pd.DataFrame()

    gene_id = df.get('locus_tag', df.get('index', pd.Series(df.index, index=df.index)))
    gene = df.get('gene', df.get('Name', gene_id))
    log2fc = pd.to_numeric(df['log2FoldChange'], errors='coerce')
    padj = pd.to_numeric(df.get('padj'), errors='coerce').fillna(1.0)

    summary = pd.DataFrame({
        'antibiotic_class': 'aminoglycoside',
        'treatment': 'tobramycin',
        'regulation': [
            classify(fc, adj) for fc, adj in zip(log2fc, padj)
        ],
        'shared_group': '',
        'gene': gene,
        'gene_id': gene_id,
        'log2FoldChange': log2fc,
        'padj': padj,
        'padj_source': 'input_or_set_to_1_if_missing',
        'signal_strength': [
            signal_strength(fc, adj) for fc, adj in zip(log2fc, padj)
        ],
    })
    return summary.dropna(subset=['log2FoldChange'])


def volcano_plot(results, title, output_path):
    import matplotlib.pyplot as plt

    df = results.dropna(subset=['log2FoldChange', 'padj']).copy()
    safe_padj = pd.to_numeric(df['padj'], errors='coerce').fillna(1.0)
    safe_padj = safe_padj.clip(lower=np.finfo(float).tiny)
    df['neg_log10_padj'] = -np.log10(safe_padj)
    df['color'] = 'grey'
    df.loc[(df['log2FoldChange'] > 2) & (df['padj'] < 0.05), 'color'] = 'red'
    df.loc[(df['log2FoldChange'] < -2) & (df['padj'] < 0.05), 'color'] = 'blue'

    fig, ax = plt.subplots(figsize=(10, 8))
    ax.scatter(df['log2FoldChange'], df['neg_log10_padj'], c=df['color'], alpha=0.5, s=10)
    ax.axhline(-np.log10(0.05), color='black', linestyle='--', linewidth=0.8)
    ax.axvline(2, color='black', linestyle='--', linewidth=0.8)
    ax.axvline(-2, color='black', linestyle='--', linewidth=0.8)
    ax.set_xlabel('log2 Fold Change')
    ax.set_ylabel('-log10 adjusted p-value')
    ax.set_title(title)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close(fig)
    print(f"Saved volcano plot to {output_path}")


def ordered_summary(df):
    if df.empty:
        return df
    df = df.sort_values(
        ['signal_strength', 'log2FoldChange'], ascending=[False, False])
    leading_columns = [
        'antibiotic_class', 'treatment', 'regulation', 'shared_group', 'gene',
        'gene_id', 'log2FoldChange', 'padj', 'signal_strength',
        'padj_source',
    ]
    return df[[col for col in leading_columns if col in df.columns] +
              [col for col in df.columns if col not in leading_columns]]


def write_dataset_outputs(df, output_dir):
    if df.empty:
        return
    summary = readable_output(ordered_summary(df))
    summary.to_csv(output_dir / 'promoter_summary.csv', index=False)

    split_tables = {}
    for regulation in ['upregulated', 'not_regulated', 'downregulated']:
        split = summary[summary['regulation'] == regulation].copy()
        split_tables[regulation] = split
        split.to_csv(output_dir / f'{regulation}_promoters.csv', index=False)

    with pd.ExcelWriter(output_dir / 'promoter_summary.xlsx', engine='openpyxl') as writer:
        summary.to_excel(writer, sheet_name='all_promoters', index=False)
        for sheet_name, split in split_tables.items():
            split.to_excel(writer, sheet_name=sheet_name, index=False)

    print(f"Readable summary written to {output_dir}")
    print(f"Rows: {len(summary)}")


imported_aminoglycoside = normalize_imported_aminoglycoside()
tobramycin = normalize_tobramycin()

write_dataset_outputs(imported_aminoglycoside, AMINOGLYCOSIDE_FINAL_DIR)
write_dataset_outputs(tobramycin, TOBRAMYCIN_FINAL_DIR)

if not args.no_plots and not tobramycin.empty:
    volcano_plot(tobramycin, 'Tobramycin vs Control',
                 TOBRAMYCIN_PLOT_DIR / 'volcano_tobramycin.png')
