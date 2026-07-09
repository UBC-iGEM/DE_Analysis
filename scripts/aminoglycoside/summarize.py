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
TOBRAMYCIN_OUTPUT_DIR = OUTPUT_ROOT / 'tobramycin'
GENTAMICIN_OUTPUT_DIR = OUTPUT_ROOT / 'gentamicin'
SHARED_OUTPUT_DIR = OUTPUT_ROOT / 'aminoglycoside_shared'
INTERMEDIATE_DIR = OUTPUT_ROOT / 'aminoglycoside' / 'intermediate'
TOBRAMYCIN_FINAL_DIR = TOBRAMYCIN_OUTPUT_DIR / 'final'
TOBRAMYCIN_PLOT_DIR = TOBRAMYCIN_OUTPUT_DIR / 'plots'
GENTAMICIN_FINAL_DIR = GENTAMICIN_OUTPUT_DIR / 'final'
SHARED_FINAL_DIR = SHARED_OUTPUT_DIR / 'final'
MPL_CACHE_DIR = OUTPUT_ROOT / '_cache' / 'matplotlib'
RAW_DIR = DATA_DIR / 'GSE228373_RAW'
TOBRAMYCIN_FINAL_DIR.mkdir(parents=True, exist_ok=True)
TOBRAMYCIN_PLOT_DIR.mkdir(parents=True, exist_ok=True)
GENTAMICIN_FINAL_DIR.mkdir(parents=True, exist_ok=True)
SHARED_FINAL_DIR.mkdir(parents=True, exist_ok=True)
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


def normalize_gentamicin():
    treated_path = RAW_DIR / 'GSM7119829_476822_S10.txt'
    untreated_path = RAW_DIR / 'GSM7119827_385274_DKP.txt'
    named_candidates_path = INTERMEDIATE_DIR / 'gentamicin_candidates_named.csv'
    candidates_path = INTERMEDIATE_DIR / 'gentamicin_candidates.csv'

    if treated_path.exists() and untreated_path.exists():
        treated = pd.read_csv(treated_path, sep='\t', names=['gene_id', 'treated'])
        untreated = pd.read_csv(untreated_path, sep='\t', names=['gene_id', 'untreated'])
        df = treated.merge(untreated, on='gene_id')
        df['treated'] = pd.to_numeric(df['treated'], errors='coerce')
        df['untreated'] = pd.to_numeric(df['untreated'], errors='coerce')
        df['log2FoldChange'] = np.log2((df['treated'] + 1) / (df['untreated'] + 1))
        df['gene'] = df['gene_id']

        if named_candidates_path.exists():
            names = pd.read_csv(named_candidates_path)[['gene_id', 'gene']].dropna()
            df = df.drop(columns=['gene']).merge(names, on='gene_id', how='left')
            df['gene'] = df['gene'].fillna(df['gene_id'])
    elif named_candidates_path.exists():
        df = pd.read_csv(named_candidates_path)
        df['log2FoldChange'] = pd.to_numeric(df['log2FC'], errors='coerce')
    elif candidates_path.exists():
        df = pd.read_csv(candidates_path)
        df['gene'] = df['gene_id']
        df['log2FoldChange'] = pd.to_numeric(df['log2FC'], errors='coerce')
    else:
        print('No gentamicin input found; skipping gentamicin summary')
        return pd.DataFrame()

    log2fc = pd.to_numeric(df['log2FoldChange'], errors='coerce')
    summary = pd.DataFrame({
        'antibiotic_class': 'aminoglycoside',
        'treatment': 'gentamicin',
        'regulation': [classify(fc, require_padj=False) for fc in log2fc],
        'shared_group': '',
        'gene': df.get('gene', df['gene_id']),
        'gene_id': df['gene_id'],
        'log2FoldChange': log2fc,
        'padj': 1.0,
        'padj_source': 'not_available_no_replicates_set_to_1',
        'signal_strength': [signal_strength(fc) for fc in log2fc],
    })
    return summary.dropna(subset=['log2FoldChange'])


def shared_aminoglycoside_promoters(tobramycin, gentamicin):
    if tobramycin.empty or gentamicin.empty:
        return pd.DataFrame()

    tob_up = tobramycin[tobramycin['regulation'] == 'upregulated']
    gen_up = gentamicin[gentamicin['regulation'] == 'upregulated']
    shared = tob_up.merge(
        gen_up,
        on='gene',
        suffixes=('_tobramycin', '_gentamicin'),
    )
    if shared.empty:
        return shared

    shared['antibiotic_class'] = 'aminoglycoside'
    shared['treatment'] = 'tobramycin_and_gentamicin'
    shared['regulation'] = 'upregulated'
    shared['shared_group'] = 'both_aminoglycosides_upregulated'
    shared['log2FoldChange'] = shared[
        ['log2FoldChange_tobramycin', 'log2FoldChange_gentamicin']
    ].min(axis=1)
    shared['signal_strength'] = shared[
        ['signal_strength_tobramycin', 'signal_strength_gentamicin']
    ].min(axis=1)
    shared['padj'] = shared['padj_tobramycin']
    shared['padj_source'] = shared['padj_source_tobramycin']
    shared['gene_id'] = shared['gene_id_tobramycin']
    return shared


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


def write_summary(df, output_path):
    if df.empty:
        return
    readable_output(df).to_csv(output_path, index=False)
    print(f"Readable summary written to {output_path}")
    print(f"Rows: {len(df)}")


tobramycin = normalize_tobramycin()
gentamicin = normalize_gentamicin()
shared = shared_aminoglycoside_promoters(tobramycin, gentamicin)

leading_columns = [
    'antibiotic_class', 'treatment', 'regulation', 'shared_group', 'gene',
    'gene_id', 'log2FoldChange', 'padj', 'signal_strength',
    'padj_source',
]

for df, output_path in [
        (tobramycin, TOBRAMYCIN_FINAL_DIR / 'promoter_summary.csv'),
        (gentamicin, GENTAMICIN_FINAL_DIR / 'promoter_summary.csv'),
        (shared, SHARED_FINAL_DIR / 'both_aminoglycosides_upregulated.csv')]:
    if not df.empty:
        df = df.sort_values(
            ['regulation', 'signal_strength'], ascending=[True, False])
        df = df[[col for col in leading_columns if col in df.columns] +
                [col for col in df.columns if col not in leading_columns]]
        write_summary(df, output_path)

if not args.no_plots and not tobramycin.empty:
    volcano_plot(tobramycin, 'Tobramycin vs Control',
                 TOBRAMYCIN_PLOT_DIR / 'volcano_tobramycin.png')
