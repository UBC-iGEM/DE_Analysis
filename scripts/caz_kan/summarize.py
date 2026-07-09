import pandas as pd
import numpy as np
from pathlib import Path
import argparse
import os

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[1]
OUTPUT_DIR = REPO_ROOT / 'outputs' / 'caz_kan'
INTERMEDIATE_DIR = OUTPUT_DIR / 'intermediate'
FINAL_DIR = OUTPUT_DIR / 'final'
CEFTAZIDIME_DIR = FINAL_DIR / 'ceftazidime'
KANAMYCIN_DIR = FINAL_DIR / 'kanamycin'
SHARED_DIR = FINAL_DIR / 'shared'
PLOT_DIR = OUTPUT_DIR / 'plots'
MPL_CACHE_DIR = REPO_ROOT / 'outputs' / '_cache' / 'matplotlib'
CEFTAZIDIME_DIR.mkdir(parents=True, exist_ok=True)
KANAMYCIN_DIR.mkdir(parents=True, exist_ok=True)
SHARED_DIR.mkdir(parents=True, exist_ok=True)
PLOT_DIR.mkdir(parents=True, exist_ok=True)
MPL_CACHE_DIR.mkdir(parents=True, exist_ok=True)
os.environ.setdefault('MPLCONFIGDIR', str(MPL_CACHE_DIR))

parser = argparse.ArgumentParser()
parser.add_argument(
    '--no-plots',
    action='store_true',
    help='Write promoter CSVs without generating volcano PNGs.',
)
args = parser.parse_args()


def require_csv(filename):
    path = INTERMEDIATE_DIR / filename
    if not path.exists():
        raise FileNotFoundError(
            f"Missing {path}. Run deseq2_caz.py and deseq2_kan.py before summarize.py."
        )
    return path

results = pd.read_csv(require_csv('results_caz.csv'), index_col=0)
results_kan = pd.read_csv(require_csv('results_kan.csv'), index_col=0)
primary = pd.read_csv(require_csv('betalactam_primary.csv'), index_col=0)
primary_kan = pd.read_csv(require_csv('kanamycin_primary.csv'), index_col=0)


def add_signal_columns(df):
    df = df.dropna(subset=['gene', 'log2FoldChange']).copy()
    df['padj'] = pd.to_numeric(df['padj'], errors='coerce').fillna(1.0)
    safe_padj = df['padj'].clip(lower=np.finfo(float).tiny)
    df['neg_log10_padj'] = -np.log10(safe_padj)
    df['signal_strength'] = df['log2FoldChange'].abs() * df['neg_log10_padj']
    return df


def readable_output(df):
    rounded = df.copy()
    number_columns = rounded.select_dtypes(include=[np.number]).columns
    rounded[number_columns] = rounded[number_columns].round(2)
    return rounded


def promoter_lists(results, output_dir):
    df = add_signal_columns(results)
    upregulated = df[(df['log2FoldChange'] > 2) & (df['padj'] < 0.05)].copy()
    downregulated = df[(df['log2FoldChange'] < -2) & (df['padj'] < 0.05)].copy()
    regulated_genes = set(upregulated['gene']).union(downregulated['gene'])
    not_regulated = df[~df['gene'].isin(regulated_genes)].copy()

    outputs = {
        'upregulated': upregulated.sort_values(
            ['signal_strength', 'log2FoldChange'], ascending=[False, False]),
        'downregulated': downregulated.sort_values(
            ['signal_strength', 'log2FoldChange'], ascending=[False, True]),
        'not_regulated': not_regulated.sort_values(
            ['signal_strength', 'log2FoldChange'], ascending=[False, False]),
    }

    for regulation, data in outputs.items():
        readable_output(data).to_csv(
            output_dir / f'{regulation}_promoters.csv')

    return outputs


def shared_upregulated_promoters(beta_lactam_up, aminoglycoside_up):
    beta = beta_lactam_up[['gene', 'log2FoldChange', 'padj', 'signal_strength']].rename(
        columns={
            'log2FoldChange': 'beta_lactam_log2FoldChange',
            'padj': 'beta_lactam_padj',
            'signal_strength': 'beta_lactam_signal_strength',
        })
    amino = aminoglycoside_up[['gene', 'log2FoldChange', 'padj', 'signal_strength']].rename(
        columns={
            'log2FoldChange': 'aminoglycoside_log2FoldChange',
            'padj': 'aminoglycoside_padj',
            'signal_strength': 'aminoglycoside_signal_strength',
        })
    both = beta.merge(amino, on='gene')
    both['shared_signal_strength'] = both[
        ['beta_lactam_signal_strength', 'aminoglycoside_signal_strength']
    ].min(axis=1)
    both = both.sort_values(
        ['shared_signal_strength', 'beta_lactam_log2FoldChange',
         'aminoglycoside_log2FoldChange'],
        ascending=[False, False, False],
    )
    readable_output(both).to_csv(
        SHARED_DIR / 'both_classes_upregulated_promoters.csv', index=False)
    return both


def write_treatment_summary(outputs, antibiotic_class, treatment, output_path):
    rows = []
    for regulation, data in outputs.items():
        summary = data.copy()
        summary['antibiotic_class'] = antibiotic_class
        summary['treatment'] = treatment
        summary['regulation'] = regulation
        summary['shared_group'] = ''
        rows.append(summary)
    summary = pd.concat(rows, sort=False)
    leading_columns = ['antibiotic_class', 'treatment', 'regulation', 'shared_group', 'gene']
    ordered_columns = leading_columns + [col for col in summary.columns if col not in leading_columns]
    summary = summary[ordered_columns]
    readable_output(summary).to_csv(output_path, index=True, index_label='gene_id')
    return summary


def write_shared_summary(both_upregulated):
    if both_upregulated.empty:
        return
    both = both_upregulated.rename(
        columns={
            'beta_lactam_log2FoldChange': 'ceftazidime_log2FoldChange',
            'beta_lactam_padj': 'ceftazidime_padj',
            'aminoglycoside_log2FoldChange': 'kanamycin_log2FoldChange',
            'aminoglycoside_padj': 'kanamycin_padj',
        }).copy()
    both['antibiotic_class'] = 'both'
    both['treatment'] = 'ceftazidime_and_kanamycin'
    both['regulation'] = 'upregulated'
    both['shared_group'] = 'both_classes_upregulated'
    leading_columns = ['antibiotic_class', 'treatment', 'regulation', 'shared_group', 'gene']
    ordered_columns = leading_columns + [col for col in both.columns if col not in leading_columns]
    readable_output(both[ordered_columns]).to_csv(
        SHARED_DIR / 'promoter_summary.csv', index=False)


beta_lactam_lists = promoter_lists(results, CEFTAZIDIME_DIR)
aminoglycoside_lists = promoter_lists(results_kan, KANAMYCIN_DIR)
both_upregulated = shared_upregulated_promoters(
    beta_lactam_lists['upregulated'], aminoglycoside_lists['upregulated'])
ceftazidime_summary = write_treatment_summary(
    beta_lactam_lists, 'beta_lactam', 'ceftazidime',
    CEFTAZIDIME_DIR / 'promoter_summary.csv')
kanamycin_summary = write_treatment_summary(
    aminoglycoside_lists, 'aminoglycoside', 'kanamycin',
    KANAMYCIN_DIR / 'promoter_summary.csv')
write_shared_summary(both_upregulated)

# Cross-reactivity
caz_genes = set(primary['gene'].dropna())
kan_genes = set(primary_kan['gene'].dropna())
overlap = caz_genes.intersection(kan_genes)
print(f"Genes upregulated by both: {len(overlap)}")
print(f"Ceftazidime results written to {CEFTAZIDIME_DIR}")
print(f"Kanamycin results written to {KANAMYCIN_DIR}")
print(f"Shared results written to {SHARED_DIR}")
print(f"Upregulated in both classes: {len(both_upregulated)}")

key_candidates = ['recA', 'recN', 'emrA', 'pspD', 'lipA', 'dgcZ']
print("\nKey candidate cross-reactivity:")
for gene in key_candidates:
    print(f"{gene}: {'CROSS-REACTIVE' if gene in kan_genes else 'specific to ceftazidime'}")

kan_specific = primary_kan[~primary_kan['gene'].isin(caz_genes)].sort_values('log2FoldChange', ascending=False)
print(f"\nKanamycin specific candidates: {len(kan_specific)}")
kan_specific.to_csv(INTERMEDIATE_DIR / 'kanamycin_specific.csv')

def volcano_plot(results, title, output_path, highlight_genes=None):
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(10, 8))
    df = results.dropna(subset=['log2FoldChange', 'padj']).copy()
    df['neg_log10_padj'] = -np.log10(df['padj'])
    df['color'] = 'grey'
    df.loc[(df['log2FoldChange'] > 2) & (df['padj'] < 0.05), 'color'] = 'red'
    df.loc[(df['log2FoldChange'] < -2) & (df['padj'] < 0.05), 'color'] = 'blue'
    ax.scatter(df['log2FoldChange'], df['neg_log10_padj'], c=df['color'], alpha=0.5, s=10)
    ax.axhline(-np.log10(0.05), color='black', linestyle='--', linewidth=0.8)
    ax.axvline(2, color='black', linestyle='--', linewidth=0.8)
    ax.axvline(-2, color='black', linestyle='--', linewidth=0.8)
    if highlight_genes:
        for gene in highlight_genes:
            gene_data = df[df['gene'] == gene]
            if not gene_data.empty:
                ax.annotate(gene, (gene_data['log2FoldChange'].values[0],
                            gene_data['neg_log10_padj'].values[0]), fontsize=8, ha='right')
    ax.set_xlabel('log2 Fold Change')
    ax.set_ylabel('-log10 adjusted p-value')
    ax.set_title(title)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close(fig)
    print(f"Saved volcano plot to {output_path}")

if not args.no_plots:
    volcano_plot(results, 'Ceftazidime vs Control', PLOT_DIR / 'volcano_ceftazidime.png',
        highlight_genes=['recA', 'recN', 'dgcZ', 'lipA', 'emrA', 'pspD'])

    volcano_plot(results_kan, 'Kanamycin vs Control', PLOT_DIR / 'volcano_kanamycin.png',
        highlight_genes=['ibpB', 'pdeL', 'rybB', 'ompG'])
