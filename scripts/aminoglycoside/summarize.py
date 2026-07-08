import numpy as np
import pandas as pd
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[1]
DATA_DIR = REPO_ROOT / 'data' / 'aminoglycoside'
OUTPUT_DIR = REPO_ROOT / 'outputs' / 'aminoglycoside'
INTERMEDIATE_DIR = OUTPUT_DIR / 'intermediate'
FINAL_DIR = OUTPUT_DIR / 'final'
RAW_DIR = DATA_DIR / 'GSE228373_RAW'
SUMMARY_PATH = FINAL_DIR / 'promoter_summary.csv'
FINAL_DIR.mkdir(parents=True, exist_ok=True)


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
    analysis_path = DATA_DIR / 'GSE224240_analysis.xlsx'
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


tobramycin = normalize_tobramycin()
gentamicin = normalize_gentamicin()
shared = shared_aminoglycoside_promoters(tobramycin, gentamicin)

summary = pd.concat([tobramycin, gentamicin, shared], sort=False)
summary = summary.sort_values(
    ['shared_group', 'treatment', 'regulation', 'signal_strength'],
    ascending=[True, True, True, False],
)

leading_columns = [
    'antibiotic_class', 'treatment', 'regulation', 'shared_group', 'gene',
    'gene_id', 'log2FoldChange', 'padj', 'signal_strength',
    'padj_source',
]
summary = summary[[col for col in leading_columns if col in summary.columns] +
                  [col for col in summary.columns if col not in leading_columns]]
readable_output(summary).to_csv(SUMMARY_PATH, index=False)

print(f"Readable aminoglycoside summary written to {SUMMARY_PATH}")
print(f"Rows: {len(summary)}")
