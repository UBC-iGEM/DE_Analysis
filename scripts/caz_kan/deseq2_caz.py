import pandas as pd
import numpy as np
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from pathlib import Path
import tarfile

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[1]
DATA_DIR = REPO_ROOT / 'data' / 'caz_kan'
OUTPUT_DIR = REPO_ROOT / 'outputs' / 'caz_kan'
INTERMEDIATE_DIR = OUTPUT_DIR / 'intermediate'
SCRATCH_DIR = OUTPUT_DIR / 'scratch'
RAW_DIR = SCRATCH_DIR / 'GSE220559_RAW'
RAW_ARCHIVE = DATA_DIR / 'GSE220559_RAW.tar'
MAPPING_PATH = INTERMEDIATE_DIR / 'gene_mapping.csv'
INTERMEDIATE_DIR.mkdir(parents=True, exist_ok=True)
SCRATCH_DIR.mkdir(parents=True, exist_ok=True)


def ensure_raw_dir():
    if not RAW_DIR.exists() and RAW_ARCHIVE.exists():
        RAW_DIR.mkdir()
        with tarfile.open(RAW_ARCHIVE) as archive:
            archive.extractall(RAW_DIR)


def gene_names_for(index):
    fallback = pd.Series(index.astype(str), index=index)
    if MAPPING_PATH.exists():
        mapping = pd.read_csv(MAPPING_PATH, index_col=0)['gene']
        return fallback.map(mapping).fillna(fallback)
    print(f"{MAPPING_PATH.name} not found; using gene IDs as promoter names")
    return fallback


def benjamini_hochberg(pvalues):
    pvalues = pd.Series(pd.to_numeric(pvalues, errors='coerce')).fillna(1.0)
    ranked = pvalues.rank(method='first').astype(int)
    adjusted = pvalues * len(pvalues) / ranked
    adjusted = adjusted.sort_values(ascending=False).cummin().sort_index()
    return adjusted.clip(upper=1.0)


def ensure_padj(results):
    results = results.copy()
    results['pvalue'] = pd.to_numeric(results['pvalue'], errors='coerce')
    results['padj'] = pd.to_numeric(results['padj'], errors='coerce')
    filled_padj = benjamini_hochberg(results['pvalue'])
    results['padj_source'] = np.where(results['padj'].notna(), 'DESeq2', 'BH_from_pvalue')
    results['padj'] = results['padj'].fillna(filled_padj).fillna(1.0)
    results.loc[results['pvalue'].isna(), 'padj_source'] = 'no_valid_pvalue_set_to_1'
    return results


ensure_raw_dir()
caz_files = sorted(RAW_DIR.glob('*CAZ*.txt*'))
wu_files = sorted(RAW_DIR.glob('*Wu*.txt*'))
if not caz_files or not wu_files:
    raise FileNotFoundError(
        f"Could not find CAZ and Wu count files in {RAW_DIR}. "
        f"Keep {RAW_ARCHIVE.name} in {DATA_DIR} or extract it to {RAW_DIR}."
    )

def load_counts(files):
    dfs = []
    for f in files:
        df = pd.read_csv(f, sep=r'\s+', index_col=0)
        df = df.iloc[:, 0:1]
        df.columns = [f.stem]
        dfs.append(df)
    return pd.concat(dfs, axis=1)

caz = load_counts(caz_files)
wu = load_counts(wu_files)
counts = pd.concat([caz, wu], axis=1).T

metadata = pd.DataFrame({
    'condition': ['ceftazidime']*3 + ['control']*3
}, index=counts.index)

dds = DeseqDataSet(counts=counts, metadata=metadata, design_factors='condition')
dds.deseq2()
stats = DeseqStats(dds, contrast=['condition', 'ceftazidime', 'control'])
stats.summary()
results = stats.results_df
results = ensure_padj(results)
results['gene'] = gene_names_for(results.index)

primary = results[(results['log2FoldChange'] > 2) & (results['padj'] < 0.05)].sort_values('log2FoldChange', ascending=False)
secondary = results[(results['log2FoldChange'] > 1) & (results['log2FoldChange'] <= 2) & (results['padj'] < 0.05)].sort_values('log2FoldChange', ascending=False)

print(f"Primary candidates: {len(primary)}")
print(primary[['gene', 'log2FoldChange', 'padj']].head(20).to_string())
print(f"\nSecondary candidates: {len(secondary)}")
print(secondary[['gene', 'log2FoldChange', 'padj']].head(20).to_string())

results.to_csv(INTERMEDIATE_DIR / 'results_caz.csv')
primary.to_csv(INTERMEDIATE_DIR / 'betalactam_primary.csv')
secondary.to_csv(INTERMEDIATE_DIR / 'betalactam_secondary.csv')
