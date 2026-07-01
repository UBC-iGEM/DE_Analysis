from pathlib import Path
import tarfile

import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

BASE_DIR = Path(__file__).resolve().parent
RAW_DIR = BASE_DIR / 'GSE220559_RAW'
RAW_TAR = BASE_DIR / 'GSE220559_RAW.tar'
MAPPING_FILE = BASE_DIR / 'gene_mapping.csv'

def discover_count_files(prefix):
    for pattern in (f'*{prefix}*.txt.gz', f'*{prefix}*.txt'):
        files = sorted(RAW_DIR.glob(pattern))
        if files:
            return files
    return []


def ensure_raw_data():
    if discover_count_files('CAZ') and discover_count_files('Wu'):
        return

    if not RAW_TAR.exists():
        raise FileNotFoundError(
            f'Could not find raw count files in {RAW_DIR} or archive {RAW_TAR}.'
        )

    RAW_DIR.mkdir(parents=True, exist_ok=True)
    with tarfile.open(RAW_TAR, 'r:*') as tar:
        tar.extractall(RAW_DIR, filter='data')


ensure_raw_data()
mapping = pd.read_csv(MAPPING_FILE, index_col=0)['gene']

caz_files = discover_count_files('CAZ')
wu_files = discover_count_files('Wu')

def load_counts(files):
    if not files:
        raise FileNotFoundError(
            f'No count files matched in {RAW_DIR}. '
            'Expected files like "*CAZ*.txt.gz" and "*Wu*.txt.gz".'
        )

    dfs = []
    for f in files:
        df = pd.read_csv(f, sep=r'\s+', index_col=0)
        df = df.iloc[:, 0:1]
        sample_name = f.name
        if sample_name.endswith('.txt.gz'):
            sample_name = sample_name[:-7]
        elif sample_name.endswith('.txt'):
            sample_name = sample_name[:-4]
        df.columns = [sample_name]
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
results['gene'] = results.index.map(mapping)

primary = results[(results['log2FoldChange'] > 2) & (results['padj'] < 0.05)].sort_values('log2FoldChange', ascending=False)
secondary = results[(results['log2FoldChange'] > 1) & (results['log2FoldChange'] <= 2) & (results['padj'] < 0.05)].sort_values('log2FoldChange', ascending=False)

print(f"Primary candidates: {len(primary)}")
print(primary[['gene', 'log2FoldChange', 'padj']].head(20).to_string())
print(f"\nSecondary candidates: {len(secondary)}")
print(secondary[['gene', 'log2FoldChange', 'padj']].head(20).to_string())

results.to_csv(BASE_DIR / 'results_caz.csv')
primary.to_csv(BASE_DIR / 'betalactam_primary.csv')
secondary.to_csv(BASE_DIR / 'betalactam_secondary.csv')
