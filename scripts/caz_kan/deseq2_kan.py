import pandas as pd
import numpy as np
from pathlib import Path
import tarfile
import os

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[1]
DATA_DIR = REPO_ROOT / 'data' / 'caz_kan'
OUTPUT_DIR = REPO_ROOT / 'outputs' / 'caz_kan'
INTERMEDIATE_DIR = OUTPUT_DIR / 'intermediate'
SCRATCH_DIR = OUTPUT_DIR / 'scratch'
RAW_DIR = SCRATCH_DIR / 'GSE220559_RAW'
RAW_ARCHIVE = DATA_DIR / 'GSE220559_RAW.tar'
EXTRACTED_RAW_DIR = DATA_DIR / 'GSE220559_RAW'
STANDARDIZED_DIR = DATA_DIR / 'standardized'
STANDARDIZED_COUNTS = STANDARDIZED_DIR / 'counts.csv'
STANDARDIZED_METADATA = STANDARDIZED_DIR / 'metadata.csv'
MAPPING_PATH = INTERMEDIATE_DIR / 'gene_mapping.csv'
MPL_CACHE_DIR = REPO_ROOT / 'outputs' / '_cache' / 'matplotlib'
INTERMEDIATE_DIR.mkdir(parents=True, exist_ok=True)
SCRATCH_DIR.mkdir(parents=True, exist_ok=True)
MPL_CACHE_DIR.mkdir(parents=True, exist_ok=True)
os.environ.setdefault('MPLCONFIGDIR', str(MPL_CACHE_DIR))

_sysconf = os.sysconf


def safe_sysconf(name):
    if name == 'SC_SEM_NSEMS_MAX':
        try:
            return _sysconf(name)
        except PermissionError:
            return 256
    return _sysconf(name)


os.sysconf = safe_sysconf

from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from pydeseq2.default_inference import DefaultInference


def ensure_raw_dir():
    if EXTRACTED_RAW_DIR.exists():
        return EXTRACTED_RAW_DIR
    if not RAW_DIR.exists() and RAW_ARCHIVE.exists():
        RAW_DIR.mkdir()
        with tarfile.open(RAW_ARCHIVE) as archive:
            archive.extractall(RAW_DIR)
    return RAW_DIR


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


def load_counts(files):
    dfs = []
    for f in files:
        df = pd.read_csv(f, sep=r'\s+', index_col=0)
        df = df.iloc[:, 0:1]
        df.columns = [f.stem]
        dfs.append(df)
    return pd.concat(dfs, axis=1)

def load_standardized_counts(treatment):
    count_table = pd.read_csv(STANDARDIZED_COUNTS)
    metadata = pd.read_csv(STANDARDIZED_METADATA)
    selected = metadata[metadata['treatment'].isin([treatment, 'control'])].copy()
    selected['condition'] = selected['treatment'].where(
        selected['treatment'] == treatment, 'control')
    sample_ids = selected['sample_id'].tolist()
    counts = count_table.set_index('gene_id')[sample_ids].T
    return counts, selected.set_index('sample_id')[['condition']]

if STANDARDIZED_COUNTS.exists() and STANDARDIZED_METADATA.exists():
    counts, metadata = load_standardized_counts('kanamycin')
else:
    raw_dir = ensure_raw_dir()
    kan_files = sorted(raw_dir.glob('*KAN*.txt*'))
    wu_files = sorted(raw_dir.glob('*Wu*.txt*'))
    if not kan_files or not wu_files:
        raise FileNotFoundError(
            f"Could not find KAN and Wu count files in {raw_dir}. "
            f"Keep {RAW_ARCHIVE.name} in {DATA_DIR}, extract it to {EXTRACTED_RAW_DIR}, "
            "or run scripts/standardize_inputs.py."
        )

    kan = load_counts(kan_files)
    wu = load_counts(wu_files)
    counts = pd.concat([kan, wu], axis=1).T
    metadata = pd.DataFrame({
        'condition': ['kanamycin']*len(kan.columns) + ['control']*len(wu.columns)
    }, index=counts.index)

inference = DefaultInference(n_cpus=1)
dds = DeseqDataSet(
    counts=counts,
    metadata=metadata,
    design_factors='condition',
    inference=inference,
    n_cpus=1,
)
dds.deseq2()
stats = DeseqStats(
    dds,
    contrast=['condition', 'kanamycin', 'control'],
    inference=inference,
    n_cpus=1,
)
stats.summary()
results_kan = stats.results_df
results_kan = ensure_padj(results_kan)
results_kan['gene'] = gene_names_for(results_kan.index)

primary_kan = results_kan[(results_kan['log2FoldChange'] > 2) & (results_kan['padj'] < 0.05)].sort_values('log2FoldChange', ascending=False)

print(f"Kanamycin primary candidates: {len(primary_kan)}")
print(primary_kan[['gene', 'log2FoldChange', 'padj']].head(20).to_string())

results_kan.to_csv(INTERMEDIATE_DIR / 'results_kan.csv')
primary_kan.to_csv(INTERMEDIATE_DIR / 'kanamycin_primary.csv')
