import pandas as pd
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[1]
DATA_DIR = REPO_ROOT / 'data' / 'aminoglycoside'
TOBRAMYCIN_DIR = REPO_ROOT / 'data' / 'tobramycin'
OUTPUT_DIR = REPO_ROOT / 'outputs' / 'aminoglycoside'
INTERMEDIATE_DIR = OUTPUT_DIR / 'intermediate'
INTERMEDIATE_DIR.mkdir(parents=True, exist_ok=True)

# Load data
df = pd.read_excel(TOBRAMYCIN_DIR / 'GSE224240_analysis.xlsx')

# Filter for significantly upregulated genes
sig = df[(df['log2FoldChange'] > 2) & (df['padj'] < 0.05)]

# Sort by fold change
sig = sig.sort_values('log2FoldChange', ascending=False)

# Preview
print(sig[['index', 'gene', 'log2FoldChange', 'padj']])

# Save
sig.to_csv(INTERMEDIATE_DIR / 'aminoglycoside_candidates.csv', index=False)
