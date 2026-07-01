from pathlib import Path
import pandas as pd
import re

DATA_ROOT = Path(__file__).parent.parent / 'data' 

# Parse GTF to extract ER3413 ID to gene name mapping
mapping = {}
with open(DATA_ROOT / 'ecoli_K12_MG1655.gtf.tar.gz', 'r') as f:
    for line in f:
        if 'old_locus_tag' in line and 'gene_id' in line:
            old = re.search(r'old_locus_tag "([^"]+)"', line)
            gene = re.search(r'gene "([^"]+)"', line)
            if old and gene:
                # Convert ER3413_1 format to ER3413_1
                old_id = old.group(1)
                mapping[old_id] = gene.group(1)

print(f"Mapped {len(mapping)} genes")

# Load candidates
candidates = pd.read_csv(Path(__file__).parent / 'gentamicin_candidates.csv')

# Map gene names
candidates['gene'] = candidates['gene_id'].map(mapping)

print(candidates[['gene_id', 'gene', 'log2FC']].head(20))
candidates.to_csv((Path(__file__).parent / 'gentamicin_candidates_named.csv'), index=False)