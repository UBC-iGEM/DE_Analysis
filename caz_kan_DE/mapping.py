from pathlib import Path
import gzip
import re
import pandas as pd

DATA_ROOT = Path(__file__).parent.parent / 'data' 

mapping = {}
with gzip.open(DATA_ROOT / 'ecoli_K12_MG1655.gtf.tar.gz', 'rt') as f:
    for line in f:
        gene_id = re.search(r'gene_id "([^"]+)"', line)
        gene = re.search(r'gene "([^"]+)"', line)
        if gene_id and gene:
            mapping[gene_id.group(1)] = gene.group(1)

pd.DataFrame(list(mapping.items()), columns=['gene_id', 'gene']).to_csv(
    'gene_mapping.csv', index=False)
print(f"Mapped {len(mapping)} genes")