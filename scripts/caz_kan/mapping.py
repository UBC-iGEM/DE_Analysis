import gzip
import re
import pandas as pd
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
EXTRACTED_RAW_DIR = DATA_DIR / 'GSE220559_RAW'
STANDARDIZED_COUNTS = DATA_DIR / 'standardized' / 'counts.csv'
GTF_PATH = DATA_DIR / 'ecoli_k12.gtf.gz'
OUTPUT_PATH = INTERMEDIATE_DIR / 'gene_mapping.csv'
INTERMEDIATE_DIR.mkdir(parents=True, exist_ok=True)
SCRATCH_DIR.mkdir(parents=True, exist_ok=True)


def ensure_raw_dir():
    if EXTRACTED_RAW_DIR.exists():
        return EXTRACTED_RAW_DIR
    if not RAW_DIR.exists() and RAW_ARCHIVE.exists():
        RAW_DIR.mkdir()
        with tarfile.open(RAW_ARCHIVE) as archive:
            archive.extractall(RAW_DIR)
    return RAW_DIR


def write_identity_mapping():
    if STANDARDIZED_COUNTS.exists():
        counts = pd.read_csv(STANDARDIZED_COUNTS)
        gene_ids = counts['gene_id'].astype(str)
        pd.DataFrame({'gene_id': gene_ids, 'gene': gene_ids}).to_csv(OUTPUT_PATH, index=False)
        print(f"{GTF_PATH.name} not found; wrote ID-only mapping for {len(gene_ids)} genes")
        return

    raw_dir = ensure_raw_dir()
    count_files = sorted(raw_dir.glob('*.txt*'))
    if not count_files:
        raise FileNotFoundError(
            f"Could not find {GTF_PATH} or raw count files in {RAW_DIR}. "
            "Add ecoli_k12.gtf.gz for gene names, or keep GSE220559_RAW.tar "
            f"in {DATA_DIR} so an ID-only mapping can be generated."
        )
    counts = pd.read_csv(count_files[0], sep=r'\s+', index_col=0)
    gene_ids = counts.index.to_series().astype(str)
    pd.DataFrame({'gene_id': gene_ids, 'gene': gene_ids}).to_csv(OUTPUT_PATH, index=False)
    print(f"{GTF_PATH.name} not found; wrote ID-only mapping for {len(gene_ids)} genes")

mapping = {}
if GTF_PATH.exists():
    with gzip.open(GTF_PATH, 'rt') as f:
        for line in f:
            gene_id = re.search(r'gene_id "([^"]+)"', line)
            gene = re.search(r'gene "([^"]+)"', line)
            if gene_id and gene:
                mapping[gene_id.group(1)] = gene.group(1)

    pd.DataFrame(list(mapping.items()), columns=['gene_id', 'gene']).to_csv(
        OUTPUT_PATH, index=False)
    print(f"Mapped {len(mapping)} genes")
else:
    write_identity_mapping()
