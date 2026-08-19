# Regulatory network analysis

This pipeline builds a typed *RegulonDB regulator → candidate* graph from all
datasets listed in `config/datasets.json`, enriches it with same-release
iModulon/PRECISE evidence, and writes a candidate table plus an interactive
HTML network. It is a heuristic screening aid: expression/activity values are
reported as evidence and burden values are proxies, not measured burden.

## Data contract and setup

Run commands from the repository root. Install the network-specific runtime:

```bash
python -m pip install -r network_analysis/requirements.txt
```

Copy `config/network_assets.example.json` to `config/network_assets.json` and
review the provider, release, license, and citation metadata. The source
manifest records where each asset comes from; `setup_data.py --download`
generates `config/network_assets.lock.json` with the exact retrieved bytes,
release metadata, upstream IDs, sizes, retrieval times, and SHA256 values.
Hashes do not need to be captured manually before the first download.

The setup validator requires:

- RegulonDB `>= 14.5.0` for regulator, sigma, and identity products. The
  default GraphQL provider queries `getDatabaseInfo`, `listAllFileNames`, and
  `getDataOfFile`; its endpoint is configurable because the documented host is
  currently a prerelease host.
- PRECISE-1K dataset release `>= 1.0`, with the model and companion expression
  file pinned to the same source revision. The iModulonDB Version 2.5.0
  platform release is recorded as contextual catalog provenance, not confused
  with the PRECISE-1K dataset release.

`pymodulon==0.2.1` is the compatibility runtime used by the recovered loader.
The model supplies iModulon matrices and metadata; gene-level expression is
loaded from `IcaData.log_tpm`, then a documented `IcaData.X`, then the
companion expression file.

Download missing files and generate the lock, verify an existing checkout, or
explicitly accept changed upstream bytes:

```bash
python network_analysis/setup_data.py --manifest config/network_assets.json --download
python network_analysis/setup_data.py --manifest config/network_assets.json
python network_analysis/setup_data.py --manifest config/network_assets.json --download --refresh-lock
```

For an offline checkout, retain both the downloaded files and the lock. The
validator never silently accepts a changed cached file.

The downloader verifies HTTPS certificates. In a managed environment with a
private certificate authority, point `SSL_CERT_FILE` at the organization’s CA
bundle; do not disable certificate verification.

The identity mapping is not optional scientifically for b-number/probe inputs:
pass the same-release mapping TSV/CSV with `--mapping`. Without it the loader
retains the source identifier and adds `identity_mapping_not_verified`; those
nodes may not join RegulonDB symbols.

## Build and candidate-direction policy

The builder reads every configured `data/<dataset>/standardized/de_results.csv`
and retains all rows as provenance. Candidate seeding is selectable:

```bash
python network_analysis/build_network.py \
  --regulator data/network_regulator_gene.tsv \
  --sigma data/network_sigma_gene.tsv \
  --mapping data/network_gene_mapping.tsv \
  --candidate-direction upregulated

python network_analysis/build_network.py \
  --regulator data/network_regulator_gene.tsv \
  --sigma data/network_sigma_gene.tsv \
  --mapping data/network_gene_mapping.tsv \
  --candidate-direction either-direction
```

`upregulated` is the default and seeds a gene if it is significant and
upregulated in any dataset. `either-direction` seeds on significant induction
or repression. Both modes retain the opposite-direction observations as
evidence; significant opposing directions are marked `conflicted` and never
count as corroboration.

Evidence tiers are transparent labels, not weighted scores:

1. `conflicted` — significant opposing directions are present;
2. `limited` — the only qualifying dataset is tobramycin;
3. `corroborated` — direction-consistent significant evidence covers every
   configured dataset in at least one class;
4. `supported` — other qualifying evidence.

The graph stores objective dataset counts, support fractions, per-class tiers,
source row IDs, padj source, identity source, and caveats. Supported class keys
are `beta_lactam`, `aminoglycoside`, `fluoroquinolone`, `polymyxin`, `cross`,
and `tf`; a duplicate dataset name is rejected to prevent double-counting.
PRECISE-1K has matched ciprofloxacin/control samples, but the selected model
does not provide matched polymyxin samples, so polymyxin expression/activity
fields remain unavailable rather than being inferred. Regulatory edges use
`activates`, `represses`, or `dual`; iModulon co-membership edges use
`co-imodulon` and are never scored as regulators.

## iModulon/PRECISE annotation

The annotation stage accepts embedded gene expression when available and falls
back to the companion expression matrix:

```bash
python network_analysis/i_modulon_analysis.py \
  --graph network_analysis/output/regulatory_network.pkl \
  --precise data/precise1k/model.json.gz \
  --expression data/precise1k/expression.csv \
  --mapping data/network_gene_mapping.tsv \
  --expression-units 'log2(TPM)' \
  --expression-normalization 'quality-controlled, uncentered PRECISE-1K expression'
```

`--expression` is optional when the serialized `IcaData` contains a usable
gene-level `log_tpm` or `X` matrix. It remains required for the selected
PRECISE-1K model because that artifact does not contain gene-level expression.
`ica.A` is iModulon activity, not gene expression.

For each candidate it preserves two separate evidence families:

- `gene_expression`: basal, induced, delta, per-background values/counts,
  units, normalization, availability, and within-class basal/induced
  percentiles from the companion matrix;
- `imodulon_activity`: the same basal/induced/delta structure from `ica.A`.

It also retains all weighted iModulon memberships, identifies the primary
iModulon separately, and reports `metabolic_burden_proxy` and
`translation_burden_proxy` with their system-category/iModulon basis. Lower
proxy values are preferable, but these are heuristic screening indicators.

## Score and visualize

```bash
python network_analysis/score_candidates.py
python network_analysis/visualize_network.py
```

The scorer filters predecessors by `node_type=regulator` and the three
regulatory edge types before computing class specificity. The visualizer
escapes all tooltip values and exposes evidence tiers, per-dataset evidence,
expression, activity, caveats, and burden proxies. Generated files belong in
the ignored `network_analysis/output/` directory and should not be committed.
