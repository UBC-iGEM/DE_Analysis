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
replace every owner-supplied value. The manifest is deliberately explicit and
must record, for each asset: URL, local path, exact version, official release
revision, retrieval date, SHA256, license, and citation. The setup validator
requires:

- RegulonDB `>= 14.5.0` for both regulator/sigma flat files and the matching
  K-12 identity mapping. A preamble declaring `14.5` is normalized to
  `14.5.0`.
- An owner-identified iModulon/PRECISE product whose exact release/version is
  `>= 2.5.0`, plus its companion gene-expression matrix from the same release.
  The model and expression assets must have matching `version` and `asset_id`.
  The public [iModulonDB update page](https://imodulondb.org/updates) describes a Version 2.5.0 database release;
  if that is the intended source, record it explicitly as an `imodulondb` asset
  and provide the corresponding downloadable model/activity and expression
  artifacts rather than silently treating the legacy PRECISE-1K JSON as that
  release.

The second requirement is not a claim that `pymodulon` itself has a 2.5.0
package release. `pymodulon==0.2.1` is the compatibility runtime used by the
recovered loader; the data-product identifier must be supplied by the data
owner. The script fails rather than guessing a provider or substituting an
unrelated package version.

Validate an already downloaded set, or download missing files from explicit
manifest URLs:

```bash
python network_analysis/setup_data.py --manifest config/network_assets.json
python network_analysis/setup_data.py --manifest config/network_assets.json --download
```

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
source row IDs, padj source, identity source, and caveats. Its class keys are
`beta_lactam`, `aminoglycoside`, `cross`, and `tf`. Regulatory edges use
`activates`, `represses`, or `dual`; iModulon co-membership edges use
`co-imodulon` and are never scored as regulators.

## iModulon/PRECISE annotation

The annotation stage requires both a model and the companion expression matrix:

```bash
python network_analysis/i_modulon_analysis.py \
  --graph network_analysis/output/regulatory_network.pkl \
  --precise data/imodulon/model.json.gz \
  --expression data/imodulon/expression.csv.gz \
  --mapping data/network_gene_mapping.tsv
```

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
