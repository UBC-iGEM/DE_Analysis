# Antibiotic-responsive promoter discovery

Identify *E. coli* genes and regulatory systems that suggest promoters worth
testing for biosensor design. Wet-lab antibiotics are **amoxicillin, cephalexin,
gentamicin, and tobramycin**: two beta-lactams and two aminoglycosides.

The primary analysis combines one GSE220559 RNA-seq benchmark, a joint PyDESeq2
count model, and same-release RegulonDB network and promoter annotations.
Ceftazidime/imipenem and kanamycin are class proxies; ciprofloxacin and polymyxin
E challenge cross-reactivity. This dataset does not directly test the four
wet-lab antibiotics.

## Project layout

```text
promoter_discovery/                 Python workflow and shared input loaders
config/                            Benchmark settings and reference manifests/locks
data/
  raw/GSE220559/                    Original RNA-seq archive
  references/regulondb/             Verified regulatory and promoter references
results/
  differential_expression/
    model/                         Joint fit, filtering, pairwise contrasts, provenance
    contrasts/<treatment>/         Matched drug/control DE tables and selected inputs
  regulatory_network/              Graph, enrichment, provenance, interactive HTML
  promoter_candidates/             Annotated candidates and class shortlists
  promoter_review/                 Gene/TU/promoter mappings and construct review
  candidate_assessment/            QC figures, stability, cross-reactivity, proposed panel
tests/                             Scientific-contract and integration regressions
```

Cached references and generated results are ignored by Git. The original
GSE220559 archive is retained. See the [data inventory](data/README.md) and
[results guide](results/README.md).

## Run the workflow

Use Python >=3.11. From the repository root:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python -m promoter_discovery.run_benchmark
python -m promoter_discovery.setup_data --download
python -m promoter_discovery.build_network
python -m promoter_discovery.score_candidates
python -m promoter_discovery.promoter_selection
python -m promoter_discovery.assess_candidates
python -m promoter_discovery.build_database
python -m promoter_discovery.visualize_network
```

Reference acquisition requires network access. Once cached, omit `--download`
to verify the files. The default settings use `config/benchmark.json`,
`config/regulondb_assets.json`, and `config/regulondb.lock.json`.

Reference sources and checksums are recorded in the manifests and locks under
`config/`. The separate PRECISE-1K manifest provides optional external context.
The primary references also include regulatory sites and the accession-pinned
MG1655 genome (U00096.3) for fragment review.
The optional custom per-dataset runner is
`python -m promoter_discovery.input_data --config <config>`; the primary
benchmark requires the joint model command above.

## View the results

Start with `beta_lactam_candidates.csv` and `aminoglycoside_candidates.csv` in
`results/promoter_candidates/`, then `results/promoter_review/promoter_review.csv`.
Open CSV files in a spreadsheet application. On macOS, view the network with:

```bash
open results/regulatory_network/regulatory_network.html
```

`results/candidate_assessment/experimental_panel.csv` proposes six promoters
per class, with no repeated promoter or connected transcription-unit group.
The same folder contains sample PCA/distances, replicate summaries and plots,
original-versus-shrunken rankings, cross-reactivity comparisons, operon support,
fragment review, and promoter context diagrams. Run the assessment after
regenerating DE, network, scoring, and promoter outputs.

```bash
python -m promoter_discovery.assess_candidates --panel-per-class 6 --upstream 200 --downstream 30 --off-target-limit 1
```

Stability measures nine post-fit screening choices within the existing fitted
gene universe, without model refits. Direct comparisons additionally use BH
correction across all ten pairs and genes. The off-target limit is a configurable
descriptive log2 effect bound (1 means twofold), not a validated specificity
threshold or an equivalence test. Operon support reports concentration in
connected annotated TUs; it does not supply correlation-adjusted p-values.
Proposed fragments extend a TSS window to cover known promoter-linked sites.
They include overlapping gene context and a genome/annotation match
check, and still require manual boundary review.

The class lists may overlap. Support labels describe measured comparisons,
not independent validation or calibrated biosensor performance. Annotated
promoter sequences need operator and fragment-boundary review before ordering;
wet-lab measurements establish transfer to the actual four antibiotics.

The generated SQLite database is written to `results/promoter_discovery.sqlite`
after candidate assessment. It indexes the DE, regulatory, promoter, candidate,
panel, QC, and fragment-review outputs without replacing the source files. It
also imports the full locked RegulonDB promoter, transcription-unit, and
interaction sets. Individual interactions retain their target type, effect,
site, evidence, and source release; the existing gene-level network remains a
summary for candidate analysis. Use
the Python API (`promoter_discovery.db_api.Database`) or the search CLI:

```bash
python -m promoter_discovery.db_cli candidates --class beta_lactam
python -m promoter_discovery.db_cli --json candidate gfcc
python -m promoter_discovery.db_cli --json panel
python -m promoter_discovery.db_cli --json interactions CRP --target-kind promoter
python -m promoter_discovery.db_cli --json promoter-paths CRP
```

The database includes the five drug/control and ten direct drug/drug contrasts.
Candidate evidence links back to individual measured genes, including names
with multiple locus tags. Candidate searches with an effect or adjusted-p-value
filter use discovery comparisons; challenge responses remain available in each
candidate's full evidence. The saved full-reference graph is required for a
complete database build.

Use `--skip-database` on `assess_candidates` when only the assessment files are
needed; `build_database` remains available for an explicit rebuild.

## Verification

```bash
python -m pytest -q
python -m promoter_discovery.setup_data
```

Tests use synthetic fixtures and do not download reference data. Cached
reference validation checks the release and file hashes. Do not edit lock
hashes to accept unexplained changes.
