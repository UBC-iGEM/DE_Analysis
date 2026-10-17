# PRISM — Promoter Response Index and Selection Map

PRISM links *E. coli* transcriptomic responses, curated regulatory interactions,
and promoter annotations in a searchable SQLite database. It helps a team
trace **why a promoter is worth testing**, inspect the evidence, and record
manual review decisions. It was developed for FAST's antibiotic-responsive
sensor design and provides a reusable example for other teams.

**Scope:** this is a research and design tool. The current example uses
antibiotic-class proxies, not measurements of the four wet-lab antibiotics
(amoxicillin, cephalexin, gentamicin, and tobramycin). It does not predict clinical
susceptibility or establish promoter fluorescence, specificity, or sensor speed.

## Try the explorer without installing dependencies

From a checkout, use Python 3.11 or 3.12:

```bash
python3 -m promoter_discovery.demo
python3 -m promoter_discovery.explorer --database results/demo.sqlite
```

Open **http://127.0.0.1:8000**. The separate demo contains two explicitly
fictional candidates and demonstrates database relationships and graph browsing.
It contains no experimental claims and refuses to overwrite an existing file.
To view the real database already generated in this repository's local results:

```bash
python3 -m promoter_discovery.explorer
```

Search by gene, class, support tier, adjusted p-value, or effect size. Select a
candidate for response tables, promoter sequences, fragment checks, and the
regulatory graph. The regulator selector keeps the graph readable; **All
regulators** gives the overview. Click a solid connection to inspect individual
reference records. Expand the assessment section for stability, cross-reactivity,
and normalized-count replicate summaries. The HTTP interface is read-only,
binds to localhost, and is intended for local use.

The explorer follows FAST's brand palette and uses bundled Geist Bold headings
and Source Serif 4 body text, so typography works offline. The typefaces retain
their SIL Open Font Licenses in `promoter_discovery/static/fonts/`.

## Reproduce the scientific example

```bash
python3.11 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements-lock.txt
python -m promoter_discovery.workflow --download
```

The workflow runs reference verification/acquisition → joint DE fit → network
integration → candidate scoring → promoter mapping → assessment and database
creation → interactive network export. With references cached, omit `--download`.
Use `--reuse-fit` to rebuild downstream outputs from an existing verified fit.
A failed stage stops the workflow. Database construction writes a temporary file
and replaces the previous database only after validation succeeds.

`requirements-lock.txt` records the validated Python 3.11 dependency versions.
`requirements.txt` keeps broader development constraints; dependency updates
should pass the tests and a fresh benchmark run before replacing the lock.

The original GSE220559 count archive is tracked; reference downloads and results
are not. The locked RegulonDB release is **14.5.0**, and the genome is MG1655
**U00096.3**. If the upstream service no longer supplies the locked bytes, obtain
the matching release from its provider or deliberately review a reference update;
do not change hashes simply to make verification pass.

For individual-stage debugging, use this order (assessment builds the database):

```bash
python -m promoter_discovery.setup_data --download
python -m promoter_discovery.run_benchmark
python -m promoter_discovery.build_network
python -m promoter_discovery.score_candidates
python -m promoter_discovery.promoter_selection
python -m promoter_discovery.assess_candidates
python -m promoter_discovery.visualize_network
```

`python -m promoter_discovery.build_database` explicitly rebuilds only the query
layer from existing outputs; it never refits DE or recalculates network statistics.
`assess_candidates --skip-database` produces only assessment files. Only load
serialized network graphs produced by a trusted pipeline run: pickle files are
executable Python serialization, not an interchange format for untrusted uploads.

## What the current example produces

| Evidence | Current benchmark |
| --- | ---: |
| Unique RNA-seq samples | 18: three per condition, including shared water controls |
| Genes retained in the joint fit | 4,176 |
| Comparisons | 5 drug/control + 10 drug/drug |
| Candidate genes | 925 |
| Proposed first screening panel | 12 promoters: six per antibiotic class |
| Candidate-linked promoter/fragment reviews | 473 |
| Annotation/genome sequence matches | 448; 25 lack coordinates |
| Full-reference curated interactions | 6,460 |
| Full-reference promoters / transcription units | 4,057 / 3,780 |
| Reference interactions without promoter context | 729 |

These are dataset-specific output counts, not accuracy or clinical-performance
metrics. Candidates are genes linked to possible promoters; 925 does not mean
925 independent promoters. The panel is a greedy, balanced first-pass screening
proposal with no repeated promoter or connected transcription-unit group. It is
not the future 2,304-construct library. Sequence matching verifies reference
consistency; it does not mean a fragment is ready to order.

## Evidence and calculation choices

- **Comparable input:** one GSE220559 study and one PyDESeq2 count model,
  `~condition`, with shared controls counted once. Genes require at least ten
  counts in three samples. Candidate seeds use adjusted p-value <0.05 and
  unshrunken log₂ fold change >2; raw and shrunken effects remain separate.
- **Explicit roles:** ceftazidime is the beta-lactam discovery proxy; imipenem
  provides within-study robustness; kanamycin is the aminoglycoside discovery
  proxy. Ciprofloxacin and polymyxin E challenge cross-reactivity. Observed drugs
  and all four planned wet-lab targets have distinct records; no one-to-one
  proxy equivalence is asserted.
- **Regulatory context:** enrichment uses all eligible measured, mapped genes,
  rather than only selected candidates. Fisher tests and BH correction provide
  exploratory gene-level evidence; correlated operon genes limit interpretation.
- **Honest graph semantics:** solid edges point to the actual curated promoter,
  TU, or gene target. Dashed edges are annotated membership. A gene/TU-target
  interaction does not establish direct regulation of an associated promoter;
  gene-target effects are not assigned to sibling genes. Missing links remain
  searchable. Individual effects, evidence, source IDs, and release are retained.
- **Assessment:** nine post-fit screens vary count/effect thresholds without
  model refits. Direct comparisons also have BH correction over all ten pairs
  and genes. An off-target bound of |log₂ effect| <1 is descriptive, not an
  equivalence test. Proposed fragments cover known promoter-linked sites and
  still require operator/boundary review and wet-lab validation.

One dataset, one strain, three replicates per condition, shared controls, and an
asymmetric number of class proxies limit generalization. Transcript abundance
can change without promoter induction. Exact-drug time/dose reporter experiments
and susceptible/resistant controls are the next biological validation step.

## Database structure and use

```text
Source counts + locked references
             ↓
DE + network + promoter + assessment artifacts (CSV/JSON/FASTA/plots)
             ↓ explicit importer; source hashes; atomic replacement
results/promoter_discovery.sqlite (schema 11)
             ↓
Python API / CLI / local read-only explorer

Versioned reviews → results/promoter_discovery.curation.sqlite
```

The relational layer stores runs, assets, genes, antibiotics, contrasts, DE
results, regulatory actors/interactions/sites, promoters, TUs and memberships,
candidates/evidence/panels, QC, replicate metrics, and assessments. Stable text
IDs, foreign keys, and join tables express many-to-many relationships. Numeric
fields support filtering; original annotation and assessment attributes are
retained in JSON columns. This uses SQLite with JSON attributes, not a separate
NoSQL or graph database.

```bash
python -m promoter_discovery.db_cli candidates --class beta_lactam --padj-max 0.05
python -m promoter_discovery.db_cli --json candidate pdhr
python -m promoter_discovery.db_cli --json panel
python -m promoter_discovery.db_cli --json interactions CRP --target-kind promoter
python -m promoter_discovery.db_cli --json promoter-paths CRP
python -m promoter_discovery.db_cli --json network pdhr
python -m promoter_discovery.db_cli --json coverage
```

```python
from promoter_discovery.db_api import Database

with Database("results/promoter_discovery.sqlite", read_only=True) as db:
    hits = db.search_candidates(antibiotic_class="beta_lactam", padj_max=0.05)
    candidate = db.get_candidate("pdhr")
    paths = db.reference_promoter_paths("CRP")
```

Queries default to the latest completed run; pass `run_id` in Python or
`--run-id` on run-aware CLI commands to select another run. Generated databases
are snapshots, not an accumulated multi-release reference archive. Rebuilding
retains old review history in the companion file, not old scientific snapshots.

### Manual curation

Measurements remain imported evidence. Reviews append status, priority, notes,
construct readiness, editor, timestamp, reason, and optional source reference.

```bash
python -m promoter_discovery.db_cli --json candidate pdhr
# Use the run ID returned above:
python -m promoter_discovery.db_cli review-candidate pdhr reviewed \
  --run-id <run-id> --editor "Your name" --reason "Checked reference evidence"
python -m promoter_discovery.db_cli --json history candidate pdhr
```

Candidate statuses: `proposed`, `reviewed`, `approved`, `excluded`. Promoter
statuses: `review`, `validated`, `ready_for_order`, `excluded`; use
`review-promoter <promoter-id> <status>` with the same audit fields.
**Keep and back up `*.curation.sqlite` files.** Rebuilds restore reviews only for
identical input content, schema, and software release; changes require fresh approval. Earlier history
remains searchable, including reviews preserved during schema upgrades.

## Repository and contribution

| Location | Purpose |
| --- | --- |
| `promoter_discovery/` | Workflow, importer, schema, queries, curation, explorer |
| `config/` | Active example configuration; source manifests and hash locks |
| `data/` | Tracked source archive; ignored cached references; [inventory](data/README.md) |
| `results/` | Ignored generated artifacts and local curation archive; [guide](results/README.md) |
| `tests/` | Synthetic scientific-contract, database, CLI, and workflow regressions |

Wiki design assets and write-up drafts stay in the ignored `wiki_assets/` folder.
The optional PRECISE-1K context module and its separate manifests remain available;
they are not part of the primary antibiotic benchmark. The custom dataset runner
is `python -m promoter_discovery.input_data --config <config>`; the importer and
assessment currently support the documented benchmark layout. General ingestion
of arbitrary studies requires an adapter and is not claimed by this release.

To contribute, use a small branch and pull request, describe the biological or
software problem, and include a focused regression for changes to scientific
meaning, imports, or queries. Do not commit generated data, local reviews, or
unreviewed reference-lock changes. Run:

```bash
python -m pytest -q
python -m promoter_discovery.setup_data
```

CI uses synthetic fixtures and does not download references. A clean workflow
run against cached, locked references is an additional release check. The future
screening layer can link construct/version → plate/well → exposure/replicate →
measurement to existing candidate and promoter IDs; this experimental layer is
not implemented yet.

## License and citation

Project-authored code and documentation use the [MIT license](LICENSE).
**Third-party data and typefaces are excluded from that grant.** RegulonDB data, the GEO source
archive, PRECISE-1K, and the GenBank genome retain their providers' terms; see
`config/*assets.json` and the [data inventory](data/README.md). Review those terms
before distributing a populated database. The fictional demo avoids redistribution
of reference records. Use [CITATION.cff](CITATION.cff) for software attribution and
also cite the original dataset and methods:

- [GSE220559 study, Microbiology Spectrum (2023)](https://pmc.ncbi.nlm.nih.gov/articles/PMC10100721/)
- [RegulonDB v12.0 resource paper, Nucleic Acids Research (2024)](https://pmc.ncbi.nlm.nih.gov/articles/PMC10767902/) — background citation; the imported release is 14.5.0.
- [PyDESeq2, Bioinformatics (2023)](https://pmc.ncbi.nlm.nih.gov/articles/PMC10502239/)

For an iGEM Software submission, the team must also host the reviewed release in
its dedicated iGEM GitLab software repository and link its wiki software page;
a GitHub repository alone does not meet the current hosting requirement.
See the [official Software award requirements](https://competition.igem.org/judging/special-prizes).
