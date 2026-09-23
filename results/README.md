# Analysis results

Generated stage folders are local and ignored by Git. Recreate them with the
commands in the [project README](../README.md).

| Folder | Contents |
| --- | --- |
| `differential_expression/model/` | Joint sample metadata, raw/normalized/variance-stabilized counts, common gene filter, direct drug contrasts, and provenance |
| `differential_expression/contrasts/<treatment>/` | Full-precision matched drug/control DE results with separate original and shrunken effects, selected raw expression tables, and metadata |
| `regulatory_network/` | Full curated graph, node/source tables, directional regulon enrichment, coverage, scoring provenance, and interactive HTML |
| `promoter_candidates/` | `annotated_candidates.csv`, `beta_lactam_candidates.csv`, and `aminoglycoside_candidates.csv` |
| `promoter_review/` | `candidate_promoter_mapping.csv` and grouped `promoter_review.csv` |
| `candidate_assessment/` | QC and panel figures, replicate summaries, screening stability, cross-reactivity, operon support, proposed fragments, and `experimental_panel.csv` |
| `promoter_discovery.sqlite` | Generated SQLite query layer linking DE, network, promoter, candidate, panel, QC, and fragment-review records |

Start with the two class shortlists, then review the promoter table. Class lists
can overlap. The benchmark uses class proxies for the four wet-lab antibiotics;
promoter annotations require construct review and experimental validation.

Open `regulatory_network/regulatory_network.html` in a browser to explore the
candidate graph. The underlying serialized graph retains the full reference
network. CSV files can be opened in a spreadsheet application.

The folder migration preserved the existing joint DE fit. Its original config,
provenance, and output checksums are recorded in
`differential_expression/model/layout_migration.json`; updated paths are recorded
in the current provenance. Subsequent model runs write fresh fit provenance.

Optional PRECISE-1K graph annotations go in `regulatory_network/precise1k/`.
The documented optional scoring command writes to `precise1k_context/`.
The custom per-dataset runner writes diagnostic outputs to `custom_analysis/`.

Search the generated database with `python -m promoter_discovery.db_cli` or use
`promoter_discovery.db_api.Database` from Python. The CSV, JSON, and FASTA files
remain the reproducible source artifacts; the SQLite file is rebuilt from them.
