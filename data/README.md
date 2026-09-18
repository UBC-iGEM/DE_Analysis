# Source data and references

`raw/GSE220559/GSE220559_RAW.tar` is the only antibiotic-expression study
retained here. The full original archive is preserved, including samples
outside the selected benchmark. Study accessions name the source folders;
treatment names identify the generated contrasts under `results/`.

The active configuration, [`config/benchmark.json`](../config/benchmark.json),
selects:

| Treatment | Benchmark role |
| --- | --- |
| Ceftazidime | Beta-lactam discovery proxy |
| Imipenem | Within-study beta-lactam robustness |
| Kanamycin | Aminoglycoside discovery proxy |
| Ciprofloxacin | Cross-reactivity challenge |
| Polymyxin E | Cross-reactivity challenge |

Shared water controls occur once in the joint count model. The wet-lab targets
remain amoxicillin, cephalexin, gentamicin, and tobramycin; this study provides
class proxies rather than exact-drug validation.

`references/regulondb/` contains the cached regulator/gene and sigma/gene
networks, gene identifiers, transcription units, promoter annotations, and the
derived canonical gene mapping. `RISet.tsv` adds curated regulatory-site positions
and evidence. `references/ecoli_mg1655/U00096.3.fasta` supplies the accession-pinned
genome for proposed fragment review. Their sources and verified hashes are recorded
in [`regulondb_assets.json`](../config/regulondb_assets.json) and
[`regulondb.lock.json`](../config/regulondb.lock.json).

Optional PRECISE-1K assets use `references/precise1k/` and the separate
`config/precise1k_assets.json` manifest. They provide external regulatory
context rather than primary antibiotic-response evidence.

All generated analysis tables are in [`results/`](../results/README.md).
Follow the [project README](../README.md) to regenerate them.
