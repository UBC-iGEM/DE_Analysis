# network_analysis/

Regulatory network pipeline: maps DE candidates onto the RegulonDB v12
*E. coli* K-12 transcriptional regulatory network, scores candidates
by regulator specificity, and produces an interactive pyvis visualisation.

## Setup

### 1. Install dependencies

```bash
pip install networkx pyvis pandas
# or, if your system Python conflicts:
pip install networkx pyvis pandas --break-system-packages
```

### 2. Download RegulonDB flat files

Go to **regulondb.ccg.unam.mx/datasets** and download these two files into
`network_analysis/data/`:

| Download page label | Save as |
|---|---|
| **NetworkRegulatorGene** — "Regulatory Interactions" → TF-gene, TF-TU etc. (this is the merged file covering ALL regulator types: TFs, sigma factors are *not* in here, but small-molecule effectors like ppGpp are) | `network_regulator_gene.tsv` |
| **NetworkSigmaGene** — sigma factor → gene interactions | `network_sigma_gene.tsv` |

The sigma-gene file is essential for RpoE→rybB and RpoH→ibpA/ibpB edges —
without it, your best aminoglycoside candidates will show `biosensor_flag = ?`
despite having clean, well-characterised upstream regulation.

### Real file format (confirmed against RegulonDB Release 14.5 / v12.0 downloads)

Both files share the same structure: a ~20-line `#`-prefixed license/citation
preamble, then **one header row that does NOT start with `#`** — it looks like
`1)colName\t2)colName...` — followed by tab-separated data.

**`network_regulator_gene.tsv`** (7 columns):
```
1  regulatorId
2  regulatorName        ← TF name, OR small-molecule effector name (e.g. "ppGpp")
3  RegulatorGeneName     ← gene encoding the regulator. EMPTY for non-protein
                            effectors — this is how the pipeline tells TF from effector.
4  regulatedGeneId
5  regulatedGeneName     ← target gene  (note: column 5, not column 2)
6  function              ← '+' activate / '-' repress / '-+' dual / '' unknown
7  confidenceLevel       ← 'C' confirmed / 'S' strong / 'W' weak / '?' unknown
```

**`network_sigma_gene.tsv`** (5 columns):
```
1  sigmaName             ← coded name, e.g. "sigma24" — NOT "RpoE"
2  regulatedGeneName     ← target gene
3  function              ← same encoding as above
4  promoterEvidence      ← bracketed evidence-code list (unused)
5  confidenceLevel        ← same encoding as above
```

Sigma codes are mapped to common names internally (standard *E. coli*
nomenclature by molecular weight):

| Code | Common name | Role |
|---|---|---|
| sigma19 | FecI | iron-citrate transport |
| sigma24 | **RpoE** | envelope/extracytoplasmic stress → **rybB, ompG** |
| sigma28 | FliA | flagellar genes |
| sigma32 | **RpoH** | heat shock / protein QC → **ibpA, ibpB** |
| sigma38 | RpoS | stationary phase / general stress |
| sigma54 | RpoN | nitrogen metabolism |
| sigma70 | RpoD | housekeeping |

### 3. Ensure DE output files exist

The pipeline reads from `../caz_kan_DE/`:
```
betalactam_primary.csv   ← from deseq2_cazgit.py
kanamycin_primary.csv    ← from deseq2_kan.py
```
Run the DE scripts first if these files are missing.

---

## Usage

Run all three steps from inside `network_analysis/`:

```bash
cd network_analysis/

# Step 1: parse RegulonDB + build graph
python build_network.py --top-n 50 --min-fc 2.0

# Step 2: score regulators and annotate candidates
python score_candidates.py

# Step 3: visualise (opens browser)
python visualize_network.py --open
```

### Key parameters

| Flag | Script | Default | Effect |
|------|--------|---------|--------|
| `--top-n` | build | 50 | Keep top-N candidates per class by log₂FC |
| `--min-fc` | build | 2.0 | Hard FC floor (2.0 = primary candidates only) |
| `--min-tf-degree` | build | 1 | Min candidates a regulator must touch to be included |
| `--min-confidence` | build | `W` | RegulonDB evidence floor: `C` > `S` > `W`. Raise to `S` for stricter, cloning-grade evidence only |
| `--no-effectors` | build | off | Drop small-molecule effectors (e.g. ppGpp) — keep only TF/sigma binding sites |
| `--cross-threshold` | score | 0.25 | Specificity margin to flag a regulator as cross-reactive |
| `--open` | viz | off | Auto-open HTML in browser after export |

### Recommended first run

```bash
# Compact, cloning-grade view: strong+confirmed evidence only, top 30 per class
python build_network.py --top-n 30 --min-tf-degree 2 --min-confidence S
python score_candidates.py
python visualize_network.py --open
```

---

## Output files

```
output/
  regulatory_network.pkl       # networkx DiGraph (for downstream use)
  node_table.csv                # flat node attribute table
  tf_specificity_scores.csv     # per-regulator class specificity + candidate counts
  annotated_candidates.csv      # candidates with regulatory flags
  scoring_summary.txt           # printed summary report
  regulatory_network.html       # interactive pyvis visualisation
```

### Key columns in `annotated_candidates.csv`

| Column | Description |
|--------|-------------|
| `gene` | Gene symbol (lowercase) |
| `group` | `caz` / `kan` / `cross` |
| `fc_caz` / `fc_kan` | log₂FC values |
| `regulators` | Comma-sep list of all regulators touching this gene |
| `best_tf` | Most class-specific regulator among them |
| `best_tf_type` | `TF` / `sigma` / `effector` — whether `best_tf` has a discrete, clonable DNA binding site (TF/sigma) or is a small molecule acting on RNAP globally (effector, e.g. ppGpp) |
| `best_tf_spec` | That regulator's specificity score (0–1) |
| `has_cross_reactive_regulator` | True if any regulator is cross-reactive |
| `regulatory_clarity` | Composite score (0–1); higher = better reporter |
| `biosensor_flag` | ✓ recommended / ⊙ possible / ✗ avoid / ? unknown |
| `note` | Flags effector-only regulation, or "no regulator in RegulonDB" |

### Key columns in `tf_specificity_scores.csv`

| Column | Description |
|--------|-------------|
| `tf` | Regulator name |
| `regulator_type` | `TF` / `sigma` / `effector` |
| `specificity_caz` / `specificity_kan` | Fraction of targets in each DE class |
| `is_cross_reactive` | True if this regulator touches both classes |
| `dominant_class` | `beta-lactam` / `aminoglycoside` / `mixed` |

---

## Visual encoding (regulatory_network.html)

- **Node colour**: beta-lactam (coral) / aminoglycoside (teal) / cross-reactive (purple) / regulator (grey)
- **Node shape**: ellipse with plain border = TF · ellipse with dashed border = sigma factor · diamond = small-molecule effector
- **Edge style**: solid green = activates · dashed red = represses · dashed orange = dual
- **Hover tooltip**: full regulatory detail including RegulonDB confidence level

---

## Known limitations

- **RegulonDB coverage is uneven.** Well-studied regulators (LexA, CpxR) have
  many annotated targets; genes like `dgcZ` may have no regulator annotated in
  v12 at all. A `biosensor_flag = ?` means "not absent — just uncharacterised",
  not "ruled out".

- **Small-molecule effectors are real regulatory information, but not clonable
  the same way.** ~300 of the 319 distinct regulators in `NetworkRegulatorGene.tsv`
  have no encoding gene (e.g. ppGpp, the stringent-response alarmone) — they act
  directly on RNA polymerase rather than via a discrete DNA binding site. The
  pipeline keeps these by default (they're scientifically informative — e.g. a
  ppGpp connection to a candidate hints at a stringent-response/translational-
  stress link relevant to aminoglycosides) but flags them via `regulator_type`
  and a caveat note, since they can't be cloned as an isolated promoter element
  the way a TF or sigma binding site can. Use `--no-effectors` to exclude them
  entirely if you only want clonable candidates in the graph.

- **The 896-gene kanamycin list.** `--top-n 50` is the default to keep the
  graph legible. Lower `--min-fc` or raise `--top-n` to explore more candidates,
  but the visualisation becomes a hairball above ~150 nodes.

- **Confidence levels matter for cloning decisions.** A `W` (weak)-confidence
  interaction is still curated evidence, but `--min-confidence S` is worth
  trying once you're narrowing toward an actual cloning shortlist.

---

## Extending the pipeline

The pickled networkx graph is the stable interchange format. Future extensions
can load it directly:

```python
import pickle, networkx as nx
G = pickle.load(open('output/regulatory_network.pkl', 'rb'))
```

Planned next steps (see project Notion):
- **iModulon integration** via `pymodulon`: add iModulon membership as a node
  attribute and co-module membership as a second edge type.
- **Evo 2 likelihood scoring**: add as an additional column in
  `annotated_candidates.csv` once the candidate list is narrowed.
- **Promoter sequence retrieval**: fetch TSS coordinates and -10/-35 elements
  from RegulonDB's PromoterSet download for candidates flagged ✓.
