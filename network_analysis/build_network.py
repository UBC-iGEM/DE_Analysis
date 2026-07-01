"""
build_network.py
────────────────
Loads RegulonDB v12 flat files and DE results, constructs a networkx DiGraph
of regulator → candidate-gene interactions, and saves it for downstream
scoring and visualisation steps.

RegulonDB v12 file formats (confirmed against real downloads, June 2026)
──────────────────────────────────────────────────────────────────────
NetworkRegulatorGene.tsv  (7 tab-separated columns, ~20-line license preamble,
                            then a "1)colName\t2)colName..." header row)
    1  regulatorId
    2  regulatorName        ← regulator (TF name, or small-molecule effector
                               name e.g. "ppGpp")
    3  RegulatorGeneName    ← gene encoding the regulator. EMPTY for non-protein
                               effectors — this is how we detect TF vs effector.
    4  regulatedGeneId
    5  regulatedGeneName    ← target gene  (NOTE: column 5, not column 2!)
    6  function             ← '+' activate / '-' repress / '-+' dual / '' unknown
    7  confidenceLevel      ← 'C' confirmed / 'S' strong / 'W' weak / '?' unknown

NetworkSigmaGene.tsv  (5 tab-separated columns, same preamble style)
    1  sigmaName            ← coded name, e.g. "sigma24" (NOT "RpoE")
    2  regulatedGeneName    ← target gene
    3  function             ← same encoding as above
    4  promoterEvidence     ← bracketed evidence-code list (unused here)
    5  confidenceLevel      ← same encoding as above

Sigma codes map to common names by molecular weight (standard E. coli
nomenclature) — see SIGMA_NAME_MAP below. sigma24 = RpoE = the envelope
stress sigma factor we need for rybB/ompG; sigma32 = RpoH for ibpA/ibpB.

Usage
-----
    python build_network.py \
        --caz       ../caz_kan_DE/betalactam_primary.csv \
        --kan       ../caz_kan_DE/kanamycin_primary.csv \
        --regulator data/network_regulator_gene.tsv \
        --sigma     data/network_sigma_gene.tsv \
        --top-n 50 \
        --out   output/

Run from inside network_analysis/.
"""

import argparse
import os
import pickle
import sys

import networkx as nx
import pandas as pd


# ── 1. Constants ───────────────────────────────────────────────────────────

EFFECT_MAP = {
    '+':  'activates',
    '-':  'represses',
    '-+': 'dual',
    '+-': 'dual',     # defensive: some RegulonDB exports use this order instead
    '?':  'unknown',
    '':   'unknown',
}

CONFIDENCE_RANK = {'C': 3, 'S': 2, 'W': 1, '?': 0, '': 0}

# Standard E. coli sigma factor naming: code (by approx. kDa) → common name.
# Confirmed codes present in RegulonDB v12 NetworkSigmaGene.tsv: 19/24/28/32/38/54/70.
SIGMA_NAME_MAP = {
    'sigma19': 'fecI',   # σ19  — iron-citrate transport
    'sigma24': 'rpoE',   # σE/σ24 — envelope/extracytoplasmic stress (rybB, ompG)
    'sigma28': 'fliA',   # σF/σ28 — flagellar genes
    'sigma32': 'rpoH',   # σH/σ32 — heat shock / protein quality control (ibpA, ibpB)
    'sigma38': 'rpoS',   # σS/σ38 — stationary phase / general stress
    'sigma54': 'rpoN',   # σN/σ54 — nitrogen metabolism
    'sigma70': 'rpoD',   # σD/σ70 — housekeeping / primary sigma
}


# ── 2. RegulonDB parsing ──────────────────────────────────────────────────

def _is_header_or_comment(parts: list[str]) -> bool:
    """
    RegulonDB files have a ~20-line '#'-prefixed license/citation preamble,
    followed by ONE column-header row that does NOT start with '#' but
    instead looks like '1)colName\\t2)colName...'. Both need skipping.
    """
    if not parts:
        return True
    first = parts[0].strip()
    if first.startswith('#') or first == '':
        return True
    # Header row pattern: "1)regulatorId", "2)regulatorName", etc.
    if len(first) >= 2 and first[0].isdigit() and ')' in first[:4]:
        return True
    return False


def _parse_regulator_gene_file(path: str) -> list[dict]:
    """Parse NetworkRegulatorGene.tsv (7-column format, target in col 5)."""
    rows = []
    with open(path, encoding='utf-8') as fh:
        for line in fh:
            line = line.rstrip('\n')
            if not line:
                continue
            parts = line.split('\t')
            if _is_header_or_comment(parts):
                continue
            if len(parts) < 7:
                continue

            regulator_name = parts[1].strip()
            regulator_gene = parts[2].strip()   # empty ⇒ non-protein effector
            target         = parts[4].strip()
            effect         = parts[5].strip()
            confidence     = parts[6].strip()

            if not regulator_name or not target:
                continue

            rows.append({
                'regulator':      regulator_name.lower(),
                'target':         target.lower(),
                'effect':         effect,
                'confidence':     confidence if confidence else '?',
                'regulator_type': 'TF' if regulator_gene else 'effector',
                'source':         'regulator-gene',
            })
    return rows


def _parse_sigma_gene_file(path: str) -> list[dict]:
    """Parse NetworkSigmaGene.tsv (5-column format, coded sigma names)."""
    rows = []
    with open(path, encoding='utf-8') as fh:
        for line in fh:
            line = line.rstrip('\n')
            if not line:
                continue
            parts = line.split('\t')
            if _is_header_or_comment(parts):
                continue
            if len(parts) < 5:
                continue

            sigma_code = parts[0].strip()
            target     = parts[1].strip()
            effect     = parts[2].strip()
            confidence = parts[4].strip()

            if not sigma_code or not target:
                continue

            regulator = SIGMA_NAME_MAP.get(sigma_code.lower(), sigma_code)

            rows.append({
                'regulator':      regulator.lower(),
                'target':         target.lower(),
                'effect':         effect,
                'confidence':     confidence if confidence else '?',
                'regulator_type': 'sigma',
                'source':         'sigma-gene',
            })
    return rows


def load_regulondb(regulator_gene_path: str,
                   sigma_gene_path: str | None = None,
                   min_confidence: str = 'W') -> pd.DataFrame:
    """
    Load and combine RegulonDB regulator-gene and sigma-gene flat files.

    min_confidence : keep interactions with confidence ≥ this level.
                     Rank order: C (confirmed) > S (strong) > W (weak) > ? (unknown).
                     Default 'W' keeps everything with any curated evidence.
                     Use 'S' to require strong-or-better evidence only.

    Returns a DataFrame with columns:
        regulator, target, effect, edge_type, confidence, regulator_type, source
    Duplicate (regulator, target) pairs — common when RegulonDB has multiple
    evidence lines for the same interaction — are collapsed to the single
    highest-confidence row.
    """
    rows = _parse_regulator_gene_file(regulator_gene_path)
    if sigma_gene_path and os.path.exists(sigma_gene_path):
        rows += _parse_sigma_gene_file(sigma_gene_path)
    else:
        if sigma_gene_path:
            print(f"[warn] sigma-gene file not found: {sigma_gene_path}")
            print("       Sigma factor edges (rpoE→rybB, rpoH→ibpA/B) will be absent.")

    df = pd.DataFrame(rows)
    if df.empty:
        sys.exit("[error] No regulatory interactions parsed. Check file paths/format.")

    df['edge_type'] = df['effect'].map(EFFECT_MAP).fillna('unknown')
    df['conf_rank'] = df['confidence'].map(CONFIDENCE_RANK).fillna(0)

    # Apply confidence floor
    min_rank = CONFIDENCE_RANK.get(min_confidence, 0)
    n_before = len(df)
    df = df[df['conf_rank'] >= min_rank].copy()
    if len(df) < n_before:
        print(f"[regulondb] confidence filter (≥{min_confidence}): "
              f"{n_before:,} → {len(df):,} rows")

    # Drop unknown-direction edges — not useful for direction-sensitive scoring
    df = df[df['edge_type'] != 'unknown'].copy()

    # Collapse duplicate (regulator, target) pairs: keep highest-confidence row
    n_before = len(df)
    df = (df.sort_values('conf_rank', ascending=False)
            .drop_duplicates(subset=['regulator', 'target'], keep='first'))
    if len(df) < n_before:
        print(f"[regulondb] collapsed {n_before - len(df):,} duplicate "
              f"(regulator, target) evidence rows")

    n_tf       = (df.drop_duplicates('regulator')['regulator_type'] == 'TF').sum()
    n_sigma    = (df.drop_duplicates('regulator')['regulator_type'] == 'sigma').sum()
    n_effector = (df.drop_duplicates('regulator')['regulator_type'] == 'effector').sum()
    print(f"[regulondb] loaded {len(df):,} interactions "
          f"({df['regulator'].nunique()} regulators → {df['target'].nunique()} genes)")
    print(f"            regulator types: {n_tf} TF, {n_sigma} sigma, {n_effector} effector")

    return df


# ── 3. DE result loading ─────────────────────────────────────────────────

def load_de_results(caz_path: str,
                    kan_path: str,
                    top_n: int | None = None,
                    min_fc: float = 0.0) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Load primary candidate CSVs produced by deseq2_caz.py / deseq2_kan.py.

    Expected columns (from pydeseq2 output):
        gene, log2FoldChange, padj  (plus baseMean, lfcSE, stat — all ignored here)

    top_n  : keep only the top-N genes by log2FC per class (guards against
              the 896-gene kanamycin list blowing up the visualisation).
    min_fc : hard lower bound on log2FC (applied before top_n).
    """
    def _read(path, label):
        df = pd.read_csv(path, index_col=0)
        df.columns = df.columns.str.strip()
        if 'gene' not in df.columns:
            raise ValueError(f"No 'gene' column in {path}. "
                             f"Columns found: {list(df.columns)}")
        df = df.dropna(subset=['gene', 'log2FoldChange'])
        df['gene'] = df['gene'].str.strip().str.lower()
        df = df[df['log2FoldChange'] >= min_fc]
        if top_n:
            df = df.nlargest(top_n, 'log2FoldChange')
        print(f"[de] {label}: {len(df)} candidates "
              f"(log2FC ≥ {min_fc:.1f}{f', top {top_n}' if top_n else ''})")
        return df

    caz = _read(caz_path, 'ceftazidime')
    kan = _read(kan_path, 'kanamycin')
    return caz, kan


# ── 4. Graph construction ────────────────────────────────────────────────

def build_graph(df_caz: pd.DataFrame,
                df_kan: pd.DataFrame,
                regulondb: pd.DataFrame,
                min_tf_degree: int = 1,
                include_effectors: bool = True) -> nx.DiGraph:
    """
    Build a directed networkx graph:
        Nodes  — candidate DE genes + their regulators (TF / sigma / effector)
        Edges  — regulator → gene interactions from RegulonDB

    Node attributes
    ───────────────
    group           : 'caz' | 'kan' | 'cross' | 'tf'
                       (NOTE: 'tf' is the umbrella group for ALL regulator
                       nodes — TFs, sigma factors, AND small-molecule
                       effectors — kept for backward compatibility with
                       scoring code that filters on group=='tf'.)
    regulator_type  : 'TF' | 'sigma' | 'effector' (only set on regulator nodes;
                       lets you distinguish a clonable TF binding site from
                       a small-molecule effector like ppGpp that has no
                       discrete DNA binding site to clone)
    fc_caz / fc_kan : log2FC values (0 for regulator nodes)
    label           : display name

    Edge attributes
    ───────────────
    edge_type  : 'activates' | 'represses' | 'dual'
    confidence : 'C' | 'S' | 'W' | '?'
    source     : 'regulator-gene' | 'sigma-gene'

    min_tf_degree     : only add a regulator node if it regulates at least
                        this many candidates. Default 1 (include all).
    include_effectors : if False, drop small-molecule effectors (e.g. ppGpp)
                        entirely — useful if you want a graph of only
                        clonable TF/sigma binding sites.
    """
    G = nx.DiGraph()

    caz_fc = dict(zip(df_caz['gene'], df_caz['log2FoldChange']))
    kan_fc = dict(zip(df_kan['gene'], df_kan['log2FoldChange']))
    caz_genes = set(caz_fc)
    kan_genes = set(kan_fc)
    cross_genes = caz_genes & kan_genes

    for gene in caz_genes | kan_genes:
        if gene in cross_genes:
            group = 'cross'
        elif gene in caz_genes:
            group = 'caz'
        else:
            group = 'kan'
        G.add_node(gene,
                   group=group,
                   regulator_type=None,
                   fc_caz=round(caz_fc.get(gene, 0.0), 3),
                   fc_kan=round(kan_fc.get(gene, 0.0), 3),
                   label=gene)

    all_candidates = caz_genes | kan_genes
    relevant = regulondb[regulondb['target'].isin(all_candidates)].copy()

    if not include_effectors:
        n_before = len(relevant)
        relevant = relevant[relevant['regulator_type'] != 'effector']
        if len(relevant) < n_before:
            print(f"[graph] excluded {n_before - len(relevant)} effector "
                  f"interactions (--no-effectors)")

    reg_degrees = relevant.groupby('regulator')['target'].nunique()
    regs_to_add = set(reg_degrees[reg_degrees >= min_tf_degree].index)

    for _, row in relevant.iterrows():
        reg = row['regulator']
        gene = row['target']
        if reg not in regs_to_add:
            continue

        if reg not in G.nodes:
            G.add_node(reg, group='tf', regulator_type=row['regulator_type'],
                      fc_caz=0.0, fc_kan=0.0, label=reg)

        if not G.has_edge(reg, gene):
            G.add_edge(reg, gene,
                       edge_type=row['edge_type'],
                       confidence=row['confidence'],
                       source=row['source'])

    n_reg = sum(1 for _, d in G.nodes(data=True) if d['group'] == 'tf')
    n_candidate = len(G) - n_reg
    n_tf_nodes       = sum(1 for _, d in G.nodes(data=True) if d.get('regulator_type') == 'TF')
    n_sigma_nodes    = sum(1 for _, d in G.nodes(data=True) if d.get('regulator_type') == 'sigma')
    n_effector_nodes = sum(1 for _, d in G.nodes(data=True) if d.get('regulator_type') == 'effector')
    print(f"[graph] {n_candidate} candidate nodes, {n_reg} regulator nodes "
          f"({n_tf_nodes} TF, {n_sigma_nodes} sigma, {n_effector_nodes} effector), "
          f"{G.number_of_edges()} edges")

    expected = ['reca', 'cpxp', 'rybb', 'ibpb', 'ibpa', 'sula', 'recn', 'dgcz', 'pdel']
    missing = [g for g in expected if g not in G.nodes]
    if missing:
        print(f"[graph] note: expected candidates absent from graph: {missing}")
        print("         (check that they appear in the DE input files and pass FC filter)")

    return G


# ── 5. Persistence ───────────────────────────────────────────────────────

def save_graph(G: nx.DiGraph, out_dir: str) -> None:
    os.makedirs(out_dir, exist_ok=True)

    pickle_path = os.path.join(out_dir, 'regulatory_network.pkl')
    with open(pickle_path, 'wb') as fh:
        pickle.dump(G, fh)
    print(f"[save] graph → {pickle_path}")

    rows = [{'gene': node, **data} for node, data in G.nodes(data=True)]
    node_df = pd.DataFrame(rows)
    node_path = os.path.join(out_dir, 'node_table.csv')
    node_df.to_csv(node_path, index=False)
    print(f"[save] node table ({len(node_df)} rows) → {node_path}")


# ── 6. Entry point ───────────────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--caz',   default='../caz_kan_DE/betalactam_primary.csv',
                   help='Ceftazidime primary candidates CSV')
    p.add_argument('--kan',   default='../caz_kan_DE/kanamycin_primary.csv',
                   help='Kanamycin primary candidates CSV')
    p.add_argument('--regulator', default='data/network_regulator_gene.tsv',
                   help='RegulonDB NetworkRegulatorGene.tsv (7-column format)')
    p.add_argument('--sigma',     default='data/network_sigma_gene.tsv',
                   help='RegulonDB NetworkSigmaGene.tsv (5-column format, optional but recommended)')
    p.add_argument('--top-n', type=int, default=50,
                   help='Keep top-N candidates per class by log2FC (default: 50)')
    p.add_argument('--min-fc', type=float, default=2.0,
                   help='Minimum log2FC threshold (default: 2.0 = primary candidates)')
    p.add_argument('--min-tf-degree', type=int, default=1,
                   help='Only include regulators touching ≥ this many candidates (default: 1)')
    p.add_argument('--min-confidence', default='W', choices=['C', 'S', 'W'],
                   help='Minimum RegulonDB confidence level to include (default: W = all evidence)')
    p.add_argument('--no-effectors', action='store_true',
                   help='Exclude small-molecule effectors (e.g. ppGpp) — keep only TF/sigma regulators')
    p.add_argument('--out',   default='output/',
                   help='Output directory')
    return p.parse_args()


def main():
    args = parse_args()

    for path, name in [(args.caz, '--caz'), (args.kan, '--kan'), (args.regulator, '--regulator')]:
        if not os.path.exists(path):
            sys.exit(f"[error] {name} file not found: {path}\n"
                     f"        Run from inside network_analysis/ and check paths.")

    regulondb = load_regulondb(args.regulator, args.sigma, min_confidence=args.min_confidence)
    df_caz, df_kan = load_de_results(args.caz, args.kan,
                                     top_n=args.top_n,
                                     min_fc=args.min_fc)
    G = build_graph(df_caz, df_kan, regulondb,
                    min_tf_degree=args.min_tf_degree,
                    include_effectors=not args.no_effectors)
    save_graph(G, args.out)
    print("\nNext: python score_candidates.py")


if __name__ == '__main__':
    main()
