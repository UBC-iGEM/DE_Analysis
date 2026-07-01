"""
score_candidates.py
───────────────────
Loads the regulatory network built by build_network.py and computes two
complementary scores:

  TF Specificity Score
  ────────────────────
  For each TF/sigma factor node, measures how class-specific its DE targets
  are. A TF with score → 1.0 for beta-lactam exclusively regulates caz
  candidates; a TF with score ≈ 0.5 co-regulates both classes equally and
  is therefore a cross-reactivity risk.

      specificity_caz(TF) = n_caz_targets / (n_caz_targets + n_kan_targets)
      specificity_kan(TF) = 1 − specificity_caz(TF)

  Candidate Regulatory Clarity Score
  ────────────────────────────────────
  For each candidate gene, summarises the regulatory context:
  - which TFs regulate it
  - the most class-specific TF among those regulators
  - whether any regulator is cross-reactive (risk flag)
  - a composite 'regulatory_clarity' score (0–1, higher = more specific)

Usage
-----
    python score_candidates.py --graph output/regulatory_network.pkl

Outputs (to output/):
    tf_specificity_scores.csv
    annotated_candidates.csv
    scoring_summary.txt
"""

import argparse
import os
import pickle
import textwrap

import networkx as nx
import pandas as pd


# ── 1. Load graph ─────────────────────────────────────────────────────────────

def load_graph(path: str) -> nx.DiGraph:
    with open(path, 'rb') as fh:
        G = pickle.load(fh)
    print(f"[load] graph: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")
    return G


# ── 2. TF specificity scoring ─────────────────────────────────────────────────

def compute_tf_specificity(G: nx.DiGraph,
                           cross_threshold: float = 0.25) -> pd.DataFrame:
    """
    For every TF/sigma node, count how many of its targets fall in each
    DE class and derive specificity scores.

    cross_threshold : a TF is flagged 'cross_reactive' if both
                      n_caz_targets ≥ 1 AND n_kan_targets ≥ 1
                      AND specificity is within this margin of 0.5.
                      Lower = stricter definition of cross-reactivity.
    """
    rows = []
    for node, data in G.nodes(data=True):
        if data.get('group') != 'tf':
            continue

        targets = list(G.successors(node))
        target_groups = {t: G.nodes[t].get('group', 'unknown') for t in targets}

        n_caz   = sum(1 for g in target_groups.values() if g == 'caz')
        n_kan   = sum(1 for g in target_groups.values() if g == 'kan')
        n_cross = sum(1 for g in target_groups.values() if g == 'cross')
        n_total = n_caz + n_kan + n_cross

        if n_total == 0:
            continue  # TF regulates candidates not in either class list (shouldn't happen)

        regulator_type = data.get('regulator_type', 'TF')

        # Specificity: how much does this TF lean toward one class?
        # Cross-list genes count toward both for penalty purposes.
        spec_caz = (n_caz + n_cross) / (n_caz + n_kan + n_cross)
        spec_kan = (n_kan + n_cross) / (n_caz + n_kan + n_cross)

        # Cross-reactivity flag: TF touches both classes
        is_cross_reactive = (n_caz + n_cross > 0) and (n_kan + n_cross > 0)

        # Dominant class label
        if n_caz > n_kan:
            dominant = 'beta-lactam'
        elif n_kan > n_caz:
            dominant = 'aminoglycoside'
        else:
            dominant = 'mixed'

        rows.append({
            'tf':                node,
            'regulator_type':    regulator_type,
            'n_caz_targets':     n_caz,
            'n_kan_targets':     n_kan,
            'n_cross_targets':   n_cross,
            'n_total_targets':   n_total,
            'specificity_caz':   round(spec_caz, 3),
            'specificity_kan':   round(spec_kan, 3),
            'is_cross_reactive': is_cross_reactive,
            'dominant_class':    dominant,
            'target_genes':      ', '.join(sorted(targets)),
        })

    df = pd.DataFrame(rows).sort_values('n_total_targets', ascending=False)
    print(f"[score] {len(df)} TF/sigma nodes scored")
    print(f"        {df['is_cross_reactive'].sum()} cross-reactive TFs flagged")
    return df


# ── 3. Candidate annotation ───────────────────────────────────────────────────

def annotate_candidates(G: nx.DiGraph,
                        tf_scores: pd.DataFrame) -> pd.DataFrame:
    """
    For each candidate gene, annotate with its regulatory context.

    Returns a DataFrame sorted within each class by log2FC (descending),
    with regulatory information appended.

    Key output columns
    ──────────────────
    gene              : gene symbol
    group             : caz / kan / cross
    fc_caz            : log2FC under ceftazidime
    fc_kan            : log2FC under kanamycin
    regulators        : comma-separated list of TFs regulating this gene
    best_tf           : most class-specific TF for this gene's group
    best_tf_spec      : that TF's specificity score
    has_cross_reactive_regulator : True if any regulator is cross-reactive
    regulatory_clarity : composite score (0–1); higher = better biosensor choice
                         penalises cross-reactive regulators, rewards high FC
    biosensor_flag    : ✓ / ⊙ / ✗ based on regulatory clarity + FC
    """
    tf_spec_map = dict(zip(tf_scores['tf'], zip(tf_scores['specificity_caz'],
                                                 tf_scores['specificity_kan'],
                                                 tf_scores['is_cross_reactive'],
                                                 tf_scores['regulator_type'])))

    rows = []
    for node, data in G.nodes(data=True):
        group = data.get('group', '')
        if group == 'tf':
            continue

        fc_caz = data.get('fc_caz', 0.0)
        fc_kan = data.get('fc_kan', 0.0)
        regulators = list(G.predecessors(node))

        if not regulators:
            # Candidate has no annotated regulator in RegulonDB — still include
            rows.append({
                'gene': node, 'group': group,
                'fc_caz': fc_caz, 'fc_kan': fc_kan,
                'regulators': '',
                'best_tf': '',
                'best_tf_spec': None,
                'has_cross_reactive_regulator': False,
                'regulatory_clarity': None,
                'biosensor_flag': '?',
                'note': 'no regulator in RegulonDB',
            })
            continue

        # Determine which specificity column to use based on this gene's class
        is_caz = group in ('caz', 'cross')
        is_kan = group in ('kan', 'cross')

        # For each regulator, get its class-relevant specificity + type
        reg_specs = []
        cross_flags = []
        for tf in regulators:
            if tf in tf_spec_map:
                spec_caz, spec_kan, is_cr, reg_type = tf_spec_map[tf]
                relevant_spec = spec_caz if group == 'caz' else spec_kan
                reg_specs.append((tf, relevant_spec, reg_type))
                cross_flags.append(is_cr)
            else:
                reg_specs.append((tf, None, 'unknown'))
                cross_flags.append(False)

        # Best TF = highest relevant specificity. On ties, prefer an actual
        # TF/sigma factor (clonable binding site) over a small-molecule
        # effector (no discrete binding site to clone).
        valid = [(tf, s, t) for tf, s, t in reg_specs if s is not None]
        type_priority = {'TF': 2, 'sigma': 2, 'effector': 0, 'unknown': 1}
        if valid:
            best_tf, best_spec, best_type = max(
                valid, key=lambda x: (x[1], type_priority.get(x[2], 0)))
        else:
            best_tf, best_spec, best_type = '', None, ''

        has_cross = any(cross_flags)

        # Regulatory clarity score (0–1)
        # = best_tf_specificity * (1 - cross_reactive_penalty)
        if best_spec is not None:
            cross_penalty = 0.4 if has_cross else 0.0
            clarity = best_spec * (1.0 - cross_penalty)
        else:
            clarity = None

        # Biosensor flag
        if clarity is None:
            flag = '?'
        elif clarity >= 0.75:
            flag = '✓'
        elif clarity >= 0.5:
            flag = '⊙'
        else:
            flag = '✗'

        note = ''
        if best_type == 'effector':
            note = 'best regulator is a small-molecule effector (no discrete binding site to clone)'

        rows.append({
            'gene': node,
            'group': group,
            'fc_caz': fc_caz,
            'fc_kan': fc_kan,
            'regulators': ', '.join(sorted(regulators)),
            'best_tf': best_tf,
            'best_tf_type': best_type,
            'best_tf_spec': round(best_spec, 3) if best_spec else None,
            'has_cross_reactive_regulator': has_cross,
            'regulatory_clarity': round(clarity, 3) if clarity else None,
            'biosensor_flag': flag,
            'note': note,
        })

    df = pd.DataFrame(rows)

    # Sort: within each group, by FC descending for the relevant antibiotic
    df['sort_fc'] = df.apply(
        lambda r: r['fc_caz'] if r['group'] == 'caz' else r['fc_kan'], axis=1)
    df = df.sort_values(['group', 'sort_fc'], ascending=[True, False])
    df = df.drop(columns='sort_fc')

    return df


# ── 4. Summary report ─────────────────────────────────────────────────────────

def print_summary(tf_scores: pd.DataFrame,
                  candidate_scores: pd.DataFrame,
                  out_dir: str) -> None:
    lines = []
    lines.append("═" * 68)
    lines.append("  REGULATORY NETWORK SCORING SUMMARY")
    lines.append("═" * 68)

    # Top TFs per class
    for class_label, col, group in [('Beta-lactam', 'specificity_caz', 'caz'),
                                     ('Aminoglycoside', 'specificity_kan', 'kan')]:
        relevant = tf_scores[
            (tf_scores['dominant_class'].str.contains(
                'beta-lactam' if group == 'caz' else 'amino') |
             (tf_scores['n_caz_targets' if group == 'caz' else 'n_kan_targets'] > 0))
        ].nlargest(5, col)
        lines.append(f"\nTop TFs for {class_label} candidates (by specificity):")
        lines.append(f"  {'TF':<14} {'spec':>6}  {'n_targets':>9}  cross-reactive?")
        lines.append("  " + "─" * 48)
        for _, row in relevant.iterrows():
            n = row['n_caz_targets'] if group == 'caz' else row['n_kan_targets']
            cr = "⚠ yes" if row['is_cross_reactive'] else "no"
            lines.append(f"  {row['tf']:<14} {row[col]:>6.3f}  {n:>9}  {cr}")

    # Candidate summary
    lines.append("\n" + "─" * 68)
    lines.append("  CANDIDATE BIOSENSOR FLAGS")
    lines.append("─" * 68)
    for group, label in [('caz', 'Beta-lactam'), ('kan', 'Aminoglycoside'), ('cross', 'Cross-reactive')]:
        sub = candidate_scores[candidate_scores['group'] == group]
        lines.append(f"\n{label} candidates ({len(sub)} total):")
        lines.append(f"  {'gene':<12} {'flag':>5}  {'fc':>6}  {'best_tf':<14}  {'clarity':>7}")
        lines.append("  " + "─" * 56)
        for _, row in sub.head(15).iterrows():
            fc = row['fc_caz'] if group == 'caz' else row['fc_kan']
            lines.append(
                f"  {row['gene']:<12} {row['biosensor_flag']:>5}  {fc:>6.2f}  "
                f"{str(row['best_tf']):<14}  "
                f"{str(row['regulatory_clarity']) if row['regulatory_clarity'] else '   n/a':>7}"
            )
        if len(sub) > 15:
            lines.append(f"  … and {len(sub) - 15} more (see annotated_candidates.csv)")

    lines.append("\n" + "═" * 68)
    report = '\n'.join(lines)
    print(report)

    report_path = os.path.join(out_dir, 'scoring_summary.txt')
    with open(report_path, 'w') as fh:
        fh.write(report + '\n')
    print(f"\n[save] summary → {report_path}")


# ── 5. Entry point ────────────────────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--graph', default='output/regulatory_network.pkl',
                   help='Pickled graph from build_network.py')
    p.add_argument('--cross-threshold', type=float, default=0.25,
                   help='Specificity margin to flag a TF as cross-reactive (default: 0.25)')
    p.add_argument('--out', default='output/',
                   help='Output directory')
    return p.parse_args()


def main():
    args = parse_args()
    if not os.path.exists(args.graph):
        import sys
        sys.exit(f"[error] graph file not found: {args.graph}\n"
                 f"        Run build_network.py first.")

    G = load_graph(args.graph)
    tf_scores = compute_tf_specificity(G, cross_threshold=args.cross_threshold)
    candidate_scores = annotate_candidates(G, tf_scores)

    os.makedirs(args.out, exist_ok=True)

    tf_path = os.path.join(args.out, 'tf_specificity_scores.csv')
    tf_scores.to_csv(tf_path, index=False)
    print(f"[save] TF scores ({len(tf_scores)} rows) → {tf_path}")

    cand_path = os.path.join(args.out, 'annotated_candidates.csv')
    candidate_scores.to_csv(cand_path, index=False)
    print(f"[save] annotated candidates ({len(candidate_scores)} rows) → {cand_path}")

    print_summary(tf_scores, candidate_scores, args.out)
    print("\nNext: python visualize_network.py")


if __name__ == '__main__':
    main()
