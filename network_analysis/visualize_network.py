"""
visualize_network.py
────────────────────
Renders the scored regulatory network as an interactive HTML file using pyvis.

Node design
───────────
  caz candidates     coral/orange   size ∝ fc_caz
  kan candidates     teal/green     size ∝ fc_kan
  cross-reactive     purple         size ∝ max(fc_caz, fc_kan)
  TF / sigma factor  grey outline   fixed size; dashed border for sigma factors

Edge design
───────────
  activates   solid green arrow
  represses   dashed red arrow
  dual        dashed orange arrow

Tooltip (on hover) includes:
  - gene name, class, log2FC values
  - biosensor flag and regulatory clarity
  - comma-separated regulators

Usage
-----
    python visualize_network.py \
        --graph  output/regulatory_network.pkl \
        --scores output/annotated_candidates.csv \
        --tf-scores output/tf_specificity_scores.csv \
        --out    output/regulatory_network.html
"""

import argparse
import os
import pickle
import webbrowser

import networkx as nx
import pandas as pd

try:
    from pyvis.network import Network
except ImportError:
    raise SystemExit(
        "[error] pyvis not installed.\n"
        "        pip install pyvis --break-system-packages\n"
        "        or:  pip install pyvis"
    )


# ── Colour palette (matches the D3 widget we built earlier) ──────────────────

COLORS = {
    'caz':   {'bg': '#FAECE7', 'border': '#993C1D', 'font': '#71290F'},
    'kan':   {'bg': '#E1F5EE', 'border': '#0F6E56', 'font': '#085041'},
    'cross': {'bg': '#EEEDFE', 'border': '#534AB7', 'font': '#3C3489'},
    'tf':    {'bg': '#F1EFE8', 'border': '#5F5E5A', 'font': '#444441'},
}

EDGE_COLORS = {
    'activates': '#1D9E75',
    'represses': '#E24B4A',
    'dual':      '#BA7517',
}


# ── Node sizing ───────────────────────────────────────────────────────────────

def node_size(group: str, fc_caz: float, fc_kan: float) -> int:
    if group == 'tf':
        return 22
    fc = max(fc_caz, fc_kan)
    # Scale: fc=2 → size 18, fc=8 → size 36, capped at 42
    return min(42, max(16, int(10 + fc * 3.2)))


# ── Tooltip HTML ─────────────────────────────────────────────────────────────

def make_tooltip(gene: str, row: pd.Series | None, tf_row: pd.Series | None) -> str:
    if tf_row is not None:
        # Regulator node (TF, sigma factor, or small-molecule effector)
        reg_type = tf_row.get('regulator_type', 'TF')
        lines = [
            f"<b>{gene}</b> ({reg_type})",
            f"Class specificity:",
            f"  beta-lactam: {tf_row.get('specificity_caz', '?'):.3f}",
            f"  aminoglycoside: {tf_row.get('specificity_kan', '?'):.3f}",
            f"Targets in DE list: {int(tf_row.get('n_total_targets', 0))}",
            f"  ceftazidime: {int(tf_row.get('n_caz_targets', 0))}",
            f"  kanamycin:   {int(tf_row.get('n_kan_targets', 0))}",
            f"Cross-reactive: {'⚠ yes' if tf_row.get('is_cross_reactive') else 'no'}",
        ]
        if reg_type == 'effector':
            lines.append("⚠ small-molecule effector — no discrete DNA binding "
                         "site to clone as an isolated promoter element")
        return '<br>'.join(lines)

    if row is None:
        return gene

    flag = row.get('biosensor_flag', '?')
    best_tf_type = row.get('best_tf_type', '')
    tf_label = f"{row.get('best_tf', '—')}" + (f" [{best_tf_type}]" if best_tf_type else "")
    lines = [
        f"<b>{gene}</b> [{row.get('group', '?').upper()}]",
        f"log₂FC  ceftazidime: {row.get('fc_caz', 0):.2f}",
        f"log₂FC  kanamycin:   {row.get('fc_kan', 0):.2f}",
        f"Regulators: {row.get('regulators', 'none')}",
        f"Best TF: {tf_label} (spec {row.get('best_tf_spec', '?')})",
        f"Cross-reactive regulator: {'⚠ yes' if row.get('has_cross_reactive_regulator') else 'no'}",
        f"Regulatory clarity: {row.get('regulatory_clarity', '?')}",
        f"Biosensor flag: {flag}",
    ]
    if row.get('note'):
        lines.append(f"Note: {row['note']}")
    return '<br>'.join(lines)


# ── Build pyvis network ───────────────────────────────────────────────────────

def build_pyvis(G: nx.DiGraph,
                candidate_scores: pd.DataFrame,
                tf_scores: pd.DataFrame) -> 'Network':
    """Construct and style the pyvis Network object."""

    net = Network(
        height='720px',
        width='100%',
        directed=True,
        bgcolor='#ffffff',
        font_color='#333333',
        notebook=False,
        cdn_resources='in_line',   # self-contained HTML, no external CDN needed
    )

    # Index score tables for O(1) lookup
    cand_idx = candidate_scores.set_index('gene') if not candidate_scores.empty else pd.DataFrame()
    tf_idx   = tf_scores.set_index('tf')          if not tf_scores.empty else pd.DataFrame()

    # ── Nodes ─────────────────────────────────────────────────────────────────
    for node, data in G.nodes(data=True):
        group = data.get('group', 'tf')
        fc_caz = data.get('fc_caz', 0.0)
        fc_kan = data.get('fc_kan', 0.0)
        regulator_type = data.get('regulator_type')  # 'TF' | 'sigma' | 'effector' | None

        col = COLORS.get(group, COLORS['tf'])
        size = node_size(group, fc_caz, fc_kan)

        # Look up scores
        cand_row = cand_idx.loc[node] if (not cand_idx.empty and node in cand_idx.index) else None
        tf_row   = tf_idx.loc[node]   if (not tf_idx.empty   and node in tf_idx.index)   else None

        tooltip = make_tooltip(node, cand_row, tf_row)

        # Label: include FC for candidates to aid visual scanning
        if group != 'tf':
            fc = fc_caz if group == 'caz' else fc_kan
            label = f"{node}\n{fc:.1f}"
        else:
            label = node

        # Shape encodes regulator subtype: sigma factors get a dashed double
        # border (distinct regulatory mechanism), small-molecule effectors
        # get a diamond (no discrete DNA binding site — not directly
        # clonable the way a TF/sigma binding site is).
        if regulator_type == 'effector':
            shape = 'diamond'
            border_width, border_dashes = 1, False
        elif regulator_type == 'sigma':
            shape = 'ellipse'
            border_width, border_dashes = 2, True
        else:
            shape = 'ellipse'
            border_width, border_dashes = 1, False

        net.add_node(
            node,
            label=label,
            title=tooltip,
            color={
                'background': col['bg'],
                'border':     col['border'],
                'highlight':  {'background': col['bg'], 'border': '#111111'},
            },
            font={'color': col['font'], 'size': 12, 'face': 'monospace'},
            size=size,
            shape=shape,
            borderWidth=border_width,
            borderWidthSelected=3,
            shapeProperties={'borderDashes': border_dashes},
        )

    # ── Edges ─────────────────────────────────────────────────────────────────
    for src, tgt, data in G.edges(data=True):
        etype = data.get('edge_type', 'activates')
        color = EDGE_COLORS.get(etype, '#B4B2A9')
        is_dashed = etype in ('represses', 'dual')

        net.add_edge(
            src, tgt,
            title=etype,
            color={'color': color, 'highlight': '#111111'},
            dashes=is_dashed,
            arrows='to',
            width=1.5,
            smooth={'type': 'curvedCW', 'roundness': 0.1},
        )

    return net


# ── Physics and interaction options ──────────────────────────────────────────

PHYSICS_OPTIONS = """{
  "physics": {
    "solver": "forceAtlas2Based",
    "forceAtlas2Based": {
      "gravitationalConstant": -80,
      "centralGravity": 0.008,
      "springLength": 140,
      "springConstant": 0.06,
      "damping": 0.4
    },
    "stabilization": {
      "enabled": true,
      "iterations": 200,
      "updateInterval": 25
    }
  },
  "interaction": {
    "hover": true,
    "tooltipDelay": 150,
    "navigationButtons": true,
    "keyboard": true
  },
  "edges": {
    "smooth": {"type": "dynamic"}
  }
}"""


# ── Legend HTML injected into the output ─────────────────────────────────────

LEGEND_HTML = """
<div style="font-family:monospace;font-size:12px;line-height:1.8;
            padding:10px 16px;border:1px solid #ddd;border-radius:6px;
            background:#fafafa;margin-bottom:8px;max-width:720px">
  <b>Node colour</b>&nbsp;
  <span style="background:#FAECE7;border:1px solid #993C1D;padding:1px 6px;border-radius:3px">
    beta-lactam specific</span>&nbsp;
  <span style="background:#E1F5EE;border:1px solid #0F6E56;padding:1px 6px;border-radius:3px">
    aminoglycoside specific</span>&nbsp;
  <span style="background:#EEEDFE;border:1px solid #534AB7;padding:1px 6px;border-radius:3px">
    cross-reactive</span>&nbsp;
  <span style="background:#F1EFE8;border:1px solid #5F5E5A;padding:1px 6px;border-radius:3px">
    regulator</span><br>
  <b>Regulator shape</b>&nbsp; ellipse, plain border = TF &nbsp;|&nbsp;
  ellipse, dashed border = sigma factor &nbsp;|&nbsp;
  diamond = small-molecule effector (e.g. ppGpp — no discrete DNA binding
  site, not directly clonable the way a TF/sigma site is)<br>
  <b>Node size</b>&nbsp; proportional to log₂FC (candidate nodes) &nbsp;|&nbsp;
  <b>Label</b>&nbsp; gene name + log₂FC&nbsp;|&nbsp;
  <b>Hover</b>&nbsp; full regulatory detail incl. RegulonDB confidence level<br>
  <b>Edge</b>&nbsp;
  <span style="color:#1D9E75">━ activates</span>&nbsp;&nbsp;
  <span style="color:#E24B4A">╌ represses</span>&nbsp;&nbsp;
  <span style="color:#BA7517">╌ dual</span><br>
  <b>Biosensor flags</b>&nbsp; ✓ recommended &nbsp; ⊙ possible &nbsp; ✗ avoid &nbsp; ? not in RegulonDB
</div>
"""


# ── Export with legend ────────────────────────────────────────────────────────

def export_html(net: 'Network', path: str, candidate_scores: pd.DataFrame) -> None:
    """Write interactive HTML, injecting a legend and candidate summary table."""
    net.set_options(PHYSICS_OPTIONS)

    # Generate raw pyvis HTML
    html = net.generate_html()

    # Build a compact candidate table to embed below the graph
    table_rows = []
    for _, row in candidate_scores.iterrows():
        fc = row['fc_caz'] if row['group'] == 'caz' else row['fc_kan']
        cross_warn = ' ⚠' if row.get('has_cross_reactive_regulator') else ''
        table_rows.append(
            f"<tr>"
            f"<td style='font-family:monospace'>{row['gene']}</td>"
            f"<td>{row['group'].upper()}</td>"
            f"<td>{fc:.2f}</td>"
            f"<td style='font-family:monospace'>{row.get('best_tf','—')}</td>"
            f"<td>{str(row.get('regulatory_clarity','?'))}{cross_warn}</td>"
            f"<td style='font-size:16px'>{row['biosensor_flag']}</td>"
            f"</tr>"
        )

    table_html = f"""
<details style="margin-top:12px;font-family:sans-serif;font-size:13px">
  <summary style="cursor:pointer;font-weight:500">
    Candidate scoring table ({len(candidate_scores)} genes) — click to expand
  </summary>
  <div style="overflow-x:auto;margin-top:8px">
  <table style="border-collapse:collapse;width:100%;font-size:12px">
    <thead style="background:#f0f0f0">
      <tr>
        <th style="text-align:left;padding:4px 8px">Gene</th>
        <th style="text-align:left;padding:4px 8px">Class</th>
        <th style="text-align:right;padding:4px 8px">log₂FC</th>
        <th style="text-align:left;padding:4px 8px">Best TF</th>
        <th style="text-align:right;padding:4px 8px">Clarity</th>
        <th style="text-align:center;padding:4px 8px">Flag</th>
      </tr>
    </thead>
    <tbody>
      {''.join(table_rows)}
    </tbody>
  </table>
  </div>
</details>
"""

    # Inject legend + table before </body>
    html = html.replace(
        '</body>',
        LEGEND_HTML + table_html + '</body>'
    )

    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, 'w', encoding='utf-8') as fh:
        fh.write(html)
    print(f"[save] interactive HTML → {path}")


# ── Entry point ───────────────────────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--graph',     default='output/regulatory_network.pkl')
    p.add_argument('--scores',    default='output/annotated_candidates.csv')
    p.add_argument('--tf-scores', default='output/tf_specificity_scores.csv')
    p.add_argument('--out',       default='output/regulatory_network.html')
    p.add_argument('--open',      action='store_true',
                   help='Open the HTML in your default browser after export')
    return p.parse_args()


def main():
    args = parse_args()
    for path, name in [(args.graph, '--graph'),
                       (args.scores, '--scores'),
                       (getattr(args, 'tf_scores'), '--tf-scores')]:
        if not os.path.exists(path):
            import sys
            sys.exit(f"[error] {name} file not found: {path}\n"
                     f"        Run build_network.py and score_candidates.py first.")

    G = pickle.load(open(args.graph, 'rb'))
    candidate_scores = pd.read_csv(args.scores)
    tf_scores        = pd.read_csv(getattr(args, 'tf_scores'))

    net = build_pyvis(G, candidate_scores, tf_scores)
    export_html(net, args.out, candidate_scores)

    if args.open:
        webbrowser.open(f'file://{os.path.abspath(args.out)}')
    else:
        print(f"\nOpen in browser:\n  {os.path.abspath(args.out)}")


if __name__ == '__main__':
    main()
