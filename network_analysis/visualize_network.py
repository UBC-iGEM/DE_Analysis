"""Render the typed regulatory network and expose all evidence families."""

from __future__ import annotations

import argparse
import html
import json
import pickle
from pathlib import Path
from typing import Any

import networkx as nx

try:
    from .dataset_registry import CLASS_REGISTRY
except ImportError:  # pragma: no cover
    from dataset_registry import CLASS_REGISTRY  # type: ignore

try:
    from pyvis.network import Network
except ImportError:  # pragma: no cover - gives a clear CLI error
    Network = None  # type: ignore[assignment,misc]


COLORS = {
    **{class_key: specification["color"] for class_key, specification in CLASS_REGISTRY.items()},
    "cross": "#805ad5",
    "tf": "#4a5568",
}
EDGE_COLORS = {"activates": "#2f855a", "represses": "#c53030", "dual": "#805ad5", "co-imodulon": "#a0aec0"}


def _value(value: Any) -> str:
    if isinstance(value, (dict, list, tuple, set)):
        return json.dumps(value, sort_keys=True, default=str)
    return "" if value is None else str(value)


def _section(title: str, value: Any) -> str:
    return f"<b>{html.escape(title)}</b>: {html.escape(_value(value))}<br>"


def build_tooltip(node: str, data: dict[str, Any]) -> str:
    """Build an escaped tooltip containing graph, evidence, and proxy data."""

    parts = [f"<b>{html.escape(str(node))}</b><br>", _section("Node type", data.get("node_type", "")), _section("Class", data.get("group", ""))]
    if data.get("node_type") == "candidate":
        for field in (
            "evidence_tier", "evidence_tier_by_class", "candidate_direction_policy", "n_datasets_observed",
            "n_datasets_significant", "direction_consistent", "dataset_evidence_json", "caveats",
            "evidence_quality_flags", "tobramycin_only", "gene_expression", "imodulon_activity",
            "imodulon_primary", "imodulons", "metabolic_burden_proxy", "translation_burden_proxy",
            "burden_proxy_system_categories", "burden_proxy_imodulons", "burden_proxy_basis",
        ):
            if field in data:
                parts.append(_section(field, data[field]))
        parts.append("<small>Expression: lower basal and higher induced values are preferred; percentiles are within-class descriptive summaries. Burden values are heuristic proxies; lower is preferable.</small>")
    else:
        parts.extend([_section("Regulator type", data.get("regulator_type", "")), _section("RegulonDB release", data.get("regulondb_release", ""))])
    return "".join(parts)


def render_graph(graph: nx.DiGraph, output_html: str | Path, title: str = "Regulatory network") -> Path:
    if Network is None:
        raise RuntimeError("pyvis is required for HTML visualization; install network_analysis/requirements.txt")
    output_path = Path(output_html)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    network = Network(height="900px", width="100%", directed=True, bgcolor="#ffffff", font_color="#1a202c", cdn_resources="in_line")
    network.barnes_hut()
    for node, data in graph.nodes(data=True):
        group = data.get("group", "tf")
        node_type = data.get("node_type", "regulator")
        network.add_node(node, label=str(node), title=build_tooltip(str(node), data), color=COLORS.get(group, COLORS["tf"]), shape="dot" if node_type == "regulator" else "ellipse")
    for source, target, data in graph.edges(data=True):
        edge_type = data.get("edge_type", "unknown")
        network.add_edge(source, target, title=html.escape(_value({"edge_type": edge_type, "source": data.get("source", "")})), color=EDGE_COLORS.get(edge_type, "#718096"), dashes=edge_type == "co-imodulon", arrows="to")
    network.set_options("""
    {"interaction":{"hover":true,"navigationButtons":true},"physics":{"stabilization":{"iterations":300}},"edges":{"smooth":{"type":"dynamic"}}}
    """)
    # pyvis.write_html() uses the Windows process encoding (often cp1252),
    # which cannot serialize names such as Greek sigma characters. Generate
    # the document in memory and own the UTF-8 write explicitly.
    output_path.write_text(network.generate_html(notebook=False), encoding="utf-8")
    return output_path


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--graph", default=root / "network_analysis" / "output" / "regulatory_network_imodulon.pkl")
    parser.add_argument("--out", default=root / "network_analysis" / "output" / "regulatory_network.html")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv)
    with Path(args.graph).open("rb") as handle:
        graph = pickle.load(handle)
    path = render_graph(graph, args.out)
    print(f"[save] visualization -> {path}")


if __name__ == "__main__":
    main()
