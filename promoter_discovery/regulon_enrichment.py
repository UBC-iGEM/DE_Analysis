"""DE-target enrichment against full regulons, independent of candidate filters."""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

try:
    from .build_network import REGULATORY_EDGE_TYPES
except ImportError:  # pragma: no cover
    from build_network import REGULATORY_EDGE_TYPES


def bh_adjust(pvalues):
    values = np.asarray(pvalues, dtype=float)
    if not len(values):
        return values
    order = np.argsort(values, kind="stable")
    adjusted = values[order] * len(values) / np.arange(1, len(values) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1].clip(0, 1)
    restored = np.empty_like(adjusted)
    restored[order] = adjusted
    return restored


ENRICHMENT_COLUMNS = ["regulator", "source_dataset", "antibiotic_class", "direction", "n_background",
                      "n_full_regulon", "n_measured_regulon", "n_responsive_background", "n_overlap",
                      "odds_ratio", "pvalue", "padj", "overlap_genes", "activation_pattern_matches",
                      "repression_pattern_matches", "ambiguous_sign_matches"]


def compute_regulon_enrichment(graph):
    observations = pd.DataFrame(graph.graph.get("source_observations", []))
    rows = []
    if observations.empty:
        return pd.DataFrame(columns=ENRICHMENT_COLUMNS)
    regulons = {}
    for regulator, target, edge in graph.edges(data=True):
        if edge.get("edge_type") in REGULATORY_EDGE_TYPES:
            regulons.setdefault(regulator, {})[target] = edge["edge_type"]
    for dataset, group in observations.groupby("source_dataset", sort=True):
        eligible = group[group.eligible_for_network & group.identity_source.eq("same_release_mapping")]
        # A gene is one unit of evidence. Opposite significant probe estimates
        # remain visible in provenance but cannot count in either directional set.
        universe = set(eligible.canonical_gene)
        directions = eligible[eligible.is_significant].groupby("canonical_gene").regulation.agg(set)
        response = {direction: {gene for gene, signs in directions.items() if signs == {direction}}
                    for direction in ("upregulated", "downregulated")}
        for regulator, edges in sorted(regulons.items()):
            targets = set(edges) & universe
            if not targets:
                continue
            for direction, responsive in response.items():
                overlap = targets & responsive
                a, b = len(overlap), len(targets - responsive)
                c, d = len(responsive - targets), len(universe - targets - responsive)
                odds, pvalue = fisher_exact([[a, b], [c, d]], alternative="greater")
                activating_effect = "activates" if direction == "upregulated" else "represses"
                repressing_effect = "represses" if direction == "upregulated" else "activates"
                rows.append({"regulator": regulator, "source_dataset": dataset,
                             "antibiotic_class": group.antibiotic_class.iloc[0], "direction": direction,
                             "n_background": len(universe), "n_full_regulon": len(edges),
                             "n_measured_regulon": len(targets), "n_responsive_background": len(responsive),
                             "n_overlap": a, "odds_ratio": odds, "pvalue": pvalue,
                             "overlap_genes": ", ".join(sorted(overlap)),
                             "activation_pattern_matches": sum(edges[g] == activating_effect for g in overlap),
                             "repression_pattern_matches": sum(edges[g] == repressing_effect for g in overlap),
                             "ambiguous_sign_matches": sum(edges[g] == "dual" for g in overlap)})
    frame = pd.DataFrame(rows)
    if frame.empty:
        return pd.DataFrame(columns=ENRICHMENT_COLUMNS)
    # One declared family covers all regulators x drug contrasts x directions.
    frame["padj"] = bh_adjust(frame.pvalue)
    return frame[ENRICHMENT_COLUMNS].sort_values(["padj", "source_dataset", "regulator", "direction"]).reset_index(drop=True)
