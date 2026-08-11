import hashlib
import json
from pathlib import Path

import networkx as nx
import pandas as pd
import pytest

from network_analysis.build_network import aggregate_candidate_evidence, build_graph, load_de_results
from network_analysis.dataset_registry import load_dataset_config, version_at_least
from network_analysis.score_candidates import annotate_candidates, compute_tf_specificity
from network_analysis.setup_data import load_manifest, sha256_file, validate_manifest
from network_analysis.i_modulon_analysis import compute_gene_expression_evidence


def _write_de(root: Path, dataset: str, rows: list[dict]) -> None:
    path = root / "data" / dataset / "standardized"
    path.mkdir(parents=True)
    pd.DataFrame(rows).to_csv(path / "de_results.csv", index=False)


def _config(root: Path) -> Path:
    config = {
        "thresholds": {"log2_fold_change": 2.0, "padj": 0.05},
        "datasets": [
            {"name": "amoxicillin", "antibiotic_class": "beta_lactam"},
            {"name": "gentamicin", "antibiotic_class": "aminoglycoside"},
            {"name": "tobramycin", "antibiotic_class": "aminoglycoside"},
        ],
    }
    path = root / "config.json"
    path.write_text(json.dumps(config), encoding="utf-8")
    return path


def test_direction_modes_preserve_conflicts_and_tobramycin_tier(tmp_path: Path):
    config_path = _config(tmp_path)
    _write_de(tmp_path, "amoxicillin", [{"gene": "geneA", "gene_id": "b0001", "log2FoldChange": 3.0, "padj": 0.01, "regulation": "upregulated"}])
    _write_de(tmp_path, "gentamicin", [{"gene": "geneA", "gene_id": "b0001", "log2FoldChange": -3.0, "padj": 0.01, "regulation": "downregulated"}])
    _write_de(tmp_path, "tobramycin", [{"gene": "geneB", "gene_id": "b0002", "log2FoldChange": 3.0, "padj": 0.01, "regulation": "upregulated"}])
    mapping = tmp_path / "mapping.tsv"
    pd.DataFrame({"canonical_gene": ["geneA", "geneB"], "locus_tag": ["b0001", "b0002"]}).to_csv(mapping, sep="\t", index=False)

    up = load_de_results(config_path, tmp_path, mapping, "upregulated")
    either = load_de_results(config_path, tmp_path, mapping, "either-direction")
    assert set(up.loc[up.candidate_seed, "canonical_gene"]) == {"genea", "geneb"}
    assert set(either.loc[either.candidate_seed, "canonical_gene"]) == {"genea", "geneb"}
    evidence = aggregate_candidate_evidence(up, load_dataset_config(config_path))
    assert evidence.set_index("canonical_gene").loc["genea", "evidence_tier"] == "conflicted"
    assert evidence.set_index("canonical_gene").loc["geneb", "evidence_tier"] == "limited"
    assert "source_row_ids" in evidence.columns


def test_co_imodulon_edges_do_not_score_as_regulators():
    graph = nx.DiGraph()
    graph.add_node("rpoe", node_type="regulator", group="tf", regulator_type="sigma")
    graph.add_node("genea", node_type="candidate", group="beta_lactam", significant_classes=["beta_lactam"])
    graph.add_node("geneb", node_type="candidate", group="aminoglycoside", significant_classes=["aminoglycoside"])
    graph.add_edge("rpoe", "genea", edge_type="activates")
    graph.add_edge("genea", "geneb", edge_type="co-imodulon")
    scores = compute_tf_specificity(graph)
    assert scores.iloc[0]["n_total_targets"] == 1
    annotated = annotate_candidates(graph, scores)
    assert annotated.set_index("gene").loc["geneb", "regulators"] == ""


def test_expression_contract_retains_basal_induced_and_backgrounds():
    expression = pd.DataFrame({"ctrl_1": [1.0], "drug_1": [8.0]}, index=["b0001"])
    groups = {"beta_lactam": {"backgrounds": [{"background": "cam", "control": ["ctrl_1"], "treated": ["drug_1"]}]}, "aminoglycoside": {"backgrounds": []}}
    evidence = compute_gene_expression_evidence(expression, "b0001", groups, "TPM", "log2")
    beta = evidence["beta_lactam"]
    assert (beta["basal"], beta["induced"], beta["delta"]) == (1.0, 8.0, 7.0)
    assert beta["available"] is True
    assert beta["per_background"][0]["n_control"] == 1
    assert beta["units"] == "TPM"


def test_manifest_checksums_versions_and_matching_release(tmp_path: Path):
    regulon = tmp_path / "regulon.tsv"
    regulon.write_text("# RegulonDB Release: 14.5\n", encoding="utf-8")
    model = tmp_path / "model.json.gz"
    expression = tmp_path / "expression.csv"
    model.write_bytes(b"model")
    expression.write_bytes(b"expression")
    digest = lambda path: sha256_file(path)
    assets = []
    for name, kind, path in (("reg", "regulondb", regulon), ("model", "imodulon", model), ("expr", "expression", expression)):
        assets.append({"name": name, "kind": kind, "url": path.as_uri(), "path": path.name, "version": "14.5.0" if kind == "regulondb" else "2.5.0", "release": "14.5.0" if kind == "regulondb" else "2.5.0", "official_revision": "r1", "retrieved_at": "2026-08-11", "sha256": digest(path), "license": "test", "citation": "test", "asset_id": "rdb" if kind == "regulondb" else "imod"})
    manifest = {"assets": assets}
    assert not validate_manifest(manifest, tmp_path)
    assets[-1]["asset_id"] = "different"
    assert any("same release/version" in error for error in validate_manifest(manifest, tmp_path))
    assert version_at_least("14.5", "14.5.0")
