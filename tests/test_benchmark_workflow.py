import json
import math
from pathlib import Path

import pandas as pd
import numpy as np
import pytest

from promoter_discovery.build_network import build_graph, load_de_results, load_regulondb
from promoter_discovery.dataset_registry import load_dataset_config
from promoter_discovery.promoter_selection import map_promoters
from promoter_discovery.regulon_enrichment import compute_regulon_enrichment
from promoter_discovery.run_benchmark import assemble_counts
from promoter_discovery.score_candidates import annotate_candidates, class_shortlist, compute_regulon_coverage
from promoter_discovery.setup_data import validate_manifest
from promoter_discovery.visualize_network import candidate_subgraph


def _fixture(tmp_path):
    config = {"thresholds": {"log2_fold_change": 2, "padj": 0.05}, "datasets": [
        {"name": "beta", "antibiotic_class": "beta_lactam", "source_study": "study", "role": "discovery"},
        {"name": "challenge", "antibiotic_class": "aminoglycoside", "source_study": "study", "role": "cross_reactivity_challenge"}]}
    path = tmp_path / "config.json"
    path.write_text(json.dumps(config))
    mapping = tmp_path / "mapping.tsv"
    pd.DataFrame({"canonical_gene": [f"gene{i}" for i in range(20)],
                  "locus_tag": [f"b{i:04}" for i in range(20)]}).to_csv(mapping, sep="\t", index=False)
    for name in ("beta", "challenge"):
        target = tmp_path / "results/differential_expression/contrasts" / name
        target.mkdir(parents=True)
        rows = [{"gene": f"b{i:04}", "gene_id": f"b{i:04}",
                 "log2FoldChange": 3 + i if name == "beta" and i < 3 else -4 if name == "challenge" and i < 3 else 6 if name == "challenge" and i == 4 else 0,
                 "padj": 0.001 if i < 3 or name == "challenge" and i == 4 else 0.9,
                 "eligible_for_network": True} for i in range(20)]
        pd.DataFrame(rows).to_csv(target / "de_results.csv", index=False)
    edges = pd.DataFrame([{"regulator": "tf", "target": f"gene{i}", "edge_type": "activates", "regulator_type": "TF"} for i in range(3)])
    return path, mapping, edges


def test_full_network_enrichment_is_invariant_to_candidate_filter(tmp_path):
    path, mapping, edges = _fixture(tmp_path)
    config = load_dataset_config(path)
    all_obs = load_de_results(path, tmp_path, mapping)
    top_obs = load_de_results(path, tmp_path, mapping, top_n=1)
    graph = build_graph(all_obs, edges, config)
    limited = build_graph(top_obs, edges, config)
    assert limited.nodes["gene0"]["node_type"] == "reference_gene"
    assert limited.number_of_edges() == 3
    assert candidate_subgraph(limited).number_of_edges() == 1
    enrichment = compute_regulon_enrichment(graph)
    pd.testing.assert_frame_equal(enrichment, compute_regulon_enrichment(limited))
    beta = enrichment[(enrichment.source_dataset == "beta") & (enrichment.direction == "upregulated")].iloc[0]
    assert beta.n_background == 20  # includes measured genes without any reference edge
    assert beta.n_overlap == beta.n_measured_regulon == 3
    assert beta.pvalue == pytest.approx(1 / math.comb(20, 3))
    assert beta.padj == pytest.approx(2 / math.comb(20, 3))  # four tests, two equally small
    assert beta.activation_pattern_matches == 3


def test_full_regulon_coverage_does_not_call_one_candidate_specific(tmp_path):
    path, mapping, edges = _fixture(tmp_path)
    edges = pd.concat([edges, pd.DataFrame([{"regulator": "tf", "target": f"outside{i}", "edge_type": "represses", "regulator_type": "TF"} for i in range(97)])], ignore_index=True)
    graph = build_graph(load_de_results(path, tmp_path, mapping, top_n=1), edges, load_dataset_config(path))
    coverage = compute_regulon_coverage(graph).iloc[0]
    assert coverage.n_total_targets == 100
    assert coverage.candidate_coverage_beta_lactam == 0.01
    scored = annotate_candidates(graph)
    assert "regulatory_clarity" not in scored
    assert scored.biosensor_flag.eq("requires_promoter_mapping_and_wetlab_validation").all()


def test_cross_class_directions_are_profiles_not_conflicts_and_challenges_do_not_seed(tmp_path):
    path, mapping, edges = _fixture(tmp_path)
    graph = build_graph(load_de_results(path, tmp_path, mapping), edges, load_dataset_config(path))
    assert "gene4" not in graph  # challenge-only induction
    assert graph.nodes["gene0"]["cross_class_direction_difference"]
    assert graph.nodes["gene0"]["evidence_tier"] == "single_comparison_supported"
    scored = annotate_candidates(graph)
    assert scored.response_profile_status.eq("responses_in_multiple_classes").all()
    baseline = scored.gene.tolist()
    for gene in baseline:
        graph.nodes[gene].update(metabolic_burden_proxy=1000, gene_expression={"unmatched": 100}, imodulon_activity={"unmatched": -999})
    assert annotate_candidates(graph).gene.tolist() == baseline


def test_unavailable_and_unmapped_genes_are_excluded_from_enrichment(tmp_path):
    path, mapping, edges = _fixture(tmp_path)
    obs = load_de_results(path, tmp_path, mapping)
    obs.loc[obs.canonical_gene.eq("gene19"), "eligible_for_network"] = False
    obs.loc[obs.canonical_gene.eq("gene18"), "identity_source"] = "source_gene"
    result = compute_regulon_enrichment(build_graph(obs, edges, load_dataset_config(path)))
    assert result.n_background.eq(18).all()


def test_opposing_regulatory_records_are_preserved_as_dual(tmp_path):
    path = tmp_path / "network.tsv"
    path.write_text("R1\tTF\ttf\tG1\tgene\t+\tS\nR1\tTF\ttf\tG1\tgene\t-\tW\n")
    frame = load_regulondb(path)
    assert len(frame) == 1
    assert frame.iloc[0].edge_type == "dual"
    assert len(frame.iloc[0].interaction_evidence) == 2


def test_promoter_mapping_groups_operon_genes_without_inventing_fragments():
    candidates = pd.DataFrame({"gene": ["a", "b", "c", "missing"]})
    tus = pd.DataFrame([{"id": "tu", "tuGenes": "a;b", "promoterId": "p", "operonName": "ab"},
                        {"id": "tu2", "tuGenes": "c", "promoterId": "", "operonName": "c"}])
    promoters = pd.DataFrame([{"id": "p", "name": "abp", "strand": "reverse", "posTSS": "100",
                              "sequence": "acgtA", "confidenceLevel": "S"}])
    mapped, review = map_promoters(candidates, tus, promoters, "14.5.0")
    assert len(review) == 1
    assert review.iloc[0].n_candidate_genes == 2
    assert review.iloc[0].candidate_genes == "a, b"
    assert "fragment_start" not in review
    assert mapped.set_index("gene").loc["missing", "mapping_status"] == "no_curated_transcription_unit"
    assert mapped.set_index("gene").loc["c", "mapping_status"] == "tu_without_promoter_annotation"


def test_count_assembly_deduplicates_controls_and_rejects_fractional_input(monkeypatch):
    from promoter_discovery import run_benchmark
    def loader(dataset):
        ids = ["ctrl1", "ctrl2", "ctrl3", *(dataset["name"] + str(i) for i in range(3))]
        expression = pd.DataFrame({"gene_id": ["g1", "g2"], **{sample: [10, 20] for sample in ids}})
        metadata = pd.DataFrame({"sample_id": ids, "condition": ["control"] * 3 + ["treated"] * 3, "source_file": ids})
        return expression, metadata
    monkeypatch.setattr(run_benchmark, "load_tar_gene_tables", loader)
    counts, metadata = assemble_counts({"datasets": [{"name": "a"}, {"name": "b"}]})
    assert len(metadata) == 9
    assert counts.shape == (2, 9)
    def fractional(dataset):
        expression, metadata = loader(dataset)
        expression["ctrl1"] = [0.5, 2]
        return expression, metadata
    monkeypatch.setattr(run_benchmark, "load_tar_gene_tables", fractional)
    with pytest.raises(ValueError, match="raw integer counts"):
        assemble_counts({"datasets": [{"name": "a"}]})


def test_primary_reference_manifest_does_not_require_precise():
    path = Path(__file__).resolve().parents[1] / "config/regulondb_assets.json"
    manifest = json.loads(path.read_text())
    assert not validate_manifest(manifest, check_files=False)


def test_each_class_shortlist_uses_its_own_effect_and_support():
    import networkx as nx
    graph = nx.DiGraph()
    for gene, tier, amino_effect, beta_effect in [("a", "multiple_comparisons_supported", 3, 12), ("b", "single_comparison_supported", 7, 3)]:
        graph.add_node(gene, node_type="candidate", significant_classes=["beta_lactam", "aminoglycoside"],
                       evidence_tier=tier, max_abs_log2_fold_change=beta_effect,
                       evidence_tier_by_class={"beta_lactam": tier, "aminoglycoside": "single_comparison_supported"},
                       dataset_evidence_json=json.dumps([{"antibiotic_class": cls, "qualifying": True,
                                                          "dataset_role": "discovery", "log2FoldChange": [effect]}
                                                         for cls, effect in [("aminoglycoside", amino_effect), ("beta_lactam", beta_effect)]]))
    candidates = annotate_candidates(graph)
    assert class_shortlist(candidates, graph, "aminoglycoside").gene.tolist() == ["b", "a"]
    assert class_shortlist(candidates, graph, "beta_lactam").gene.tolist() == ["a", "b"]


def test_benchmark_config_is_one_study_with_explicit_proxy_roles():
    path = Path(__file__).resolve().parents[1] / "config/benchmark.json"
    config = load_dataset_config(path)
    assert {d["source_study"] for d in config["datasets"]} == {"GSE220559"}
    assert {d["de_method"] for d in config["datasets"]} == {"pydeseq2"}
    assert set(config["benchmark"]["wetlab_antibiotics"]) == {"amoxicillin", "cephalexin", "gentamicin", "tobramycin"}
    assert sum(d["role"] == "cross_reactivity_challenge" for d in config["datasets"]) == 2


def test_joint_count_model_exports_matched_and_direct_evidence_with_provenance(tmp_path, monkeypatch):
    from promoter_discovery import run_benchmark
    rng = np.random.default_rng(42)
    genes = [f"b{i:04}" for i in range(150)]
    ids = [f"ctrl{i}" for i in range(3)] + [f"beta{i}" for i in range(3)] + [f"amino{i}" for i in range(3)]
    means = np.full((150, 9), 80.0)
    means[:10, 3:6] *= 12
    means[10:20, 6:9] *= 12
    master = pd.DataFrame(rng.negative_binomial(40, 40 / (40 + means)), index=genes, columns=ids)
    def loader(dataset):
        selected = ids[:3] + [dataset["name"] + str(i) for i in range(3)]
        expression = master[selected].reset_index(names="gene_id")
        metadata = pd.DataFrame({"sample_id": selected, "condition": ["control"] * 3 + ["treated"] * 3, "source_file": selected})
        return expression, metadata
    monkeypatch.setattr(run_benchmark, "load_tar_gene_tables", loader)
    archive = tmp_path / "counts-source.tsv"
    master.to_csv(archive, sep="\t")
    config = {"thresholds": {"log2_fold_change": 2.0, "padj": 0.05},
              "benchmark": {"min_count": 10, "min_samples": 3},
              "datasets": [{"name": name, "antibiotic_class": cls, "de_method": "pydeseq2",
                            "source_study": "synthetic", "role": "discovery", "archive": str(archive),
                            "treatment": name, "comparison": name + "_vs_water"}
                           for name, cls in [("beta", "beta_lactam"), ("amino", "aminoglycoside")]]}
    path = tmp_path / "config.json"
    path.write_text(json.dumps(config))
    out = tmp_path / "results/differential_expression/model"
    run_benchmark.run(path, tmp_path, out)
    provenance = json.loads((out / "benchmark_provenance.json").read_text())
    assert provenance["n_samples"] == 9
    assert provenance["n_genes_retained"] == 150
    result = pd.read_csv(tmp_path / "results/differential_expression/contrasts/beta/de_results.csv").set_index("gene_id")
    assert result.loc[genes[:10], "regulation"].eq("upregulated").sum() >= 8
    assert result.loc[genes[:10], "log2FoldChange"].median() > 3
    assert result.lfcSE.gt(0).all()
    assert result.shrunk_log2FoldChange.notna().all()
    assert result.shrunk_lfcSE.gt(0).all()
    assert (out / "joint_vst_counts.csv").exists()
    assert pd.read_csv(out / "joint_raw_counts.csv", index_col=0).shape == master.shape
    assert result.lfc_ci_low.le(result.log2FoldChange).all()
    pairs = pd.read_csv(out / "drug_pairwise_contrasts.csv").set_index("gene_id")
    assert pairs.loc[genes[:10], "log2FoldChange"].median() > 3
    assert pairs.numerator_dataset.eq("beta").all()
    observations = load_de_results(path, tmp_path)
    assert len(observations) == 300
    assert observations.attrs["benchmark_provenance"]["pairwise_sha256"] == provenance["pairwise_sha256"]
    export = tmp_path / "results/differential_expression/contrasts/beta/de_results.csv"
    export.write_text(export.read_text() + "\n")
    with pytest.raises(ValueError, match="not verified by joint-fit provenance"):
        load_de_results(path, tmp_path)
