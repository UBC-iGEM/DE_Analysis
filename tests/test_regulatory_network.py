import json
import ssl
from pathlib import Path
from types import SimpleNamespace

import networkx as nx
import pandas as pd
import pytest

from promoter_discovery.build_network import aggregate_candidate_evidence, build_graph, load_de_results, load_identity_mapping, save_graph
from promoter_discovery.dataset_registry import dataset_caveats, load_dataset_config, version_at_least
from promoter_discovery.score_candidates import annotate_candidates, compute_regulon_coverage
from promoter_discovery import setup_data, score_candidates, i_modulon_analysis
from promoter_discovery.setup_data import normalize_gene_product_mapping, validate_manifest
from promoter_discovery.i_modulon_analysis import compute_gene_expression_evidence, load_embedded_expression, validate_expression_alignment
from promoter_discovery.visualize_network import _edge_arrows


def test_dual_role_candidates_survive_annotation_and_default_scoring(tmp_path, monkeypatch):
    config_path = _config(tmp_path)
    for dataset in ("amoxicillin", "gentamicin", "tobramycin"):
        _write_de(tmp_path, dataset, [
            {"gene": gene, "gene_id": locus, "log2FoldChange": 3.0, "padj": 0.01}
            for gene, locus in (("geneA", "b0001"), ("geneB", "b0002"))
        ])
    config = load_dataset_config(config_path)
    observations = load_de_results(config_path, tmp_path)
    interactions = pd.DataFrame([
        {"regulator": "genea", "target": target, "edge_type": "activates", "regulator_type": "TF"}
        for target in ("genea", "geneb")
    ])
    graph = build_graph(observations, interactions, config)
    assert graph.nodes["genea"]["node_type"] == "candidate"
    assert graph.nodes["genea"]["is_regulator"] is True
    assert graph.nodes["genea"]["group"] == "cross"
    assert graph.nodes["genea"]["n_datasets_significant"] == 3
    output = tmp_path / "results" / "regulatory_network"
    save_graph(graph, output)

    model = SimpleNamespace(
        M=pd.DataFrame({"stress": [1.0, 2.0]}, index=["b0001", "b0002"]),
        A=pd.DataFrame({"ctrl": [0.0], "drug": [5.0]}, index=["stress"]),
        gene_table=pd.DataFrame({"gene_name": ["genea", "geneb"]}, index=["b0001", "b0002"]),
        imodulon_table=pd.DataFrame({"system_category": ["Metabolism"]}, index=["stress"]),
    )
    groups = {key: {"backgrounds": []} for key in i_modulon_analysis.CLASS_REGISTRY}
    groups["beta_lactam"]["backgrounds"] = [
        {"background": "MH", "control": ["ctrl"], "treated": ["drug"]}
    ]
    expression = pd.DataFrame({"ctrl": [1.0, 2.0], "drug": [4.0, 6.0]}, index=model.M.index)
    graph, candidates, discovered, _ = i_modulon_analysis.annotate_graph(
        graph, model, i_modulon_analysis.build_gene_lookup(model), model.M.gt(0), {},
        groups, expression, "log2(TPM)", "fixture",
    )
    assert set(candidates["gene"]) == {"genea", "geneb"}
    assert graph.edges["genea", "geneb"]["edge_type"] == "activates"
    annotated_path = output / "regulatory_network_imodulon.pkl"
    i_modulon_analysis.save_outputs(graph, candidates, pd.DataFrame(), output, annotated_path, discovered, [])
    # Exercise the CLI defaults against both raw and annotated graph files.
    monkeypatch.setattr(score_candidates, "__file__", str(tmp_path / "promoter_discovery" / "score_candidates.py"))
    assert score_candidates.parse_args([]).graph == output / "regulatory_network.pkl"
    score_candidates.main(["--graph", str(annotated_path)])
    scored = pd.read_csv(tmp_path / "results/promoter_candidates/annotated_candidates.csv").set_index("gene")
    assert set(scored.index) == {"genea", "geneb"}
    assert scored.loc["genea", "regulators"] == "genea"
    assert scored.loc["geneb", "regulators"] == "genea"
    assert json.loads(scored.loc["genea", "gene_expression_delta"])["beta_lactam"] == 3.0
    assert scored.loc["genea", "metabolic_burden_proxy"] == 1.0
    assert json.loads(scored.loc["geneb", "imodulon_activity_delta"])["beta_lactam"] == 5.0
    regulators = pd.read_csv(output / "regulon_candidate_coverage.csv").set_index("regulator")
    assert regulators.loc["genea", "n_total_targets"] == 2


def test_ssl_context_verifies_certificates_and_hostnames(monkeypatch: pytest.MonkeyPatch):
    monkeypatch.delenv("SSL_CERT_FILE", raising=False)
    context = setup_data._ssl_context()
    assert context.verify_mode == ssl.CERT_REQUIRED
    assert context.check_hostname is True


def test_ssl_context_honors_configured_ca_bundle(monkeypatch: pytest.MonkeyPatch):
    calls: list[str | None] = []

    def create_default_context(*, cafile: str | None = None) -> ssl.SSLContext:
        calls.append(cafile)
        return ssl.SSLContext(ssl.PROTOCOL_TLS_CLIENT)

    monkeypatch.setenv("SSL_CERT_FILE", "/path/to/project-ca.pem")
    monkeypatch.setattr(setup_data.ssl, "create_default_context", create_default_context)
    setup_data._ssl_context()
    assert calls == ["/path/to/project-ca.pem"]


def _write_de(root: Path, dataset: str, rows: list[dict]) -> None:
    path = root / "results/differential_expression/contrasts" / dataset
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
    assert aggregate_candidate_evidence(up, load_dataset_config(config_path)).iloc[0]["candidate_direction_policy"] == "upregulated"
    assert aggregate_candidate_evidence(either, load_dataset_config(config_path)).iloc[0]["candidate_direction_policy"] == "either-direction"
    evidence = aggregate_candidate_evidence(up, load_dataset_config(config_path))
    assert evidence.set_index("canonical_gene").loc["genea", "evidence_tier"] == "single_comparison_supported"
    assert evidence.set_index("canonical_gene").loc["geneb", "evidence_tier"] == "single_comparison_supported"
    assert "source_row_ids" in evidence.columns


def test_candidate_normalization_uses_configured_thresholds(tmp_path: Path):
    config = {
        "thresholds": {"log2_fold_change": 2.0, "padj": 0.1},
        "datasets": [
            {"name": "amoxicillin", "antibiotic_class": "beta_lactam"},
        ],
    }
    config_path = tmp_path / "config.json"
    config_path.write_text(json.dumps(config), encoding="utf-8")
    _write_de(
        tmp_path,
        "amoxicillin",
        [
            {
                "gene": "geneA",
                "gene_id": "b0001",
                "log2FoldChange": 3.0,
                "padj": 0.08,
                "regulation": "not_regulated",
            },
            {
                "gene": "geneB",
                "gene_id": "b0002",
                "log2FoldChange": 2.0,
                "padj": 0.01,
                "regulation": "not_regulated",
            },
        ],
    )

    observations = load_de_results(config_path, tmp_path)
    indexed = observations.set_index("canonical_gene")
    assert indexed.loc["genea", "regulation"] == "upregulated"
    assert bool(indexed.loc["genea", "candidate_seed"]) is True
    assert indexed.loc["geneb", "regulation"] == "not_regulated"
    assert bool(indexed.loc["geneb", "candidate_seed"]) is False


def test_top_n_uses_strongest_signal_per_gene_and_class(tmp_path: Path):
    config = {
        "thresholds": {"log2_fold_change": 2.0, "padj": 0.05},
        "datasets": [
            {"name": "beta", "antibiotic_class": "beta_lactam"},
            {"name": "amino", "antibiotic_class": "aminoglycoside"},
        ],
    }
    config_path = tmp_path / "config.json"
    config_path.write_text(json.dumps(config), encoding="utf-8")
    _write_de(
        tmp_path,
        "beta",
        [
            {"gene": "geneA", "gene_id": "b0001", "log2FoldChange": 3.0, "padj": 0.01, "regulation": "upregulated"},
            {"gene": "geneB", "gene_id": "b0002", "log2FoldChange": 5.0, "padj": 0.01, "regulation": "upregulated"},
        ],
    )
    _write_de(
        tmp_path,
        "amino",
        [
            {"gene": "geneA", "gene_id": "b0001", "log2FoldChange": 10.0, "padj": 0.01, "regulation": "upregulated"},
            {"gene": "geneC", "gene_id": "b0003", "log2FoldChange": 9.0, "padj": 0.01, "regulation": "upregulated"},
        ],
    )

    observations = load_de_results(config_path, tmp_path, top_n=1)
    selected = set(observations.loc[observations["candidate_seed"], "canonical_gene"])
    assert selected == {"genea", "geneb"}


def test_co_imodulon_edges_do_not_score_as_regulators():
    graph = nx.DiGraph()
    graph.add_node("rpoe", node_type="regulator", group="tf", regulator_type="sigma")
    graph.add_node("genea", node_type="candidate", group="beta_lactam", significant_classes=["beta_lactam"])
    graph.add_node("geneb", node_type="candidate", group="aminoglycoside", significant_classes=["aminoglycoside"])
    graph.add_edge("rpoe", "genea", edge_type="activates")
    graph.add_edge("genea", "geneb", edge_type="co-imodulon")
    scores = compute_regulon_coverage(graph)
    assert scores.iloc[0]["n_total_targets"] == 1
    annotated = annotate_candidates(graph)
    assert annotated.set_index("gene").loc["geneb", "regulators"] == ""


def test_candidate_scores_rank_strong_evidence_before_conflicts():
    graph = nx.DiGraph()
    for tier in ("single_comparison_supported", "within_class_direction_difference", "multiple_comparisons_supported"):
        graph.add_node(
            tier,
            node_type="candidate",
            group="beta_lactam",
            significant_classes=["beta_lactam"],
            evidence_tier=tier,
            max_abs_log2_fold_change=3.0,
        )

    annotated = annotate_candidates(graph)
    assert annotated["gene"].tolist() == [
        "multiple_comparisons_supported",
        "single_comparison_supported",
        "within_class_direction_difference",
    ]


def test_only_regulatory_edges_render_with_arrows():
    assert _edge_arrows("activates") == "to"
    assert _edge_arrows("represses") == "to"
    assert _edge_arrows("dual") == "to"
    assert _edge_arrows("co-imodulon") == ""


def test_extended_antibiotic_classes_are_supported_and_scored(tmp_path: Path):
    config = {
        "thresholds": {"log2_fold_change": 2.0, "padj": 0.05},
        "datasets": [
            {"name": "ciprofloxacin", "antibiotic_class": "fluoroquinolone"},
            {"name": "polymixinE", "antibiotic_class": "polymyxin"},
        ],
    }
    config_path = tmp_path / "config.json"
    config_path.write_text(json.dumps(config), encoding="utf-8")
    assert load_dataset_config(config_path)["datasets"][0]["antibiotic_class"] == "fluoroquinolone"

    graph = nx.DiGraph()
    graph.add_node("rpoe", node_type="regulator", group="tf", regulator_type="sigma")
    graph.add_node("gyrA", node_type="candidate", group="fluoroquinolone", significant_classes=["fluoroquinolone"])
    graph.add_node("pmrA", node_type="candidate", group="polymyxin", significant_classes=["polymyxin"])
    graph.add_edge("rpoe", "gyrA", edge_type="activates")
    graph.add_edge("rpoe", "pmrA", edge_type="activates")
    scores = compute_regulon_coverage(graph)
    score = scores.iloc[0]
    assert score["n_fluoroquinolone_targets"] == 1
    assert score["n_polymyxin_targets"] == 1
    assert bool(score["has_candidates_in_multiple_classes"])


def test_dataset_config_rejects_duplicate_names(tmp_path: Path):
    path = tmp_path / "config.json"
    path.write_text(json.dumps({
        "thresholds": {"log2_fold_change": 2.0, "padj": 0.05},
        "datasets": [
            {"name": "gentamicin", "antibiotic_class": "aminoglycoside"},
            {"name": "gentamicin", "antibiotic_class": "aminoglycoside"},
        ],
    }), encoding="utf-8")
    with pytest.raises(ValueError, match="Duplicate dataset name"):
        load_dataset_config(path)


def test_current_amoxicillin_comparisons_have_provenance_caveats():
    resistant = dataset_caveats("amoxicillin_resistant_vs_wt")
    exposed = dataset_caveats("amoxicillin_resistant_amox_vs_wt_amox")
    assert resistant and "resistance background" in resistant[0]
    assert exposed and "not an antibiotic-versus-control contrast" in exposed[0]


def test_expression_contract_retains_basal_induced_and_backgrounds():
    expression = pd.DataFrame({"ctrl_1": [1.0], "drug_1": [8.0]}, index=["b0001"])
    groups = {"beta_lactam": {"backgrounds": [{"background": "cam", "control": ["ctrl_1"], "treated": ["drug_1"]}]}, "aminoglycoside": {"backgrounds": []}}
    evidence = compute_gene_expression_evidence(expression, "b0001", groups, "TPM", "log2")
    beta = evidence["beta_lactam"]
    assert (beta["basal"], beta["induced"], beta["delta"]) == (1.0, 8.0, 7.0)
    assert beta["available"] is True
    assert beta["per_background"][0]["n_control"] == 1
    assert beta["units"] == "TPM"


def test_embedded_expression_precedes_external_fallback_and_aligns_to_ica():
    ica = SimpleNamespace(
        A=pd.DataFrame([[1.0]], index=["imod1"], columns=["sample1"]),
        M=pd.DataFrame([[0.5]], index=["b0001"], columns=["imod1"]),
        log_tpm=pd.DataFrame([[4.0]], index=["b0001"], columns=["sample1"]),
        X=pd.DataFrame([[99.0]], index=["b0001"], columns=["sample1"]),
    )
    expression, metadata = load_embedded_expression(ica)
    assert metadata["source"] == "IcaData.log_tpm"
    assert metadata["units"] == "log2(TPM)"
    assert expression.loc["b0001", "sample1"] == 4.0
    validate_expression_alignment(expression, ica)

    ica.A.columns = ["unmatched"]
    with pytest.raises(ValueError, match="no sample IDs"):
        validate_expression_alignment(expression, ica)


def test_manifest_versions_and_matching_release_without_hand_entered_hashes(tmp_path: Path):
    regulon = tmp_path / "regulon.tsv"
    regulon.write_text("# RegulonDB Release: 14.5\n", encoding="utf-8")
    model = tmp_path / "model.json.gz"
    expression = tmp_path / "expression.csv"
    model.write_bytes(b"model")
    expression.write_bytes(b"expression")
    assets = []
    for name, kind, path, release, revision in (
        ("reg", "regulondb", regulon, "14.5.0", "RegulonDB-14.5.0"),
        ("model", "imodulon", model, "1.0", "v1.0"),
        ("expr", "expression", expression, "1.0", "v1.0"),
    ):
        assets.append({"name": name, "kind": kind, "provider": "file", "url": path.as_uri(), "path": path.name, "dataset_release": release, "source_revision": revision, "license": "test", "citation": "test", "asset_id": "rdb" if kind == "regulondb" else "imod"})
    manifest = {"minimum_versions": {"regulondb": "14.5.0", "precise1k": "1.0"}, "assets": assets}
    assert not validate_manifest(manifest, tmp_path)
    assets[-1]["source_revision"] = "v2.0"
    assert any("share dataset_release" in error for error in validate_manifest(manifest, tmp_path))
    assert version_at_least("14.5", "14.5.0")


def test_normalize_regulondb_gene_product_mapping(tmp_path: Path):
    raw = tmp_path / "GeneProductAllIdentifiersSet.tsv"
    raw.write_text(
        "1)geneId\t2)geneName\t6)geneSynonyms\t7)otherDbsGeneIds\n"
        "RDB1\talr\talr5,EG10001\t[REFSEQ:b4053][OU-MICROARRAY:b4053]\n",
        encoding="utf-8",
    )
    output = tmp_path / "mapping.tsv"
    result = normalize_gene_product_mapping(raw, output, "14.5.0")
    frame = pd.read_csv(output, sep="\t")
    assert result["rows"] == 1
    assert frame.iloc[0]["canonical_gene"] == "alr"
    assert frame.iloc[0]["canonical_locus_tag"] == "b4053"
    assert "b4053" in frame.iloc[0]["aliases"]
    lookup = load_identity_mapping(output)
    assert lookup["alr5"] == ("alr", "b4053")


def test_graphql_download_writes_raw_products_lock_and_mapping(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    responses = {
        "NetworkRegulatorGene": {"_id": "r1", "fileName": "NetworkRegulatorGene", "content": "1)regulatorId\t2)regulatorName\t3)RegulatorGeneName\t4)regulatedGeneId\t5)regulatedGeneName\t6)function\t7)confidenceLevel\nR1\tRpoD\trpoD\tG1\talr\t+\tW\n", "columnsDetails": "# columns"},
        "NetworkSigmaGene": {"_id": "s1", "fileName": "NetworkSigmaGene", "content": "1)sigmaName\t2)regulatedGeneName\t3)function\t4)promoterEvidence\t5)confidenceLevel\nsigma70\talr\t+\tE\tW\n", "columnsDetails": "# columns"},
        "GeneProductAllIdentifiersSet": {"_id": "m1", "fileName": "GeneProductAllIdentifiersSet", "content": "1)geneId\t2)geneName\t6)geneSynonyms\t7)otherDbsGeneIds\nG1\talr\talr5\t[REFSEQ:b4053]\n", "columnsDetails": "# columns"},
    }

    def fake_graphql(endpoint: str, query: str, variables: dict | None = None) -> dict:
        if "getDatabaseInfo" in query:
            return {"getDatabaseInfo": [{"regulonDBVersion": "14.5.0", "releaseDate": "2026-01-28"}], "listAllFileNames": sorted(responses)}
        return {"getDataOfFile": responses[variables["fileName"]]}

    monkeypatch.setattr(setup_data, "_graphql_request", fake_graphql)
    manifest = {
        "minimum_versions": {"regulondb": "14.5.0", "precise1k": "1.0"},
        "regulondb": {"endpoint": "https://example.test/graphql"},
        "assets": [
            {"name": "reg", "kind": "regulondb", "provider": "regulondb_graphql", "remote_name": "NetworkRegulatorGene", "path": "data/reg.tsv"},
            {"name": "sigma", "kind": "regulondb", "provider": "regulondb_graphql", "remote_name": "NetworkSigmaGene", "path": "data/sigma.tsv"},
            {"name": "mapping", "kind": "regulondb_mapping", "provider": "regulondb_graphql", "remote_name": "GeneProductAllIdentifiersSet", "path": "data/raw_mapping.tsv", "derived_path": "data/references/regulondb/network_gene_mapping.tsv"},
        ],
    }
    lock_path = tmp_path / "precise1k.lock.json"
    (tmp_path / "data").mkdir()
    (tmp_path / "data/model.json.gz").write_bytes(b"model")
    (tmp_path / "data/expression.csv").write_bytes(b"gene\tctrl\n")
    manifest["assets"].extend([
        {"name": "model", "kind": "imodulon", "provider": "file", "url": (tmp_path / "data/model.json.gz").as_uri(), "path": "data/model.json.gz", "dataset_release": "1.0", "source_revision": "v1.0"},
        {"name": "expression", "kind": "expression", "provider": "file", "url": (tmp_path / "data/expression.csv").as_uri(), "path": "data/expression.csv", "dataset_release": "1.0", "source_revision": "v1.0"},
    ])
    lock = setup_data.download_assets(manifest, tmp_path, lock_path)
    assert lock["assets"]["reg"]["remote_id"] == "r1"
    assert lock["assets"]["reg"]["source_url"] == "https://example.test/graphql"
    assert (tmp_path / "data/references/regulondb/network_gene_mapping.tsv").exists()
    assert "b4053" in (tmp_path / "data/references/regulondb/network_gene_mapping.tsv").read_text(encoding="utf-8")
    assert lock["assets"]["mapping:derived"]["path"] == "data/references/regulondb/network_gene_mapping.tsv"
    assert not validate_manifest(manifest, tmp_path, lock=lock, require_lock=True)
    cached_lock = setup_data.download_assets(manifest, tmp_path, lock_path)
    assert cached_lock["assets"]["reg"]["remote_id"] == "r1"
    (tmp_path / "data/references/regulondb/network_gene_mapping.tsv").write_text("modified\n", encoding="utf-8")
    assert any(
        "SHA256 mismatch against lock for mapping:derived" in error
        for error in validate_manifest(manifest, tmp_path, lock=cached_lock, require_lock=True)
    )
