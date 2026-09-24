"""Small, self-contained pipeline outputs for database integration tests."""

import csv
import json
import pickle

import networkx as nx
import pytest


def _write_csv(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


@pytest.fixture
def database_root(tmp_path):
    root = tmp_path / "project"
    config = {
        "datasets": [
            {"name": "benchmark_ceftazidime", "treatment": "ceftazidime", "output_directory": "ceftazidime", "role": "discovery"},
            {"name": "benchmark_kanamycin", "treatment": "kanamycin", "output_directory": "kanamycin", "role": "cross_reactivity_challenge"},
        ]
    }
    for relative, data in (
        ("config/benchmark.json", config),
        ("config/regulondb.lock.json", {"assets": {}}),
        ("results/differential_expression/model/benchmark_provenance.json", {}),
        ("results/regulatory_network/scoring_provenance.json", {"regulondb_release": "14.5.0"}),
    ):
        path = root / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(data))

    for observed in ("ceftazidime", "kanamycin"):
        _write_csv(
            root / f"results/differential_expression/contrasts/{observed}/de_results.csv",
            [
                {"gene_id": gene, "comparison": f"{observed}_vs_water", "baseMean": 20,
                 "log2FoldChange": effect, "lfcSE": 0.4, "pvalue": 0.001,
                 "padj": padj, "lfc_ci_low": effect - 0.8, "lfc_ci_high": effect + 0.8,
                 "regulation": "upregulated" if padj < 0.05 else "not_regulated",
                 "eligible_for_network": "True"}
                for gene, effect, padj in (("b0001", 3 if observed == "ceftazidime" else 5, 0.01 if observed == "ceftazidime" else 0.001),
                                          ("b0002", 2.5 if observed == "ceftazidime" else 0.2, 0.02 if observed == "ceftazidime" else 0.8))
            ],
        )
    _write_csv(
        root / "results/differential_expression/model/drug_pairwise_contrasts.csv",
        [{"gene_id": gene, "numerator_dataset": "benchmark_ceftazidime", "denominator_dataset": "benchmark_kanamycin",
          "baseMean": 20, "log2FoldChange": 1.2, "lfcSE": 0.3, "pvalue": 0.02, "padj": 0.04,
          "lfc_ci_low": 0.6, "lfc_ci_high": 1.8, "regulation": "not_regulated", "eligible_for_network": "True"}
         for gene in ("b0001", "b0002")],
    )
    _write_csv(
        root / "results/regulatory_network/source_observations.csv",
        [
            {"canonical_gene": "genea", "source_gene_id": gene,
             "source_dataset": dataset, "dataset_role": role,
             "log2FoldChange": effect, "padj": padj,
             "lfc_ci_low": effect - 0.8, "lfc_ci_high": effect + 0.8,
             "regulation": "upregulated"}
            for gene, dataset, role, effect, padj in (
                ("b0001", "benchmark_ceftazidime", "discovery", 3, 0.01),
                ("b0002", "benchmark_ceftazidime", "discovery", 2.5, 0.02),
                ("b0001", "benchmark_kanamycin", "cross_reactivity_challenge", 5, 0.001),
            )
        ],
    )
    graph = nx.DiGraph(network_scope="full_reference")
    graph.add_node("crp", node_type="regulator", label="CRP", is_regulator=True)
    graph.add_node("genea", node_type="candidate")
    graph.add_edge("crp", "genea", edge_type="activates", confidence="S", interaction_evidence=[{"source": "RegulonDB"}])
    path = root / "results/regulatory_network/regulatory_network.pkl"
    with path.open("wb") as handle:
        pickle.dump(graph, handle)

    _write_csv(root / "results/promoter_candidates/annotated_candidates.csv", [
        {"gene": "genea", "group": "beta_lactam", "evidence_tier": "multiple_comparisons_supported",
         "max_abs_log2_fold_change": 3, "response_profile_json": "{}", "evidence_quality_flags": "[]"}
    ])
    _write_csv(root / "results/promoter_review/promoter_review.csv", [
        {"promoter_id": "P1", "tu_ids": "T1", "operon_name": "operonA", "regulondb_release": "14.5.0",
         "promoter_name": "pA", "strand": "forward", "tss": 100, "sigma_factor": "sigma70",
         "annotated_sequence": "ACGT", "mapping_status": "curated_promoter_mapped", "candidate_genes": "genea"}
    ])
    _write_csv(root / "results/promoter_review/candidate_promoter_mapping.csv", [
        {"gene": "genea", "tu_id": "T1", "operon_name": "operonA", "tu_genes": "GeneA;GeneB",
         "promoter_id": "P1", "regulondb_release": "14.5.0"}
    ])
    assessment = root / "results/candidate_assessment"
    _write_csv(assessment / "experimental_panel.csv", [
        {"gene": "genea", "promoter_id": "P1", "antibiotic_class": "beta_lactam",
         "selection_order": 1, "selection_basis": "test", "transcription_group": "T1"}
    ])
    _write_csv(assessment / "sample_qc.csv", [
        {"sample_id": "S1", "source_file": "S1.txt", "condition": "benchmark_ceftazidime",
         "library_counts": 1000, "genes_detected": 2, "PC1": 0.1, "PC2": 0.2}
    ])
    _write_csv(assessment / "fragment_review.csv", [
        {"promoter_id": "P1", "genome_accession": "U00096.3", "fragment_start": 90,
         "fragment_end": 110, "strand": "forward", "fragment_sequence": "ACGT",
         "reference_sequence_matches": "True", "n_curated_sites": 1,
         "overlapping_genes": "genea", "construct_status": "review"}
    ])
    _write_csv(assessment / "operon_support.csv", [
        {"regulator": "crp", "source_dataset": "benchmark_ceftazidime", "direction": "upregulated",
         "overlap_genes": "genea", "dominant_group_gene_fraction": 1,
         "odds_ratio": 2, "padj": 0.04}
    ])
    _write_csv(assessment / "regulatory_sites.csv", [
        {"id": "site-1", "promoterID": "P1", "promoterName": "pA", "regulatorId": "crp",
         "regulatorName": "CRP", "tss": 100, "sigmaF": "sigma70", "tfrsLeft": 95,
         "tfrsRight": 98, "strand": "forward", "tfrsSeq": "ACGT",
         "riFunction": "activator", "tfrsEvidence": "curated", "confidenceLevel": "S"}
    ])
    return root
