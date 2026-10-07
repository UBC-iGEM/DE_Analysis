import numpy as np
import pandas as pd
import pytest

from promoter_discovery.assess_candidates import (
    candidate_stability, cross_reactivity, fragment_review,
    operon_support, oriented_sequence, select_panel, transcription_groups,
)


def test_direct_comparison_reverses_intervals_and_keeps_null_inconclusive():
    config = {"thresholds": {"padj": 0.05, "log2_fold_change": 2}, "datasets": [
        {"name": "a", "antibiotic_class": "beta_lactam", "role": "discovery"},
        {"name": "b", "antibiotic_class": "other", "role": "cross_reactivity_challenge"}]}
    observations = pd.DataFrame([
        {"canonical_gene": "g", "source_gene_id": "id", "source_dataset": "a", "dataset_role": "discovery", "antibiotic_class": "beta_lactam",
         "eligible_for_network": True, "padj": 0.001, "log2FoldChange": 4, "lfc_ci_low": 3, "lfc_ci_high": 5},
        {"canonical_gene": "g", "source_gene_id": "id", "source_dataset": "b", "dataset_role": "cross_reactivity_challenge", "antibiotic_class": "other",
         "eligible_for_network": True, "padj": 0.9, "log2FoldChange": 0, "lfc_ci_low": -2, "lfc_ci_high": 2}])
    pairs = pd.DataFrame([
        {"gene_id": "id", "numerator_dataset": "b", "denominator_dataset": "a", "log2FoldChange": -4,
         "lfc_ci_low": -5, "lfc_ci_high": -3, "eligible_for_network": True, "pvalue": 0.001, "padj": 0.001},
        {"gene_id": "other_gene", "numerator_dataset": "b", "denominator_dataset": "a", "log2FoldChange": 0,
         "lfc_ci_low": -1, "lfc_ci_high": 1, "eligible_for_network": True, "pvalue": 0.9, "padj": 0.9}])
    row = cross_reactivity(observations, pairs, config, 1).iloc[0]
    assert (row.direct_log2_difference, row.direct_ci_low, row.direct_ci_high) == (4, 3, 5)
    assert row.family_padj == pytest.approx(0.002)  # includes the non-candidate gene
    assert row.direct_status == "discovery_response_stronger"
    assert not row.counterpart_interval_within_bound  # non-significance does not establish a bounded effect
    pairs.loc[0, "eligible_for_network"] = False
    assert cross_reactivity(observations, pairs, config, 1).iloc[0].direct_status == "inconclusive"


def test_stability_reports_fixed_model_screens_and_class_specific_rank():
    rows = [{"canonical_gene": g, "source_gene_id": g, "dataset_role": "discovery", "antibiotic_class": "beta_lactam",
             "eligible_for_network": True, "padj": 0.001, "log2FoldChange": fc, "shrunk_log2FoldChange": shrunk}
            for g, fc, shrunk in [("stable", 3, 2.8), ("fragile", 5, 2)]]
    raw = pd.DataFrame([[100] * 3, [10] * 3], index=["stable", "fragile"])
    config = {"thresholds": {"padj": 0.05, "log2_fold_change": 2}, "benchmark": {"min_count": 10, "min_samples": 3}}
    result = candidate_stability(pd.DataFrame(rows), raw, config).set_index("gene")
    assert result.loc["stable", "screen_fraction"] == 1
    assert result.loc["fragile", "screen_fraction"] == pytest.approx(1 / 3)
    assert result.loc["fragile", "raw_rank"] == 1
    assert result.loc["stable", "shrunk_rank"] == 1


def test_tu_context_collapses_overlapping_units_without_altering_enrichment():
    tus = pd.DataFrame({"tuGenes": ["a;b", "b;c", "d", "a;b"]})
    groups = transcription_groups(tus)
    result = operon_support(pd.DataFrame([{"overlap_genes": "a, b, c, d, unannotated", "padj": 0.01}]), groups).iloc[0]
    assert result.n_transcription_groups == 3
    assert result.n_genes_without_tu == 1
    assert result.dominant_group_gene_fraction == pytest.approx(3 / 5)
    assert result.padj == 0.01


@pytest.mark.parametrize("strand", ["forward", "reverse"])
def test_fragments_include_sites_outside_window_and_use_transcription_orientation(strand):
    genome = "ACGTCAGT" * 30
    tss = 100
    a, b = (tss - 3, tss + 2) if strand == "forward" else (tss - 2, tss + 3)
    annotated = oriented_sequence(genome, a, b, strand).lower()
    annotated = annotated[:3] + annotated[3].upper() + annotated[4:]
    promoters = pd.DataFrame([{"id": "p", "posTSS": str(tss), "strand": strand, "sequence": annotated,
                               "boxMinus10pos": "", "boxMinus35pos": ""}])
    review = pd.DataFrame([{"promoter_id": "p", "promoter_name": "test", "candidate_genes": "g", "sigma_factor": "sigma70", "promoter_confidence": "S"}])
    sites = pd.DataFrame([{"promoterID": "p", "tfrsLeft": "75", "tfrsRight": "130", "tfrsID": "site", "regulatorName": "TF"}])
    genes = pd.DataFrame([{"geneName": "nearby", "leftEndPos": "125", "rightEndPos": "150"}])
    result, _ = fragment_review(review, promoters, sites, genes, genome, 10, 5)
    row = result.iloc[0]
    assert (row.fragment_start, row.fragment_end) == (75, 130)
    assert row.fragment_sequence == oriented_sequence(genome, 75, 130, strand)
    assert row.reference_sequence_matches
    assert row.n_sites_outside_initial_window == 1
    assert row.overlapping_genes == "nearby"
    promoters.loc[0, "sequence"] = "nnnAnn"
    assert not fragment_review(review, promoters, sites, genes, genome, 10, 5)[0].iloc[0].reference_sequence_matches
    promoters.loc[0, "posTSS"] = ""
    missing = fragment_review(review, promoters, sites, genes, genome, 10, 5)[0]
    assert len(missing) == 1 and not missing.iloc[0].reference_sequence_matches
    assert missing.iloc[0].construct_status.startswith("missing_tss_or_strand")


def test_panel_balances_classes_without_reusing_promoters_or_tu_groups():
    rows = []
    for cls in ["beta_lactam", "aminoglycoside"]:
        for gene, promoter, group, effect, match in [("a", "p1", "ab", 5, True), ("b", "p2", "ab", 4, True),
                                                    ("c", "p3", "c", 3, True), ("d", "p4", "d", 10, False)]:
            rows.append({"gene": gene, "promoter_id": promoter, "transcription_group": group, "antibiotic_class": cls,
                         "shrunk_effect": effect, "reference_sequence_matches": match, "n_supporting_comparisons": 1,
                         "screen_fraction": 1, "confidence_rank": 2, "regulators": "tf", "direct_stronger_fraction": 0.5,
                         "challenge_max_abs_log2FC": 1})
    result = select_panel(pd.DataFrame(rows), 1)
    assert result.antibiotic_class.tolist() == ["beta_lactam", "aminoglycoside"]
    assert result.gene.tolist() == ["a", "c"]
    assert result.promoter_id.is_unique and result.transcription_group.is_unique
    assert not result.gene.eq("d").any()
