import numpy as np
import pandas as pd
import pytest
from promoter_discovery import input_data as run_analysis

from promoter_discovery.input_data import (
    benjamini_hochberg,
    classify,
    finalize_results,
    median_ratio_size_factors,
    load_excel_de_results,
)


def test_benjamini_hochberg_preserves_index_and_monotonic_adjustment():
    pvalues = pd.Series([0.04, 0.01, 0.03], index=["gene_c", "gene_a", "gene_b"])
    adjusted = benjamini_hochberg(pvalues)

    assert adjusted.index.tolist() == pvalues.index.tolist()
    assert np.allclose(adjusted.to_numpy(), [0.04, 0.03, 0.04])


def test_classify_uses_strict_configured_thresholds():
    assert classify(2.01, 0.049, 2.0, 0.05) == "upregulated"
    assert classify(-2.01, 0.049, 2.0, 0.05) == "downregulated"
    assert classify(2.0, 0.01, 2.0, 0.05) == "not_regulated"
    assert classify(3.0, 0.05, 2.0, 0.05) == "not_regulated"


def test_finalize_results_computes_adjustment_signal_and_order():
    raw = pd.DataFrame(
        {
            "gene": ["geneA", "geneB"],
            "gene_id": ["b0001", "b0002"],
            "log2FoldChange": [3.0, -1.0],
            "pvalue": [0.01, 0.5],
        }
    )
    dataset = {
        "name": "example",
        "antibiotic_class": "beta_lactam",
        "treatment": "example-drug",
        "comparison": "treated_vs_control",
    }

    result = finalize_results(
        raw,
        dataset,
        {"log2_fold_change": 2.0, "padj": 0.05},
    )

    assert result["gene"].tolist() == ["geneA", "geneB"]
    assert result["regulation"].tolist() == ["upregulated", "not_regulated"]
    assert result["padj_source"].tolist() == ["BH_from_pvalue", "BH_from_pvalue"]
    assert result["comparison"].tolist() == ["treated_vs_control"] * 2
    assert result.iloc[0]["signal_strength"] > result.iloc[1]["signal_strength"]


def test_median_ratio_size_factors_are_scale_invariant():
    counts = pd.DataFrame(
        {
            "sample_a": [10, 20, 30],
            "sample_b": [20, 40, 60],
        }
    )
    factors = median_ratio_size_factors(counts)
    normalized = counts.div(factors, axis=1)

    assert np.allclose(normalized["sample_a"], normalized["sample_b"])


@pytest.fixture
def excel_dataset(tmp_path):
    workbook = tmp_path / "analysis.xlsx"
    with pd.ExcelWriter(workbook, engine="openpyxl") as writer:
        pd.DataFrame({
            "gene_id": ["b0001"], "gene": ["geneA"],
            "log2FoldChange": [3.0], "pvalue": [0.001], "padj": [0.01],
        }).to_excel(writer, sheet_name="MH vs tob", index=False)
        pd.DataFrame({
            "locus": ["b0001"], "EcMH1_S1": [10], "Ectob1_S7": [80],
            "NCMMH1_S73": [20], "NCMtob1_S76": [30],
        }).to_excel(writer, sheet_name="counts", index=False)
    return {
        "name": "tobramycin", "antibiotic_class": "aminoglycoside",
        "input_type": "excel_de_results",
        "treatment": "tobramycin", "workbook": str(workbook),
        "de_sheet": "MH vs tob", "counts_sheet": "counts",
        "control_samples": ["EcMH1_S1"], "treated_samples": ["Ectob1_S7"],
    }


def test_excel_counts_match_explicit_de_comparison(excel_dataset, tmp_path, monkeypatch):
    results, counts, metadata = load_excel_de_results(excel_dataset)
    assert counts.columns.tolist() == ["gene_id", "EcMH1_S1", "Ectob1_S7"]
    assert metadata["sample_id"].tolist() == ["EcMH1_S1", "Ectob1_S7"]
    assert metadata["condition"].tolist() == ["control", "treated"]
    assert metadata["sample_group"].tolist() == ["control", "tobramycin"]
    assert results.loc[0, "log2FoldChange"] == 3.0
    assert results.loc[0, "padj"] == 0.01
    monkeypatch.setattr(run_analysis, "REPO_ROOT", tmp_path)
    run_analysis.write_standardized(excel_dataset, expression=counts, metadata=metadata, results=results)
    output = tmp_path / "results/differential_expression/contrasts/tobramycin"
    assert pd.read_csv(output / "counts.csv").equals(counts)
    assert pd.read_csv(output / "expression.csv").equals(counts)
    assert pd.read_excel(output / "standardized_inputs.xlsx", sheet_name="metadata")["condition"].tolist() == ["control", "treated"]


@pytest.mark.parametrize("samples, message", [
    ([], "requires explicit"),
    (["EcMH1_S1"], "duplicate or overlapping"),
    (["missing_sample"], "Missing sample columns"),
])
def test_excel_counts_reject_invalid_sample_assignments(excel_dataset, samples, message):
    excel_dataset["treated_samples"] = samples
    with pytest.raises(ValueError, match=message):
        load_excel_de_results(excel_dataset)
