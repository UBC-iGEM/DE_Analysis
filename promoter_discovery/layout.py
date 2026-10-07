"""Shared project layout for source data and generated analysis stages."""
from pathlib import Path

DE_MODEL = Path("results/differential_expression/model")
DE_CONTRASTS = Path("results/differential_expression/contrasts")
REGULATORY_NETWORK = Path("results/regulatory_network")
PROMOTER_CANDIDATES = Path("results/promoter_candidates")
PROMOTER_REVIEW = Path("results/promoter_review")
REGULONDB = Path("data/references/regulondb")


def contrast_directory(root, dataset):
    """Keep dataset evidence IDs stable while naming directories by comparison."""
    return Path(root) / DE_CONTRASTS / dataset.get("output_directory", dataset["name"])
