"""Shared dataset and evidence-contract definitions.

This module is intentionally dependency-light.  It is imported by the network
builder, iModulon annotator, scorer, and tests so that dataset names, class
keys, caveats, and threshold validation cannot drift between pipeline stages.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any


CLASS_REGISTRY: dict[str, dict[str, Any]] = {
    "beta_lactam": {"label": "Beta-lactam", "datasets": ("amoxicillin", "ceftazidime")},
    "aminoglycoside": {"label": "Aminoglycoside", "datasets": ("gentamicin", "tobramycin")},
}

DATASET_CAVEATS: dict[str, str] = {
    "amoxicillin": "Resistant-mutant/fluoxetine comparison; not a direct wild-type antibiotic exposure.",
    "ceftazidime": "Gene identifiers are supplied as b-number locus tags and require same-release mapping.",
    "gentamicin": "Probe-to-gene annotation contains duplicate symbols; collapse is deterministic and provenance-preserving.",
    "tobramycin": "Independent aminoglycoside dataset; evidence is limited when this is the only qualifying dataset.",
}

REQUIRED_DE_COLUMNS = {"gene", "log2FoldChange", "padj"}
VALID_CLASSES = set(CLASS_REGISTRY) | {"cross", "tf"}


def normalize_version(value: str | None) -> tuple[int, ...] | None:
    """Return a comparable numeric version tuple, accepting ``14.5``."""

    if value is None:
        return None
    text = str(value).strip().lstrip("vV")
    pieces = text.split(".")
    if not pieces or any(not piece.isdigit() for piece in pieces):
        return None
    return tuple(int(piece) for piece in pieces)


def version_at_least(value: str | None, minimum: str) -> bool:
    actual = normalize_version(value)
    wanted = normalize_version(minimum)
    if actual is None or wanted is None:
        return False
    padded_actual = actual + (0,) * max(0, len(wanted) - len(actual))
    padded_wanted = wanted + (0,) * max(0, len(actual) - len(wanted))
    return padded_actual >= padded_wanted


def load_dataset_config(path: str | Path) -> dict[str, Any]:
    """Load and validate the repository dataset configuration."""

    config_path = Path(path)
    with config_path.open(encoding="utf-8") as handle:
        config = json.load(handle)
    thresholds = config.get("thresholds")
    datasets = config.get("datasets")
    if not isinstance(thresholds, dict) or not isinstance(datasets, list) or not datasets:
        raise ValueError(f"{config_path} must contain non-empty thresholds and datasets")
    for key in ("log2_fold_change", "padj"):
        if key not in thresholds:
            raise ValueError(f"{config_path} thresholds missing {key!r}")
    names: set[str] = set()
    for dataset in datasets:
        if not isinstance(dataset, dict):
            raise ValueError("Each configured dataset must be an object")
        for key in ("name", "antibiotic_class"):
            if not dataset.get(key):
                raise ValueError(f"Dataset entry missing {key!r}")
        name = str(dataset["name"])
        if name in names:
            raise ValueError(f"Duplicate dataset name: {name}")
        names.add(name)
        if dataset["antibiotic_class"] not in CLASS_REGISTRY:
            raise ValueError(
                f"Dataset {name} has unsupported antibiotic_class "
                f"{dataset['antibiotic_class']!r}; expected {sorted(CLASS_REGISTRY)}"
            )
    config["_path"] = str(config_path.resolve())
    return config


def configured_datasets_by_class(config: dict[str, Any]) -> dict[str, list[str]]:
    result = {key: [] for key in CLASS_REGISTRY}
    for dataset in config["datasets"]:
        result[dataset["antibiotic_class"]].append(dataset["name"])
    return result


def dataset_caveats(name: str) -> list[str]:
    caveat = DATASET_CAVEATS.get(name)
    return [caveat] if caveat else []
