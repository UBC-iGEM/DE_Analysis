"""Import documented biological context and source identifier associations."""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path
import re


def read_rows(path: Path) -> list[dict]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def insert_identifiers(connection, root: Path, run_id: str, release: str) -> None:
    asset = "results/regulatory_network/source_observations.csv"
    for index, row in enumerate(read_rows(root / asset)):
        gene = row["canonical_gene"]
        connection.execute("INSERT OR IGNORE INTO genes(gene_id,canonical_name) VALUES (?,?)", (gene, gene))
        for namespace, field in (("gene_symbol", "canonical_gene"), ("source_gene_id", "source_gene_id"),
                                 ("locus_tag", "canonical_locus_tag")):
            identifier = row.get(field)
            if identifier:
                connection.execute(
                    "INSERT INTO gene_identifiers(run_id,gene_id,namespace,identifier,source_asset_id,"
                    "source_row,source_release,mapping_basis) VALUES (?,?,?,?,?,?,?,?)",
                    (run_id, gene, namespace, identifier, asset, index, release,
                     row.get("identity_source") or "imported source-observation association"),
                )
    # A symbol may refer to more than one locus. Never select one arbitrarily.
    connection.execute(
        "UPDATE genes SET locus_tag=(SELECT MIN(identifier) FROM gene_identifiers i "
        "WHERE i.gene_id=genes.gene_id AND i.run_id=? AND i.namespace='locus_tag' "
        "HAVING COUNT(DISTINCT identifier)=1) WHERE gene_id IN "
        "(SELECT gene_id FROM gene_identifiers WHERE run_id=?)", (run_id, run_id),
    )


def insert_sample_metadata(connection, root: Path, run_id: str, config: dict) -> None:
    context_asset = "config/study_context.json"
    metadata_asset = "results/differential_expression/model/joint_metadata.csv"
    context = json.loads((root / context_asset).read_text())
    study = config.get("benchmark", {}).get("study", "GSE220559")
    if context["study_id"] != study:
        raise ValueError("Study context does not match the benchmark study")
    for field in ("source_reference", "source_section", "curation_note"):
        if not isinstance(context.get(field), str) or not context[field].strip():
            raise ValueError(f"Study context requires {field}")
    metadata = read_rows(root / metadata_asset)
    # Older pandas exports used an unnamed sample index. Preserve its values.
    for row in metadata:
        if "sample_id" not in row and "" in row:
            row["sample_id"] = row.pop("")
    qc = read_rows(root / "results/candidate_assessment/sample_qc.csv")
    by_sample = {row["sample_id"]: row for row in metadata}
    qc_ids = [row["sample_id"] for row in qc]
    if len(by_sample) != len(metadata) or len(set(qc_ids)) != len(qc_ids) or set(qc_ids) != set(by_sample):
        raise ValueError("Sample metadata and QC must contain the same unique samples")
    datasets = {item["name"]: item for item in config["datasets"]}
    for row in qc:
        sample = row["sample_id"]
        original = by_sample[sample]
        condition = original["condition"]
        if row["condition"] != condition or row.get("source_file") != original.get("source_file"):
            raise ValueError(f"Sample metadata conflicts with QC: {sample}")
        if condition != "water" and condition not in datasets:
            raise ValueError(f"Unknown sample condition: {condition}")
        if condition not in context["conditions"]:
            raise ValueError(f"Study context is missing condition: {condition}")
        exposure = context["conditions"][condition]
        for field in ("concentration", "exposure_minutes"):
            number = exposure.get(field)
            if number is not None and (type(number) not in (int, float) or not math.isfinite(number) or number < 0):
                raise ValueError(f"{condition}: {field} must be a finite nonnegative number or null")
        unit = exposure.get("concentration_unit")
        if exposure.get("concentration") is None:
            if unit is not None:
                raise ValueError(f"{condition}: concentration and unit must both be unknown or both supplied")
        elif not isinstance(unit, str) or not unit.strip():
            raise ValueError(f"{condition}: numeric concentration requires a unit")
        expected_drug = None if condition == "water" else datasets[condition]["treatment"]
        if exposure.get("antibiotic_id") != expected_drug:
            raise ValueError(f"Study context compound conflicts with sample condition: {condition}")
        replicate = original.get("replicate") or None
        replicate_basis = "sample metadata" if replicate is not None else "unavailable"
        if replicate is None:
            choices = list(datasets.values()) if condition == "water" else [datasets[condition]]
            labels = set()
            for dataset in choices:
                pattern = dataset.get("sample_name_regex")
                match = re.search(pattern, original.get("source_file") or sample) if pattern else None
                if match and match.groupdict().get("replicate"):
                    labels.add(match.group("replicate"))
            if len(labels) > 1:
                raise ValueError(f"Conflicting replicate labels for {sample}")
            if labels:
                replicate = labels.pop()
                replicate_basis = "replicate label parsed from source filename with configured sample_name_regex"
        role = "control" if condition == "water" else datasets[condition].get("role")
        connection.execute(
            "UPDATE samples SET replicate=?,antibiotic_role=?,metadata_json=? WHERE study_id=? AND sample_id=?",
            (replicate, role, json.dumps({"original": original, "source_asset": metadata_asset,
                                         "replicate_basis": replicate_basis, "role_source": "config/benchmark.json"}),
             study, sample),
        )
        connection.execute(
            "INSERT INTO sample_exposures(run_id,study_id,sample_id,antibiotic_id,medium,growth_phase,"
            "concentration,concentration_unit,concentration_basis,exposure_minutes,source_asset_id,"
            "source_reference,curation_note) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)",
            (run_id, study, sample, expected_drug, context.get("medium"), context.get("growth_phase"),
             exposure.get("concentration"), exposure.get("concentration_unit"), exposure.get("concentration_basis"),
             exposure.get("exposure_minutes"), context_asset, context["source_reference"],
             context["source_section"] + ": " + context["curation_note"]),
        )
