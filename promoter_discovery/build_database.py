"""Build the generated SQLite query layer from existing pipeline artifacts."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import pickle
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable

from .database import DEFAULT_DATABASE, close_database, open_database, transaction
from .curation import archive_path, preserve_existing_reviews, restore_reviews


TARGETS = {
    "ceftazidime": ("cephalexin", "beta_lactam", "proxy"),
    "imipenem": ("amoxicillin", "beta_lactam", "proxy"),
    "kanamycin": ("gentamicin", "aminoglycoside", "proxy"),
    "ciprofloxacin": ("cross_reactivity", "other", "challenge"),
    "polymyxin_e": ("cross_reactivity", "other", "challenge"),
}
REFERENCE_FILES = (
    "data/references/regulondb/TUSet.tsv",
    "data/references/regulondb/PromoterSet.tsv",
    "data/references/regulondb/RISet.tsv",
)


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _read_reference(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        reader.fieldnames = [name.split(")", 1)[-1].strip() for name in reader.fieldnames or []]
        return [{key: value.strip() for key, value in row.items()} for row in reader]


def _number(value: str | None, default=None):
    if value in (None, "", "NA", "NaN", "nan"):
        return default
    try:
        number = float(value)
        return int(number) if number.is_integer() else number
    except (TypeError, ValueError):
        return default


def _boolean(value: str | None, default=None):
    if value in (None, ""):
        return default
    return str(value).strip().lower() in {"1", "true", "yes"}


def _strand(value: str | None):
    normalized = (value or "").strip().lower()
    return {"forward": "+", "reverse": "-", "+": "+", "-": "-"}.get(normalized, "?") if normalized else None


def _json(value: str | None):
    if value in (None, ""):
        return None
    try:
        json.loads(value)
        return value
    except (TypeError, json.JSONDecodeError):
        return json.dumps(value)


def _file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _required_files(root: Path, config: dict) -> list[Path]:
    results = root / "results"
    files = [
        root / "config/benchmark.json",
        root / "config/regulondb.lock.json",
        results / "differential_expression/model/benchmark_provenance.json",
        results / "differential_expression/model/drug_pairwise_contrasts.csv",
        results / "regulatory_network/regulatory_network.pkl",
        results / "regulatory_network/source_observations.csv",
        results / "regulatory_network/scoring_provenance.json",
        results / "candidate_assessment" / "experimental_panel.csv",
        results / "candidate_assessment" / "fragment_review.csv",
        results / "candidate_assessment" / "operon_support.csv",
        results / "candidate_assessment" / "regulatory_sites.csv",
        results / "candidate_assessment" / "sample_qc.csv",
        results / "promoter_candidates" / "annotated_candidates.csv",
        results / "promoter_review" / "promoter_review.csv",
        results / "promoter_review" / "candidate_promoter_mapping.csv",
    ]
    files.extend(
        results / "differential_expression/contrasts" / dataset["output_directory"] / "de_results.csv"
        for dataset in config["datasets"]
    )
    files.extend(root / relative for relative in REFERENCE_FILES)
    missing = [path for path in files if not path.exists()]
    if missing:
        raise FileNotFoundError("Database inputs missing: " + ", ".join(map(str, missing)))
    return files


def _verify_reference_files(root: Path, release: str) -> None:
    lock = json.loads((root / "config/regulondb.lock.json").read_text())
    records = {item.get("path"): item for item in lock.get("assets", {}).values()}
    for relative in REFERENCE_FILES:
        record = records.get(relative)
        if (not record or record.get("regulondb_release") != release
                or record.get("sha256") != _file_sha256(root / relative)):
            raise ValueError(f"Reference is not verified by the {release} lock: {relative}")


def _input_hashes(root: Path, files: Iterable[Path]) -> tuple[str, dict[str, str]]:
    hashes = {path.relative_to(root).as_posix(): _file_sha256(path) for path in files}
    combined = hashlib.sha256()
    for path, digest in sorted(hashes.items()):
        combined.update(f"{path}:{digest}\n".encode())
    return combined.hexdigest(), hashes


def _insert_run(connection, run_id: str, input_hash: str, config_hash: str, release: str) -> None:
    connection.execute(
        "INSERT INTO analysis_runs(run_id, created_at, status, config_sha256, input_sha256, reference_release) "
        "VALUES (?, ?, 'building', ?, ?, ?)",
        (run_id, datetime.now(timezone.utc).isoformat(), config_hash, input_hash, release),
    )


def _insert_assets(connection, root: Path, run_id: str, files: list[Path], hashes: dict[str, str], release: str) -> None:
    for path in files:
        relative = path.relative_to(root).as_posix()
        connection.execute(
            "INSERT INTO assets(asset_id, run_id, path, kind, source_release, sha256, size_bytes) "
            "VALUES (?, ?, ?, ?, ?, ?, ?)",
            (relative, run_id, relative, path.suffix.lstrip(".") or "file", release,
             hashes[relative], path.stat().st_size),
        )


def _insert_antibiotics(connection) -> None:
    for observed, (target, antibiotic_class, proxy_role) in TARGETS.items():
        connection.execute(
            "INSERT OR IGNORE INTO antibiotics(antibiotic_id, target_name, observed_compound, antibiotic_class, proxy_role) "
            "VALUES (?, ?, ?, ?, ?)",
            (observed, target, observed, antibiotic_class, proxy_role),
        )


def _insert_genes(connection, rows: Iterable[dict[str, str]]) -> None:
    for row in rows:
        gene_id = row.get("gene_id") or row.get("gene")
        if not gene_id:
            continue
        connection.execute(
            "INSERT INTO genes(gene_id, canonical_name, locus_tag) VALUES (?, ?, ?) "
            "ON CONFLICT(gene_id) DO UPDATE SET canonical_name=COALESCE(excluded.canonical_name, genes.canonical_name), locus_tag=COALESCE(excluded.locus_tag, genes.locus_tag)",
            (gene_id, row.get("gene") or gene_id, row.get("canonical_locus_tag") or row.get("source_locus_tag")),
        )


def _insert_de(connection, root: Path, run_id: str, config: dict) -> dict[str, str]:
    contrast_root = root / "results/differential_expression/contrasts"
    contrast_ids = {}
    for dataset in config["datasets"]:
        observed = dataset["treatment"]
        path = contrast_root / dataset["output_directory"] / "de_results.csv"
        rows = _read_csv(path)
        if not rows:
            continue
        first = rows[0]
        comparison_id = f"{run_id}:{observed}"
        contrast_ids[dataset["name"]] = comparison_id
        target, antibiotic_class, proxy_role = TARGETS.get(observed, (observed, "other", "control"))
        connection.execute(
            "INSERT INTO contrasts(contrast_id, run_id, antibiotic_id, name, numerator, denominator, comparison_type, antibiotic_class, proxy_role) "
            "VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)",
            (comparison_id, run_id, observed, first.get("comparison") or f"{observed}_vs_water", observed, "water", "drug_control", antibiotic_class, proxy_role),
        )
        _insert_genes(connection, rows)
        for index, row in enumerate(rows):
            connection.execute(
                "INSERT INTO de_results(run_id, contrast_id, gene_id, base_mean, log2_fold_change, log2_fold_change_se, "
                "shrunk_log2_fold_change, shrunk_lfc_se, p_value, padj, ci_low, ci_high, regulation, eligible, source_row_id) "
                "VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
                (run_id, comparison_id, row.get("gene_id") or row.get("gene"), _number(row.get("baseMean")),
                 _number(row.get("log2FoldChange")), _number(row.get("lfcSE")), _number(row.get("shrunk_log2FoldChange")),
                 _number(row.get("shrunk_lfcSE")), _number(row.get("pvalue")), _number(row.get("padj")),
                 _number(row.get("lfc_ci_low")), _number(row.get("lfc_ci_high")), row.get("regulation"),
                 int(bool(_boolean(row.get("eligible_for_network"), True))), f"{observed}:{index}"),
            )
    return contrast_ids


def _insert_pairwise(connection, root: Path, run_id: str, config: dict) -> None:
    datasets = {dataset["name"]: dataset for dataset in config["datasets"]}
    rows = _read_csv(root / "results/differential_expression/model/drug_pairwise_contrasts.csv")
    for index, row in enumerate(rows):
        numerator = row["numerator_dataset"]
        denominator = row["denominator_dataset"]
        if numerator not in datasets or denominator not in datasets:
            raise ValueError(f"Unknown direct-comparison dataset: {numerator}, {denominator}")
        contrast_id = f"{run_id}:pair:{numerator}:{denominator}"
        connection.execute(
            "INSERT OR IGNORE INTO contrasts(contrast_id, run_id, name, numerator, denominator, comparison_type, proxy_role) "
            "VALUES (?, ?, ?, ?, ?, 'drug_drug', 'direct_comparison')",
            (contrast_id, run_id, f"{numerator}_vs_{denominator}", numerator, denominator),
        )
        gene_id = row["gene_id"]
        connection.execute("INSERT OR IGNORE INTO genes(gene_id) VALUES (?)", (gene_id,))
        connection.execute(
            "INSERT INTO de_results(run_id, contrast_id, gene_id, base_mean, log2_fold_change, "
            "log2_fold_change_se, p_value, padj, ci_low, ci_high, regulation, eligible, source_row_id) "
            "VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
            (run_id, contrast_id, gene_id, _number(row.get("baseMean")),
             _number(row.get("log2FoldChange")), _number(row.get("lfcSE")),
             _number(row.get("pvalue")), _number(row.get("padj")),
             _number(row.get("lfc_ci_low")), _number(row.get("lfc_ci_high")),
             row.get("regulation"), int(bool(_boolean(row.get("eligible_for_network"), True))),
             f"pairwise:{index}"),
        )


def _insert_candidate_evidence(connection, root: Path, run_id: str, contrast_ids: dict[str, str]) -> None:
    candidate_ids = {
        row[0] for row in connection.execute(
            "SELECT candidate_id FROM candidates WHERE run_id = ?", (run_id,)
        )
    }
    for row in _read_csv(root / "results/regulatory_network/source_observations.csv"):
        candidate_id = row["canonical_gene"]
        if candidate_id not in candidate_ids:
            continue
        source_dataset = row["source_dataset"]
        if source_dataset not in contrast_ids:
            raise ValueError(f"Unknown candidate evidence dataset: {source_dataset}")
        connection.execute(
            "INSERT INTO candidate_evidence(run_id, candidate_id, contrast_id, source_gene_id, "
            "effect, padj, ci_low, ci_high, direction, evidence_role) "
            "VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
            (run_id, candidate_id, contrast_ids[source_dataset], row["source_gene_id"],
             _number(row.get("log2FoldChange")), _number(row.get("padj")),
             _number(row.get("lfc_ci_low")), _number(row.get("lfc_ci_high")),
             row.get("regulation"), row.get("dataset_role")),
        )


def _insert_samples(connection, root: Path, run_id: str) -> None:
    study_id = "GSE220559"
    connection.execute("INSERT OR IGNORE INTO studies(study_id, accession, organism, strain, source_archive) VALUES (?, ?, ?, ?, ?)", (study_id, study_id, "Escherichia coli", "K-12", "GSE220559_RAW.tar"))
    for row in _read_csv(root / "results/candidate_assessment/sample_qc.csv"):
        sample_id = row["sample_id"]
        connection.execute("INSERT OR IGNORE INTO samples(sample_id, study_id, source_file, condition) VALUES (?, ?, ?, ?)", (sample_id, study_id, row.get("source_file"), row.get("condition") or "unknown"))
        connection.execute("INSERT INTO sample_qc(run_id, study_id, sample_id, library_size, genes_detected, pc1, pc2) VALUES (?, ?, ?, ?, ?, ?, ?)", (run_id, study_id, sample_id, _number(row.get("library_counts")), _number(row.get("genes_detected")), _number(row.get("PC1")), _number(row.get("PC2"))))


def _insert_reference_network(connection, root: Path, run_id: str, release: str) -> None:
    """Keep every same-release interaction and its promoter/TU context."""
    reference = root / "data/references/regulondb"
    for row in _read_reference(reference / "PromoterSet.tsv"):
        connection.execute(
            "INSERT INTO promoters(promoter_id, name, genome_accession, strand, tss, sigma_factor, "
            "sequence, annotation_status, source_release, annotation_json) "
            "VALUES (?, ?, 'U00096.3', ?, ?, ?, ?, 'reference', ?, ?)",
            (row["id"], row.get("name"), _strand(row.get("strand")), _number(row.get("posTSS")),
             row.get("sigmaFactor"), row.get("sequence"), release, json.dumps(row)),
        )
    for row in _read_reference(reference / "TUSet.tsv"):
        tu_id = row["id"]
        connection.execute(
            "INSERT INTO transcription_units(tu_id, name, source_release, annotation_json) VALUES (?, ?, ?, ?)",
            (tu_id, row.get("name"), release, json.dumps(row)),
        )
        for order, name in enumerate((row.get("tuGenes") or "").split(";")):
            gene_id = name.strip().lower()
            if gene_id:
                connection.execute("INSERT OR IGNORE INTO genes(gene_id, canonical_name) VALUES (?, ?)", (gene_id, gene_id))
                connection.execute("INSERT OR IGNORE INTO tu_genes(tu_id, gene_id, gene_order) VALUES (?, ?, ?)", (tu_id, gene_id, order))
        promoter_id = row.get("promoterId")
        if promoter_id:
            connection.execute("INSERT INTO tu_promoters(tu_id, promoter_id) VALUES (?, ?)", (tu_id, promoter_id))
            connection.execute("UPDATE promoters SET tu_id = COALESCE(tu_id, ?) WHERE promoter_id = ?", (tu_id, promoter_id))

    for row in _read_reference(reference / "RISet.tsv"):
        interaction_type = row["type"]
        actor_type, target_kind = interaction_type.split("-", 1)
        target_kind = target_kind.lower()
        actor_id = row["regulatorId"]
        actor_name = row["regulatorName"]
        connection.execute(
            "INSERT OR IGNORE INTO regulatory_actors(actor_id, name, actor_type, source_release) VALUES (?, ?, ?, ?)",
            (actor_id, actor_name, actor_type, release),
        )
        target = row.get("targetTuOrGene", "")
        source_target_id, _, target_name = target.partition(":")
        promoter_id = row.get("promoterID") if target_kind == "promoter" else None
        if promoter_id:
            source_target_id = promoter_id
        tu_id = source_target_id if target_kind == "tu" else None
        gene_id = (target_name or source_target_id).lower() if target_kind == "gene" else None
        if gene_id:
            connection.execute("INSERT OR IGNORE INTO genes(gene_id, canonical_name) VALUES (?, ?)", (gene_id, gene_id))
        site_id = row.get("tfrsID") or None
        if site_id:
            regulator_id = actor_name.lower() if actor_type in {"TF", "regulator"} else None
            if regulator_id:
                connection.execute(
                    "INSERT OR IGNORE INTO regulators(regulator_id, name, regulator_type, source_release) VALUES (?, ?, ?, ?)",
                    (regulator_id, actor_name, actor_type, release),
                )
            connection.execute(
                "INSERT OR IGNORE INTO regulatory_sites(site_id, promoter_id, regulator_id, start, end, strand, "
                "sequence, function, evidence, confidence, source_release, annotation_json) "
                "VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
                (site_id, row.get("promoterID") or None, regulator_id, _number(row.get("tfrsLeft")),
                 _number(row.get("tfrsRight")), _strand(row.get("strand")), row.get("tfrsSeq"),
                 row.get("riFunction"), row.get("tfrsEvidence"), row.get("confidenceLevel"), release, json.dumps(row)),
            )
        connection.execute(
            "INSERT INTO regulatory_interactions(interaction_id, run_id, actor_id, interaction_type, target_kind, "
            "source_target_id, promoter_id, tu_id, gene_id, site_id, effect, conformation, confidence, "
            "site_evidence, interaction_evidence, evidence_category, pmids, source_release, source_json) "
            "VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
            (row["id"], run_id, actor_id, interaction_type, target_kind, source_target_id or None,
             promoter_id, tu_id, gene_id, site_id, row.get("riFunction"), row.get("cnfName"),
             row.get("confidenceLevel"), row.get("tfrsEvidence"), row.get("riEvidence"),
             row.get("riEvCategory"), row.get("riPMIDS"), release, json.dumps(row)),
        )


def _insert_promoters(connection, root: Path, run_id: str) -> None:
    rows = _read_csv(root / "results/promoter_review/promoter_review.csv")
    for row in rows:
        promoter_id = row.get("promoter_id")
        if not promoter_id:
            continue
        tu_ids = [item for item in (row.get("tu_ids") or "").split(";") if item]
        for tu_id in tu_ids:
            connection.execute("INSERT OR IGNORE INTO transcription_units(tu_id, name, source_release) VALUES (?, ?, ?)", (tu_id, row.get("operon_name"), row.get("regulondb_release")))
        connection.execute(
            "INSERT INTO promoters(promoter_id, tu_id, name, start, end, strand, tss, sigma_factor, sequence, annotation_status, source_release) "
            "VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?) "
            "ON CONFLICT(promoter_id) DO UPDATE SET "
            "tu_id=COALESCE(excluded.tu_id, promoters.tu_id), "
            "name=COALESCE(excluded.name, promoters.name), "
            "strand=COALESCE(excluded.strand, promoters.strand), "
            "tss=COALESCE(excluded.tss, promoters.tss), "
            "sigma_factor=COALESCE(excluded.sigma_factor, promoters.sigma_factor), "
            "sequence=COALESCE(excluded.sequence, promoters.sequence), "
            "annotation_status=excluded.annotation_status",
            (promoter_id, tu_ids[0] if tu_ids else None, row.get("promoter_name"), None, None,
             _strand(row.get("strand")), _number(row.get("tss")), row.get("sigma_factor"),
             row.get("annotated_sequence"), row.get("mapping_status") or "review", row.get("regulondb_release")),
        )
        for tu_id in tu_ids:
            connection.execute("INSERT OR IGNORE INTO tu_promoters(tu_id, promoter_id) VALUES (?, ?)", (tu_id, promoter_id))
        for name in (row.get("candidate_genes") or "").split(";"):
            gene_id = name.strip().lower()
            if gene_id:
                connection.execute("INSERT OR IGNORE INTO genes(gene_id, canonical_name) VALUES (?, ?)", (gene_id, gene_id))
                connection.execute("INSERT OR IGNORE INTO promoter_genes(promoter_id, gene_id, relationship) VALUES (?, ?, 'candidate')", (promoter_id, gene_id))

    for row in _read_csv(root / "results/promoter_review/candidate_promoter_mapping.csv"):
        candidate_id = row["gene"].strip().lower()
        tu_id = row.get("tu_id") or None
        promoter_id = row.get("promoter_id") or None
        if tu_id:
            connection.execute("INSERT OR IGNORE INTO transcription_units(tu_id, name, source_release) VALUES (?, ?, ?)", (tu_id, row.get("operon_name"), row.get("regulondb_release")))
            for order, name in enumerate((row.get("tu_genes") or "").split(";")):
                gene_id = name.strip().lower()
                if gene_id:
                    connection.execute("INSERT OR IGNORE INTO genes(gene_id, canonical_name) VALUES (?, ?)", (gene_id, gene_id))
                    connection.execute("INSERT OR IGNORE INTO tu_genes(tu_id, gene_id, gene_order) VALUES (?, ?, ?)", (tu_id, gene_id, order))
            connection.execute("INSERT OR IGNORE INTO candidate_tus(run_id, candidate_id, tu_id) VALUES (?, ?, ?)", (run_id, candidate_id, tu_id))
        if promoter_id:
            if tu_id:
                connection.execute("INSERT OR IGNORE INTO tu_promoters(tu_id, promoter_id) VALUES (?, ?)", (tu_id, promoter_id))
            connection.execute("INSERT OR IGNORE INTO candidate_promoters(run_id, candidate_id, promoter_id) VALUES (?, ?, ?)", (run_id, candidate_id, promoter_id))
            connection.execute("INSERT OR IGNORE INTO promoter_genes(promoter_id, gene_id, relationship) VALUES (?, ?, 'candidate')", (promoter_id, candidate_id))


def _insert_candidates(connection, root: Path, run_id: str) -> None:
    rows = _read_csv(root / "results/promoter_candidates/annotated_candidates.csv")
    _insert_genes(connection, rows)
    for row in rows:
        candidate_id = row.get("gene") or row.get("gene_id")
        if not candidate_id:
            continue
        connection.execute("INSERT OR IGNORE INTO candidates(candidate_id, run_id, gene_id, promoter_id, antibiotic_class, support_tier, ranking_score, response_summary_json, review_flags_json) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)", (candidate_id, run_id, candidate_id, row.get("promoter_id") or None, row.get("shortlist_class") or row.get("group") or "unknown", row.get("class_support_tier") or row.get("evidence_tier"), _number(row.get("class_max_abs_log2_fold_change") or row.get("max_abs_log2_fold_change")), _json(row.get("response_profile_json")), _json(row.get("evidence_quality_flags"))))


def _insert_panel(connection, root: Path, run_id: str) -> None:
    panel_path = root / "results/candidate_assessment/experimental_panel.csv"
    for row in _read_csv(panel_path):
        candidate_id = row.get("gene")
        if not candidate_id:
            continue
        connection.execute("INSERT INTO candidate_panels(run_id, panel_id, candidate_id, promoter_id, antibiotic_class, selection_order, selection_rationale, transcription_group) VALUES (?, ?, ?, ?, ?, ?, ?, ?)", (run_id, "experimental_panel", candidate_id, row.get("promoter_id") or None, row.get("antibiotic_class") or "unknown", _number(row.get("selection_order")), row.get("selection_basis"), row.get("transcription_group") or None))


def _insert_assessment(connection, root: Path, run_id: str) -> None:
    for row in _read_csv(root / "results/candidate_assessment/fragment_review.csv"):
        connection.execute("INSERT OR IGNORE INTO fragment_reviews(fragment_id, run_id, promoter_id, genome_accession, start, end, strand, sequence, sequence_match, missing_coordinates, regulatory_site_count, overlapping_genes_json, review_status) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)", (row.get("promoter_id"), run_id, row.get("promoter_id"), row.get("genome_accession"), _number(row.get("fragment_start")), _number(row.get("fragment_end")), _strand(row.get("strand")), row.get("fragment_sequence"), int(bool(_boolean(row.get("reference_sequence_matches"), False))), int(_number(row.get("fragment_start")) is None or _number(row.get("fragment_end")) is None), _number(row.get("n_curated_sites")), _json(row.get("overlapping_genes")), row.get("construct_status") or "review"))
    for row in _read_csv(root / "results/candidate_assessment/operon_support.csv"):
        support_id = f"{row.get('regulator')}:{row.get('source_dataset')}:{row.get('direction')}"
        connection.execute("INSERT OR IGNORE INTO operon_support(support_id, run_id, regulated_set, dominant_group_fraction, enrichment, padj, direction, metadata_json) VALUES (?, ?, ?, ?, ?, ?, ?, ?)", (support_id, run_id, row.get("overlap_genes"), _number(row.get("dominant_group_gene_fraction")), _number(row.get("odds_ratio")), _number(row.get("padj")), row.get("direction"), json.dumps(row)))


def _insert_regulatory_sites(connection, root: Path) -> None:
    path = root / "results/candidate_assessment/regulatory_sites.csv"
    if not path.exists():
        return
    for row in _read_csv(path):
        promoter_id = row.get("promoterID") or None
        regulator_id = row.get("regulatorId") or row.get("regulatorName") or None
        if promoter_id:
            connection.execute("INSERT OR IGNORE INTO promoters(promoter_id, name, tss, sigma_factor, annotation_status) VALUES (?, ?, ?, ?, ?)", (promoter_id, row.get("promoterName"), _number(row.get("tss")), row.get("sigmaF"), "site_annotation"))
        if regulator_id:
            connection.execute("INSERT OR IGNORE INTO regulators(regulator_id, name) VALUES (?, ?)", (regulator_id, row.get("regulatorName") or regulator_id))
        connection.execute("INSERT OR IGNORE INTO regulatory_sites(site_id, promoter_id, regulator_id, start, end, strand, sequence, function, evidence, confidence) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)", (row.get("id"), promoter_id, regulator_id, _number(row.get("tfrsLeft")), _number(row.get("tfrsRight")), _strand(row.get("strand")), row.get("tfrsSeq"), row.get("riFunction"), row.get("tfrsEvidence"), row.get("confidenceLevel")))


def _insert_network(connection, root: Path, run_id: str, release: str) -> None:
    """Import every regulatory edge from the saved full-reference graph."""
    path = root / "results/regulatory_network/regulatory_network.pkl"
    try:
        with path.open("rb") as handle:
            graph = pickle.load(handle)
    except (ImportError, ModuleNotFoundError) as exc:
        raise RuntimeError("Network import requires networkx; install requirements.txt") from exc
    if graph.graph.get("network_scope") != "full_reference":
        raise ValueError("Expected the full-reference regulatory graph")
    regulator_nodes = set()
    for node, attrs in graph.nodes(data=True):
        node = str(node)
        node_type = str(attrs.get("node_type", "")).lower()
        connection.execute("INSERT OR IGNORE INTO genes(gene_id, canonical_name, locus_tag) VALUES (?, ?, ?)", (node, attrs.get("canonical_gene") or node, attrs.get("canonical_locus_tag") or None))
        if attrs.get("is_regulator") or node_type in {"regulator", "sigma"}:
            regulator_nodes.add(node)
            connection.execute("INSERT OR IGNORE INTO regulators(regulator_id, name, regulator_type, source_release) VALUES (?, ?, ?, ?)", (node, attrs.get("label") or attrs.get("canonical_gene") or node, attrs.get("regulator_type") or node_type or "regulator", release))
    for source, target, attrs in graph.edges(data=True):
        source, target = str(source), str(target)
        edge_type = attrs.get("edge_type")
        if edge_type not in {"activates", "represses", "dual"}:
            continue
        if source not in regulator_nodes:
            raise ValueError(f"Regulatory edge source {source!r} is not a regulator")
        edge_id = hashlib.sha256(f"{run_id}:{source}:{target}:{edge_type}".encode()).hexdigest()[:24]
        connection.execute("INSERT INTO regulatory_edges(edge_id, run_id, regulator_id, target_gene_id, edge_type, effect, evidence, confidence, source_release, metadata_json) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)", (edge_id, run_id, source, target, edge_type, edge_type, json.dumps(attrs.get("interaction_evidence", []), default=str), attrs.get("confidence"), release, json.dumps(attrs, default=str)))
    count = connection.execute("SELECT COUNT(*) FROM regulatory_edges WHERE run_id = ?", (run_id,)).fetchone()[0]
    if not count:
        raise ValueError("The full-reference graph contains no regulatory edges")


def build_database(root: str | Path = ".", database_path: str | Path = DEFAULT_DATABASE) -> Path:
    root = Path(root).resolve()
    database_path = Path(database_path)
    if not database_path.is_absolute():
        database_path = root / database_path
    config = json.loads((root / "config/benchmark.json").read_text())
    inputs = _required_files(root, config)
    release = json.loads((root / "results/regulatory_network/scoring_provenance.json").read_text())["regulondb_release"]
    _verify_reference_files(root, release)
    input_hash, hashes = _input_hashes(root, inputs)
    run_id = f"run-{input_hash[:16]}"
    database_path.parent.mkdir(parents=True, exist_ok=True)
    temp_handle, temp_name = tempfile.mkstemp(prefix="promoter-discovery-", suffix=".sqlite", dir=database_path.parent)
    os.close(temp_handle)
    temp_path = Path(temp_name)
    connection = None
    try:
        connection = open_database(temp_path)
        with transaction(connection):
            _insert_run(connection, run_id, input_hash, hashes["config/benchmark.json"], release)
            _insert_assets(connection, root, run_id, inputs, hashes, release)
            _insert_antibiotics(connection)
            contrast_ids = _insert_de(connection, root, run_id, config)
            _insert_pairwise(connection, root, run_id, config)
            _insert_samples(connection, root, run_id)
            _insert_reference_network(connection, root, run_id, release)
            _insert_candidates(connection, root, run_id)
            _insert_promoters(connection, root, run_id)
            _insert_panel(connection, root, run_id)
            _insert_regulatory_sites(connection, root)
            _insert_candidate_evidence(connection, root, run_id, contrast_ids)
            _insert_network(connection, root, run_id, release)
            _insert_assessment(connection, root, run_id)
            preserve_existing_reviews(database_path)
            restore_reviews(connection, archive_path(database_path), run_id)
            if connection.execute("PRAGMA foreign_key_check").fetchone():
                raise ValueError("Database contains broken foreign-key relationships")
            coverage = _reference_coverage(connection, root)
            connection.execute("UPDATE analysis_runs SET status='complete', summary_json=? WHERE run_id=?", (json.dumps({"input_files": len(inputs), "reference_coverage": coverage}), run_id))
        connection.execute("PRAGMA journal_mode = DELETE")
        close_database(connection)
        connection = None
        os.replace(temp_path, database_path)
    except Exception:
        if connection is not None:
            connection.close()
        temp_path.unlink(missing_ok=True)
        Path(f"{temp_path}-wal").unlink(missing_ok=True)
        Path(f"{temp_path}-shm").unlink(missing_ok=True)
        raise
    return database_path


def _reference_coverage(connection, root: Path) -> dict:
    """Check complete source import and report unresolved contextual links."""
    counts = {}
    for filename, table in (("RISet.tsv", "regulatory_interactions"),
                            ("PromoterSet.tsv", "promoters"), ("TUSet.tsv", "transcription_units")):
        source_ids = {row["id"] for row in _read_reference(root / "data/references/regulondb" / filename)}
        key = {"regulatory_interactions": "interaction_id", "promoters": "promoter_id",
               "transcription_units": "tu_id"}[table]
        imported_ids = {row[0] for row in connection.execute(f"SELECT {key} FROM {table}")}
        if not source_ids <= imported_ids:
            raise ValueError(f"Incomplete reference import: {filename}")
        counts[table] = len(source_ids)
    paths = connection.execute("SELECT COUNT(DISTINCT interaction_id) FROM reference_regulatory_paths").fetchone()[0]
    if paths != counts["regulatory_interactions"]:
        raise ValueError("Regulatory paths do not cover every reference interaction")
    counts["interactions_by_target"] = dict(connection.execute(
        "SELECT target_kind, COUNT(*) FROM regulatory_interactions GROUP BY target_kind"))
    counts["interactions_without_promoter"] = connection.execute(
        "SELECT COUNT(*) FROM regulatory_interactions ri WHERE NOT EXISTS "
        "(SELECT 1 FROM reference_regulatory_paths p WHERE p.interaction_id=ri.interaction_id AND p.promoter_id IS NOT NULL)"
    ).fetchone()[0]
    return counts


def main(argv=None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path("."))
    parser.add_argument("--database", type=Path, default=DEFAULT_DATABASE)
    args = parser.parse_args(argv)
    print(build_database(args.root, args.database))


if __name__ == "__main__":
    main()
