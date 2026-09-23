"""Build the generated SQLite query layer from existing pipeline artifacts."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable

from .database import DEFAULT_DATABASE, close_database, open_database, transaction


TARGETS = {
    "ceftazidime": ("cephalexin", "beta_lactam", "proxy"),
    "imipenem": ("amoxicillin", "beta_lactam", "proxy"),
    "kanamycin": ("gentamicin", "aminoglycoside", "proxy"),
    "ciprofloxacin": ("cross_reactivity", "other", "challenge"),
    "polymyxin_e": ("cross_reactivity", "other", "challenge"),
}


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


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


def _sha256(paths: Iterable[Path]) -> str:
    digest = hashlib.sha256()
    for path in sorted(paths):
        digest.update(str(path).encode())
        digest.update(path.read_bytes())
    return digest.hexdigest()


def _required_files(root: Path) -> list[Path]:
    results = root / "results"
    files = [
        results / "candidate_assessment" / "experimental_panel.csv",
        results / "candidate_assessment" / "fragment_review.csv",
        results / "candidate_assessment" / "operon_support.csv",
        results / "candidate_assessment" / "sample_qc.csv",
        results / "promoter_candidates" / "annotated_candidates.csv",
        results / "promoter_review" / "promoter_review.csv",
        results / "promoter_review" / "candidate_promoter_mapping.csv",
    ]
    files.extend((results / "differential_expression" / "contrasts").glob("*/de_results.csv"))
    missing = [path for path in files if not path.exists()]
    if missing:
        raise FileNotFoundError("Database inputs missing: " + ", ".join(map(str, missing)))
    return files


def _insert_run(connection, run_id: str, input_hash: str) -> None:
    connection.execute(
        "INSERT INTO analysis_runs(run_id, created_at, status, input_sha256, reference_release) "
        "VALUES (?, ?, 'building', ?, ?)",
        (run_id, datetime.now(timezone.utc).isoformat(), input_hash, "RegulonDB 14.5.0"),
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


def _insert_de(connection, root: Path, run_id: str) -> None:
    contrast_root = root / "results/differential_expression/contrasts"
    for path in sorted(contrast_root.glob("*/de_results.csv")):
        observed = path.parent.name
        rows = _read_csv(path)
        if not rows:
            continue
        first = rows[0]
        comparison_id = f"{run_id}:{observed}"
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


def _insert_samples(connection, root: Path, run_id: str) -> None:
    study_id = "GSE220559"
    connection.execute("INSERT OR IGNORE INTO studies(study_id, accession, organism, strain, source_archive) VALUES (?, ?, ?, ?, ?)", (study_id, study_id, "Escherichia coli", "K-12", "GSE220559_RAW.tar"))
    for row in _read_csv(root / "results/candidate_assessment/sample_qc.csv"):
        sample_id = row["sample_id"]
        connection.execute("INSERT OR IGNORE INTO samples(sample_id, study_id, source_file, condition) VALUES (?, ?, ?, ?)", (sample_id, study_id, row.get("source_file"), row.get("condition") or "unknown"))
        connection.execute("INSERT INTO sample_qc(run_id, study_id, sample_id, library_size, genes_detected, pc1, pc2) VALUES (?, ?, ?, ?, ?, ?, ?)", (run_id, study_id, sample_id, _number(row.get("library_counts")), _number(row.get("genes_detected")), _number(row.get("PC1")), _number(row.get("PC2"))))


def _insert_promoters(connection, root: Path) -> None:
    rows = _read_csv(root / "results/promoter_review/promoter_review.csv")
    for row in rows:
        promoter_id = row.get("promoter_id")
        if not promoter_id:
            continue
        tu_id = (row.get("tu_ids") or "").split(";")[0] or None
        if tu_id:
            connection.execute("INSERT OR IGNORE INTO transcription_units(tu_id, name, source_release) VALUES (?, ?, ?)", (tu_id, row.get("operon_name"), row.get("regulondb_release")))
        connection.execute("INSERT OR IGNORE INTO promoters(promoter_id, tu_id, name, start, end, strand, tss, sigma_factor, sequence, annotation_status, source_release) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)", (promoter_id, tu_id, row.get("promoter_name"), None, None, _strand(row.get("strand")), _number(row.get("tss")), row.get("sigma_factor"), row.get("annotated_sequence"), row.get("mapping_status") or "review", row.get("regulondb_release")))
        for gene_id in (row.get("candidate_genes") or row.get("tu_genes") or "").split(";"):
            if gene_id:
                resolved = connection.execute("SELECT gene_id FROM genes WHERE gene_id = ? OR canonical_name = ? LIMIT 1", (gene_id, gene_id)).fetchone()
                resolved_id = resolved[0] if resolved else gene_id
                connection.execute("INSERT OR IGNORE INTO genes(gene_id, canonical_name) VALUES (?, ?)", (resolved_id, gene_id))
                connection.execute("INSERT OR IGNORE INTO promoter_genes(promoter_id, gene_id, relationship) VALUES (?, ?, ?)", (promoter_id, resolved_id, "candidate_or_tu"))
        if tu_id:
            connection.execute("INSERT OR IGNORE INTO tu_promoters(tu_id, promoter_id) VALUES (?, ?)", (tu_id, promoter_id))


def _insert_candidates(connection, root: Path, run_id: str) -> None:
    rows = _read_csv(root / "results/promoter_candidates/annotated_candidates.csv")
    _insert_genes(connection, rows)
    for row in rows:
        candidate_id = row.get("gene") or row.get("gene_id")
        if not candidate_id:
            continue
        connection.execute("INSERT OR IGNORE INTO candidates(candidate_id, run_id, gene_id, promoter_id, antibiotic_class, support_tier, ranking_score, response_summary_json, review_flags_json) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)", (candidate_id, run_id, candidate_id, row.get("promoter_id") or None, row.get("shortlist_class") or row.get("group") or "unknown", row.get("class_support_tier") or row.get("evidence_tier"), _number(row.get("class_max_abs_log2_fold_change") or row.get("max_abs_log2_fold_change")), _json(row.get("response_profile_json")), _json(row.get("evidence_quality_flags"))))
    panel_path = root / "results/candidate_assessment/experimental_panel.csv"
    if panel_path.exists():
        for row in _read_csv(panel_path):
            candidate_id = row.get("gene")
            if not candidate_id:
                continue
            connection.execute("INSERT OR IGNORE INTO candidate_panels(run_id, panel_id, candidate_id, promoter_id, antibiotic_class, selection_order, selection_rationale, transcription_group) VALUES (?, ?, ?, ?, ?, ?, ?, ?)", (run_id, "experimental_panel", candidate_id, row.get("promoter_id") or None, row.get("antibiotic_class") or "unknown", _number(row.get("selection_order")), row.get("selection_basis"), row.get("transcription_group") or None))


def _insert_assessment(connection, root: Path, run_id: str) -> None:
    for row in _read_csv(root / "results/candidate_assessment/fragment_review.csv"):
        connection.execute("INSERT OR IGNORE INTO fragment_reviews(fragment_id, run_id, promoter_id, genome_accession, start, end, strand, sequence, sequence_match, missing_coordinates, regulatory_site_count, overlapping_genes_json, review_status) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)", (row.get("promoter_id"), run_id, row.get("promoter_id"), row.get("genome_accession"), _number(row.get("fragment_start")), _number(row.get("fragment_end")), _strand(row.get("strand")), row.get("fragment_sequence"), int(bool(_boolean(row.get("reference_sequence_matches"), False))), int(_number(row.get("fragment_start")) is None or _number(row.get("fragment_end")) is None), _number(row.get("n_curated_sites")), _json(row.get("overlapping_genes")), row.get("construct_status") or "review"))
    for row in _read_csv(root / "results/candidate_assessment/operon_support.csv"):
        support_id = f"{row.get('regulator')}:{row.get('source_dataset')}:{row.get('direction')}"
        connection.execute("INSERT OR IGNORE INTO operon_support(support_id, run_id, regulated_set, dominant_group_fraction, enrichment, padj, direction, metadata_json) VALUES (?, ?, ?, ?, ?, ?, ?, ?)", (support_id, run_id, row.get("overlap_genes"), _number(row.get("dominant_group_gene_fraction")), _number(row.get("odds_ratio")), _number(row.get("padj")), row.get("direction"), json.dumps(row)))


def build_database(root: str | Path = ".", database_path: str | Path = DEFAULT_DATABASE) -> Path:
    root = Path(root).resolve()
    database_path = Path(database_path)
    if not database_path.is_absolute():
        database_path = root / database_path
    inputs = _required_files(root)
    input_hash = _sha256(inputs)
    run_id = f"run-{input_hash[:16]}"
    database_path.parent.mkdir(parents=True, exist_ok=True)
    temp_handle, temp_name = tempfile.mkstemp(prefix="promoter-discovery-", suffix=".sqlite", dir=database_path.parent)
    os.close(temp_handle)
    temp_path = Path(temp_name)
    connection = None
    try:
        connection = open_database(temp_path)
        with transaction(connection):
            _insert_run(connection, run_id, input_hash)
            _insert_antibiotics(connection)
            _insert_de(connection, root, run_id)
            _insert_samples(connection, root, run_id)
            _insert_promoters(connection, root)
            _insert_candidates(connection, root, run_id)
            _insert_assessment(connection, root, run_id)
            connection.execute("UPDATE analysis_runs SET status='complete', summary_json=? WHERE run_id=?", (json.dumps({"input_files": len(inputs)}), run_id))
        close_database(connection)
        connection = None
        os.replace(temp_path, database_path)
    except Exception:
        if connection is not None:
            connection.close()
        temp_path.unlink(missing_ok=True)
        raise
    return database_path


def main(argv=None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path("."))
    parser.add_argument("--database", type=Path, default=DEFAULT_DATABASE)
    args = parser.parse_args(argv)
    print(build_database(args.root, args.database))


if __name__ == "__main__":
    main()
