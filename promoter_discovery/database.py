"""SQLite connection and schema helpers for integrated project results.

The database is a generated query layer over the pipeline's reproducible files.
Scientific result tables are intentionally added in later schema migrations;
this module establishes the shared metadata and reference-data foundation.
"""

from __future__ import annotations

import sqlite3
from contextlib import contextmanager
from pathlib import Path
from typing import Iterator


SCHEMA_VERSION = 5
DEFAULT_DATABASE = Path("results/promoter_discovery.sqlite")


SCHEMA_SQL = """
CREATE TABLE IF NOT EXISTS database_metadata (
    key TEXT PRIMARY KEY,
    value TEXT NOT NULL
);

CREATE TABLE IF NOT EXISTS analysis_runs (
    run_id TEXT PRIMARY KEY,
    created_at TEXT NOT NULL,
    status TEXT NOT NULL CHECK (status IN ('building', 'complete', 'failed')),
    pipeline_version TEXT,
    git_revision TEXT,
    config_sha256 TEXT,
    input_sha256 TEXT,
    reference_release TEXT,
    summary_json TEXT
);

CREATE TABLE IF NOT EXISTS assets (
    asset_id TEXT PRIMARY KEY,
    run_id TEXT REFERENCES analysis_runs(run_id),
    path TEXT NOT NULL,
    kind TEXT NOT NULL,
    source_release TEXT,
    sha256 TEXT NOT NULL,
    size_bytes INTEGER,
    UNIQUE (run_id, path, sha256)
);

CREATE TABLE IF NOT EXISTS studies (
    study_id TEXT PRIMARY KEY,
    accession TEXT NOT NULL UNIQUE,
    organism TEXT,
    strain TEXT,
    design TEXT,
    source_archive TEXT
);

CREATE TABLE IF NOT EXISTS samples (
    sample_id TEXT NOT NULL,
    study_id TEXT NOT NULL REFERENCES studies(study_id),
    source_file TEXT,
    condition TEXT NOT NULL,
    antibiotic_role TEXT,
    replicate TEXT,
    metadata_json TEXT,
    PRIMARY KEY (study_id, sample_id)
);

CREATE TABLE IF NOT EXISTS antibiotics (
    antibiotic_id TEXT PRIMARY KEY,
    target_name TEXT NOT NULL,
    observed_compound TEXT,
    antibiotic_class TEXT NOT NULL CHECK (antibiotic_class IN ('beta_lactam', 'aminoglycoside', 'other')),
    proxy_role TEXT NOT NULL CHECK (proxy_role IN ('target', 'proxy', 'challenge', 'control')),
    notes TEXT,
    UNIQUE (target_name, observed_compound)
);

CREATE TABLE IF NOT EXISTS genes (
    gene_id TEXT PRIMARY KEY,
    canonical_name TEXT,
    locus_tag TEXT,
    aliases_json TEXT,
    regulondb_id TEXT UNIQUE
);

CREATE TABLE IF NOT EXISTS regulators (
    regulator_id TEXT PRIMARY KEY,
    name TEXT,
    regulator_type TEXT,
    source_release TEXT
);

CREATE TABLE IF NOT EXISTS transcription_units (
    tu_id TEXT PRIMARY KEY,
    name TEXT,
    source_release TEXT,
    annotation_json TEXT
);

CREATE TABLE IF NOT EXISTS promoters (
    promoter_id TEXT PRIMARY KEY,
    tu_id TEXT REFERENCES transcription_units(tu_id),
    name TEXT,
    genome_accession TEXT,
    start INTEGER,
    end INTEGER,
    strand TEXT CHECK (strand IN ('+', '-', '?') OR strand IS NULL),
    tss INTEGER,
    sigma_factor TEXT,
    sequence TEXT,
    annotation_status TEXT NOT NULL DEFAULT 'review',
    source_release TEXT,
    annotation_json TEXT,
    CHECK (start IS NULL OR start >= 1),
    CHECK (end IS NULL OR end >= 1)
);

CREATE TABLE IF NOT EXISTS regulatory_sites (
    site_id TEXT PRIMARY KEY,
    promoter_id TEXT REFERENCES promoters(promoter_id),
    regulator_id TEXT REFERENCES regulators(regulator_id),
    start INTEGER,
    end INTEGER,
    strand TEXT CHECK (strand IN ('+', '-', '?') OR strand IS NULL),
    sequence TEXT,
    function TEXT,
    evidence TEXT,
    confidence TEXT,
    source_release TEXT,
    annotation_json TEXT
);

CREATE TABLE IF NOT EXISTS contrasts (
    contrast_id TEXT PRIMARY KEY,
    run_id TEXT NOT NULL REFERENCES analysis_runs(run_id),
    antibiotic_id TEXT REFERENCES antibiotics(antibiotic_id),
    name TEXT NOT NULL,
    numerator TEXT NOT NULL,
    denominator TEXT NOT NULL,
    comparison_type TEXT NOT NULL CHECK (comparison_type IN ('drug_control', 'drug_drug')),
    antibiotic_class TEXT,
    proxy_role TEXT,
    effect_cutoff REAL,
    padj_cutoff REAL,
    metadata_json TEXT,
    UNIQUE (run_id, name)
);

CREATE TABLE IF NOT EXISTS de_results (
    run_id TEXT NOT NULL REFERENCES analysis_runs(run_id),
    contrast_id TEXT NOT NULL REFERENCES contrasts(contrast_id),
    gene_id TEXT NOT NULL REFERENCES genes(gene_id),
    base_mean REAL,
    log2_fold_change REAL,
    log2_fold_change_se REAL,
    shrunk_log2_fold_change REAL,
    shrunk_lfc_se REAL,
    p_value REAL,
    padj REAL,
    ci_low REAL,
    ci_high REAL,
    regulation TEXT,
    eligible INTEGER NOT NULL DEFAULT 1 CHECK (eligible IN (0, 1)),
    source_row_id TEXT,
    PRIMARY KEY (run_id, contrast_id, gene_id)
);

CREATE TABLE IF NOT EXISTS regulatory_edges (
    edge_id TEXT PRIMARY KEY,
    run_id TEXT REFERENCES analysis_runs(run_id),
    regulator_id TEXT NOT NULL REFERENCES regulators(regulator_id),
    target_gene_id TEXT NOT NULL REFERENCES genes(gene_id),
    edge_type TEXT NOT NULL,
    effect TEXT,
    evidence TEXT,
    confidence TEXT,
    source_release TEXT,
    metadata_json TEXT,
    UNIQUE (run_id, regulator_id, target_gene_id, edge_type, source_release)
);

CREATE TABLE IF NOT EXISTS tu_genes (
    tu_id TEXT NOT NULL REFERENCES transcription_units(tu_id),
    gene_id TEXT NOT NULL REFERENCES genes(gene_id),
    gene_order INTEGER,
    PRIMARY KEY (tu_id, gene_id)
);

CREATE TABLE IF NOT EXISTS tu_promoters (
    tu_id TEXT NOT NULL REFERENCES transcription_units(tu_id),
    promoter_id TEXT NOT NULL REFERENCES promoters(promoter_id),
    PRIMARY KEY (tu_id, promoter_id)
);

CREATE TABLE IF NOT EXISTS promoter_genes (
    promoter_id TEXT NOT NULL REFERENCES promoters(promoter_id),
    gene_id TEXT NOT NULL REFERENCES genes(gene_id),
    relationship TEXT,
    PRIMARY KEY (promoter_id, gene_id)
);

CREATE TABLE IF NOT EXISTS candidates (
    candidate_id TEXT NOT NULL,
    run_id TEXT NOT NULL REFERENCES analysis_runs(run_id),
    gene_id TEXT NOT NULL REFERENCES genes(gene_id),
    promoter_id TEXT REFERENCES promoters(promoter_id),
    antibiotic_class TEXT NOT NULL,
    support_tier TEXT,
    ranking_score REAL,
    response_summary_json TEXT,
    status TEXT NOT NULL DEFAULT 'proposed',
    review_flags_json TEXT,
    PRIMARY KEY (run_id, candidate_id),
    UNIQUE (run_id, gene_id, promoter_id)
);

CREATE TABLE IF NOT EXISTS candidate_evidence (
    run_id TEXT NOT NULL REFERENCES analysis_runs(run_id),
    candidate_id TEXT NOT NULL,
    contrast_id TEXT NOT NULL REFERENCES contrasts(contrast_id),
    effect REAL,
    padj REAL,
    ci_low REAL,
    ci_high REAL,
    direction TEXT,
    evidence_role TEXT,
    PRIMARY KEY (run_id, candidate_id, contrast_id),
    FOREIGN KEY (run_id, candidate_id) REFERENCES candidates(run_id, candidate_id)
);

CREATE TABLE IF NOT EXISTS candidate_panels (
    run_id TEXT NOT NULL REFERENCES analysis_runs(run_id),
    panel_id TEXT NOT NULL,
    candidate_id TEXT NOT NULL,
    promoter_id TEXT REFERENCES promoters(promoter_id),
    antibiotic_class TEXT NOT NULL,
    selection_order INTEGER,
    selection_rationale TEXT,
    transcription_group TEXT,
    PRIMARY KEY (run_id, panel_id, candidate_id),
    FOREIGN KEY (run_id, candidate_id) REFERENCES candidates(run_id, candidate_id),
    UNIQUE (run_id, panel_id, promoter_id),
    UNIQUE (run_id, panel_id, transcription_group)
);

CREATE TABLE IF NOT EXISTS sample_qc (
    run_id TEXT NOT NULL REFERENCES analysis_runs(run_id),
    study_id TEXT NOT NULL,
    sample_id TEXT NOT NULL,
    library_size REAL,
    genes_detected INTEGER,
    pc1 REAL,
    pc2 REAL,
    distance_summary_json TEXT,
    PRIMARY KEY (run_id, study_id, sample_id),
    FOREIGN KEY (study_id, sample_id) REFERENCES samples(study_id, sample_id)
);

CREATE TABLE IF NOT EXISTS replicate_metrics (
    run_id TEXT NOT NULL REFERENCES analysis_runs(run_id),
    candidate_id TEXT,
    contrast_id TEXT REFERENCES contrasts(contrast_id),
    condition TEXT,
    replicate_count INTEGER,
    mean_value REAL,
    standard_deviation REAL,
    coefficient_of_variation REAL,
    metrics_json TEXT
);

CREATE TABLE IF NOT EXISTS operon_support (
    support_id TEXT PRIMARY KEY,
    run_id TEXT NOT NULL REFERENCES analysis_runs(run_id),
    tu_id TEXT REFERENCES transcription_units(tu_id),
    regulated_set TEXT,
    dominant_group_fraction REAL,
    enrichment REAL,
    padj REAL,
    direction TEXT,
    metadata_json TEXT
);

CREATE TABLE IF NOT EXISTS fragment_reviews (
    fragment_id TEXT PRIMARY KEY,
    run_id TEXT NOT NULL REFERENCES analysis_runs(run_id),
    promoter_id TEXT REFERENCES promoters(promoter_id),
    genome_accession TEXT,
    start INTEGER,
    end INTEGER,
    strand TEXT CHECK (strand IN ('+', '-', '?') OR strand IS NULL),
    sequence TEXT,
    sequence_match INTEGER CHECK (sequence_match IN (0, 1) OR sequence_match IS NULL),
    missing_coordinates INTEGER NOT NULL DEFAULT 0 CHECK (missing_coordinates IN (0, 1)),
    regulatory_site_count INTEGER,
    overlapping_genes_json TEXT,
    review_status TEXT NOT NULL DEFAULT 'review',
    metadata_json TEXT
);

CREATE INDEX IF NOT EXISTS idx_assets_run ON assets(run_id);
CREATE INDEX IF NOT EXISTS idx_samples_condition ON samples(condition);
CREATE INDEX IF NOT EXISTS idx_promoters_tu ON promoters(tu_id);
CREATE INDEX IF NOT EXISTS idx_sites_promoter ON regulatory_sites(promoter_id);
CREATE INDEX IF NOT EXISTS idx_sites_regulator ON regulatory_sites(regulator_id);
CREATE INDEX IF NOT EXISTS idx_contrasts_run ON contrasts(run_id);
CREATE INDEX IF NOT EXISTS idx_de_gene ON de_results(gene_id, contrast_id);
CREATE INDEX IF NOT EXISTS idx_de_padj ON de_results(contrast_id, padj);
CREATE INDEX IF NOT EXISTS idx_edges_regulator ON regulatory_edges(regulator_id);
CREATE INDEX IF NOT EXISTS idx_edges_target ON regulatory_edges(target_gene_id);
CREATE INDEX IF NOT EXISTS idx_tu_genes_gene ON tu_genes(gene_id);
CREATE INDEX IF NOT EXISTS idx_tu_promoters_promoter ON tu_promoters(promoter_id);
CREATE INDEX IF NOT EXISTS idx_candidates_class ON candidates(antibiotic_class, support_tier);
CREATE INDEX IF NOT EXISTS idx_candidates_gene ON candidates(gene_id);
CREATE INDEX IF NOT EXISTS idx_evidence_contrast ON candidate_evidence(contrast_id, padj);
CREATE INDEX IF NOT EXISTS idx_panel_class ON candidate_panels(panel_id, antibiotic_class);
CREATE INDEX IF NOT EXISTS idx_qc_sample ON sample_qc(sample_id);
CREATE INDEX IF NOT EXISTS idx_replicates_candidate ON replicate_metrics(candidate_id, contrast_id);
CREATE INDEX IF NOT EXISTS idx_operon_run ON operon_support(run_id, tu_id);
CREATE INDEX IF NOT EXISTS idx_fragments_promoter ON fragment_reviews(promoter_id);
"""


def open_database(path: str | Path = DEFAULT_DATABASE) -> sqlite3.Connection:
    """Open a SQLite database with project safety settings enabled."""
    database_path = Path(path)
    database_path.parent.mkdir(parents=True, exist_ok=True)
    connection = sqlite3.connect(database_path)
    connection.row_factory = sqlite3.Row
    connection.execute("PRAGMA foreign_keys = ON")
    connection.execute("PRAGMA journal_mode = WAL")
    initialize_schema(connection)
    return connection


def initialize_schema(connection: sqlite3.Connection) -> None:
    """Create the current schema without deleting existing records."""
    connection.executescript(SCHEMA_SQL)
    connection.execute(
        "INSERT OR REPLACE INTO database_metadata(key, value) VALUES (?, ?)",
        ("schema_version", str(SCHEMA_VERSION)),
    )
    connection.commit()


def schema_version(connection: sqlite3.Connection) -> int:
    """Return the stored schema version, or zero for an uninitialized DB."""
    row = connection.execute(
        "SELECT value FROM database_metadata WHERE key = 'schema_version'"
    ).fetchone()
    return int(row["value"]) if row else 0


@contextmanager
def transaction(connection: sqlite3.Connection) -> Iterator[sqlite3.Connection]:
    """Run a group of writes atomically and roll back on failure."""
    try:
        yield connection
    except Exception:
        connection.rollback()
        raise
    else:
        connection.commit()


def close_database(connection: sqlite3.Connection) -> None:
    """Commit pending work and close a database connection."""
    connection.commit()
    connection.close()
