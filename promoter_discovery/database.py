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


SCHEMA_VERSION = 2
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

CREATE INDEX IF NOT EXISTS idx_assets_run ON assets(run_id);
CREATE INDEX IF NOT EXISTS idx_samples_condition ON samples(condition);
CREATE INDEX IF NOT EXISTS idx_promoters_tu ON promoters(tu_id);
CREATE INDEX IF NOT EXISTS idx_sites_promoter ON regulatory_sites(promoter_id);
CREATE INDEX IF NOT EXISTS idx_sites_regulator ON regulatory_sites(regulator_id);
CREATE INDEX IF NOT EXISTS idx_contrasts_run ON contrasts(run_id);
CREATE INDEX IF NOT EXISTS idx_de_gene ON de_results(gene_id, contrast_id);
CREATE INDEX IF NOT EXISTS idx_de_padj ON de_results(contrast_id, padj);
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
