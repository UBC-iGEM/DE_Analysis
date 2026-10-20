"""Content fingerprints and retained evidence snapshots for local releases."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import sqlite3
import subprocess
import tempfile


def source_provenance() -> dict:
    """Fingerprint the executing source, including uncommitted Python changes.

    Git is descriptive; content hashes determine identity even in an exported
    checkout without Git. Documentation-only commits do not invalidate reviews.
    """
    root = Path(__file__).resolve().parents[1]
    paths = sorted((root / "promoter_discovery").glob("*.py"))
    for name in ("requirements.txt", "requirements-lock.txt"):
        path = root / name
        if path.is_file():
            paths.append(path)
    manifest = {}
    for path in paths:
        relative_path = path.relative_to(root).as_posix()
        manifest[relative_path] = hashlib.sha256(path.read_bytes()).hexdigest()
    digest = hashlib.sha256(json.dumps(manifest, sort_keys=True).encode()).hexdigest()
    revision, dirty = None, None
    try:
        revision = subprocess.run(
            ["git", "rev-parse", "HEAD"], cwd=root, check=True, capture_output=True,
            text=True, timeout=5,
        ).stdout.strip()
        status = subprocess.run(
            ["git", "status", "--porcelain", "--", "promoter_discovery", "requirements.txt", "requirements-lock.txt"],
            cwd=root, check=True, capture_output=True, text=True, timeout=5,
        ).stdout
        dirty = bool(status.strip())
    except (OSError, subprocess.SubprocessError):
        pass
    return {"git_revision": revision, "code_dirty": dirty,
            "code_sha256": digest, "code_manifest": manifest}


def preserve_snapshot(path: Path, new_run_id: str, new_schema: int) -> Path | None:
    """Back up a superseded scientific snapshot, including committed WAL data.

    Failure aborts replacement of the live database. A content-addressed name
    avoids overwriting an earlier snapshot or creating duplicate backups.
    """
    if not path.exists():
        return None
    source = sqlite3.connect(path.resolve().as_uri() + "?mode=ro", uri=True)
    temporary = None
    try:
        old_schema = source.execute("SELECT value FROM database_metadata WHERE key='schema_version'").fetchone()
        runs = [row[0] for row in source.execute("SELECT run_id FROM analysis_runs ORDER BY run_id")]
        if old_schema and int(old_schema[0]) == new_schema and runs == [new_run_id]:
            return None
        folder = path.parent / "database_snapshots" / path.stem
        folder.mkdir(parents=True, exist_ok=True)
        descriptor, name = tempfile.mkstemp(prefix="snapshot-", suffix=".sqlite", dir=folder)
        os.close(descriptor)
        temporary = Path(name)
        destination = sqlite3.connect(temporary)
        try:
            source.backup(destination)
            destination.execute("PRAGMA journal_mode=DELETE")
            if destination.execute("PRAGMA integrity_check").fetchone()[0] != "ok":
                raise ValueError("Previous database snapshot failed integrity validation")
        finally:
            destination.close()
        digest = hashlib.sha256(temporary.read_bytes()).hexdigest()
        target = folder / f"sha256-{digest}.sqlite"
        if target.exists():
            if hashlib.sha256(target.read_bytes()).hexdigest() != digest:
                raise ValueError(f"Existing evidence snapshot has changed: {target}")
        else:
            os.replace(temporary, target)
        return target
    finally:
        source.close()
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def prepare_replacement(path: Path) -> None:
    """Remove WAL state through SQLite before replacing its main file.

    Rebuilds require a quiescent database. Fail if a connection prevents the
    journal-mode transition; never unlink journal files manually.
    """
    if not path.exists():
        return
    connection = sqlite3.connect(path.resolve().as_uri() + "?mode=rw", uri=True, timeout=1)
    try:
        if connection.execute("PRAGMA journal_mode").fetchone()[0] == "wal":
            busy, _, _ = connection.execute("PRAGMA wal_checkpoint(TRUNCATE)").fetchone()
            if busy:
                raise RuntimeError("Database is in use; stop curation/explorer connections before rebuilding")
            if connection.execute("PRAGMA journal_mode=DELETE").fetchone()[0] != "delete":
                raise RuntimeError("Could not clear the previous database's WAL state")
    finally:
        connection.close()
