"""Persist user reviews independently of the rebuildable scientific database."""

from __future__ import annotations

import json
import sqlite3
from pathlib import Path


ARCHIVE_SCHEMA = """CREATE TABLE IF NOT EXISTS review_history (
    edit_id TEXT PRIMARY KEY,
    run_id TEXT NOT NULL,
    reference_release TEXT,
    entity_type TEXT NOT NULL,
    entity_id TEXT NOT NULL,
    edited_at TEXT NOT NULL,
    edit_json TEXT NOT NULL,
    review_json TEXT NOT NULL
)"""


def archive_path(database_path: str | Path) -> Path:
    return Path(database_path).with_suffix(".curation.sqlite")


def attach_archive(connection: sqlite3.Connection, path: Path) -> None:
    if "curation_archive" not in {row[1] for row in connection.execute("PRAGMA database_list")}:
        connection.execute("ATTACH DATABASE ? AS curation_archive", (str(path),))
        connection.execute(ARCHIVE_SCHEMA.replace("review_history", "curation_archive.review_history", 1))


def archive_review(connection: sqlite3.Connection, edit: dict, review: dict) -> None:
    release = connection.execute(
        "SELECT reference_release FROM analysis_runs WHERE run_id = ?", (edit["run_id"],)
    ).fetchone()[0]
    connection.execute(
        "INSERT OR IGNORE INTO curation_archive.review_history VALUES (?, ?, ?, ?, ?, ?, ?, ?)",
        (edit["edit_id"], edit["run_id"], release, edit["entity_type"], edit["entity_id"],
         edit["edited_at"], json.dumps(edit), json.dumps(review)),
    )


def preserve_existing_reviews(database_path: Path) -> None:
    """Capture reviews from older generated databases, including schema 9."""
    if not database_path.exists():
        return
    old = sqlite3.connect(f"{database_path.as_uri()}?mode=ro", uri=True)
    old.row_factory = sqlite3.Row
    try:
        tables = {row[0] for row in old.execute("SELECT name FROM sqlite_master WHERE type='table'")}
        if "curation_edits" not in tables:
            return
        edits = [dict(row) for row in old.execute("SELECT * FROM curation_edits")]
        if not edits:
            return
        with sqlite3.connect(archive_path(database_path)) as archive:
            archive.execute(ARCHIVE_SCHEMA)
            for edit in edits:
                table = "candidate_reviews" if edit["entity_type"] == "candidate" else "promoter_reviews"
                review = old.execute(f"SELECT * FROM {table} WHERE edit_id = ?", (edit["edit_id"],)).fetchone()
                if review is None or edit["run_id"] is None:
                    raise ValueError(f"Cannot preserve incomplete review {edit['edit_id']}")
                release = old.execute("SELECT reference_release FROM analysis_runs WHERE run_id = ?", (edit["run_id"],)).fetchone()[0]
                archive.execute(
                    "INSERT OR IGNORE INTO review_history VALUES (?, ?, ?, ?, ?, ?, ?, ?)",
                    (edit["edit_id"], edit["run_id"], release, edit["entity_type"], edit["entity_id"],
                     edit["edited_at"], json.dumps(edit), json.dumps(dict(review))),
                )
    finally:
        old.close()


def restore_reviews(connection: sqlite3.Connection, path: Path, run_id: str) -> None:
    """Replay matching-run annotations; retain older-run history in the archive."""
    if not path.exists():
        return
    with sqlite3.connect(f"{path.as_uri()}?mode=ro", uri=True) as archive:
        records = archive.execute(
            "SELECT edit_json, review_json FROM review_history WHERE run_id=? ORDER BY edited_at, edit_id",
            (run_id,),
        ).fetchall()
    for edit_json, review_json in records:
        edit, review = json.loads(edit_json), json.loads(review_json)
        for table, values in (("curation_edits", edit),
                              ("candidate_reviews" if edit["entity_type"] == "candidate" else "promoter_reviews", review)):
            columns = list(values)
            connection.execute(
                f"INSERT INTO {table} ({', '.join(columns)}) VALUES ({', '.join('?' for _ in columns)})",
                [values[column] for column in columns],
            )


def review_history(path: Path, entity_type: str, entity_id: str) -> list[dict]:
    if not path.exists():
        return []
    with sqlite3.connect(f"{path.as_uri()}?mode=ro", uri=True) as archive:
        rows = archive.execute(
            "SELECT edit_json, review_json, reference_release FROM review_history "
            "WHERE entity_type=? AND entity_id=? ORDER BY edited_at, edit_id", (entity_type, entity_id)
        ).fetchall()
    return [{**json.loads(edit), "review": json.loads(review), "reference_release": release}
            for edit, review, release in rows]
