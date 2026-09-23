"""Small read-oriented Python API for the generated project database."""

from __future__ import annotations

import sqlite3
import hashlib
from datetime import datetime, timezone
from pathlib import Path

from .database import DEFAULT_DATABASE, open_database


class Database:
    """Convenience queries that return dictionaries suitable for scripts or JSON."""

    def __init__(self, path: str | Path = DEFAULT_DATABASE):
        self.connection = open_database(path)

    def close(self) -> None:
        self.connection.close()

    def _rows(self, query: str, parameters=()) -> list[dict]:
        return [dict(row) for row in self.connection.execute(query, parameters).fetchall()]

    def search_candidates(
        self,
        antibiotic_class: str | None = None,
        support_tier: str | None = None,
        padj_max: float | None = None,
        effect_min: float | None = None,
        run_id: str | None = None,
    ) -> list[dict]:
        clauses = []
        parameters: list[object] = []
        if antibiotic_class:
            clauses.append("c.antibiotic_class = ?")
            parameters.append(antibiotic_class)
        if support_tier:
            clauses.append("c.support_tier = ?")
            parameters.append(support_tier)
        if run_id:
            clauses.append("c.run_id = ?")
            parameters.append(run_id)
        if padj_max is not None or effect_min is not None:
            clauses.append("EXISTS (SELECT 1 FROM candidate_evidence e WHERE e.run_id = c.run_id AND e.candidate_id = c.candidate_id")
            if padj_max is not None:
                clauses.append("e.padj <= ?")
                parameters.append(padj_max)
            if effect_min is not None:
                clauses.append("ABS(e.effect) >= ?")
                parameters.append(effect_min)
            clauses.append(")")
        where = " WHERE " + " AND ".join(clauses) if clauses else ""
        return self._rows("SELECT * FROM candidate_summary c" + where + " ORDER BY c.ranking_score DESC NULLS LAST", parameters)

    def get_candidate(self, candidate_id: str, run_id: str | None = None) -> dict | None:
        rows = self._rows("SELECT * FROM candidate_summary WHERE candidate_id = ?" + (" AND run_id = ?" if run_id else "") + " LIMIT 1", (candidate_id, run_id) if run_id else (candidate_id,))
        if not rows:
            return None
        result = rows[0]
        result["evidence"] = self._rows("SELECT * FROM candidate_evidence WHERE run_id = ? AND candidate_id = ? ORDER BY padj", (result["run_id"], candidate_id))
        result["regulators"] = self._rows("SELECT DISTINCT r.* FROM regulatory_edges e JOIN regulators r ON r.regulator_id = e.regulator_id WHERE e.run_id = ? AND e.target_gene_id = ?", (result["run_id"], result["gene_id"]))
        return result

    def search_regulator(self, regulator: str) -> list[dict]:
        return self._rows("SELECT * FROM regulator_candidate_paths WHERE regulator_id = ? OR regulator_name = ? ORDER BY antibiotic_class, candidate_id", (regulator, regulator))

    def panel(self, panel_id: str = "experimental_panel", run_id: str | None = None) -> list[dict]:
        query = "SELECT * FROM panel_review WHERE panel_id = ?"
        parameters: list[object] = [panel_id]
        if run_id:
            query += " AND run_id = ?"
            parameters.append(run_id)
        return self._rows(query + " ORDER BY selection_order", parameters)

    def sample_qc(self, run_id: str | None = None) -> list[dict]:
        if run_id:
            return self._rows("SELECT * FROM sample_qc_summary WHERE run_id = ? ORDER BY sample_id", (run_id,))
        return self._rows("SELECT * FROM sample_qc_summary ORDER BY run_id, sample_id")

    def record_candidate_review(
        self,
        run_id: str,
        candidate_id: str,
        review_status: str,
        editor: str,
        reason: str,
        priority: int | None = None,
        construct_ready: bool | None = None,
        notes: str | None = None,
        source_reference: str | None = None,
    ) -> str:
        """Append a validated candidate review and return its edit ID."""
        valid = {"proposed", "reviewed", "approved", "excluded"}
        if review_status not in valid:
            raise ValueError(f"review_status must be one of {sorted(valid)}")
        exists = self.connection.execute("SELECT 1 FROM candidates WHERE run_id = ? AND candidate_id = ?", (run_id, candidate_id)).fetchone()
        if not exists:
            raise ValueError(f"Unknown candidate {candidate_id!r} for run {run_id!r}")
        edited_at = datetime.now(timezone.utc).isoformat()
        edit_id = "edit-" + hashlib.sha256(f"{run_id}:{candidate_id}:{edited_at}:{editor}".encode()).hexdigest()[:20]
        previous = self.connection.execute("SELECT review_status FROM candidate_reviews WHERE run_id = ? AND candidate_id = ? ORDER BY edit_id DESC LIMIT 1", (run_id, candidate_id)).fetchone()
        with self.connection:
            self.connection.execute("INSERT INTO curation_edits(edit_id, run_id, entity_type, entity_id, field_name, old_value, new_value, editor, edited_at, reason, source_reference) VALUES (?, ?, 'candidate', ?, 'review_status', ?, ?, ?, ?, ?, ?)", (edit_id, run_id, candidate_id, previous[0] if previous else None, review_status, editor, edited_at, reason, source_reference))
            self.connection.execute("INSERT INTO candidate_reviews(edit_id, run_id, candidate_id, review_status, priority, construct_ready, notes) VALUES (?, ?, ?, ?, ?, ?, ?)", (edit_id, run_id, candidate_id, review_status, priority, None if construct_ready is None else int(construct_ready), notes))
        return edit_id
