"""Small read-oriented Python API for the generated project database."""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from uuid import uuid4

from .database import DEFAULT_DATABASE, open_database
from .curation import archive_path, archive_review, attach_archive, review_history


class Database:
    """Convenience queries that return dictionaries suitable for scripts or JSON."""

    def __init__(self, path: str | Path = DEFAULT_DATABASE):
        self.path = Path(path).resolve()
        if not self.path.is_file():
            raise FileNotFoundError(f"Database not found: {self.path}. Run python -m promoter_discovery.build_database")
        self.connection = open_database(self.path)

    def close(self) -> None:
        self.connection.close()

    def _rows(self, query: str, parameters=()) -> list[dict]:
        return [dict(row) for row in self.connection.execute(query, parameters).fetchall()]

    def _run_id(self, run_id: str | None = None) -> str | None:
        if run_id is not None:
            return run_id
        row = self.connection.execute(
            "SELECT run_id FROM analysis_runs WHERE status='complete' ORDER BY created_at DESC, run_id DESC LIMIT 1"
        ).fetchone()
        return row[0] if row else None

    def runs(self) -> list[dict]:
        return self._rows("SELECT * FROM analysis_runs ORDER BY created_at DESC, run_id DESC")

    def reference_coverage(self, run_id: str | None = None) -> dict:
        row = self.connection.execute("SELECT summary_json FROM analysis_runs WHERE run_id=?", (self._run_id(run_id),)).fetchone()
        return json.loads(row[0] or "{}").get("reference_coverage", {}) if row else {}

    def search_candidates(
        self,
        antibiotic_class: str | None = None,
        support_tier: str | None = None,
        padj_max: float | None = None,
        effect_min: float | None = None,
        run_id: str | None = None,
        query: str | None = None,
    ) -> list[dict]:
        if antibiotic_class is not None and antibiotic_class not in {"beta_lactam", "aminoglycoside"}:
            raise ValueError("class must be beta_lactam or aminoglycoside")
        if padj_max is not None and not 0 <= padj_max <= 1:
            raise ValueError("padj_max must be between 0 and 1")
        if effect_min is not None and effect_min < 0:
            raise ValueError("effect_min must be nonnegative")
        clauses = []
        parameters: list[object] = []
        if antibiotic_class:
            clauses.append("c.antibiotic_class = ?")
            parameters.append(antibiotic_class)
        if support_tier:
            clauses.append("c.support_tier = ?")
            parameters.append(support_tier)
        clauses.append("c.run_id = ?")
        parameters.append(self._run_id(run_id))
        if query:
            clauses.append("(instr(lower(c.candidate_id), lower(?)) > 0 OR instr(lower(c.canonical_name), lower(?)) > 0)")
            parameters.extend((query, query))
        if padj_max is not None or effect_min is not None:
            evidence_clauses = [
                "e.run_id = c.run_id",
                "e.candidate_id = c.candidate_id",
                "e.evidence_role = 'discovery'",
            ]
            if padj_max is not None:
                evidence_clauses.append("e.padj <= ?")
                parameters.append(padj_max)
            if effect_min is not None:
                evidence_clauses.append("ABS(e.effect) >= ?")
                parameters.append(effect_min)
            clauses.append(
                "EXISTS (SELECT 1 FROM candidate_evidence e WHERE "
                + " AND ".join(evidence_clauses) + ")"
            )
        where = " WHERE " + " AND ".join(clauses) if clauses else ""
        return self._rows("SELECT * FROM candidate_summary c" + where + " ORDER BY c.ranking_score DESC NULLS LAST", parameters)

    def get_candidate(self, candidate_id: str, run_id: str | None = None) -> dict | None:
        rows = self._rows("SELECT * FROM candidate_summary WHERE lower(candidate_id) = lower(?) AND run_id = ?", (candidate_id, self._run_id(run_id)))
        if not rows:
            return None
        result = rows[0]
        candidate_id = result["candidate_id"]
        result["evidence"] = self._rows("SELECT * FROM candidate_evidence WHERE run_id = ? AND candidate_id = ? ORDER BY padj", (result["run_id"], candidate_id))
        result["regulators"] = self._rows("SELECT DISTINCT r.* FROM regulatory_edges e JOIN regulators r ON r.regulator_id = e.regulator_id WHERE e.run_id = ? AND e.target_gene_id = ?", (result["run_id"], result["gene_id"]))
        result["transcription_units"] = self._rows(
            "SELECT tu.* FROM candidate_tus ct JOIN transcription_units tu ON tu.tu_id = ct.tu_id "
            "WHERE ct.run_id = ? AND ct.candidate_id = ? ORDER BY tu.tu_id",
            (result["run_id"], candidate_id),
        )
        result["promoters"] = self._rows(
            "SELECT p.* FROM candidate_promoters cp JOIN promoters p ON p.promoter_id = cp.promoter_id "
            "WHERE cp.run_id = ? AND cp.candidate_id = ? ORDER BY p.promoter_id",
            (result["run_id"], candidate_id),
        )
        result["regulatory_paths"] = self.candidate_paths(candidate_id, result["run_id"])
        return result

    def search_regulator(self, regulator: str, run_id: str | None = None) -> list[dict]:
        return self._rows("SELECT * FROM regulator_candidate_paths WHERE (lower(regulator_id) = lower(?) OR lower(regulator_name) = lower(?)) AND run_id=? ORDER BY antibiotic_class, candidate_id", (regulator, regulator, self._run_id(run_id)))

    def reference_interactions(
        self, actor: str, target_kind: str | None = None, run_id: str | None = None
    ) -> list[dict]:
        """Return individual curated interactions for an actor name or source ID."""
        if target_kind is not None and target_kind not in {"promoter", "tu", "gene"}:
            raise ValueError("target_kind must be promoter, tu, or gene")
        query = (
            "SELECT ri.interaction_id, ri.run_id, a.actor_id, a.name AS actor_name, a.actor_type, "
            "ri.interaction_type, ri.target_kind, "
            "COALESCE(ri.promoter_id, ri.tu_id, ri.gene_id) AS target_id, "
            "ri.site_id, ri.effect, ri.conformation, ri.confidence, ri.site_evidence, "
            "ri.interaction_evidence, ri.evidence_category, ri.pmids, ri.source_release "
            "FROM regulatory_interactions ri JOIN regulatory_actors a ON a.actor_id = ri.actor_id "
            "WHERE (lower(a.actor_id) = lower(?) OR lower(a.name) = lower(?))"
        )
        parameters: list[object] = [actor, actor]
        if target_kind:
            query += " AND ri.target_kind = ?"
            parameters.append(target_kind)
        query += " AND ri.run_id = ?"
        parameters.append(self._run_id(run_id))
        return self._rows(query + " ORDER BY ri.target_kind, target_id, ri.interaction_id", parameters)

    def reference_promoter_paths(self, actor: str, run_id: str | None = None) -> list[dict]:
        """Return all target types, with direct regulation distinguished from context."""
        query = "SELECT * FROM reference_promoter_paths WHERE (lower(actor_id) = lower(?) OR lower(actor_name) = lower(?))"
        parameters: list[object] = [actor, actor]
        query += " AND run_id = ?"
        parameters.append(self._run_id(run_id))
        return self._rows(query + " ORDER BY promoter_id, tu_id, gene_id, interaction_id", parameters)

    def candidate_paths(self, candidate_id: str, run_id: str | None = None) -> list[dict]:
        """Gene-target effects apply to that gene, never to siblings in its TU."""
        run_id = self._run_id(run_id)
        candidate = self.connection.execute(
            "SELECT candidate_id, gene_id FROM candidates WHERE lower(candidate_id)=lower(?) AND run_id=?",
            (candidate_id, run_id),
        ).fetchone()
        if candidate is None:
            return []
        return self._rows(
            "SELECT DISTINCT p.* FROM reference_regulatory_paths p WHERE p.run_id=? "
            "AND p.interaction_id IN (SELECT ri.interaction_id FROM regulatory_interactions ri WHERE ri.run_id=? AND "
            "(ri.gene_id=? OR ri.promoter_id IN (SELECT cp.promoter_id FROM candidate_promoters cp WHERE cp.run_id=? "
            "AND cp.candidate_id=? UNION SELECT tp.promoter_id FROM tu_promoters tp JOIN tu_genes tg "
            "ON tg.tu_id=tp.tu_id WHERE tg.gene_id=?) OR ri.tu_id IN (SELECT ct.tu_id FROM candidate_tus ct "
            "WHERE ct.run_id=? AND ct.candidate_id=? UNION SELECT tu_id FROM tu_genes WHERE gene_id=?))) "
            "ORDER BY p.actor_name, p.interaction_id, p.promoter_id, p.tu_id, p.gene_id",
            (run_id, run_id, candidate["gene_id"], run_id, candidate["candidate_id"], candidate["gene_id"],
             run_id, candidate["candidate_id"], candidate["gene_id"]),
        )

    def get_promoter(self, promoter_id: str, run_id: str | None = None) -> dict | None:
        run_id = self._run_id(run_id)
        rows = self._rows("SELECT * FROM promoters WHERE promoter_id=?", (promoter_id,))
        if not rows:
            return None
        result = rows[0]
        result["regulatory_paths"] = self._rows(
            "SELECT * FROM reference_regulatory_paths WHERE promoter_id=? AND run_id=? ORDER BY actor_name, interaction_id",
            (promoter_id, run_id),
        )
        result["fragments"] = self._rows("SELECT * FROM fragment_reviews WHERE promoter_id=? AND run_id=?", (promoter_id, run_id))
        result["review"] = self._latest_review("promoter", promoter_id, run_id)
        return result

    def candidate_network(self, candidate_id: str, run_id: str | None = None) -> dict | None:
        """Export a local graph with evidence on each individual regulation edge."""
        candidate = self.get_candidate(candidate_id, run_id)
        if candidate is None:
            return None
        nodes, edges = {}, {}

        def node(kind, identifier, label=None):
            key = f"{kind}:{identifier}"
            nodes[key] = {"id": key, "entity_id": identifier, "kind": kind, "label": label or identifier}
            return key

        gene = node("gene", candidate["gene_id"], candidate["canonical_name"])
        for path in candidate["regulatory_paths"]:
            actor = node("actor", path["actor_id"], path["actor_name"])
            target = node(path["target_kind"], path["target_id"])
            edges[path["interaction_id"]] = {
                "id": path["interaction_id"], "source": actor, "target": target,
                "kind": "regulation", "effect": path["effect"], "confidence": path["confidence"],
                "evidence": path["interaction_evidence"], "site_evidence": path["site_evidence"],
                "site_id": path["site_id"], "conformation": path["conformation"],
                "pmids": path["pmids"], "source_release": path["source_release"],
                "target_kind": path["target_kind"],
            }
            promoter = node("promoter", path["promoter_id"], path["promoter_name"]) if path["promoter_id"] else None
            tu = node("tu", path["tu_id"]) if path["tu_id"] else None
            for source, target, relation in ((promoter, tu, "promoter–TU association"),
                                              (tu, gene if path["gene_id"] == candidate["gene_id"] else None, "TU membership")):
                if source and target:
                    key = f"context:{source}:{target}"
                    edges[key] = {"id": key, "source": source, "target": target,
                                  "kind": "context", "label": relation}
        return {"candidate_id": candidate["candidate_id"], "run_id": candidate["run_id"],
                "nodes": list(nodes.values()), "edges": list(edges.values())}

    def panel(self, panel_id: str = "experimental_panel", run_id: str | None = None) -> list[dict]:
        query = "SELECT * FROM panel_review WHERE panel_id = ?"
        parameters: list[object] = [panel_id]
        query += " AND run_id = ?"
        parameters.append(self._run_id(run_id))
        return self._rows(query + " ORDER BY selection_order", parameters)

    def sample_qc(self, run_id: str | None = None) -> list[dict]:
        return self._rows("SELECT * FROM sample_qc_summary WHERE run_id = ? ORDER BY sample_id", (self._run_id(run_id),))

    def _latest_review(self, entity_type: str, entity_id: str, run_id: str | None) -> dict | None:
        table = "candidate_reviews" if entity_type == "candidate" else "promoter_reviews"
        rows = self._rows(
            f"SELECT r.* FROM {table} r JOIN curation_edits e ON e.edit_id=r.edit_id "
            "WHERE e.entity_id=? AND e.run_id=? ORDER BY e.edited_at DESC, e.edit_id DESC LIMIT 1",
            (entity_id, run_id),
        )
        return rows[0] if rows else None

    def review_history(self, entity_type: str, entity_id: str) -> list[dict]:
        if entity_type not in {"candidate", "promoter"}:
            raise ValueError("entity_type must be candidate or promoter")
        return review_history(archive_path(self.path), entity_type, entity_id)

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
        if priority is not None and (type(priority) is not int or priority < 0):
            raise ValueError("priority must be a nonnegative integer")
        return self._record_review("candidate", candidate_id, run_id, review_status, editor, reason,
                                   construct_ready, notes, source_reference, priority)

    def record_promoter_review(
        self, run_id: str, promoter_id: str, review_status: str, editor: str, reason: str,
        construct_ready: bool | None = None, notes: str | None = None, source_reference: str | None = None,
    ) -> str:
        return self._record_review("promoter", promoter_id, run_id, review_status, editor, reason,
                                   construct_ready, notes, source_reference)

    def _record_review(self, entity_type, entity_id, run_id, status, editor, reason,
                       construct_ready, notes, source_reference, priority=None) -> str:
        valid = ({"proposed", "reviewed", "approved", "excluded"} if entity_type == "candidate"
                 else {"review", "validated", "ready_for_order", "excluded"})
        if status not in valid:
            raise ValueError(f"review_status must be one of {sorted(valid)}")
        if not editor.strip() or not reason.strip():
            raise ValueError("editor and reason are required")
        if construct_ready is not None and type(construct_ready) is not bool:
            raise ValueError("construct_ready must be a boolean")
        if not self.connection.execute("SELECT 1 FROM analysis_runs WHERE run_id=? AND status='complete'", (run_id,)).fetchone():
            raise ValueError(f"Unknown completed run {run_id!r}")
        if entity_type == "candidate":
            exists = self.connection.execute("SELECT 1 FROM candidates WHERE run_id=? AND candidate_id=?", (run_id, entity_id)).fetchone()
        else:
            exists = self.connection.execute("SELECT 1 FROM promoters WHERE promoter_id=?", (entity_id,)).fetchone()
        if not exists:
            raise ValueError(f"Unknown {entity_type} {entity_id!r} for run {run_id!r}")
        edited_at = datetime.now(timezone.utc).isoformat()
        edit_id = "edit-" + uuid4().hex
        previous = self._latest_review(entity_type, entity_id, run_id)
        review = {"edit_id": edit_id, f"{entity_type}_id": entity_id, "review_status": status,
                  "construct_ready": None if construct_ready is None else int(construct_ready), "notes": notes}
        if entity_type == "candidate":
            review.update(run_id=run_id, priority=priority)
        edit = {"edit_id": edit_id, "run_id": run_id, "entity_type": entity_type, "entity_id": entity_id,
                "field_name": "review", "old_value": json.dumps(previous) if previous else None,
                "new_value": json.dumps(review), "editor": editor, "edited_at": edited_at,
                "reason": reason, "source_reference": source_reference}
        attach_archive(self.connection, archive_path(self.path))
        with self.connection:
            for table, values in (("curation_edits", edit), (f"{entity_type}_reviews", review)):
                self.connection.execute(
                    f"INSERT INTO {table} ({', '.join(values)}) VALUES ({', '.join('?' for _ in values)})",
                    list(values.values()),
                )
            archive_review(self.connection, edit, review)
        return edit_id
