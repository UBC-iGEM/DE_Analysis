"""Command-line searches for the generated promoter-discovery database."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from .database import DEFAULT_DATABASE
from .db_api import Database


def _print_rows(rows: list[dict], as_json: bool) -> None:
    if as_json:
        print(json.dumps(rows, indent=2, sort_keys=True, default=str))
        return
    if not rows:
        print("No matching records.")
        return
    columns = list(rows[0])
    print("\t".join(columns))
    for row in rows:
        print("\t".join("" if row[column] is None else str(row[column]) for column in columns))


def main(argv=None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--database", type=Path, default=DEFAULT_DATABASE)
    parser.add_argument("--json", action="store_true", dest="as_json", help="write JSON instead of a table")
    subparsers = parser.add_subparsers(dest="command", required=True)

    candidates = subparsers.add_parser("candidates", help="search candidate promoters")
    candidates.add_argument("--class", dest="antibiotic_class")
    candidates.add_argument("--support-tier")
    candidates.add_argument("--padj-max", type=float)
    candidates.add_argument("--effect-min", type=float)
    candidates.add_argument("--run-id")
    candidates.add_argument("--query", help="match a candidate ID or gene name")

    candidate = subparsers.add_parser("candidate", help="show one candidate and its evidence")
    candidate.add_argument("candidate_id")
    candidate.add_argument("--run-id")

    regulator = subparsers.add_parser("regulator", help="show regulator-to-candidate paths")
    regulator.add_argument("regulator")
    regulator.add_argument("--run-id")

    interactions = subparsers.add_parser("interactions", help="show curated reference interactions for an actor")
    interactions.add_argument("actor")
    interactions.add_argument("--target-kind", choices=("promoter", "tu", "gene"))
    interactions.add_argument("--run-id")

    paths = subparsers.add_parser("promoter-paths", help="show reference actor-to-promoter-to-gene paths")
    paths.add_argument("actor")
    paths.add_argument("--run-id")

    panel = subparsers.add_parser("panel", help="show a selected panel")
    panel.add_argument("--panel-id", default="experimental_panel")
    panel.add_argument("--run-id")

    qc = subparsers.add_parser("qc", help="show sample quality-control records")
    qc.add_argument("--run-id")

    coverage = subparsers.add_parser("coverage", help="show reference import and mapping coverage")
    coverage.add_argument("--run-id")
    network = subparsers.add_parser("network", help="export a candidate graph with direct and context edges")
    network.add_argument("candidate_id")
    network.add_argument("--run-id")
    promoter = subparsers.add_parser("promoter", help="show a promoter, fragment checks, and review")
    promoter.add_argument("promoter_id")
    promoter.add_argument("--run-id")
    history = subparsers.add_parser("history", help="show persistent reviews across analysis runs")
    history.add_argument("entity_type", choices=("candidate", "promoter"))
    history.add_argument("entity_id")
    for kind in ("candidate", "promoter"):
        review = subparsers.add_parser(f"review-{kind}", help=f"append a versioned {kind} review")
        review.add_argument("entity_id")
        review.add_argument("status")
        review.add_argument("--run-id", required=True)
        review.add_argument("--editor", required=True)
        review.add_argument("--reason", required=True)
        review.add_argument("--notes")
        review.add_argument("--source-reference")
        review.add_argument("--construct-ready", action=argparse.BooleanOptionalAction, default=None)
        if kind == "candidate":
            review.add_argument("--priority", type=int)

    args = parser.parse_args(argv)
    try:
        database = Database(args.database)
    except (FileNotFoundError, ValueError) as error:
        parser.error(str(error))
    try:
        if args.command == "candidates":
            rows = database.search_candidates(args.antibiotic_class, args.support_tier, args.padj_max, args.effect_min, args.run_id, args.query)
        elif args.command == "candidate":
            result = database.get_candidate(args.candidate_id, args.run_id)
            rows = [result] if result else []
        elif args.command == "regulator":
            rows = database.search_regulator(args.regulator, args.run_id)
        elif args.command == "interactions":
            rows = database.reference_interactions(args.actor, args.target_kind, args.run_id)
        elif args.command == "promoter-paths":
            rows = database.reference_promoter_paths(args.actor, args.run_id)
        elif args.command == "panel":
            rows = database.panel(args.panel_id, args.run_id)
        elif args.command == "coverage":
            rows = [database.reference_coverage(args.run_id)]
        elif args.command == "network":
            result = database.candidate_network(args.candidate_id, args.run_id)
            rows = [result] if result else []
        elif args.command == "promoter":
            result = database.get_promoter(args.promoter_id, args.run_id)
            rows = [result] if result else []
        elif args.command == "history":
            rows = database.review_history(args.entity_type, args.entity_id)
        elif args.command.startswith("review-"):
            kwargs = {"construct_ready": args.construct_ready, "notes": args.notes,
                      "source_reference": args.source_reference}
            if args.command == "review-candidate":
                kwargs["priority"] = args.priority
                edit_id = database.record_candidate_review(args.run_id, args.entity_id, args.status, args.editor, args.reason, **kwargs)
            else:
                edit_id = database.record_promoter_review(args.run_id, args.entity_id, args.status, args.editor, args.reason, **kwargs)
            rows = [{"edit_id": edit_id}]
        else:
            rows = database.sample_qc(args.run_id)
        _print_rows(rows, args.as_json)
    except ValueError as error:
        parser.error(str(error))
    finally:
        database.close()


if __name__ == "__main__":
    main()
