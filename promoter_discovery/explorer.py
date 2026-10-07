"""Serve the local candidate explorer: python -m promoter_discovery.explorer."""

from __future__ import annotations

import argparse
import json
import sqlite3
from http.server import BaseHTTPRequestHandler, HTTPServer
from pathlib import Path
from urllib.parse import parse_qs, unquote, urlsplit

from .database import DEFAULT_DATABASE
from .db_api import Database


def query_database(database: Database, path: str, parameters: dict[str, list[str]]):
    """Shared JSON routes for the browser and focused integration tests."""
    def value(key):
        return parameters.get(key, [None])[0] or None

    run_id = value("run_id")
    if path == "/api/candidates":
        return database.search_candidates(
            antibiotic_class=value("class"), support_tier=value("support_tier"),
            padj_max=float(value("padj_max")) if value("padj_max") else None,
            effect_min=float(value("effect_min")) if value("effect_min") else None,
            run_id=run_id, query=value("query"),
        )
    if path == "/api/runs":
        return database.runs()
    if path == "/api/coverage":
        return database.reference_coverage(run_id)
    if path == "/api/panel":
        return database.panel(run_id=run_id)
    for prefix, method in (("/api/candidate/", database.get_candidate),
                           ("/api/network/", database.candidate_network),
                           ("/api/promoter/", database.get_promoter),
                           ("/api/regulator/", database.reference_promoter_paths)):
        if path.startswith(prefix):
            result = method(unquote(path[len(prefix):]), run_id)
            if result is None:
                raise LookupError("Record not found")
            return result
    raise LookupError("Route not found")


def create_server(database_path: str | Path, port: int = 8000) -> HTTPServer:
    database_path = Path(database_path).resolve()
    database = Database(database_path)
    database.close()
    page = Path(__file__).with_name("explorer.html").read_bytes()

    class Handler(BaseHTTPRequestHandler):
        def do_GET(self):
            url = urlsplit(self.path)
            if url.path == "/":
                self.respond(200, page, "text/html; charset=utf-8")
                return
            if not url.path.startswith("/api/"):
                self.respond(404, b"Not found", "text/plain; charset=utf-8")
                return
            database = None
            try:
                database = Database(database_path)
                payload = query_database(database, url.path, parse_qs(url.query))
                status = 200
            except LookupError as error:
                status, payload = 404, {"error": str(error)}
            except (ValueError, FileNotFoundError) as error:
                status, payload = 400, {"error": str(error)}
            except sqlite3.Error:
                status, payload = 500, {"error": "Database query failed. Check the server terminal."}
            finally:
                if database is not None:
                    database.close()
            self.respond(status, json.dumps(payload, allow_nan=False).encode(), "application/json; charset=utf-8")

        def respond(self, status, body, content_type):
            self.send_response(status)
            self.send_header("Content-Type", content_type)
            self.send_header("Content-Length", str(len(body)))
            self.send_header("Cache-Control", "no-store")
            self.end_headers()
            self.wfile.write(body)

    return HTTPServer(("127.0.0.1", port), Handler)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--database", type=Path, default=DEFAULT_DATABASE)
    parser.add_argument("--port", type=int, default=8000)
    args = parser.parse_args(argv)
    try:
        server = create_server(args.database, args.port)
    except (OSError, ValueError) as error:
        parser.error(str(error))
    print(f"Promoter explorer: http://127.0.0.1:{server.server_port}", flush=True)
    print(f"Database: {args.database.resolve()}", flush=True)
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        server.server_close()


if __name__ == "__main__":
    main()
