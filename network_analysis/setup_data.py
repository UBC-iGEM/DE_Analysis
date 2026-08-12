"""Acquire and validate versioned regulatory-network analysis assets.

The manifest describes *where* an asset comes from.  A generated lock file
records the exact bytes retrieved, including release metadata and SHA-256,
so subsequent runs can verify the cached inputs without requiring hand-copied
hashes in the source manifest.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import ssl
import tempfile
import time
import urllib.error
import urllib.request
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

try:
    from .dataset_registry import normalize_version, version_at_least
except ImportError:  # pragma: no cover
    from dataset_registry import normalize_version, version_at_least  # type: ignore


MINIMUM_REGULONDB = "14.5.0"
MINIMUM_PRECISE1K = "1.0"
DEFAULT_REGULONDB_ENDPOINT = "https://regulondb-prerelease.ccg.unam.mx/graphql"
GRAPHQL_TIMEOUT_SECONDS = 45
GRAPHQL_RETRIES = 2
LOCK_SCHEMA_VERSION = 1
GRAPHQL_QUERY = """query($fileName: String!) {
  getDataOfFile(fileName: $fileName) {
    _id
    fileName
    content
    columnsDetails
  }
}"""


def _ssl_context() -> ssl.SSLContext:
    """Build a verified context, honoring an explicitly configured CA bundle."""

    cafile = os.environ.get("SSL_CERT_FILE")
    if cafile:
        return ssl.create_default_context(cafile=cafile)
    try:
        import certifi
    except ImportError:  # pragma: no cover - certifi is normally transitive
        return ssl.create_default_context()
    return ssl.create_default_context(cafile=certifi.where())


def sha256_file(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def _asset_release(asset: dict[str, Any]) -> str | None:
    value = asset.get("dataset_release") or asset.get("release")
    if isinstance(value, str):
        return value.strip() or None
    return None


def _asset_revision(asset: dict[str, Any]) -> str | None:
    for key in ("source_revision", "official_revision"):
        value = asset.get(key)
        if value:
            return str(value).strip()
    return None


def _asset_provider(asset: dict[str, Any]) -> str:
    provider = asset.get("provider")
    if provider:
        return str(provider).strip().lower()
    return "http"


def _asset_expected_hash(asset: dict[str, Any]) -> str | None:
    value = asset.get("expected_sha256") or asset.get("sha256")
    if not value or str(value).startswith("REPLACE_WITH"):
        return None
    return str(value).lower()


def load_manifest(path: str | Path) -> dict[str, Any]:
    with Path(path).open(encoding="utf-8") as handle:
        manifest = json.load(handle)
    if not isinstance(manifest, dict) or not isinstance(manifest.get("assets"), list) or not manifest["assets"]:
        raise ValueError("Asset manifest must contain a non-empty assets list")
    names: set[str] = set()
    for asset in manifest["assets"]:
        if not isinstance(asset, dict):
            raise ValueError("Each asset manifest entry must be an object")
        missing = {field for field in ("name", "kind", "path") if not asset.get(field)}
        if missing:
            raise ValueError(f"Asset {asset.get('name', '<unnamed>')} missing manifest fields: {sorted(missing)}")
        if asset["name"] in names:
            raise ValueError(f"Duplicate asset name: {asset['name']}")
        names.add(str(asset["name"]))
        provider = _asset_provider(asset)
        if provider in {"regulondb_graphql", "graphql"}:
            if not asset.get("remote_name"):
                raise ValueError(f"GraphQL asset {asset['name']} needs remote_name")
        elif provider in {"http", "https", "file"}:
            if not str(asset.get("url", "")).startswith(("file://", "https://", "http://")):
                raise ValueError(f"HTTP asset {asset['name']} needs an explicit http(s) or file URL")
        else:
            raise ValueError(f"Unsupported asset provider for {asset['name']}: {provider}")
        expected = _asset_expected_hash(asset)
        if expected is not None and not re.fullmatch(r"[0-9a-f]{64}", expected):
            raise ValueError(f"Asset {asset['name']} expected_sha256 must be a 64-character hex digest")
        for key in ("dataset_release", "source_revision"):
            if asset.get(key) is not None and not str(asset[key]).strip():
                raise ValueError(f"Asset {asset['name']} has an empty {key}")
    return manifest


def load_lock(path: str | Path) -> dict[str, Any] | None:
    lock_path = Path(path)
    if not lock_path.exists():
        return None
    with lock_path.open(encoding="utf-8") as handle:
        lock = json.load(handle)
    if not isinstance(lock, dict) or lock.get("schema_version") != LOCK_SCHEMA_VERSION:
        raise ValueError(f"Unsupported asset lock schema: {lock_path}")
    if not isinstance(lock.get("assets"), dict):
        raise ValueError(f"Asset lock must contain an assets object: {lock_path}")
    return lock


def _atomic_write_bytes(destination: Path, content: bytes) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary_name = tempfile.mkstemp(prefix=f".{destination.name}.", dir=destination.parent)
    temporary = Path(temporary_name)
    try:
        with os.fdopen(fd, "wb") as handle:
            handle.write(content)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, destination)
    finally:
        temporary.unlink(missing_ok=True)


def _atomic_write_json(destination: Path, payload: dict[str, Any]) -> None:
    encoded = (json.dumps(payload, indent=2, ensure_ascii=False) + "\n").encode("utf-8")
    _atomic_write_bytes(destination, encoded)


def _graphql_request(endpoint: str, query: str, variables: dict[str, Any] | None = None) -> dict[str, Any]:
    payload = json.dumps({"query": query, "variables": variables or {}}).encode("utf-8")
    request = urllib.request.Request(endpoint, data=payload, headers={"Content-Type": "application/json", "Accept": "application/json"}, method="POST")
    last_error: Exception | None = None
    for attempt in range(GRAPHQL_RETRIES + 1):
        try:
            with urllib.request.urlopen(request, timeout=GRAPHQL_TIMEOUT_SECONDS, context=_ssl_context()) as response:
                result = json.loads(response.read().decode("utf-8"))
            if result.get("errors"):
                raise ValueError(f"RegulonDB GraphQL errors: {result['errors']}")
            data = result.get("data")
            if not isinstance(data, dict):
                raise ValueError("RegulonDB GraphQL response has no data object")
            return data
        except (OSError, urllib.error.URLError, TimeoutError, json.JSONDecodeError, ValueError) as error:
            last_error = error
            if attempt < GRAPHQL_RETRIES:
                time.sleep(0.5 * (attempt + 1))
    raise RuntimeError(f"RegulonDB GraphQL request failed after retries: {last_error}") from last_error


def _select_regulondb_release(info: Any, minimum: str) -> dict[str, Any]:
    if not isinstance(info, list):
        raise ValueError("RegulonDB getDatabaseInfo did not return a list")
    candidates = [row for row in info if isinstance(row, dict) and normalize_version(str(row.get("regulonDBVersion", "")))]
    if not candidates:
        raise ValueError("RegulonDB getDatabaseInfo returned no numeric release")
    selected = max(candidates, key=lambda row: normalize_version(str(row["regulonDBVersion"])) or ())
    release = str(selected["regulonDBVersion"])
    if not version_at_least(release, minimum):
        raise ValueError(f"RegulonDB release {release} is below required {minimum}")
    return selected


def _regulondb_context(manifest: dict[str, Any], assets: Iterable[dict[str, Any]]) -> dict[str, Any]:
    endpoints = {str(asset.get("endpoint") or manifest.get("regulondb", {}).get("endpoint") or DEFAULT_REGULONDB_ENDPOINT) for asset in assets}
    if len(endpoints) != 1:
        raise ValueError("All RegulonDB GraphQL assets must use one endpoint per acquisition")
    endpoint = endpoints.pop()
    data = _graphql_request(endpoint, "{ getDatabaseInfo { regulonDBVersion ecocycVersion lcVersion releaseDate note } listAllFileNames }")
    release_info = _select_regulondb_release(data.get("getDatabaseInfo"), str(manifest.get("minimum_versions", {}).get("regulondb", MINIMUM_REGULONDB)))
    available = set(data.get("listAllFileNames") or [])
    required = {str(asset["remote_name"]) for asset in assets}
    missing = required - available
    if missing:
        raise ValueError(f"RegulonDB GraphQL products are unavailable: {sorted(missing)}")
    return {"endpoint": endpoint, "release": release_info, "available_products": sorted(available)}


def _fetch_http(url: str) -> bytes:
    with urllib.request.urlopen(url, timeout=GRAPHQL_TIMEOUT_SECONDS, context=_ssl_context()) as response:
        return response.read()


def _record_asset(asset: dict[str, Any], destination: Path, retrieved_at: str, **extra: Any) -> dict[str, Any]:
    record: dict[str, Any] = {
        "path": str(asset["path"]),
        "provider": _asset_provider(asset),
        "source_url": asset.get("url"),
        "source_revision": _asset_revision(asset),
        "dataset_release": _asset_release(asset),
        "sha256": sha256_file(destination),
        "size": destination.stat().st_size,
        "retrieved_at": retrieved_at,
    }
    record.update({key: value for key, value in extra.items() if value is not None})
    expected = _asset_expected_hash(asset)
    if expected and record["sha256"].lower() != expected:
        raise ValueError(f"SHA256 mismatch for {asset['name']}: {destination}")
    return record


def _mapping_tokens(value: str) -> list[str]:
    return [token.strip() for token in re.split(r"[|;,]", value or "") if token.strip()]


def normalize_gene_product_mapping(raw_path: str | Path, destination: str | Path, source_release: str) -> dict[str, Any]:
    """Derive the stable mapping contract consumed by build_network.py."""

    raw = Path(raw_path)
    rows: list[dict[str, str]] = []
    with raw.open(encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header: list[str] | None = None
        for values in reader:
            if not values or all(not value.strip() for value in values):
                continue
            cleaned = [re.sub(r"^\d+\)", "", value.strip().lstrip("\ufeff")).lower() for value in values]
            if header is None:
                header = cleaned
                continue
            row = dict(zip(header, values, strict=False))
            gene = row.get("genename", "").strip()
            gene_id = row.get("geneid", "").strip()
            synonyms = _mapping_tokens(row.get("genesynonyms", ""))
            external = row.get("otherdbsgeneids", "")
            b_numbers = re.findall(r"\bb\d{4}\b", external, re.IGNORECASE)
            b_numbers = list(dict.fromkeys(value.lower() for value in b_numbers))
            status = "mapped" if len(b_numbers) == 1 else "ambiguous" if len(b_numbers) > 1 else "missing_locus_tag"
            aliases = list(dict.fromkeys([gene.lower(), *[value.lower() for value in synonyms], *b_numbers, gene_id.lower()]))
            rows.append({
                "canonical_gene": gene.lower(),
                "canonical_locus_tag": b_numbers[0] if len(b_numbers) == 1 else "",
                "aliases": "|".join(value for value in aliases if value),
                "regulondb_gene_id": gene_id,
                "identity_mapping_status": status,
                "source_release": source_release,
            })
    if not rows:
        raise ValueError(f"No gene rows parsed from RegulonDB mapping file: {raw}")
    output = Path(destination)
    output.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary_name = tempfile.mkstemp(prefix=f".{output.name}.", dir=output.parent, text=True)
    temporary = Path(temporary_name)
    try:
        with os.fdopen(fd, "w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)
        os.replace(temporary, output)
    finally:
        temporary.unlink(missing_ok=True)
    return {"path": str(output), "sha256": sha256_file(output), "size": output.stat().st_size, "rows": len(rows), "source_release": source_release}


def download_assets(manifest: dict[str, Any], root: str | Path | None = None, lock_path: str | Path | None = None, refresh_lock: bool = False) -> dict[str, Any]:
    root_path = Path(root or ".")
    lock_file = Path(lock_path) if lock_path else root_path / "config" / "network_assets.lock.json"
    prior = load_lock(lock_file) or {"schema_version": LOCK_SCHEMA_VERSION, "assets": {}}
    assets = [asset for asset in manifest["assets"]]
    asset_names = {str(asset["name"]) for asset in assets}
    lock_assets: dict[str, Any] = {name: record for name, record in prior.get("assets", {}).items() if name in asset_names}
    regulondb_assets = [asset for asset in assets if _asset_provider(asset) in {"regulondb_graphql", "graphql"}]
    context = _regulondb_context(manifest, regulondb_assets) if regulondb_assets else None
    release = str(context["release"]["regulonDBVersion"]) if context else None
    for asset in assets:
        destination = root_path / asset["path"]
        prior_record = lock_assets.get(asset["name"], {})
        extra = {key: prior_record[key] for key in ("remote_name", "remote_id", "columns_details", "regulondb_release") if key in prior_record}
        if destination.exists() and prior_record and not refresh_lock:
            observed = sha256_file(destination)
            if observed.lower() != str(prior_record.get("sha256", "")).lower():
                raise ValueError(f"Existing file checksum does not match lock for {asset['name']}: {destination}")
            print(f"[setup] verified cached {asset['name']}: {destination}")
        elif destination.exists() and not refresh_lock and not prior_record:
            print(f"[setup] adopting existing untracked asset {asset['name']}: {destination}")
        else:
            provider = _asset_provider(asset)
            if provider in {"regulondb_graphql", "graphql"}:
                response = _graphql_request(str(context["endpoint"]), GRAPHQL_QUERY, {"fileName": asset["remote_name"]})
                remote = response.get("getDataOfFile")
                if not isinstance(remote, dict) or remote.get("fileName") != asset["remote_name"] or not isinstance(remote.get("content"), str):
                    raise ValueError(f"RegulonDB GraphQL returned malformed product for {asset['name']}")
                _atomic_write_bytes(destination, remote["content"].encode("utf-8"))
                extra = {
                    "source_url": context["endpoint"],
                    "remote_name": remote.get("fileName"),
                    "remote_id": remote.get("_id"),
                    "columns_details": remote.get("columnsDetails"),
                    "regulondb_release": release,
                }
            else:
                url = str(asset["url"])
                print(f"[setup] downloading {asset['name']} -> {destination}")
                _atomic_write_bytes(destination, _fetch_http(url))
                extra = {}
        record = _record_asset(asset, destination, prior_record.get("retrieved_at") or _now(), **extra)
        lock_assets[asset["name"]] = record
    mapping_asset = next((asset for asset in assets if asset.get("derived_path")), None)
    if mapping_asset:
        raw_path = root_path / mapping_asset["path"]
        derived = normalize_gene_product_mapping(raw_path, root_path / mapping_asset["derived_path"], release or str(_asset_release(mapping_asset) or "unknown"))
        lock_assets[f"{mapping_asset['name']}:derived"] = {"source_asset": mapping_asset["name"], **derived, "retrieved_at": _now()}
    lock = {"schema_version": LOCK_SCHEMA_VERSION, "generated_at": _now(), "assets": lock_assets}
    _atomic_write_json(lock_file, lock)
    return lock


def validate_manifest(manifest: dict[str, Any], root: str | Path | None = None, check_files: bool = True, lock: dict[str, Any] | None = None, require_lock: bool = False) -> list[str]:
    """Validate source metadata, local assets, lock hashes, and release pairing."""

    root_path = Path(root or ".")
    errors: list[str] = []
    min_regulondb = str(manifest.get("minimum_versions", {}).get("regulondb", MINIMUM_REGULONDB))
    min_precise = str(manifest.get("minimum_versions", {}).get("precise1k", MINIMUM_PRECISE1K))
    assets = manifest["assets"]
    regulondb_assets = [asset for asset in assets if str(asset["kind"]).lower().startswith("regulondb")]
    precise_assets = [asset for asset in assets if str(asset["kind"]).lower() in {"imodulon", "imodulondb", "expression", "imodulon_expression"}]
    lock_assets = (lock or {}).get("assets", {})
    release_values = {str((lock_assets.get(asset["name"], {}) or {}).get("regulondb_release") or _asset_release(asset) or "") for asset in regulondb_assets}
    release_values.discard("")
    if release_values and not all(version_at_least(value, min_regulondb) for value in release_values):
        errors.append(f"RegulonDB release must be >= {min_regulondb}: {sorted(release_values)}")
    for asset in assets:
        destination = root_path / asset["path"]
        if check_files:
            if not destination.exists():
                errors.append(f"Missing local asset {asset['name']}: {destination}")
                continue
            observed = sha256_file(destination)
            expected = _asset_expected_hash(asset)
            if expected and observed.lower() != expected:
                errors.append(f"SHA256 mismatch for {asset['name']}: {destination}")
            record = lock_assets.get(asset["name"])
            if require_lock and not record:
                errors.append(f"Missing lock entry for {asset['name']}")
            elif record:
                if observed.lower() != str(record.get("sha256", "")).lower():
                    errors.append(f"SHA256 mismatch against lock for {asset['name']}: {destination}")
                if destination.stat().st_size != int(record.get("size", destination.stat().st_size)):
                    errors.append(f"Size mismatch against lock for {asset['name']}: {destination}")
    model = next((asset for asset in precise_assets if str(asset["kind"]).lower() in {"imodulon", "imodulondb"}), None)
    expression = next((asset for asset in precise_assets if str(asset["kind"]).lower() in {"expression", "imodulon_expression"}), None)
    if model is None or expression is None:
        errors.append("Manifest must identify both an iModulon model and a companion expression fallback")
    else:
        model_release = _asset_release(model)
        expression_release = _asset_release(expression)
        if not version_at_least(model_release, min_precise) or not version_at_least(expression_release, min_precise):
            errors.append(f"PRECISE-1K dataset release must be >= {min_precise}")
        if (model_release, _asset_revision(model), model.get("asset_id")) != (expression_release, _asset_revision(expression), expression.get("asset_id")):
            errors.append("iModulon model and expression assets must share dataset_release, source_revision, and asset_id")
    return errors


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", default=root / "config" / "network_assets.json")
    parser.add_argument("--root", default=root)
    parser.add_argument("--lock", default=root / "config" / "network_assets.lock.json")
    parser.add_argument("--download", action="store_true", help="Download missing assets and generate/update the lock")
    parser.add_argument("--refresh-lock", action="store_true", help="Replace cached assets and accept new observed hashes")
    parser.add_argument("--no-file-check", action="store_true", help="Validate source metadata without requiring local files")
    parser.add_argument("--no-lock", action="store_true", help="Do not require a generated lock entry")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv)
    manifest = load_manifest(args.manifest)
    if args.download:
        lock = download_assets(manifest, args.root, args.lock, args.refresh_lock)
    else:
        lock = load_lock(args.lock)
    errors = validate_manifest(manifest, args.root, check_files=not args.no_file_check, lock=lock, require_lock=not args.no_lock and not args.no_file_check)
    if errors:
        raise SystemExit("[setup] validation failed:\n- " + "\n- ".join(errors))
    print(f"[setup] validated {len(manifest['assets'])} assets")


if __name__ == "__main__":
    main()
