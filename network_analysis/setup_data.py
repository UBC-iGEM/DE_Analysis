"""Validate/download versioned network-analysis assets.

The upstream product carrying the requested ``iModulon >=2.5.0`` identifier is
intentionally supplied by the owner in the manifest; this script does not
guess a provider or silently substitute a Python package version.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import urllib.request
from pathlib import Path
from typing import Any

try:
    from .dataset_registry import normalize_version, version_at_least
except ImportError:  # pragma: no cover
    from dataset_registry import normalize_version, version_at_least  # type: ignore


MINIMUM_REGULONDB = "14.5.0"
MINIMUM_IMODULON = "2.5.0"
REQUIRED_ASSET_FIELDS = {"name", "kind", "url", "path", "version", "sha256", "license", "citation", "retrieved_at", "official_revision"}


def sha256_file(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _asset_version(asset: dict[str, Any]) -> str | None:
    return str(asset.get("version") or asset.get("release") or "").strip() or None


def load_manifest(path: str | Path) -> dict[str, Any]:
    with Path(path).open(encoding="utf-8") as handle:
        manifest = json.load(handle)
    if not isinstance(manifest.get("assets"), list) or not manifest["assets"]:
        raise ValueError("Asset manifest must contain a non-empty assets list")
    names: set[str] = set()
    for asset in manifest["assets"]:
        missing = {field for field in REQUIRED_ASSET_FIELDS if not asset.get(field)}
        if missing:
            raise ValueError(f"Asset {asset.get('name', '<unnamed>')} missing manifest fields: {sorted(missing)}")
        if asset["name"] in names:
            raise ValueError(f"Duplicate asset name: {asset['name']}")
        names.add(asset["name"])
        if not re.fullmatch(r"[0-9a-fA-F]{64}", str(asset["sha256"])):
            raise ValueError(f"Asset {asset['name']} sha256 must be a 64-character hex digest")
        if normalize_version(_asset_version(asset)) is None:
            raise ValueError(f"Asset {asset['name']} needs a numeric release/version")
    return manifest


def _regulondb_release(path: Path) -> str | None:
    with path.open(encoding="utf-8", errors="replace") as handle:
        for _ in range(80):
            line = handle.readline()
            if not line:
                break
            match = re.search(r"RegulonDB\s+Release\s*:\s*v?([0-9]+(?:\.[0-9]+)*)", line, re.IGNORECASE)
            if match:
                return match.group(1)
    return None


def validate_manifest(manifest: dict[str, Any], root: str | Path | None = None, check_files: bool = True) -> list[str]:
    """Validate provenance, minimum versions, checksums, and iModulon pairing."""

    root_path = Path(root or ".")
    errors: list[str] = []
    regulondb_assets = [asset for asset in manifest["assets"] if str(asset["kind"]).lower().startswith("regulondb")]
    imodulon_assets = [asset for asset in manifest["assets"] if str(asset["kind"]).lower() in {"imodulon", "imodulondb", "expression", "imodulon_expression"}]
    for asset in regulondb_assets:
        version = _asset_version(asset)
        if not version_at_least(version, MINIMUM_REGULONDB):
            errors.append(f"{asset['name']} release {version!r} is below RegulonDB {MINIMUM_REGULONDB}")
        if check_files:
            path = root_path / asset["path"]
            if not path.exists():
                errors.append(f"Missing local asset {asset['name']}: {path}")
            elif sha256_file(path).lower() != str(asset["sha256"]).lower():
                errors.append(f"SHA256 mismatch for {asset['name']}: {path}")
            elif str(asset["kind"]).lower() == "regulondb":
                observed = _regulondb_release(path)
                if observed and not version_at_least(observed, MINIMUM_REGULONDB):
                    errors.append(f"{asset['name']} file declares RegulonDB {observed}, below {MINIMUM_REGULONDB}")
    for asset in imodulon_assets:
        version = _asset_version(asset)
        if not version_at_least(version, MINIMUM_IMODULON):
            errors.append(f"{asset['name']} release/version {version!r} is below iModulon {MINIMUM_IMODULON}")
        if check_files:
            path = root_path / asset["path"]
            if not path.exists():
                errors.append(f"Missing local asset {asset['name']}: {path}")
            elif sha256_file(path).lower() != str(asset["sha256"]).lower():
                errors.append(f"SHA256 mismatch for {asset['name']}: {path}")
    model = next((asset for asset in imodulon_assets if str(asset["kind"]).lower() in {"imodulon", "imodulondb"}), None)
    expression = next((asset for asset in imodulon_assets if str(asset["kind"]).lower() in {"expression", "imodulon_expression"}), None)
    if model is None or expression is None:
        errors.append("Manifest must identify both an iModulon model and its companion expression asset")
    elif (_asset_version(model), model.get("asset_id")) != (_asset_version(expression), expression.get("asset_id")):
        errors.append("iModulon model and expression assets must share the same release/version and asset_id")
    return errors


def download_assets(manifest: dict[str, Any], root: str | Path | None = None) -> None:
    root_path = Path(root or ".")
    for asset in manifest["assets"]:
        destination = root_path / asset["path"]
        destination.parent.mkdir(parents=True, exist_ok=True)
        if destination.exists():
            observed = sha256_file(destination)
            if observed.lower() == str(asset["sha256"]).lower():
                print(f"[setup] verified existing {asset['name']}: {destination}")
                continue
            raise ValueError(f"Existing file checksum does not match manifest for {asset['name']}: {destination}")
        if not str(asset["url"]).startswith(("file://", "https://", "http://")):
            raise ValueError(f"Asset {asset['name']} URL must be an explicit http(s) or file URL")
        print(f"[setup] downloading {asset['name']} -> {destination}")
        urllib.request.urlretrieve(asset["url"], destination)
        observed = sha256_file(destination)
        if observed.lower() != str(asset["sha256"]).lower():
            raise ValueError(f"Downloaded checksum does not match manifest for {asset['name']}")


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", default=root / "config" / "network_assets.json")
    parser.add_argument("--root", default=root)
    parser.add_argument("--download", action="store_true", help="Download missing assets, then validate them")
    parser.add_argument("--no-file-check", action="store_true", help="Validate manifest metadata without requiring files")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv)
    manifest = load_manifest(args.manifest)
    if args.download:
        download_assets(manifest, args.root)
    errors = validate_manifest(manifest, args.root, check_files=not args.no_file_check)
    if errors:
        raise SystemExit("[setup] validation failed:\n- " + "\n- ".join(errors))
    print(f"[setup] validated {len(manifest['assets'])} assets")


if __name__ == "__main__":
    main()
