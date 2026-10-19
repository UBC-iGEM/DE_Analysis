"""Run the benchmark-to-database workflow in order from the repository root."""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def run_workflow(*, download: bool = False, reuse_fit: bool = False) -> None:
    """Stop at the first failed stage; assessment builds the database exactly once."""
    stages = [("setup_data", ["--download"] if download else [])]
    if not reuse_fit:
        stages.append(("run_benchmark", []))
    stages.extend((stage, []) for stage in (
        "build_network", "score_candidates", "promoter_selection", "assess_candidates", "visualize_network",
    ))
    for stage, arguments in stages:
        print(f"[workflow] {stage}", flush=True)
        subprocess.run([sys.executable, "-m", f"promoter_discovery.{stage}", *arguments],
                       cwd=ROOT, check=True)
    print("[workflow] Complete. Run python -m promoter_discovery.explorer", flush=True)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--download", action="store_true", help="fetch locked reference assets if needed")
    parser.add_argument("--reuse-fit", action="store_true", help="reuse the existing joint fit; downstream provenance is still checked")
    args = parser.parse_args(argv)
    try:
        run_workflow(download=args.download, reuse_fit=args.reuse_fit)
    except subprocess.CalledProcessError as error:
        parser.exit(error.returncode, f"Workflow stopped: {error.cmd[2]} failed. See the stage output above.\n")


if __name__ == "__main__":
    main()
