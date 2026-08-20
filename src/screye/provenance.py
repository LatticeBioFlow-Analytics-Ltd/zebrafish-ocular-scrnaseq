"""Run provenance capture.

Every run writes a manifest recording exactly which package versions, input
files and parameters produced the outputs sitting next to it. Without this, a
figure in a manuscript cannot be traced back to the code that made it.
"""

from __future__ import annotations

import hashlib
import json
import logging
import platform
import subprocess
from dataclasses import asdict
from datetime import UTC, datetime
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path

from .config import Config

log = logging.getLogger(__name__)

TRACKED_PACKAGES = (
    "scanpy", "anndata", "numpy", "pandas", "scipy",
    "scikit-learn", "leidenalg", "igraph", "matplotlib", "h5py",
)


def _package_versions() -> dict[str, str]:
    versions = {}
    for package in TRACKED_PACKAGES:
        try:
            versions[package] = version(package)
        except PackageNotFoundError:
            versions[package] = "not installed"
    return versions


def _git_commit() -> str | None:
    """Return the current commit, marked dirty if the tree has uncommitted changes."""
    try:
        commit = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            capture_output=True, text=True, timeout=5, check=True,
        ).stdout.strip()
        dirty = subprocess.run(
            ["git", "status", "--porcelain"],
            capture_output=True, text=True, timeout=5, check=True,
        ).stdout.strip()
        return f"{commit}-dirty" if dirty else commit
    except (subprocess.SubprocessError, FileNotFoundError):
        return None


def _sha256(path: Path, chunk: int = 1 << 20) -> str:
    """Checksum an input file so a changed input cannot pass unnoticed."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(chunk):
            digest.update(block)
    return digest.hexdigest()


def write_manifest(cfg: Config, path: Path, extra: dict | None = None) -> dict:
    """Write `run_manifest.json` alongside the results and return its contents."""
    manifest = {
        "timestamp_utc": datetime.now(UTC).isoformat(timespec="seconds"),
        "git_commit": _git_commit(),
        "python": platform.python_version(),
        "platform": platform.platform(),
        "packages": _package_versions(),
        "seed": cfg.seed,
        "inputs": [
            {
                "name": s.name,
                "timepoint": s.timepoint,
                "path": str(s.h5),
                "bytes": s.h5.stat().st_size,
                "sha256": _sha256(s.h5),
            }
            for s in cfg.samples
        ],
        "parameters": {"qc": asdict(cfg.qc), "cluster": asdict(cfg.cluster)},
    }
    if extra:
        manifest["results"] = extra

    path.write_text(json.dumps(manifest, indent=2, default=str))
    log.info("Provenance written to %s", path)
    return manifest
