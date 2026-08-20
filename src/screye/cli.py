"""Command-line entry point.

    screye                    # run with config/config.yaml
    screye --dry-run          # resolve and validate everything, run nothing
    screye -c other.yaml -v   # alternative config, debug logging

Run it from anywhere: paths in the config are resolved against the repository
root, located by walking up for `pyproject.toml`.
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

from .config import Config
from .pipeline import run

log = logging.getLogger("screye")


def find_repo_root(start: Path | None = None) -> Path:
    """Walk upwards until `pyproject.toml` is found.

    Anchoring on a marker file rather than the working directory means the CLI,
    the notebook and the test suite all resolve the same paths from wherever
    they happen to be invoked.
    """
    here = (start or Path.cwd()).resolve()
    for candidate in (here, *here.parents):
        if (candidate / "pyproject.toml").exists():
            return candidate
    raise FileNotFoundError(
        f"No pyproject.toml found at or above {here}. "
        "Run screye from inside the repository, or pass --root."
    )


def describe(cfg: Config) -> str:
    """Human-readable summary of exactly what a run would do."""
    lines = ["", "Inputs"]
    for spec in cfg.samples:
        size = f"{spec.h5.stat().st_size / 1e6:,.1f} MB" if spec.h5.exists() else "MISSING"
        mark = " " if spec.h5.exists() else "!"
        lines.append(f"  {mark} {spec.name:<6} {spec.timepoint:<6} {size:>12}  {spec.h5}")

    lines += [
        "",
        "Parameters",
        f"    mitochondrial prefix   {cfg.qc.mito_prefix!r}",
        f"    outlier threshold      {cfg.qc.n_mads} MADs",
        f"    max mitochondrial %    {cfg.qc.max_pct_mito}",
        f"    doublet detection      {'on' if cfg.qc.run_doublet_detection else 'off'}",
        (f"    HVGs / PCs / neighbors {cfg.cluster.n_top_genes} / "
         f"{cfg.cluster.n_pcs} / {cfg.cluster.n_neighbors}"),
        (f"    leiden resolutions     {list(cfg.cluster.leiden_resolutions)} "
         f"(primary {cfg.cluster.primary_resolution})"),
        f"    t-SNE                  {'on' if cfg.cluster.compute_tsne else 'off'}",
        ("    batch correction       "
         f"{'ON via ' + cfg.cluster.batch_key if cfg.cluster.batch_correct else 'off'}"),
        f"    random seed            {cfg.seed}",
        "",
        "Outputs",
        f"    {cfg.outdir}",
        "",
    ]
    return "\n".join(lines)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="screye",
        description="scRNA-seq QC, clustering and annotation across timepoints.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Parameters live in config/config.yaml; nothing is hard-coded in the source.",
    )
    parser.add_argument("-c", "--config", default="config/config.yaml",
                        help="config file, relative to the repo root (default: %(default)s)")
    parser.add_argument("--root", type=Path, default=None,
                        help="repository root; auto-detected from pyproject.toml if omitted")
    parser.add_argument("-n", "--dry-run", action="store_true",
                        help="validate inputs and print the plan without running")
    parser.add_argument("-v", "--verbose", action="store_true", help="debug logging")
    args = parser.parse_args(argv)

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s | %(levelname)-7s | %(message)s",
        datefmt="%H:%M:%S",
    )

    try:
        root = (args.root or find_repo_root()).resolve()
        config_path = Path(args.config)
        if not config_path.is_absolute():
            config_path = root / config_path

        cfg = Config.from_yaml(config_path).resolve(root)
        log.info("Repository root: %s", root)
        log.info("Configuration:   %s", config_path)

        if args.dry_run:
            cfg.validate_inputs()
            print(describe(cfg))
            log.info("Dry run: configuration is valid and all inputs are present.")
            return 0

        print(describe(cfg))
        run(cfg)

    except (FileNotFoundError, ValueError) as exc:
        log.error("%s", exc)
        return 1
    except KeyboardInterrupt:
        log.warning("Interrupted.")
        return 130
    return 0


if __name__ == "__main__":
    sys.exit(main())
