"""Generate two small CellRanger-format .h5 files for demonstration.

The real day-5 and day-8 libraries are unpublished, so this script synthesises
data with the same structure: a 10x HDF5 matrix with planted cell populations
expressing the marker panel in `config/markers_compartment.yaml`, plus mitochondrial genes,
low-quality cells and doublets so that every QC branch is actually exercised.

Running it makes the pipeline reproducible end-to-end from a clean clone with
no downloads and no access to the original data.

    python tests/fixtures/make_synthetic_h5.py
"""

from __future__ import annotations

import argparse
from pathlib import Path

import h5py
import numpy as np
import scipy.sparse as sp
import yaml

RNG_SEED = 0
N_BACKGROUND_GENES = 1200


def write_10x_h5(path: Path, matrix: sp.csc_matrix, barcodes: list[str],
                 genes: list[str], genome: str = "GRCz11") -> None:
    """Write a CellRanger v3 style HDF5 file (genes x cells, CSC)."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(path, "w") as handle:
        grp = handle.create_group("matrix")
        grp.create_dataset("data", data=matrix.data.astype(np.int32))
        grp.create_dataset("indices", data=matrix.indices.astype(np.int64))
        grp.create_dataset("indptr", data=matrix.indptr.astype(np.int64))
        grp.create_dataset("shape", data=np.array(matrix.shape, dtype=np.int32))
        grp.create_dataset("barcodes", data=np.array(barcodes, dtype="S"))

        features = grp.create_group("features")
        n = len(genes)
        features.create_dataset("id", data=np.array(
            [f"ENSDARG{i:011d}" for i in range(n)], dtype="S"))
        features.create_dataset("name", data=np.array(genes, dtype="S"))
        features.create_dataset("feature_type", data=np.array(
            ["Gene Expression"] * n, dtype="S"))
        features.create_dataset("genome", data=np.array([genome] * n, dtype="S"))
        features.create_dataset("_all_tag_keys", data=np.array([b"genome"]))


def simulate(markers: dict[str, list[str]], n_cells: int, seed: int,
             population_weights: dict[str, float]) -> tuple[sp.csc_matrix, list[str], list[str]]:
    """Simulate counts with planted, marker-defined populations."""
    rng = np.random.default_rng(seed)

    marker_genes = [g for genes in markers.values() for g in genes]
    mito_genes = [f"mt-{n}" for n in ("nd1", "nd2", "co1", "co2", "cyb", "atp6")]
    background = [f"gene{i:04d}" for i in range(N_BACKGROUND_GENES)]
    genes = marker_genes + mito_genes + background
    gene_index = {g: i for i, g in enumerate(genes)}

    types = list(markers)
    probs = np.array([population_weights.get(t, 1.0) for t in types], dtype=float)
    probs /= probs.sum()
    assignments = rng.choice(len(types), size=n_cells, p=probs)

    counts = rng.negative_binomial(2, 0.4, size=(n_cells, len(genes))).astype(np.int32)

    for cell, type_idx in enumerate(assignments):
        for gene in markers[types[type_idx]]:
            counts[cell, gene_index[gene]] += rng.negative_binomial(40, 0.35)

    # Depth variation, so MAD-based filtering has something to act on.
    depth = rng.lognormal(0.0, 0.35, size=n_cells)[:, None]
    counts = (counts * depth).astype(np.int32)

    # Stressed cells: high mitochondrial fraction, few genes.
    n_bad = max(1, n_cells // 20)
    bad = rng.choice(n_cells, n_bad, replace=False)
    counts[np.ix_(bad, [gene_index[g] for g in mito_genes])] += 400
    counts[bad[: n_bad // 2], : len(marker_genes)] //= 8

    # Doublets: sum two random cells.
    n_doublet = max(1, n_cells // 25)
    a = rng.choice(n_cells, n_doublet, replace=False)
    b = rng.choice(n_cells, n_doublet, replace=False)
    counts[a] += counts[b]

    barcodes = [f"{''.join(rng.choice(list('ACGT'), 16))}-1" for _ in range(n_cells)]
    # CellRanger stores genes x cells.
    return sp.csc_matrix(counts.T), barcodes, genes


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--n-cells", type=int, default=1500)
    parser.add_argument("--markers", default="config/markers_compartment.yaml")
    parser.add_argument("--outdir", default="data")
    args = parser.parse_args()

    markers = yaml.safe_load(Path(args.markers).read_text())["cell_types"]
    outdir = Path(args.outdir)

    # Different composition per timepoint: later stage is more differentiated
    # and less proliferative, which is the contrast the pipeline should surface.
    compositions = {
        "day5": {"proliferating": 3.0, "ocular_photoreceptor": 0.5,
                 "neural_progenitor_radial_glia": 2.0, "neuron_differentiated": 1.0},
        "day8": {"proliferating": 0.4, "ocular_photoreceptor": 2.0,
                 "neural_progenitor_radial_glia": 0.6, "neuron_differentiated": 2.5},
    }

    for i, (sample, weights) in enumerate(compositions.items()):
        matrix, barcodes, genes = simulate(markers, args.n_cells, RNG_SEED + i, weights)
        path = outdir / f"{sample}_filtered_feature_bc_matrix.h5"
        write_10x_h5(path, matrix, barcodes, genes)
        print(f"wrote {path}  ({matrix.shape[1]} cells x {matrix.shape[0]} genes)")


if __name__ == "__main__":
    main()
