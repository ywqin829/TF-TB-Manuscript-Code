#!/usr/bin/env python
"""Extract candidate TF expression from the in vivo GSE326212 BAL atlas.

The full atlas is sparse and large. This script streams the CSR matrix in
row chunks and aggregates expression by sample and coarse cell type. It does
not treat individual cells as biological replicates.
"""

from pathlib import Path
import sys

import anndata as ad
import h5py
import numpy as np
import pandas as pd

PROJECT = Path(__file__).resolve().parents[1]
FULL = Path(
    r"C:\Users\yw_qi\Downloads\tb_host_immune\GSE326212"
    r"\GSE326212_human_tb_bal_adata.h5ad"
)
SUBSET = Path(
    r"C:\Users\yw_qi\Downloads\tb_host_immune\GSE326212"
    r"\outputs\03_subset_annotated.h5ad"
)
OUT = PROJECT / "results" / "gse326212_candidate_tf_pseudobulk.tsv"

CANDIDATES = [
    "MAF", "SNAPC4", "TBX21", "ARNTL", "STAT4",
    "MEF2B", "RUNX3", "ZBED1", "NFATC2", "PAX5",
]
CHUNK = 20_000


def decode(values):
    return [x.decode("utf-8") if isinstance(x, bytes) else str(x) for x in values]


def main():
    if not FULL.exists() or not SUBSET.exists():
        raise FileNotFoundError("GSE326212 input data are not present.")

    subset = ad.read_h5ad(SUBSET, backed="r")
    obs = subset.obs[
        ["sample", "subject_id", "group_norm", "celltype_coarse",
         "predicted_doublet"]
    ].copy()
    obs.index = obs.index.astype(str)
    obs["celltype_coarse"] = obs["celltype_coarse"].astype("string")
    obs["sample"] = obs["sample"].astype("string")
    keep = (
        ~obs["predicted_doublet"].fillna(False).astype(bool)
        & obs["celltype_coarse"].notna()
        & obs["sample"].notna()
    ).to_numpy()

    key = (
        obs["sample"].fillna("") + "|" +
        obs["celltype_coarse"].fillna("")
    )
    key_cat = pd.Categorical(key)
    key_codes = key_cat.codes.astype(np.int64)
    n_keys = len(key_cat.categories)

    with h5py.File(FULL, "r") as f:
        genes = decode(f["var"]["_index"][:])
        gene_to_index = {g: i for i, g in enumerate(genes)}
        present = [g for g in CANDIDATES if g in gene_to_index]
        missing = sorted(set(CANDIDATES) - set(present))
        target_indices = np.asarray(
            [gene_to_index[g] for g in present], dtype=np.int64
        )

        # A lookup table turns each sparse column index into target-column
        # position (or -1 when the gene is not requested).
        lookup = np.full(len(genes), -1, dtype=np.int32)
        target_position = {int(v): i for i, v in enumerate(target_indices)}
        for v, i in target_position.items():
            lookup[v] = i

        indptr = f["X"]["indptr"][:].astype(np.int64)
        n_cells = len(indptr) - 1
        if n_cells != len(obs):
            raise ValueError(
                f"Full/subset cell mismatch: {n_cells} vs {len(obs)}"
            )

        sums = np.zeros((n_keys, len(present)), dtype=np.float64)
        nz = np.zeros((n_keys, len(present)), dtype=np.float64)
        counts = np.zeros(n_keys, dtype=np.int64)

        for start in range(0, n_cells, CHUNK):
            end = min(start + CHUNK, n_cells)
            lo = int(indptr[start])
            hi = int(indptr[end])
            ix = f["X"]["indices"][lo:hi]
            data = f["X"]["data"][lo:hi]
            local_ptr = indptr[start:end + 1] - lo
            row_ids = np.repeat(
                np.arange(end - start, dtype=np.int64),
                np.diff(local_ptr),
            )
            target_col = lookup[ix.astype(np.int64)]
            hit = target_col >= 0
            dense = np.zeros(
                (end - start, len(present)), dtype=np.float64
            )
            dense[row_ids[hit], target_col[hit]] = data[hit]

            valid = keep[start:end]
            codes = key_codes[start:end][valid]
            if codes.size:
                np.add.at(sums, codes, dense[valid])
                np.add.at(nz, codes, (dense[valid] > 0).astype(float))
                np.add.at(counts, codes, 1)

            if start % (CHUNK * 5) == 0:
                print(f"processed {end:,}/{n_cells:,} cells", flush=True)

    key_meta = key_cat.categories.to_series(index=range(n_keys))
    key_meta = key_meta.str.split("|", n=1, expand=True, regex=False)
    key_meta.columns = ["sample", "celltype_coarse"]
    key_meta["sample"] = key_meta["sample"].astype(str)

    sample_lookup = obs.drop_duplicates("sample").set_index("sample")[
        ["subject_id", "group_norm"]
    ]
    rows = []
    for i in range(n_keys):
        if counts[i] == 0:
            continue
        sample = key_meta.loc[i, "sample"]
        subject = sample_lookup.loc[sample, "subject_id"] if sample in sample_lookup.index else pd.NA
        group = sample_lookup.loc[sample, "group_norm"] if sample in sample_lookup.index else pd.NA
        for j, tf in enumerate(present):
            rows.append({
                "sample": sample,
                "subject_id": subject,
                "group_norm": group,
                "celltype_coarse": key_meta.loc[i, "celltype_coarse"],
                "n_cells": int(counts[i]),
                "TF": tf,
                "mean_expression": sums[i, j] / counts[i],
                "percent_expressing": 100.0 * nz[i, j] / counts[i],
            })

    result = pd.DataFrame(rows).sort_values(
        ["celltype_coarse", "TF", "group_norm", "sample"]
    )
    OUT.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(OUT, sep="\t", index=False)

    print(f"Wrote {OUT}")
    print(f"cells after QC: {keep.sum():,}; retained TF genes: {present}")
    if missing:
        print(f"Missing TF genes: {missing}", file=sys.stderr)


if __name__ == "__main__":
    main()
