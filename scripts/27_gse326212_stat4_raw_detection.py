#!/usr/bin/env python
"""Extract STAT4 raw-count detection metrics from GSE326212."""

from pathlib import Path

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
OUT = PROJECT / "results" / "gse326212_stat4_raw_detection.tsv"
CHUNK = 20_000


def decode(values):
    return [x.decode("utf-8") if isinstance(x, bytes) else str(x) for x in values]


def main():
    subset = ad.read_h5ad(SUBSET, backed="r")
    obs = subset.obs[
        ["sample", "subject_id", "group_norm", "celltype_coarse", "predicted_doublet"]
    ].copy()
    keep = (
        ~obs["predicted_doublet"].fillna(False).astype(bool)
        & obs["celltype_coarse"].notna()
        & obs["sample"].notna()
    ).to_numpy()

    key = obs["sample"].astype(str) + "|" + obs["celltype_coarse"].astype(str)
    key_cat = pd.Categorical(key)
    codes = key_cat.codes.astype(np.int64)
    n_keys = len(key_cat.categories)

    with h5py.File(FULL, "r") as f:
        genes = decode(f["var"]["_index"][:])
        stat4_idx = genes.index("STAT4")
        ptr = f["layers"]["counts"]["indptr"][:].astype(np.int64)
        n_cells = len(ptr) - 1

        total = np.zeros(n_keys, dtype=np.float64)
        nz1 = np.zeros(n_keys, dtype=np.float64)
        nz2 = np.zeros(n_keys, dtype=np.float64)
        counts = np.zeros(n_keys, dtype=np.int64)

        for start in range(0, n_cells, CHUNK):
            end = min(start + CHUNK, n_cells)
            lo = int(ptr[start])
            hi = int(ptr[end])
            ix = f["layers"]["counts"]["indices"][lo:hi]
            data = f["layers"]["counts"]["data"][lo:hi]
            local_ptr = ptr[start:end + 1] - lo
            row_ids = np.repeat(
                np.arange(end - start, dtype=np.int64),
                np.diff(local_ptr),
            )
            hit = ix == stat4_idx
            values = np.zeros(end - start, dtype=np.float64)
            values[row_ids[hit]] = data[hit]

            valid = keep[start:end]
            group = codes[start:end][valid]
            v = values[valid]
            if group.size:
                np.add.at(total, group, v)
                np.add.at(nz1, group, v >= 1)
                np.add.at(nz2, group, v >= 2)
                np.add.at(counts, group, 1)

    sample_lookup = obs.drop_duplicates("sample").set_index("sample")[
        ["subject_id", "group_norm"]
    ]
    rows = []
    for i, category in enumerate(key_cat.categories):
        if counts[i] == 0:
            continue
        sample, celltype = category.split("|", 1)
        rows.append({
            "sample": sample,
            "subject_id": sample_lookup.loc[sample, "subject_id"],
            "group_norm": sample_lookup.loc[sample, "group_norm"],
            "celltype_coarse": celltype,
            "n_cells": int(counts[i]),
            "mean_stat4_count": total[i] / counts[i],
            "percent_count_ge1": 100 * nz1[i] / counts[i],
            "percent_count_ge2": 100 * nz2[i] / counts[i],
        })

    pd.DataFrame(rows).sort_values(
        ["celltype_coarse", "group_norm", "sample"]
    ).to_csv(OUT, sep="\t", index=False)
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
