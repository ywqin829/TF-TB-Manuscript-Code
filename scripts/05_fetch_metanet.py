#!/usr/bin/env python3
"""Extract selected METANet networks without downloading the 15 GB compendium."""

from __future__ import annotations

import csv
import gzip
import io
import json
import math
import struct
import time
import zlib
from pathlib import Path

import pandas as pd
import requests


PROJECT = Path(__file__).resolve().parents[1]
RESULTS = PROJECT / "results"
URL = (
    "https://zenodo.org/api/records/20614479/files/"
    "METANet_compendium.zip/content"
)
RECORD_API = "https://zenodo.org/api/records/20614479"
CHUNK = 2_000_000

TARGET_GENES = {
    "TBX21",
    "PAX5",
    "NFATC2",
    "IFNG",
    "CD40LG",
    "EBF1",
    "CD19",
    "CD79A",
    "BLK",
    "FCER2",
    "POU2AF1",
    "IKZF3",
}
CORE_TFS = {"TBX21", "PAX5", "NFATC2"}


def request_range(url: str, start: int, end: int, retries: int = 3) -> bytes:
    headers = {"Range": f"bytes={start}-{end}"}
    for attempt in range(retries):
        try:
            response = requests.get(url, headers=headers, timeout=120)
            response.raise_for_status()
            if response.status_code != 206:
                raise RuntimeError(f"Range request returned {response.status_code}")
            return response.content
        except Exception:
            if attempt + 1 == retries:
                raise
            time.sleep(2 ** attempt)
    raise RuntimeError("Unreachable")


def parse_zip64_extra(extra: bytes, values: dict[str, int]) -> dict[str, int]:
    pos = 0
    while pos + 4 <= len(extra):
        tag, length = struct.unpack_from("<2H", extra, pos)
        data = extra[pos + 4 : pos + 4 + length]
        pos += 4 + length
        if tag != 1:
            continue
        cursor = 0
        for key in ("uncompressed_size", "compressed_size", "local_header_offset"):
            if values[key] == 0xFFFFFFFF:
                values[key] = struct.unpack_from("<Q", data, cursor)[0]
                cursor += 8
    return values


def central_directory(url: str, total_size: int) -> list[dict]:
    tail_start = max(0, total_size - CHUNK)
    tail = request_range(url, tail_start, total_size - 1)
    eocd = tail.rfind(b"PK\x05\x06")
    if eocd < 0:
        raise RuntimeError("ZIP EOCD not found")
    entries = struct.unpack_from("<4s4H2LH", tail, eocd)
    central_size = entries[5]
    central_offset = entries[6]

    if central_offset == 0xFFFFFFFF or central_size == 0xFFFFFFFF:
        locator = tail.rfind(b"PK\x06\x07", 0, eocd)
        if locator < 0:
            raise RuntimeError("ZIP64 locator not found")
        zip64_offset = struct.unpack_from("<Q", tail, locator + 8)[0]
        zip64_data = request_range(url, zip64_offset, zip64_offset + 200)
        marker = zip64_data.find(b"PK\x06\x06")
        if marker < 0:
            raise RuntimeError("ZIP64 EOCD not found")
        zip64 = struct.unpack_from("<4sQ2H2L4Q", zip64_data, marker)
        central_size = zip64[8]
        central_offset = zip64[9]

    central = request_range(url, central_offset, central_offset + central_size - 1)
    rows = []
    pos = 0
    while pos < len(central):
        (
            _sig,
            _made,
            _needed,
            _flags,
            compression,
            _mtime,
            _mdate,
            _crc,
            comp_size,
            uncomp_size,
            name_len,
            extra_len,
            comment_len,
            _disk,
            _internal,
            _external,
            local_offset,
        ) = struct.unpack_from("<4s6H3L5H2L", central, pos)
        name = central[pos + 46 : pos + 46 + name_len].decode("utf-8", "replace")
        extra = central[
            pos + 46 + name_len : pos + 46 + name_len + extra_len
        ]
        values = parse_zip64_extra(
            extra,
            {
                "uncompressed_size": uncomp_size,
                "compressed_size": comp_size,
                "local_header_offset": local_offset,
            },
        )
        rows.append(
            {
                "name": name,
                "compression": compression,
                **values,
            }
        )
        pos += 46 + name_len + extra_len + comment_len
    return rows


def write_member_gzip(url: str, row: dict, destination: Path) -> None:
    offset = row["local_header_offset"]
    header = request_range(url, offset, offset + 29)
    (
        _sig,
        _version,
        _flags,
        _compression,
        _mtime,
        _mdate,
        _crc,
        _comp_size,
        _uncomp_size,
        filename_len,
        extra_len,
    ) = struct.unpack("<4s5H3L2H", header)
    data_start = offset + 30 + filename_len + extra_len
    data_end = data_start + row["compressed_size"] - 1
    response = requests.get(
        url,
        headers={"Range": f"bytes={data_start}-{data_end}"},
        stream=True,
        timeout=(30, 300),
    )
    response.raise_for_status()
    if response.status_code != 206:
        raise RuntimeError(f"Range request returned {response.status_code}")

    temporary = destination.with_suffix(destination.suffix + ".part")
    try:
        with gzip.open(temporary, "wb") as output:
            if row["compression"] == 8:
                decompressor = zlib.decompressobj(-15)
                for chunk in response.iter_content(chunk_size=1_048_576):
                    if chunk:
                        output.write(decompressor.decompress(chunk))
                output.write(decompressor.flush())
            elif row["compression"] == 0:
                for chunk in response.iter_content(chunk_size=1_048_576):
                    if chunk:
                        output.write(chunk)
            else:
                raise RuntimeError(
                    f"Unsupported compression method {row['compression']}"
                )
        temporary.replace(destination)
    finally:
        if temporary.exists():
            temporary.unlink()


def read_network(path: Path) -> pd.DataFrame:
    columns = list(pd.read_csv(path, sep="\t", nrows=0).columns)
    if "TF" in columns and "GENE" in columns and "pred_prob" in columns:
        frame = pd.read_csv(
            path,
            sep="\t",
            usecols=["TF", "GENE", "pred_prob"],
            dtype={"TF": str, "GENE": str, "pred_prob": float},
            low_memory=False,
        ).rename(columns={"TF": "source", "GENE": "target", "pred_prob": "score"})
    else:
        frame = pd.read_csv(
            path,
            sep="\t",
            header=None,
            names=["source", "target", "score"],
            dtype={"source": str, "target": str, "score": float},
            low_memory=False,
        )
    frame["source"] = frame["source"].astype(str).str.split(".").str[0]
    frame["target"] = frame["target"].astype(str).str.split(".").str[0]
    frame["score"] = pd.to_numeric(frame["score"], errors="coerce")
    return frame.dropna(subset=["source", "target", "score"]).drop_duplicates()


def map_symbols_to_ids(symbols: list[str]) -> dict[str, str]:
    if not symbols:
        return {}
    response = requests.post(
        "https://mygene.info/v3/query",
        data={
            "q": ",".join(symbols),
            "scopes": "symbol",
            "fields": "symbol,ensembl.gene",
            "species": "human",
        },
        timeout=120,
    )
    response.raise_for_status()
    result = response.json()
    mapping = {}
    for record in result:
        symbol = record.get("symbol")
        ensembl = record.get("ensembl")
        if isinstance(ensembl, dict):
            ensembl = [ensembl]
        if not isinstance(ensembl, list):
            continue
        for item in ensembl:
            gene_id = item.get("gene") if item else None
            if symbol and gene_id:
                mapping[gene_id] = symbol
    return mapping


def map_ids_to_symbols(gene_ids: list[str]) -> dict[str, str]:
    mapping = {}
    for offset in range(0, len(gene_ids), 800):
        batch = gene_ids[offset : offset + 800]
        response = requests.post(
            "https://mygene.info/v3/query",
            data={
                "q": ",".join(batch),
                "scopes": "ensembl.gene",
                "fields": "symbol,ensembl.gene",
                "species": "human",
            },
            timeout=120,
        )
        response.raise_for_status()
        for record in response.json():
            symbol = record.get("symbol")
            ensembl = record.get("ensembl")
            if isinstance(ensembl, dict):
                ensembl = [ensembl]
            if not isinstance(ensembl, list):
                continue
            for item in ensembl:
                gene_id = item.get("gene") if item else None
                if symbol and gene_id:
                    mapping[gene_id] = symbol
    return mapping


def main() -> None:
    RESULTS.mkdir(parents=True, exist_ok=True)
    metadata = requests.get(RECORD_API, timeout=60).json()
    compendium = next(
        item for item in metadata["files"] if item["key"] == "METANet_compendium.zip"
    )
    total_size = compendium["size"]
    directory = central_directory(URL, total_size)

    selected = {
        "whole_blood_full": next(
            row
            for row in directory
            if row["name"] == "METANet_compendium/Networks/METANets/whole_blood.tsv"
        ),
        "whole_blood_tissue_specific": next(
            row
            for row in directory
            if row["name"]
            == "METANet_compendium/Networks/Tissue_Specific_METANets/whole_blood.tsv"
        ),
        "lung_tissue_specific": next(
            row
            for row in directory
            if row["name"]
            == "METANet_compendium/Networks/Tissue_Specific_METANets/lung.tsv"
        ),
    }

    networks = {}
    manifest = []
    for label, row in selected.items():
        network_path = RESULTS / f"metanet_{label}.tsv.gz"
        if not network_path.exists():
            write_member_gzip(URL, row, network_path)
        networks[label] = read_network(network_path)
        manifest.append(
            {
                "network": label,
                "member": row["name"],
                "compressed_bytes": row["compressed_size"],
                "uncompressed_bytes": row["uncompressed_size"],
                "edges": len(networks[label]),
            }
        )

    pd.DataFrame(manifest).to_csv(
        RESULTS / "metanet_manifest.tsv", sep="\t", index=False
    )

    id_to_symbol = map_symbols_to_ids(sorted(TARGET_GENES))
    target_ids = {
        symbol: gene_id
        for gene_id, symbol in id_to_symbol.items()
        if symbol in TARGET_GENES
    }
    id_symbols = {value: key for key, value in target_ids.items()}

    direct_rows = []
    upstream_rows = []
    overlap_rows = []

    module_files = {
        "Restoration_Core": RESULTS / "restoration_core_genes.txt",
        "Disease_Driver": RESULTS / "disease_driver_genes.txt",
    }
    module_symbols = {
        label: [line.strip() for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]
        for label, path in module_files.items()
    }
    module_symbol_map = {
        label: map_symbols_to_ids(genes) for label, genes in module_symbols.items()
    }

    for network_name, network in networks.items():
        target_universe = set(network["target"])
        for tf in sorted(CORE_TFS):
            tf_id = target_ids.get(tf)
            if tf_id is None:
                continue
            outgoing = network.loc[network["source"] == tf_id].sort_values(
                "score", ascending=False, ignore_index=True
            )
            outgoing_targets = set(outgoing["target"])
            for target in sorted(TARGET_GENES):
                target_id = target_ids.get(target)
                if target_id is None:
                    continue
                hit = outgoing.loc[outgoing["target"] == target_id]
                direct_rows.append(
                    {
                        "network": network_name,
                        "TF": tf,
                        "target": target,
                        "edge_present": not hit.empty,
                        "score": float(hit["score"].iloc[0]) if not hit.empty else math.nan,
                        "rank_within_TF": int(hit.index[0] + 1) if not hit.empty else math.nan,
                        "outgoing_edges": len(outgoing),
                    }
                )

            for module_name, symbol_map in module_symbol_map.items():
                module_ids = set(symbol_map) & target_universe
                overlap = outgoing_targets & module_ids
                universe_n = len(target_universe)
                regulon_n = len(outgoing_targets & target_universe)
                module_n = len(module_ids)
                overlap_n = len(overlap)
                if universe_n == 0 or regulon_n == 0 or module_n == 0:
                    p_value = math.nan
                else:
                    p_value = hypergeom_sf(
                        overlap_n - 1, universe_n, module_n, regulon_n
                    )
                overlap_rows.append(
                    {
                        "network": network_name,
                        "TF": tf,
                        "module": module_name,
                        "overlap": overlap_n,
                        "TF_targets": regulon_n,
                        "module_genes_in_network": module_n,
                        "network_genes": universe_n,
                        "fold_enrichment": (
                            (overlap_n + 0.5) / (regulon_n + 1)
                        )
                        / ((module_n + 0.5) / (universe_n + 1)),
                        "p_value": p_value,
                    }
                )

        core_target_ids = {
            target_ids[tf]: tf for tf in CORE_TFS if tf in target_ids
        }
        incoming = network.loc[network["target"].isin(core_target_ids)].copy()
        if not incoming.empty:
            source_ids = incoming["source"].drop_duplicates().tolist()
            source_map = map_ids_to_symbols(source_ids)
            incoming["source_symbol"] = incoming["source"].map(source_map)
            incoming["target_symbol"] = incoming["target"].map(core_target_ids)
            incoming = incoming.sort_values(
                ["target_symbol", "score"], ascending=[True, False]
            )
            for target_symbol, group in incoming.groupby("target_symbol"):
                for rank, (_, row) in enumerate(group.head(15).iterrows(), start=1):
                    upstream_rows.append(
                        {
                            "network": network_name,
                            "target_TF": target_symbol,
                            "candidate_upstream_TF": row["source_symbol"],
                            "score": row["score"],
                            "rank": rank,
                            "incoming_edges": len(group),
                        }
                    )

    direct = pd.DataFrame(direct_rows)
    upstream = pd.DataFrame(upstream_rows)
    overlap = pd.DataFrame(overlap_rows)
    if not overlap.empty:
        overlap["FDR"] = bh_adjust(overlap["p_value"].tolist())
        overlap = overlap.sort_values(["module", "FDR", "p_value"])

    direct.to_csv(RESULTS / "metanet_direct_edges.tsv", sep="\t", index=False)
    upstream.to_csv(RESULTS / "metanet_upstream_candidates.tsv", sep="\t", index=False)
    overlap.to_csv(RESULTS / "metanet_module_overlap.tsv", sep="\t", index=False)


def hypergeom_sf(k: int, population: int, successes: int, draws: int) -> float:
    """Upper-tail hypergeometric probability via log-combinations."""
    if k < 0:
        return 1.0
    lower = max(0, draws - (population - successes))
    upper = min(draws, successes)
    if k >= upper:
        return 0.0
    denominator = log_comb(population, draws)
    probability = 0.0
    for observed in range(k + 1, upper + 1):
        probability += math.exp(
            log_comb(successes, observed)
            + log_comb(population - successes, draws - observed)
            - denominator
        )
    return min(1.0, probability)


def log_comb(n: int, k: int) -> float:
    if k < 0 or k > n:
        return -math.inf
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)


def bh_adjust(p_values: list[float]) -> list[float]:
    indexed = [
        (value, index)
        for index, value in enumerate(p_values)
        if not math.isnan(value)
    ]
    adjusted = [math.nan] * len(p_values)
    if not indexed:
        return adjusted
    indexed.sort()
    m = len(indexed)
    running = 1.0
    for reverse_rank, (value, original_index) in enumerate(reversed(indexed), start=1):
        rank = m - reverse_rank + 1
        running = min(running, value * m / rank)
        adjusted[original_index] = min(1.0, running)
    return adjusted


if __name__ == "__main__":
    main()
