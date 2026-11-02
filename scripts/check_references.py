#!/usr/bin/env python3
"""Cross-check the existing numbered references against Crossref."""

from __future__ import annotations

import re
import time
from pathlib import Path

import docx
import pandas as pd
import requests


PROJECT = Path(__file__).resolve().parents[1]
MANUSCRIPT = Path(
    r"C:\Users\yw_qi\OneDrive\TF-Manuscript-inflammation research(2)(1).docx"
)


def main() -> None:
    document = docx.Document(MANUSCRIPT)
    references = []
    for paragraph in document.paragraphs:
        match = re.match(r"^(\d+)\.\s+(.*)$", paragraph.text.strip())
        if match:
            references.append((int(match.group(1)), match.group(2)))

    rows = []
    for number, reference in references:
        try:
            response = requests.get(
                "https://api.crossref.org/works",
                params={"query.bibliographic": reference, "rows": 1},
                timeout=30,
            )
            response.raise_for_status()
            item = response.json()["message"]["items"][0]
            rows.append(
                {
                    "reference_number": number,
                    "original_reference": reference,
                    "crossref_score": item.get("score"),
                    "matched_title": (item.get("title") or [""])[0],
                    "doi": item.get("DOI", ""),
                }
            )
        except Exception as exc:
            rows.append(
                {
                    "reference_number": number,
                    "original_reference": reference,
                    "crossref_score": None,
                    "matched_title": f"ERROR: {exc}",
                    "doi": "",
                }
            )
        time.sleep(0.2)

    output = PROJECT / "results" / "reference_crossref_candidates.tsv"
    pd.DataFrame(rows).to_csv(output, sep="\t", index=False)


if __name__ == "__main__":
    main()
