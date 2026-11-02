#!/usr/bin/env python3
import os
"""Inventory the canonical review and upload outputs using the ASCII project path."""

import csv
from pathlib import Path

root = Path(os.environ.get("TF_TB_PROJECT_ROOT", Path.cwd()))
rows = []


def add(path, purpose):
    if not path.is_file():
        raise FileNotFoundError(path)
    rows.append((str(path), str(path.parent), purpose))


for path in sorted((root / "figures").glob("*.pdf")):
    purpose = "Main manuscript figure" if not path.name.startswith("FigureS") else "Supplementary result figure"
    add(path, purpose)
for path in sorted((root / "scripts" / "figure_sources").glob("*.R")):
    add(path, "R entry point for the matching integrated PDF figure")
for name in ("README.md",):
    add(root / "scripts" / "figure_sources" / name, "Figure source reproduction instructions")
for pattern in ("gse31348_cibersort*.tsv", "gse89403_tbri_bcell_time_adjusted*.tsv"):
    for path in sorted((root / "results").glob(pattern)):
        add(path, "Sensitivity analysis result or retained model input")
for pattern in ("gse31348_cibersort*.txt", "gse89403_tbri_bcell*.txt"):
    for path in sorted((root / "logs").glob(pattern)):
        add(path, "Sensitivity analysis method and result summary")
for name in (
    "06_figures.R", "08_build_manuscript.py", "21_sync_markdown.py",
    "22_build_submission_package.py", "37_tbri_bcell_time_adjusted.R",
    "38_gse31348_cibersort_lm22_sensitivity.R",
    "39_export_integrated_figure_sources.py", "40_package_figure_sources.py",
    "41_submission_preflight.py", "42_final_delivery_manifest.py",
):
    add(root / "scripts" / name, "Canonical generation, packaging or verification script")
for name, purpose in (
    ("Manuscript_Inflammation_Research.docx", "Manuscript with embedded Table 1"),
    ("Title_Page.docx", "Title page"),
    ("Cover_Letter.docx", "Cover letter"),
    ("Supplementary_Information.docx", "Captions and embedded supplementary figures"),
    ("Figure_Editorial_Audit.md", "Definitive main and supplementary figure map"),
    ("BioRender_Figure1_Prompt.md", "BioRender design prompt for Figure 1"),
    ("Submission_Preflight_20261031.txt", "Three-part exact-file submission preflight"),
    ("SUBMISSION_GUIDE.md", "Upload field guide"),
):
    add(root / "docs" / name, purpose)
add(root / "TF_TB_Figure_R_Sources_20261031.zip", "R source bundle prepared for GitHub review")
add(root / "TF_TB_Submission_20261030.zip", "15-file prepared manuscript upload archive")

manifest = root / "docs" / "Final_Delivery_Manifest_20261031.tsv"
with manifest.open("w", encoding="utf-8", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t")
    writer.writerow(("Absolute_path", "Parent_output_directory", "Purpose"))
    writer.writerows(rows)
print(f"Wrote {len(rows)} output paths: {manifest}")

