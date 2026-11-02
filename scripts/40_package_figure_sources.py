#!/usr/bin/env python3
"""Package the canonical R figure sources for review and GitHub preparation."""

from pathlib import Path
from zipfile import ZipFile, ZIP_DEFLATED

root = Path(__file__).resolve().parents[1]
archive_path = root / "TF_TB_Figure_R_Sources_20261031.zip"
source_dir = root / "scripts" / "figure_sources"

files = sorted(source_dir.glob("*.R")) + [source_dir / "README.md"]
files += [root / "scripts" / name for name in (
    "00_config.R", "06_figures.R", "20_gse326212_invivo_figure.R",
    "37_tbri_bcell_time_adjusted.R",
    "38_gse31348_cibersort_lm22_sensitivity.R", "40_editorial_core_validation.R",
    "39_export_integrated_figure_sources.py",
)]
files.append(root / "inputs" / "figure_source_context" / "Fig3F_GSE192483_3TF_UMAP.pdf")

missing = [str(path) for path in files if not path.is_file()]
if missing:
    raise FileNotFoundError("Missing figure source inputs: " + ", ".join(missing))

with ZipFile(archive_path, "w", ZIP_DEFLATED) as archive:
    for path in files:
        archive.write(path, path.relative_to(root).as_posix())

print(f"Packaged {len(files)} source files: {archive_path}")
