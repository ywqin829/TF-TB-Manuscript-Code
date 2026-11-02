#!/usr/bin/env python3
"""Normalize openxlsx output for strict XLSX readers while retaining every table cell.

Run after scripts/07_export_tables.R. The R export has stale worksheet dimensions
and references absent drawing parts; Excel may repair these silently, but Python
readers either truncate sheets or reject the file. Rewriting through openpyxl
also keeps the public supplementary workbook portable.
"""
from pathlib import Path
from tempfile import NamedTemporaryFile
import os

from openpyxl import Workbook, load_workbook
from openpyxl.styles import Alignment, Font, PatternFill
from openpyxl.utils import get_column_letter


root = Path(__file__).resolve().parents[1]
path = root / "supplement" / "Supplementary_Tables.xlsx"
source = load_workbook(path, read_only=True, data_only=True)
target = Workbook()
target.remove(target.active)
header_fill = PatternFill("solid", fgColor="1F4E79")
for name in source.sheetnames:
    sheet = source[name]
    sheet.reset_dimensions()  # openxlsx declares A1 even when data extend farther.
    output = target.create_sheet(name)
    for row in sheet.iter_rows(values_only=True):
        output.append(row)
    output.freeze_panes = "A2"
    if output.max_column:
        output.auto_filter.ref = output.dimensions
    for cell in output[1]:
        cell.fill = header_fill
        cell.font = Font(color="FFFFFF", bold=True)
        cell.alignment = Alignment(wrap_text=True)
    for col in output.columns:
        width = min(55, max(12, max(len(str(c.value or "")) for c in list(col)[:100]) + 2))
        output.column_dimensions[get_column_letter(col[0].column)].width = width
source.close()
with NamedTemporaryFile(dir=path.parent, suffix=".xlsx", delete=False) as handle:
    temporary = Path(handle.name)
try:
    target.save(temporary)
    os.replace(temporary, path)
finally:
    temporary.unlink(missing_ok=True)
print(f"Finalized {len(target.sheetnames)} sheets: {path}")
