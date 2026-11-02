#!/usr/bin/env python3
import os
"""Check the exact prepared submission files for internal notes and wrong files."""

from pathlib import Path
import re

import fitz
from docx import Document

root = Path(os.environ.get("TF_TB_PROJECT_ROOT", Path.cwd()))
package = root / "submission_package_20261030"
report = root / "docs" / "Submission_Preflight_20261031.txt"
expected = {
    "Manuscript_Inflammation_Research.docx", "Title_Page.docx",
    "Cover_Letter.docx", "Supplementary_Information.docx",
    *(f"Figures/Figure{i}_{stem}.pdf" for i, stem in [
        (1, "StudyDesign"), (2, "Discovery"), (3, "SingleCellContext"),
        (4, "OutcomeModels"), (5, "TBRI_OutcomeAssociation"),
        (6, "CompositionConfounding"), (7, "InVivo_Candidate_TF_Context"),
    ]),
    "Supplementary/FigureS1_METANet.pdf",
    "Supplementary/FigureS2_CellTypeSpecificMtbResponse.pdf",
    "Supplementary/FigureS3_TBRI_Bcell_TimeAdjusted.pdf",
    "Supplementary/FigureS4_GSE31348_CIBERSORT_LM22_Bcell_Trajectory.pdf",
}
dialogue = re.compile(
    r"system prompt|\buser\s*:|\bassistant\s*:|conversation log|as an AI|"
    r"以下是修改后的|我来为你|请将以下内容|希望这对你有帮助|"
    r"\b(?:Claude|Codex|ChatGPT)\b", re.I,
)
working = re.compile(
    r"\b(?:TODO|TBD|FIXME|placeholder)\b|\[INSERT|\[FILL IN|"
    r"待办|待补|待定|待确认|工作稿|内部备注|草稿备注|"
    r"update before submission|at submission|to be completed|citation needed", re.I,
)


def docx_lines(path):
    doc = Document(path)
    for p in doc.paragraphs:
        yield p.text
    for table in doc.tables:
        for row in table.rows:
            for cell in row.cells:
                yield cell.text
    for section in doc.sections:
        for part in (section.header, section.footer):
            for p in part.paragraphs:
                yield p.text
            for table in part.tables:
                for row in table.rows:
                    for cell in row.cells:
                        yield cell.text


actual = {p.relative_to(package).as_posix() for p in package.rglob("*") if p.is_file()}
problems = []
for name in sorted(expected - actual):
    problems.append(f"MISSING: {name}")
for name in sorted(actual - expected):
    problems.append(f"UNEXPECTED: {name}")
for name in sorted(actual & expected):
    path = package / name
    if path.stat().st_size == 0:
        problems.append(f"EMPTY: {name}")
        continue
    if path.suffix.lower() == ".docx":
        lines = docx_lines(path)
    else:
        with fitz.open(path) as pdf:
            lines = [page.get_text() for page in pdf]
    for number, text in enumerate(lines, 1):
        for label, pattern in (("INTERNAL_DIALOGUE", dialogue), ("WORKING_NOTE", working)):
            match = pattern.search(text)
            if match:
                excerpt = " ".join(text[max(0, match.start()-55):match.end()+90].split())
                problems.append(f"{label}: {name}, unit {number}: {excerpt}")

status = "PASS" if not problems else "BLOCKED"
report.write_text(
    f"Submission preflight: {status}\n"
    f"Exact prepared package: {package}\n"
    f"Expected files: {len(expected)}; actual files: {len(actual)}\n"
    "Checks: internal prompt/dialogue; unfinished notes or markers; missing, stale or unexpected upload files.\n"
    + ("\n".join(problems) + "\n" if problems else "No hits.\n"),
    encoding="utf-8",
)
print(report.read_text(encoding="utf-8"))
raise SystemExit(0 if status == "PASS" else 1)


