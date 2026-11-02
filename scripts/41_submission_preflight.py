#!/usr/bin/env python3
"""Run the three submission-package checks required before upload."""
import os
import re
from pathlib import Path

import fitz
from docx import Document

root = Path(os.environ.get("TF_TB_PROJECT_ROOT", Path.cwd()))
package = root / "submission_package_20261102"
report = root / "docs" / "Submission_Preflight_20261102.txt"
expected = {
    "Manuscript_Inflammation_Research.docx",
    "Title_Page.docx",
    "Cover_Letter.docx",
    *(f"Figures/Figure{i}_{stem}.pdf" for i, stem in [
        (1, "StudyDesign"), (2, "Discovery"), (3, "CompositionConfounding"),
        (4, "SingleCellContext"), (5, "OutcomeModels"),
        (6, "TBRI_OutcomeAssociation"),
    ]),
    *(f"Supplementary/FigureS{i}_{stem}.pdf" for i, stem in [
        (1, "METANet"), (2, "CellTypeSpecificMtbResponse"),
        (3, "TBRI_Bcell_TimeAdjusted"),
        (4, "GSE31348_CIBERSORT_LM22_Bcell_Trajectory"),
        (5, "InVivo_Candidate_TF_Context"),
    ]),
    "Supplementary_Code.zip",
}

dialogue_terms = [
    "system" + " prompt", r"\b" + "user" + r"\\s*:",
    r"\b" + "assistant" + r"\\s*:", "conversation" + " log",
    "as an " + "AI", "以下是修改后的", "我来为你", "请将以下内容",
    r"\\b(?:" + "|".join(["Claude", "Codex", "ChatGPT"]) + r")\\b",
]
unfinished_terms = [
    r"\\b(?:" + "|".join(["TO" + "DO", "T" + "BD", "FIX" + "ME"]) + r")\\b",
    "\\\\" + "[" + "INSERT", "\\\\" + "[" + "FILL IN", "place" + "holder",
    "待办", "待补", "待定", "待确认", "工作稿", "内部备注", "草稿备注",
    "update before submission", "at submission", "to be completed", "citation needed",
]
dialogue = re.compile("|".join(dialogue_terms), re.I)
unfinished = re.compile("|".join(unfinished_terms), re.I)


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


actual = {p.relative_to(package).as_posix() for p in package.rglob("*") if p.is_file()} if package.exists() else set()
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
    elif path.suffix.lower() == ".pdf":
        with fitz.open(path) as pdf:
            lines = [page.get_text() for page in pdf]
    else:
        continue
    for number, text in enumerate(lines, 1):
        for label, pattern in (("INTERNAL_DIALOGUE", dialogue), ("UNFINISHED_MARKER", unfinished)):
            match = pattern.search(text)
            if match:
                excerpt = " ".join(text[max(0, match.start()-55):match.end()+90].split())
                problems.append(f"{label}: {name}, unit {number}: {excerpt}")

status = "PASS" if not problems else "BLOCKED"
report.parent.mkdir(parents=True, exist_ok=True)
report.write_text(
    f"Submission preflight: {status}\n"
    f"Exact prepared package: {package}\n"
    f"Expected files: {len(expected)}; actual files: {len(actual)}\n"
    "Checks: internal dialogue; unfinished markers; missing, stale or unexpected upload files.\n"
    + ("\n".join(problems) + "\n" if problems else "No hits.\n"),
    encoding="utf-8",
)
print(report.read_text(encoding="utf-8"))
raise SystemExit(0 if status == "PASS" else 1)
