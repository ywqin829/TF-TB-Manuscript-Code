#!/usr/bin/env python3
"""Run the mandatory manuscript-no-internal-notes gate on the upload set."""

from __future__ import annotations

import json
import re
import zipfile
from pathlib import Path


PROJECT = Path(__file__).resolve().parents[1]
PACKAGE = PROJECT / "submission_package_20261030"
REPORT = PROJECT / "docs" / "PRE_SUBMISSION_GATE.md"

EXPECTED = {
    "Cover_Letter.docx",
    "Manuscript_Inflammation_Research.docx",
    "Supplementary_Information.docx",
    "Title_Page.docx",
    "Figures/Figure1_StudyDesign.pdf",
    "Figures/Figure2_Discovery.pdf",
    "Figures/Figure3_SingleCellContext.pdf",
    "Figures/Figure4_OutcomeModels.pdf",
    "Figures/Figure5_TBRI_OutcomeAssociation.pdf",
    "Figures/Figure6_CompositionConfounding.pdf",
    "Figures/Figure7_InVivo_Candidate_TF_Context.pdf",
    "Supplementary/FigureS1_METANet.pdf",
    "Supplementary/FigureS2_GSE218065_descriptive.pdf",
    "Supplementary/FigureS3_CellTypeSpecificMtbResponse.pdf",
    "Supplementary/FigureS4_TBRI_Calibration_DCA.pdf",
    "Supplementary/FigureS5_TBRI_Bcell_TimeAdjusted.pdf",
    "Supplementary/FigureS6_GSE31348_CIBERSORT_LM22_Bcell_Trajectory.pdf",
}

INTERNAL_PATTERNS = [
    r"system prompt",
    r"\buser:",
    r"\bassistant:",
    r"conversation log",
    r"as an AI",
    r"以下是修改后的",
    r"我来为你",
    r"请将以下内容",
    r"希望这对你有帮助",
]

WORKING_PATTERNS = [
    r"\bTODO\b",
    r"\bTBD\b",
    r"\bFIXME\b",
    r"\bplaceholder\b",
    r"\[INSERT",
    r"\[FILL IN",
    r"待办",
    r"待补",
    r"待定",
    r"待确认",
    r"工作稿",
    r"内部备注",
    r"草稿备注",
    r"update before submission",
    r"at submission",
    r"to be completed",
    r"citation needed",
]


def xml_text_from_office(path: Path) -> str:
    chunks = []
    with zipfile.ZipFile(path) as archive:
        for name in archive.namelist():
            if not name.endswith(".xml"):
                continue
            data = archive.read(name).decode("utf-8", errors="ignore")
            chunks.append(re.sub(r"<[^>]+>", " ", data))
    return "\n".join(chunks)


def code_text_from_zip(path: Path) -> str:
    chunks = []
    with zipfile.ZipFile(path) as archive:
        for name in archive.namelist():
            if name.endswith("/"):
                continue
            if Path(name).suffix.lower() in {
                ".r", ".py", ".ps1", ".md", ".txt", ".tsv", ".csv", ".json", ".yml", ".yaml"
            }:
                chunks.append(archive.read(name).decode("utf-8", errors="ignore"))
    return "\n".join(chunks)


def matches(patterns: list[str], text: str) -> list[str]:
    return [pattern for pattern in patterns if re.search(pattern, text, flags=re.I)]


def main() -> None:
    actual = {
        str(path.relative_to(PACKAGE)).replace("\\", "/")
        for path in PACKAGE.rglob("*")
        if path.is_file()
    }

    missing = sorted(EXPECTED - actual)
    extra = sorted(actual - EXPECTED)

    stale_pattern = re.compile(r"(\.bak|\.tmp|~$|_DO_NOT_UPLOAD|old|superseded|draft)", re.I)
    stale = sorted(name for name in actual if stale_pattern.search(name))

    internal_hits: list[dict] = []
    working_hits: list[dict] = []

    for relative in sorted(EXPECTED & actual):
        path = PACKAGE / relative
        suffix = path.suffix.lower()
        if suffix in {".docx", ".xlsx"}:
            text = xml_text_from_office(path)
        elif suffix == ".zip":
            text = code_text_from_zip(path)
        elif suffix in {".txt", ".md", ".csv", ".tsv"}:
            text = path.read_text(encoding="utf-8", errors="ignore")
        else:
            text = ""

        for pattern in matches(INTERNAL_PATTERNS, text):
            internal_hits.append({"file": relative, "pattern": pattern})
        for pattern in matches(WORKING_PATTERNS, text):
            working_hits.append({"file": relative, "pattern": pattern})

    passed = not missing and not extra and not stale and not internal_hits and not working_hits

    lines = [
        "# Pre-submission Gate",
        "",
        f"Result: {'PASS' if passed else 'BLOCKED'}",
        "",
        "## 1. Internal prompts or dialogue",
        "",
        "PASS" if not internal_hits else json.dumps(internal_hits, indent=2),
        "",
        "## 2. TODO or working notes",
        "",
        "PASS" if not working_hits else json.dumps(working_hits, indent=2),
        "",
        "## 3. Wrong, stale, or missing upload files",
        "",
        f"Missing: {missing if missing else 'none'}",
        "",
        f"Unexpected: {extra if extra else 'none'}",
        "",
        f"Stale names: {stale if stale else 'none'}",
        "",
        f"Upload files checked: {len(actual)}",
    ]
    REPORT.write_text("\n".join(lines) + "\n", encoding="utf-8")

    if not passed:
        raise SystemExit("Pre-submission gate blocked")


if __name__ == "__main__":
    main()
