#!/usr/bin/env python3
"""Automated structural and integrity checks for the submission package."""

from __future__ import annotations

import re
from pathlib import Path

import docx
from openpyxl import load_workbook
from PIL import Image


PROJECT = Path(__file__).resolve().parents[1]
DOCS = PROJECT / "docs"
FIGS = PROJECT / "figures"
SUPP = PROJECT / "supplement"

checks: list[tuple[str, bool, str]] = []


def check(name: str, condition: bool, detail: str = "") -> None:
    checks.append((name, bool(condition), detail))


def document_text(path: Path) -> str:
    document = docx.Document(path)
    return "\n".join(paragraph.text for paragraph in document.paragraphs)


MANUSCRIPT = DOCS / "Manuscript_Inflammation_Research.docx"
TITLE_PAGE = DOCS / "Title_Page.docx"
COVER = DOCS / "Cover_Letter.docx"
SUPPLEMENT_TEXT = DOCS / "Supplementary_Information.docx"

for path in [MANUSCRIPT, TITLE_PAGE, COVER, SUPPLEMENT_TEXT]:
    check(f"Exists: {path.name}", path.exists(), str(path))

if MANUSCRIPT.exists():
    manuscript = docx.Document(MANUSCRIPT)
    paragraphs = [paragraph.text for paragraph in manuscript.paragraphs]
    full_text = "\n".join(paragraphs)

    abstract_start = paragraphs.index("Abstract")
    keyword_index = next(
        index for index, text in enumerate(paragraphs) if text.startswith("Keywords:")
    )
    abstract = " ".join(paragraphs[abstract_start + 1 : keyword_index])
    abstract_words = len(re.findall(r"\b[\w-]+\b", abstract))
    check(
        "Abstract length 150-200 words",
        150 <= abstract_words <= 200,
        f"{abstract_words} words",
    )

    narrative_text = full_text.split("\nReferences\n", 1)[0]
    narrative_word_count = len(re.findall(r"\b[\w-]+\b", narrative_text))
    check(
        "Narrative text under approximately 6,000 words excluding references",
        narrative_word_count < 6000,
        f"{narrative_word_count} words including title, abstract, body and legends",
    )

    check(
        "Seven main figures embedded",
        len(manuscript.inline_shapes) == 7,
        f"{len(manuscript.inline_shapes)} inline images",
    )
    check(
        "Main cohort table contains three cohorts",
        bool(manuscript.tables)
        and len(manuscript.tables[0].rows) == 4
        and len(manuscript.tables[0].columns) == 7,
        f"{len(manuscript.tables)} table(s)",
    )

    body = full_text.split("\nReferences\n", 1)[0]
    citation_groups = re.findall(r"\[([0-9,\s\-]+)\]", body)
    citation_numbers = {
        int(number)
        for group in citation_groups
        for number in re.findall(r"\d+", group)
    }
    reference_numbers = {
        int(match.group(1))
        for match in re.finditer(r"(?m)^(\d+)\.\s", full_text)
    }
    max_group_size = max(
        (len(re.findall(r"\d+", group)) for group in citation_groups),
        default=0,
    )
    check(
        "Citation groups contain at most three references",
        max_group_size <= 3,
        f"largest group = {max_group_size}",
    )
    range_hits = re.findall(r"[\[(]\s*\d+\s*[-–]\s*\d+\s*[\])]", body)
    check(
        "No numeric citation ranges",
        not range_hits,
        f"hits = {range_hits}",
    )
    paragraphs_with_stacked_groups = []
    for paragraph in paragraphs:
        groups = re.findall(r"\[([0-9,\s\-–]+)\]", paragraph)
        stacked = [
            group for group in groups
            if len(re.findall(r"\d+", group)) > 1
        ]
        if len(stacked) >= 2:
            paragraphs_with_stacked_groups.append(stacked)
    check(
        "No paragraph contains two or more stacked citation groups",
        not paragraphs_with_stacked_groups,
        f"hits = {paragraphs_with_stacked_groups}",
    )
    check(
        "Every reference is cited",
        reference_numbers.issubset(citation_numbers),
        f"uncited = {sorted(reference_numbers - citation_numbers)}",
    )
    check(
        "Every citation has a reference entry",
        citation_numbers.issubset(reference_numbers),
        f"missing = {sorted(citation_numbers - reference_numbers)}",
    )

    forbidden = [
        "TODO",
        "FIXME",
        "work note",
        "working dialogue",
        "internal prompt",
        "selection leakage",
        "the previous analysis was wrong",
        "do not submit",
    ]
    hits = [token for token in forbidden if token.lower() in full_text.lower()]
    check(
        "No internal notes or placeholders in manuscript",
        not hits,
        f"hits = {hits}",
    )

    check(
        "No em dash used as a sentence connector",
        "\u2014" not in body,
        f"count = {body.count(chr(0x2014))}",
    )
    banned_ai_terms = [
        "innovative", "pioneering", "revolutionary paradigm",
        "transformative framework", "remarkable", "unprecedented",
        "breakthrough performance", "general-purpose", "at its essence",
        "pave the way for", "profound challenges", "stems from",
    ]
    banned_hits = {
        term: body.lower().count(term)
        for term in banned_ai_terms
        if body.lower().count(term) >= 3
    }
    check(
        "Banned AI-tone vocabulary below threshold",
        not banned_hits,
        f"hits = {banned_hits}",
    )

for path in [COVER, SUPPLEMENT_TEXT, TITLE_PAGE]:
    if not path.exists():
        continue
    text = document_text(path)
    forbidden = ["TODO", "FIXME", "working dialogue", "internal prompt", "do not submit"]
    hits = [token for token in forbidden if token.lower() in text.lower()]
    check(
        f"No internal notes or placeholders in {path.name}",
        not hits,
        f"hits = {hits}",
    )

required_figures = [
    "Figure1_StudyDesign",
    "Figure2_Discovery",
    "Figure3_SingleCellContext",
    "Figure4_OutcomeModels",
    "Figure5_TBRI_OutcomeAssociation",
    "Figure6_CompositionConfounding",
    "Figure7_InVivo_Candidate_TF_Context",
    "FigureS3_CellTypeSpecificMtbResponse",
    "FigureS1_METANet",
    "FigureS2_GSE218065_descriptive",
    "FigureS4_TBRI_Calibration_DCA",
    "FigureS_time_adjusted_tbri_bcell_GSE89403",
    "FigureS_GSE31348_CIBERSORT_LM22_Bcell_trajectory",
]
for stem in required_figures:
    for suffix in [".pdf"]:
        path = FIGS / f"{stem}{suffix}"
        check(f"Exists: {stem}{suffix}", path.exists(), str(path))
        if path.exists() and path.stat().st_size < 1_000:
            check(f"Non-empty: {stem}{suffix}", False, f"{path.stat().st_size} bytes")

workbook_path = SUPP / "Supplementary_Tables.xlsx"
check("Supplementary workbook exists", workbook_path.exists(), str(workbook_path))
if workbook_path.exists():
    workbook = load_workbook(workbook_path, read_only=True)
    check(
        "Supplementary workbook has 51 result sheets",
        len(workbook.sheetnames) == 51,
        ", ".join(workbook.sheetnames),
    )

required_paths = [
    PROJECT / "results" / "analysis_evidence_map.tsv",
    PROJECT / "logs" / "internal_analysis_amendments.tsv",
    PROJECT / "logs" / "gse40553_discovery_summary.txt",
    PROJECT / "logs" / "gse89403_outcome_summary.txt",
]
for path in required_paths:
    check(f"Traceability file exists: {path.name}", path.exists(), str(path))

passed = sum(condition for _, condition, _ in checks)
failed = len(checks) - passed

lines = [
    "# Submission Package QC",
    "",
    f"Checks: {len(checks)}",
    f"Passed: {passed}",
    f"Failed: {failed}",
    "",
    "| Check | Status | Detail |",
    "|---|---:|---|",
]
for name, condition, detail in checks:
    status = "PASS" if condition else "FAIL"
    detail = detail.replace("|", "\\|")
    lines.append(f"| {name} | {status} | {detail} |")

(DOCS / "QC_REPORT.md").write_text("\n".join(lines) + "\n", encoding="utf-8")

if failed:
    raise SystemExit(f"{failed} QC checks failed")
