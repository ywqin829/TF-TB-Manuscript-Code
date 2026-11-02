#!/usr/bin/env python
"""Regenerate searchable Markdown from the authoritative DOCX/XLSX outputs."""

from pathlib import Path
import re

from docx import Document
from markitdown import MarkItDown

PROJECT = Path(__file__).resolve().parents[1]
DOCS = PROJECT / "docs"
MD = DOCS / "markdown"
SUPPLEMENT = PROJECT / "supplement"


def convert(source: Path, destination: Path) -> None:
    converter = MarkItDown()
    result = converter.convert_local(source)
    destination.write_text(result.markdown, encoding="utf-8")


def replace_image_placeholders(
    markdown_path: Path,
    replacements: list[tuple[str, str]],
) -> None:
    text = markdown_path.read_text(encoding="utf-8")
    pattern = re.compile(r"!\[[^\]]*\]\(data:image/[^)]*\)")
    matches = pattern.findall(text)
    if len(matches) != len(replacements):
        raise ValueError(
            f"{markdown_path.name}: expected {len(replacements)} image placeholders, "
            f"found {len(matches)}"
        )
    iterator = iter(replacements)
    text = pattern.sub(lambda _: next(iterator), text)
    markdown_path.write_text(text, encoding="utf-8")


def build_abstract() -> None:
    document = Document(DOCS / "Manuscript_Inflammation_Research.docx")
    title = ""
    abstract = ""
    keywords = ""
    in_abstract = False
    for paragraph in document.paragraphs:
        text = paragraph.text.strip()
        if paragraph.style.name == "Title" and not title:
            title = text
            continue
        if text == "Abstract":
            in_abstract = True
            continue
        if in_abstract and text.startswith("1. Introduction"):
            break
        if in_abstract and text.startswith("Keywords:"):
            keywords = text
            continue
        if in_abstract and text and not abstract and text != "Abstract":
            abstract = text

    content = f"# {title}\n\n## Abstract\n\n{abstract}"
    if keywords:
        content += f"\n\n**{keywords}**"
    (MD / "Abstract.md").write_text(content + "\n", encoding="utf-8")


def build_index() -> None:
    root = PROJECT.as_posix()
    content = f"""# Markdown Documents

The DOCX files in `docs/` are the authoritative submission artifacts.
The analysis workbook is an internal support file and is excluded from the upload set.
These Markdown files preserve headings, prose, citations,
and table text for indexing and reuse.

## Submission documents

- [Manuscript_Inflammation_Research.md](<{root}/docs/markdown/Manuscript_Inflammation_Research.md>)  
  Full revised manuscript, including title, abstract, body, figure legends,
  and references.
- [Table1_Cohort_Characteristics.md](<{root}/docs/Table1_Cohort_Characteristics.md>)  
  Verified cohort sources, time points, treatment, HIV status and outcome definitions.
- [Abstract.md](<{root}/docs/markdown/Abstract.md>)  
  Standalone abstract and keywords.
- [Cover_Letter.md](<{root}/docs/markdown/Cover_Letter.md>)  
  Markdown version of the cover letter.
- [Title_Page.md](<{root}/docs/markdown/Title_Page.md>)  
  Title page and affiliations.
- [Supplementary_Information.md](<{root}/docs/markdown/Supplementary_Information.md>)  
  Supplementary figure descriptions.

## Analysis and delivery documents

- [Supplementary_Tables.md](<{root}/docs/markdown/Supplementary_Tables.md>)  
  Internal searchable analysis workbook (51 result sheets; not a submission supplement).
- [InVivo_SingleCell_Dataset_Assessment.md](<{root}/docs/markdown/InVivo_SingleCell_Dataset_Assessment.md>)  
  Fitness-for-purpose assessment of the available human TB single-cell datasets.
- [Expert_Discussion_Publication_Assessment.md](<{root}/docs/Expert_Discussion_Publication_Assessment.md>)  
  Integrated scientific, reviewer, journal-fit and publication-strategy report.
- [Complete_Results_Report.md](<{root}/docs/Complete_Results_Report.md>)  
  Full Chinese result compendium covering all analysis domains and supplementary tables.
- [External_Review_Response_20261029.md](<{root}/docs/External_Review_Response_20261029.md>)  
  Point-by-point response to two rounds of technical and interpretive review.
- [CITATION_ATTRIBUTION_AUDIT.md](<{root}/docs/CITATION_ATTRIBUTION_AUDIT.md>)  
  Before/after record of citation groups split into claim-specific attribution.
- [ANALYSIS_SUMMARY.md](<{root}/docs/ANALYSIS_SUMMARY.md>)  
  Current findings and scientific interpretation.
- [PRE_SUBMISSION_REVIEW.md](<{root}/docs/PRE_SUBMISSION_REVIEW.md>)  
  Editorial and integrity assessment.
- [SUBMISSION_GUIDE.md](<{root}/docs/SUBMISSION_GUIDE.md>)  
  Upload mapping and submission instructions.
- [QC_REPORT.md](<{root}/docs/QC_REPORT.md>)  
  Automated document, citation, and figure checks.
- [PRE_SUBMISSION_GATE.md](<{root}/docs/PRE_SUBMISSION_GATE.md>)  
  Required internal-note and upload-file gate.
"""
    (MD / "INDEX.md").write_text(content, encoding="utf-8")


def main() -> None:
    MD.mkdir(parents=True, exist_ok=True)
    convert(
        DOCS / "Manuscript_Inflammation_Research.docx",
        MD / "Manuscript_Inflammation_Research.md",
    )
    replace_image_placeholders(
        MD / "Manuscript_Inflammation_Research.md",
        [
            ("[Figure 1](<../../figures/Figure1_StudyDesign.pdf>)"),
            ("[Figure 2](<../../figures/Figure2_Discovery.pdf>)"),
            ("[Figure 3](<../../figures/Figure3_CompositionConfounding.pdf>)"),
            ("[Figure 4](<../../figures/Figure4_SingleCellContext.pdf>)"),
            ("[Figure 5](<../../figures/Figure5_OutcomeModels.pdf>)"),
            ("[Figure 6](<../../figures/Figure6_TBRI_OutcomeAssociation.pdf>)"),
        ],
    )
    convert(DOCS / "Cover_Letter.docx", MD / "Cover_Letter.md")
    convert(DOCS / "Title_Page.docx", MD / "Title_Page.md")
    convert(
        DOCS / "Supplementary_Information.docx",
        MD / "Supplementary_Information.md",
    )
    replace_image_placeholders(
        MD / "Supplementary_Information.md",
        [
            ("[Figure S1](<../../figures/FigureS1_METANet.pdf>)"),
            ("[Figure S2](<../../figures/FigureS2_CellTypeSpecificMtbResponse.pdf>)"),
            ("[Figure S3](<../../figures/FigureS3_TBRI_Bcell_TimeAdjusted.pdf>)"),
            ("[Figure S4](<../../figures/FigureS4_GSE31348_CIBERSORT_LM22_Bcell_Trajectory.pdf>)"),
            ("[Figure S5](<../../figures/FigureS5_InVivo_Candidate_TF_Context.pdf>)"),
        ],
    )
    convert(
        SUPPLEMENT / "Supplementary_Tables.xlsx",
        MD / "Supplementary_Tables.md",
    )
    build_abstract()
    build_index()


if __name__ == "__main__":
    main()
