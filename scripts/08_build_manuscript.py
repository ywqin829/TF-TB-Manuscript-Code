#!/usr/bin/env python3
"""Build the revised manuscript, title page, cover letter, and supplement."""

from __future__ import annotations

import csv
import json
import re
import atexit
import tempfile
from pathlib import Path

import docx
import fitz
from docx import Document
from docx.enum.section import WD_SECTION
from docx.enum.table import WD_TABLE_ALIGNMENT
from docx.enum.text import WD_ALIGN_PARAGRAPH, WD_BREAK
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Cm, Inches, Pt


PROJECT = Path(__file__).resolve().parents[1]
FIGS = PROJECT / "figures"
SUPP = PROJECT / "supplement"
DOCS = PROJECT / "docs"
RESULTS = PROJECT / "results"
ORIGINAL = Path(
    r"C:\Users\yw_qi\OneDrive\TF-Manuscript-inflammation research(2)(1).docx"
)

_figure_cache = tempfile.TemporaryDirectory(prefix="tf_tb_docx_figures_")
atexit.register(_figure_cache.cleanup)

# Main-manuscript DOCX files are text-only submission layouts. Figure files are
# uploaded separately, while their legends are collected and placed after the
# reference list.
FIGURE_LEGENDS: list[list[tuple[str, bool, bool]]] = []


def figure_embed(filename: str) -> Path:
    """Rasterize a PDF only in a temporary directory for DOCX embedding."""
    final_stems = {
        "Figure3_SingleCellContext": "Figure4_SingleCellContext",
        "Figure4_OutcomeModels": "Figure5_OutcomeModels",
        "Figure5_TBRI_OutcomeAssociation": "Figure6_TBRI_OutcomeAssociation",
        "Figure6_CompositionConfounding": "Figure3_CompositionConfounding",
        "Figure7_InVivo_Candidate_TF_Context": "FigureS5_InVivo_Candidate_TF_Context",
    }
    stem = Path(filename).stem
    pdf = FIGS / (final_stems.get(stem, stem) + ".pdf")
    png = Path(_figure_cache.name) / Path(filename).with_suffix(".png")
    if not png.exists():
        with fitz.open(pdf) as document:
            document[0].get_pixmap(dpi=240, alpha=False).save(png)
    return png

TITLE = (
    "Immune cell composition shapes the interpretation of blood transcriptional "
    "recovery during tuberculosis treatment"
)
SHORT_TITLE = "Composition-aware TB recovery"

AUTHORS = [
    "Jianhao Xu",
    "Cuiqi Yang",
    "Yang Lu",
    "Xinran Ding",
    "Jiahui Tang",
    "Liangqiong Chen",
    "Feifan Xue",
    "Yongwei Qin",
]
AFFILIATIONS = [
    "Department of Laboratory Medicine, Affiliated Jianhu Hospital of Xinglin College, "
    "Nantong University, Jianhu 224700, China",
    "Department of Pathogen Biology, Medical College, Nantong University, "
    "No. 19 Qixiu Road, Nantong 226001, China",
    "The Intensive Care Unit, Affiliated Hospital of Nantong University, "
    "Nantong 226001, China",
    "Department of Medical Laboratory, Haian People's Hospital, Affiliated Haian Hospital "
    "of Nantong University, Haian 226600, China",
    "Department of Clinical Laboratory, Affiliated Nantong Hospital of Shanghai University, "
    "Nantong, China",
]


def set_cell_shading(cell, fill: str) -> None:
    tc_pr = cell._tc.get_or_add_tcPr()
    shd = OxmlElement("w:shd")
    shd.set(qn("w:fill"), fill)
    tc_pr.append(shd)


def set_page_number(paragraph) -> None:
    paragraph.alignment = WD_ALIGN_PARAGRAPH.CENTER
    run = paragraph.add_run()
    begin = OxmlElement("w:fldChar")
    begin.set(qn("w:fldCharType"), "begin")
    instr = OxmlElement("w:instrText")
    instr.set(qn("xml:space"), "preserve")
    instr.text = " PAGE "
    end = OxmlElement("w:fldChar")
    end.set(qn("w:fldCharType"), "end")
    run._r.extend([begin, instr, end])


def configure_document(document: Document, double_space: bool = True) -> None:
    section = document.sections[0]
    section.top_margin = Cm(2.0)
    section.bottom_margin = Cm(2.0)
    section.left_margin = Cm(2.2)
    section.right_margin = Cm(2.2)

    styles = document.styles
    normal = styles["Normal"]
    normal.font.name = "Times New Roman"
    normal._element.rPr.rFonts.set(qn("w:eastAsia"), "宋体")
    normal.font.size = Pt(11)

    for style_name in ["Normal", "Heading 1", "Heading 2", "Title"]:
        style = styles[style_name]
        style.font.name = "Times New Roman"
        style._element.rPr.rFonts.set(qn("w:eastAsia"), "宋体")
        if double_space:
            style.paragraph_format.line_spacing = 2.0
        style.paragraph_format.space_after = Pt(0)

    styles["Heading 1"].font.size = Pt(13)
    styles["Heading 1"].font.bold = True
    styles["Heading 2"].font.size = Pt(11.5)
    styles["Heading 2"].font.bold = True
    styles["Title"].font.size = Pt(16)
    styles["Title"].font.bold = True

    footer = section.footer.paragraphs[0]
    set_page_number(footer)


def add_runs(paragraph, parts: list[tuple[str, bool, bool]]) -> None:
    for text, bold, italic in parts:
        run = paragraph.add_run(text)
        run.bold = bold
        run.italic = italic


def add_body(document: Document, text: str) -> None:
    paragraph = document.add_paragraph(text)
    paragraph.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY


def add_heading(document: Document, text: str, level: int = 1) -> None:
    document.add_heading(text, level=level)


def add_figure(document: Document, filename: str, legend_parts: list[tuple[str, bool, bool]]) -> None:
    FIGURE_LEGENDS.append(legend_parts)
    # Keep two empty paragraphs at the original call site so the deterministic
    # editorial restructuring remains stable, but do not embed the figure.
    document.add_paragraph()
    document.add_paragraph()


def format_table(table) -> None:
    table.style = "Table Grid"
    table.alignment = WD_TABLE_ALIGNMENT.CENTER
    for row_index, row in enumerate(table.rows):
        for cell in row.cells:
            cell.vertical_alignment = 1
            for paragraph in cell.paragraphs:
                paragraph.paragraph_format.line_spacing = 1.0
                paragraph.paragraph_format.space_after = Pt(2)
                for run in paragraph.runs:
                    run.font.name = "Times New Roman"
                    run.font.size = Pt(9)
            if row_index == 0:
                set_cell_shading(cell, "D9EAF7")
                for paragraph in cell.paragraphs:
                    for run in paragraph.runs:
                        run.bold = True


def add_cohort_characteristics_table(document: Document) -> None:
    path = RESULTS / "cohort_characteristics.tsv"
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    headers = [
        "Dataset",
        "Analysis role",
        "Platform; country",
        "Participants / samples",
        "Sampling and treatment",
        "HIV and susceptibility",
        "Outcome definition",
    ]
    caption = document.add_paragraph()
    caption.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY
    add_runs(
        caption,
        [
            ("Table 1. ", True, False),
            (
                "Cohorts and analytical roles. GSE40553 is the Bloom et al. 2012 "
                "treatment-response study [39]; GSE31348 is the Cliff et al. 2013 "
                "study [8]; GSE89403 is the Thompson et al. 2017 outcome study [9]. "
                "All three source cohorts consisted of HIV-negative participants. "
                "Counts describe the samples retained in the current analysis.",
                False,
                False,
            ),
        ],
    )
    table = document.add_table(rows=1, cols=len(headers))
    for index, header in enumerate(headers):
        table.rows[0].cells[index].text = header

    for row in rows:
        values = [
            f"{row['Dataset']} ({row['Source']}; PMID {row['PMID']})",
            row["Role"],
            f"{row['Platform']}; {row['Country']}",
            f"{row['Current_participants']} / {row['Current_samples']}",
            f"{row['Timepoint_sample_counts']}; {row['Treatment']}",
            f"HIV: {row['HIV_status']}. Susceptibility: {row['Drug_susceptibility']}.",
            row["Outcome_definition"],
        ]
        cells = table.add_row().cells
        for index, value in enumerate(values):
            cells[index].text = value

    format_table(table)
    for row in table.rows:
        for cell in row.cells:
            for paragraph in cell.paragraphs:
                for run in paragraph.runs:
                    run.font.size = Pt(7.5)


def parse_original_references() -> dict[int, str]:
    document = docx.Document(ORIGINAL)
    references = {}
    for paragraph in document.paragraphs:
        match = re.match(r"^(\d+)\.\s+(.*)$", paragraph.text.strip())
        if match:
            references[int(match.group(1))] = match.group(2)
    return references


ORIGINAL_REFS = parse_original_references()

NEW_REFS = {
    19: (
        "Szabo SJ, Kim ST, Costa GL, Zhang X, Fathman CG, Glimcher LH. A novel transcription "
        "factor, T-bet, directs Th1 lineage commitment. Cell. 2000;100(6):655-669. "
        "https://doi.org/10.1016/S0092-8674(00)80702-3"
    ),
    20: (
        "Lazarevic V, Glimcher LH, Lord GM. T-bet: a bridge between innate and adaptive "
        "immunity. Nat Rev Immunol. 2013;13(11):777-789. https://doi.org/10.1038/nri3536"
    ),
    37: (
        "Jung WJ, Acharya S, Ruskin DP, Liao S, Akbary Moghaddam V, Erdenebaatar Z, "
        "Brent MR. METANet: a supervised ensemble learning framework for reconstructing "
        "direct and functional tissue-specific transcription factor networks. "
        "Bioinformatics. 2026;42(10):btag738. https://doi.org/10.1093/bioinformatics/btag738"
    ),
    38: (
        "Wang L, Ma H, Wen Z, Niu L, Chen X, Liu H, et al. Single-cell RNA-sequencing "
        "reveals heterogeneity and intercellular crosstalk in human tuberculosis lung. "
        "J Infect. 2023;87(5):373-384. https://doi.org/10.1016/j.jinf.2023.09.004"
    ),
    39: (
        "Garcia-Alonso L, Holland CH, Ibrahim MM, Turei D, Saez-Rodriguez J. Benchmark "
        "and integration of resources for the estimation of human transcription factor "
        "activities. Genome Res. 2019;29(8):1363-1375. "
        "https://doi.org/10.1101/gr.240663.118"
    ),
    40: (
        "Badia-i-Mompel P, Velez Santiago J, Braunger J, Geiss C, Dimitrov D, Muller-Dott S, "
        "et al. decoupleR: ensemble of computational methods to infer biological activities "
        "from omics data. Bioinform Adv. 2022;2(1):vbac016. "
        "https://doi.org/10.1093/bioadv/vbac016"
    ),
    41: (
        "Hanzelmann S, Castelo R, Guinney J. GSVA: gene set variation analysis for "
        "microarray and RNA-Seq data. BMC Bioinformatics. 2013;14:7. "
        "https://doi.org/10.1186/1471-2105-14-7"
    ),
    42: (
        "Bates D, Machler M, Bolker B, Walker S. Fitting linear mixed-effects models using "
        "lme4. J Stat Softw. 2015;67(1):1-48. https://doi.org/10.18637/jss.v067.i01"
    ),
    43: (
        "Robin X, Turck N, Hainard A, Tiberti N, Lisacek F, Sanchez JC, Muller M. pROC: "
        "an open-source package for R and S+ to analyze and compare ROC curves. "
        "BMC Bioinformatics. 2011;12:77. https://doi.org/10.1186/1471-2105-12-77"
    ),
    44: (
        "Vickers AJ, Elkin EB. Decision curve analysis: a novel method for evaluating "
        "prediction models. Med Decis Making. 2006;26(6):565-574. "
        "https://doi.org/10.1177/0272989X06295361"
    ),
    45: (
        "Finotello F, Mayer C, Plattner C, Laschober G, Rieder D, Hackl H, et al. Molecular "
        "and pharmacological modulators of the tumor immune contexture revealed by "
        "deconvolution of RNA-seq data. Genome Med. 2019;11:34. "
        "https://doi.org/10.1186/s13073-019-0638-6"
    ),
    46: (
        "Sturm G, Finotello F, Petitprez F, Zhang JD, Baumbach J, Fridman WH, et al. "
        "Comprehensive evaluation of transcriptome-based cell-type quantification methods "
        "for immuno-oncology. Bioinformatics. 2019;35(14):i436-i445. "
        "https://doi.org/10.1093/bioinformatics/btz363"
    ),
    47: (
        "Oelen R, de Vries DH, Brugge H, Gordon MG, Vochteloo M, Ye CJ, et al. "
        "Single-cell RNA-sequencing of peripheral blood mononuclear cells reveals "
        "widespread, context-specific gene expression regulation upon pathogenic exposure. "
        "Nat Commun. 2022;13:2893. https://doi.org/10.1038/s41467-022-30893-5"
    ),
    48: (
        "Branchett WJ, Kim JW, Shields J, Chakravarty P, Lee J, Novsarka I, et al. "
        "Airway immune signatures of protection and disease progression in recent human "
        "tuberculosis household contacts. Nat Immunol. 2026;27(8):1577-1589. "
        "https://doi.org/10.1038/s41590-026-02544-0"
    ),
    49: (
        "Kearns K, Tippalagama R, Chawla A, Greenbaum J, De Silva AD, Gunasinghe W, et al. "
        "Enrichment of CD4(-)CD8(-) cytotoxic Vdelta1/3 T cells in pulmonary tuberculosis "
        "disease. iScience. 2026;29(9):117241. https://doi.org/10.1016/j.isci.2026.117241"
    ),
    50: (
        "Raju B, Hoshino Y, Belitskaya-Levy I, Dawson R, Ress S, Gold JA, et al. Gene "
        "expression profiles of bronchoalveolar cells in pulmonary TB. Tuberculosis (Edinb). "
        "2008;88(1):39-51. https://doi.org/10.1016/j.tube.2007.07.003"
    ),
    51: (
        "Joosten SA, van Meijgaarden KE, Del Nonno F, Baiocchini A, Petrone L, Vanini V, et al. "
        "Patients with tuberculosis have a dysfunctional circulating B-cell compartment, which "
        "normalizes following successful treatment. PLoS Pathog. 2016;12(6):e1005687. "
        "https://doi.org/10.1371/journal.ppat.1005687"
    ),
    52: (
        "Ritchie ME, Phipson B, Wu D, Hu Y, Law CW, Shi W, Smyth GK. limma powers differential "
        "expression analyses for RNA-sequencing and microarray studies. Nucleic Acids Res. "
        "2015;43(7):e47. https://doi.org/10.1093/nar/gkv007"
    ),
    53: (
        "Bloom CI, Graham CM, Berry MP, Wilkinson KA, Oni T, Rozakeas F, et al. Detectable "
        "changes in the blood transcriptome are present after two weeks of antituberculosis "
        "therapy. PLoS One. 2012;7(10):e46191. https://doi.org/10.1371/journal.pone.0046191"
    ),
}

ALL_REFERENCES = {**ORIGINAL_REFS, **NEW_REFS}
CITED_REFERENCE_NUMBERS = [1, 2, 3, 4, 5, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 24, 25, 31, 33, 34, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 51, 52, 53]
REFERENCE_RENUMBER = {old: new for new, old in enumerate(CITED_REFERENCE_NUMBERS, start=1)}
REFERENCES = {
    REFERENCE_RENUMBER[old]: text
    for old, text in ALL_REFERENCES.items()
    if old in REFERENCE_RENUMBER
}


def reference_text(number: int) -> str:
    return REFERENCES[number].replace("\t", " ").strip()


def set_submission_metadata(document: Document) -> None:
    """Set human-facing DOCX metadata and remove generator defaults."""
    props = document.core_properties
    props.author = "yongwei Qin"
    props.last_modified_by = "yongwei Qin"
    props.comments = ""
    props.subject = ""


def build_title_page() -> None:
    document = Document()
    set_submission_metadata(document)
    configure_document(document, double_space=True)
    title = document.add_paragraph(style="Title")
    title.alignment = WD_ALIGN_PARAGRAPH.CENTER
    title.add_run(TITLE)

    subtitle = document.add_paragraph()
    subtitle.alignment = WD_ALIGN_PARAGRAPH.CENTER
    run = subtitle.add_run(f"Running title: {SHORT_TITLE}")
    run.italic = True
    run.font.size = Pt(10)

    document.add_paragraph()
    author_p = document.add_paragraph()
    author_p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    author_p.add_run(", ".join(AUTHORS)).bold = True

    for index, affiliation in enumerate(AFFILIATIONS, start=1):
        p = document.add_paragraph()
        p.paragraph_format.line_spacing = 1.15
        p.add_run(f"{index}. {affiliation}").font.size = Pt(9)

    document.add_paragraph()
    corr = document.add_paragraph()
    corr.add_run("Corresponding author\n").bold = True
    corr.add_run(
        "Yongwei Qin, Department of Pathogen Biology, Medical College, Nantong University, "
        "No. 19 Qixiu Road, Nantong 226001, China\n"
        "Email: yw_qin@hotmail.com"
    )

    metadata = document.add_table(rows=5, cols=2)
    format_table(metadata)
    rows = [
        ("Article type", "Original Research"),
        ("Target journal", "Inflammation Research"),
        ("Main figures", "7"),
        ("Supplementary figures", "4"),
        ("Supplementary tables", "None"),
    ]
    for row, (label, value) in zip(metadata.rows, rows):
        row.cells[0].text = label
        row.cells[1].text = value

    document.save(DOCS / "Title_Page.docx")


def build_cover_letter() -> None:
    document = Document()
    set_submission_metadata(document)
    configure_document(document, double_space=False)
    date = document.add_paragraph("31 October 2026")
    document.add_paragraph("Editor-in-Chief")
    document.add_paragraph("Inflammation Research")
    document.add_paragraph("")
    p = document.add_paragraph()
    p.add_run("Re: ").bold = True
    p.add_run(TITLE)

    salutation = document.add_paragraph("Dear Editor,")
    salutation.paragraph_format.space_after = Pt(8)

    paragraphs = [
        (
            "We submit our original research article for consideration in Inflammation Research. "
            "The study asks how much of the peripheral transcriptional recovery observed during "
            "tuberculosis treatment reflects changing leukocyte composition rather than an "
            "intrinsic transcription-factor program."
        ),
        (
            "The central advance is not the identification of a new master regulator. Using "
            "patient-blocked models, we show that treatment is reproducibly accompanied by "
            "B-cell-associated transcriptional recovery, whereas a systematic screen of 271 "
            "TFs found no regulator that retained all cross-cohort criteria after adjustment "
            "for estimated immune-cell composition. This reframes an apparent TBX21/PAX5 "
            "regulatory signature as a largely composition-linked recovery program."
        ),
        (
            "The analysis uses three separate longitudinal datasets with distinct roles: GSE40553 "
            "for discovery, GSE31348 for treatment-program review and score development, and "
            "GSE89403 for outcome association. It includes a South-Africa-only sensitivity "
            "analysis, matched hg19 TPM reconstruction, and separation of between- and "
            "within-participant correlations. The retrospectively selected four-gene TB "
            "Resolution Index is intentionally reported as a B-cell-associated proxy without incremental outcome value: a marker-only "
            "B-cell score that excludes CD79A and EBF1 remains strongly correlated with it."
        ),
        (
            "The manuscript is relevant to inflammation research because it distinguishes "
            "inflammatory contraction and immune reconstitution from intracellular regulatory "
            "activity measured indirectly from whole blood. METANet is used only as a candidate "
            "edge prior in the supplement, and patient-derived BAL single-cell data provide "
            "cross-sectional disease-state context without being treated as a treatment trajectory."
        ),
        (
            "The manuscript is not under consideration elsewhere. The authors declare no "
            "competing interests. The source cohorts are public and de-identified, and no new "
            "human-participant recruitment was performed."
        ),
        "We hope the manuscript will be suitable for peer review.",
    ]
    for text in paragraphs:
        p = document.add_paragraph(text)
        p.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY
        p.paragraph_format.space_after = Pt(8)

    document.add_paragraph("")
    document.add_paragraph("Yours sincerely,")
    document.add_paragraph("Yongwei Qin")
    document.add_paragraph("On behalf of all authors")
    from editorial_structure import revise_cover
    revise_cover(document)
    document.save(DOCS / "Cover_Letter.docx")


def build_main_manuscript() -> None:
    document = Document()
    set_submission_metadata(document)
    configure_document(document, double_space=True)

    title = document.add_paragraph(style="Title")
    title.alignment = WD_ALIGN_PARAGRAPH.CENTER
    title.add_run(TITLE)

    short = document.add_paragraph()
    short.alignment = WD_ALIGN_PARAGRAPH.CENTER
    short.add_run(f"Running title: {SHORT_TITLE}").italic = True

    add_heading(document, "Abstract", 1)
    abstract_sections = [
        (
            "Objective and design: ",
            "To separate transcription-factor recovery from immune-cell composition and "
            "evaluate a cohort-standardized score."
        ),
        (
            "Material or subjects: ",
            "Discovery: GSE40553 (74 participants); development: GSE31348 (27); outcome: "
            "GSE89403 (82); in vivo context: GSE326212 BAL single cells."
        ),
        (
            "Treatment: ",
            "Longitudinal tuberculosis cohorts received standard multidrug therapy; untreated "
            "latent-infection and household-contact groups were comparators."
        ),
        (
            "Methods: ",
            "Donor-blocked limma, DoRothEA/VIPER, quanTIseq, and patient-level mixed models "
            "were used; TBRI was selected retrospectively from 16,644 gene combinations "
            "ranked in GSE31348."
        ),
        (
            "Results: ",
            "We identified 171 restoration-core and 175 disease-driver genes. Among 271 screened "
            "TFs, MAF met all criteria before but not after composition adjustment. "
            "Early outcome AUCs were 0.68-0.73 through Week 4, whereas Week 24 reflected "
            "concurrent recovery state. TBRI "
            "showed a pooled correlation with B-cell fraction (rho = 0.90; between-participant "
            "0.87, within-participant 0.83) and added no significant early discrimination."
        ),
        (
            "Conclusions: ",
            "B-cell-associated transcriptional recovery was reproducible during treatment, but "
            "its outcome association was exploratory and not robust after multiplicity "
            "correction. A substantial component of inferred TF recovery tracks cell composition; no "
            "single composition-independent TF axis or early-monitoring utility is established."
        ),
    ]
    p = document.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY
    for label, text in abstract_sections:
        p.add_run(label).bold = True
        p.add_run(text + " ")

    keywords = document.add_paragraph()
    keywords.add_run("Keywords: ").bold = True
    keywords.add_run(
        "tuberculosis; treatment response; B cells; cell composition; STAT4; "
        "single-cell RNA sequencing; biomarker"
    )

    add_heading(document, "1. Introduction", 1)
    intro = [
        (
            "Tuberculosis remains a major infectious-disease burden, and current treatment "
            "requires prolonged multidrug therapy [2]. Drug-resistant and otherwise difficult-to-treat "
            "forms remain a major clinical challenge [1]. Although most drug-susceptible disease "
            "can be cured, cumulative toxicity and treatment interruption continue to complicate "
            "care [2,3]. Hepatotoxicity is a major safety concern [4], and shortened regimens have "
            "not eliminated the risk of relapse [5]. Host-directed strategies and better "
            "treatment monitoring therefore remain important complements to pathogen-directed "
            "therapy [6]. A useful host marker should reflect the transition from active disease "
            "to durable immune control rather than transient bacterial clearance."
        ),
        (
            "Blood transcriptomic studies have shown that successful treatment produces "
            "coordinated changes across innate and adaptive immune modules. Defective T-cell "
            "exhaustion responses have been reported in pulmonary tuberculosis [7]. Treatment follows "
            "distinct kinetic phases, including early inflammatory contraction and later "
            "humoral re-emergence [8]. Blood RNA signatures can also forecast treatment outcome "
            "[9]. Interferon-inducible and modular blood programs have been described across "
            "independent cohorts [15,16], and host RNA states can identify active disease [13]. "
            "Although these patterns are reproducible, their interpretation has two unresolved problems. "
            "First, whole-blood module recovery may reflect a change in leukocyte composition "
            "rather than coordinated regulation within a fixed cell population. Second, "
            "regulon-level transcription-factor activity is often inferred from the same "
            "transcriptome and can be over-interpreted as a causal upstream mechanism."
        ),
        (
            "B cells contribute to tuberculosis immunity through antibody production, "
            "antigen presentation, cytokine activity, and organization of pulmonary lymphoid "
            "structures [10,11]."
        ),
        (
            "Previous experimental and clinical evidence suggests that B cells and antibody "
            "responses may contribute to host immunity against M. tuberculosis [12,14]."
        ),
        (
            "Circulating B-cell dysfunction has also been reported to normalize after "
            "successful treatment [51]."
        ),
        (
            "At the transcriptional level, PAX5 and EBF1 are core B-cell identity regulators "
            "[20,21]."
        ),
        (
            "TBX21 encodes T-bet, a central regulator of type 1 immunity [18,19]. It is biologically plausible "
            "that recovery of T-cell help and B-cell programs occurs together, but it does not "
            "follow that TBX21 directly drives PAX5 in the same cell or that either TF is "
            "independent of changing cell proportions. Lung parenchyma-homing CD4 T cells also "
            "contribute to Mtb control [22]."
        ),
        (
            "We therefore reanalysed public longitudinal treatment cohorts using the patient "
            "as the unit of inference. The question was whether TF-regulon recovery added "
            "information beyond changing B-cell abundance, not whether B-cell recovery occurs. "
            "We defined a donor-robust reversal program, screened TFs for composition robustness, "
            "evaluated association with treatment failure, and evaluated a compact cohort-standardized score in "
            "an independent outcome-linked cohort. METANet served as a prior, and patient-derived "
            "single-cell data supplied in vivo context."
        ),
    ]
    for paragraph_text in intro:
        add_body(document, paragraph_text)

    add_heading(document, "2. Materials and methods", 1)
    add_heading(document, "2.1 Study design and data sources", 2)
    methods_paragraphs = [
        (
            "This secondary analysis used three longitudinal blood cohorts: GSE40553 for "
            "discovery (74 participants, 176 samples) [53], GSE31348 for treatment-program "
            "review and model development (27, 135) [9], and GSE89403 for outcome association "
            "(82, 312; 74 cured and 8 not cured) [10]. "
            "GSE192483 supplied 11 high- or low-18F-FDG-avidity lung tissues from six patients "
            "[38]. GSE326212 supplied in vivo bronchoalveolar lavage (BAL) single-cell "
            "disease-state data [48]. GSE218065 was a descriptive one-donor stimulation example."
        ),
        (
            "GSE40553 retained 176 of 204 samples with an explicit LTB/TB time label (74 "
            "participants; seven baseline samples from the UK and the remainder from South Africa). "
            "GSE31348 retained all 135 samples from 27 participants. GSE89403 began with 914 GSM "
            "records; 624 eligible GSM records mapped to 312 subject-time IDs from 82 participants "
            "in the processed analysis matrix. Each eligible ID had two GEO records. The matrix "
            "contains one expression column per ID; whether the source preprocessing selected or "
            "combined technical replicates could not be independently reconstructed. GSM mappings, "
            "full flow, and batch counts are provided in the reproducibility archive."
        ),
        (
            "Dataset roles differed by analysis. GSE31348 supported cross-cohort review of "
            "treatment-associated expression programs, TBRI development, and one condition of "
            "the exploratory 271-TF screen. GSE89403 was not used to select TBRI genes or "
            "direction, but it contributed to that TF screen and supplied the outcome-association "
            "analysis; TF findings are therefore screen-derived, not independent validation. "
            "GSE326212 was reserved for cross-sectional in vivo cell-state contextualization "
            "and did not select the blood-based index."
        ),
    ]
    for paragraph_text in methods_paragraphs:
        add_body(document, paragraph_text)
    # Placeholder retains the original paragraph indexing used by the
    # deterministic editorial restructuring; the actual table is appended
    # after the figure legends below.
    document.add_paragraph()
    add_figure(
        document,
        "Figure1_StudyDesign.png",
        [
            ("Fig. 1 ", True, False),
            (
                "Study design and dataset roles. Discovery, model development, external "
                "outcome analysis, cross-sectional single-cell context, and the orthogonal METANet regulatory "
                "prior are separated. Sample and participant numbers are shown for each cohort. "
                "All longitudinal inference uses the patient as the biological unit.",
                False,
                False,
            ),
        ],
    )

    add_heading(document, "2.2 Preprocessing and longitudinal differential expression", 2)
    add_body(
        document,
        "GSE40553 expression values were quantile-normalised, log2-scaled when required, "
        "filtered at the 25th percentile, and collapsed from probes to symbols by selecting "
        "the highest-expression probe per gene. Participant identifiers were parsed from "
        "sample titles. GSE31348 was processed with RMA and probe-to-symbol mapping. GSE89403 "
        "used the deposited log2 expression matrix, with probes or transcripts aggregated "
        "within symbols. Expression matrices and metadata were checked for exact sample order "
        "before all models."
    )
    add_body(
        document,
        "Differential expression used limma [52]. For GSE40553, a six-group design was fitted "
        "with a within-participant consensus correlation estimated by duplicateCorrelation. "
        "The two primary contrasts were active TB versus latent TB (TB_0M versus LTB) and "
        "6 months versus baseline (TB_6M versus TB_0M). Genes were declared significant at "
        "Benjamini-Hochberg FDR < 0.05 and absolute log2 fold change > 1. The restoration-core "
        "set required lower expression in active disease and higher expression after treatment. "
        "The disease-driver set required the opposite direction."
    )

    add_heading(document, "2.3 Module scores, enrichment, and longitudinal summaries", 2)
    add_body(
        document,
        "Single-sample module scores were calculated with ssGSEA using GSVA [29]. GO biological-"
        "process over-representation used clusterProfiler and org.Hs.eg.db with the expressed "
        "transcriptome as background. Module trajectories in GSE40553 and GSE31348 were "
        "summarised with linear mixed models containing participant random intercepts [30]. "
        "For GSE31348, Week 0 was the reference and p values were adjusted within each feature."
    )

    add_heading(document, "2.4 TF-regulon activity and candidate prioritisation", 2)
    add_body(
        document,
        "TF-regulon activity was inferred using DoRothEA confidence levels A-C and the VIPER "
        "method implemented in decoupleR [27,28]. To avoid treating repeated samples as "
        "independent observations, the activity matrix was analysed with donor-blocked limma "
        "models in GSE40553 and with participant random intercepts in GSE31348 and GSE89403. "
        "TF activity is an inferred regulon property, not a direct protein-activity measurement."
    )
    add_body(
        document,
        "For prioritisation, DoRothEA A-D target sets were tested for overlap with the "
        "restoration-core genes using the hypergeometric distribution. We report unadjusted "
        "and FDR-adjusted enrichment. We also calculated Spearman correlations between changes "
        "in TF activity and changes in module score among participants sampled at both baseline "
        "and 6 months. These donor-level correlations are internal concordance analyses, not "
        "independent validation."
    )
    add_body(
        document,
        "After the prior-guided analyses, we applied a post hoc uniform screen to all 271 TFs "
        "with inferable regulon activity using four directional criteria: "
        "higher activity after treatment in GSE40553, higher activity at Week 26 in GSE31348, "
        "lower activity in poor outcomes in GSE89403, and positive paired-donor correlation with "
        "the restoration module. Each criterion required FDR < 0.05. The numeric TF effects in "
        "the three cohorts are cohort-specific model coefficients and are not a common effect "
        "scale. Screen-derived candidates "
        "were then refitted with the original design plus three principal components of "
        "quanTIseq-estimated cell fractions. A candidate was considered composition-robust only "
        "if the original direction and FDR < 0.05 were retained in every cohort used for that "
        "screen. Because GSE89403 contributed one screening criterion and was also the outcome "
        "cohort, it does not independently validate the screened TFs."
    )

    add_heading(document, "2.5 Development of the TB Resolution Index", 2)
    add_body(
        document,
        "TBRI was developed retrospectively. An initial LASSO set from GSE40553 contained 19 "
        "genes; all 2-5-gene combinations from this set (16,644 candidate combinations) were "
        "trained in GSE40553 and the top 400 were ranked by treatment-state discrimination in "
        "GSE31348. CD79A, EBF1, FCGR1A and GK were retained as a compact subset, and signs were "
        "assigned to match the treatment-response direction defined in GSE40553. The simple "
        "score z(CD79A) + z(EBF1) - z(FCGR1A) - z(GK) was then analysed as a cohort-standardized "
        "signature. Because GSE31348 contributed to subset selection, apparent treatment-state "
        "discrimination in this cohort is subject to selection-related optimism and is not considered "
        "independent validation. A GSE31348-mean/SD sensitivity "
        "score was applied to GSE89403. No same-cohort GSE31348 AUC is reported as performance; "
        "GSE89403 AUCs measure outcome discrimination."
    )

    add_heading(document, "2.6 External outcome association in GSE89403", 2)
    add_body(
        document,
        "GSE89403 contains longitudinal samples at diagnosis, Week 1, Week 4, and Week 24. "
        "TBRI and TF-regulon values were evaluated with patient-level mixed models of the form "
        "value ~ group x time + (1 | patient). Contrasts between successful treatment and poor "
        "outcome were estimated at each time. We additionally refitted key models with a B-cell "
        "marker module (MS4A1, CD19, CD79A, CD79B) or T/NK marker module (CD3D, CD3E, NKG7, "
        "GNLY) as a compositional sensitivity analysis."
    )
    add_body(
        document,
        "Outcome discrimination was evaluated at each timepoint with pROC and DeLong confidence "
        "intervals [31]. To avoid apparent calibration from in-sample predictions, logistic "
        "recalibration at each timepoint used leave-one-out probabilities. Calibration slope, "
        "calibration intercept, Brier score, and decision-curve net benefit were treated as "
        "exploratory because only eight participants had poor outcomes [32]."
    )

    add_heading(document, "2.7 Cell-composition analysis", 2)
    add_body(
        document,
        "quanTIseq estimated cell fractions [33,34]. Microarrays were explicitly exponentiated; "
        "GSE89403 TPM was rebuilt from official hg19-aligned raw counts and hg19 gene lengths. Logit fractions "
        "were summarized by three PCs. Estimated fractions were verified to be non-negative and "
        "to sum approximately to one and were treated as composition-associated sensitivity "
        "estimates, not absolute cell counts. Original and screen-derived TF models were compared with "
        "and without PCs by likelihood-ratio test, with incremental marginal R2 and adjusted "
        "contrasts. In GSE89403, B-cell fraction alone was compared with B-cell fraction plus "
        "TBRI by DeLong test, likelihood-ratio test, and delta AUC. TBRI-B-cell correlation was "
        "also separated into a between-participant correlation of participant means and a "
        "within-participant repeated-measures correlation controlling for participant identity. "
        "A marker-only B-cell score (MS4A1, CD19 and CD79B, excluding CD79A and EBF1) was used "
        "to check whether the association depended on genes present in TBRI. For gene-level "
        "attribution, donor-blocked limma models for all restoration-core and disease-driver "
        "genes were compared with and without composition PCs; attenuation metrics were retained "
        "as an exploratory check because small cohort-specific effects make ratios unstable."
    )
    add_body(
        document,
        "We further separated each B-cell exposure into its participant mean and within-participant "
        "deviation, then fitted TBRI against both terms, categorical treatment time, and a participant "
        "random intercept. The primary exposure was logit-transformed quanTIseq B-cell fraction; raw "
        "fraction and a three-gene B-cell marker score excluding TBRI genes were sensitivity analyses. "
        "For GSE31348, we also estimated naive plus memory B-cell fractions using the original "
        "LM22/CIBERSORT method path through IOBR (microarray setting, 200 permutations), then fitted "
        "a participant-intercept time model and compared estimates with quanTIseq. These methods "
        "share the same expression profiles and do not constitute independent cell-count validation."
    )

    add_heading(document, "2.8 Single-cell context and orthogonal TF-network prior", 2)
    add_body(
        document,
        "GSE192483 was descriptive. GSE326212 contains BAL from active tuberculosis, "
        "infection-controlled contacts, and progressors. Among 894,159 cells, doublets were "
        "excluded; expression was aggregated by sample and cell type when at least 20 cells were "
        "available. The main set included 52 participants (19 IGRA-negative, 12 IGRA-positive, "
        "21 active TB), one sample each; models adjusted sex and smoking."
    )
    add_body(
        document,
        "EGAS00001005376 (358,360 cells, 120 donors) was pseudobulked by donor, timepoint and "
        "cell type and analyzed with donor random intercepts [35]. It is supplementary acute-"
        "stimulation evidence, not a treatment model (Supplementary Fig. S2)."
    )
    add_body(
        document,
        "METANet whole-blood and lung maps were queried for TBX21, PAX5, NFATC2, and MAF target "
        "edges and for candidate upstream TFs [25]. Module-overlap enrichment used the tissue-specific "
        "maps and was FDR-adjusted. METANet predictions are derived from GTEx expression, "
        "cis-regulatory activity, and aggregated ChIP-seq labels; they are context-informed "
        "candidate edges, not tuberculosis-specific causal evidence."
    )

    add_heading(document, "2.9 Statistical analysis", 2)
    add_body(
        document,
        "Analyses used R 4.5 and are reproducible from the supplied scripts. Two-sided tests "
        "and a type I error rate of 0.05 were used. Patient was the biological unit for all "
        "longitudinal inference. Exact P values, effect estimates, 95% confidence intervals, "
        "and FDR values for the principal findings are reported in the text and figures. HIV status was reported as "
        "negative in all three source cohorts. "
        "Body mass, diabetes and detailed regimen covariates were not harmonised across cohorts."
    )

    add_heading(document, "3. Results", 1)
    add_heading(
        document,
        "3.1 Donor-blocked analysis identifies an inflammatory-to-B-cell transcriptional reversal",
        2,
    )
    add_body(
        document,
        "In GSE40553, active TB differed extensively from latent infection and treatment "
        "reversed many changes (Fig. 2A, B). Donor-blocked models identified 171 "
        "restoration-core and 175 disease-driver genes. The duplicateCorrelation consensus "
        "parameter was 0.075; genome-wide disease and treatment log2 fold changes were strongly "
        "opposite (Pearson r = -0.753, Spearman rho = -0.687 across 21,186 genes). These "
        "contrasts share the baseline group, so this correlation alone cannot establish biological reversal. "
        "With both gene sets fixed before external testing, 162/171 restoration-core genes "
        "increased and 156/163 measurable disease-driver genes decreased from Week 0 to Week 26 "
        "in GSE31348 (expression-matched permutation P < 0.0002 for each set; 27 paired patients). "
        "Twelve disease-driver genes were absent from the external expression matrix. This "
        "cohort has no latent-TB comparator, so it tests direction rather than full normalization. Excluding "
        "seven UK baseline samples, 161/171 restoration-core (94.2%) and "
        "161/175 disease-driver (92.0%) genes were recovered from 28 South African active-TB "
        "baseline and 38 LTB samples. This South-Africa-only restoration module retained B-cell "
        "proliferation and receptor signalling as top GO processes and increased from Week 0 to "
        "Week 26 in GSE31348. In the sensitivity GSVA rerun, its effect retained 92.1% of the "
        "original module effect (P = 7.3 x 10^-30). SSGSEA normalized scores from separate fits "
        "use a run-specific multiplier; the original score in this rerun correlated perfectly "
        "with the primary score but was multiplied by 2.079, so only effects within the same fit "
        "are compared. Candidate TF directions were unchanged in this population subset: "
        "TBX21, PAX5, MAF and NFATC2 all increased after treatment (effects 2.39, 2.17, 1.56 and "
        "1.36; FDR < 0.001). Restoration-core "
        "expression increased during treatment (6-month estimate 0.434), the disease-driver "
        "score decreased (-0.485), and both remained displaced from baseline at 12 months "
        "(Fig. 2C)."
    )
    add_body(
        document,
        "GO enrichment separated the modules biologically. The restoration core was enriched "
        "for B-cell proliferation, B-cell receptor signalling, and B-cell activation. The "
        "disease-driver module was dominated by responses to lipopolysaccharide, bacterial "
        "molecules, and virus. These results reproduce the broad inflammatory-to-adaptive "
        "transition described previously [8], but they also show that the adaptive component "
        "is strongly represented by B-cell processes."
    )
    add_figure(
        document,
        "Figure2_Discovery.png",
        [
            ("Fig. 2 ", True, False),
            (
                "Discovery of a donor-robust treatment-reversal program. A, Active TB versus "
                "latent TB. B, Month 6 versus baseline. C, Module scores across treatment; "
                "points show means and error bars show 95% confidence intervals. D, Changes in "
                "TF-regulon activity from baseline to Month 6. E, Spearman correlations between "
                "paired participant changes in TF activity and the restoration-core score (n = 24 paired participants); "
                "correlation P values are Benjamini-Hochberg FDR-adjusted across the TF correlations. "
                "Panels A and B use the same 21,186-gene universe, with the discovery-defined module labels applied to the treatment contrast in panel B. "
                "F, Target-set enrichment of treatment-responsive TFs in the restoration core. "
                "All expression models use participant-blocked inference.",
                False,
                False,
            ),
        ],
    )

    add_heading(
        document,
        "3.2 Initial TF candidates track recovery, but a systematic screen does not establish a master regulator",
        2,
    )
    add_body(
        document,
        "Regulon activity showed coordinated reversal during treatment. PAX5, TBX21, and NFATC2 "
        "inferred activity scores all increased from baseline to Month 6 (PAX5 model estimate "
        "1.72, FDR 9.4 x 10^-6; TBX21 2.18, FDR 1.7 x 10^-6; NFATC2 1.31, FDR 1.5 x 10^-6). "
        "These are differences in VIPER regulon-activity scores, not log2 expression fold "
        "changes. Interferon and NF-kB regulons including IRF1, STAT2, STAT1, and NFKB1 "
        "decreased (Fig. 2D)."
    )
    add_body(
        document,
        "There were 35 baseline and 25 Month 6 samples; 24 participants had both, leaving 11 "
        "baseline-only and one Month 6-only sample. All 24 paired participants were South African. "
        "Their change in TBX21 "
        "activity showed the strongest correlation with restoration-core change (Spearman "
        "rho = 0.84, FDR = 1.0 x 10^-4). PAX5 showed a weaker correlation (rho = 0.41, nominal "
        "P = 0.044, FDR = 0.18), and NFATC2 was not significant after FDR correction (rho = 0.36, "
        "FDR = 0.24; Fig. 2E)."
    )
    add_body(
        document,
        "Target-set enrichment was less conclusive. PAX5 targets were enriched in the "
        "restoration core (5 of 55 regulon genes; fold enrichment 9.5; FDR = 0.0099). TBX21 "
        "showed nominal enrichment (3 of 76; fold 4.1; unadjusted P = 0.036), but this did not "
        "survive correction across TFs (FDR = 0.24). NFATC2 also had nominal enrichment without "
        "FDR support (FDR = 0.069; Fig. 2F). These analyses identify internal correlates, but "
        "they were initiated from a small prior-guided candidate set and do not by themselves "
        "establish an upstream regulator."
    )
    add_body(
        document,
        "The systematic screen changed this interpretation. MAF was the only one of 271 TFs "
        "meeting all four criteria: higher after treatment in GSE40553 (1.29, FDR = 0.0013), "
        "higher at Week 26 in GSE31348 (0.11, FDR = 0.0010), lower in poor outcomes in GSE89403 "
        "(-0.26, FDR = 0.0055), and positively correlated with donor restoration (rho = 0.67, "
        "FDR = 0.014). TBX21 met three criteria; PAX5 and NFATC2 met two. After composition "
        "adjustment, MAF retained the GSE89403 outcome association (FDR = 0.046), but that was "
        "within one of the screening cohorts and does not constitute independent confirmation; "
        "it also lost treatment and Week-26 support, and no TF met all three cohort criteria. "
        "METANet ranked "
        "MAF among predicted upstream inputs to PAX5, providing a prior rather than validation. "
        "The data therefore support a composition-linked recovery program rather than a single "
        "composition-independent master regulator."
    )

    add_heading(
        document,
        "3.3 Cross-cohort treatment check confirms module direction but is not independent validation",
        2,
    )
    add_body(
        document,
        "In GSE31348, the restoration-core score increased from Week 0 to Week 26 (mixed-model "
        "estimate 0.203, P = 8.2 x 10^-29), while the disease-driver score decreased (-0.331, "
        "P = 2.3 x 10^-27). Directional change was already detectable at Week 1, although the "
        "largest separation occurred at the end of treatment. TBX21 and NFATC2 activity also "
        "increased over the treatment course, and PAX5 activity returned above baseline by "
        "Week 26."
    )
    add_body(
        document,
        "Because the four genes came from combinations ranked in GSE31348, no same-cohort AUC "
        "is reported as treatment-state discrimination performance. Both the earlier value near "
        "0.99 and the simpler score value of 0.757 were generated after selection in that cohort "
        "and are not independent performance estimates. The four genes were manually retained "
        "from the top-ranked combinations as a compact, interpretable set rather than declared "
        "the uniquely best model. GSE89403, which did not select the TBRI genes or signs, is the "
        "only cohort used for outcome association."
    )

    add_heading(document, "3.4 Lesional single-cell data place the candidate regulators in distinct immune compartments", 2)
    add_body(
        document,
        "In 11 lung tissues from six tuberculosis patients (GSE192483), PAX5 was confined "
        "largely to B-lineage cells, while TBX21 and NFATC2 were detected in T and NK "
        "compartments (Fig. 3A,B). The analysis supports lineage "
        "plausibility and the interpretation that the bulk PAX5 signal may contain a B-cell "
        "composition component. Because GSE192483 does not contain an independent healthy-control "
        "arm, no disease-versus-health cell-level claim is made. The one-donor GSE218065 "
        "stimulation library is not used as patient-level evidence."
    )
    add_figure(
        document,
        "Figure3_SingleCellContext.png",
        [
            ("Fig. 3 ", True, False),
            (
                "Candidate-regulator localization in tuberculous lung tissue (GSE192483). "
                "A, Patient-tissue-level median expression and detection of five candidate "
                "TFs across annotated cell classes in 11 tissues from six patients. "
                "B, PAX5, TBX21 and NFATC2 expression on the UMAP. Both "
                "panels are descriptive; no cell-level inferential P values are reported.",
                False,
                False,
            ),
        ],
    )

    add_heading(
        document,
        "3.5 Patient-level models assess outcome association and compositional effects",
        2,
    )
    add_body(
        document,
        "In GSE89403, TBX21 activity was lower in poor outcomes (unadjusted group P = 0.027), "
        "but the group-by-time interaction was not significant (P = 0.38). After adjustment for "
        "the T/NK marker module, the time effect attenuated (P = 0.10) and the group term no "
        "longer survived correction across features (FDR = 0.15; Fig. 4B). PAX5 activity was "
        "also lower in poor outcomes before adjustment (P = 0.020), but the difference "
        "disappeared after adjustment for the B-cell marker module (P = 0.81; Fig. 4A)."
    )
    add_body(
        document,
        "The B-cell marker module showed lower values in poor outcomes with a nominal "
        "group-by-time interaction (P = 0.036), but this did not survive correction across "
        "features (FDR = 0.14; Fig. 4D). The PAX5 target module showed a similar unadjusted interaction "
        "(P = 0.040; FDR = 0.20), which attenuated after B-cell adjustment (P = 0.83; Fig. 4C). The TBX21 "
        "target module did not separate outcome groups (interaction P = 0.41; Fig. 4E). "
        "The corresponding TBRI trajectory is shown in Fig. 4F. B-cell-associated recovery was reproducible during treatment, but its outcome "
        "association remained exploratory."
    )
    add_figure(
        document,
        "Figure4_OutcomeModels.png",
        [
            ("Fig. 4 ", True, False),
            (
                "Patient-level longitudinal models in GSE89403. A-E, Model-estimated trajectories "
                "for PAX5 activity, TBX21 activity, the PAX5 target module, the B-cell marker "
                "module, and the TBX21 target module. F, TBRI trajectory. Points show estimated "
                "marginal means and error bars show 95% confidence intervals from linear mixed "
                "models with participant random intercepts. Blue denotes successful treatment "
                "and red denotes poor outcome. Group-by-time values shown in the panels are nominal "
                "(unadjusted) P values; FDR-adjusted values are reported in the Results and supporting tables.",
                False,
                False,
            ),
        ],
    )

    add_heading(
        document,
        "3.6 External outcome association is weak and does not validate monitoring utility",
        2,
    )
    add_body(
        document,
        "In GSE89403, the cohort-standardized TBRI produced AUCs of 0.680 at baseline (95% CI 0.470-0.890), "
        "0.726 at Week 1 (0.528-0.925), 0.679 at Week 4 (0.430-0.927), and 0.933 at Week 24 "
        "(0.870-0.997; Fig. 5B). Week 1 showed nominal discrimination, but only six poor-outcome "
        "participants had samples at that timepoint; the baseline and Week 4 intervals included 0.5. "
        "These results do not establish robust early predictive utility. The Week 24 estimate reflects "
        "concordance with the treatment state already reached rather than early prediction. A "
        "sensitivity score using GSE31348 means and SDs directly in GSE89403 "
        "produced similar AUCs (0.70, 0.73, 0.69 and 0.94), although cross-platform calibration "
        "is not deployment-grade. Score distributions are shown in Fig. 5A and the Week 1/Week 4 "
        "ROC curves in Fig. 5C."
    )
    add_body(
        document,
        "Exploratory calibration estimates are reported in the derived results archive; the small "
        "poor-outcome group does not support a clinical calibration claim."
    )
    add_figure(
        document,
        "Figure5_TBRI_OutcomeAssociation.png",
        [
            ("Fig. 5 ", True, False),
            (
                "External outcome association of the cohort-standardized TB Resolution Index. "
                "A, Score distribution by treatment outcome and time. B, Outcome discrimination "
                "with 95% DeLong intervals. C, Week 1 and Week 4 ROC curves, with AUCs, 95% DeLong intervals and poor-outcome sample counts.",
                False,
                False,
            ),
        ],
    )

    add_heading(
        document,
        "3.7 Immune-cell composition is associated with a substantial component of inferred TF recovery",
        2,
    )
    add_body(
        document,
        "After explicit non-log/TPM input correction, quanTIseq decomposed all three bulk cohorts into reference immune-cell fractions. "
        "Adding three composition principal components significantly improved every tested "
        "PAX5, TBX21, and NFATC2 model after FDR correction (Fig. 6A, B). The incremental "
        "marginal R2 for TBX21 was 0.17 in GSE40553, 0.20 in GSE31348, and 0.15 in GSE89403. "
        "For PAX5, the corresponding increments were 0.12, 0.11, and 0.11. These effects were "
        "smaller for NFATC2, with increments of 0.04, 0.10, and 0.05."
    )
    add_body(
        document,
        "Preselected TF treatment contrasts in GSE31348 also shrank after composition adjustment: "
        "PAX5 changed from 0.280 (95% CI 0.180 to 0.381) to 0.015 (-0.123 to 0.153), "
        "TBX21 from 0.257 (0.157 to 0.358) to -0.005 (-0.134 to 0.123), and MAF "
        "from 0.111 (0.050 to 0.172) to 0.029 (-0.060 to 0.119). Effects and 95% CIs for "
        "all three cohorts are provided in the reproducibility archive. These changes do not "
        "identify a causal mediation effect. CD79A and CD19 overlap the PAX5 target set; in "
        "GSE89403, an outcome contrast adjusted for the nonoverlapping MS4A1 and CD79B marker "
        "score was -0.042 (95% CI -0.213 to 0.129)."
    )
    add_body(
        document,
        "An exploratory gene-level donor-blocked analysis showed that adding composition PCs "
        "changed the treatment coefficient for many module genes: the composition joint FDR was "
        "below 0.05 for 74.3% of restoration-core genes in GSE40553 and 83.6% in GSE31348; the "
        "median attenuation was 26.1% and 57.1%, respectively. However, per-gene attenuation "
        "magnitudes did not correlate across cohorts (restoration rho = 0.06; disease-driver "
        "rho = -0.03), including after requiring an absolute base effect of at least 0.5 log2 "
        "units. We therefore do not infer a stable gene-level cell-intrinsic fraction; the "
        "composition interpretation remains at the program level."
    )
    add_body(
        document,
        "The external index tracked estimated B-cell composition. Across 312 GSE89403 samples, the "
        "pooled descriptive Spearman correlation between TBRI and B-cell fraction was 0.903. "
        "Among 82 participants, the between-participant rho was 0.871 and the within-participant "
        "repeated-measures rho was 0.832 (P = 1.9 x 10^-60), showing both stable "
        "between-person differences and coordinated within-person change (Fig. 3C). B-cell "
        "fractions were lower in poor outcomes (Fig. 3D). A marker-only B-cell score excluding "
        "CD79A and EBF1 (MS4A1, CD19 and CD79B) remained correlated with TBRI (pooled rho = 0.826; "
        "between-participant rho = 0.857; within-participant rho = 0.849), and only EBF1 among "
        "the four TBRI genes appears in the fixed quanTIseq 153-gene signature. Thus "
        "the association is not solely caused by exact score-gene sharing, although TBRI remains "
        "a B-cell-associated proxy rather than an independent biomarker. B-cell fraction alone "
        "achieved outcome AUCs of 0.694, 0.760, 0.654 and 0.889 at baseline, Week 1, Week 4 and "
        "Week 24. Adding TBRI changed AUC by -0.002, -0.007, +0.031 and +0.038; Week 24 was "
        "borderline by DeLong test (P = 0.047) but not by likelihood-ratio test (P = 0.051), and "
        "early timepoints showed no improvement (DeLong P = 0.948, 0.507 and 0.481; Fig. 3E, F)."
    )
    add_body(
        document,
        "After adjustment for categorical treatment time and participant intercepts, the "
        "within-participant logit B-cell fraction remained associated with TBRI in 312 records "
        "from 82 participants (coefficient 1.816, 95% Wald CI 1.548-2.084; P = 2.42 x 10^-30). "
        "The association also persisted for raw fraction and the three-gene marker score, which "
        "excludes TBRI genes (Supplementary Fig. S3). In GSE31348, LM22/CIBERSORT estimated a "
        "Week 26 versus baseline total B-cell fraction increase of 0.00779 (95% Wald CI "
        "0.00327-0.01231; P = 0.00101) in 135 records from 27 participants; its agreement with "
        "quanTIseq was moderate (Spearman rho = 0.493; Supplementary Fig. S4). In a paired "
        "27-patient comparison, the estimated Week 26 B-cell fraction increased by 0.00779 "
        "(95% CI 0.00140-0.01418) with CIBERSORT and 0.02597 (0.01264-0.03931) with "
        "quanTIseq; estimated neutrophil fraction changes were -0.12320 (95% CI -0.15036 "
        "to -0.09604) and -0.12624 (-0.17142 to -0.08107), respectively. These are relative fractions, not "
        "absolute cell counts. The mixed-model and paired-comparison confidence intervals were obtained from "
        "different inferential procedures and are therefore reported separately. The time-adjusted "
        "association supports coordinated change but cannot establish causation or independent "
        "biomarker value."
    )
    add_body(
        document,
        "This result reframes the biological interpretation. The treatment-associated B-cell "
        "program is reproducible, but the TBRI does not provide clear incremental information "
        "beyond the estimated B-cell composition signal. Likewise, the apparent coordinated "
        "TBX21-PAX5 regulon reversal attenuated after adjustment for estimated immune-cell composition. The "
        "screen-derived MAF candidate retained its outcome association after adjustment but lost "
        "treatment-response support. The remaining "
        "associations support cell-state hypotheses, but no TF can be interpreted as a "
        "demonstrated intracellular master regulator."
    )
    add_figure(
        document,
        "Figure6_CompositionConfounding.png",
        [
            ("Fig. 6 ", True, False),
            (
                "Estimated immune-cell composition influences inferred TF activity and RNA biomarker interpretation. "
                "A, Marginal R2 with and without composition principal components. "
                "B, likelihood-ratio tests for adding composition principal components, shown as Benjamini-Hochberg FDR-adjusted values. C, Correlation between TBRI and quanTIseq-estimated "
                "B-cell fraction in GSE89403. D, Estimated B-cell fraction by treatment outcome. "
                "E, Outcome AUC for B-cell fraction alone versus B-cell fraction plus TBRI. "
                "F, Incremental AUC of TBRI over B-cell fraction. All models retain patient-level "
                "inference where repeated samples are present. Panel C shows the pooled "
                "descriptive correlation; participant-level between- and within-participant "
                "estimates are reported in the Results.",
                False,
                False,
            ),
        ],
    )

    add_heading(
        document,
        "3.8 In vivo airway single-cell data contextualize the candidate TFs",
        2,
    )
    add_body(
        document,
        "GSE326212 provided independent in vivo BAL data from active tuberculosis, "
        "IGRA-negative contacts, IGRA-positive contacts, and progressors. After doublet "
        "exclusion, 893,829 cells were aggregated into patient-sample and coarse-cell-type "
        "pseudobulk units. The resulting cell-type map separated B-lineage PAX5 and MEF2B, "
        "T/NK-associated TBX21 and RUNX3, and broadly expressed MAF and STAT4. These patterns "
        "are descriptive and were not generated by treating cells as independent donors."
    )
    add_body(
        document,
        "In 19 IGRA-negative and 21 active-TB participants, STAT4 pseudobulk expression was "
        "higher in PTB MoAM than IGRA-negative contacts (0.0232, SE 0.0040, FDR = 1.2 x 10^-4) "
        "and IGRA-positive contacts (0.0225, SE 0.0050, FDR = 0.011). It reflected a higher "
        "STAT4-positive fraction (+3.82 points, SE 0.76, P = 8.4 x 10^-6); expression per "
        "positive cell did not differ (0.27, P = 0.12). Raw counts confirmed count>=1 and "
        "mean-count differences (P = 0.018 and 0.028); count>=2 attenuated (0.11). The effect "
        "persisted after project/cell-number adjustment, in both neutrophil strata and all 15 "
        "leave-one-project-out models. MAF was lower in PTB CD4 T cells (FDR = 0.092). This "
        "agrees with prior BAL Th1 evidence [50] and does not support a broad TBX21-PAX5 axis."
    )
    add_body(
        document,
        "This analysis addresses the in vivo context that acute stimulation cannot provide, but "
        "it does not create a longitudinal treatment result. GSE326212 is cross-sectional; age "
        "was unavailable and project was partly confounded with disease group. The acute "
        "EGAS00001005376 stimulation experiment is therefore retained only as supplementary "
        "cell-state evidence, where it is explicitly not used to contradict clinical recovery."
        " (Supplementary Fig. S2)."
    )
    add_figure(
        document,
        "Figure7_InVivo_Candidate_TF_Context.png",
        [
            ("Fig. 7 ", True, False),
            (
                "In vivo airway single-cell context for candidate TFs (GSE326212). "
                "A, Patient-sample-level pseudobulk expression across coarse BAL cell types. "
                "Colour indicates relative expression scaled within each TF, and dot size "
                "indicates the median percentage of expressing cells. B, Sample-level contrasts "
                "between active tuberculosis and IGRA-negative contacts. Colour shows the effect "
                "in standard-error units; asterisks indicate BH-adjusted P values (* < 0.05, "
                "** < 0.01, *** < 0.001). Cells were aggregated by sample and cell "
                "type before all comparisons. This is cross-sectional in vivo evidence, not a "
                "treatment trajectory. Raw-count threshold sensitivity is reported in the Supplementary Information.",
                False,
                False,
            ),
        ],
    )

    add_heading(document, "4. Discussion", 1)
    discussion = [
        (
            "Earlier longitudinal studies established that blood transcriptional profiles change "
            "during effective tuberculosis therapy, with an early decline in inflammatory programs "
            "followed by a slower return of lymphocyte-associated expression [8,40]. The question "
            "raised by those observations is therefore not simply whether the blood transcriptome "
            "responds to treatment, but what biological process a returning signature represents. "
            "Our reanalysis suggests that apparent recovery in whole blood should be interpreted as "
            "an integrated change in the sampled immune system, rather than automatically as "
            "restoration of gene regulation within individual cells. A reproducible expression "
            "change can be informative about disease state while remaining ambiguous about the "
            "cellular process that generated it."
        ),
        (
            "The opposing trajectories of innate and B-cell-associated transcripts provide a useful "
            "example. Active tuberculosis is characterized by prominent interferon-inducible and "
            "myeloid blood programs [16,17], whereas later B-cell-related expression during therapy "
            "was already described by Cliff and colleagues [9]. In the present analysis, a reduction "
            "in estimated neutrophil representation accompanied a smaller rise in relative B-cell "
            "fraction. Because cell fractions are compositional, the apparent increase of one "
            "population can partly reflect contraction of another; neither deconvolution method "
            "measures absolute B-cell numbers. This does not make the B-cell signal biologically "
            "trivial. Independent cellular studies have described abnormalities in circulating B-cell "
            "phenotype and function during active tuberculosis and improvement after treatment [51]. "
            "Rather, our findings show why B-cell-associated RNA recovery and B-cell numerical or "
            "functional reconstitution must be kept distinct."
        ),
        (
            "This distinction also refines the meaning of treatment-associated normalization. An early "
            "fall in inflammatory transcripts need not coincide with slower restoration of adaptive "
            "immune features, and neither necessarily marks complete resolution of tissue disease. "
            "The staged responses described by Bloom and Cliff are consistent with different parts of "
            "the host response evolving on different timescales [9,53]. Thompson and colleagues "
            "further related blood RNA patterns to treatment outcome and pulmonary imaging, including "
            "evidence that inflammatory abnormalities may persist despite microbiological cure [9]. "
            "Taken together, these studies argue against treating a single blood score as a universal "
            "measure of bacterial clearance, immune competence and recovery of lung pathology."
        ),
        (
            "The consequences are especially important for inferred transcription-factor activity. "
            "Regulon methods estimate activity from target-gene expression rather than directly "
            "measuring transcription-factor action within a defined cell type [27,28]. When the "
            "relative abundance of target-expressing populations changes, a bulk regulon score may "
            "shift without a corresponding change in intracellular regulation. In our data, "
            "treatment-associated PAX5, TBX21 and MAF effects were attenuated after accounting for estimated "
            "immune-cell composition. PAX5 and EBF1 are established B-lineage regulators [20,21], "
            "whereas TBX21 is associated with type 1 immune differentiation [18,19]. Attenuation "
            "after adjustment is therefore a sensitivity analysis for composition, not a causal estimate "
            "of how much recovery is mediated by cell abundance. Because composition estimates and "
            "regulon scores are derived from overlapping expression information, adjustment can remove "
            "genuine biological variation as well as confounding."
        ),
        (
            "The same interpretive distinction applies to blood RNA biomarkers. A score can "
            "discriminate treatment states because it summarizes a broad change in immune-cell mixture, "
            "yet contribute little information about future outcome once that mixture is taken into "
            "account. TBRI illustrates this problem: its close association with estimated B-cell "
            "fraction, including after accounting for treatment time, makes it a useful representation "
            "of the recovering blood expression state, but not an independently established measure of "
            "treatment success. The limited number of poor outcomes prevents a definitive assessment of "
            "incremental prognostic value. Previous studies identified blood RNA signatures associated "
            "with treatment outcomes [9,17], whereas signatures developed for active-disease diagnosis "
            "or future disease progression address different clinical questions [23,24]. A clinically "
            "useful treatment-monitoring marker must therefore demonstrate added value over clinical, "
            "microbiological and cell-composition information, not merely track treatment time."
        ),
        (
            "There is a broader implication for studies of inflammation based on mixed-cell specimens. "
            "Longitudinal expression changes can arise from altered leukocyte trafficking, shifts in "
            "population abundance, changes in the activation state of existing cells, or combinations "
            "of these processes. Composition-aware analysis can reveal when an apparent regulatory "
            "signature is sensitive to the cellular mixture, but it should not automatically remove every "
            "immune-cell-associated signal. If redistribution of leukocyte populations is itself part of "
            "the treatment response, indiscriminate adjustment risks subtracting the biology under study."
        ),
        (
            "The available data also constrain the level of inference. The discovery contrasts share a "
            "baseline group, although directional testing in a separate cohort reduces reliance on this "
            "internally correlated comparison. The two deconvolution approaches use the same bulk "
            "transcriptomes and provide algorithmic sensitivity analyses rather than orthogonal cell "
            "counts. GSE89403 includes only eight poor outcomes, and the lung and airway single-cell "
            "datasets provide cross-sectional localization rather than treatment-linked evidence of "
            "intracellular recovery. The cohorts were HIV-negative, and clinical covariates were not "
            "consistently available. These limitations restrict claims about mechanism and clinical "
            "transportability without negating the reproducible treatment-associated expression pattern."
        ),
        (
            "The next step is a more discriminating validation strategy rather than another search for "
            "a universal recovery gene or master regulator. Paired absolute leukocyte counts and flow-"
            "cytometric B-cell measurements, combined with expression profiling of purified populations "
            "or donor-linked longitudinal single-cell data, would separate population turnover from "
            "within-cell transcriptional change. Larger prospective cohorts with serial microbiological "
            "measurements and, when feasible, imaging could then establish whether composition-sensitive "
            "and composition-independent RNA features provide complementary clinical information."
        ),
        (
            "In conclusion, tuberculosis treatment-associated blood transcriptional recovery is "
            "reproducible, but its mechanistic and clinical meaning depends on distinguishing changes "
            "in the cells present from changes occurring within those cells. The present analysis does "
            "not establish absolute B-cell reconstitution, a composition-independent master regulator, "
            "or independent early prognostic value for TBRI. It does provide a composition-aware "
            "framework for interpreting treatment-associated blood RNA changes and for designing the "
            "cell-resolved and prospective studies needed to test their biological and clinical value."
        ),
    ]
    for paragraph_text in discussion:
        add_body(document, paragraph_text)

    add_heading(document, "Funding", 1)
    add_body(
        document,
        "This work was supported by the China Postdoctoral Science Foundation (2020M671561) "
        "and the Jiangsu Postdoctoral Science Foundation (2020Z130)."
    )

    add_heading(document, "Author contributions", 1)
    add_body(
        document,
        "Jianhao Xu: Conceptualization, Methodology, Formal analysis, Data curation, Writing "
        "- original draft. Cuiqi Yang: Data curation, Validation, Formal analysis. Yang Lu: "
        "Formal analysis, Data curation. Xinran Ding: Investigation, Data curation. Jiahui Tang: "
        "Investigation. Liangqiong Chen: Investigation. Feifan Xue: Investigation, Validation. "
        "Yongwei Qin: Conceptualization, Supervision, Resources, Writing - review and editing."
    )

    add_heading(document, "Ethics approval", 1)
    add_body(
        document,
        "This study used publicly available, de-identified datasets. No new human participants "
        "or identifiable private information were collected. Ethical approval was granted by "
        "the original studies."
    )

    add_heading(document, "Consent to participate", 1)
    add_body(
        document,
        "Not applicable. The present study did not recruit participants."
    )

    add_heading(document, "Data availability", 1)
    add_body(
        document,
        "Source transcriptomic data are available from NCBI GEO under GSE40553, GSE31348, "
        "GSE89403, GSE192483, GSE218065, and GSE326212. The acute-stimulation EGA resource is "
        "EGAS00001005376. The reanalysis scripts, processed result tables, and figure-generation "
        "code and derived result tables accompany this submission as a reproducibility archive. "
        "No new raw sequencing data were generated."
    )

    add_heading(document, "Code availability", 1)
    add_body(
        document,
        "The R and Python scripts used to reproduce the reported analyses and figures are "
        "included in the accompanying Supplementary_Code.zip archive and prepared for public "
        "repository deposition. The workflow uses fixed random seeds and records "
        "software versions in the supplied log files."
    )

    # Preserve paragraph positions used by the editorial restructuring without
    # leaving an internal AI declaration or implementation detail in the DOCX.
    document.add_paragraph()
    document.add_paragraph()

    add_heading(document, "Competing interests", 1)
    add_body(document, "The authors declare no competing interests.")

    add_heading(document, "References", 1)
    for number in sorted(REFERENCES):
        paragraph = document.add_paragraph()
        paragraph.paragraph_format.line_spacing = 1.15
        paragraph.paragraph_format.left_indent = Cm(0.5)
        paragraph.paragraph_format.first_line_indent = Cm(-0.5)
        paragraph.add_run(f"{number}. ").bold = True
        paragraph.add_run(reference_text(number))

    add_heading(document, "Figure legends", 1)
    # Keep the submission legend block in figure-number order even when the
    # editorial restructuring moves a figure's Results section earlier or later.
    def legend_sort_key(legend_parts: list[tuple[str, bool, bool]]) -> tuple[int, int]:
        import re
        text = "".join(part[0] for part in legend_parts)
        match = re.search(r"Fig\.\s+(S?)(\d+)", text)
        if not match:
            return (2, 999)
        # editorial_structure.revise_main remaps the original figure calls
        # (1, 2, 3, 4, 5, 6, 7) to the final submission labels
        # (1, 2, 4, 5, 6, 3, S5). Sort by the final labels before that remap.
        original_number = int(match.group(2))
        final_order = {1: 0, 2: 1, 6: 2, 3: 3, 4: 4, 5: 5, 7: 6}
        return (0, final_order.get(original_number, 999))

    for legend_parts in sorted(FIGURE_LEGENDS, key=legend_sort_key):
        # Figure 7 was remapped to supplementary Figure S5. Keep it with the
        # other supplementary legends below rather than mixing it into the
        # primary Figure 1-6 legend block.
        legend_text = "".join(part[0] for part in legend_parts)
        # The source legend is still labelled Figure 7 at this stage; the
        # editorial pass remaps that label to Supplementary Figure S5 later.
        if "Fig. S5" in legend_text or "Figure 7" in legend_text or "Fig. 7" in legend_text:
            continue
        paragraph = document.add_paragraph()
        paragraph.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY
        add_runs(paragraph, legend_parts)

    # Supplementary figures are supplied as separate PDF files, but their
    # legends belong in the manuscript after the references. Keep only figures
    # that are explicitly cited in the main text.
    supplementary_legends = {
        "S1": (
            "Figure S1. Orthogonal METANet analysis. A, Predicted edge rank for canonical TBX21 "
            "and PAX5 targets in the whole-blood network and absence from the corresponding "
            "tissue-specific network. Lower rank percentiles indicate stronger predicted edges. "
            "B, Tissue-specific target-set overlap with the discovery modules. C, Top predicted "
            "upstream TFs of TBX21 and PAX5 in the whole-blood tissue-specific map. METANet is an "
            "orthogonal prior and does not establish disease-specific causality [25]."
        ),
        "S2": (
            "Figure S2. Acute Mtb-stimulation pseudobulk analysis. The analysis uses donor-aware "
            "linear mixed models and does not treat cells as replicates. It is retained as "
            "supplementary cell-state evidence only. Acute in vitro stimulation is not a treatment "
            "model and is not interpreted as contradicting the in vivo longitudinal recovery "
            "analyses. Asterisks indicate FDR * < 0.05, ** < 0.01, *** < 0.001."
        ),
        "S3": (
            "Figure S3. Time-adjusted within-participant association of TBRI with "
            "quanTIseq-estimated B-cell fraction in GSE89403. The plot displays residualized "
            "values after adjusting both axes for categorical treatment time and participant. "
            "The partial Pearson r describes the residual association; the displayed P value "
            "comes from a mixed model with participant mean B-cell level, categorical treatment "
            "time, and participant random intercept (82 participants, 312 records). Both "
            "quantities derive from the same bulk transcriptomes."
        ),
        "S4": (
            "Figure S4. Paired Week 26-minus-baseline changes in estimated B-cell "
            "and neutrophil fractions in GSE31348. Points show mean within-patient changes "
            "and bars show 95% paired t confidence intervals (27 participants). CIBERSORT "
            "and quanTIseq use the same bulk expression profiles; neither measures absolute "
            "cell counts. The simultaneous neutrophil decline limits interpretation of the "
            "relative B-cell increase."
        ),
        "S5": (
            "Figure S5. In vivo airway single-cell context for candidate TFs (GSE326212). "
            "A, Patient-sample-level pseudobulk expression across coarse BAL cell types. "
            "Colour indicates relative expression scaled within each TF, and dot size indicates "
            "the median percentage of expressing cells. B, Sample-level contrasts between active "
            "tuberculosis and IGRA-negative contacts. Colour shows the effect in standard-error "
            "units; asterisks indicate BH-adjusted P values (* < 0.05, ** < 0.01, *** < 0.001). "
            "Cells were aggregated by sample and cell type before all comparisons. This is "
            "cross-sectional in vivo evidence, not a treatment trajectory. Raw-count threshold "
            "sensitivity is reported in the Supplementary Information."
        ),
    }
    # The editorial restructuring below inserts the final Results order and
    # carries the citations for S1-S5 into the finished manuscript. Keep only
    # those five figures, all of which are explicitly cited in that final text.
    cited_supplement = ["S1", "S2", "S3", "S4", "S5"]
    if cited_supplement:
        add_heading(document, "Supplementary figure legends", 1)
        for key in sorted(cited_supplement, key=lambda x: int(x[1:])):
            add_body(document, supplementary_legends[key])

    add_heading(document, "Table", 1)
    add_cohort_characteristics_table(document)

    from editorial_structure import revise_main
    revise_main(document)
    document.save(DOCS / "Manuscript_Inflammation_Research.docx")


def build_supplementary_information() -> None:
    document = Document()
    configure_document(document, double_space=True)
    title = document.add_paragraph(style="Title")
    title.alignment = WD_ALIGN_PARAGRAPH.CENTER
    title.add_run("Supplementary Information")
    subtitle = document.add_paragraph()
    subtitle.alignment = WD_ALIGN_PARAGRAPH.CENTER
    subtitle.add_run(TITLE)

    add_heading(document, "Contents", 1)
    contents = [
        "Figure S1. Orthogonal METANet regulatory-network analysis.",
        "Figure S2. Acute Mtb-stimulation pseudobulk analysis.",
        "Figure S3. Time-adjusted within-participant TBRI and B-cell association.",
        "Figure S4. Paired B-cell and neutrophil fraction changes in GSE31348.",
    ]
    for line in contents:
        document.add_paragraph(line, style="List Bullet")

    add_heading(document, "Figure S1", 1)
    add_body(
        document,
        "Figure S1. Orthogonal METANet analysis. A, Predicted edge rank for canonical TBX21 "
        "and PAX5 targets in the whole-blood network and absence from the corresponding "
        "tissue-specific network. Lower rank percentiles indicate stronger predicted edges. "
        "B, Tissue-specific target-set overlap with the discovery modules. C, Top predicted "
        "upstream TFs of TBX21 and PAX5 in the whole-blood tissue-specific map. METANet is an "
        "orthogonal prior and does not establish disease-specific causality [25]."
    )

    add_heading(document, "Figure S2", 1)
    add_body(
        document,
        "Figure S2. Acute Mtb-stimulation pseudobulk analysis. The analysis uses donor-aware "
        "linear mixed models and does not treat cells as replicates. It is retained as "
        "supplementary cell-state evidence only. Acute in vitro stimulation is not a treatment "
        "model and is not interpreted as contradicting the in vivo longitudinal recovery "
        "analyses. Asterisks indicate FDR * < 0.05, ** < 0.01, *** < 0.001."
    )

    add_heading(document, "Figure S3", 1)
    add_body(document, "Figure S3. Time-adjusted within-participant association of TBRI with "
             "quanTIseq-estimated B-cell fraction in GSE89403. The plot displays residualized "
             "values after adjusting both axes for categorical treatment time and participant. "
             "The partial Pearson r describes the residual association; the displayed P value "
             "comes from a mixed model with participant mean B-cell level, categorical treatment "
             "time, and participant random intercept (82 participants, 312 records). Both "
             "quantities derive from the same bulk transcriptomes.")

    add_heading(document, "Figure S4", 1)
    add_body(document, "Figure S4. Paired Week 26-minus-baseline changes in estimated B-cell "
             "and neutrophil fractions in GSE31348. Points show mean within-patient changes "
             "and bars show 95% paired t confidence intervals (27 participants). CIBERSORT "
             "and quanTIseq use the same bulk expression profiles; neither measures absolute "
             "cell counts. The simultaneous neutrophil decline limits interpretation of the "
             "relative B-cell increase.")

    from editorial_structure import revise_supplement
    revise_supplement(document, figure_embed, add_body, add_heading)
    document.save(DOCS / "Supplementary_Information.docx")


def main() -> None:
    DOCS.mkdir(parents=True, exist_ok=True)
    SUPP.mkdir(parents=True, exist_ok=True)
    build_title_page()
    build_cover_letter()
    build_main_manuscript()


if __name__ == "__main__":
    main()



