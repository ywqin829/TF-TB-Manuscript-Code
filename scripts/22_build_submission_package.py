#!/usr/bin/env python
"""Assemble the upload-ready submission package from authoritative outputs."""

from pathlib import Path
import shutil
import zipfile

PROJECT = Path(__file__).resolve().parents[1]
DOCS = PROJECT / "docs"
FIGS = PROJECT / "figures"
SUPP = PROJECT / "supplement"
PACKAGE = PROJECT / "submission_package_20261031"
ARCHIVE = PROJECT / "TF_TB_Submission_20261031.zip"
ZENODO_ARCHIVE = PROJECT / "Zenodo_Derived_Results.zip"

DERIVED_RESULTS = [
    "cohort_inclusion_flow.tsv", "cohort_timepoint_counts.tsv",
    "cohort_batch_group_counts.tsv", "gse89403_gsm_subject_time_map.tsv",
    "gse31348_fixed_gene_direction.tsv", "gse31348_fixed_gene_direction_summary.tsv",
    "gse31348_paired_cell_fraction_changes.tsv", "gse31348_paired_cell_fraction_effects.tsv",
    "gse31348_cibersort_lm22_legacy_settings.tsv",
    "gse31348_quantiseq_explicit_nonlog.tsv", "key_tf_composition_effects.tsv",
    "pax5_bcell_marker_target_overlap.tsv", "pax5_nonoverlap_marker_sensitivity.tsv",
    "candidate_composition_robustness.tsv", "unbiased_tf_cross_cohort_screen.tsv",
    "fixed_core_genes_discovery_evidence.tsv",
]


def copy_file(source: Path, relative_destination: str) -> None:
    if not source.exists():
        raise FileNotFoundError(source)
    destination = PACKAGE / relative_destination
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, destination)


def build_code_archive() -> None:
    destination = PACKAGE / "Supplementary_Code.zip"
    excluded_script_prefixes = (
        "09_", "10_", "11_", "21_", "22_"
    )
    include_scripts = [
        path for path in sorted((PROJECT / "scripts").glob("*"))
        if path.is_file()
        and not path.name.startswith(excluded_script_prefixes)
    ]
    include_logs = sorted((PROJECT / "logs").glob("*.txt"))
    metadata = [
        PROJECT / "results" / "analysis_evidence_map.tsv",
        PROJECT / "results" / "dataset_manifest.tsv",
        PROJECT / "results" / "sample_counts_by_dataset.tsv",
    ]

    with zipfile.ZipFile(
        destination, "w", compression=zipfile.ZIP_DEFLATED
    ) as archive:
        archive.write(PROJECT / "README.md", "SupplementaryCode/README.md")
        for script in include_scripts:
            if script.is_file():
                archive.write(
                    script,
                    f"SupplementaryCode/scripts/{script.name}",
                )
        for script in sorted((PROJECT / "scripts" / "figure_sources").glob("*")):
            if script.is_file():
                archive.write(script, f"SupplementaryCode/scripts/figure_sources/{script.name}")
        for log in include_logs:
            archive.write(log, f"SupplementaryCode/logs/{log.name}")
        for item in metadata:
            if item.exists():
                archive.write(
                    item,
                    f"SupplementaryCode/metadata/{item.name}",
                )
        for name in DERIVED_RESULTS:
            source = PROJECT / "results" / name
            if not source.exists():
                raise FileNotFoundError(source)
            archive.write(source, f"SupplementaryCode/derived_results/{name}")


def build_zenodo_archive() -> None:
    """Stage derived tables separately from figures and manuscript files."""
    note = (
        "Derived results for the TB transcriptional recovery reanalysis.\n"
        "Original expression data remain in GEO. These tables contain fixed gene-set "
        "validation, paired deconvolution effects, TF contrasts and cohort provenance.\n"
        "GSE89403 has two eligible GSM records per subject-time ID; the supplied "
        "processed expression matrix has one column per ID. Its upstream replicate "
        "combination or selection could not be independently reconstructed.\n"
    )
    with zipfile.ZipFile(ZENODO_ARCHIVE, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        archive.writestr("README.txt", note)
        for name in DERIVED_RESULTS:
            source = PROJECT / "results" / name
            if not source.exists():
                raise FileNotFoundError(source)
            archive.write(source, f"derived_results/{name}")


def main() -> None:
    if PACKAGE.resolve().parent != PROJECT.resolve():
        raise RuntimeError(f"Refusing to clean an unexpected package path: {PACKAGE}")
    PACKAGE.mkdir(parents=True, exist_ok=True)
    (PACKAGE / "Figures").mkdir(parents=True, exist_ok=True)
    (PACKAGE / "Supplementary").mkdir(parents=True, exist_ok=True)
    for folder in (PACKAGE / "Figures", PACKAGE / "Supplementary"):
        for old_file in folder.iterdir():
            if old_file.is_file():
                old_file.unlink()

    copy_file(DOCS / "Manuscript_Inflammation_Research.docx", "Manuscript_Inflammation_Research.docx")
    copy_file(DOCS / "Title_Page.docx", "Title_Page.docx")
    copy_file(DOCS / "Cover_Letter.docx", "Cover_Letter.docx")

    for number, stem in [
        (1, "Figure1_StudyDesign"),
        (2, "Figure2_Discovery"),
        (3, "Figure3_CompositionConfounding"),
        (4, "Figure4_SingleCellContext"),
        (5, "Figure5_OutcomeModels"),
        (6, "Figure6_TBRI_OutcomeAssociation"),
    ]:
        suffix = ".pdf"
        copy_file(FIGS / f"{stem}{suffix}", f"Figures/{stem}{suffix}")

    copy_file(
        FIGS / "FigureS1_METANet.pdf",
        "Supplementary/FigureS1_METANet.pdf",
    )
    copy_file(
        FIGS / "FigureS2_CellTypeSpecificMtbResponse.pdf",
        "Supplementary/FigureS2_CellTypeSpecificMtbResponse.pdf",
    )
    copy_file(
        FIGS / "FigureS3_TBRI_Bcell_TimeAdjusted.pdf",
        "Supplementary/FigureS3_TBRI_Bcell_TimeAdjusted.pdf",
    )
    copy_file(
        FIGS / "FigureS4_GSE31348_CIBERSORT_LM22_Bcell_Trajectory.pdf",
        "Supplementary/FigureS4_GSE31348_CIBERSORT_LM22_Bcell_Trajectory.pdf",
    )

    copy_file(FIGS / "FigureS5_InVivo_Candidate_TF_Context.pdf", "Supplementary/FigureS5_InVivo_Candidate_TF_Context.pdf")

    build_zenodo_archive()
    build_code_archive()

    if ARCHIVE.exists():
        ARCHIVE.unlink()
    with zipfile.ZipFile(
        ARCHIVE, "w", compression=zipfile.ZIP_DEFLATED
    ) as archive:
        for path in sorted(PACKAGE.rglob("*")):
            if path.is_file():
                archive.write(
                    path,
                    path.relative_to(PACKAGE).as_posix(),
                )


if __name__ == "__main__":
    main()
