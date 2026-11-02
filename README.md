# Composition-aware tuberculosis treatment transcriptome reanalysis

This repository contains the analysis and figure-generation code for the manuscript **Immune cell composition shapes the interpretation of blood transcriptional recovery during tuberculosis treatment**.

## Scope

The study is a secondary analysis of public transcriptomic datasets. No new samples or raw sequencing data were generated.

- Discovery and longitudinal recovery: **GSE40553**
- Independent directional validation: **GSE31348**
- Longitudinal outcome and TBRI/composition analyses: **GSE89403**
- Additional single-cell context: **GSE192483**, **GSE218065**, **GSE326212**
- Acute stimulation supplementary context: **EGAS00001005376**

## Repository layout

- scripts/: preprocessing, statistical analyses, sensitivity analyses, and manuscript figure sources
- scripts/figure_sources/: R source files for the integrated main and supplementary figures
- results/, figures/, and tables/ are generated locally and are not committed to this code repository
- No participant-level data are distributed here; obtain source data from the accession records above

## Reproducibility

The project-specific runner is scripts/run_all.ps1. Before running it, set TF_TB_PROJECT_ROOT and the input environment variables consumed by scripts/00_config.R:

- TF_TB_GSE40553_RDS
- TF_TB_GSE31348_RDS
- TF_TB_GSE89403_RDS

Additional raw and single-cell inputs are documented in the corresponding scripts. Some analyses require external reference resources such as LM22/CIBERSORT or NicheNet databases; those are downloaded or configured separately and are not redistributed here.

The analysis uses R and Python. Exact package versions and run logs should be archived with the manuscript submission or a DOI-backed data record.

## Data availability

All source datasets are publicly available through NCBI GEO or EGA under the accession numbers listed above. Derived result tables are generated locally from the supplied scripts and are provided in the manuscript submission archive; source data remain available from the accession records above.

## Code availability

Public repository: https://github.com/ywqin829/TF-TB-Manuscript-Code

## Legacy files

The original top-level scripts are retained for provenance. The supported reproducibility workflow is the code under `scripts/`; legacy scripts should not be used for the current manuscript without review.


