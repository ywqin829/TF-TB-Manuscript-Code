# Decoding the Ignition-Effector Cascade: TBX21-PAX5 Axis in TB Resolution

**Manuscript**: Respiratory Medicine (2025)

## Study Overview

This repository contains the complete analysis code for the manuscript describing the **Ignition-Effector Cascade** — a mechanism where T/NK cells (TBX21/NFATC2) activate first and drive B-cell (PAX5/EBF1) restoration during successful TB treatment.

## Figure Organization

| Script | Figure | Content | Dataset |
|--------|--------|---------|---------|
| `01_Figure2_GSE40553.R` | Fig 2A-B | Volcano plots (Disease vs Recovery) | GSE40553 |
| | Fig 2C | Venn diagrams (Transcriptional Reversal) | GSE40553 |
| | Fig 2D | Mfuzz temporal clustering | GSE40553 |
| | Fig 2E | GO enrichment (Reactome/GO) | GSE40553 |
| | Fig 2F | VIPER TF activity heatmap | GSE40553 |
| | Fig 2G | Fisher TF enrichment bubble plot | GSE40553 |
| `02_Figure3_SingleCell.R` | Fig 3A-C | PBMC scRNA: TF DotPlot/VlnPlot/FeaturePlot | GSE218065 |
| | Fig 3D-F | Lung scRNA: TF pseudobulk analysis | GSE192483 |
| | Fig 3G-I | MTB PBMC: temporal dynamics + CellChat + Monocle3 | EGAS00001005376 |
| | Fig 3J | Sankey network: Cell → TF → Function | GSE40553 |
| | Fig 3K | 9-gene boxplots across timepoints | GSE40553 |
| `03_Figure4_CIBERSORT_NicheNet.R` | Fig 4A | CIBERSORT immune cell dynamics | GSE31348 |
| | Fig 4B | PAX5/TBX21/NFATC2 longitudinal expression | GSE31348 |
| | Fig 4C | NicheNet ligand-receptor prediction | GSE31348 |
| | Fig 4D | scRNA pseudobulk validation (lung B/Plasma) | GSE192483 |
| `04_Figure5_GSE89403.R` | Fig 5A-F | PAX5/TBX21/NFATC2 activity + target score validation | GSE89403 |
| `05_Figure6_HealingScore.R` | Fig 6A | LASSO coefficient profile | GSE40553 |
| | Fig 6B | LASSO cross-validation | GSE40553 |
| | Fig 6C | Lollipop plot (19 LASSO genes) | GSE40553 |
| | Fig 6D | Diagnostic ROC (GSE31348, AUC=0.992) | GSE31348 |
| | Fig 6E-G | Healing Score: Dynamic ROC, Calibration, DCA | GSE89403 |

## Data Dependencies

| Dataset | Accession | Type | Download |
|---------|-----------|------|----------|
| Discovery cohort | GSE40553 | Illumina HumanHT-12 v4 array | GEO: GSE40553 |
| Validation cohort 1 | GSE31348 | Affymetrix HG-U133 Plus 2.0 | GEO: GSE31348 |
| Outcome cohort | GSE89403 | RNA-seq log2 expression | GEO: GSE89403 |
| Additional validation | GSE19435 | Illumina array | GEO: GSE19435 |
| PBMC scRNA-seq | GSE218065 | 10x Genomics | GEO: GSE218065 |
| Lung scRNA-seq | GSE192483 | 10x Genomics | GEO: GSE192483 |
| MTB timecourse scRNA | EGAS00001005376 | 10x Genomics | EGA |

### NicheNet Dependency
Fig 4C requires external NicheNet database files:
- `ligand_target_matrix.rds`
- `lr_network.rds`
- `weighted_networks.rds`

Download from: https://github.com/saeyslab/nichenetr

## Output Directory

All figures are saved to `output/` as PDF files.

## Session Info

Analysis performed with R 4.x, Seurat v5, and Bioconductor packages.
