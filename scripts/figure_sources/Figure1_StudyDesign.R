# Run from the project root: Rscript scripts/figure_sources/Figure1_StudyDesign.R
# Generated from scripts/06_figures.R by scripts/39_export_integrated_figure_sources.py.
# Uses existing analysis tables and source-context panels; no upstream analysis is rerun.

source(file.path("scripts", "00_config.R"))

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
  library(forcats)
  library(patchwork)
  library(cowplot)
  library(grid)
  library(scales)
  library(emmeans)
  library(lme4)
  library(pROC)
  library(pdftools)
  library(magick)
})

discovery <- readRDS(out_path("results", "gse40553_discovery_results.rds"))
development <- readRDS(out_path("results", "gse31348_development_results.rds"))
outcome <- readRDS(out_path("results", "gse89403_outcome_results.rds"))

disease <- read_tsv(out_path("results", "gse40553_disease_0m_vs_ltb.tsv"), show_col_types = FALSE)
treatment <- read_tsv(out_path("results", "gse40553_treatment_6m_vs_0m.tsv"), show_col_types = FALSE)
tf_table <- read_tsv(out_path("results", "gse40553_tf_activity_subject_blocked.tsv"), show_col_types = FALSE)
donor_cor <- read_tsv(out_path("results", "gse40553_paired_donor_correlations.tsv"), show_col_types = FALSE)
metanet_direct <- read_tsv(out_path("results", "metanet_direct_edges.tsv"), show_col_types = FALSE)
metanet_overlap <- read_tsv(out_path("results", "metanet_module_overlap.tsv"), show_col_types = FALSE)
metanet_upstream <- read_tsv(out_path("results", "metanet_upstream_candidates.tsv"), show_col_types = FALSE)
composition_models <- read_tsv(out_path("results", "tf_activity_composition_models.tsv"), show_col_types = FALSE)
tbri_increment <- read_tsv(out_path("results", "tbri_incremental_value.tsv"), show_col_types = FALSE)
comp_89403 <- read_tsv(out_path("results", "gse89403_quantiseq.tsv"), show_col_types = FALSE)
pb_long <- read_tsv(out_path("results", "mtb_pbmc_pseudobulk_long.tsv"), show_col_types = FALSE)
pb_models <- read_tsv(out_path("results", "mtb_pbmc_celltype_time_models.tsv"), show_col_types = FALSE)

figure_dir <- out_path("figures")

panel_tag_theme <- theme(
  plot.tag = element_text(face = "bold", size = 20),
  plot.tag.position = c(0.01, 0.99)
)

volcano_plot <- function(data, title, module_colors) {
  data <- data %>%
    mutate(
      neglog10_fdr = pmin(-log10(pmax(adj.P.Val, 1e-50)), 25),
      Direction = case_when(
        adj.P.Val < 0.05 & logFC > 1 ~ "Higher after comparison",
        adj.P.Val < 0.05 & logFC < -1 ~ "Lower after comparison",
        TRUE ~ "Not significant"
      )
    )
  label_genes <- c("PAX5", "TBX21", "NFATC2", "EBF1", "CD79A", "CD19",
                   "IFNG", "CD40LG", "IRF1", "STAT1")
  ggplot(data, aes(logFC, neglog10_fdr)) +
    geom_point(aes(colour = Module, alpha = Direction == "Not significant"),
               size = 0.7) +
    scale_alpha_manual(values = c(`TRUE` = 0.22, `FALSE` = 0.85), guide = "none") +
    scale_colour_manual(values = module_colors, drop = FALSE) +
    geom_vline(xintercept = c(-1, 1), linetype = 2, colour = "grey60") +
    geom_hline(yintercept = -log10(0.05), linetype = 2, colour = "grey60") +
    ggrepel::geom_text_repel(
      data = data %>% filter(Gene %in% label_genes),
      aes(label = Gene),
      size = 4,
      min.segment.length = 0,
      max.overlaps = 30,
      seed = 20261028
    ) +
    coord_cartesian(xlim = c(-4, 4)) +
    labs(x = "log2 fold change", y = expression(-log[10]~FDR), title = title) +
    theme_tf(12) +
    theme(legend.position = "bottom", legend.title = element_blank())
}

module_colors <- c(
  "Restoration core" = "#2C7FB8",
  "Disease driver" = "#D95F0E",
  "Other" = "#BDBDBD"
)



figure1_boxes <- tribble(
  ~stage, ~xmin, ~xmax, ~ymin, ~ymax, ~label, ~fill,
  "Discovery", 0.5, 2.5, 5.3, 7.0,
  "GSE40553\n74 patients | 176 samples\nTranscriptional reversal\nTF-regulon activity",
  "#DEEBF7",
  "Development", 3.0, 5.0, 5.3, 7.0,
  "GSE31348\n27 patients | 135 samples\nModule dynamics\nScore development",
  "#FEE0D2",
  "External outcome", 5.5, 7.5, 5.3, 7.0,
  "GSE89403\n82 patients | 312 samples\nOutcome-linked\nTBRI association",
  "#E5F5E0",
  "In vivo airway", 0.4, 2.7, 2.0, 3.7,
  "GSE326212\n17 libraries | 77 samples\nCross-sectional BAL scRNA-seq\nDisease-state context",
  "#FFF7BC",
  "Lesional context", 3.0, 5.2, 2.0, 3.7,
  "GSE192483\n6 patients | 11 tissues\nCell-lineage localization\nand lesion gradients",
  "#FDD0A2",
  "Regulatory prior", 5.5, 7.6, 2.0, 3.7,
  "METANet\nWhole blood + lung\nOrthogonal TF-target\ncandidate edges",
  "#EDE7F6"
)

figure1_arrows <- tribble(
  ~x, ~y, ~xend, ~yend,
  2.5, 6.15, 3.0, 6.15,
  5.0, 6.15, 5.5, 6.15,
  2.0, 5.3, 1.6, 3.7,
  4.5, 5.3, 4.1, 3.7,
  6.3, 5.3, 6.3, 3.7
)

p1 <- ggplot() +
  geom_rect(
    data = figure1_boxes,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill),
    colour = "black", linewidth = 0.8
  ) +
  geom_text(
    data = figure1_boxes,
    aes((xmin + xmax) / 2, (ymin + ymax) / 2, label = label),
    size = 4.4, lineheight = 1.0
  ) +
  geom_segment(
    data = figure1_arrows,
    aes(x = x, y = y, xend = xend, yend = yend),
    arrow = arrow(length = unit(0.18, "inches"), type = "closed"),
    linewidth = 1
  ) +
  annotate(
    "text", x = 4.0, y = 0.9,
    label = "All longitudinal inference uses patients as the biological unit",
    fontface = "bold", size = 5.2
  ) +
  scale_fill_identity() +
  coord_cartesian(xlim = c(0, 8), ylim = c(0.5, 7.5), clip = "off") +
  theme_void(base_size = 14)

save_plot_pair(p1, "Figure1_StudyDesign", 12, 7.5)

