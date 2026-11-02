# Run from the project root: Rscript scripts/figure_sources/FigureS1_METANet.R
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



# Supplementary figure: METANet orthogonal regulatory prior.
direct_plot <- metanet_direct %>%
  filter(network %in% c("whole_blood_full", "whole_blood_tissue_specific")) %>%
  filter(TF %in% c("TBX21", "PAX5")) %>%
  mutate(
    rank_percentile = if_else(
      edge_present,
      rank_within_TF / outgoing_edges,
      NA_real_
    ),
    network = recode(
      network,
      whole_blood_full = "Whole blood",
      whole_blood_tissue_specific = "Whole blood\ntissue-specific"
    )
  )

ps1a <- ggplot(direct_plot, aes(target, network)) +
  geom_tile(aes(fill = rank_percentile), colour = "white", linewidth = 0.8) +
  geom_text(
    data = direct_plot %>% filter(is.na(rank_percentile)),
    label = "absent", size = 3.4, colour = "grey35"
  ) +
  facet_wrap(~TF, ncol = 1) +
  scale_fill_gradientn(
    colours = c("#4575B4", "#FFFFBF", "#D73027"),
    na.value = "#F1F1F1",
    trans = "reverse",
    name = "Target-rank\npercentile"
  ) +
  labs(x = NULL, y = NULL, title = "Canonical edges in METANet") +
  theme_tf(10) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1))

overlap_plot <- metanet_overlap %>%
  filter(network != "whole_blood_full", is.finite(p_value)) %>%
  mutate(
    contrast = paste(TF, module, sep = " | "),
    contrast = factor(contrast, levels = rev(unique(contrast))),
    network = recode(
      network,
      whole_blood_tissue_specific = "Whole blood",
      lung_tissue_specific = "Lung"
    )
  )

ps1b <- ggplot(overlap_plot, aes(fold_enrichment, contrast, colour = network)) +
  geom_vline(xintercept = 1, linetype = 2, colour = "grey60") +
  geom_point(aes(size = overlap), position = position_dodge(width = 0.35)) +
  scale_x_log10() +
  labs(x = "Fold enrichment of discovery-module genes", y = NULL,
       title = "METANet target-set overlap", colour = "Network", size = "Genes") +
  theme_tf(10)

upstream_plot <- metanet_upstream %>%
  filter(network == "whole_blood_tissue_specific", target_TF %in% c("PAX5", "TBX21")) %>%
  group_by(target_TF) %>%
  slice_head(n = 8) %>%
  ungroup() %>%
  mutate(
    candidate_upstream_TF = factor(
      candidate_upstream_TF,
      levels = rev(unique(candidate_upstream_TF))
    )
  )

ps1c <- ggplot(upstream_plot, aes(score, candidate_upstream_TF, colour = target_TF)) +
  geom_segment(aes(x = 0.45, xend = score, yend = candidate_upstream_TF),
               linewidth = 0.7) +
  geom_point(size = 2.8) +
  facet_wrap(~target_TF, scales = "free_y", ncol = 2) +
  labs(x = "Predicted edge score", y = NULL,
       title = "Whole-blood upstream candidates", colour = "Target TF") +
  theme_tf(10) +
  theme(legend.position = "none")

supp_figure1 <- ps1a / (ps1b | ps1c) +
  plot_layout(heights = c(1.2, 1)) +
  plot_annotation(tag_levels = "A") &
  panel_tag_theme

save_plot_pair(supp_figure1, "FigureS1_METANet", 15, 11)

