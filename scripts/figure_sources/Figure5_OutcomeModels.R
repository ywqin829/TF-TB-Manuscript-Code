# Run from the project root: Rscript scripts/figure_sources/Figure5_OutcomeModels.R
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



trajectory_panel <- function(feature, title, p_interaction) {
  fit <- outcome$models[[feature]]$unadjusted$fit
  emm <- as.data.frame(emmeans(fit, ~ Group | TimePoint))
  emm$TimePoint <- factor(
    emm$TimePoint,
    levels = c("0_Baseline", "1_Week1", "2_Week4", "3_Week24"),
    labels = c("0M", "W1", "W4", "W24")
  )
  emm$Group <- factor(
    emm$Group,
    levels = c("Successful_Treatment", "Poor_Outcome"),
    labels = c("Success", "Poor outcome")
  )
  ggplot(emm, aes(TimePoint, emmean, colour = Group, group = Group)) +
    geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                  width = 0.13, linewidth = 0.8,
                  position = position_dodge(width = 0.18)) +
    geom_line(linewidth = 1.1, position = position_dodge(width = 0.18)) +
    geom_point(size = 2.7, position = position_dodge(width = 0.18)) +
    scale_colour_manual(values = c("Success" = "#1976D2", "Poor outcome" = "#D32F2F")) +
    labs(
      x = NULL, y = "Model-estimated value", colour = "Treatment outcome", title = title,
      subtitle = paste0("Group x time P = ", format.pval(p_interaction, digits = 2))
    ) +
    theme_tf(11) +
    theme(legend.position = "bottom")
}

interaction_p <- function(feature) {
  outcome$model_anova %>%
    filter(
      feature == !!feature,
      term == "Group:TimePoint",
      adjustment == "unadjusted"
    ) %>%
    pull(`Pr(>F)`) %>%
    first()
}

feature_names <- c(
  PAX5 = "PAX5 activity",
  TBX21 = "TBX21 activity",
  PAX5_target_module = "PAX5 target module",
  B_cell_marker_module = "B-cell marker module",
  TBX21_target_module = "TBX21 target module"
)

trajectory_plots <- lapply(names(feature_names), function(feature) {
  trajectory_panel(feature, feature_names[[feature]], interaction_p(feature))
})

score_fit_emm <- as.data.frame(emmeans(outcome$score_mixed, ~ Group | TimePoint))
score_fit_emm$TimePoint <- factor(
  score_fit_emm$TimePoint,
  levels = c("0_Baseline", "1_Week1", "2_Week4", "3_Week24"),
  labels = c("0M", "W1", "W4", "W24")
)
score_fit_emm$Group <- factor(
  score_fit_emm$Group,
  levels = c("Successful_Treatment", "Poor_Outcome"),
  labels = c("Success", "Poor outcome")
)
score_interaction <- outcome$score_anova %>%
  filter(term == "Group:TimePoint") %>%
  pull(`Pr(>F)`) %>%
  first()

p4f <- ggplot(score_fit_emm, aes(TimePoint, emmean, colour = Group, group = Group)) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                width = 0.13, linewidth = 0.8,
                position = position_dodge(width = 0.18)) +
  geom_line(linewidth = 1.1, position = position_dodge(width = 0.18)) +
  geom_point(size = 2.7, position = position_dodge(width = 0.18)) +
  scale_colour_manual(values = c("Success" = "#1976D2", "Poor outcome" = "#D32F2F")) +
  labs(
    x = NULL, y = "Model-estimated score", colour = "Treatment outcome", title = "TB Resolution Index",
    subtitle = paste0("Group x time P = ", format.pval(score_interaction, digits = 2))
  ) +
  theme_tf(11) +
  theme(legend.position = "bottom")

figure4_legend <- cowplot::get_legend(
  p4f + theme(legend.position = "bottom", legend.justification = "center")
)
figure4_panels <- wrap_plots(
  lapply(c(trajectory_plots, list(p4f)), function(p) p + theme(legend.position = "none")),
  ncol = 3
) + plot_annotation(tag_levels = "A") & panel_tag_theme
figure4 <- cowplot::plot_grid(
  figure4_panels, figure4_legend, ncol = 1, rel_heights = c(1, 0.07)
)

save_plot_pair(figure4, "Figure5_OutcomeModels", 15, 10)

