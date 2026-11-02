# Run from the project root: Rscript scripts/figure_sources/Figure6_CompositionConfounding.R
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



# Figure 6: Does cell composition explain the apparent TF recovery?
composition_key <- composition_models %>%
  filter(TF %in% c("PAX5", "TBX21", "NFATC2")) %>%
  mutate(
    dataset = factor(dataset, levels = c("GSE40553", "GSE31348", "GSE89403")),
    TF = factor(TF, levels = c("PAX5", "TBX21", "NFATC2"))
  )

p6a <- ggplot(composition_key, aes(y = TF)) +
  geom_segment(
    aes(x = R2_base, xend = R2_with_composition, yend = TF),
    colour = "grey65", linewidth = 1
  ) +
  geom_point(
    aes(x = R2_base, colour = "Base model"),
    size = 3.1
  ) +
  geom_point(
    aes(x = R2_with_composition, colour = "With composition PCs"),
    size = 3.1
  ) +
  facet_wrap(~dataset, nrow = 1) +
  scale_colour_manual(
    values = c(
      "Base model" = "#4393C3",
      "With composition PCs" = "#B2182B"
    )
  ) +
  scale_x_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.08))) +
  labs(
    x = "Marginal R2",
    y = NULL,
    colour = NULL,
    title = "Composition PCs increase model fit"
  ) +
  theme_tf(11) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold")
  )

p6b <- ggplot(
  composition_key,
  aes(dataset, TF, fill = -log10(pmax(FDR, 1e-30)))
) +
  geom_tile(colour = "white", linewidth = 1) +
  geom_text(
    aes(label = if_else(FDR < 0.001, "<0.001", sprintf("%.3f", FDR))),
    size = 3.8, fontface = "bold"
  ) +
  scale_fill_gradientn(
    colours = c("#F7FBFC", "#FEB24C", "#B2182B"),
    name = "-log10 FDR"
  ) +
  labs(x = NULL, y = NULL, title = "Composition model likelihood-ratio test") +
  theme_tf(11) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1),
        legend.position = "right")

score_comp <- outcome$score_data %>%
  mutate(SampleID = as.character(SampleID)) %>%
  left_join(comp_89403, by = "SampleID") %>%
  rename(B_cell_fraction = `B cell`)

tbri_vs_b <- cor.test(
  score_comp$HealingScore,
  score_comp$B_cell_fraction,
  method = "spearman",
  exact = FALSE
)

p6c <- ggplot(
  score_comp,
  aes(B_cell_fraction * 100, HealingScore, colour = Group)
) +
  geom_point(size = 1.8, alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 1.1) +
  scale_colour_manual(
    values = c("Successful_Treatment" = "#1976D2",
               "Poor_Outcome" = "#D32F2F"),
    labels = c("Success", "Poor outcome")
  ) +
  labs(
    x = "Estimated B-cell fraction (%)",
    y = "TB Resolution Index",
    title = "TBRI is strongly correlated with B-cell fraction",
    subtitle = paste0("Spearman rho = ", sprintf("%.3f", tbri_vs_b$estimate))
  ) +
  theme_tf(11) +
  theme(legend.position = "bottom", legend.title = element_blank())

comp_trajectory <- score_comp %>%
  group_by(Group, TimePoint) %>%
  summarise(
    mean = mean(B_cell_fraction),
    se = sd(B_cell_fraction) / sqrt(n()),
    .groups = "drop"
  ) %>%
  mutate(
    TimePoint = factor(
      TimePoint,
      levels = c("0_Baseline", "1_Week1", "2_Week4", "3_Week24"),
      labels = c("0M", "W1", "W4", "W24")
    )
  )

p6d <- ggplot(comp_trajectory, aes(TimePoint, mean, colour = Group, group = Group)) +
  geom_errorbar(aes(ymin = mean - 1.96 * se, ymax = mean + 1.96 * se),
                width = 0.12, linewidth = 0.8,
                position = position_dodge(width = 0.18)) +
  geom_line(linewidth = 1.1, position = position_dodge(width = 0.18)) +
  geom_point(size = 2.8, position = position_dodge(width = 0.18)) +
  scale_colour_manual(
    values = c("Successful_Treatment" = "#1976D2",
               "Poor_Outcome" = "#D32F2F"),
    labels = c("Success", "Poor outcome")
  ) +
  labs(x = NULL, y = "Estimated B-cell fraction",
       title = "B-cell fraction during treatment") +
  theme_tf(11) +
  theme(legend.position = "bottom", legend.title = element_blank())

auc_long <- tbri_increment %>%
  select(TimePoint, AUC_Bcell_only, AUC_Bcell_plus_TBRI) %>%
  pivot_longer(-TimePoint, names_to = "Model", values_to = "AUC") %>%
  mutate(
    TimePoint = factor(
      TimePoint,
      levels = c("0_Baseline", "1_Week1", "2_Week4", "3_Week24"),
      labels = c("0M", "W1", "W4", "W24")
    ),
    Model = recode(
      Model,
      AUC_Bcell_only = "B-cell fraction",
      AUC_Bcell_plus_TBRI = "B-cell fraction + TBRI"
    )
  )

p6e <- ggplot(auc_long, aes(TimePoint, AUC, colour = Model, group = Model)) +
  geom_hline(yintercept = 0.5, linetype = 2, colour = "grey60") +
  geom_line(linewidth = 1.1) +
  geom_point(size = 3) +
  scale_colour_manual(values = c(
    "B-cell fraction" = "#4393C3",
    "B-cell fraction + TBRI" = "#B2182B"
  )) +
  labs(x = NULL, y = "Outcome AUC",
       title = "Incremental discrimination") +
  theme_tf(11) +
  theme(legend.position = "bottom", legend.title = element_blank())

p6f <- tbri_increment %>%
  mutate(
    TimePoint = factor(
      TimePoint,
      levels = c("0_Baseline", "1_Week1", "2_Week4", "3_Week24"),
      labels = c("0M", "W1", "W4", "W24")
    )
  ) %>%
  ggplot(aes(TimePoint, delta_AUC)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey60") +
  geom_col(fill = "#8073AC", width = 0.65) +
  geom_text(aes(label = sprintf("%+.3f", delta_AUC)),
            vjust = -0.35, size = 3.8, fontface = "bold") +
  scale_y_continuous(limits = c(-0.012, 0.055), breaks = c(0, 0.02, 0.04)) +
  labs(x = NULL, y = "Delta AUC", title = "TBRI gain over B-cell fraction") +
  theme_tf(11)

figure6 <- (p6a | p6b) / (p6c | p6d) / (p6e | p6f) +
  plot_annotation(tag_levels = "A") &
  panel_tag_theme

save_plot_pair(figure6, "Figure6_CompositionConfounding", 13, 13)

