# Run from the project root: Rscript scripts/figure_sources/FigureS4_TBRI_Calibration_DCA.R
# Generated from scripts/06_figures.R; uses existing model outputs.

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


score_auc <- read_tsv(out_path("results", "gse89403_score_external_auc.tsv"), show_col_types = FALSE) %>%
  mutate(
    TimePoint = factor(
      TimePoint,
      levels = c("0_Baseline", "1_Week1", "2_Week4", "3_Week24"),
      labels = c("0M", "W1", "W4", "W24")
    )
  )

p5a <- ggplot(outcome$score_data, aes(TimePoint, HealingScore, fill = Group)) +
  geom_boxplot(width = 0.62, outlier.shape = NA, alpha = 0.75,
               position = position_dodge(width = 0.72)) +
  geom_point(
    aes(colour = Group),
    position = position_jitterdodge(jitter.width = 0.14, dodge.width = 0.72),
    size = 1.4, alpha = 0.55, show.legend = FALSE
  ) +
  scale_fill_manual(values = c(
    "Successful_Treatment" = "#1976D2",
    "Poor_Outcome" = "#D32F2F"
  ), labels = c("Success", "Poor outcome")) +
  scale_colour_manual(values = c(
    "Successful_Treatment" = "#1976D2",
    "Poor_Outcome" = "#D32F2F"
  )) +
  labs(x = NULL, y = "TB Resolution Score", title = "Score distribution") +
  theme_tf(11) +
  theme(legend.position = "bottom", legend.title = element_blank())

p5b <- ggplot(score_auc, aes(TimePoint, AUC)) +
  geom_hline(yintercept = 0.5, linetype = 2, colour = "grey55") +
  geom_errorbar(aes(ymin = lower_95, ymax = upper_95), width = 0.12,
                linewidth = 0.9, colour = "#2166AC") +
  geom_point(size = 3.2, colour = "#2166AC") +
  geom_text(aes(label = sprintf("%.2f", AUC)), nudge_y = 0.07, size = 3.8) +
  coord_cartesian(ylim = c(0.35, 1.02)) +
  labs(x = NULL, y = "AUC (95% CI)", title = "Outcome discrimination") +
  theme_tf(11)

week1 <- outcome$score_predictions %>% filter(TimePoint == "1_Week1")
week1_roc <- roc(
  factor(week1$Group, levels = c("Poor_Outcome", "Successful_Treatment")),
  week1$HealingScore,
  levels = c("Poor_Outcome", "Successful_Treatment"),
  direction = "<",
  quiet = TRUE
)
roc_df <- data.frame(
  FPR = 1 - week1_roc$specificities,
  TPR = week1_roc$sensitivities
)

p5c <- ggplot(roc_df, aes(FPR, TPR)) +
  geom_abline(linetype = 2, colour = "grey55") +
  geom_path(linewidth = 1.2, colour = "#B2182B") +
  annotate("text", x = 0.64, y = 0.22,
           label = paste0("AUC = ", sprintf("%.3f", as.numeric(auc(week1_roc)))),
           size = 4.3, fontface = "bold") +
  coord_equal() +
  labs(x = "1 - specificity", y = "Sensitivity", title = "Week 1 outcome ROC") +
  theme_tf(11)

calibration_data <- week1 %>%
  mutate(
    outcome_poor = as.integer(Group == "Poor_Outcome"),
    bin = ntile(loo_probability_poor, 5)
  ) %>%
  group_by(bin) %>%
  summarise(
    predicted = mean(loo_probability_poor),
    observed = mean(outcome_poor),
    n = n(),
    .groups = "drop"
  )

p5d <- ggplot(calibration_data, aes(predicted, observed)) +
  geom_abline(linetype = 2, colour = "grey55") +
  geom_line(linewidth = 1, colour = "#2166AC") +
  geom_point(aes(size = n), colour = "#2166AC") +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(x = "Predicted probability of poor outcome",
       y = "Observed proportion", title = "Week 1 calibration", size = "n") +
  theme_tf(11) +
  theme(legend.position = "none")

dca_week1 <- outcome$dca %>%
  filter(TimePoint == "1_Week1") %>%
  mutate(
    Strategy = case_when(
      variable == "loo_probability_poor" ~ "TBRI",
      variable == "all" ~ "Treat all",
      TRUE ~ "Treat none"
    )
  )

p5e <- ggplot(dca_week1, aes(threshold, net_benefit, colour = Strategy)) +
  geom_line(linewidth = 1.1) +
  scale_colour_manual(values = c(
    "TBRI" = "#B2182B",
    "Treat all" = "#4393C3",
    "Treat none" = "grey45"
  )) +
  coord_cartesian(ylim = c(-0.03, 0.12)) +
  labs(x = "Threshold probability", y = "Net benefit",
       title = "Week 1 decision curve") +
  theme_tf(11) +
  theme(legend.position = "bottom", legend.title = element_blank())

component_data <- outcome$score_data %>%
  group_by(TimePoint) %>%
  summarise(
    n_success = sum(Group == "Successful_Treatment"),
    n_poor = sum(Group == "Poor_Outcome"),
    median_score = median(HealingScore),
    .groups = "drop"
  )

p5f <- ggplot(component_data, aes(TimePoint, median_score, group = 1)) +
  geom_line(linewidth = 1.1, colour = "#6A51A3") +
  geom_point(size = 3, colour = "#6A51A3") +
  geom_text(aes(label = paste0("S", n_success, "/P", n_poor)),
            nudge_y = 0.22, size = 3.4) +
  labs(x = NULL, y = "Median score", title = "Whole-cohort score trajectory") +
  theme_tf(11)

supp_figure4 <- plot_grid(
  p5d, p5e,
  ncol = 2,
  labels = c("A", "B"),
  label_size = 20
)
save_plot_pair(supp_figure4, "FigureS4_TBRI_Calibration_DCA", 15, 7)

