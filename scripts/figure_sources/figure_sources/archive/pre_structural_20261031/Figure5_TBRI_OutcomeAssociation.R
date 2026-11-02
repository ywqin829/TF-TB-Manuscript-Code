# Run from the project root: Rscript scripts/figure_sources/Figure5_TBRI_OutcomeAssociation.R
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



score_auc <- read_tsv(out_path("results", "gse89403_score_external_auc.tsv"), show_col_types = FALSE) %>%
  mutate(
    TimePoint = factor(
      TimePoint,
      levels = c("0_Baseline", "1_Week1", "2_Week4", "3_Week24"),
      labels = c("Baseline", "Week 1", "Week 4", "Week 24")
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
  theme_tf(11) +
  scale_x_discrete(labels = c(
    "0_Baseline" = "Baseline",
    "1_Week1" = "Week 1",
    "2_Week4" = "Week 4",
    "3_Week24" = "Week 24"
  ))

week1 <- outcome$score_predictions %>% filter(TimePoint == "1_Week1")
# Display both available early-treatment visits with fixed ROC orientation.
roc_visits <- c("1_Week1", "2_Week4")
roc_labels <- character(2)
roc_parts <- vector("list", 2)
auc_reference <- readr::read_tsv(out_path("results", "gse89403_score_external_auc.tsv"), show_col_types = FALSE)
for (i in seq_along(roc_visits)) {
  visit <- outcome$score_predictions %>% filter(TimePoint == roc_visits[i])
  visit_roc <- roc(visit$Group, visit$HealingScore,
                   levels = c("Poor_Outcome", "Successful_Treatment"),
                   direction = "<", quiet = TRUE)
  visit_ci <- as.numeric(ci.auc(visit_roc, method = "delong"))
  visit_auc <- as.numeric(auc(visit_roc))
  ref <- auc_reference %>% filter(TimePoint == roc_visits[i])
  stopifnot(nrow(ref) == 1L, abs(visit_auc - ref$AUC) < 1e-8,
            abs(visit_ci[1] - ref$lower_95) < 1e-8,
            abs(visit_ci[3] - ref$upper_95) < 1e-8)
  roc_labels[i] <- sprintf("Week %s: AUC %.3f (95%% CI %.3f-%.3f); poor outcome n = %d",
                           c("1", "4")[i], visit_auc, visit_ci[1], visit_ci[3],
                           sum(visit$Group == "Poor_Outcome"))
  roc_parts[[i]] <- data.frame(FPR = 1 - visit_roc$specificities,
                              TPR = visit_roc$sensitivities, Visit = c("Week 1", "Week 4")[i])
}
roc_df <- bind_rows(roc_parts)
p5c <- ggplot(roc_df, aes(FPR, TPR, colour = Visit, linetype = Visit)) +
  geom_abline(linetype = 2, colour = "grey65") +
  geom_path(linewidth = 1.1) +
  scale_colour_manual(values = c("Week 1" = "#B2182B", "Week 4" = "#2166AC")) +
  scale_linetype_manual(values = c("Week 1" = "solid", "Week 4" = "longdash")) +
  coord_equal() +
  labs(x = "1 - specificity", y = "Sensitivity", title = "Early-treatment outcome ROC",
       colour = NULL, linetype = NULL,
       caption = paste(roc_labels, collapse = "\n")) +
  theme_tf(11) +
  theme(legend.position = "bottom", legend.margin = margin(0, 0, 0, 0),
        plot.caption = element_text(size = 10, hjust = 0, lineheight = 1.15))

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

# Three complementary panels: distributions, AUC estimates, and early ROC curves.
# Fix each plotting region to a square; canvas space must not stretch panels.
p5a <- p5a +
  scale_x_discrete(labels = c("0_Baseline" = "Baseline", "1_Week1" = "Week 1",
                              "2_Week4" = "Week 4", "3_Week24" = "Week 24")) +
  theme(aspect.ratio = 1)
p5b <- p5b + theme(aspect.ratio = 1)
p5c <- p5c + labs(caption = paste(sub("; poor outcome", "\nPoor outcome", roc_labels), collapse = "\n"))
figure5 <- (p5a | p5b | p5c) +
  plot_layout(widths = c(1, 1, 1)) +
  plot_annotation(tag_levels = "A") &
  panel_tag_theme

save_plot_pair(figure5, "Figure5_TBRI_OutcomeAssociation", 15, 6)

