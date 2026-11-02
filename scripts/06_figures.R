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

p2a <- volcano_plot(disease, "Active TB vs latent TB", module_colors)
p2b <- volcano_plot(treatment, "6 months vs baseline", module_colors)

module_long <- discovery$module_scores %>%
  pivot_longer(
    cols = c(Restoration_Core, Disease_Driver),
    names_to = "Module",
    values_to = "Score"
  ) %>%
  mutate(
    Module = recode(
      Module,
      Restoration_Core = "Restoration core",
      Disease_Driver = "Disease driver"
    ),
    Group = factor(
      Group,
      levels = c("LTB", "TB_0M", "TB_0.5M", "TB_2M", "TB_6M", "TB_12M")
    )
  )

module_summary <- module_long %>%
  group_by(Module, Group) %>%
  summarise(
    mean = mean(Score),
    se = sd(Score) / sqrt(n()),
    .groups = "drop"
  )

module_panel <- function(module_name) {
  ggplot(
    filter(module_summary, Module == module_name),
    aes(Group, mean, group = 1, colour = Module)
  ) +
    geom_errorbar(aes(ymin = mean - 1.96 * se, ymax = mean + 1.96 * se),
                  width = 0.15, linewidth = 0.8) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 3) +
    scale_colour_manual(values = module_colors[module_name], guide = "none") +
    labs(x = NULL, y = "ssGSEA score", title = module_name) +
    theme_tf(11) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1))
}

p2c <- module_panel("Restoration core") | module_panel("Disease driver")

tf_plot_data <- tf_table %>%
  filter(is.finite(treatment_logFC)) %>%
  arrange(treatment_FDR) %>%
  slice_head(n = 18) %>%
  mutate(
    TF = factor(TF, levels = rev(TF)),
    Enriched = if_else(significant_target_enrichment, "Targets enriched", "Not enriched")
  )

p2d <- ggplot(tf_plot_data, aes(treatment_logFC, TF, colour = Enriched)) +
  geom_vline(xintercept = 0, colour = "grey50", linewidth = 0.5) +
  geom_segment(aes(x = 0, xend = treatment_logFC, yend = TF), linewidth = 0.8) +
  geom_point(size = 3.2) +
  scale_colour_manual(values = c(
    "Targets enriched" = "#D73027",
    "Not enriched" = "#4575B4"
  )) +
  labs(x = "VIPER TF activity-score change, 6M vs baseline", y = NULL,
       title = "Regulon activity reversal") +
  theme_tf(11) +
  theme(legend.position = "bottom", legend.title = element_blank())

donor_plot_data <- donor_cor %>%
  filter(is.finite(rho_restoration)) %>%
  arrange(desc(rho_restoration)) %>%
  slice_head(n = 18) %>%
  mutate(TF = factor(TF, levels = rev(TF)))

p2e <- ggplot(donor_plot_data, aes(rho_restoration, TF)) +
  geom_vline(xintercept = 0, colour = "grey50", linewidth = 0.5) +
  geom_segment(aes(x = 0, xend = rho_restoration, yend = TF), linewidth = 0.8,
               colour = "#2C7FB8") +
  geom_point(aes(colour = FDR_restoration < 0.05), size = 3.2) +
  scale_colour_manual(values = c(`TRUE` = "#D73027", `FALSE` = "grey55"),
                      labels = c(`TRUE` = "FDR < 0.05", `FALSE` = "FDR ≥ 0.05")) +
  labs(x = "Spearman rho, paired donor changes", y = NULL,
       title = "TF activity and restoration program") +
  theme_tf(11) +
  theme(legend.position = "bottom", legend.title = element_blank())

fisher_plot_data <- tf_table %>%
  filter(significant_activity, is.finite(fold_enrichment), is.finite(FDR)) %>%
  arrange(FDR) %>%
  slice_head(n = 12) %>%
  mutate(
    TF = factor(TF, levels = rev(TF)),
    Enrichment = if_else(FDR < 0.05, "FDR < 0.05", "FDR ≥ 0.05")
  )

p2f <- ggplot(fisher_plot_data, aes(fold_enrichment, TF)) +
  geom_segment(aes(x = 1, xend = fold_enrichment, yend = TF),
               colour = "grey65", linewidth = 0.7) +
  geom_point(aes(size = overlap, colour = Enrichment)) +
  scale_x_log10() +
  scale_colour_manual(values = c("FDR < 0.05" = "#D73027", "FDR ≥ 0.05" = "#74ADD1")) +
  labs(x = "Fold enrichment in restoration core", y = NULL,
       title = "TF target-set enrichment", size = "Overlap") +
  theme_tf(11) +
  theme(legend.position = "right")

figure2 <- (p2a | p2b) / (p2c | p2d) / (p2e | p2f) +
  plot_layout(guides = "keep") +
  plot_annotation(tag_levels = "A") &
  panel_tag_theme

save_plot_pair(figure2, "Figure2_Discovery", 14, 13)

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

render_pdf_panel <- function(path, trim = TRUE) {
  png_path <- tempfile(fileext = ".png")
  pdftools::pdf_convert(path, format = "png", dpi = 600,
                        filenames = png_path, verbose = FALSE)
  image <- magick::image_read(png_path)
  if (trim) image <- magick::image_trim(image)
  image <- magick::image_background(image, "white")
  image <- magick::image_flatten(image)
  magick::image_border(image, "white", "40x40")
}

lung_tf <- read_tsv(out_path("results", "gse192483_candidate_tf_localization.tsv"),
                    show_col_types = FALSE) %>%
  filter(TF %in% c("PAX5", "TBX21", "NFATC2", "MAF", "STAT4"),
         n_patients >= 2) %>%
  mutate(TF = factor(TF, levels = rev(c("PAX5", "TBX21", "NFATC2", "MAF", "STAT4"))),
         CellType = factor(CellType, levels = sort(unique(CellType))))
p3a <- ggplot(lung_tf,
              aes(CellType, TF, size = median_percent_expressing,
                  colour = median_mean_expression)) +
  geom_point() +
  scale_size_continuous(range = c(0.5, 7)) +
  scale_colour_gradient(low = "#2C7FB8", high = "#D95F0E") +
  guides(colour = guide_colourbar(title.position = "top",
                                 barwidth = unit(4, "cm"), barheight = unit(0.35, "cm")),
         size = guide_legend(title.position = "top", nrow = 1)) +
  labs(x = "Lung cell type", y = NULL, size = "% expressing",
       colour = "Mean expression") +
  theme_tf(13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
        axis.text.y = element_text(size = 12, face = "bold"),
        legend.position = "bottom",
        legend.box = "horizontal",
        legend.spacing.x = unit(0.4, "cm"))
p3b_img <- render_pdf_panel(file.path(
  out_path("inputs", "figure_source_context"), "Fig3F_GSE192483_3TF_UMAP.pdf"
))
p3b <- ggdraw() + draw_image(p3b_img)
figure3 <- plot_grid(p3a, p3b, ncol = 1, labels = c("A", "B"),
                     label_size = 20, rel_heights = c(1, 0.75))
ggsave(file.path(figure_dir, "Figure3_SingleCellContext.pdf"),
       figure3, width = 15, height = 9.5, device = cairo_pdf, bg = "white")

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

save_plot_pair(figure4, "Figure4_OutcomeModels", 15, 10)

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

# Figure 7: donor-aware, cell-type-specific response to Mtb stimulation.
pb_summary <- pb_long %>%
  group_by(gene, cell_type, timepoint) %>%
  summarise(
    mean = mean(expression),
    se = sd(expression) / sqrt(n()),
    n_donors = n_distinct(assignment),
    .groups = "drop"
  ) %>%
  mutate(
    timepoint = factor(timepoint, levels = c("UT", "3hMTB", "24hMTB"))
  )

pb_effect <- pb_models %>%
  filter(
    (gene == "TBX21" & cell_type %in% c("NK", "CD8T", "CD4T", "B")) |
      (gene %in% c("PAX5", "EBF1") & cell_type == "B") |
      (gene %in% c("NFATC2", "IFNG", "CD40LG") &
         cell_type %in% c("NK", "CD8T", "CD4T"))
  ) %>%
  mutate(
    comparison = factor(
      contrast,
      levels = c("3hMTB - UT", "24hMTB - UT"),
      labels = c("3 h", "24 h")
    ),
    feature = paste(gene, cell_type, sep = " | "),
    significance = case_when(
      FDR < 0.001 ~ "***",
      FDR < 0.01 ~ "**",
      FDR < 0.05 ~ "*",
      TRUE ~ ""
    )
  )

p7a <- pb_effect %>%
  mutate(
    gene = factor(gene, levels = c("TBX21", "PAX5", "EBF1", "NFATC2", "IFNG", "CD40LG")),
    cell_time = factor(
      paste(cell_type, comparison, sep = " · "),
      levels = rev(as.vector(outer(c("B", "CD4T", "CD8T", "NK"),
                                    c("3 h", "24 h"), paste, sep = " · ")))
    )
  ) %>%
  ggplot(aes(gene, cell_time, fill = estimate)) +
  geom_tile(width = 0.92, height = 0.92, colour = "white", linewidth = 0.5) +
  geom_text(aes(label = significance), size = 4, fontface = "bold") +
  scale_fill_gradient2(
    low = "#4575B4", mid = "white", high = "#D73027",
    midpoint = 0, name = "Change versus unstimulated",
    breaks = c(-0.2, 0, 0.2)
  ) +
  guides(fill = guide_colorbar(title.position = "top",
                               barwidth = grid::unit(4, "cm"),
                               barheight = grid::unit(0.3, "cm"))) +
  labs(x = NULL, y = NULL, title = "Mtb response by cell state") +
  theme_tf(11) +
  coord_fixed(ratio = 1) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom",
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 10))

celltype_gene_plot <- function(gene_name, cell_types) {
  data <- pb_summary %>%
    filter(gene == gene_name, cell_type %in% cell_types)
  p <- ggplot(
    data,
    aes(timepoint, mean, colour = cell_type, group = cell_type)
  ) +
    geom_errorbar(
      aes(ymin = mean - 1.96 * se, ymax = mean + 1.96 * se),
      width = 0.12, linewidth = 0.7,
      position = position_dodge(width = 0.16)
    ) +
    geom_line(linewidth = 1.1, position = position_dodge(width = 0.16)) +
    geom_point(size = 2.7, position = position_dodge(width = 0.16)) +
    labs(x = NULL, y = "Pseudobulk expression", title = gene_name) +
    theme_tf(10) +
    theme(legend.position = "bottom", legend.title = element_blank())
  p
}

p7b <- celltype_gene_plot("TBX21", c("NK", "CD8T", "CD4T", "B"))
p7c <- celltype_gene_plot("PAX5", c("B"))
p7d <- celltype_gene_plot("NFATC2", c("NK", "CD8T", "CD4T"))
p7e <- celltype_gene_plot("IFNG", c("NK", "CD8T", "CD4T"))
p7f <- celltype_gene_plot("CD40LG", c("NK", "CD8T", "CD4T"))

figure7 <- (p7a | p7b | p7c) / (p7d | p7e | p7f) +
  plot_layout(heights = c(1.2, 1)) +
  plot_annotation(tag_levels = "A") &
  panel_tag_theme

save_plot_pair(figure7, "FigureS2_CellTypeSpecificMtbResponse", 15, 11)

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

# Supplementary figure: descriptive one-donor stimulation example.
pS2a_img <- render_pdf_panel(file.path(
  out_path("inputs", "figure_source_context"), "Fig3A_GSE218065_scRNA_Dotplot_annotation.pdf"
))
pS2b_img <- render_pdf_panel(file.path(
  out_path("inputs", "figure_source_context"), "Fig3B_GSE218065_PBMC_Mtb_auto_annotation_umap.pdf"
))
pS2c_img <- render_pdf_panel(file.path(
  out_path("inputs", "figure_source_context"), "Fig3C_GSE218065_6TF-scRNA-UMAP.pdf"
))

pS2a <- ggdraw() + draw_image(pS2a_img)
pS2b <- ggdraw() + draw_image(pS2b_img)
pS2c <- ggdraw() + draw_image(pS2c_img)

supp_figure2 <- plot_grid(
  plot_grid(pS2a, labels = "A", label_size = 20),
  plot_grid(pS2b, pS2c, ncol = 2, rel_widths = c(0.9, 1.1),
            labels = c("B", "C"), label_size = 20),
  ncol = 1, rel_heights = c(0.43, 0.57)
)
ggsave(file.path(figure_dir, "archive", "FigureS2_GSE218065_descriptive.pdf"),
       supp_figure2, width = 15, height = 11, device = cairo_pdf, bg = "white")

writeLines(c(
  "Main figures generated:",
  "Figure1_StudyDesign",
  "Figure2_Discovery",
  "Figure3_SingleCellContext",
  "Figure4_OutcomeModels",
  "Figure5_TBRI_OutcomeAssociation",
  "Figure6_CompositionConfounding",
  "Supplementary figures generated:",
  "FigureS1_METANet",
  "FigureS2_CellTypeSpecificMtbResponse"
), out_path("logs", "figure_generation.txt"))
