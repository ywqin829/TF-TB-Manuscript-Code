# Run from the project root: Rscript scripts/figure_sources/Figure2_Discovery.R
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

