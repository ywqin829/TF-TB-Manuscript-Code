# Run from the project root: Rscript scripts/figure_sources/FigureS3_CellTypeSpecificMtbResponse.R
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

save_plot_pair(figure7, "FigureS3_CellTypeSpecificMtbResponse", 15, 11)

