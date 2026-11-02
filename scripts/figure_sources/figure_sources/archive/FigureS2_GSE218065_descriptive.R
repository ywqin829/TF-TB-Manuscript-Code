# Run from the project root: Rscript scripts/figure_sources/FigureS2_GSE218065_descriptive.R
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
ggsave(file.path(figure_dir, "FigureS2_GSE218065_descriptive.pdf"),
       supp_figure2, width = 15, height = 11, device = cairo_pdf, bg = "white")

