source(file.path("scripts", "00_config.R"))

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(patchwork)
})

pseudo <- read_tsv(
  out_path("results", "gse326212_candidate_tf_pseudobulk.tsv"),
  show_col_types = FALSE
)
models <- read_tsv(
  out_path("results", "gse326212_candidate_tf_models.tsv"),
  show_col_types = FALSE
)
lopo <- read_tsv(
  out_path("results", "gse326212_stat4_leave_one_project_out.tsv"),
  show_col_types = FALSE
)

tf_order <- c(
  "MAF", "SNAPC4", "TBX21", "STAT4", "MEF2B",
  "RUNX3", "ZBED1", "NFATC2", "PAX5"
)
celltype_order <- c(
  "Alveolar_mac", "Mo-derived_AM", "Monocyte", "Neutrophil",
  "cDC1", "CD4_T", "CD8_T", "gdT_NK_ILC", "Treg", "B"
)

dot_data <- pseudo %>%
  filter(
    n_cells >= 20,
    TF %in% tf_order,
    celltype_coarse %in% celltype_order
  ) %>%
  group_by(TF, celltype_coarse) %>%
  summarise(
    median_expression = median(mean_expression),
    median_percent_expressing = median(percent_expressing),
    n_patient_samples = n(),
    .groups = "drop"
  ) %>%
  group_by(TF) %>%
  mutate(
    relative_expression = as.numeric(scale(median_expression))
  ) %>%
  ungroup() %>%
  mutate(
    TF = factor(TF, levels = rev(tf_order)),
    celltype_coarse = factor(celltype_coarse, levels = celltype_order)
  )

p7a <- ggplot(
  dot_data,
  aes(
    x = celltype_coarse,
    y = TF,
    size = median_percent_expressing,
    colour = relative_expression
  )
) +
  geom_point() +
  scale_size_continuous(range = c(0.3, 6.5), name = "% expressing") +
  scale_colour_gradient2(
    low = "#2166AC",
    mid = "#F7F7F7",
    high = "#B2182B",
    midpoint = 0,
    name = "Relative\nexpression"
  ) +
  labs(
    x = NULL,
    y = NULL,
    title = "In vivo BAL cell-type context"
  ) +
  theme_tf(11) +
  coord_fixed(ratio = 1) +
  theme(
    axis.text.x = element_text(angle = 55, hjust = 1, size = 10),
    panel.grid.major = element_blank(),
    legend.position = "bottom",
    legend.box = "vertical"
  )

effect_data <- models %>%
  filter(
    comparison == "PTB_vs_IGRAneg",
    TF %in% tf_order,
    celltype_coarse %in% celltype_order,
    !is.na(p_value)
  ) %>%
  mutate(
    t_effect = estimate / std_error,
    significance = case_when(
      FDR < 0.001 ~ "***",
      FDR < 0.01 ~ "**",
      FDR < 0.05 ~ "*",
      TRUE ~ ""
    ),
    TF = factor(TF, levels = rev(tf_order)),
    celltype_coarse = factor(celltype_coarse, levels = celltype_order)
  )

p7b <- ggplot(
  effect_data,
  aes(celltype_coarse, TF, fill = t_effect)
) +
  geom_tile(colour = "white", linewidth = 0.8) +
  geom_text(aes(label = significance), size = 3.8, fontface = "bold") +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-4, 4),
    oob = scales::squish,
    name = "PTB vs IGRA-\n(SE units)"
  ) +
  labs(
    x = NULL,
    y = NULL,
    title = "PTB vs IGRA-negative contacts"
  ) +
  theme_tf(11) +
  coord_fixed(ratio = 1) +
  theme(
    axis.text.x = element_text(angle = 55, hjust = 1, size = 10),
    legend.position = "bottom"
  )

full_stat4_estimate <- models %>%
  filter(comparison == "PTB_vs_IGRAneg", TF == "STAT4", celltype_coarse == "Mo-derived_AM") %>%
  pull(estimate)
stopifnot(length(full_stat4_estimate) == 1L, is.finite(full_stat4_estimate))

p7c <- ggplot(
  lopo,
  aes(x = estimate, y = reorder(excluded_project, estimate))
) +
  geom_vline(xintercept = full_stat4_estimate, linetype = 2, colour = "grey40") +
  geom_point(colour = "#B2182B", size = 2.8) +
  scale_x_continuous(limits = c(0.0215, 0.0245), breaks = c(0.022, 0.023, 0.024), labels = scales::number_format(accuracy = 0.001)) +
  labs(
    x = "Effect estimate",
    y = NULL,
    title = "STAT4: project sensitivity",
    subtitle = "Dashed line: all-project estimate"
  ) +
  theme_tf(9) +
  theme(
    axis.text.y = element_text(size = 7),
    plot.title = element_text(size = 11, hjust = 0),
    plot.subtitle = element_text(size = 8)
  )

figure <- (p7a | p7b) +
  plot_layout(widths = c(1, 1)) +
  plot_annotation(tag_levels = "A")

save_plot_pair(figure, "Figure7_InVivo_Candidate_TF_Context", 14, 7.7)

writeLines(c(
  "Figure 7 uses patient-sample-level GSE326212 BAL pseudobulk; cells are not replicates.",
  "Panel A shows median expression and detection across patient samples by coarse cell type; relative expression is scaled within each TF.",
  "Panel B shows PTB vs IGRA-negative sample-level model effects in standard-error units; stars are BH-adjusted P values.",
  "Leave-one-project-out estimates for STAT4 remain in results/gse326212_stat4_leave_one_project_out.tsv and are summarized in the Results, not displayed as a main-figure panel.",
  "Covariates: sex and smoking. Project/batch remains a limitation.",
  "This is in vivo disease-state evidence, not a longitudinal treatment analysis."
), out_path("logs", "gse326212_invivo_figure_summary.txt"))
