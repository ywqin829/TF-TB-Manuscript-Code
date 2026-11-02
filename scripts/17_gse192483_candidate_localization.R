source(file.path("scripts", "00_config.R"))

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
  library(ggplot2)
})

seurat_path <- Sys.getenv("GSE192483_RDS", unset = "")
if (!nzchar(seurat_path)) {
  fallback <- "data/GSE192483_3-Celltype/sce_celltype.rds"
  if (file.exists(fallback)) seurat_path <- fallback
}
if (!file.exists(seurat_path)) {
  stop(
    "Set GSE192483_RDS to the GSE192483 Seurat object. ",
    "Current value: ", seurat_path
  )
}

candidates <- c(
  "MAF", "SNAPC4", "TBX21", "ARNTL", "STAT4",
  "MEF2B", "RUNX3", "ZBED1", "NFATC2", "PAX5"
)

obj <- readRDS(seurat_path)
available <- intersect(candidates, rownames(obj))
if (length(available) == 0) {
  stop("None of the candidate TFs are present in the GSE192483 object.")
}

expr <- GetAssayData(obj, assay = "RNA", layer = "data")[
  available, , drop = FALSE
]
meta <- obj@meta.data[, c(
  "Patient", "Sample_title", "Tissue_Type", "celltype"
)]
meta$CellID <- rownames(meta)

sample_key <- interaction(
  meta$Patient,
  meta$Sample_title,
  meta$Tissue_Type,
  meta$celltype,
  drop = TRUE
)
cell_groups <- split(seq_len(nrow(meta)), sample_key)

sample_level <- bind_rows(lapply(cell_groups, function(idx) {
  if (length(idx) < 20) return(NULL)
  first <- meta[idx[1], , drop = FALSE]
  bind_rows(lapply(available, function(tf) {
    values <- as.numeric(expr[tf, idx])
    data.frame(
      Patient = first$Patient,
      Sample = first$Sample_title,
      Tissue_Type = first$Tissue_Type,
      CellType = first$celltype,
      TF = tf,
      n_cells = length(idx),
      mean_log_normalized_expression = mean(values),
      percent_expressing = 100 * mean(values > 0)
    )
  }))
}))

celltype_summary <- sample_level %>%
  group_by(TF, CellType) %>%
  summarise(
    n_patients = n_distinct(Patient),
    n_patient_tissue_groups = n(),
    median_mean_expression = median(mean_log_normalized_expression),
    q1_mean_expression = quantile(mean_log_normalized_expression, 0.25),
    q3_mean_expression = quantile(mean_log_normalized_expression, 0.75),
    median_percent_expressing = median(percent_expressing),
    .groups = "drop"
  ) %>%
  arrange(TF, desc(median_mean_expression))

top_compartment <- celltype_summary %>%
  group_by(TF) %>%
  slice_max(median_mean_expression, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(TF, top_in_vivo_compartment = CellType)

summary_table <- celltype_summary %>%
  left_join(top_compartment, by = "TF") %>%
  arrange(TF, desc(median_mean_expression))

write_tsv(
  sample_level,
  out_path("results", "gse192483_candidate_tf_patient_celltype.tsv")
)
write_tsv(
  summary_table,
  out_path("results", "gse192483_candidate_tf_localization.tsv")
)

plot_data <- celltype_summary %>%
  filter(n_patients >= 2) %>%
  mutate(
    TF = factor(TF, levels = rev(available)),
    CellType = factor(CellType, levels = sort(unique(CellType)))
  )

p <- ggplot(
  plot_data,
  aes(
    x = CellType,
    y = TF,
    size = median_percent_expressing,
    color = median_mean_expression
  )
) +
  geom_point() +
  scale_size_continuous(range = c(0.3, 6)) +
  scale_color_viridis_c(option = "C") +
  labs(
    x = NULL,
    y = NULL,
    size = "% expressing",
    color = "Mean expression"
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.text.y = element_text(face = "bold", color = "black"),
    legend.position = "right"
  )

dir.create(out_path("figures", "archive"), recursive = TRUE, showWarnings = FALSE)
ggsave(out_path("figures", "archive", "Figure7_GSE192483_Candidate_TF_Localization.pdf"),
       plot = p, width = 11, height = 6.8, units = "in", device = cairo_pdf)

writeLines(c(
  paste0("Input Seurat object: ", seurat_path),
  paste0("Candidate TFs present: ", paste(available, collapse = ", ")),
  "Donor/tissue/cell-type groups required at least 20 cells.",
  "Summary statistics are medians across patient-tissue groups; cells are not treated as replicates.",
  "GSE192483 is in vivo cross-sectional lung tissue, not a longitudinal treatment cohort."
), out_path("logs", "gse192483_candidate_localization_summary.txt"))
