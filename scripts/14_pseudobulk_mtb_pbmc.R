source(file.path("scripts", "00_config.R"))

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
  library(lme4)
  library(lmerTest)
  library(emmeans)
})

set.seed(20261028)

mtb_path <- "D:/RDS/1M-scBloodNL/Mtb/mtb_final.rds"
if (!file.exists(mtb_path)) stop("Missing local single-cell object: ", mtb_path)

obj <- readRDS(mtb_path)
genes <- intersect(
  c("TBX21", "PAX5", "EBF1", "NFATC2", "IFNG", "CD40LG",
    "CD19", "CD79A", "RUNX3", "TCF7", "PRDM1", "XBP1"),
  rownames(obj)
)

meta <- obj@meta.data[, c("assignment", "timepoint", "cell_type_lowerres")]
meta$assignment <- factor(meta$assignment)
meta$timepoint <- factor(meta$timepoint, levels = c("UT", "3hMTB", "24hMTB"))
meta$cell_type_lowerres <- factor(meta$cell_type_lowerres)

expr <- LayerData(obj, assay = "RNA", layer = "data")[genes, , drop = FALSE]
expr_dense <- as.matrix(expr)
group_id <- paste(meta$assignment, meta$timepoint, meta$cell_type_lowerres, sep = "|")
group_levels <- unique(group_id)
group_n <- tabulate(match(group_id, group_levels), nbins = length(group_levels))
names(group_n) <- group_levels
sums <- rowsum(t(expr_dense), group = group_id, reorder = TRUE)
pseudobulk <- sums / group_n[rownames(sums)]

parts <- strsplit(rownames(pseudobulk), "\\|", fixed = FALSE)
pb_meta <- data.frame(
  assignment = factor(vapply(parts, `[`, character(1), 1)),
  timepoint = factor(vapply(parts, `[`, character(1), 2),
                     levels = c("UT", "3hMTB", "24hMTB")),
  cell_type = factor(vapply(parts, `[`, character(1), 3)),
  n_cells = as.numeric(group_n[rownames(pseudobulk)]),
  stringsAsFactors = FALSE
)

pseudobulk_long <- as.data.frame(pseudobulk) %>%
  rownames_to_column("group_id") %>%
  bind_cols(pb_meta) %>%
  filter(n_cells >= 20, cell_type %in% c("NK", "CD4T", "CD8T", "B")) %>%
  pivot_longer(all_of(genes), names_to = "gene", values_to = "expression")

write_tsv(
  pseudobulk_long,
  out_path("results", "mtb_pbmc_pseudobulk_long.tsv")
)

fit_gene <- function(x, gene_name, cell_type_name) {
  x <- x %>%
    filter(complete.cases(expression, assignment, timepoint)) %>%
    droplevels()
  if (n_distinct(x$timepoint) < 2 || n_distinct(x$assignment) < 4) {
    return(NULL)
  }
  fit <- lmer(
    expression ~ timepoint + (1 | assignment),
    data = x,
    REML = TRUE,
    control = lmerControl(optimizer = "bobyqa")
  )
  emm <- emmeans(fit, ~ timepoint)
  con <- as.data.frame(emmeans::contrast(
    emm, method = "trt.vs.ctrl", ref = 1, adjust = "none"
  ))
  con$gene <- gene_name
  con$cell_type <- cell_type_name
  con$n_observations <- nrow(x)
  con$n_donors <- n_distinct(x$assignment)
  con
}

keys <- pseudobulk_long %>% distinct(gene, cell_type)
models <- bind_rows(lapply(seq_len(nrow(keys)), function(i) {
  gene_name <- keys$gene[i]
  cell_type_name <- keys$cell_type[i]
  fit_gene(
    filter(pseudobulk_long, gene == gene_name, cell_type == cell_type_name),
    gene_name,
    cell_type_name
  )
})) %>%
  group_by(gene) %>%
  mutate(FDR = p.adjust(p.value, method = "BH")) %>%
  ungroup() %>%
  arrange(gene, cell_type, contrast)

write_tsv(
  models,
  out_path("results", "mtb_pbmc_celltype_time_models.tsv")
)

# Percent positive cells for descriptive support.
positive_long <- lapply(genes, function(gene_name) {
  positive <- as.numeric(expr_dense[gene_name, ] > 0)
  sums <- rowsum(positive, group = group_id, reorder = TRUE)
  values <- as.numeric(sums[, 1]) / group_n[rownames(sums)]
  data.frame(
    group_id = rownames(sums),
    gene = gene_name,
    pct_positive = values,
    stringsAsFactors = FALSE
  )
}) %>%
  bind_rows() %>%
  left_join(
    data.frame(group_id = rownames(pseudobulk), pb_meta,
               stringsAsFactors = FALSE),
    by = "group_id"
  ) %>%
  filter(n_cells >= 20, cell_type %in% c("NK", "CD4T", "CD8T", "B"))

write_tsv(
  positive_long,
  out_path("results", "mtb_pbmc_celltype_percent_positive.tsv")
)

writeLines(c(
  paste0("Cells analysed: ", ncol(obj)),
  paste0("Donors detected: ", n_distinct(meta$assignment)),
  paste0("Cell types analysed: ", paste(levels(droplevels(meta$cell_type_lowerres)), collapse = ", ")),
  paste0("Genes analysed: ", paste(genes, collapse = ", ")),
  "Pseudobulk groups required at least 20 cells.",
  "Mixed models included donor random intercepts."
), out_path("logs", "mtb_pseudobulk_summary.txt"))
