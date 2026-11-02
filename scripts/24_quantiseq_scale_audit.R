source(file.path("scripts", "00_config.R"))

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
  library(GenomicFeatures)
  library(TxDb.Hsapiens.UCSC.hg19.knownGene)
  library(org.Hs.eg.db)
  library(immunedeconv)
})

inputs <- readRDS(out_path("results", "analysis_inputs.rds"))
raw_counts_path <- out_path(
  "data", "GSE89403_official",
  "GSE89403_rawCounts_GeneNames_AllSamples.csv.gz"
)
if (!file.exists(raw_counts_path)) {
  stop("Missing GSE89403 official raw-count matrix: ", raw_counts_path)
}

matrix_to_nonlog <- function(expr, dataset) {
  if (any(!is.finite(expr))) {
    stop(dataset, " contains non-finite expression values.")
  }
  if (max(expr, na.rm = TRUE) < 50) {
    return(2^expr)
  }
  expr
}

run_quantiseq <- function(expr_nonlog, arrays, dataset) {
  result <- immunedeconv::deconvolute(
    expr_nonlog,
    method = "quantiseq",
    tumor = FALSE,
    arrays = arrays
  )
  result <- as.data.frame(result)
  rownames(result) <- result$cell_type
  result$cell_type <- NULL
  as.matrix(result)
}

comp_40553_explicit <- run_quantiseq(
  matrix_to_nonlog(inputs$gse40553$expr_symbol, "GSE40553"),
  arrays = TRUE,
  dataset = "GSE40553"
)
comp_31348_explicit <- run_quantiseq(
  matrix_to_nonlog(inputs$gse31348$expression, "GSE31348"),
  arrays = TRUE,
  dataset = "GSE31348"
)

raw <- fread(raw_counts_path)
if (!all(c("V1", "symbol") %in% names(raw))) {
  stop("Unexpected raw-count columns.")
}
counts <- as.matrix(raw[, -(1:2)])
storage.mode(counts) <- "double"
rownames(counts) <- raw$V1

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
exons_by_gene <- exonsBy(txdb, by = "gene")
entrez_lengths <- sum(width(reduce(exons_by_gene)))

ensembl_to_entrez <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = rownames(counts),
  keytype = "ENSEMBL",
  columns = "ENTREZID"
) %>%
  filter(!is.na(ENTREZID), !duplicated(ENSEMBL)) %>%
  mutate(gene_length = entrez_lengths[ENTREZID]) %>%
  filter(!is.na(gene_length), gene_length > 0)

valid_genes <- intersect(ensembl_to_entrez$ENSEMBL, rownames(counts))
lengths <- entrez_lengths[ensembl_to_entrez$ENTREZID[
  match(valid_genes, ensembl_to_entrez$ENSEMBL)
]]
names(lengths) <- valid_genes

counts_valid <- counts[valid_genes, , drop = FALSE]
rate <- sweep(counts_valid, 1, lengths, "/")
tpm <- sweep(rate, 2, colSums(rate), "/") * 1e6

meta_89403 <- inputs$gse89403$metadata
selected_samples <- intersect(rownames(meta_89403), colnames(tpm))
if (length(selected_samples) != nrow(meta_89403)) {
  stop("TPM matrix does not contain all selected GSE89403 samples.")
}
tpm <- tpm[, selected_samples, drop = FALSE]

# Collapse duplicate probes/transcripts consistently by symbol after TPM.
symbol_for_gene <- raw$symbol[match(rownames(tpm), raw$V1)]
tpm_symbol <- rowsum(tpm, group = symbol_for_gene, reorder = TRUE)
tpm_symbol <- tpm_symbol[rownames(tpm_symbol) != "", , drop = FALSE]

comp_89403_tpm <- run_quantiseq(
  tpm_symbol,
  arrays = FALSE,
  dataset = "GSE89403"
)

load_old_comp <- function(stem) {
  read.delim(
    out_path("results", paste0(stem, "_quantiseq.tsv")),
    check.names = FALSE
  ) %>%
    column_to_rownames("SampleID") %>%
    as.matrix() %>%
    t()
}

comp_40553_old <- load_old_comp("gse40553")
comp_31348_old <- load_old_comp("gse31348")
comp_89403_old <- load_old_comp("gse89403")

compare_comp <- function(old, new, dataset) {
  common_samples <- intersect(colnames(old), colnames(new))
  common_cells <- intersect(rownames(old), rownames(new))
  bind_rows(lapply(common_cells, function(cell) {
    rho <- cor(
      old[cell, common_samples],
      new[cell, common_samples],
      method = "spearman",
      use = "pairwise.complete.obs"
    )
    data.frame(
      dataset = dataset,
      cell_type = cell,
      spearman_old_vs_nonlog = rho,
      median_abs_difference = median(
        abs(old[cell, common_samples] - new[cell, common_samples]),
        na.rm = TRUE
      )
    )
  }))
}

comparison <- bind_rows(
  compare_comp(comp_40553_old, comp_40553_explicit, "GSE40553"),
  compare_comp(comp_31348_old, comp_31348_explicit, "GSE31348"),
  compare_comp(comp_89403_old, comp_89403_tpm, "GSE89403")
)

write_tsv(
  comparison,
  out_path("results", "quantiseq_scale_correction_comparison.tsv")
)
write_tsv(
  as.data.frame(t(comp_40553_explicit)) %>% rownames_to_column("SampleID"),
  out_path("results", "gse40553_quantiseq_explicit_nonlog.tsv")
)
write_tsv(
  as.data.frame(t(comp_31348_explicit)) %>% rownames_to_column("SampleID"),
  out_path("results", "gse31348_quantiseq_explicit_nonlog.tsv")
)
write_tsv(
  as.data.frame(t(comp_89403_tpm)) %>% rownames_to_column("SampleID"),
  out_path("results", "gse89403_quantiseq_tpm.tsv")
)

writeLines(c(
  "quanTIseq was run explicitly on non-log expression.",
  "Microarray matrices were transformed with 2^x before deconvolution.",
  "GSE89403 TPM was recomputed from official hg19-aligned raw counts and hg19 gene lengths.",
  paste0("GSE89403 samples aligned: ", ncol(tpm_symbol)),
  paste0("GSE89403 genes with valid Ensembl-to-Entrez length: ", length(valid_genes)),
  "The comparison table reports old-vs-corrected Spearman correlation and median absolute fraction difference for each cell type."
), out_path("logs", "quantiseq_scale_audit_summary.txt"))
