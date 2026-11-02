source(file.path("scripts", "00_config.R"))

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
  library(GenomicFeatures)
  library(TxDb.Hsapiens.UCSC.hg19.knownGene)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(org.Hs.eg.db)
  library(immunedeconv)
})

raw_path <- out_path(
  "data", "GSE89403_official",
  "GSE89403_rawCounts_GeneNames_AllSamples.csv.gz"
)
raw <- fread(raw_path)
counts <- as.matrix(raw[, -(1:2)])
storage.mode(counts) <- "double"
rownames(counts) <- raw$V1

make_tpm <- function(txdb) {
  gene_width <- sum(width(reduce(exonsBy(txdb, by = "gene"))))
  mapping <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = rownames(counts),
    keytype = "ENSEMBL",
    columns = "ENTREZID"
  ) %>%
    filter(!is.na(ENTREZID), !duplicated(ENSEMBL)) %>%
    mutate(gene_length = gene_width[ENTREZID]) %>%
    filter(!is.na(gene_length), gene_length > 0)

  genes <- intersect(mapping$ENSEMBL, rownames(counts))
  lengths <- gene_width[mapping$ENTREZID[match(genes, mapping$ENSEMBL)]]
  names(lengths) <- genes
  rate <- sweep(counts[genes, , drop = FALSE], 1, lengths, "/")
  tpm <- sweep(rate, 2, colSums(rate), "/") * 1e6
  symbols <- raw$symbol[match(rownames(tpm), raw$V1)]
  tpm_symbol <- rowsum(tpm, group = symbols, reorder = TRUE)
  tpm_symbol[rownames(tpm_symbol) != "", , drop = FALSE]
}

run_q <- function(tpm) {
  x <- immunedeconv::deconvolute(
    tpm, method = "quantiseq", tumor = FALSE, arrays = FALSE
  ) %>%
    as.data.frame()
  rownames(x) <- x$cell_type
  x$cell_type <- NULL
  as.matrix(x)
}

inputs <- readRDS(out_path("results", "analysis_inputs.rds"))
samples <- rownames(inputs$gse89403$metadata)
hg19 <- make_tpm(TxDb.Hsapiens.UCSC.hg19.knownGene)[, samples, drop = FALSE]
hg38 <- make_tpm(TxDb.Hsapiens.UCSC.hg38.knownGene)[, samples, drop = FALSE]
q19 <- run_q(hg19)
q38 <- run_q(hg38)

compare <- bind_rows(lapply(intersect(rownames(q19), rownames(q38)), function(cell) {
  data.frame(
    cell_type = cell,
    spearman_hg19_vs_hg38 = cor(
      q19[cell, samples], q38[cell, samples],
      method = "spearman", use = "pairwise.complete.obs"
    ),
    median_abs_fraction_difference = median(
      abs(q19[cell, samples] - q38[cell, samples]),
      na.rm = TRUE
    )
  )
}))

write_tsv(compare, out_path("results", "gse89403_hg19_hg38_quantiseq_comparison.tsv"))
write_tsv(
  as.data.frame(t(q38)) %>% rownames_to_column("SampleID"),
  out_path("results", "gse89403_quantiseq_hg38_sensitivity.tsv")
)

writeLines(c(
  "Raw GSE89403 gene counts were generated after STAR alignment to hg19.",
  "The primary TPM matrix uses hg19 gene lengths; hg38 is retained only as an annotation-mismatch sensitivity analysis.",
  paste0("HG19/HG38 cell types compared: ", nrow(compare)),
  paste0("Minimum Spearman correlation: ", round(min(compare$spearman_hg19_vs_hg38, na.rm = TRUE), 4)),
  paste0("Maximum median absolute fraction difference: ",
         signif(max(compare$median_abs_fraction_difference, na.rm = TRUE), 4))
), out_path("logs", "gse89403_hg19_hg38_sensitivity_summary.txt"))
