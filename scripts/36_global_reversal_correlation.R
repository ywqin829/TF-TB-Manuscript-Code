source(file.path("scripts", "00_config.R"))

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
})

disease <- read_tsv(
  out_path("results", "gse40553_disease_0m_vs_ltb.tsv"),
  show_col_types = FALSE
) %>%
  select(Gene, disease_logFC = logFC)

treatment <- read_tsv(
  out_path("results", "gse40553_treatment_6m_vs_0m.tsv"),
  show_col_types = FALSE
) %>%
  select(Gene, treatment_logFC = logFC)

paired <- inner_join(disease, treatment, by = "Gene")
pearson <- cor.test(
  paired$disease_logFC,
  paired$treatment_logFC,
  method = "pearson"
)
spearman <- cor.test(
  paired$disease_logFC,
  paired$treatment_logFC,
  method = "spearman",
  exact = FALSE
)

result <- tibble(
  Comparison = "GSE40553 active TB vs LTB logFC versus Month 6 vs baseline logFC",
  N_genes = nrow(paired),
  Pearson_r = unname(pearson$estimate),
  Pearson_P = pearson$p.value,
  Spearman_rho = unname(spearman$estimate),
  Spearman_P = spearman$p.value,
  Interpretation = "Negative genome-wide coupling supports bidirectional treatment reversal"
)

write_tsv(
  result,
  out_path("results", "gse40553_global_reversal_correlation.tsv")
)
writeLines(c(
  paste0("Genes compared: ", nrow(paired)),
  paste0("Pearson r: ", signif(unname(pearson$estimate), 4)),
  paste0("Spearman rho: ", signif(unname(spearman$estimate), 4)),
  "This is the genome-wide disease/treatment logFC correlation.",
  "It is distinct from the duplicateCorrelation consensus parameter of 0.075."
), out_path("logs", "global_reversal_correlation_summary.txt"))
