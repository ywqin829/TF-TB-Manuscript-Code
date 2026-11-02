#!/usr/bin/env Rscript
# Publish the fixed discovery gene contrasts with independent validation results.
source("scripts/00_config.R")
suppressPackageStartupMessages({library(dplyr);library(readr);library(tibble)})
x <- readRDS(out_path("results","gse40553_discovery_results.rds"))
to_gene <- function(z,prefix) {
  z <- as.data.frame(z)
  z$Gene <- rownames(z)
  as_tibble(z) %>% rename_with(~paste0(prefix,"_",.x),-Gene)
}
discovery <- full_join(to_gene(x$disease,"Disease"),to_gene(x$treatment6m,"Treatment6m"),by="Gene") %>%
  mutate(RestorationCore=Gene %in% x$restoration_core,
         DiseaseDriver=Gene %in% x$disease_driver) %>%
  filter(RestorationCore | DiseaseDriver)
write_tsv(discovery,out_path("results","fixed_core_genes_discovery_evidence.tsv"))
