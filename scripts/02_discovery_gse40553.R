source(file.path("scripts", "00_config.R"))

suppressPackageStartupMessages({
  library(limma)
  library(GSVA)
  library(decoupleR)
  library(dorothea)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
  library(clusterProfiler)
  library(org.Hs.eg.db)
})

set.seed(20261028)

analysis_inputs <- readRDS(out_path("results", "analysis_inputs.rds"))
obj <- analysis_inputs$gse40553
expr <- obj$expr_symbol
meta <- obj$pheno

meta$Group <- factor(
  meta$Group,
  levels = c("LTB", "TB_0M", "TB_0.5M", "TB_2M", "TB_6M", "TB_12M")
)
meta$SubjectID <- factor(meta$SubjectID)

design <- model.matrix(~ 0 + Group, data = meta)
colnames(design) <- levels(meta$Group)

dc_gene <- duplicateCorrelation(expr, design, block = meta$SubjectID)
fit_gene <- lmFit(
  expr,
  design,
  block = meta$SubjectID,
  correlation = dc_gene$consensus.correlation
)
contrast_gene <- makeContrasts(
  Disease = TB_0M - LTB,
  Treatment6M = TB_6M - TB_0M,
  levels = design
)
fit_gene_contrast <- eBayes(contrasts.fit(fit_gene, contrast_gene))

disease <- topTable(
  fit_gene_contrast, coef = "Disease", number = Inf, adjust.method = "BH", sort.by = "P"
)
treatment6m <- topTable(
  fit_gene_contrast, coef = "Treatment6M", number = Inf, adjust.method = "BH", sort.by = "P"
)

restoration_core <- intersect(
  rownames(disease)[disease$adj.P.Val < 0.05 & disease$logFC < -1],
  rownames(treatment6m)[treatment6m$adj.P.Val < 0.05 & treatment6m$logFC > 1]
)
disease_driver <- intersect(
  rownames(disease)[disease$adj.P.Val < 0.05 & disease$logFC > 1],
  rownames(treatment6m)[treatment6m$adj.P.Val < 0.05 & treatment6m$logFC < -1]
)

if (length(restoration_core) < 20 || length(disease_driver) < 20) {
  stop("Discovery module definitions failed; inspect differential-expression output.")
}

disease_out <- disease %>%
  rownames_to_column("Gene") %>%
  mutate(Module = case_when(
    Gene %in% restoration_core ~ "Restoration core",
    Gene %in% disease_driver ~ "Disease driver",
    TRUE ~ "Other"
  ))

treatment_out <- treatment6m %>%
  rownames_to_column("Gene") %>%
  mutate(Module = case_when(
    Gene %in% restoration_core ~ "Restoration core",
    Gene %in% disease_driver ~ "Disease driver",
    TRUE ~ "Other"
  ))

write_tsv(disease_out, out_path("results", "gse40553_disease_0m_vs_ltb.tsv"))
write_tsv(treatment_out, out_path("results", "gse40553_treatment_6m_vs_0m.tsv"))
writeLines(restoration_core, out_path("results", "restoration_core_genes.txt"))
writeLines(disease_driver, out_path("results", "disease_driver_genes.txt"))

module_list <- list(
  Restoration_Core = restoration_core,
  Disease_Driver = disease_driver
)

gsva_params <- ssgseaParam(expr, module_list, minSize = 5)
module_scores <- t(gsva(gsva_params, verbose = FALSE)) %>%
  as.data.frame() %>%
  rownames_to_column("SampleID") %>%
  left_join(
    meta %>% rownames_to_column("SampleID") %>%
      dplyr::select(SampleID, SubjectID, Group),
    by = "SampleID"
  )

write_tsv(module_scores, out_path("results", "gse40553_module_scores.tsv"))

data(dorothea_hs, package = "dorothea")
regulons_abc <- dorothea_hs %>%
  filter(confidence %in% c("A", "B", "C"))

tf_cache <- out_path("results", "gse40553_tf_activity_cache.rds")
if (file.exists(tf_cache)) {
  tf_activity <- readRDS(tf_cache)
} else {
  tf_long <- decoupleR::run_viper(
    mat = expr,
    network = regulons_abc,
    .source = "tf",
    .target = "target",
    .mor = "mor",
    minsize = 5,
    verbose = FALSE
  )

  tf_activity <- tf_long %>%
    filter(statistic == "viper") %>%
    pivot_wider(id_cols = source, names_from = condition, values_from = score) %>%
    column_to_rownames("source") %>%
    as.matrix()

  tf_activity <- tf_activity[, rownames(meta), drop = FALSE]
  tf_activity[is.na(tf_activity)] <- 0
  saveRDS(tf_activity, tf_cache)
}

dc_tf <- duplicateCorrelation(tf_activity, design, block = meta$SubjectID)
fit_tf <- lmFit(
  tf_activity,
  design,
  block = meta$SubjectID,
  correlation = dc_tf$consensus.correlation
)
fit_tf_contrast <- eBayes(contrasts.fit(fit_tf, contrast_gene))

tf_disease <- topTable(
  fit_tf_contrast, coef = "Disease", number = Inf, adjust.method = "BH", sort.by = "P"
)
tf_treatment <- topTable(
  fit_tf_contrast, coef = "Treatment6M", number = Inf, adjust.method = "BH", sort.by = "P"
)

regulons_all <- dorothea_hs %>%
  filter(confidence %in% c("A", "B", "C", "D")) %>%
  distinct(tf, target)

expressed_universe <- intersect(rownames(expr), unique(regulons_all$target))
recovery_in_universe <- intersect(restoration_core, expressed_universe)
fisher_rows <- lapply(split(regulons_all, regulons_all$tf), function(x) {
  targets <- intersect(unique(x$target), expressed_universe)
  overlap <- intersect(targets, recovery_in_universe)
  m <- length(targets)
  k <- length(overlap)
  if (m < 5 || k < 2) return(NULL)
  p <- phyper(k - 1, m, length(expressed_universe) - m,
              length(recovery_in_universe), lower.tail = FALSE)
  data.frame(
    TF = x$tf[1],
    overlap = k,
    regulon_size_in_background = m,
    p_value = p
  )
})

fisher_tbl <- bind_rows(fisher_rows) %>%
  mutate(
    FDR = p.adjust(p_value, method = "BH"),
    fold_enrichment = (overlap / length(recovery_in_universe)) /
      (regulon_size_in_background / length(expressed_universe))
  ) %>%
  arrange(p_value)

tf_treatment_df <- data.frame(
  TF = rownames(tf_treatment),
  treatment_logFC = tf_treatment$logFC,
  treatment_aveExpr = tf_treatment$AveExpr,
  treatment_t = tf_treatment$t,
  treatment_p = tf_treatment$P.Value,
  treatment_FDR = tf_treatment$adj.P.Val
)

tf_disease_df <- data.frame(
  TF = rownames(tf_disease),
  disease_logFC = tf_disease$logFC,
  disease_p = tf_disease$P.Value,
  disease_FDR = tf_disease$adj.P.Val
)

tf_activity_out <- tf_treatment_df %>%
  left_join(tf_disease_df, by = "TF") %>%
  left_join(fisher_tbl, by = "TF") %>%
  arrange(treatment_FDR, p_value) %>%
  mutate(
    significant_activity = treatment_FDR < 0.05 & abs(treatment_logFC) > 1,
    significant_target_enrichment = FDR < 0.05
  )

write_tsv(tf_activity_out, out_path("results", "gse40553_tf_activity_subject_blocked.tsv"))
write_tsv(fisher_tbl, out_path("results", "gse40553_tf_target_fisher.tsv"))

elite_drivers <- tf_activity_out %>%
  filter(significant_activity, significant_target_enrichment)
write_tsv(elite_drivers, out_path("results", "gse40553_elite_tf_candidates.tsv"))

paired_ids <- module_scores %>%
  dplyr::select(SubjectID, Group) %>%
  distinct() %>%
  count(SubjectID, name = "n_timepoints") %>%
  filter(n_timepoints >= 2) %>%
  pull(SubjectID)

paired_module <- module_scores %>%
  filter(SubjectID %in% paired_ids) %>%
  dplyr::select(SampleID, SubjectID, Group, Restoration_Core, Disease_Driver)

tf_paired <- tf_activity[, paired_module$SampleID, drop = FALSE] %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("SampleID") %>%
  left_join(paired_module, by = "SampleID")

donor_delta_cor <- bind_rows(lapply(
  unique(paired_module$SubjectID),
  function(id) {
    x <- tf_paired %>% filter(SubjectID == id)
    b <- x %>% filter(Group == "TB_0M")
    e <- x %>% filter(Group == "TB_6M")
    if (nrow(b) != 1 || nrow(e) != 1) return(NULL)
    tf_cols <- rownames(tf_activity)
    delta <- e[tf_cols] - b[tf_cols]
    data.frame(
      SubjectID = id,
      TF = tf_cols,
      delta_tf_activity = as.numeric(delta[1, ]),
      delta_restoration = e$Restoration_Core - b$Restoration_Core,
      delta_disease = e$Disease_Driver - b$Disease_Driver
    )
  }
))

donor_correlations <- donor_delta_cor %>%
  group_by(TF) %>%
  summarise(
    n_paired_donors = sum(complete.cases(delta_tf_activity,
                                         delta_restoration,
                                         delta_disease)),
    rho_restoration = suppressWarnings(cor(
      delta_tf_activity, delta_restoration,
      method = "spearman", use = "complete.obs"
    )),
    rho_disease = suppressWarnings(cor(
      delta_tf_activity, delta_disease,
      method = "spearman", use = "complete.obs"
    )),
    .groups = "drop"
  ) %>%
  arrange(desc(rho_restoration))

write_tsv(donor_delta_cor, out_path("results", "gse40553_paired_donor_tf_deltas.tsv"))
write_tsv(donor_correlations, out_path("results", "gse40553_paired_donor_correlations.tsv"))

go_table <- function(genes, label) {
  if (length(genes) < 10) return(tibble())
  ego <- enrichGO(
    gene = genes,
    OrgDb = org.Hs.eg.db,
    keyType = "SYMBOL",
    ont = "BP",
    universe = rownames(expr),
    pAdjustMethod = "BH",
    pvalueCutoff = 0.1,
    qvalueCutoff = 0.2,
    readable = TRUE
  )
  as.data.frame(ego) %>%
    slice_head(n = 15) %>%
    mutate(Module = label)
}

go_results <- bind_rows(
  go_table(restoration_core, "Restoration core"),
  go_table(disease_driver, "Disease driver")
)
write_tsv(go_results, out_path("results", "gse40553_module_go.tsv"))

saveRDS(
  list(
    design = design,
    gene_consensus_correlation = dc_gene$consensus.correlation,
    tf_consensus_correlation = dc_tf$consensus.correlation,
    disease = disease,
    treatment6m = treatment6m,
    restoration_core = restoration_core,
    disease_driver = disease_driver,
    module_scores = module_scores,
    tf_activity = tf_activity,
    tf_activity_table = tf_activity_out,
    fisher = fisher_tbl,
    elite_drivers = elite_drivers,
    donor_correlations = donor_correlations,
    go = go_results
  ),
  out_path("results", "gse40553_discovery_results.rds")
)

writeLines(c(
  paste0("Restoration core genes: ", length(restoration_core)),
  paste0("Disease driver genes: ", length(disease_driver)),
  paste0("Gene-level within-subject correlation: ", round(dc_gene$consensus.correlation, 4)),
  paste0("TF-level within-subject correlation: ", round(dc_tf$consensus.correlation, 4)),
  paste0("Elite TF candidates: ", nrow(elite_drivers))
), out_path("logs", "gse40553_discovery_summary.txt"))
