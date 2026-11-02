source(file.path("scripts", "00_config.R"))

suppressPackageStartupMessages({
  library(limma)
  library(GSVA)
  library(lme4)
  library(lmerTest)
  library(emmeans)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
  library(clusterProfiler)
  library(org.Hs.eg.db)
})

inputs <- readRDS(out_path("results", "analysis_inputs.rds"))
discovery <- readRDS(out_path("results", "gse40553_discovery_results.rds"))
obj <- inputs$gse40553
expr <- obj$expr_symbol
meta <- obj$pheno

meta$Group <- factor(
  meta$Group,
  levels = c("LTB", "TB_0M", "TB_0.5M", "TB_2M", "TB_6M", "TB_12M")
)
meta$SubjectID <- factor(meta$SubjectID)
region_column <- "geographical location:ch1"
south_africa <- grepl("South Africa", meta[[region_column]], fixed = TRUE)
meta_sa <- meta[south_africa, , drop = FALSE]
expr_sa <- expr[, rownames(meta_sa), drop = FALSE]

design <- model.matrix(~ 0 + Group, data = meta_sa)
colnames(design) <- levels(meta_sa$Group)
dc <- duplicateCorrelation(expr_sa, design, block = meta_sa$SubjectID)
fit <- lmFit(
  expr_sa,
  design,
  block = meta_sa$SubjectID,
  correlation = dc$consensus.correlation
)
contrast <- makeContrasts(
  Disease = TB_0M - LTB,
  Treatment6M = TB_6M - TB_0M,
  levels = design
)
fit <- eBayes(contrasts.fit(fit, contrast))
disease_sa <- topTable(fit, coef = "Disease", number = Inf, adjust.method = "BH")
treatment_sa <- topTable(fit, coef = "Treatment6M", number = Inf, adjust.method = "BH")

restoration_sa <- intersect(
  rownames(disease_sa)[disease_sa$adj.P.Val < 0.05 & disease_sa$logFC < -1],
  rownames(treatment_sa)[treatment_sa$adj.P.Val < 0.05 & treatment_sa$logFC > 1]
)
disease_driver_sa <- intersect(
  rownames(disease_sa)[disease_sa$adj.P.Val < 0.05 & disease_sa$logFC > 1],
  rownames(treatment_sa)[treatment_sa$adj.P.Val < 0.05 & treatment_sa$logFC < -1]
)

tf_sa <- discovery$tf_activity[, rownames(meta_sa), drop = FALSE]
tf_dc <- duplicateCorrelation(tf_sa, design, block = meta_sa$SubjectID)
tf_fit <- lmFit(
  tf_sa,
  design,
  block = meta_sa$SubjectID,
  correlation = tf_dc$consensus.correlation
)
tf_fit <- eBayes(contrasts.fit(tf_fit, contrast))
candidate_tfs <- c("TBX21", "PAX5", "MAF", "NFATC2")
candidate_tf_contrasts <- bind_rows(lapply(candidate_tfs, function(tf) {
  disease_row <- topTable(
    tf_fit,
    coef = "Disease",
    number = Inf,
    sort.by = "none"
  )[tf, , drop = FALSE]
  treatment_row <- topTable(
    tf_fit,
    coef = "Treatment6M",
    number = Inf,
    sort.by = "none"
  )[tf, , drop = FALSE]
  data.frame(
    TF = tf,
    disease_effect_SA = disease_row$logFC,
    disease_FDR_SA = disease_row$adj.P.Val,
    treatment_effect_SA = treatment_row$logFC,
    treatment_FDR_SA = treatment_row$adj.P.Val
  )
}))

overlap_summary <- bind_rows(
  data.frame(
    Module = "Restoration_Core",
    Original_genes = length(discovery$restoration_core),
    SA_only_genes = length(restoration_sa),
    Overlap = length(intersect(discovery$restoration_core, restoration_sa)),
    Percent_original_retained = 100 *
      length(intersect(discovery$restoration_core, restoration_sa)) /
      length(discovery$restoration_core),
    Jaccard = length(intersect(discovery$restoration_core, restoration_sa)) /
      length(union(discovery$restoration_core, restoration_sa))
  ),
  data.frame(
    Module = "Disease_Driver",
    Original_genes = length(discovery$disease_driver),
    SA_only_genes = length(disease_driver_sa),
    Overlap = length(intersect(discovery$disease_driver, disease_driver_sa)),
    Percent_original_retained = 100 *
      length(intersect(discovery$disease_driver, disease_driver_sa)) /
      length(discovery$disease_driver),
    Jaccard = length(intersect(discovery$disease_driver, disease_driver_sa)) /
      length(union(discovery$disease_driver, disease_driver_sa))
  )
)

module_list <- list(
  Original_Restoration = discovery$restoration_core,
  SouthAfrica_Restoration = restoration_sa
)
module_list <- lapply(module_list, intersect, rownames(inputs$gse31348$expression))
module_list <- module_list[lengths(module_list) >= 5]

if (length(module_list) > 0) {
  gsva_scores <- t(gsva(
    ssgseaParam(inputs$gse31348$expression, module_list, minSize = 5),
    verbose = FALSE
  )) %>%
    as.data.frame() %>%
    rownames_to_column("SampleID")

  meta_31348 <- inputs$gse31348$metadata
  meta_31348$SampleID <- rownames(meta_31348)
  meta_31348$TimePoint <- factor(
    meta_31348$TimePoint,
    levels = c("Week 0", "Week 1", "Week 2", "Week 4", "Week 26")
  )
  score_data <- gsva_scores %>% left_join(meta_31348, by = "SampleID")
  write_tsv(
    score_data,
    out_path("results", "gse40553_south_africa_gse31348_module_scores.tsv")
  )
  primary_development <- readRDS(
    out_path("results", "gse31348_development_results.rds")
  )$module_scores
  score_scale_match <- score_data %>%
    dplyr::select(SampleID, Original_Restoration) %>%
    inner_join(
      primary_development %>% dplyr::select(SampleID, Restoration_Core),
      by = "SampleID"
    )
  scale_factor <- mean(
    score_scale_match$Original_Restoration /
      score_scale_match$Restoration_Core
  )
  score_spearman <- cor(
    score_scale_match$Original_Restoration,
    score_scale_match$Restoration_Core,
    method = "spearman"
  )
  score_models <- bind_rows(lapply(
    setdiff(names(gsva_scores), "SampleID"),
    function(feature) {
      fit <- lmer(
        as.formula(paste0(feature, " ~ TimePoint + (1 | SubjectID)")),
        data = score_data
      )
      emm <- emmeans(fit, ~ TimePoint)
      tab <- as.data.frame(contrast(
        emm, method = "trt.vs.ctrl", ref = 1, adjust = "none"
      )) %>%
        filter(contrast == "Week 26 - Week 0") %>%
        mutate(
          module = feature,
          overall_time_p = anova(fit)$`Pr(>F)`[1],
          gsva_scale_factor_vs_primary = scale_factor,
          gsva_spearman_vs_primary = score_spearman
        )
      tab
    }
  ))
} else {
  score_models <- data.frame()
}

go_table <- function(genes, module_label) {
  if (length(genes) < 10) return(tibble())
  enrichGO(
    gene = genes,
    universe = rownames(expr_sa),
    OrgDb = org.Hs.eg.db,
    keyType = "SYMBOL",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.1,
    qvalueCutoff = 0.2,
    readable = TRUE
  ) %>%
    as.data.frame() %>%
    slice_head(n = 15) %>%
    mutate(Module = module_label)
}

go_restoration_sa <- go_table(
  restoration_sa,
  "South Africa restoration core"
)
go_driver_sa <- go_table(
  disease_driver_sa,
  "South Africa disease driver"
)
go_union <- enrichGO(
  gene = unique(c(restoration_sa, disease_driver_sa)),
  universe = rownames(expr_sa),
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pAdjustMethod = "BH",
  qvalueCutoff = 0.2,
  readable = TRUE
) %>%
  as.data.frame()

write_tsv(
  disease_sa %>% rownames_to_column("Gene") %>% mutate(Module = "SA_disease"),
  out_path("results", "gse40553_south_africa_disease.tsv")
)
write_tsv(
  treatment_sa %>% rownames_to_column("Gene") %>% mutate(Module = "SA_treatment"),
  out_path("results", "gse40553_south_africa_treatment.tsv")
)
write.table(
  go_union,
  out_path("results", "gse40553_south_africa_module_go.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write_tsv(
  go_restoration_sa,
  out_path("results", "gse40553_south_africa_restoration_go.tsv")
)
write_tsv(
  go_driver_sa,
  out_path("results", "gse40553_south_africa_disease_driver_go.tsv")
)
writeLines(
  restoration_sa,
  out_path("results", "gse40553_south_africa_restoration_core.txt")
)
writeLines(
  disease_driver_sa,
  out_path("results", "gse40553_south_africa_disease_driver.txt")
)
write_tsv(
  overlap_summary,
  out_path("results", "gse40553_sa_vs_original_module_overlap.tsv")
)
write_tsv(
  score_models,
  out_path("results", "gse40553_sa_core_gse31348_models.tsv")
)
write_tsv(
  candidate_tf_contrasts,
  out_path("results", "gse40553_south_africa_candidate_tf_contrasts.tsv")
)

writeLines(c(
  paste0("South Africa samples used: ", nrow(meta_sa)),
  paste0("South Africa active-TB baseline: ", sum(meta_sa$Group == "TB_0M")),
  paste0("South Africa LTB: ", sum(meta_sa$Group == "LTB")),
  paste0("Original restoration core: ", length(discovery$restoration_core)),
  paste0("South-Africa-only restoration core: ", length(restoration_sa)),
  paste0("Overlap: ", overlap_summary$Overlap[1]),
  paste0("GSE31348 SA-core Week26 effect available: ", nrow(score_models)),
  paste0(
    "GSE31348 SA-core Week26 estimate: ",
    signif(
      score_models$estimate[score_models$module == "SouthAfrica_Restoration"],
      4
    )
  ),
  paste0(
    "GSE31348 SA-core Week26 P: ",
    signif(
      score_models$p.value[score_models$module == "SouthAfrica_Restoration"],
      4
    )
  ),
  paste0("SA restoration GO terms (top 15): ", nrow(go_restoration_sa)),
  paste0("Top SA restoration GO term: ",
         if (nrow(go_restoration_sa) > 0) go_restoration_sa$Description[1] else "none"),
  "UK baseline samples were excluded from this sensitivity analysis."
), out_path("logs", "gse40553_south_africa_sensitivity_summary.txt"))
