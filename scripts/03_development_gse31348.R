source(file.path("scripts", "00_config.R"))

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
  library(GSVA)
  library(decoupleR)
  library(dorothea)
  library(lme4)
  library(lmerTest)
  library(emmeans)
  library(pROC)
})

set.seed(20261028)

analysis_inputs <- readRDS(out_path("results", "analysis_inputs.rds"))
discovery <- readRDS(out_path("results", "gse40553_discovery_results.rds"))
obj <- analysis_inputs$gse31348

expr <- obj$expression
meta <- obj$metadata
meta$TimePoint <- factor(
  meta$TimePoint,
  levels = c("Week 0", "Week 1", "Week 2", "Week 4", "Week 26")
)
meta$SubjectID <- factor(meta$SubjectID)

module_list <- list(
  Restoration_Core = discovery$restoration_core,
  Disease_Driver = discovery$disease_driver
)
module_list <- lapply(module_list, function(x) intersect(x, rownames(expr)))

module_scores <- t(gsva(
  ssgseaParam(expr, module_list, minSize = 5),
  verbose = FALSE
)) %>%
  as.data.frame() %>%
  rownames_to_column("SampleID") %>%
  left_join(
    meta %>% rownames_to_column("SampleID") %>%
      select(SampleID, SubjectID, TimePoint),
    by = "SampleID"
  )

data(dorothea_hs, package = "dorothea")
regulons_abc <- dorothea_hs %>%
  filter(confidence %in% c("A", "B", "C"))

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

fit_mixed_series <- function(values, feature_name) {
  d <- data.frame(
    value = as.numeric(values),
    SubjectID = meta$SubjectID,
    TimePoint = meta$TimePoint
  )
  fit <- lmer(value ~ TimePoint + (1 | SubjectID), data = d, REML = TRUE)
  emm <- emmeans(fit, ~ TimePoint)
  contrast <- contrast(emm, method = "trt.vs.ctrl", ref = 1, adjust = "none")
  as.data.frame(contrast) %>%
    mutate(
      feature = feature_name,
      overall_time_p = anova(fit)$`Pr(>F)`[1]
    )
}

module_models <- bind_rows(lapply(
  c("Restoration_Core", "Disease_Driver"),
  function(feature) fit_mixed_series(module_scores[[feature]], feature)
)) %>%
  group_by(feature) %>%
  mutate(FDR = p.adjust(p.value, method = "BH")) %>%
  ungroup()

key_tfs <- c("PAX5", "TBX21", "NFATC2")
available_key_tfs <- intersect(key_tfs, rownames(tf_activity))
tf_models <- bind_rows(lapply(
  available_key_tfs,
  function(feature) fit_mixed_series(tf_activity[feature, ], feature)
)) %>%
  group_by(feature) %>%
  mutate(FDR = p.adjust(p.value, method = "BH")) %>%
  ungroup()

score_genes <- c("CD79A", "EBF1", "FCGR1A", "GK")
missing_score_genes <- setdiff(score_genes, rownames(expr))
if (length(missing_score_genes) > 0) {
  stop("Missing score genes in GSE31348: ", paste(missing_score_genes, collapse = ", "))
}

expr_z <- t(scale(t(expr[score_genes, , drop = FALSE])))
healing_score <- (expr_z["CD79A", ] + expr_z["EBF1", ]) -
  (expr_z["FCGR1A", ] + expr_z["GK", ])

score_data <- data.frame(
  SampleID = colnames(expr),
  SubjectID = meta$SubjectID,
  TimePoint = meta$TimePoint,
  HealingScore = as.numeric(healing_score)
)

score_fit <- lmer(
  HealingScore ~ TimePoint + (1 | SubjectID),
  data = score_data,
  REML = TRUE
)
score_emm <- emmeans(score_fit, ~ TimePoint)
score_contrasts <- as.data.frame(contrast(
  score_emm, method = "trt.vs.ctrl", ref = 1, adjust = "BH"
)) %>%
  mutate(overall_time_p = anova(score_fit)$`Pr(>F)`[1])

score_roc <- roc(
  response = factor(
    ifelse(score_data$TimePoint == "Week 0", "Baseline", "Follow-up"),
    levels = c("Baseline", "Follow-up")
  ),
  predictor = score_data$HealingScore,
  subset = score_data$TimePoint %in% c("Week 0", "Week 26"),
  direction = "<",
  quiet = TRUE
)
score_auc <- data.frame(
  comparison = "Week 26 vs Week 0 (development, not validation)",
  AUC = as.numeric(auc(score_roc)),
  lower_95 = as.numeric(ci.auc(score_roc))[1],
  upper_95 = as.numeric(ci.auc(score_roc))[3]
)

write_tsv(module_scores, out_path("results", "gse31348_module_scores.tsv"))
write_tsv(module_models, out_path("results", "gse31348_module_mixed_models.tsv"))
write_tsv(tf_models, out_path("results", "gse31348_key_tf_mixed_models.tsv"))
write_tsv(score_data, out_path("results", "gse31348_healing_score.tsv"))
write_tsv(score_contrasts, out_path("results", "gse31348_score_time_contrasts.tsv"))
write_tsv(score_auc, out_path("results", "gse31348_score_development_auc.tsv"))

saveRDS(
  list(
    module_scores = module_scores,
    module_models = module_models,
    tf_activity = tf_activity,
    tf_models = tf_models,
    score_data = score_data,
    score_fit = score_fit,
    score_contrasts = score_contrasts,
    score_auc = score_auc,
    score_genes = score_genes
  ),
  out_path("results", "gse31348_development_results.rds")
)

writeLines(c(
  "GSE31348 is analyzed as a model-development cohort.",
  paste0("Participants: ", n_distinct(meta$SubjectID)),
  paste0("Samples: ", nrow(meta)),
  paste0("Cohort-standardized TBRI treatment-state AUC (Week 26 vs Week 0): ",
         round(score_auc$AUC, 4)),
  "This is treatment-state separation, not outcome prediction or independent validation."
), out_path("logs", "gse31348_development_summary.txt"))
