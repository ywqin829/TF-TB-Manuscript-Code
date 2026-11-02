source(file.path("scripts", "00_config.R"))

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
  library(decoupleR)
  library(dorothea)
  library(lme4)
  library(lmerTest)
  library(emmeans)
  library(pROC)
  library(dcurves)
  library(rms)
})

set.seed(20261028)

analysis_inputs <- readRDS(out_path("results", "analysis_inputs.rds"))
obj <- analysis_inputs$gse89403
expr <- obj$expression
meta <- obj$metadata

meta$TimePoint <- factor(
  meta$TimePoint,
  levels = c("0_Baseline", "1_Week1", "2_Week4", "3_Week24")
)
meta$Group <- factor(
  meta$Group,
  levels = c("Successful_Treatment", "Poor_Outcome")
)
meta$SubjectID <- factor(meta$SubjectID)

data(dorothea_hs, package = "dorothea")
regulons_abc <- dorothea_hs %>%
  filter(confidence %in% c("A", "B", "C"))

tf_cache <- out_path("results", "gse89403_tf_activity_cache.rds")
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

expr_z <- t(scale(t(expr)))

module_genes <- list(
  PAX5_target_module = c("CD19", "CD79A", "FCER2", "POU2AF1", "BLK"),
  TBX21_target_module = c("IFNG", "CD40LG"),
  B_cell_marker_module = c("MS4A1", "CD19", "CD79A", "CD79B"),
  TNK_marker_module = c("CD3D", "CD3E", "NKG7", "GNLY")
)
module_genes <- lapply(module_genes, intersect, rownames(expr_z))

module_scores <- as.data.frame(lapply(module_genes, function(g) {
  colMeans(expr_z[g, , drop = FALSE])
}))
module_scores$SampleID <- rownames(meta)

key_tfs <- c("PAX5", "TBX21", "NFATC2")
available_key_tfs <- intersect(key_tfs, rownames(tf_activity))

feature_matrix <- rbind(
  tf_activity[available_key_tfs, , drop = FALSE],
  t(as.matrix(module_scores[, names(module_genes), drop = FALSE]))
)

fit_group_time <- function(values, feature_name, adjust_covariate = NULL) {
  d <- data.frame(
    value = as.numeric(values),
    Group = meta$Group,
    TimePoint = meta$TimePoint,
    SubjectID = meta$SubjectID
  )
  if (!is.null(adjust_covariate)) d$covariate <- as.numeric(adjust_covariate)

  formula_txt <- if (is.null(adjust_covariate)) {
    "value ~ Group * TimePoint + (1 | SubjectID)"
  } else {
    "value ~ Group * TimePoint + covariate + (1 | SubjectID)"
  }

  fit <- lmer(
    as.formula(formula_txt),
    data = d,
    REML = TRUE,
    control = lmerControl(optimizer = "bobyqa")
  )
  aov <- as.data.frame(anova(fit))
  aov$term <- rownames(aov)
  aov$feature <- feature_name
  aov$adjustment <- if (is.null(adjust_covariate)) "unadjusted" else "cell-marker adjusted"
  rownames(aov) <- NULL

  emm <- emmeans(fit, ~ Group | TimePoint)
  between <- as.data.frame(emmeans::contrast(
    emm, method = "revpairwise", adjust = "none"
  ))
  between$feature <- feature_name
  between$comparison <- "Poor_Outcome vs Successful_Treatment"
  between$adjustment <- if (is.null(adjust_covariate)) "unadjusted" else "cell-marker adjusted"

  list(fit = fit, anova = aov, between = between)
}

model_specs <- list(
  PAX5 = list(
    values = tf_activity["PAX5", ],
    covariate = module_scores$B_cell_marker_module
  ),
  TBX21 = list(
    values = tf_activity["TBX21", ],
    covariate = module_scores$TNK_marker_module
  ),
  NFATC2 = list(
    values = tf_activity["NFATC2", ],
    covariate = module_scores$TNK_marker_module
  ),
  PAX5_target_module = list(
    values = module_scores$PAX5_target_module,
    covariate = module_scores$B_cell_marker_module
  ),
  TBX21_target_module = list(
    values = module_scores$TBX21_target_module,
    covariate = module_scores$TNK_marker_module
  ),
  B_cell_marker_module = list(
    values = module_scores$B_cell_marker_module,
    covariate = NULL
  ),
  TNK_marker_module = list(
    values = module_scores$TNK_marker_module,
    covariate = NULL
  )
)

models <- lapply(names(model_specs), function(feature) {
  spec <- model_specs[[feature]]
  out <- list(unadjusted = fit_group_time(spec$values, feature, NULL))
  if (!is.null(spec$covariate)) {
    out$adjusted <- fit_group_time(spec$values, feature, spec$covariate)
  }
  out
})
names(models) <- names(model_specs)

model_anova <- bind_rows(lapply(models, function(x) {
  bind_rows(x$unadjusted$anova, x$adjusted$anova)
})) %>%
  group_by(term, adjustment) %>%
  mutate(FDR = p.adjust(`Pr(>F)`, method = "BH")) %>%
  ungroup()

model_between <- bind_rows(lapply(models, function(x) {
  bind_rows(x$unadjusted$between, x$adjusted$between)
})) %>%
  group_by(feature, adjustment) %>%
  mutate(FDR = p.adjust(p.value, method = "BH")) %>%
  ungroup()

score_genes <- c("CD79A", "EBF1", "FCGR1A", "GK")
missing_score_genes <- setdiff(score_genes, rownames(expr))
if (length(missing_score_genes) > 0) {
  stop("Missing score genes in GSE89403: ", paste(missing_score_genes, collapse = ", "))
}

score_expr_z <- t(scale(t(expr[score_genes, , drop = FALSE])))
healing_score <- (score_expr_z["CD79A", ] + score_expr_z["EBF1", ]) -
  (score_expr_z["FCGR1A", ] + score_expr_z["GK", ])

score_data <- data.frame(
  SampleID = rownames(meta),
  SubjectID = meta$SubjectID,
  Group = meta$Group,
  TimePoint = meta$TimePoint,
  HealingScore = as.numeric(healing_score)
)

score_mixed <- lmer(
  HealingScore ~ Group * TimePoint + (1 | SubjectID),
  data = score_data,
  REML = TRUE,
  control = lmerControl(optimizer = "bobyqa")
)
score_anova <- as.data.frame(anova(score_mixed))
score_anova$term <- rownames(score_anova)
rownames(score_anova) <- NULL

score_emm <- emmeans(score_mixed, ~ Group | TimePoint)
score_between_time <- as.data.frame(emmeans::contrast(
  score_emm, method = "revpairwise", adjust = "BH"
)) %>%
  mutate(comparison = "Poor_Outcome vs Successful_Treatment")

eval_score_time <- function(x) {
  y <- factor(
    as.character(x$Group),
    levels = c("Poor_Outcome", "Successful_Treatment")
  )
  r <- roc(
    response = y,
    predictor = x$HealingScore,
    levels = c("Poor_Outcome", "Successful_Treatment"),
    direction = "<",
    quiet = TRUE
  )
  ci <- as.numeric(ci.auc(r, method = "delong"))
  data.frame(
    TimePoint = unique(x$TimePoint)[1],
    n_success = sum(x$Group == "Successful_Treatment"),
    n_poor = sum(x$Group == "Poor_Outcome"),
    AUC = as.numeric(auc(r)),
    lower_95 = ci[1],
    upper_95 = ci[3]
  )
}

score_auc <- score_data %>%
  group_split(TimePoint) %>%
  lapply(eval_score_time) %>%
  bind_rows()

loo_logistic_prob <- function(x) {
  y <- as.integer(x$Group == "Poor_Outcome")
  pred <- numeric(length(y))
  for (i in seq_along(y)) {
    train <- x[-i, , drop = FALSE]
    fit <- suppressWarnings(glm(
      y[-i] ~ HealingScore,
      data = train,
      family = binomial()
    ))
    pred[i] <- predict(
      fit,
      newdata = x[i, , drop = FALSE],
      type = "response"
    )
  }
  x$loo_probability_poor <- pred
  x
}

score_predictions <- score_data %>%
  group_split(TimePoint) %>%
  lapply(loo_logistic_prob) %>%
  bind_rows()

calibration_metrics <- score_predictions %>%
  group_by(TimePoint) %>%
  group_modify(~ {
    y <- as.integer(.x$Group == "Poor_Outcome")
    p <- pmin(pmax(.x$loo_probability_poor, 1e-6), 1 - 1e-6)
    lp <- qlogis(p)
    cal <- suppressWarnings(glm(y ~ lp, family = binomial()))
    data.frame(
      calibration_intercept = unname(coef(cal)[1]),
      calibration_slope = unname(coef(cal)[2]),
      brier_score = mean((p - y)^2)
    )
  }) %>%
  ungroup()

dca_results <- score_predictions %>%
  group_split(TimePoint) %>%
  lapply(function(x) {
    x$outcome_poor <- as.integer(x$Group == "Poor_Outcome")
    out <- dca(
      outcome_poor ~ loo_probability_poor,
      data = x,
      thresholds = seq(0.02, 0.80, by = 0.02)
    )$dca
    out$TimePoint <- unique(x$TimePoint)[1]
    out
  }) %>%
  bind_rows()

bootstrap_auc <- function(x, B = 2000) {
  set.seed(20261028)
  y <- factor(
    as.character(x$Group),
    levels = c("Poor_Outcome", "Successful_Treatment")
  )
  values <- replicate(B, {
    idx <- sample(seq_len(nrow(x)), replace = TRUE)
    if (length(unique(y[idx])) < 2) return(NA_real_)
    as.numeric(auc(roc(
      y[idx], x$HealingScore[idx],
      levels = levels(y), direction = "<", quiet = TRUE
    )))
  })
  data.frame(
    TimePoint = unique(x$TimePoint)[1],
    bootstrap_lower_95 = quantile(values, 0.025, na.rm = TRUE),
    bootstrap_upper_95 = quantile(values, 0.975, na.rm = TRUE)
  )
}

bootstrap_ci <- score_data %>%
  group_split(TimePoint) %>%
  lapply(bootstrap_auc) %>%
  bind_rows()

score_auc <- score_auc %>%
  left_join(bootstrap_ci, by = "TimePoint") %>%
  left_join(calibration_metrics, by = "TimePoint")

write_tsv(module_scores, out_path("results", "gse89403_module_scores.tsv"))
write_tsv(model_anova, out_path("results", "gse89403_mixed_model_anova.tsv"))
write_tsv(model_between, out_path("results", "gse89403_group_contrasts_by_time.tsv"))
write_tsv(score_data, out_path("results", "gse89403_healing_score.tsv"))
write_tsv(score_anova, out_path("results", "gse89403_score_mixed_model_anova.tsv"))
write_tsv(score_between_time, out_path("results", "gse89403_score_group_contrasts.tsv"))
write_tsv(score_auc, out_path("results", "gse89403_score_external_auc.tsv"))
write_tsv(score_predictions, out_path("results", "gse89403_loo_score_predictions.tsv"))
write_tsv(dca_results, out_path("results", "gse89403_score_dca.tsv"))

saveRDS(
  list(
    tf_activity = tf_activity,
    module_scores = module_scores,
    models = models,
    model_anova = model_anova,
    model_between = model_between,
    score_data = score_data,
    score_mixed = score_mixed,
    score_anova = score_anova,
    score_between_time = score_between_time,
    score_predictions = score_predictions,
    score_auc = score_auc,
    dca = dca_results
  ),
  out_path("results", "gse89403_outcome_results.rds")
)

writeLines(c(
  "All longitudinal models use patient-level random intercepts.",
  "The TBX21 variable is read directly from tf_activity['TBX21', ].",
  "The cohort-standardized TBRI is associated with outcome in GSE89403; this is not external validation of deployment-ready fixed parameters.",
  "Calibration and DCA use leave-one-out probabilities within each time point.",
  "The poor-outcome group contains only 8 patients; interval estimates reflect this limitation."
), out_path("logs", "gse89403_outcome_summary.txt"))
