source(file.path("scripts", "00_config.R"))

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
  library(lme4)
  library(lmerTest)
  library(pROC)
})

set.seed(20261028)

inputs <- readRDS(out_path("results", "analysis_inputs.rds"))
discovery <- readRDS(out_path("results", "gse40553_discovery_results.rds"))
development <- readRDS(out_path("results", "gse31348_development_results.rds"))
outcome <- readRDS(out_path("results", "gse89403_outcome_results.rds"))

load_corrected_composition <- function(filename) {
  path <- out_path("results", filename)
  if (!file.exists(path)) {
    stop(
      "Corrected quanTIseq output not found: ", path,
      ". Run scripts/24_quantiseq_scale_audit.R first."
    )
  }
  read.delim(path, check.names = FALSE) %>%
    column_to_rownames("SampleID") %>%
    as.matrix() %>%
    t()
}

comp_40553 <- load_corrected_composition(
  "gse40553_quantiseq_explicit_nonlog.tsv"
)
comp_31348 <- load_corrected_composition(
  "gse31348_quantiseq_explicit_nonlog.tsv"
)
comp_89403 <- load_corrected_composition("gse89403_quantiseq_tpm.tsv")

align_composition <- function(comp, meta) {
  missing <- setdiff(rownames(meta), colnames(comp))
  if (length(missing) > 0) {
    stop("Composition matrix is missing samples: ", paste(missing, collapse = ", "))
  }
  comp[, rownames(meta), drop = FALSE]
}

make_comp_pcs <- function(comp, n_pc = 3) {
  view <- t(comp)
  view <- log(pmin(pmax(view, 1e-4), 1 - 1e-4) /
                (1 - pmin(pmax(view, 1e-4), 1 - 1e-4)))
  view <- apply(view, 2, function(z) {
    finite <- is.finite(z)
    if (!any(finite)) return(rep(0, length(z)))
    z[!finite] <- stats::median(z[finite], na.rm = TRUE)
    z
  })
  keep <- apply(view, 2, stats::var, na.rm = TRUE) > 0
  view <- view[, keep, drop = FALSE]
  if (ncol(view) < 1) stop("No non-constant cell-type columns for composition PCA.")
  view <- scale(view)
  pca <- prcomp(view, center = TRUE, scale. = FALSE)
  n_pc <- min(n_pc, ncol(pca$x))
  scores <- pca$x[, seq_len(n_pc), drop = FALSE]
  colnames(scores) <- paste0("composition_PC", seq_len(n_pc))
  scores
}

fit_variance_models <- function(values, meta, formula_base, formula_full) {
  d <- cbind(
    data.frame(value = as.numeric(values)),
    meta
  )
  # Both formulas are data-column names, never free-form user input.
  base_fit <- lmer(as.formula(formula_base), data = d, REML = FALSE)
  full_fit <- lmer(as.formula(formula_full), data = d, REML = FALSE)
  lr <- anova(base_fit, full_fit)
  r2 <- function(fit) {
    fixed <- var(predict(fit, re.form = NA))
    vc <- as.data.frame(VarCorr(fit))
    random <- sum(vc$vcov)
    residual <- sigma(fit)^2
    fixed / (fixed + random + residual)
  }
  list(
    base = base_fit,
    full = full_fit,
    lr_p = lr$`Pr(>Chisq)`[2],
    r2_base = r2(base_fit),
    r2_full = r2(full_fit),
    delta_r2 = r2(full_fit) - r2(base_fit)
  )
}

run_tf_composition <- function(tf_matrix, meta, comp, dataset, formula_base, formula_full) {
  meta <- meta[colnames(tf_matrix), , drop = FALSE]
  stopifnot(identical(rownames(meta), colnames(tf_matrix)))
  pcs <- make_comp_pcs(align_composition(comp, meta))
  model_meta <- cbind(
    meta,
    as.data.frame(pcs)
  )
  candidates <- intersect(
    c(
      "MAF", "SNAPC4", "TBX21", "ARNTL", "STAT4", "MEF2B", "RUNX3",
      "DUX4", "ZBED1", "NFATC2", "PAX5", "CEBPG", "AHR", "BACH2",
      "SPIB", "RBPJ"
    ),
    rownames(tf_matrix)
  )
  bind_rows(lapply(candidates, function(tf) {
    result <- fit_variance_models(
      tf_matrix[tf, ],
      model_meta,
      formula_base,
      formula_full
    )
    data.frame(
      dataset = dataset,
      TF = tf,
      n_samples = nrow(model_meta),
      n_subjects = n_distinct(model_meta$SubjectID),
      LRT_p = result$lr_p,
      R2_base = result$r2_base,
      R2_with_composition = result$r2_full,
      delta_R2 = result$delta_r2
    )
  }))
}

meta_40553 <- inputs$gse40553$pheno
meta_40553$Group <- factor(
  meta_40553$Group,
  levels = c("LTB", "TB_0M", "TB_0.5M", "TB_2M", "TB_6M", "TB_12M")
)
meta_40553$SubjectID <- factor(meta_40553$SubjectID)
comp40553_pcs <- make_comp_pcs(comp_40553)
meta_40553_pcs <- cbind(meta_40553, as.data.frame(comp40553_pcs))

meta_31348 <- inputs$gse31348$metadata
meta_31348$TimePoint <- factor(
  meta_31348$TimePoint,
  levels = c("Week 0", "Week 1", "Week 2", "Week 4", "Week 26")
)
meta_31348$SubjectID <- factor(meta_31348$SubjectID)
comp31348_pcs <- make_comp_pcs(comp_31348)
meta_31348_pcs <- cbind(meta_31348, as.data.frame(comp31348_pcs))

meta_89403 <- inputs$gse89403$metadata
meta_89403$Group <- factor(
  meta_89403$Group,
  levels = c("Successful_Treatment", "Poor_Outcome")
)
meta_89403$TimePoint <- factor(
  meta_89403$TimePoint,
  levels = c("0_Baseline", "1_Week1", "2_Week4", "3_Week24")
)
meta_89403$SubjectID <- factor(meta_89403$SubjectID)
comp89403_pcs <- make_comp_pcs(comp_89403)
meta_89403_pcs <- cbind(meta_89403, as.data.frame(comp89403_pcs))

tf_models <- bind_rows(
  run_tf_composition(
    discovery$tf_activity, meta_40553, comp_40553, "GSE40553",
    "value ~ Group + (1 | SubjectID)",
    "value ~ Group + composition_PC1 + composition_PC2 + composition_PC3 + (1 | SubjectID)"
  ),
  run_tf_composition(
    development$tf_activity, meta_31348, comp_31348, "GSE31348",
    "value ~ TimePoint + (1 | SubjectID)",
    "value ~ TimePoint + composition_PC1 + composition_PC2 + composition_PC3 + (1 | SubjectID)"
  ),
  run_tf_composition(
    outcome$tf_activity, meta_89403, comp_89403, "GSE89403",
    "value ~ Group * TimePoint + (1 | SubjectID)",
    "value ~ Group * TimePoint + composition_PC1 + composition_PC2 + composition_PC3 + (1 | SubjectID)"
  )
) %>%
  group_by(dataset) %>%
  mutate(FDR = p.adjust(LRT_p, method = "BH")) %>%
  ungroup() %>%
  arrange(dataset, LRT_p)

# Incremental outcome value of TBRI over the estimated B-cell fraction.
score_external <- read.delim(
  out_path("results", "gse89403_healing_score.tsv"),
  check.names = FALSE
)
score_external$SampleID <- as.character(score_external$SampleID)
score_external$TimePoint <- factor(
  score_external$TimePoint,
  levels = c("0_Baseline", "1_Week1", "2_Week4", "3_Week24")
)

comp_89403_df <- as.data.frame(t(comp_89403)) %>%
  rownames_to_column("SampleID")
score_comp <- score_external %>%
  left_join(comp_89403_df, by = "SampleID") %>%
  mutate(
    B_fraction = `B cell`,
    outcome_poor = as.integer(Group == "Poor_Outcome")
  )

incremental_value <- bind_rows(lapply(
  levels(score_comp$TimePoint),
  function(tp) {
    x <- filter(score_comp, TimePoint == tp)
    if (length(unique(x$outcome_poor)) < 2) return(NULL)
    base <- glm(outcome_poor ~ B_fraction, data = x, family = binomial())
    full <- glm(outcome_poor ~ B_fraction + HealingScore, data = x, family = binomial())
    # The outcome is modelled as probability of poor outcome. Flip the
    # probability so that higher values mean favorable outcome, matching
    # the direction used for clinical recovery.
    base_prob <- 1 - predict(base, type = "response")
    full_prob <- 1 - predict(full, type = "response")
    y <- factor(x$Group, levels = c("Poor_Outcome", "Successful_Treatment"))
    base_roc <- roc(y, base_prob, levels = levels(y), direction = "<", quiet = TRUE)
    full_roc <- roc(y, full_prob, levels = levels(y), direction = "<", quiet = TRUE)
    test <- pROC::roc.test(base_roc, full_roc, method = "delong")
    data.frame(
      TimePoint = tp,
      n = nrow(x),
      n_poor = sum(x$outcome_poor),
      AUC_Bcell_only = as.numeric(auc(base_roc)),
      AUC_Bcell_plus_TBRI = as.numeric(auc(full_roc)),
      delta_AUC = as.numeric(auc(full_roc) - auc(base_roc)),
      DeLong_p = test$p.value,
      LRT_p = anova(base, full, test = "Chisq")$`Pr(>Chi)`[2]
    )
  }
))

write_tsv(tf_models, out_path("results", "tf_activity_composition_models.tsv"))
write_tsv(incremental_value, out_path("results", "tbri_incremental_value.tsv"))
write_tsv(
  as.data.frame(t(comp_40553)) %>% rownames_to_column("SampleID"),
  out_path("results", "gse40553_quantiseq.tsv")
)
write_tsv(
  as.data.frame(t(comp_31348)) %>% rownames_to_column("SampleID"),
  out_path("results", "gse31348_quantiseq.tsv")
)
write_tsv(
  as.data.frame(t(comp_89403)) %>% rownames_to_column("SampleID"),
  out_path("results", "gse89403_quantiseq.tsv")
)

writeLines(c(
  "Composition analysis used explicitly non-log microarray input and raw-count-derived TPM for GSE89403.",
  "TF activity was modelled with composition principal components as covariates.",
  "TBRI incremental value was tested against a model containing only the estimated B-cell fraction.",
  paste0("TF rows:", nrow(tf_models)),
  paste0("Incremental-value rows:", nrow(incremental_value))
), out_path("logs", "composition_analysis_summary.txt"))
