source(file.path("scripts", "00_config.R"))

suppressPackageStartupMessages({
  library(limma)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
})

inputs <- readRDS(out_path("results", "analysis_inputs.rds"))
discovery <- readRDS(out_path("results", "gse40553_discovery_results.rds"))
development <- readRDS(out_path("results", "gse31348_development_results.rds"))
outcome <- readRDS(out_path("results", "gse89403_outcome_results.rds"))

screen <- read_tsv(
  out_path("results", "unbiased_tf_cross_cohort_screen.tsv"),
  show_col_types = FALSE
)

top_candidates <- screen %>%
  arrange(desc(evidence_count),
          pmax(GSE40553_treatment_FDR, GSE31348_week26_FDR,
               GSE89403_outcome_FDR, donor_FDR_restoration, na.rm = TRUE)) %>%
  slice_head(n = 10) %>%
  pull(TF)

load_composition <- function(stem) {
  read.delim(
    out_path("results", paste0(stem, "_quantiseq.tsv")),
    check.names = FALSE
  ) %>%
    column_to_rownames("SampleID") %>%
    as.matrix()
}

comp_40553 <- load_composition("gse40553")
comp_31348 <- load_composition("gse31348")
comp_89403 <- load_composition("gse89403")

make_comp_pcs <- function(comp, n_pc = 3) {
  view <- log(pmin(pmax(comp, 1e-4), 1 - 1e-4) /
                (1 - pmin(pmax(comp, 1e-4), 1 - 1e-4)))
  view <- apply(view, 2, function(z) {
    finite <- is.finite(z)
    if (!any(finite)) return(rep(0, length(z)))
    z[!finite] <- median(z[finite], na.rm = TRUE)
    z
  })
  view <- view[, apply(view, 2, var, na.rm = TRUE) > 0, drop = FALSE]
  pca <- prcomp(scale(view), center = TRUE, scale. = FALSE)
  n_pc <- min(n_pc, ncol(pca$x))
  scores <- pca$x[, seq_len(n_pc), drop = FALSE]
  colnames(scores) <- paste0("composition_PC", seq_len(n_pc))
  scores
}

extract_adjusted <- function(
    tf_matrix, meta, comp, dataset,
    contrast_matrix, contrast_name) {
  message("Fitting composition-adjusted contrast for ", dataset)
  meta <- meta[colnames(tf_matrix), , drop = FALSE]
  comp <- comp[rownames(meta), , drop = FALSE]
  pcs <- make_comp_pcs(comp)

  base_terms <- switch(
    dataset,
    GSE40553 = "~ 0 + Group",
    GSE31348 = "~ 0 + TimePoint",
    GSE89403 = "~ 0 + GroupTime"
  )
  design <- model.matrix(
    as.formula(paste(base_terms, "+ composition_PC1 + composition_PC2 + composition_PC3")),
    data = cbind(meta, as.data.frame(pcs))
  )
  colnames(design) <- make.names(colnames(design))

  full_contrast <- matrix(
    0,
    nrow = ncol(design),
    ncol = ncol(contrast_matrix),
    dimnames = list(colnames(design), colnames(contrast_matrix))
  )
  shared_terms <- intersect(rownames(contrast_matrix), colnames(design))
  full_contrast[shared_terms, ] <- contrast_matrix[shared_terms, , drop = FALSE]

  dc <- duplicateCorrelation(tf_matrix, design, block = meta$SubjectID)
  fit <- lmFit(
    tf_matrix,
    design,
    block = meta$SubjectID,
    correlation = dc$consensus.correlation
  )
  fit <- eBayes(contrasts.fit(fit, full_contrast))

  topTable(
    fit,
    coef = contrast_name,
    number = Inf,
    adjust.method = "BH"
  ) %>%
    rownames_to_column("TF") %>%
    transmute(
      TF,
      !!paste0(dataset, "_adjusted_logFC") := logFC,
      !!paste0(dataset, "_adjusted_p") := P.Value,
      !!paste0(dataset, "_adjusted_FDR") := adj.P.Val
    )
}

meta_40553 <- inputs$gse40553$pheno
meta_40553$SubjectID <- factor(meta_40553$SubjectID)
meta_40553$Group <- factor(
  meta_40553$Group,
  levels = c("LTB", "TB_0M", "TB_0.5M", "TB_2M", "TB_6M", "TB_12M")
)

meta_31348 <- inputs$gse31348$metadata
meta_31348$SubjectID <- factor(meta_31348$SubjectID)
meta_31348$TimePoint <- factor(
  meta_31348$TimePoint,
  levels = c("Week 0", "Week 1", "Week 2", "Week 4", "Week 26")
)

meta_89403 <- inputs$gse89403$metadata
meta_89403$SubjectID <- factor(meta_89403$SubjectID)
meta_89403$Group <- factor(
  meta_89403$Group,
  levels = c("Successful_Treatment", "Poor_Outcome")
)
meta_89403$TimePoint <- factor(
  meta_89403$TimePoint,
  levels = c("0_Baseline", "1_Week1", "2_Week4", "3_Week24")
)
meta_89403$GroupTime <- interaction(meta_89403$Group, meta_89403$TimePoint, drop = TRUE)

design_40553 <- model.matrix(
  ~ 0 + Group,
  data = meta_40553
)
colnames(design_40553) <- make.names(colnames(design_40553))
contrast_40553 <- makeContrasts(
  Treatment6m = GroupTB_6M - GroupTB_0M,
  levels = design_40553
)

design_31348 <- model.matrix(
  ~ 0 + TimePoint,
  data = meta_31348
)
colnames(design_31348) <- make.names(colnames(design_31348))
contrast_31348 <- makeContrasts(
  Week26 = TimePointWeek.26 - TimePointWeek.0,
  levels = design_31348
)

design_89403 <- model.matrix(
  ~ 0 + GroupTime,
  data = meta_89403
)
colnames(design_89403) <- make.names(colnames(design_89403))
poor_cols <- grep("^GroupTimePoor_Outcome\\.", colnames(design_89403), value = TRUE)
success_cols <- grep("^GroupTimeSuccessful_Treatment\\.", colnames(design_89403), value = TRUE)
contrast_89403 <- matrix(
  0,
  nrow = ncol(design_89403),
  ncol = 1,
  dimnames = list(colnames(design_89403), "OutcomeDifference")
)
contrast_89403[poor_cols, 1] <- 1 / length(poor_cols)
contrast_89403[success_cols, 1] <- -1 / length(success_cols)

adjusted <- list(
  extract_adjusted(
    discovery$tf_activity,
    meta_40553,
    comp_40553,
    "GSE40553",
    contrast_40553,
    "Treatment6m"
  ),
  extract_adjusted(
    development$tf_activity,
    meta_31348,
    comp_31348,
    "GSE31348",
    contrast_31348,
    "Week26"
  ),
  extract_adjusted(
    outcome$tf_activity,
    meta_89403,
    comp_89403,
    "GSE89403",
    contrast_89403,
    "OutcomeDifference"
  )
) %>%
  purrr::reduce(full_join, by = "TF")

robustness <- screen %>%
  select(
    TF, evidence_count,
    GSE40553_treatment_logFC, GSE40553_treatment_FDR,
    GSE31348_week26_logFC, GSE31348_week26_FDR,
    GSE89403_poor_vs_success_logFC, GSE89403_outcome_FDR,
    donor_rho_restoration, donor_FDR_restoration
  ) %>%
  left_join(adjusted, by = "TF") %>%
  mutate(
    adjusted_discovery_same_direction =
      sign(GSE40553_adjusted_logFC) == sign(GSE40553_treatment_logFC),
    adjusted_validation_same_direction =
      sign(GSE31348_adjusted_logFC) == sign(GSE31348_week26_logFC),
    adjusted_outcome_same_direction =
      sign(GSE89403_adjusted_logFC) == sign(GSE89403_poor_vs_success_logFC),
    adjusted_discovery_supported =
      adjusted_discovery_same_direction & GSE40553_adjusted_FDR < 0.05,
    adjusted_validation_supported =
      adjusted_validation_same_direction & GSE31348_adjusted_FDR < 0.05,
    adjusted_outcome_supported =
      adjusted_outcome_same_direction & GSE89403_adjusted_FDR < 0.05,
    adjusted_evidence_count = rowSums(
      cbind(
        adjusted_discovery_supported,
        adjusted_validation_supported,
        adjusted_outcome_supported
      ),
      na.rm = TRUE
    )
  ) %>%
  arrange(desc(adjusted_evidence_count), evidence_count, TF)

write_tsv(
  robustness,
  out_path("results", "candidate_composition_robustness.tsv")
)

writeLines(c(
  paste0("Screen-derived top candidates: ", paste(top_candidates, collapse = ", ")),
  "Models include the original design, three quanTIseq composition PCs, and donor-blocked inference.",
  "Support requires the same effect direction and adjusted FDR < 0.05 in each cohort.",
  "Composition PCs are summaries of estimated cell fractions; residual confounding remains possible."
), out_path("logs", "candidate_composition_robustness_summary.txt"))
