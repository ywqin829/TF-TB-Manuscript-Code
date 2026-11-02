source(file.path("scripts", "00_config.R"))

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
  library(limma)
  library(quantiseqr)
})

inputs <- readRDS(out_path("results", "analysis_inputs.rds"))
discovery <- readRDS(out_path("results", "gse40553_discovery_results.rds"))

load_composition <- function(filename) {
  read.delim(
    out_path("results", filename),
    check.names = FALSE
  ) %>%
    column_to_rownames("SampleID") %>%
    as.matrix() %>%
    t()
}

make_comp_pcs <- function(comp, meta, n_pc = 3) {
  missing <- setdiff(rownames(meta), colnames(comp))
  if (length(missing) > 0) {
    stop("Composition matrix is missing samples: ", paste(missing, collapse = ", "))
  }
  values <- t(comp[, rownames(meta), drop = FALSE])
  values <- pmin(pmax(values, 1e-4), 1 - 1e-4)
  values <- log(values / (1 - values))
  keep <- apply(values, 2, var, na.rm = TRUE) > 0
  values <- values[, keep, drop = FALSE]
  pca <- prcomp(scale(values), center = TRUE, scale. = FALSE)
  n_pc <- min(n_pc, ncol(pca$x))
  scores <- pca$x[, seq_len(n_pc), drop = FALSE]
  colnames(scores) <- paste0("composition_PC", seq_len(n_pc))
  scores
}

fit_gene_attenuation <- function(
    expression,
    meta,
    genes,
    module_map,
    dataset,
    factor_name,
    test_level,
    reference_level,
    composition) {
  genes <- intersect(genes, rownames(expression))
  model_meta <- meta[colnames(expression), , drop = FALSE]
  pcs <- make_comp_pcs(composition, model_meta)
  model_meta <- cbind(model_meta, as.data.frame(pcs))

  formula_base <- as.formula(paste0("~ 0 + ", factor_name))
  formula_full <- as.formula(paste0(
    "~ 0 + ", factor_name,
    " + composition_PC1 + composition_PC2 + composition_PC3"
  ))
  design_base <- model.matrix(formula_base, data = model_meta)
  design_full <- model.matrix(formula_full, data = model_meta)
  colnames(design_base) <- make.names(colnames(design_base), unique = TRUE)
  colnames(design_full) <- make.names(colnames(design_full), unique = TRUE)

  test_coefficient <- make.names(paste0(factor_name, test_level))
  reference_coefficient <- make.names(paste0(factor_name, reference_level))
  treatment_expression <- paste0(
    test_coefficient, " - ", reference_coefficient
  )
  base_contrast <- makeContrasts(
    contrasts = treatment_expression,
    levels = design_base
  )
  colnames(base_contrast) <- "Treatment"
  full_contrast <- makeContrasts(
    contrasts = treatment_expression,
    levels = design_full
  )
  colnames(full_contrast) <- "Treatment"

  expression_subset <- expression[genes, rownames(model_meta), drop = FALSE]
  block <- model_meta$SubjectID

  dc_base <- duplicateCorrelation(expression_subset, design_base, block = block)
  fit_base <- lmFit(
    expression_subset,
    design_base,
    block = block,
    correlation = dc_base$consensus.correlation
  )
  fit_base <- eBayes(contrasts.fit(fit_base, base_contrast))
  base_table <- topTable(
    fit_base,
    coef = "Treatment",
    number = Inf,
    sort.by = "none"
  )

  dc_full <- duplicateCorrelation(expression_subset, design_full, block = block)
  fit_full <- lmFit(
    expression_subset,
    design_full,
    block = block,
    correlation = dc_full$consensus.correlation
  )
  fit_full <- eBayes(fit_full)
  pc_columns <- grep("^composition_PC", colnames(design_full), value = TRUE)
  composition_test <- topTable(
    fit_full,
    coef = pc_columns,
    number = Inf,
    sort.by = "none"
  )
  fit_full_contrast <- eBayes(contrasts.fit(fit_full, full_contrast))
  full_table <- topTable(
    fit_full_contrast,
    coef = "Treatment",
    number = Inf,
    sort.by = "none"
  )

  stopifnot(identical(rownames(base_table), genes))
  stopifnot(identical(rownames(full_table), genes))
  retained <- full_table$logFC / base_table$logFC
  data.frame(
    dataset = dataset,
    gene = genes,
    module = unname(module_map[genes]),
    effect_base = base_table$logFC,
    effect_adjusted = full_table$logFC,
    effect_retained_fraction = retained,
    composition_attenuation = 1 - retained,
    sign_consistent = sign(base_table$logFC) == sign(full_table$logFC),
    treatment_FDR_base = base_table$adj.P.Val,
    treatment_FDR_adjusted = full_table$adj.P.Val,
    composition_joint_p = composition_test$P.Value,
    composition_joint_FDR = composition_test$adj.P.Val
  )
}

comp_40553 <- load_composition("gse40553_quantiseq_explicit_nonlog.tsv")
comp_31348 <- load_composition("gse31348_quantiseq_explicit_nonlog.tsv")

module_genes <- c(
  setNames(rep("Restoration", length(discovery$restoration_core)),
           discovery$restoration_core),
  setNames(rep("Disease_driver", length(discovery$disease_driver)),
           discovery$disease_driver)
)

meta_40553 <- inputs$gse40553$pheno
meta_40553$Group <- factor(
  meta_40553$Group,
  levels = c("LTB", "TB_0M", "TB_0.5M", "TB_2M", "TB_6M", "TB_12M")
)
meta_40553$SubjectID <- factor(meta_40553$SubjectID)

meta_31348 <- inputs$gse31348$metadata
meta_31348$TimePoint <- factor(
  meta_31348$TimePoint,
  levels = c("Week 0", "Week 1", "Week 2", "Week 4", "Week 26")
)
meta_31348$SubjectID <- factor(meta_31348$SubjectID)

attenuation <- bind_rows(
  fit_gene_attenuation(
    inputs$gse40553$expr_symbol,
    meta_40553,
    names(module_genes),
    module_genes,
    "GSE40553",
    "Group",
    "TB_6M",
    "TB_0M",
    comp_40553
  ),
  fit_gene_attenuation(
    inputs$gse31348$expression,
    meta_31348,
    names(module_genes),
    module_genes,
    "GSE31348",
    "TimePoint",
    "Week 26",
    "Week 0",
    comp_31348
  )
)

summary_table <- attenuation %>%
  group_by(dataset, module) %>%
  summarise(
    n_genes = n(),
    median_attenuation = median(composition_attenuation, na.rm = TRUE),
    q25_attenuation = quantile(composition_attenuation, 0.25, na.rm = TRUE),
    q75_attenuation = quantile(composition_attenuation, 0.75, na.rm = TRUE),
    prop_attenuation_ge_25pct = mean(composition_attenuation >= 0.25, na.rm = TRUE),
    prop_attenuation_ge_50pct = mean(composition_attenuation >= 0.50, na.rm = TRUE),
    prop_sign_reversed = mean(!sign_consistent, na.rm = TRUE),
    prop_composition_FDR_05 = mean(composition_joint_FDR < 0.05, na.rm = TRUE),
    .groups = "drop"
  )

cross_cohort <- attenuation %>%
  select(dataset, gene, module, composition_attenuation) %>%
  pivot_wider(names_from = dataset, values_from = composition_attenuation) %>%
  filter(is.finite(GSE40553), is.finite(GSE31348))

cross_cohort_tests <- cross_cohort %>%
  group_by(module) %>%
  summarise(
    n_genes = n(),
    spearman_attenuation = cor(
      GSE40553, GSE31348,
      method = "spearman",
      use = "complete.obs"
    ),
    .groups = "drop"
  )

effect_floor_sensitivity <- bind_rows(lapply(
  c(0.25, 0.5, 1.0),
  function(effect_floor) {
    filtered <- attenuation %>%
      filter(abs(effect_base) >= effect_floor)
    filtered %>%
      group_by(dataset, module) %>%
      summarise(
        effect_floor = effect_floor,
        n_genes = n(),
        median_attenuation = median(composition_attenuation, na.rm = TRUE),
        prop_attenuation_ge_25pct = mean(
          composition_attenuation >= 0.25,
          na.rm = TRUE
        ),
        prop_attenuation_ge_50pct = mean(
          composition_attenuation >= 0.50,
          na.rm = TRUE
        ),
        prop_sign_reversed = mean(!sign_consistent, na.rm = TRUE),
        prop_composition_FDR_05 = mean(
          composition_joint_FDR < 0.05,
          na.rm = TRUE
        ),
        .groups = "drop"
      )
  }
))

cross_cohort_by_floor <- bind_rows(lapply(
  c(0.25, 0.5, 1.0),
  function(effect_floor) {
    attenuation %>%
      filter(abs(effect_base) >= effect_floor) %>%
      select(dataset, gene, module, composition_attenuation) %>%
      pivot_wider(names_from = dataset, values_from = composition_attenuation) %>%
      filter(is.finite(GSE40553), is.finite(GSE31348)) %>%
      group_by(module) %>%
      summarise(
        effect_floor = effect_floor,
        n_genes = n(),
        spearman_attenuation = cor(
          GSE40553, GSE31348,
          method = "spearman",
          use = "complete.obs"
        ),
        .groups = "drop"
      )
  }
))

# Independent B-cell marker check in GSE89403. These markers exclude
# CD79A and EBF1 so that the check is not driven by TBRI's B-lineage genes.
score <- read_tsv(
  out_path("results", "gse89403_healing_score.tsv"),
  show_col_types = FALSE
)
independent_markers <- c("MS4A1", "CD19", "CD79B")
expression_89403 <- inputs$gse89403$expression
missing_markers <- setdiff(independent_markers, rownames(expression_89403))
if (length(missing_markers) > 0) {
  stop("Missing independent B-cell markers: ", paste(missing_markers, collapse = ", "))
}
marker_z <- t(scale(t(expression_89403[independent_markers, , drop = FALSE])))
marker_score <- colMeans(marker_z, na.rm = TRUE)

bc <- score %>%
  select(SampleID, SubjectID, HealingScore) %>%
  mutate(
    Independent_Bcell_marker = marker_score[SampleID],
    SubjectID = factor(SubjectID)
  ) %>%
  filter(is.finite(HealingScore), is.finite(Independent_Bcell_marker))

pooled <- cor.test(
  bc$HealingScore,
  bc$Independent_Bcell_marker,
  method = "spearman",
  exact = FALSE
)
means <- bc %>%
  group_by(SubjectID) %>%
  summarise(
    HealingScore = mean(HealingScore),
    Independent_Bcell_marker = mean(Independent_Bcell_marker),
    .groups = "drop"
  )
between <- cor.test(
  means$HealingScore,
  means$Independent_Bcell_marker,
  method = "spearman",
  exact = FALSE
)
within_fit <- lm(
  HealingScore ~ Independent_Bcell_marker + SubjectID,
  data = bc
)
within_coef <- summary(within_fit)$coefficients["Independent_Bcell_marker", ]
within_df <- within_fit$df.residual
within_t <- unname(within_coef["t value"])
within_rho <- sign(within_t) * sqrt(within_t^2 / (within_t^2 + within_df))

marker_correlation <- bind_rows(
  tibble(
    Component = "Pooled samples",
    Method = "Spearman correlation",
    N_samples = nrow(bc),
    N_participants = n_distinct(bc$SubjectID),
    Estimate = unname(pooled$estimate),
    P_value = pooled$p.value
  ),
  tibble(
    Component = "Between participants",
    Method = "Spearman correlation of participant means",
    N_samples = nrow(means),
    N_participants = nrow(means),
    Estimate = unname(between$estimate),
    P_value = between$p.value
  ),
  tibble(
    Component = "Within participants",
    Method = "Repeated-measures correlation controlling for participant",
    N_samples = nrow(bc),
    N_participants = n_distinct(bc$SubjectID),
    Estimate = within_rho,
    P_value = unname(within_coef["Pr(>|t|)"])
  )
)

signature_path <- system.file(
  "extdata", "TIL10_signature.txt",
  package = "quantiseqr"
)
signature <- read.delim(signature_path, check.names = FALSE)
score_genes <- c("CD79A", "EBF1", "FCGR1A", "GK")
signature_overlap <- tibble(
  Gene = score_genes,
  In_quanTIseq_fixed_signature = score_genes %in% signature$ID,
  B_cell_weight = signature$B.cells[match(score_genes, signature$ID)]
)

fraction_diagnostics <- bind_rows(lapply(
  c("GSE40553", "GSE31348", "GSE89403"),
  function(dataset) {
    filename <- switch(
      dataset,
      GSE40553 = "gse40553_quantiseq_explicit_nonlog.tsv",
      GSE31348 = "gse31348_quantiseq_explicit_nonlog.tsv",
      GSE89403 = "gse89403_quantiseq_tpm.tsv"
    )
    fractions <- read.delim(
      out_path("results", filename),
      check.names = FALSE
    )
    values <- as.matrix(fractions[, setdiff(names(fractions), "SampleID")])
    data.frame(
      dataset = dataset,
      n_samples = nrow(fractions),
      min_fraction = min(values, na.rm = TRUE),
      max_fraction = max(values, na.rm = TRUE),
      min_fraction_sum = min(rowSums(values), na.rm = TRUE),
      median_fraction_sum = median(rowSums(values), na.rm = TRUE),
      max_fraction_sum = max(rowSums(values), na.rm = TRUE),
      max_abs_sum_error = max(abs(rowSums(values) - 1), na.rm = TRUE),
      n_negative = sum(values < -1e-8, na.rm = TRUE)
    )
  }
))

write_tsv(
  attenuation,
  out_path("results", "gene_composition_attenuation.tsv")
)
write_tsv(
  summary_table,
  out_path("results", "gene_composition_attenuation_summary.tsv")
)
write_tsv(
  cross_cohort_tests,
  out_path("results", "gene_composition_attenuation_cross_cohort.tsv")
)
write_tsv(
  effect_floor_sensitivity,
  out_path("results", "gene_composition_attenuation_effect_floor.tsv")
)
write_tsv(
  cross_cohort_by_floor,
  out_path("results", "gene_composition_attenuation_cross_cohort_by_floor.tsv")
)
write_tsv(
  marker_correlation,
  out_path("results", "tbri_independent_bcell_marker_correlation.tsv")
)
write_tsv(
  signature_overlap,
  out_path("results", "tbri_quantiseq_signature_overlap.tsv")
)
write_tsv(
  fraction_diagnostics,
  out_path("results", "quantiseq_fraction_diagnostics.tsv")
)

writeLines(c(
  "Gene-level limma models quantify attenuation after adding three composition PCs.",
  "Both models use patient-blocked duplicateCorrelation inference.",
  paste0("Genes modelled: ", nrow(attenuation)),
  paste0(
    "GSE40553 restoration median attenuation: ",
    signif(summary_table$median_attenuation[
      summary_table$dataset == "GSE40553" &
        summary_table$module == "Restoration"
    ], 4)
  ),
  paste0(
    "GSE31348 restoration median attenuation: ",
    signif(summary_table$median_attenuation[
      summary_table$dataset == "GSE31348" &
        summary_table$module == "Restoration"
    ], 4)
  ),
  paste0("Independent B-cell marker pooled rho: ", signif(unname(pooled$estimate), 4)),
  paste0("Independent B-cell marker between-participant rho: ", signif(unname(between$estimate), 4)),
  paste0("Independent B-cell marker within-participant rho: ", signif(within_rho, 4)),
  "Effect-floor sensitivity used absolute base treatment effects of 0.25, 0.5 and 1.0 log2 units.",
  "These attenuation estimates are composition-associated, not causal cell-intrinsic fractions."
), out_path("logs", "composition_framework_summary.txt"))
