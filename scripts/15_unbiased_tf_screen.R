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

tf40553 <- discovery$tf_activity
tf31348 <- development$tf_activity
tf89403 <- outcome$tf_activity

meta40553 <- inputs$gse40553$pheno
meta40553$SubjectID <- factor(meta40553$SubjectID)
meta40553$Group <- factor(
  meta40553$Group,
  levels = c("LTB", "TB_0M", "TB_0.5M", "TB_2M", "TB_6M", "TB_12M")
)

meta31348 <- inputs$gse31348$metadata
meta31348$SubjectID <- factor(meta31348$SubjectID)
meta31348$TimePoint <- factor(
  meta31348$TimePoint,
  levels = c("Week 0", "Week 1", "Week 2", "Week 4", "Week 26")
)

meta89403 <- inputs$gse89403$metadata
meta89403$SubjectID <- factor(meta89403$SubjectID)
meta89403$Group <- factor(
  meta89403$Group,
  levels = c("Successful_Treatment", "Poor_Outcome")
)
meta89403$TimePoint <- factor(
  meta89403$TimePoint,
  levels = c("0_Baseline", "1_Week1", "2_Week4", "3_Week24")
)

# Donor-blocked treatment response in the validation cohort.
design31348 <- model.matrix(~ 0 + TimePoint, data = meta31348)
colnames(design31348) <- make.names(levels(meta31348$TimePoint))
dc31348 <- duplicateCorrelation(tf31348, design31348, block = meta31348$SubjectID)
fit31348 <- lmFit(
  tf31348,
  design31348,
  block = meta31348$SubjectID,
  correlation = dc31348$consensus.correlation
)
con31348 <- makeContrasts(Week26 = Week.26 - Week.0, levels = design31348)
fit31348 <- eBayes(contrasts.fit(fit31348, con31348))
res31348 <- topTable(fit31348, coef = "Week26", number = Inf, adjust.method = "BH") %>%
  rownames_to_column("TF") %>%
  transmute(
    TF,
    GSE31348_week26_logFC = logFC,
    GSE31348_week26_FDR = adj.P.Val
  )

# Donor-blocked average outcome contrast in GSE89403.
meta89403$GroupTime <- interaction(meta89403$Group, meta89403$TimePoint, drop = TRUE)
design89403 <- model.matrix(~ 0 + GroupTime, data = meta89403)
colnames(design89403) <- levels(meta89403$GroupTime)
dc89403 <- duplicateCorrelation(tf89403, design89403, block = meta89403$SubjectID)
fit89403 <- lmFit(
  tf89403,
  design89403,
  block = meta89403$SubjectID,
  correlation = dc89403$consensus.correlation
)
poor_cols <- grep("^Poor_Outcome\\.", colnames(design89403), value = TRUE)
success_cols <- grep("^Successful_Treatment\\.", colnames(design89403), value = TRUE)
outcome_contrast <- matrix(
  0,
  nrow = ncol(design89403),
  ncol = 1,
  dimnames = list(colnames(design89403), "OutcomeDifference")
)
outcome_contrast[poor_cols, 1] <- 1 / length(poor_cols)
outcome_contrast[success_cols, 1] <- -1 / length(success_cols)
fit89403 <- eBayes(contrasts.fit(fit89403, outcome_contrast))
res89403 <- topTable(fit89403, coef = "OutcomeDifference", number = Inf,
                     adjust.method = "BH") %>%
  rownames_to_column("TF") %>%
  transmute(
    TF,
    GSE89403_poor_vs_success_logFC = logFC,
    GSE89403_outcome_FDR = adj.P.Val
  )

screen <- discovery$tf_activity_table %>%
  select(TF, GSE40553_treatment_logFC = treatment_logFC,
         GSE40553_treatment_FDR = treatment_FDR) %>%
  left_join(res31348, by = "TF") %>%
  left_join(res89403, by = "TF") %>%
  left_join(
    read_tsv(out_path("results", "gse40553_paired_donor_correlations.tsv"),
             show_col_types = FALSE) %>%
      select(TF, donor_rho_restoration = rho_restoration,
             donor_FDR_restoration = FDR_restoration),
    by = "TF"
  ) %>%
  mutate(
    discovery_positive = GSE40553_treatment_logFC > 0 & GSE40553_treatment_FDR < 0.05,
    validation_positive = GSE31348_week26_logFC > 0 & GSE31348_week26_FDR < 0.05,
    outcome_lower_in_poor = GSE89403_poor_vs_success_logFC < 0 & GSE89403_outcome_FDR < 0.05,
    donor_positive = donor_rho_restoration > 0 & donor_FDR_restoration < 0.05,
    evidence_count = rowSums(cbind(
      discovery_positive, validation_positive,
      outcome_lower_in_poor, donor_positive
    ), na.rm = TRUE)
  ) %>%
  arrange(desc(evidence_count),
          pmax(GSE40553_treatment_FDR, GSE31348_week26_FDR,
               GSE89403_outcome_FDR, donor_FDR_restoration, na.rm = TRUE))

write_tsv(screen, out_path("results", "unbiased_tf_cross_cohort_screen.tsv"))

writeLines(c(
  paste0("TFs screened: ", nrow(screen)),
  paste0("TFs with all four evidence directions: ",
         sum(screen$evidence_count == 4, na.rm = TRUE)),
  paste0("TFs with at least three evidence directions: ",
         sum(screen$evidence_count >= 3, na.rm = TRUE)),
  "Significance uses FDR < 0.05 in each cohort."
), out_path("logs", "unbiased_tf_screen_summary.txt"))
