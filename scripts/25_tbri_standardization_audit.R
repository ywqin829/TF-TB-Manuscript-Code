source(file.path("scripts", "00_config.R"))

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(pROC)
})

inputs <- readRDS(out_path("results", "analysis_inputs.rds"))
dev <- readRDS(out_path("results", "gse31348_development_results.rds"))
outcome <- readRDS(out_path("results", "gse89403_outcome_results.rds"))

score_genes <- c("CD79A", "EBF1", "FCGR1A", "GK")

train_expr <- inputs$gse31348$expression[score_genes, , drop = FALSE]
train_means <- rowMeans(train_expr, na.rm = TRUE)
train_sds <- apply(train_expr, 1, sd, na.rm = TRUE)
if (any(!is.finite(train_sds)) || any(train_sds == 0)) {
  stop("Non-finite or zero training SD for TBRI genes.")
}

eval_expr <- inputs$gse89403$expression[score_genes, , drop = FALSE]
fixed_z <- sweep(eval_expr, 1, train_means, "-")
fixed_z <- sweep(fixed_z, 1, train_sds, "/")
fixed_score <- (fixed_z["CD79A", ] + fixed_z["EBF1", ]) -
  (fixed_z["FCGR1A", ] + fixed_z["GK", ])

current <- outcome$score_data
fixed <- data.frame(
  SampleID = names(fixed_score),
  FixedScore = as.numeric(fixed_score)
)
comparison <- current %>%
  left_join(fixed, by = "SampleID")

auc_task <- function(data, predictor, task) {
  bind_rows(lapply(split(data, data$TimePoint), function(x) {
    y <- factor(
      as.character(x$Group),
      levels = c("Poor_Outcome", "Successful_Treatment")
    )
    r <- roc(
      y, x[[predictor]],
      levels = levels(y), direction = "<", quiet = TRUE
    )
    data.frame(
      task = task,
      predictor = predictor,
      TimePoint = as.character(x$TimePoint[1]),
      n_success = sum(x$Group == "Successful_Treatment"),
      n_poor = sum(x$Group == "Poor_Outcome"),
      AUC = as.numeric(auc(r)),
      lower_95 = as.numeric(ci.auc(r, method = "delong"))[1],
      upper_95 = as.numeric(ci.auc(r, method = "delong"))[3]
    )
  }))
}

auc_current <- auc_task(
  comparison, "HealingScore",
  "Outcome discrimination using cohort-internal z standardization"
)
auc_fixed <- auc_task(
  comparison, "FixedScore",
  "Outcome discrimination using GSE31348-frozen mean and SD"
)

baseline <- comparison %>%
  select(SubjectID, TimePoint, HealingScore) %>%
  filter(TimePoint == "0_Baseline") %>%
  select(SubjectID, BaselineScore = HealingScore)
delta_data <- comparison %>%
  left_join(baseline, by = "SubjectID") %>%
  filter(!is.na(BaselineScore)) %>%
  mutate(DeltaScore = HealingScore - BaselineScore)

auc_delta <- auc_task(
  delta_data, "DeltaScore",
  "Outcome discrimination using within-participant baseline change"
) %>%
  filter(TimePoint != "0_Baseline")

auc_table <- bind_rows(auc_current, auc_fixed, auc_delta)

write_tsv(
  comparison,
  out_path("results", "gse89403_tbri_standardization_comparison.tsv")
)
write_tsv(
  auc_table,
  out_path("results", "gse89403_tbri_standardization_auc.tsv")
)

writeLines(c(
  "Current TBRI uses gene-wise z scores calculated within each cohort.",
  "A frozen-parameter sensitivity score used GSE31348 gene means and SDs applied directly to GSE89403.",
  "Within-participant change was calculated only for participants with a baseline sample.",
  "The development GSE31348 AUC is treatment-state separation (Week 26 vs baseline), not outcome prediction.",
  paste0("GSE89403 samples with baseline for delta score: ",
         n_distinct(delta_data$SubjectID[delta_data$TimePoint == "1_Week1"]),
         " at Week 1, ",
         n_distinct(delta_data$SubjectID[delta_data$TimePoint == "2_Week4"]),
         " at Week 4, and ",
         n_distinct(delta_data$SubjectID[delta_data$TimePoint == "3_Week24"]),
         " at Week 24.")
), out_path("logs", "tbri_standardization_audit_summary.txt"))
