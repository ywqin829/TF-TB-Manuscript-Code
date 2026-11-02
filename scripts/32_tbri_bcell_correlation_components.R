source(file.path("scripts", "00_config.R"))

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
})

score <- read_tsv(
  out_path("results", "gse89403_healing_score.tsv"),
  show_col_types = FALSE
)
composition <- read_tsv(
  out_path("results", "gse89403_quantiseq_tpm.tsv"),
  show_col_types = FALSE
)

dat <- score %>%
  select(SampleID, SubjectID, TimePoint, HealingScore) %>%
  inner_join(
    composition %>% select(SampleID, B_cell = `B cell`),
    by = "SampleID"
  ) %>%
  filter(
    is.finite(HealingScore),
    is.finite(B_cell),
    !is.na(SubjectID)
  ) %>%
  mutate(SubjectID = factor(SubjectID))

if (nrow(dat) != nrow(score)) {
  stop("TBRI and B-cell data did not align for all GSE89403 score samples.")
}
if (n_distinct(dat$SubjectID) != 82) {
  stop("Expected 82 GSE89403 participants in TBRI-B-cell analysis.")
}

pooled_test <- cor.test(
  dat$HealingScore,
  dat$B_cell,
  method = "spearman",
  exact = FALSE
)

subject_means <- dat %>%
  group_by(SubjectID) %>%
  summarise(
    HealingScore = mean(HealingScore),
    B_cell = mean(B_cell),
    .groups = "drop"
  )
between_test <- cor.test(
  subject_means$HealingScore,
  subject_means$B_cell,
  method = "spearman",
  exact = FALSE
)

# Repeated-measures correlation is the partial correlation of TBRI and the
# B-cell fraction after conditioning on participant identity. This is
# algebraically the same fixed-effects ANCOVA definition used by rmcorr.
within_fit <- lm(
  HealingScore ~ B_cell + SubjectID,
  data = dat
)
within_coef <- summary(within_fit)$coefficients["B_cell", ]
within_df <- within_fit$df.residual
within_t <- unname(within_coef["t value"])
within_r <- sign(within_t) * sqrt(within_t^2 / (within_t^2 + within_df))
within_p <- unname(within_coef["Pr(>|t|)"])

results <- bind_rows(
  tibble(
    Component = "Pooled samples",
    Method = "Spearman correlation",
    N_samples = nrow(dat),
    N_participants = n_distinct(dat$SubjectID),
    Estimate = unname(pooled_test$estimate),
    P_value = pooled_test$p.value
  ),
  tibble(
    Component = "Between participants",
    Method = "Spearman correlation of participant means",
    N_samples = nrow(subject_means),
    N_participants = nrow(subject_means),
    Estimate = unname(between_test$estimate),
    P_value = between_test$p.value
  ),
  tibble(
    Component = "Within participants",
    Method = "Repeated-measures correlation controlling for participant",
    N_samples = nrow(dat),
    N_participants = n_distinct(dat$SubjectID),
    Estimate = within_r,
    P_value = within_p
  )
)

write_tsv(
  results,
  out_path("results", "gse89403_tbri_bcell_correlation_components.tsv")
)
writeLines(c(
  paste0("GSE89403 TBRI-B-cell records aligned: ", nrow(dat)),
  paste0("Participants: ", n_distinct(dat$SubjectID)),
  paste0("Pooled Spearman rho: ", signif(unname(pooled_test$estimate), 4)),
  paste0("Between-participant Spearman rho: ", signif(unname(between_test$estimate), 4)),
  paste0("Within-participant repeated-measures rho: ", signif(within_r, 4)),
  paste0("Within-participant P: ", signif(within_p, 4)),
  "The pooled coefficient treats repeated samples as correlated observations and is descriptive, not an independent-patient correlation."
), out_path("logs", "gse89403_tbri_bcell_correlation_components_summary.txt"))
