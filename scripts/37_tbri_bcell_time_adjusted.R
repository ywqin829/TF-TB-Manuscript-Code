#!/usr/bin/env Rscript
# Time-adjusted within/between participant analysis of the GSE89403 TBRI--B-cell association.
# This script intentionally reuses the final retained 312 participant-time records and existing
# quanTIseq output; it does not reprocess RNA-seq data or refit the TBRI gene-selection model.

source("scripts/00_config.R")
suppressPackageStartupMessages({
  library(dplyr)
  library(lme4)
  library(lmerTest)
  library(ggplot2)
  library(readr)
})

score <- read_tsv(out_path("results", "gse89403_healing_score.tsv"), show_col_types = FALSE)
comp <- read_tsv(out_path("results", "gse89403_quantiseq_tpm.tsv"), show_col_types = FALSE)
inputs <- readRDS(out_path("results", "analysis_inputs.rds"))

stopifnot(all(c("SampleID", "SubjectID", "TimePoint", "HealingScore") %in% names(score)))
stopifnot(all(c("SampleID", "B cell") %in% names(comp)))

# Independent marker score deliberately excludes the two B-lineage genes in TBRI (CD79A, EBF1).
# Genes are standardized across the fixed final GSE89403 analysis set, as for the existing score.
independent_markers <- c("MS4A1", "CD19", "CD79B")
expr <- inputs$gse89403$expression
missing_markers <- setdiff(independent_markers, rownames(expr))
if (length(missing_markers) > 0L) stop("Missing marker genes: ", paste(missing_markers, collapse = ", "))
expr_z <- t(scale(t(expr)))
marker_score <- colMeans(expr_z[independent_markers, , drop = FALSE])
marker_df <- data.frame(SampleID = names(marker_score), Independent_B_marker_score = unname(marker_score))

dat <- score %>%
  inner_join(comp %>% select(SampleID, B_cell_fraction = `B cell`), by = "SampleID") %>%
  inner_join(marker_df, by = "SampleID") %>%
  mutate(
    SubjectID = factor(SubjectID),
    TimePoint = factor(TimePoint, levels = c("0_Baseline", "1_Week1", "2_Week4", "3_Week24")),
    # quanTIseq proportions include zeros.  The continuity correction is fixed from this analysis set.
    B_cell_logit = qlogis(pmin(pmax(B_cell_fraction, 1 / (2 * n())), 1 - 1 / (2 * n())))
  )

if (nrow(dat) != 312L || n_distinct(dat$SubjectID) != 82L) {
  stop("Expected the frozen final GSE89403 set (312 records from 82 participants); got ",
       nrow(dat), " records from ", n_distinct(dat$SubjectID), " participants.")
}
if (anyNA(dat$TimePoint)) stop("Unexpected GSE89403 time point label.")

add_within_between <- function(data, variable) {
  data %>%
    group_by(SubjectID) %>%
    mutate(
      exposure_between = mean(.data[[variable]], na.rm = TRUE),
      exposure_within = .data[[variable]] - exposure_between
    ) %>%
    ungroup()
}

fit_exposure <- function(data, variable, label, scale_label) {
  model_data <- add_within_between(data, variable)
  null_fit <- lmer(
    HealingScore ~ exposure_between + TimePoint + (1 | SubjectID),
    data = model_data, REML = FALSE
  )
  full_fit <- lmer(
    HealingScore ~ exposure_within + exposure_between + TimePoint + (1 | SubjectID),
    data = model_data, REML = FALSE
  )
  lrt <- anova(null_fit, full_fit)
  lrt_p <- lrt$`Pr(>Chisq)`[2]
  co <- as.data.frame(coef(summary(full_fit)))
  get_row <- function(term, component) {
    r <- co[term, , drop = FALSE]
    data.frame(
      Exposure = label,
      Scale = scale_label,
      Component = component,
      Estimate = r$Estimate,
      StdError = r$`Std. Error`,
      DF = r$df,
      t_value = r$`t value`,
      P_value = r$`Pr(>|t|)`,
      Wald95_low = r$Estimate - 1.96 * r$`Std. Error`,
      Wald95_high = r$Estimate + 1.96 * r$`Std. Error`,
      Within_component_LRT_P = lrt_p,
      N_records = nrow(model_data),
      N_participants = n_distinct(model_data$SubjectID),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }
  list(
    coefficients = bind_rows(
      get_row("exposure_within", "Within-participant deviation"),
      get_row("exposure_between", "Participant mean")
    ),
    comparison = data.frame(
      Exposure = label,
      Scale = scale_label,
      Null_logLik = as.numeric(logLik(null_fit)),
      Full_logLik = as.numeric(logLik(full_fit)),
      LRT_chisq = lrt$Chisq[2],
      LRT_df = lrt$Df[2],
      LRT_P = lrt_p,
      N_records = nrow(model_data),
      N_participants = n_distinct(model_data$SubjectID),
      stringsAsFactors = FALSE
    ),
    model_data = model_data,
    model = full_fit
  )
}

# Primary composition analysis uses a bounded logit transformation. Raw fraction is a prespecified
# scale sensitivity check. The independent marker score is a biological sensitivity check.
fits <- list(
  fit_exposure(dat, "B_cell_logit", "quanTIseq B-cell fraction", "logit fraction (primary)"),
  fit_exposure(dat, "B_cell_fraction", "quanTIseq B-cell fraction", "raw fraction (sensitivity)"),
  fit_exposure(dat, "Independent_B_marker_score", "Independent B-cell marker score", "mean expression z-score")
)

coef_table <- bind_rows(lapply(fits, `[[`, "coefficients"))
comparison_table <- bind_rows(lapply(fits, `[[`, "comparison"))
write_tsv(coef_table, out_path("results", "gse89403_tbri_bcell_time_adjusted.tsv"))
write_tsv(comparison_table, out_path("results", "gse89403_tbri_bcell_time_adjusted_model_comparisons.tsv"))

# Partial correlation is descriptive; mixed-model inference remains the primary test.
# Residualize both axes on the same treatment-time and participant terms.
plot_dat <- fits[[1]]$model_data %>%
  mutate(
    TBRI_time_adjusted = residuals(lm(HealingScore ~ TimePoint + SubjectID, data = .)),
    B_cell_time_adjusted = residuals(lm(B_cell_logit ~ TimePoint + SubjectID, data = .))
  )
partial_r <- cor(plot_dat$B_cell_time_adjusted, plot_dat$TBRI_time_adjusted)
primary_p <- coef_table$P_value[
  coef_table$Exposure == "quanTIseq B-cell fraction" &
    coef_table$Scale == "logit fraction (primary)" &
    coef_table$Component == "Within-participant deviation"
]
stopifnot(length(primary_p) == 1L, is.finite(primary_p), is.finite(partial_r))
write_tsv(
  tibble(
    Descriptive_partial_Pearson_r = partial_r,
    Mixed_model_within_P = primary_p,
    N_records = nrow(plot_dat),
    N_participants = n_distinct(plot_dat$SubjectID)
  ),
  out_path("results", "gse89403_tbri_bcell_time_adjusted_figure_stats.tsv")
)
p <- ggplot(plot_dat, aes(B_cell_time_adjusted, TBRI_time_adjusted)) +
  geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey75") +
  geom_vline(xintercept = 0, linewidth = 0.25, colour = "grey75") +
  geom_point(alpha = 0.62, size = 1.65, colour = "#176B87") +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.75, colour = "#B7472A") +
  labs(
    x = "Logit quanTIseq B-cell fraction residual\nafter treatment time and participant adjustment",
    y = "TBRI residual after treatment time\nand participant adjustment",
    title = "Time-adjusted TBRI–B-cell association",
    subtitle = sprintf("Partial Pearson r = %.3f (descriptive); mixed-model P = %.2e", partial_r, primary_p)
  ) +
  theme_classic(base_size = 13) +
  theme(plot.subtitle = element_text(size = 11.5, margin = margin(b = 8)))
ggsave(out_path("figures", "FigureS3_TBRI_Bcell_TimeAdjusted.pdf"), p, width = 7.0, height = 4.8)

write_tsv(
  dat %>% select(SampleID, SubjectID, TimePoint, HealingScore, B_cell_fraction, B_cell_logit, Independent_B_marker_score),
  out_path("results", "gse89403_tbri_bcell_time_adjusted_input.tsv")
)

message("Completed time-adjusted TBRI--B-cell analysis using ", nrow(dat), " records from ", n_distinct(dat$SubjectID), " participants.")
