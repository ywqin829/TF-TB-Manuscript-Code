#!/usr/bin/env Rscript
# LM22/CIBERSORT sensitivity analysis for GSE31348.
# Reuses the preprocessed microarray expression matrix and follows the legacy author's IOBR call
# (arrays=TRUE, LM22, 200 permutations).  It is a same-bulk-data algorithmic sensitivity analysis,
# not an independent experimental validation and not CIBERSORTx.

source("scripts/00_config.R")
suppressPackageStartupMessages({
  library(IOBR)
  library(dplyr)
  library(readr)
  library(tibble)
  library(lme4)
  library(lmerTest)
  library(ggplot2)
})

gse31348 <- readRDS(paths$gse31348)
expr <- as.matrix(gse31348$expression)
meta <- as.data.frame(gse31348$metadata) %>%
  rename(SampleID = geo_accession) %>%
  mutate(
    SubjectID = factor(SubjectID),
    TimePoint = factor(TimePoint, levels = c("Week 0", "Week 1", "Week 2", "Week 4", "Week 26"))
  )

if (ncol(expr) != 135L || nrow(meta) != 135L || n_distinct(meta$SubjectID) != 27L) {
  stop("Expected frozen GSE31348 matrix (135 records; 27 patients); found ", ncol(expr), " records and ", n_distinct(meta$SubjectID), " patients.")
}
if (!identical(colnames(expr), meta$SampleID)) stop("Expression and metadata sample order do not agree.")

# The legacy GSE31348 script used IOBR::deconvo_cibersort with arrays=TRUE and perm=200.
# IOBR's CIBERSORT backend supplies LM22 internally; output columns retain the CIBERSORT suffix.
set.seed(20261028)
ciber <- IOBR::deconvo_cibersort(expr, arrays = TRUE, perm = 200, seed = 20261028)
ciber <- as.data.frame(ciber)
if (!"ID" %in% names(ciber)) stop("Unexpected IOBR CIBERSORT result: missing sample ID column.")

normalised_names <- tolower(gsub("[^a-z0-9]+", "", names(ciber), ignore.case = TRUE))
naive_col <- names(ciber)[normalised_names == "bcellsnaivecibersort"]
memory_col <- names(ciber)[normalised_names == "bcellsmemorycibersort"]
if (length(naive_col) != 1L || length(memory_col) != 1L) {
  stop("Could not identify LM22 naive and memory B-cell columns. Available columns: ",
       paste(names(ciber), collapse = ", "))
}

ciber_out <- ciber %>%
  rename(SampleID = ID) %>%
  left_join(meta, by = "SampleID") %>%
  mutate(
    CIBERSORT_naive_B = .data[[naive_col]],
    CIBERSORT_memory_B = .data[[memory_col]],
    CIBERSORT_total_B = CIBERSORT_naive_B + CIBERSORT_memory_B
  )
if (anyNA(ciber_out$SubjectID)) stop("CIBERSORT results contain samples absent from metadata.")

quantiseq <- read_tsv(out_path("results", "gse31348_quantiseq_explicit_nonlog.tsv"), show_col_types = FALSE) %>%
  select(SampleID, quanTIseq_B_cell = `B cell`)
compare <- ciber_out %>% inner_join(quantiseq, by = "SampleID")
rho <- cor(compare$CIBERSORT_total_B, compare$quanTIseq_B_cell, method = "spearman")

# Participant is the biological unit: a random intercept protects against treating five repeated
# microarray samples per participant as independent observations.
b_model <- lmer(CIBERSORT_total_B ~ TimePoint + (1 | SubjectID), data = ciber_out, REML = FALSE)
b_coef <- as.data.frame(coef(summary(b_model))) %>%
  rownames_to_column("Term") %>%
  filter(Term == "(Intercept)" | grepl("^TimePoint", Term)) %>%
  transmute(
    Term,
    Estimate,
    StdError = `Std. Error`,
    DF = df,
    t_value = `t value`,
    P_value = `Pr(>|t|)`,
    Wald95_low = Estimate - 1.96 * StdError,
    Wald95_high = Estimate + 1.96 * StdError,
    N_records = nrow(ciber_out),
    N_participants = n_distinct(ciber_out$SubjectID)
  )

summary_table <- ciber_out %>%
  group_by(TimePoint) %>%
  summarise(
    N_records = n(),
    N_participants = n_distinct(SubjectID),
    median_CIBERSORT_total_B = median(CIBERSORT_total_B),
    IQR_low = quantile(CIBERSORT_total_B, 0.25),
    IQR_high = quantile(CIBERSORT_total_B, 0.75),
    .groups = "drop"
  ) %>%
  mutate(
    method = "CIBERSORT (LM22 via IOBR)",
    comparison_with_quanTIseq = paste0("Spearman rho across 135 matched samples = ", sprintf("%.3f", rho))
  )

write_tsv(ciber_out, out_path("results", "gse31348_cibersort_lm22_legacy_settings.tsv"))
write_tsv(summary_table, out_path("results", "gse31348_cibersort_lm22_bcell_trajectory.tsv"))
write_tsv(b_coef, out_path("results", "gse31348_cibersort_lm22_bcell_mixed_model.tsv"))
write_tsv(
  data.frame(
    Measure = "CIBERSORT LM22 total B cells (naive + memory) versus quanTIseq B cell",
    Spearman_rho = rho,
    N_records = nrow(compare),
    N_participants = n_distinct(compare$SubjectID),
    Interpretation = "Same bulk expression data; algorithmic sensitivity only, not independent validation."
  ),
  out_path("results", "gse31348_cibersort_quantiseq_bcell_agreement.tsv")
)

p <- ggplot(ciber_out, aes(TimePoint, CIBERSORT_total_B)) +
  geom_boxplot(aes(fill = TimePoint), width = 0.65, outlier.shape = NA, alpha = 0.75) +
  geom_line(aes(group = SubjectID), colour = "grey45", alpha = 0.28, linewidth = 0.35) +
  geom_point(aes(colour = TimePoint), position = position_jitter(width = 0.08, height = 0), alpha = 0.65, size = 1.4) +
  scale_fill_brewer(palette = "Blues", guide = "none") +
  scale_colour_brewer(palette = "Blues", guide = "none") +
  labs(
    title = "GSE31348 LM22/CIBERSORT sensitivity analysis",
    subtitle = paste0("Total B cells = naïve + memory\nCIBERSORT–quanTIseq Spearman rho = ", sprintf("%.2f", rho)),
    x = "Treatment time", y = "Estimated total B-cell fraction"
  ) +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 11.5))
ggsave(out_path("figures", "FigureS4_GSE31348_CIBERSORT_LM22_Bcell_Trajectory.pdf"), p, width = 7.0, height = 4.8)

writeLines(c(
  "GSE31348 CIBERSORT sensitivity analysis.",
  "Method: IOBR::deconvo_cibersort, arrays=TRUE, internal LM22, perm=200, seed=20261028.",
  "Input: existing frozen preprocessed GSE31348 microarray expression matrix (135 repeated records from 27 participants).",
  "B-cell measure: LM22 naive B + memory B fractions.",
  paste0("CIBERSORT--quanTIseq total B-cell Spearman rho: ", sprintf("%.6f", rho)),
  "Interpretation boundary: both algorithms use the same bulk transcriptome; agreement is algorithmic robustness, not independent cell-count validation."
), out_path("logs", "gse31348_cibersort_lm22_sensitivity_summary.txt"))

message("Completed GSE31348 CIBERSORT/LM22 sensitivity analysis for ", nrow(ciber_out), " repeated records from ", n_distinct(ciber_out$SubjectID), " participants.")
