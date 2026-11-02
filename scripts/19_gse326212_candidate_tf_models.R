source(file.path("scripts", "00_config.R"))

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
})

pseudo <- read_tsv(
  out_path("results", "gse326212_candidate_tf_pseudobulk.tsv"),
  show_col_types = FALSE
)
meta <- read_csv(
  "C:/Users/yw_qi/Downloads/tb_host_immune/GSE326212/outputs/reanalysis/metadata_authoritative.csv",
  show_col_types = FALSE
)
sample_covariates <- read_csv(
  "C:/Users/yw_qi/Downloads/tb_host_immune/GSE326212/outputs/02_obs_meta.csv",
  show_col_types = FALSE
) %>%
  select(sample, sex, smoking, project)

model_data <- pseudo %>%
  filter(n_cells >= 20) %>%
  left_join(
    meta %>%
      select(sample, subject_id, group_norm, in_main_subject_set),
    by = c("sample", "subject_id", "group_norm")
  ) %>%
  left_join(sample_covariates, by = "sample") %>%
  filter(
    !is.na(in_main_subject_set),
    in_main_subject_set,
    group_norm %in% c("PTB", "IGRAneg", "IGRApos")
  ) %>%
  mutate(
    group_norm = factor(group_norm, levels = c("IGRAneg", "IGRApos", "PTB"))
  )

subject_unit_audit <- model_data %>%
  distinct(sample, subject_id, group_norm) %>%
  group_by(group_norm) %>%
  summarise(
    samples = n(),
    unique_subjects = n_distinct(subject_id),
    max_samples_per_subject = max(table(subject_id)),
    .groups = "drop"
  )
write_tsv(
  subject_unit_audit,
  out_path("results", "gse326212_subject_unit_audit.tsv")
)

fit_one <- function(d, comparison) {
  if (n_distinct(d$group_norm) < 2 || nrow(d) < 8) return(NULL)
  d$group_norm <- droplevels(d$group_norm)
  covariates <- c("sex", "smoking")[
    vapply(d[c("sex", "smoking")], nlevels, integer(1)) > 1
  ]
  fit <- lm(
    reformulate(c("group_norm", covariates), response = "mean_expression"),
    data = d
  )
  coef_table <- summary(fit)$coefficients
  term <- grep("^group_norm", rownames(coef_table), value = TRUE)
  if (length(term) == 0) return(NULL)
  term <- term[1]
  data.frame(
    comparison = comparison,
    TF = unique(d$TF),
    celltype_coarse = unique(d$celltype_coarse),
    n_IGRAneg = sum(d$group_norm == "IGRAneg"),
    n_IGRApos = sum(d$group_norm == "IGRApos"),
    n_PTB = sum(d$group_norm == "PTB"),
    coefficient = term,
    estimate = coef_table[term, "Estimate"],
    std_error = coef_table[term, "Std. Error"],
    p_value = coef_table[term, "Pr(>|t|)"]
  )
}

comparisons <- list(
  PTB_vs_IGRAneg = model_data %>% filter(group_norm %in% c("PTB", "IGRAneg")),
  PTB_vs_IGRApos = model_data %>% filter(group_norm %in% c("PTB", "IGRApos"))
)

results <- bind_rows(lapply(names(comparisons), function(comp) {
  bind_rows(lapply(split(
    comparisons[[comp]],
    interaction(comparisons[[comp]]$TF,
                comparisons[[comp]]$celltype_coarse,
                drop = TRUE)
  ), fit_one, comparison = comp))
})) %>%
  group_by(comparison) %>%
  mutate(FDR = p.adjust(p_value, method = "BH")) %>%
  ungroup() %>%
  arrange(comparison, p_value, desc(abs(estimate)))

write_tsv(
  results,
  out_path("results", "gse326212_candidate_tf_models.tsv")
)

writeLines(c(
  paste0("Model rows: ", nrow(results)),
  "Unit of inference: patient sample x coarse cell type; cells are not replicates.",
  "Models compare PTB with IGRA-negative or IGRA-positive BAL samples using donor-level pseudobulk.",
  "Covariates: sex and smoking. Project/batch could not be fully separated from group and remains a limitation.",
  "This is cross-sectional in vivo disease-state evidence, not a treatment-response analysis."
), out_path("logs", "gse326212_candidate_tf_models_summary.txt"))
