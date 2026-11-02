source(file.path("scripts", "00_config.R"))

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
})

raw <- read_tsv(
  out_path("results", "gse326212_stat4_raw_detection.tsv"),
  show_col_types = FALSE
)
meta <- read_csv(
  "C:/Users/yw_qi/Downloads/tb_host_immune/GSE326212/outputs/reanalysis/metadata_authoritative.csv",
  show_col_types = FALSE
) %>%
  select(sample, group_norm, in_main_subject_set)
covariates <- read_csv(
  "C:/Users/yw_qi/Downloads/tb_host_immune/GSE326212/outputs/02_obs_meta.csv",
  show_col_types = FALSE
) %>%
  select(sample, sex, smoking, project)

dat <- raw %>%
  filter(celltype_coarse == "Mo-derived_AM") %>%
  left_join(meta, by = c("sample", "group_norm")) %>%
  left_join(covariates, by = "sample") %>%
  filter(
    in_main_subject_set,
    group_norm %in% c("IGRAneg", "PTB")
  ) %>%
  mutate(
    group_norm = factor(group_norm, levels = c("IGRAneg", "PTB"))
  )

fit_one <- function(outcome) {
  fit <- lm(
    reformulate(
      c("group_norm", "sex", "smoking", "project"),
      response = outcome
    ),
    data = dat
  )
  co <- summary(fit)$coefficients
  ci <- confint(fit, parm = "group_normPTB")
  data.frame(
    outcome = outcome,
    n_IGRAneg = sum(dat$group_norm == "IGRAneg"),
    n_PTB = sum(dat$group_norm == "PTB"),
    estimate = co["group_normPTB", "Estimate"],
    std_error = co["group_normPTB", "Std. Error"],
    lower_95 = ci[1],
    upper_95 = ci[2],
    p_value = co["group_normPTB", "Pr(>|t|)"]
  )
}

dat$log1p_mean_stat4_count <- log1p(dat$mean_stat4_count)
results <- bind_rows(lapply(
  c("log1p_mean_stat4_count", "percent_count_ge1", "percent_count_ge2"),
  fit_one
))

write_tsv(
  results,
  out_path("results", "gse326212_stat4_raw_threshold_models.tsv")
)
