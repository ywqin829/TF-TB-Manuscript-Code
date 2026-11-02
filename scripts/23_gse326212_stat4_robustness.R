source(file.path("scripts", "00_config.R"))

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
  library(lme4)
  library(lmerTest)
})

pseudo <- read_tsv(
  out_path("results", "gse326212_candidate_tf_pseudobulk.tsv"),
  show_col_types = FALSE
)
meta <- read_csv(
  "C:/Users/yw_qi/Downloads/tb_host_immune/GSE326212/outputs/reanalysis/metadata_authoritative.csv",
  show_col_types = FALSE
)
covariates <- read_csv(
  "C:/Users/yw_qi/Downloads/tb_host_immune/GSE326212/outputs/02_obs_meta.csv",
  show_col_types = FALSE
) %>%
  select(sample, sex, smoking, project)

base <- pseudo %>%
  filter(TF == "STAT4", n_cells >= 20) %>%
  left_join(
    meta %>% select(sample, subject_id, group_norm, in_main_subject_set),
    by = c("sample", "subject_id", "group_norm")
  ) %>%
  left_join(covariates, by = "sample") %>%
  filter(
    !is.na(in_main_subject_set),
    in_main_subject_set,
    group_norm %in% c("IGRAneg", "IGRApos", "PTB")
  ) %>%
  mutate(
    group_norm = factor(group_norm, levels = c("IGRAneg", "IGRApos", "PTB")),
    mean_positive = if_else(
      percent_expressing > 0,
      mean_expression / (percent_expressing / 100),
      0
    )
  )

extract_lm <- function(data, outcome, label, covariates = "sex + smoking") {
  formula <- as.formula(
    paste0(outcome, " ~ group_norm + ", covariates)
  )
  fit <- lm(formula, data = data)
  co <- summary(fit)$coefficients
  term <- grep("^group_normPTB$", rownames(co), value = TRUE)
  if (!length(term)) return(NULL)
  ci <- confint(fit, parm = term)
  data.frame(
    analysis = label,
    outcome = outcome,
    n_IGRAneg = sum(data$group_norm == "IGRAneg"),
    n_PTB = sum(data$group_norm == "PTB"),
    estimate = co[term, "Estimate"],
    std_error = co[term, "Std. Error"],
    lower_95 = ci[1],
    upper_95 = ci[2],
    p_value = co[term, "Pr(>|t|)"]
  )
}

moam <- base %>%
  filter(celltype_coarse == "Mo-derived_AM")

models <- bind_rows(
  extract_lm(moam, "mean_expression", "MoAM group only", covariates = "1"),
  extract_lm(moam, "mean_expression", "MoAM sex + smoking"),
  extract_lm(
    moam, "mean_expression", "MoAM sex + smoking + n_cells",
    covariates = "sex + smoking + scale(n_cells)"
  ),
  extract_lm(moam, "percent_expressing", "MoAM detection rate"),
  extract_lm(moam, "mean_positive", "MoAM mean among positive cells")
)

meta_neut <- read_csv(
  "C:/Users/yw_qi/Downloads/tb_host_immune/GSE326212/outputs/reanalysis/metadata_authoritative.csv",
  show_col_types = FALSE
) %>%
  select(sample, neut_split)
moam <- moam %>%
  select(-any_of("neut_split")) %>%
  left_join(meta_neut, by = "sample") %>%
  mutate(
    group_neut = case_when(
      group_norm == "IGRAneg" ~ "IGRAneg",
      group_norm == "PTB" & neut_split == "neut_lo" ~ "PTB_neut_lo",
      group_norm == "PTB" & neut_split == "neut_hi" ~ "PTB_neut_hi",
      TRUE ~ NA_character_
    ),
    group_neut = factor(
      group_neut,
      levels = c("IGRAneg", "PTB_neut_lo", "PTB_neut_hi")
    )
  )

neut_fit <- lm(
  mean_expression ~ group_neut + sex + smoking,
  data = filter(moam, !is.na(group_neut))
)
neut_co <- summary(neut_fit)$coefficients
for (term in c("group_neutPTB_neut_lo", "group_neutPTB_neut_hi")) {
  ci <- confint(neut_fit, parm = term)
  models <- bind_rows(
    models,
    data.frame(
      analysis = paste("MoAM", sub("group_neut", "", term), "vs IGRAneg"),
      outcome = "mean_expression",
      n_IGRAneg = sum(moam$group_neut == "IGRAneg", na.rm = TRUE),
      n_PTB = sum(moam$group_neut == sub("group_neut", "", term), na.rm = TRUE),
      estimate = neut_co[term, "Estimate"],
      std_error = neut_co[term, "Std. Error"],
      lower_95 = ci[1],
      upper_95 = ci[2],
      p_value = neut_co[term, "Pr(>|t|)"]
    )
  )
}

project_fit <- lm(
  mean_expression ~ group_norm + sex + smoking + factor(project),
  data = moam
)
project_co <- summary(project_fit)$coefficients
project_ci <- confint(project_fit, parm = "group_normPTB")
models <- bind_rows(
  models,
  data.frame(
    analysis = "MoAM project fixed effects",
    outcome = "mean_expression",
    n_IGRAneg = sum(moam$group_norm == "IGRAneg"),
    n_PTB = sum(moam$group_norm == "PTB"),
    estimate = project_co["group_normPTB", "Estimate"],
    std_error = project_co["group_normPTB", "Std. Error"],
    lower_95 = project_ci[1],
    upper_95 = project_ci[2],
    p_value = project_co["group_normPTB", "Pr(>|t|)"]
  )
)

mixed_fit <- lmer(
  mean_expression ~ group_norm + sex + smoking + (1 | project),
  data = moam,
  REML = TRUE
)
mixed_co <- summary(mixed_fit)$coefficients
models <- bind_rows(
  models,
  data.frame(
    analysis = "MoAM project random intercept",
    outcome = "mean_expression",
    n_IGRAneg = sum(moam$group_norm == "IGRAneg"),
    n_PTB = sum(moam$group_norm == "PTB"),
    estimate = mixed_co["group_normPTB", "Estimate"],
    std_error = mixed_co["group_normPTB", "Std. Error"],
    lower_95 = mixed_co["group_normPTB", "Estimate"] -
      1.96 * mixed_co["group_normPTB", "Std. Error"],
    upper_95 = mixed_co["group_normPTB", "Estimate"] +
      1.96 * mixed_co["group_normPTB", "Std. Error"],
    p_value = mixed_co["group_normPTB", "Pr(>|t|)"]
  )
)

projects <- unique(moam$project)
lopo <- bind_rows(lapply(projects, function(excluded) {
  data <- filter(moam, project != excluded)
  if (n_distinct(data$group_norm) < 2 || nrow(data) < 8) return(NULL)
  fit <- lm(mean_expression ~ group_norm + sex + smoking, data = data)
  co <- summary(fit)$coefficients
  data.frame(
    excluded_project = excluded,
    n_IGRAneg = sum(data$group_norm == "IGRAneg"),
    n_PTB = sum(data$group_norm == "PTB"),
    estimate = co["group_normPTB", "Estimate"],
    p_value = co["group_normPTB", "Pr(>|t|)"]
  )
}))

myeloid_types <- c(
  "Alveolar_mac", "Mo-derived_AM", "Monocyte", "Neutrophil", "cDC1"
)
myeloid <- base %>%
  filter(celltype_coarse %in% myeloid_types) %>%
  mutate(celltype_coarse = factor(celltype_coarse, levels = myeloid_types))

myeloid_fit <- lmer(
  mean_expression ~ group_norm + celltype_coarse + sex + smoking +
    (1 | sample),
  data = myeloid
)
myeloid_co <- summary(myeloid_fit)$coefficients
myeloid_effect <- data.frame(
  analysis = "Myeloid compartments, sample random intercept",
  outcome = "mean_expression",
  n_IGRAneg_samples = n_distinct(myeloid$sample[myeloid$group_norm == "IGRAneg"]),
  n_PTB_samples = n_distinct(myeloid$sample[myeloid$group_norm == "PTB"]),
  estimate = myeloid_co["group_normPTB", "Estimate"],
  std_error = myeloid_co["group_normPTB", "Std. Error"],
  p_value = myeloid_co["group_normPTB", "Pr(>|t|)"]
)

write_tsv(
  models,
  out_path("results", "gse326212_stat4_model_robustness.tsv")
)
write_tsv(
  lopo,
  out_path("results", "gse326212_stat4_leave_one_project_out.tsv")
)
write_tsv(
  myeloid_effect,
  out_path("results", "gse326212_stat4_myeloid_composite.tsv")
)
write_tsv(
  moam,
  out_path("results", "gse326212_stat4_moam_samples.tsv")
)

writeLines(c(
  paste0("Mo-derived AM samples: ", nrow(moam)),
  paste0("Projects represented in Mo-derived AM: ", paste(projects, collapse = ", ")),
  paste0("Leave-one-project-out models: ", nrow(lopo)),
  paste0("Leave-one-project-out positive estimates: ", sum(lopo$estimate > 0, na.rm = TRUE), "/", nrow(lopo)),
  paste0("Leave-one-project-out nominal p<0.05: ", sum(lopo$p_value < 0.05, na.rm = TRUE), "/", nrow(lopo)),
  "Project is partially confounded with disease group; project-adjusted and leave-one-project-out models are sensitivity analyses, not proof of independence."
), out_path("logs", "gse326212_stat4_robustness_summary.txt"))
