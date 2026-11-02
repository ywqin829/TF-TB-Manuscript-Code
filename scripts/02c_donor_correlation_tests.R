source(file.path("scripts", "00_config.R"))

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

donor_delta <- read_tsv(
  out_path("results", "gse40553_paired_donor_tf_deltas.tsv"),
  show_col_types = FALSE
)

test_cor <- function(x, module_delta) {
  complete <- complete.cases(x$delta_tf_activity, x[[module_delta]])
  if (sum(complete) < 4) return(c(rho = NA, p = NA))
  fit <- cor.test(
    x$delta_tf_activity[complete],
    x[[module_delta]][complete],
    method = "spearman",
    exact = FALSE
  )
  c(rho = unname(fit$estimate), p = fit$p.value)
}

results <- donor_delta %>%
  group_by(TF) %>%
  group_modify(~ {
    restoration <- test_cor(.x, "delta_restoration")
    disease <- test_cor(.x, "delta_disease")
    data.frame(
      n_paired_donors = sum(complete.cases(
        .x$delta_tf_activity, .x$delta_restoration
      )),
      rho_restoration = restoration["rho"],
      p_restoration = restoration["p"],
      rho_disease = disease["rho"],
      p_disease = disease["p"]
    )
  }) %>%
  ungroup() %>%
  mutate(
    FDR_restoration = p.adjust(p_restoration, method = "BH"),
    FDR_disease = p.adjust(p_disease, method = "BH")
  ) %>%
  arrange(desc(rho_restoration))

write_tsv(
  results,
  out_path("results", "gse40553_paired_donor_correlations.tsv")
)
