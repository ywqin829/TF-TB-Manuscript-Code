source(file.path("scripts", "00_config.R"))

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
  library(lme4)
  library(lmerTest)
})

set.seed(20261028)

inputs <- readRDS(out_path("results", "analysis_inputs.rds"))
discovery <- readRDS(out_path("results", "gse40553_discovery_results.rds"))
development <- readRDS(out_path("results", "gse31348_development_results.rds"))
outcome <- readRDS(out_path("results", "gse89403_outcome_results.rds"))
metanet_upstream <- read_tsv(
  out_path("results", "metanet_upstream_candidates.tsv"),
  show_col_types = FALSE
)
composition <- read_tsv(
  out_path("results", "tf_activity_composition_models.tsv"),
  show_col_types = FALSE
)

candidates <- metanet_upstream %>%
  filter(network == "whole_blood_tissue_specific", rank <= 15) %>%
  pull(candidate_upstream_TF) %>%
  unique() %>%
  intersect(
    Reduce(
      intersect,
      list(
        rownames(discovery$tf_activity),
        rownames(development$tf_activity),
        rownames(outcome$tf_activity)
      )
    )
  )

extract_effect <- function(matrix, meta, formula, term_pattern, label) {
  rows <- lapply(intersect(candidates, rownames(matrix)), function(tf) {
    d <- cbind(data.frame(value = as.numeric(matrix[tf, rownames(meta)])), meta)
    fit <- lmer(as.formula(formula), data = d, REML = TRUE)
    tab <- as.data.frame(anova(fit))
    tab$term <- rownames(tab)
    p <- tab$`Pr(>F)`[grepl(term_pattern, tab$term)][1]
    data.frame(
      TF = tf,
      cohort = label,
      p_value = p
    )
  })
  bind_rows(rows)
}

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

meta_89403 <- inputs$gse89403$metadata
meta_89403$Group <- factor(
  meta_89403$Group,
  levels = c("Successful_Treatment", "Poor_Outcome")
)
meta_89403$TimePoint <- factor(
  meta_89403$TimePoint,
  levels = c("0_Baseline", "1_Week1", "2_Week4", "3_Week24")
)
meta_89403$SubjectID <- factor(meta_89403$SubjectID)

effects <- bind_rows(
  extract_effect(
    discovery$tf_activity, meta_40553,
    "value ~ Group + (1 | SubjectID)", "^Group$", "GSE40553"
  ),
  extract_effect(
    development$tf_activity, meta_31348,
    "value ~ TimePoint + (1 | SubjectID)", "^TimePoint$", "GSE31348"
  ),
  extract_effect(
    outcome$tf_activity, meta_89403,
    "value ~ Group + (1 | SubjectID)", "^Group$", "GSE89403"
  )
) %>%
  group_by(cohort) %>%
  mutate(FDR = p.adjust(p_value, method = "BH")) %>%
  ungroup() %>%
  pivot_wider(names_from = cohort, values_from = c(p_value, FDR))

direction <- data.frame(
  TF = intersect(candidates, rownames(discovery$tf_activity)),
  GSE40553_treatment_logFC = as.numeric(
    discovery$tf_activity_table[match(
      intersect(candidates, rownames(discovery$tf_activity)),
      discovery$tf_activity_table$TF
    ), "treatment_logFC"]
  )
)

discovery_sig <- discovery$tf_activity_table %>%
  select(TF, treatment_logFC, treatment_FDR) %>%
  rename(
    GSE40553_treatment_logFC = treatment_logFC,
    GSE40553_treatment_FDR = treatment_FDR
  )

candidate_table <- metanet_upstream %>%
  filter(network == "whole_blood_tissue_specific", rank <= 15) %>%
  group_by(candidate_upstream_TF) %>%
  summarise(
    METANet_targets = paste(sort(unique(target_TF)), collapse = ";"),
    best_METANet_rank = min(rank),
    best_METANet_score = max(score),
    .groups = "drop"
  ) %>%
  rename(TF = candidate_upstream_TF) %>%
  filter(TF %in% candidates) %>%
  left_join(discovery_sig, by = "TF") %>%
  left_join(effects, by = "TF") %>%
  left_join(
    composition %>%
      select(dataset, TF, delta_R2) %>%
      pivot_wider(names_from = dataset, values_from = delta_R2,
                  names_prefix = "composition_deltaR2_"),
    by = "TF"
  ) %>%
  mutate(
    direction_consistent = case_when(
      GSE40553_treatment_logFC > 0 &
        FDR_GSE31348 < 0.05 &
        FDR_GSE89403 < 0.05 ~ "yes",
      TRUE ~ "not established"
    )
  ) %>%
  arrange(FDR_GSE89403, FDR_GSE31348, FDR_GSE40553)

write_tsv(candidate_table, out_path("results", "metanet_candidate_cross_cohort.tsv"))

writeLines(c(
  paste0("Candidate TFs tested: ", nrow(candidate_table)),
  paste0("Direction-consistent candidates: ",
         sum(candidate_table$direction_consistent == "yes")),
  "Cross-cohort presence is candidate prioritisation, not causal validation."
), out_path("logs", "metanet_candidate_validation_summary.txt"))
