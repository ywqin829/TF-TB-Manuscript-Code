source(file.path("scripts", "00_config.R"))

suppressPackageStartupMessages({
  library(dplyr)
  library(lme4)
  library(lmerTest)
  library(emmeans)
  library(readr)
})

module_scores <- read_tsv(
  out_path("results", "gse40553_module_scores.tsv"),
  show_col_types = FALSE
) %>%
  mutate(
    Group = factor(
      Group,
      levels = c("LTB", "TB_0M", "TB_0.5M", "TB_2M", "TB_6M", "TB_12M")
    ),
    SubjectID = factor(SubjectID)
  )

tb_scores <- module_scores %>%
  filter(Group != "LTB") %>%
  droplevels()

fit_module <- function(feature) {
  fit <- lmer(
    as.formula(paste(feature, "~ Group + (1 | SubjectID)")),
    data = tb_scores,
    REML = TRUE
  )
  emm <- emmeans(fit, ~ Group)
  contrasts <- as.data.frame(emmeans::contrast(
    emm,
    method = "trt.vs.ctrl",
    ref = 1,
    adjust = "BH"
  ))
  contrasts$feature <- feature
  contrasts$overall_time_p <- anova(fit)$`Pr(>F)`[1]
  contrasts
}

module_models <- bind_rows(
  fit_module("Restoration_Core"),
  fit_module("Disease_Driver")
)

write_tsv(
  module_models,
  out_path("results", "gse40553_module_mixed_models.tsv")
)
