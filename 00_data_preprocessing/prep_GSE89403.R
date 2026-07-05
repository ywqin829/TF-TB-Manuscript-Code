# ============================================================
# prep_GSE89403.R — Outcome Clinical Cohort
# ============================================================
# Input:  log2Exp.csv, TBCatalysis.txt (from GEO GSE89403)
# Output: GSE89403_processed.rds (mat_final, meta_final, tf_acts)
# ============================================================

library(dplyr)
library(tibble)
library(data.table)
library(GEOquery)
library(decoupleR)
library(dorothea)

# ── 1. Load expression ──────────────────────────────────────
log2_exp <- fread("log2Exp.csv", data.table = FALSE)
mat_final <- log2_exp %>%
  tibble::column_to_rownames("V1") %>%
  as.matrix()

# ── 2. Load metadata ────────────────────────────────────────
gset <- getGEO(filename = "GSE89403_series_matrix.txt.gz")
pdata <- pData(gset)

meta_final <- pdata %>%
  dplyr::select(geo_accession, `subject:ch1`, `time:ch1`, `treatmentresult:ch1`) %>%
  rename(SubjectID = `subject:ch1`, TimeRaw = `time:ch1`, OutcomeRaw = `treatmentresult:ch1`) %>%
  mutate(
    TimePoint = case_when(
      TimeRaw == "DX" ~ "0_Baseline",
      TimeRaw == "day_7" ~ "1_Week1",
      TimeRaw == "week_4" ~ "2_Week4",
      TimeRaw == "week_24" ~ "3_Week24",
      TRUE ~ "Excluded"
    ),
    Group = case_when(
      OutcomeRaw == "Definite Cure" ~ "Successful_Treatment",
      OutcomeRaw == "Not Cured" ~ "Poor_Outcome",
      TRUE ~ "Excluded"
    )
  ) %>%
  filter(Group != "Excluded", TimePoint != "Excluded") %>%
  filter(geo_accession %in% colnames(mat_final))

# Align
mat_final <- mat_final[, meta_final$geo_accession]
cat(sprintf("Final matrix: %d genes x %d samples\n", nrow(mat_final), ncol(mat_final)))

# ── 3. TF activity (DoRothEA + decoupleR) ──────────────────
data(dorothea_hs, package = "dorothea")
regulons <- dorothea_hs %>% filter(confidence %in% c("A", "B", "C"))

tf_acts <- decoupleR::run_wmean(
  mat = mat_final,
  network = regulons,
  .source = "tf",
  .target = "target",
  .mor = "mor",
  minsize = 5,
  verbose = FALSE
)

tf_matrix <- tf_acts %>%
  filter(statistic == "norm_wmean") %>%
  pivot_wider(id_cols = source, names_from = condition, values_from = score) %>%
  column_to_rownames("source") %>%
  as.matrix()

# ── 4. Save ─────────────────────────────────────────────────
saveRDS(list(
  mat_final = mat_final,
  meta_final = meta_final,
  tf_matrix = tf_matrix,
  tf_acts = tf_acts
), file = "GSE89403_processed.rds")

cat("Done: GSE89403_processed.rds saved.\n")
