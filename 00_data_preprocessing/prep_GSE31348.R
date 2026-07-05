# ============================================================
# prep_GSE31348.R — Validation Cohort 1 (Affymetrix)
# ============================================================
# Input:  GSE31348_RAW/*.CEL.gz (from GEO)
# Output: GSE31348_processed.rds (expression + metadata)
# ============================================================

library(affy)
library(hgu133plus2.db)
library(AnnotationDbi)
library(dplyr)
library(GEOquery)
library(stringr)

# ── 1. Read CEL files & RMA ─────────────────────────────────
raw_data <- ReadAffy(celfile.path = "GSE31348_RAW")
eset <- rma(raw_data)
expr_matrix <- exprs(eset)

# Clean column names
colnames(expr_matrix) <- str_extract(colnames(expr_matrix), "GSM[0-9]+")

# ── 2. Load metadata ────────────────────────────────────────
gset <- getGEO(filename = "GSE31348_series_matrix.txt.gz")
pdat <- pData(gset)

pdat_clean <- pdat %>%
  dplyr::select(geo_accession, `patient_id:ch1`, `time point:ch1`) %>%
  dplyr::rename(SubjectID = `patient_id:ch1`, TimePoint = `time point:ch1`)

rownames(pdat_clean) <- pdat_clean$geo_accession
pdat_sorted <- pdat_clean[colnames(expr_matrix), ]
pdat_sorted$TimePoint <- factor(pdat_sorted$TimePoint,
                                levels = c("Week 0", "Week 1", "Week 2", "Week 4", "Week 26"))
stopifnot(identical(colnames(expr_matrix), rownames(pdat_sorted)))

# ── 3. Probe → Gene ────────────────────────────────────────
gene_symbols <- mapIds(hgu133plus2.db,
                       keys = rownames(expr_matrix),
                       column = "SYMBOL", keytype = "PROBEID",
                       multiVals = "first")

expr_df <- as.data.frame(expr_matrix)
expr_df$Symbol <- gene_symbols
expr_df <- na.omit(expr_df)

expr_final <- expr_df %>%
  dplyr::group_by(Symbol) %>%
  dplyr::summarise_all(mean) %>%
  tibble::column_to_rownames("Symbol") %>%
  as.matrix()

cat(sprintf("Final matrix: %d genes x %d samples\n", nrow(expr_final), ncol(expr_final)))

# ── 4. Save ─────────────────────────────────────────────────
saveRDS(list(
  expression = expr_final,
  metadata = pdat_sorted
), file = "GSE31348_processed.rds")

cat("Done: GSE31348_processed.rds saved.\n")
