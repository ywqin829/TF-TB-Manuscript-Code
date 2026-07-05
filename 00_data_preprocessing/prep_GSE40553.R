# ============================================================
# prep_GSE40553.R — Discovery Cohort Data Processing
# ============================================================
# Input:  GSE40553_series_matrix.txt.gz (from GEO)
# Output: GSE40553_processed.rds (list with exprs_symbol, pdata,
#         genes_core_recovery, genes_core_disease, DEG results)
# ============================================================

library(GEOquery)
library(limma)
library(dplyr)
library(AnnotationDbi)
library(illuminaHumanv4.db)
set.seed(42)

options(timeout = 90000)

# ── 1. Load GEO ─────────────────────────────────────────────
gset <- getGEO(filename = "GSE40553_series_matrix.txt.gz")
pdata <- pData(gset)
exprs_data <- exprs(gset)

# ── 2. Define groups ────────────────────────────────────────
illness <- pdata$characteristics_ch1
status  <- pdata$characteristics_ch1.4
pdata$Group <- NA
pdata$Group[grepl("LTB", illness) & grepl("untreated latent TB", status)] <- "LTB"
pdata$Group[grepl("PTB", illness) & grepl("pre-treatment", status)] <- "TB_0M"
pdata$Group[grepl("PTB", illness) & grepl("2 weeks post treatment", status)] <- "TB_0.5M"
pdata$Group[grepl("PTB", illness) & grepl("2 months post treatment", status)] <- "TB_2M"
pdata$Group[grepl("PTB", illness) & grepl("6 months post treatment", status)] <- "TB_6M"
pdata$Group[grepl("PTB", illness) & grepl("12 months post treatment", status)] <- "TB_12M"
pdata$Group <- factor(pdata$Group, levels = c("LTB", "TB_0M", "TB_0.5M", "TB_2M", "TB_6M", "TB_12M"))

# ── 3. Remove NA group samples ──────────────────────────────
valid_samples <- rownames(pdata)[!is.na(pdata$Group)]
pdata <- pdata[valid_samples, ]
exprs_data <- exprs_data[, valid_samples]
exprs_data <- exprs_data[, match(rownames(pdata), colnames(exprs_data))]
stopifnot(all(colnames(exprs_data) == rownames(pdata)))

# ── 4. NA filtering and filling ─────────────────────────────
expr_matrix <- exprs_data
keep_na <- rowSums(is.na(expr_matrix)) < ncol(expr_matrix) * 0.3
expr_matrix <- expr_matrix[keep_na, ]
for (i in seq_len(ncol(expr_matrix))) {
  idx <- is.na(expr_matrix[, i])
  if (any(idx)) expr_matrix[idx, i] <- min(expr_matrix[, i], na.rm = TRUE)
}

# ── 5. Log2 check ───────────────────────────────────────────
if (min(expr_matrix) <= 0) expr_matrix <- expr_matrix + abs(min(expr_matrix)) + 1
if (max(expr_matrix) > 100) expr_matrix <- log2(expr_matrix + 1)

# ── 6. Quantile normalization ───────────────────────────────
expr_matrix <- normalizeBetweenArrays(expr_matrix, method = "quantile")

# ── 7. Expression filtering ─────────────────────────────────
threshold_val  <- quantile(expr_matrix, 0.25)
min_group_size <- min(table(pdata$Group))
keep_probe <- rowSums(expr_matrix > threshold_val) >= min_group_size
expr_matrix <- expr_matrix[keep_probe, ]

# ── 8. Probe → Gene (illuminaHumanv4) ──────────────────────
anno <- AnnotationDbi::select(
  illuminaHumanv4.db,
  keys = rownames(expr_matrix),
  columns = c("SYMBOL"),
  keytype = "PROBEID"
) %>% filter(!is.na(SYMBOL)) %>% distinct(PROBEID, .keep_all = TRUE)

expr_df <- as.data.frame(expr_matrix) %>% mutate(PROBEID = rownames(.))
merged <- inner_join(anno, expr_df, by = "PROBEID")
expr_cols <- colnames(expr_matrix)

best_probe <- merged %>%
  group_by(SYMBOL) %>%
  slice_max(order_by = rowMeans(across(all_of(expr_cols))), n = 1) %>%
  ungroup()

exprs_symbol <- as.matrix(best_probe[, expr_cols])
rownames(exprs_symbol) <- best_probe$SYMBOL

cat(sprintf("Final matrix: %d genes x %d samples\n", nrow(exprs_symbol), ncol(exprs_symbol)))

# ── 9. DEG analysis (limma) ─────────────────────────────────
design <- model.matrix(~0 + pdata$Group)
colnames(design) <- levels(pdata$Group)
contrasts <- makeContrasts(
  Disease = TB_0M - LTB,
  Recovery = TB_6M - TB_0M,
  levels = design
)

fit <- lmFit(exprs_symbol, design)
fit2 <- contrasts.fit(fit, contrasts)
fit2 <- eBayes(fit2)

Disease <- topTable(fit2, coef = "Disease", number = Inf, adjust.method = "BH")
Recovery <- topTable(fit2, coef = "Recovery", number = Inf, adjust.method = "BH")

# ── 10. Define Transcriptional Reversal genes ──────────────
# Restoration Core: down in disease (logFC<0) AND up after treatment (logFC>0)
genes_disease_down <- Disease %>% filter(logFC < 0, adj.P.Val < 0.05) %>% rownames()
genes_disease_up   <- Disease %>% filter(logFC > 0, adj.P.Val < 0.05) %>% rownames()
genes_treat_up     <- Recovery %>% filter(logFC > 0, adj.P.Val < 0.05) %>% rownames()
genes_treat_down   <- Recovery %>% filter(logFC < 0, adj.P.Val < 0.05) %>% rownames()

genes_core_recovery <- intersect(genes_disease_down, genes_treat_up)   # 169
genes_core_disease  <- intersect(genes_disease_up, genes_treat_down)     # 171

cat(sprintf("Restoration Core: %d genes\n", length(genes_core_recovery)))
cat(sprintf("Disease Drivers: %d genes\n", length(genes_core_disease)))

# ── 11. Save ────────────────────────────────────────────────
saveRDS(list(
  exprs_symbol = exprs_symbol,
  pdata = pdata,
  Disease = Disease,
  Recovery = Recovery,
  genes_core_recovery = genes_core_recovery,
  genes_core_disease = genes_core_disease,
  genes_disease_down = genes_disease_down,
  genes_disease_up = genes_disease_up,
  genes_treat_up = genes_treat_up,
  genes_treat_down = genes_treat_down
), file = "GSE40553_processed.rds")

cat("Done: GSE40553_processed.rds saved.\n")
