# ============================================================
# prep_GSE19435.R — Additional Validation Cohort
# ============================================================
# Input:  GSE19435_series_matrix.txt.gz (from GEO)
# Output: GSE19435_processed.rds
# ============================================================

library(GEOquery)
library(limma)
library(dplyr)

options(timeout = 90000)

# ── 1. Load GEO ─────────────────────────────────────────────
gset <- getGEO(filename = "GSE19435_series_matrix.txt.gz")
pdata <- pData(gset)
exprs_data <- exprs(gset)

# ── 2. Define groups ────────────────────────────────────────
# GSE19435: HC (healthy), TB_0M (pretreatment), TB_2M, TB_12M
illness <- pdata$characteristics_ch1
pdata$Group <- NA
pdata$Group[grepl("healthy", illness, ignore.case = TRUE)] <- "HC"
pdata$Group[grepl("PTB|pulmonary tuberculosis", illness, ignore.case = TRUE) &
              grepl("0 month|pre", illness, ignore.case = TRUE)] <- "TB_0M"
pdata$Group[grepl("PTB|pulmonary tuberculosis", illness, ignore.case = TRUE) &
              grepl("2 month", illness, ignore.case = TRUE)] <- "TB_2M"
pdata$Group[grepl("PTB|pulmonary tuberculosis", illness, ignore.case = TRUE) &
              grepl("12 month", illness, ignore.case = TRUE)] <- "TB_12M"

pdata$Group <- factor(pdata$Group, levels = c("HC", "TB_0M", "TB_2M", "TB_12M"))
valid_samples <- rownames(pdata)[!is.na(pdata$Group)]
pdata <- pdata[valid_samples, ]
exprs_data <- exprs_data[, valid_samples]
exprs_data <- exprs_data[, match(rownames(pdata), colnames(exprs_data))]

# ── 3. Normalization (quantile) ─────────────────────────────
if (min(exprs_data) <= 0) exprs_data <- exprs_data + abs(min(exprs_data)) + 1
if (max(exprs_data) > 100) exprs_data <- log2(exprs_data + 1)

# ── 4. DEG ──────────────────────────────────────────────────
design <- model.matrix(~0 + pdata$Group)
colnames(design) <- levels(pdata$Group)
fit <- lmFit(exprs_data, design)
contrasts <- makeContrasts(TB_0M_vs_HC = TB_0M - HC, TB_12M_vs_TB_0M = TB_12M - TB_0M, levels = design)
fit2 <- contrasts.fit(fit, contrasts)
fit2 <- eBayes(fit2)

de_hc <- topTable(fit2, coef = "TB_0M_vs_HC", number = Inf)
de_tx <- topTable(fit2, coef = "TB_12M_vs_TB_0M", number = Inf)

saveRDS(list(expression = exprs_data, metadata = pdata, de_hc = de_hc, de_tx = de_tx),
        file = "GSE19435_processed.rds")
cat("Done: GSE19435_processed.rds saved.\n")
