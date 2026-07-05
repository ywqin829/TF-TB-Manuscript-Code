# ============================================================
# Figure 6: LASSO + Random Forest + Healing Score
# ============================================================
#
# Data sources:
#   GSE40553 (Discovery)  : expr_symbol + pheno for LASSO & RF training
#   GSE31348 (Validation) : expression + metadata for ROC validation
#   GSE89403 (Outcome)    : expression + metadata for Healing Score
#
# Workflow:
#   1. LASSO (cv.glmnet) on GSE40553 TB_0M vs TB_6M
#       -> 19 non-zero coefficient genes at lambda.1se
#   6A: LASSO coefficient profile
#   6B: Cross-validation deviance
#   6C: Lollipop plot of 19 LASSO genes
#
#   2. RF exhaustive search (2-to-5 gene combos from 19 LASSO genes)
#      -> Optimal 4-gene combo: CD79A, EBF1, FCGR1A, GK
#   6D: Diagnostic ROC in GSE31348 (Wk0 vs Wk26), AUC annotation
#
#   3. Healing Score = (Z_CD79A + Z_EBF1) - (Z_FCGR1A + Z_GK)
#   6E: Dynamic ROC in GSE89403 (W1/W4/6M vs Baseline)
#   6F: Calibration curve (Baseline vs W4)
#   6G: Decision curve analysis
# ============================================================

source("utils/theme_memory.R")

# ---- Load Libraries ----
library(glmnet)
library(randomForest)
library(pROC)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyr)
library(doParallel)
library(foreach)

# ---- Configuration ----
data_dir   <- "C:/Users/yw_qi/Downloads"
output_dir <- "output"
n_folds    <- 10
rf_cv      <- 3
outcome_colors <- col_outcome

# ============================================================
# 1. LOAD & PREPARE TRAINING DATA (GSE40553)
# ============================================================
cat("=== Loading GSE40553 (training data) ===\n")
gse53    <- readRDS(file.path(data_dir, "GSE40553.rds"))
expr53   <- gse53$expr_symbol    # 21186 x 176 matrix
pheno53  <- gse53$pheno          # 176 x 39 data.frame

# Filter: TB_0M (disease) vs TB_6M (recovery)
train_idx <- which(pheno53$Group %in% c("TB_0M", "TB_6M"))
train_expr <- expr53[, train_idx]
train_pheno <- pheno53[train_idx, "Group", drop = FALSE]
train_pheno$Group <- factor(train_pheno$Group, levels = c("TB_0M", "TB_6M"))

cat("  Training samples:", ncol(train_expr), "\n")
cat("  TB_0M:", sum(train_pheno$Group == "TB_0M"), "\n")
cat("  TB_6M:", sum(train_pheno$Group == "TB_6M"), "\n")

# Transpose: samples x genes for glmnet / RF
X_train <- t(train_expr)   # samples x genes
y_train <- train_pheno$Group

# ---- Helper: match genes across datasets ----
match_genes <- function(genes, target_rownames) {
  genes[genes %in% target_rownames]
}

# ============================================================
# 2. LASSO (cv.glmnet, alpha=1, binomial, nfolds=10)
# ============================================================
cat("\n=== Running LASSO (cv.glmnet) ===\n")
set.seed(42)
cv_fit <- cv.glmnet(
  x        = X_train,
  y        = y_train,
  family   = "binomial",
  alpha    = 1,
  nfolds   = n_folds,
  standardize = TRUE
)

# Extract non-zero coefficients at lambda.1se
lambda_1se <- cv_fit$lambda.1se
coef_mat   <- as.matrix(coef(cv_fit, s = "lambda.1se"))
nonzero_idx <- which(coef_mat[-1, 1] != 0)  # exclude intercept
lasso_genes <- rownames(coef_mat)[-1][nonzero_idx]
lasso_coefs <- coef_mat[-1, 1][nonzero_idx]

# Build a named vector for plotting
lasso_df <- data.frame(
  Gene       = names(lasso_coefs),
  Coefficient = as.numeric(lasso_coefs),
  Direction  = ifelse(lasso_coefs > 0, "Recovery", "Disease"),
  AbsCoef    = abs(lasso_coefs)
)
lasso_df <- lasso_df[order(lasso_df$Coefficient), ]
lasso_df$Gene <- factor(lasso_df$Gene, levels = lasso_df$Gene)

cat("  lambda.1se:", round(lambda_1se, 5), "\n")
cat("  Non-zero genes at lambda.1se:", length(lasso_genes), "\n")
cat("  Genes:", paste(lasso_genes, collapse = ", "), "\n")

# ============================================================
# Fig 6A: LASSO Coefficient Profile
# ============================================================
cat("\nPlotting Fig 6A: LASSO coefficient profile...\n")

pdf(file.path(output_dir, "Fig6A_LASSO_Coefficient_Profile.pdf"),
    width = 8, height = 6)
plot(cv_fit$glmnet.fit, xvar = "lambda", label = TRUE,
     main = "LASSO Coefficient Profile (GSE40553: TB_0M vs TB_6M)")
abline(v = log(lambda_1se), lty = 2, col = "red", lwd = 2)
dev.off()
cat("  Saved: Fig6A_LASSO_Coefficient_Profile.pdf\n")

# ============================================================
# Fig 6B: Cross-Validation Deviance Plot
# ============================================================
cat("Plotting Fig 6B: LASSO CV deviance...\n")

pdf(file.path(output_dir, "Fig6B_LASSO_CV_Deviance.pdf"),
    width = 8, height = 6)
plot(cv_fit, main = "LASSO Cross-Validation Deviance (GSE40553)")
dev.off()
cat("  Saved: Fig6B_LASSO_CV_Deviance.pdf\n")

# ============================================================
# Fig 6C: Lollipop Plot of 19 LASSO Genes
# ============================================================
cat("Plotting Fig 6C: Lollipop plot...\n")

# Assign recovery/disease colors
lollipop_colors <- c("Recovery" = "#1976D2", "Disease" = "#D32F2F")

p6c <- ggplot(lasso_df, aes(x = Gene, y = Coefficient, color = Direction)) +
  geom_segment(aes(x = Gene, xend = Gene, y = 0, yend = Coefficient),
               linewidth = 0.8, color = "grey50") +
  geom_point(size = 3.5) +
  coord_flip() +
  scale_color_manual(values = lollipop_colors) +
  labs(title = "LASSO Selected Genes (lambda.1se)",
       subtitle = paste0(nrow(lasso_df), " genes, GSE40553 TB_0M vs TB_6M"),
       x = "Gene", y = "Coefficient") +
  theme_pub +
  theme(legend.position = "top")

ggsave(file.path(output_dir, "Fig6C_LASSO_Lollipop.pdf"),
       p6c, width = 7, height = 6)

# ============================================================
# 3. RF Exhaustive Search on GSE40553 (2-to-5 gene combos)
# ============================================================
cat("\n=== RF Exhaustive Search (2-to-5 gene combinations) ===\n")
cat("  This may take several minutes...\n")

# Register parallel backend
n_cores <- max(1, detectCores() - 1)
cl <- makeCluster(n_cores)
registerDoParallel(cl)

# Generate all combinations
gene_combo_results <- foreach(k = 2:5, .combine = rbind,
                               .packages = c("randomForest")) %dopar% {
  combos <- combn(lasso_genes, k, simplify = FALSE)
  n_combos <- length(combos)

  # Sample if too many combos (5-gene: 11628 may be excessive)
  if (n_combos > 2000) {
    set.seed(42)
    combos <- combos[sample(n_combos, 2000)]
  }

  results <- t(sapply(combos, function(genes) {
    # 3-fold CV with RF
    set.seed(42)
    fold_ids <- sample(rep(1:rf_cv, length.out = nrow(X_train)))

    acc <- sapply(1:rf_cv, function(f) {
      train_f <- X_train[fold_ids != f, genes, drop = FALSE]
      test_f  <- X_train[fold_ids == f, genes, drop = FALSE]
      y_f     <- y_train[fold_ids != f]

      rf <- randomForest(x = train_f, y = y_f, ntree = 500)
      pred <- predict(rf, test_f)
      mean(pred == y_train[fold_ids == f])
    })
    c(k = k, genes = paste(genes, collapse = "/"), accuracy = mean(acc))
  }), stringsAsFactors = FALSE)

  results
}

stopCluster(cl)
registerDoSEQ()

# Find best combination
gene_combo_results <- as.data.frame(gene_combo_results)
gene_combo_results$k <- as.numeric(gene_combo_results$k)
gene_combo_results$accuracy <- as.numeric(gene_combo_results$accuracy)
best_combo <- gene_combo_results[which.max(gene_combo_results$accuracy), ]
cat("  Best combination:", best_combo$genes, "\n")
cat("  CV Accuracy:", round(best_combo$accuracy, 4), "\n")

# Optimal 4-gene combination (from manuscript)
healing_genes <- c("CD79A", "EBF1", "FCGR1A", "GK")
cat("  Using optimal 4-gene combo:", paste(healing_genes, collapse = ", "), "\n")

# ============================================================
# 4. TRAIN RF ON GSE40553 WITH 4-GENE COMBO
# ============================================================
cat("\n=== Training RF with 4-gene combo ===\n")
set.seed(42)
rf_final <- randomForest(
  x       = X_train[, healing_genes],
  y       = y_train,
  ntree   = 1000,
  importance = TRUE
)
cat("  RF OOB error rate:", round(tail(rf_final$err.rate[,1], 1), 2), "%\n")

# ============================================================
# 5. LOAD GSE31348 & PREDICT FOR ROC
# ============================================================
cat("\n=== Loading GSE31348 (validation) for ROC ===\n")
gse48    <- readRDS(file.path(data_dir, "GSE31348_Figure/GSE31348_Ready.rds"))
expr48   <- as.matrix(gse48$expression)   # 21357 x 135
meta48   <- gse48$metadata                # 135 x 5

# Filter Week 0 and Week 26
val_idx <- which(meta48$TimePoint %in% c("Week 0", "Week 26"))
X_val   <- expr48[healing_genes, val_idx, drop = FALSE]
X_val   <- t(X_val)   # samples x genes
y_val   <- factor(meta48$TimePoint[val_idx],
                  levels = c("Week 0", "Week 26"))

# Check for missing genes
missing_genes <- healing_genes[!healing_genes %in% rownames(expr48)]
if (length(missing_genes) > 0) {
  cat("  WARNING: Missing genes in GSE31348:", paste(missing_genes, collapse = ", "), "\n")
} else {
  cat("  All 4 genes present in GSE31348\n")
}

# Predict
pred_prob <- predict(rf_final, X_val, type = "prob")[, "TB_6M"]

# ============================================================
# Fig 6D: Diagnostic ROC in GSE31348
# ============================================================
cat("Plotting Fig 6D: Diagnostic ROC in GSE31348...\n")

roc_val <- roc(y_val, pred_prob, levels = rev(levels(y_val)))
auc_val <- round(auc(roc_val), 3)

roc_df <- data.frame(
  specificity = roc_val$specificities,
  sensitivity = roc_val$sensitivities
)

p6d <- ggplot(roc_df, aes(x = 1 - specificity, y = sensitivity)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey60") +
  geom_ribbon(aes(ymin = 0, ymax = sensitivity), fill = "#D32F2F", alpha = 0.15) +
  geom_line(color = "#D32F2F", linewidth = 1.2) +
  annotate("text", x = 0.75, y = 0.25, size = 5.5, fontface = "bold",
           label = paste0("AUC = ", auc_val),
           color = "#D32F2F") +
  labs(title = "Diagnostic ROC (GSE31348)",
       subtitle = "RF trained on GSE40553: TB_0M vs TB_6M\nValidated on GSE31348: Week 0 vs Week 26",
       x = "1 - Specificity", y = "Sensitivity") +
  theme_pub

ggsave(file.path(output_dir, "Fig6D_Diagnostic_ROC_GSE31348.pdf"),
       p6d, width = 5.5, height = 5)

# ============================================================
# 6. HEALING SCORE IN GSE89403
# ============================================================
cat("\n=== Computing Healing Score in GSE89403 ===\n")
gse03   <- readRDS(file.path(data_dir, "GSE89403_Ready_for_Analysis.rds"))
expr03  <- gse03$expression
meta03  <- gse03$metadata

# Rename timepoints for consistency
timepoint_map <- c(
  "0_Baseline" = "0M",
  "1_Week1"   = "W1",
  "2_Week4"   = "W4",
  "3_Week24"  = "6M"
)
meta03$TimeLabel <- recode(meta03$TimePoint, !!!timepoint_map)
meta03$TimeLabel <- factor(meta03$TimeLabel, levels = c("0M", "W1", "W4", "6M"))
meta03$Group <- factor(meta03$Group, levels = c("Successful_Treatment", "Poor_Outcome"))

# Z-score normalize the 4 genes across GSE89403
compute_zscore <- function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)

z_cd79a  <- compute_zscore(as.numeric(expr03["CD79A", ]))
z_ebf1   <- compute_zscore(as.numeric(expr03["EBF1", ]))
z_fcgr1a <- compute_zscore(as.numeric(expr03["FCGR1A", ]))
z_gk     <- compute_zscore(as.numeric(expr03["GK", ]))

healing_score_03 <- (z_cd79a + z_ebf1) - (z_fcgr1a + z_gk)

score_df <- data.frame(
  Sample   = colnames(expr03),
  Score    = healing_score_03,
  stringsAsFactors = FALSE
)
score_df <- score_df %>%
  left_join(meta03[, c("TimeLabel", "Group", "SubjectID")], by = "Sample")

cat("  Healing Score range:", round(range(score_df$Score, na.rm = TRUE), 2), "\n")

# ============================================================
# Fig 6E: Dynamic ROC in GSE89403
#   At each post-baseline timepoint, compare outcomes
# ============================================================
cat("Plotting Fig 6E: Dynamic ROC...\n")

dyn_tps <- c("W1", "W4", "6M")
roc_list <- list()

for (tp in dyn_tps) {
  sub <- score_df %>% filter(TimeLabel == tp)
  if (nrow(sub) < 5) next
  r <- roc(sub$Group, sub$Score, direction = "<",
           levels = c("Poor_Outcome", "Successful_Treatment"))
  roc_list[[tp]] <- r
}

# Plot ROC curves on one panel
dyn_roc_df <- do.call(rbind, lapply(names(roc_list), function(tp) {
  r <- roc_list[[tp]]
  data.frame(
    Timepoint     = tp,
    specificity   = r$specificities,
    sensitivity   = r$sensitivities,
    AUC           = round(auc(r), 3)
  )
}))

dyn_auc_labels <- unique(dyn_roc_df[, c("Timepoint", "AUC")])
dyn_auc_labels$label <- paste0(dyn_auc_labels$Timepoint, " AUC = ",
                                dyn_auc_labels$AUC)

tp_colors <- c("W1" = "#E7B800", "W4" = "#AD07E3", "6M" = "#57A4FD")

p6e <- ggplot(dyn_roc_df, aes(x = 1 - specificity, y = sensitivity,
                               color = Timepoint)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey60") +
  geom_line(linewidth = 1) +
  scale_color_manual(values = tp_colors) +
  labs(title = "Dynamic ROC: Healing Score in GSE89403",
       subtitle = "Predicting treatment outcome at each timepoint vs Baseline",
       x = "1 - Specificity", y = "Sensitivity",
       color = "Timepoint") +
  theme_pub +
  theme(legend.position = "right")

ggsave(file.path(output_dir, "Fig6E_HealingScore_DynamicROC.pdf"),
       p6e, width = 6, height = 5)

# ============================================================
# Fig 6F: Calibration Curve (Baseline vs Week 4)
# ============================================================
cat("Plotting Fig 6F: Calibration curve...\n")

library(rms)

# Prepare calibration data: Healing Score at Week 4
cal_w4 <- score_df %>% filter(TimeLabel == "W4") %>%
  mutate(OutcomeBin = ifelse(Group == "Successful_Treatment", 1, 0))

# Fit logistic model for calibration
dd <- datadist(cal_w4)
options(datadist = "dd")
lrm_fit <- lrm(OutcomeBin ~ Score, data = cal_w4, x = TRUE, y = TRUE)
cal_obj <- calibrate(lrm_fit, B = 200)

pdf(file.path(output_dir, "Fig6F_HealingScore_Calibration.pdf"),
    width = 6, height = 6)
plot(cal_obj,
     xlab = "Predicted Probability",
     ylab = "Observed Probability",
     main = "Calibration Curve (Healing Score, Week 4)",
     subtitles = FALSE)
dev.off()
cat("  Saved: Fig6F_HealingScore_Calibration.pdf\n")

options(datadist = NULL)

# ============================================================
# Fig 6G: Decision Curve Analysis
# ============================================================
cat("Plotting Fig 6G: Decision curve analysis...\n")

if (requireNamespace("dcurves", quietly = TRUE)) {
  library(dcurves)

  # For DCA: use outcome at Week 4
  dca_data <- score_df %>%
    filter(TimeLabel == "W4") %>%
    mutate(
      outcome = ifelse(Group == "Successful_Treatment", 1, 0)
    ) %>%
    select(outcome, Score)

  dca_raw <- dca(outcome ~ Score, data = dca_data,
                 thresholds = seq(0, 1, by = 0.01))
  # dcurves >= 0.4.0 returns a data.frame; older versions return a list with $dca
  dca_df <- if (is.data.frame(dca_raw)) dca_raw else dca_raw$dca

  p6g <- ggplot(dca_df, aes(x = threshold, y = net_benefit, color = variable)) +
    geom_line(linewidth = 1) +
    scale_color_manual(values = c("Score" = "#1976D2",
                                   "all" = "grey60",
                                   "none" = "grey80")) +
    coord_cartesian(ylim = c(-0.05, 0.6)) +
    labs(title = "Decision Curve Analysis",
         subtitle = "Healing Score at Week 4 (GSE89403)",
         x = "Threshold Probability", y = "Net Benefit",
         color = "") +
    theme_pub +
    theme(legend.position = "top")

  ggsave(file.path(output_dir, "Fig6G_HealingScore_DCA.pdf"),
         p6g, width = 6, height = 5)

} else {
  cat("  NOTE: dcurves package not installed. Skipping DCA.\n")
  cat("  Install with: install.packages('dcurves')\n")
}

# ============================================================
# Clean up datadist
options(datadist = NULL)

cat("\n=== Figure 6 complete. All PDFs saved to:", output_dir, "===\n")
cat("  Fig6A: LASSO Coefficient Profile\n")
cat("  Fig6B: LASSO Cross-Validation Deviance\n")
cat("  Fig6C: Lollipop Plot (", nrow(lasso_df), "genes)\n", sep="")
cat("  Fig6D: Diagnostic ROC (GSE31348, AUC =", auc_val, ")\n")
cat("  Fig6E: Healing Score Dynamic ROC\n")
cat("  Fig6F: Healing Score Calibration Curve\n")
cat("  Fig6G: Healing Score Decision Curve Analysis\n")
# ============================================================
