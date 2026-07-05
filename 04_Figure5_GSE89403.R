# ============================================================
# Figure 5: GSE89403 Outcome Validation — PAX5 & NFATC2
# ============================================================
# Dataset: GSE89403 (TB treatment outcome cohort)
#   - expression : 15879 genes x 312 samples (log2 expression)
#   - metadata   : 312 x 7 (TimePoint, Group, SubjectID, etc.)
#
# Timepoint mapping: 0_Baseline -> 0M, 1_Week1 -> W1,
#                    2_Week4 -> W4, 3_Week24 -> 6M
# Outcome groups   : Successful_Treatment, Poor_Outcome
#
# Fig 5A-C: PAX5 validation
#   5A: PAX5 activity (expression) trajectory across timepoints
#   5B: PAX5 target score (CD19/CD79A/FCER2/POU2AF1/BLK mean) trajectory
#   5C: PAX5 activity boxplot at 6M endpoint
#
# Fig 5D-F: NFATC2 (TBX21 surrogate) validation
#   5D: NFATC2 activity trajectory
#   5E: NFATC2 target score (IFNG/CD40LG mean) trajectory
#   5F: NFATC2 activity boxplot at 6M endpoint
# ============================================================

source("utils/theme_memory.R")

library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyr)

# ---- Configuration ----
data_dir   <- "C:/Users/yw_qi/Downloads"
output_dir <- "output"
outcome_colors <- col_outcome  # from theme_memory.R

# ---- Load Data ----
gse <- readRDS(file.path(data_dir, "GSE89403_Ready_for_Analysis.rds"))
expr_mat <- gse$expression   # 15879 x 312 matrix
meta     <- gse$metadata     # 312 x 7 data.frame

# ---- Prepare metadata: rename timepoints ----
timepoint_map <- c(
  "0_Baseline" = "0M",
  "1_Week1"   = "W1",
  "2_Week4"   = "W4",
  "3_Week24"  = "6M"
)
meta <- meta %>%
  mutate(
    TimeLabel  = recode(TimePoint, !!!timepoint_map),
    TimeLabel  = factor(TimeLabel, levels = c("0M", "W1", "W4", "6M")),
    Group      = factor(Group, levels = c("Successful_Treatment", "Poor_Outcome"))
  )
rownames(meta) <- rownames(gse$metadata)

# ---- Helper: extract gene expression for a sample set ----
get_expr <- function(gene, samples) {
  as.numeric(expr_mat[gene, samples])
}

# ---- Helper: compute target score (mean of multiple genes) ----
target_score <- function(genes, samples) {
  present <- genes[genes %in% rownames(expr_mat)]
  if (length(present) == 0) return(rep(NA, length(samples)))
  colMeans(expr_mat[present, samples, drop = FALSE], na.rm = TRUE)
}

# ============================================================
# Fig 5A: PAX5 Activity Trajectory
# ============================================================
cat("Plotting Fig 5A: PAX5 activity trajectory...\n")

pax5_expr <- data.frame(
  Sample   = colnames(expr_mat),
  Activity = get_expr("PAX5", colnames(expr_mat)),
  stringsAsFactors = FALSE
)
pax5_expr <- pax5_expr %>%
  left_join(meta[, c("TimeLabel", "Group")], by = "Sample")

p5a <- ggplot(pax5_expr, aes(x = TimeLabel, y = Activity, color = Group)) +
  geom_jitter(width = 0.15, size = 1.8, alpha = 0.5) +
  geom_smooth(aes(group = Group), method = "loess", se = TRUE,
              alpha = 0.15, linewidth = 0.8, span = 0.9) +
  stat_summary(aes(group = Group), fun = mean, geom = "line",
               linewidth = 1.2, alpha = 0.9) +
  stat_summary(aes(group = Group), fun = mean, geom = "point",
               size = 3, alpha = 0.9) +
  scale_color_manual(values = outcome_colors) +
  labs(title = "PAX5 Activity Trajectory",
       subtitle = "GSE89403",
       x = "Timepoint", y = "PAX5 Expression (log2)",
       color = "Outcome") +
  theme_pub +
  theme(legend.position = "top")

ggsave(file.path(output_dir, "Fig5A_PAX5_Activity_Trajectory.pdf"),
       p5a, width = 6, height = 5)

# ============================================================
# Fig 5B: PAX5 Target Score Trajectory
# ============================================================
cat("Plotting Fig 5B: PAX5 target score trajectory...\n")

pax5_targets <- c("CD19", "CD79A", "FCER2", "POU2AF1", "BLK")
pax5_ts <- data.frame(
  Sample      = colnames(expr_mat),
  TargetScore = target_score(pax5_targets, colnames(expr_mat)),
  stringsAsFactors = FALSE
)
pax5_ts <- pax5_ts %>%
  left_join(meta[, c("TimeLabel", "Group")], by = "Sample")

p5b <- ggplot(pax5_ts, aes(x = TimeLabel, y = TargetScore, color = Group)) +
  geom_jitter(width = 0.15, size = 1.8, alpha = 0.5) +
  geom_smooth(aes(group = Group), method = "loess", se = TRUE,
              alpha = 0.15, linewidth = 0.8, span = 0.9) +
  stat_summary(aes(group = Group), fun = mean, geom = "line",
               linewidth = 1.2, alpha = 0.9) +
  stat_summary(aes(group = Group), fun = mean, geom = "point",
               size = 3, alpha = 0.9) +
  scale_color_manual(values = outcome_colors) +
  labs(title = "PAX5 Target Score Trajectory",
       subtitle = "CD19 / CD79A / FCER2 / POU2AF1 / BLK",
       x = "Timepoint", y = "Mean Target Expression (log2)",
       color = "Outcome") +
  theme_pub +
  theme(legend.position = "top")

ggsave(file.path(output_dir, "Fig5B_PAX5_TargetScore_Trajectory.pdf"),
       p5b, width = 6, height = 5)

# ============================================================
# Fig 5C: PAX5 Activity Boxplot at 6M
# ============================================================
cat("Plotting Fig 5C: PAX5 activity boxplot at 6M...\n")

pax5_6m <- pax5_expr %>% filter(TimeLabel == "6M")

p5c <- ggplot(pax5_6m, aes(x = Group, y = Activity, fill = Group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.5,
               color = "black", linewidth = 0.6) +
  geom_jitter(width = 0.12, size = 2, alpha = 0.6, color = "gray30") +
  scale_fill_manual(values = outcome_colors) +
  stat_compare_means(method = "wilcox.test", label = "p.format",
                     label.x.npc = "center", size = 4.5) +
  labs(title = "PAX5 Activity at 6M Endpoint",
       subtitle = "GSE89403",
       x = "", y = "PAX5 Expression (log2)") +
  theme_pub +
  theme(legend.position = "none")

ggsave(file.path(output_dir, "Fig5C_PAX5_Activity_Boxplot_6M.pdf"),
       p5c, width = 4.5, height = 5)

# ============================================================
# Fig 5D: NFATC2 Activity Trajectory  (TBX21 surrogate)
# ============================================================
cat("Plotting Fig 5D: NFATC2 activity trajectory...\n")

nfatc2_expr <- data.frame(
  Sample   = colnames(expr_mat),
  Activity = get_expr("NFATC2", colnames(expr_mat)),
  stringsAsFactors = FALSE
)
nfatc2_expr <- nfatc2_expr %>%
  left_join(meta[, c("TimeLabel", "Group")], by = "Sample")

p5d <- ggplot(nfatc2_expr, aes(x = TimeLabel, y = Activity, color = Group)) +
  geom_jitter(width = 0.15, size = 1.8, alpha = 0.5) +
  geom_smooth(aes(group = Group), method = "loess", se = TRUE,
              alpha = 0.15, linewidth = 0.8, span = 0.9) +
  stat_summary(aes(group = Group), fun = mean, geom = "line",
               linewidth = 1.2, alpha = 0.9) +
  stat_summary(aes(group = Group), fun = mean, geom = "point",
               size = 3, alpha = 0.9) +
  scale_color_manual(values = outcome_colors) +
  labs(title = "NFATC2 Activity Trajectory",
       subtitle = "GSE89403 (TBX21 surrogate)",
       x = "Timepoint", y = "NFATC2 Expression (log2)",
       color = "Outcome") +
  theme_pub +
  theme(legend.position = "top")

ggsave(file.path(output_dir, "Fig5D_NFATC2_Activity_Trajectory.pdf"),
       p5d, width = 6, height = 5)

# ============================================================
# Fig 5E: NFATC2 Target Score Trajectory
# ============================================================
cat("Plotting Fig 5E: NFATC2 target score trajectory...\n")

nfatc2_targets <- c("IFNG", "CD40LG")
nfatc2_ts <- data.frame(
  Sample      = colnames(expr_mat),
  TargetScore = target_score(nfatc2_targets, colnames(expr_mat)),
  stringsAsFactors = FALSE
)
nfatc2_ts <- nfatc2_ts %>%
  left_join(meta[, c("TimeLabel", "Group")], by = "Sample")

# Handle any NA rows
nfatc2_ts <- nfatc2_ts %>% filter(!is.na(TargetScore))

p5e <- ggplot(nfatc2_ts, aes(x = TimeLabel, y = TargetScore, color = Group)) +
  geom_jitter(width = 0.15, size = 1.8, alpha = 0.5) +
  geom_smooth(aes(group = Group), method = "loess", se = TRUE,
              alpha = 0.15, linewidth = 0.8, span = 0.9) +
  stat_summary(aes(group = Group), fun = mean, geom = "line",
               linewidth = 1.2, alpha = 0.9) +
  stat_summary(aes(group = Group), fun = mean, geom = "point",
               size = 3, alpha = 0.9) +
  scale_color_manual(values = outcome_colors) +
  labs(title = "NFATC2 Target Score Trajectory",
       subtitle = "IFNG / CD40LG",
       x = "Timepoint", y = "Mean Target Expression (log2)",
       color = "Outcome") +
  theme_pub +
  theme(legend.position = "top")

ggsave(file.path(output_dir, "Fig5E_NFATC2_TargetScore_Trajectory.pdf"),
       p5e, width = 6, height = 5)

# ============================================================
# Fig 5F: NFATC2 Activity Boxplot at 6M
# ============================================================
cat("Plotting Fig 5F: NFATC2 activity boxplot at 6M...\n")

nfatc2_6m <- nfatc2_expr %>% filter(TimeLabel == "6M")

p5f <- ggplot(nfatc2_6m, aes(x = Group, y = Activity, fill = Group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.5,
               color = "black", linewidth = 0.6) +
  geom_jitter(width = 0.12, size = 2, alpha = 0.6, color = "gray30") +
  scale_fill_manual(values = outcome_colors) +
  stat_compare_means(method = "wilcox.test", label = "p.format",
                     label.x.npc = "center", size = 4.5) +
  labs(title = "NFATC2 Activity at 6M Endpoint",
       subtitle = "GSE89403",
       x = "", y = "NFATC2 Expression (log2)") +
  theme_pub +
  theme(legend.position = "none")

ggsave(file.path(output_dir, "Fig5F_NFATC2_Activity_Boxplot_6M.pdf"),
       p5f, width = 4.5, height = 5)

# ============================================================
cat("\nFigure 5 complete. All PDFs saved to:", output_dir, "\n")
# ============================================================
