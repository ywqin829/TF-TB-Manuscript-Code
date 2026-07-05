# ============================================================
# Shared Theme & Utility Functions
# ============================================================
# Source this file at the top of each figure script:
#   source("utils/theme_memory.R")
# ============================================================

library(ggplot2)
library(ggpubr)

# ── Publication Theme ───────────────────────────────────────
theme_pub <- theme_classic() +
  theme(
    strip.background = element_blank(),
    strip.text        = element_text(face = "bold", size = 12, color = "black"),
    axis.text.x       = element_text(color = "black", size = 12, face = "bold"),
    axis.text.y       = element_text(color = "black", size = 12, face = "bold"),
    axis.title        = element_text(face = "bold", size = 14),
    plot.title        = element_text(face = "bold", hjust = 0.5, size = 16),
    plot.subtitle     = element_text(hjust = 0.5, size = 14),
    legend.title      = element_text(face = "bold", size = 14),
    legend.text       = element_text(size = 12)
  )

# ── Color Palettes ──────────────────────────────────────────
# Timepoint colors (GSE40553: 6 stages)
col_6tp <- c(
  "LTB"    = "#F7ABE8",
  "TB_0M"  = "#FC4E07",
  "TB_0.5M" = "#E7B800",
  "TB_2M"  = "#AD07E3",
  "TB_6M"  = "#57A4FD",
  "TB_12M" = "#1B9E77"
)

# Treatment outcome colors
col_outcome <- c(
  "Successful_Treatment" = "#1976D2",
  "Poor_Outcome"         = "#D32F2F"
)

# GSE31348 timepoints (Blues palette)
col_31348 <- c(
  "Week 0"  = "#F7FBFF",
  "Week 1"  = "#DEEBF7",
  "Week 2"  = "#9ECAE1",
  "Week 4"  = "#6BAED6",
  "Week 26" = "#2171B5"
)

# MTB PBMC timepoints
col_mtb <- c("UT" = "#E69F00", "3hMTB" = "#D55E00", "24hMTB" = "#CC79A7")

# ── Statistical Comparisons ─────────────────────────────────
# GSE40553: compare each timepoint to TB_0M
comp_0M <- list(
  c("TB_0M", "LTB"),
  c("TB_0M", "TB_6M")
)

# GSE31348: compare each to Week 0
comp_w0 <- list(
  c("Week 0", "Week 1"),
  c("Week 0", "Week 2"),
  c("Week 0", "Week 4"),
  c("Week 0", "Week 26")
)

# ── Helper: Paired Wilcox Test with TIGIT-style plot ────────
# For pseudobulk-level paired comparisons across timepoints
plot_paired_timecourse <- function(data, gene_name, y_label = "Average Expression") {
  pd <- position_jitter(width = 0.1, seed = 42)
  ggplot(data, aes(x = timepoint, y = Mean_Expression)) +
    geom_boxplot(aes(fill = timepoint), outlier.shape = NA, alpha = 0.8,
                 width = 0.6, color = "black", linewidth = 0.6) +
    geom_line(aes(group = Donor), color = "gray75", alpha = 0.7,
              linetype = "11", linewidth = 0.35, position = pd) +
    geom_point(size = 2.5, shape = 21, fill = "white", color = "gray10",
               stroke = 1, position = pd) +
    facet_wrap(~ cell_type, scales = "free_y", ncol = 3) +
    scale_fill_manual(values = c("UT"="#D9D9D9", "3hMTB"="#4C78A8",
                                  "24hMTB"="#E45756")) +
    stat_compare_means(comparisons = list(c("UT", "3hMTB"), c("UT", "24hMTB")),
                       paired = TRUE, method = "wilcox.test",
                       label = "p.signif", tip.length = 0.01, size = 5) +
    labs(title = paste0(gene_name, " Expression Dynamics"),
         y = y_label, x = "") +
    theme_classic() +
    theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
          axis.text.x = element_text(size = 12, face = "bold", color = "black",
                                      angle = 45, hjust = 1),
          axis.text.y = element_text(size = 11, color = "black"),
          axis.title.y = element_text(size = 14, face = "bold"),
          strip.text = element_text(size = 14, face = "bold", color = "black"),
          strip.background = element_blank(),
          legend.position = "none")
}
