# ================================================================
# 01_Figure2_GSE40553.R
# Figure 2: Discovery cohort (GSE40553) —
#   Transcriptional reversal, temporal clustering, TF activity
# ================================================================
# Datasets:  GSE40553 (Illumina HumanHT-12 v4, 176 samples)
# Panels:    A — Volcano (TB_0M vs LTB)
#            B — Volcano (TB_6M vs TB_0M)
#            C — Venn (Transcriptional Reversal)
#            D — Mfuzz temporal clustering
#            E — GO BP enrichment (compareCluster)
#            F — VIPER TF activity heatmap
#            G — Fisher TF enrichment bubble plot
# ================================================================

# ---- 0. Setup ----
source("utils/theme_memory.R")

library(EnhancedVolcano)
library(limma)
library(ggVennDiagram)
library(patchwork)
library(Mfuzz)
library(clusterProfiler)
library(org.Hs.eg.db)
library(decoupleR)
library(dorothea)
library(pheatmap)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(ggrepel)
library(matrixStats)
library(RColorBrewer)

dir.create("output", showWarnings = FALSE)

# ================================================================
# DATA LOADING
# ================================================================
# The preprocessed GSE40553 RDS contains:
#   exprs_symbol : 21186 genes x 176 samples — gene-level expression
#   pheno        : 176 x 39 metadata (contains "Group" column)
#   group        : factor(6 levels: LTB/TB_0M/TB_0.5M/TB_2M/TB_6M/TB_12M)
# All align by position (index 1:176).

rds_path <- file.path("..", "GSE40553.rds")
if (!file.exists(rds_path)) {
  rds_path <- "C:/Users/yw_qi/Downloads/GSE40553.rds"
}
obj <- readRDS(rds_path)

exprs_symbol <- obj$expr_symbol
pdata <- obj$pheno
group <- obj$group
rm(obj)

# Sanity check
stopifnot(ncol(exprs_symbol) == length(group))
stopifnot(all(as.character(pdata$Group) == as.character(group)))

cat(sprintf("Loaded GSE40553: %d genes x %d samples\n", nrow(exprs_symbol), ncol(exprs_symbol)))
cat("Group distribution:\n")
print(table(group))


# ================================================================
# ---- RUN ONCE: DEG computation with limma ----
# Uncomment and run this block once to generate deg files.
# After running, subsequent figure generation loads the saved CSV.
# ================================================================
if (FALSE) {
  # ---- limma: Disease (TB_0M vs LTB) ----
  design_disease <- model.matrix(~ 0 + factor(group, levels = c("LTB", "TB_0M")))
  colnames(design_disease) <- c("LTB", "TB_0M")
  keep_d <- group %in% c("LTB", "TB_0M")
  fit_d <- lmFit(exprs_symbol[, keep_d], design_disease[keep_d, ])
  cont_d <- makeContrasts(TB_0M_vs_LTB = TB_0M - LTB, levels = design_disease)
  fit_d <- contrasts.fit(fit_d, cont_d)
  fit_d <- eBayes(fit_d)
  deg_disease <- topTable(fit_d, coef = 1, number = Inf)
  write.csv(deg_disease, "output/DEG_disease_TB0M_vs_LTB.csv", row.names = TRUE)

  # ---- limma: Recovery (TB_6M vs TB_0M) ----
  design_recovery <- model.matrix(~ 0 + factor(group, levels = c("TB_0M", "TB_6M")))
  colnames(design_recovery) <- c("TB_0M", "TB_6M")
  keep_r <- group %in% c("TB_0M", "TB_6M")
  fit_r <- lmFit(exprs_symbol[, keep_r], design_recovery[keep_r, ])
  cont_r <- makeContrasts(TB_6M_vs_TB0M = TB_6M - TB_0M, levels = design_recovery)
  fit_r <- contrasts.fit(fit_r, cont_r)
  fit_r <- eBayes(fit_r)
  deg_recovery <- topTable(fit_r, coef = 1, number = Inf)
  write.csv(deg_recovery, "output/DEG_recovery_TB6M_vs_TB0M.csv", row.names = TRUE)

  cat("DEG files saved to output/\n")
}

# ---- Load pre-computed DEGs (or create them inline) ----
if (file.exists("output/DEG_disease_TB0M_vs_LTB.csv")) {
  deg_disease <- read.csv("output/DEG_disease_TB0M_vs_LTB.csv", row.names = 1)
} else {
  # Inline fallback: run limma on the fly
  design_disease <- model.matrix(~ 0 + factor(group, levels = c("LTB", "TB_0M")))
  colnames(design_disease) <- c("LTB", "TB_0M")
  keep_d <- group %in% c("LTB", "TB_0M")
  fit_d <- lmFit(exprs_symbol[, keep_d], design_disease[keep_d, ])
  cont_d <- makeContrasts(TB_0M_vs_LTB = TB_0M - LTB, levels = design_disease)
  fit_d <- contrasts.fit(fit_d, cont_d)
  fit_d <- eBayes(fit_d)
  deg_disease <- topTable(fit_d, coef = 1, number = Inf)
  write.csv(deg_disease, "output/DEG_disease_TB0M_vs_LTB.csv", row.names = TRUE)
}

if (file.exists("output/DEG_recovery_TB6M_vs_TB0M.csv")) {
  deg_recovery <- read.csv("output/DEG_recovery_TB6M_vs_TB0M.csv", row.names = 1)
} else {
  design_recovery <- model.matrix(~ 0 + factor(group, levels = c("TB_0M", "TB_6M")))
  colnames(design_recovery) <- c("TB_0M", "TB_6M")
  keep_r <- group %in% c("TB_0M", "TB_6M")
  fit_r <- lmFit(exprs_symbol[, keep_r], design_recovery[keep_r, ])
  cont_r <- makeContrasts(TB_6M_vs_TB0M = TB_6M - TB_0M, levels = design_recovery)
  fit_r <- contrasts.fit(fit_r, cont_r)
  fit_r <- eBayes(fit_r)
  deg_recovery <- topTable(fit_r, coef = 1, number = Inf)
  write.csv(deg_recovery, "output/DEG_recovery_TB6M_vs_TB0M.csv", row.names = TRUE)
}

# ---- Define transcriptional reversal gene sets ----
# Filter: at least one group mean > 0.1 to avoid low-expression artifacts
mean_ltb <- rowMeans(exprs_symbol[, group == "LTB"])
mean_tb0 <- rowMeans(exprs_symbol[, group == "TB_0M"])
mean_tb6 <- rowMeans(exprs_symbol[, group == "TB_6M"])
keep_expr <- mean_tb0 > 0.1 | mean_ltb > 0.1 | mean_tb6 > 0.1

# DEG significance thresholds
sig_disease <- with(deg_disease, adj.P.Val < 0.05 & abs(logFC) > 1)
sig_recovery <- with(deg_recovery, adj.P.Val < 0.05 & abs(logFC) > 1)

# Restoration Core: Down in TB_0M (logFC < 0) ∩ Up after Tx (logFC > 0)
genes_down_in_TB <- rownames(deg_disease)[sig_disease & deg_disease$logFC < 0]
genes_up_after_Tx <- rownames(deg_recovery)[sig_recovery & deg_recovery$logFC > 0]
genes_core_recovery <- intersect(genes_down_in_TB, genes_up_after_Tx)
genes_core_recovery <- intersect(genes_core_recovery, rownames(exprs_symbol)[keep_expr])

# Disease Drivers: Up in TB_0M (logFC > 0) ∩ Down after Tx (logFC < 0)
genes_up_in_TB <- rownames(deg_disease)[sig_disease & deg_disease$logFC > 0]
genes_down_after_Tx <- rownames(deg_recovery)[sig_recovery & deg_recovery$logFC < 0]
genes_core_disease <- intersect(genes_up_in_TB, genes_down_after_Tx)
genes_core_disease <- intersect(genes_core_disease, rownames(exprs_symbol)[keep_expr])

cat(sprintf("Restoration Core (Down in TB ∩ Up after Tx): %d genes\n",
            length(genes_core_recovery)))
cat(sprintf("Disease Drivers (Up in TB ∩ Down after Tx):  %d genes\n",
            length(genes_core_disease)))

# Save gene lists for downstream use
writeLines(genes_core_recovery, "output/genes_core_recovery.txt")
writeLines(genes_core_disease, "output/genes_core_disease.txt")


# ================================================================
# Fig 2A — Volcano: TB_0M vs LTB (Disease hit)
# ================================================================
cat("Generating Fig 2A: Volcano — TB_0M vs LTB\n")

# Label top 20 genes by absolute fold change (within significant)
deg_disease$Gene <- rownames(deg_disease)
deg_disease$Label <- ""
top_n <- 20
top_idx <- order(abs(deg_disease$logFC) * (deg_disease$adj.P.Val < 0.05),
                 decreasing = TRUE)[1:min(top_n, nrow(deg_disease))]
top_idx <- top_idx[!is.na(top_idx)]
deg_disease$Label[top_idx] <- deg_disease$Gene[top_idx]

p2a <- EnhancedVolcano(deg_disease,
                       lab = deg_disease$Label,
                       x = "logFC",
                       y = "adj.P.Val",
                       title = "Disease: TB_0M vs LTB",
                       subtitle = "",
                       caption = "",
                       pCutoff = 0.05,
                       FCcutoff = 1,
                       labSize = 4.5,
                       pointSize = 0.6,
                       colAlpha = 0.8,
                       legendPosition = "top",
                       legendLabels = c("NS", expression(log[2]~FC), "FDR",
                                        expression(FDR~+~log[2]~FC)),
                       drawConnectors = TRUE,
                       widthConnectors = 0.5,
                       colConnectors = "grey50",
                       max.overlaps = Inf) +
  theme_pub

ggsave("output/Fig2A_volcano_TB0M_vs_LTB.pdf", p2a, width = 7, height = 7)
cat("  Saved: output/Fig2A_volcano_TB0M_vs_LTB.pdf\n")


# ================================================================
# Fig 2B — Volcano: TB_6M vs TB_0M (Recovery/Follow-up)
# ================================================================
cat("Generating Fig 2B: Volcano — TB_6M vs TB_0M\n")

deg_recovery$Gene <- rownames(deg_recovery)
deg_recovery$Label <- ""
top_idx <- order(abs(deg_recovery$logFC) * (deg_recovery$adj.P.Val < 0.05),
                 decreasing = TRUE)[1:min(top_n, nrow(deg_recovery))]
top_idx <- top_idx[!is.na(top_idx)]
deg_recovery$Label[top_idx] <- deg_recovery$Gene[top_idx]

p2b <- EnhancedVolcano(deg_recovery,
                       lab = deg_recovery$Label,
                       x = "logFC",
                       y = "adj.P.Val",
                       title = "Recovery: TB_6M vs TB_0M",
                       subtitle = "",
                       caption = "",
                       pCutoff = 0.05,
                       FCcutoff = 1,
                       labSize = 4.5,
                       pointSize = 0.6,
                       colAlpha = 0.8,
                       legendPosition = "top",
                       legendLabels = c("NS", expression(log[2]~FC), "FDR",
                                        expression(FDR~+~log[2]~FC)),
                       drawConnectors = TRUE,
                       widthConnectors = 0.5,
                       colConnectors = "grey50",
                       max.overlaps = Inf) +
  theme_pub

ggsave("output/Fig2B_volcano_TB6M_vs_TB0M.pdf", p2b, width = 7, height = 7)
cat("  Saved: output/Fig2B_volcano_TB6M_vs_TB0M.pdf\n")


# ================================================================
# Fig 2C — Venn diagrams: Transcriptional Reversal
# ================================================================
cat("Generating Fig 2C: Venn diagrams\n")

# ── Helper: all expressed genes (background) ──
background_genes <- rownames(exprs_symbol)[keep_expr]

# ── Left: Restoration Core (Down in TB ∩ Up after Tx) ──
# List A: genes significantly DOWN in TB_0M vs LTB
down_in_TB <- rownames(deg_disease)[sig_disease & deg_disease$logFC < 0]
down_in_TB <- intersect(down_in_TB, background_genes)
# List B: genes significantly UP after treatment (TB_6M vs TB_0M)
up_after_Tx <- rownames(deg_recovery)[sig_recovery & deg_recovery$logFC > 0]
up_after_Tx <- intersect(up_after_Tx, background_genes)

venn_recovery <- list(
  "Down in TB"  = down_in_TB,
  "Up after Tx" = up_after_Tx
)

p2c_left <- ggVennDiagram(venn_recovery,
                          set_size = 4.5,
                          label_alpha = 0,
                          label_percent_digit = 1) +
  scale_fill_gradient(low = "#F0F8FF", high = "#2E86C1") +
  labs(title = "Restoration Core") +
  theme_pub +
  theme(legend.position = "none",
        plot.margin = margin(10, 10, 10, 10))

# ── Right: Disease Drivers (Up in TB ∩ Down after Tx) ──
# List C: genes significantly UP in TB_0M vs LTB
up_in_TB <- rownames(deg_disease)[sig_disease & deg_disease$logFC > 0]
up_in_TB <- intersect(up_in_TB, background_genes)
# List D: genes significantly DOWN after treatment (TB_6M vs TB_0M)
down_after_Tx <- rownames(deg_recovery)[sig_recovery & deg_recovery$logFC < 0]
down_after_Tx <- intersect(down_after_Tx, background_genes)

venn_disease <- list(
  "Up in TB"       = up_in_TB,
  "Down after Tx"  = down_after_Tx
)

p2c_right <- ggVennDiagram(venn_disease,
                           set_size = 4.5,
                           label_alpha = 0,
                           label_percent_digit = 1) +
  scale_fill_gradient(low = "#F0F8FF", high = "#E64B35") +
  labs(title = "Disease Drivers") +
  theme_pub +
  theme(legend.position = "none",
        plot.margin = margin(10, 10, 10, 10))

p2c <- p2c_left + p2c_right +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 16))

ggsave("output/Fig2C_venn_transcriptional_reversal.pdf", p2c, width = 10, height = 5)
cat("  Saved: output/Fig2C_venn_transcriptional_reversal.pdf\n")


# ================================================================
# Fig 2D — Mfuzz temporal clustering (Restoration Core)
# ================================================================
cat("Generating Fig 2D: Mfuzz temporal clustering\n")

# ---- Median expression per timepoint ----
timepoints <- c("LTB", "TB_0M", "TB_0.5M", "TB_2M", "TB_6M", "TB_12M")
expr_long <- lapply(timepoints, function(tp) {
  cols <- which(group == tp)
  if (length(cols) == 0) return(NULL)
  rowMedians(exprs_symbol[genes_core_recovery, cols, drop = FALSE])
})
names(expr_long) <- timepoints
expr_long <- Filter(Negate(is.null), expr_long)
expr_mat <- do.call(cbind, expr_long)
colnames(expr_mat) <- names(expr_long)
rownames(expr_mat) <- genes_core_recovery

# ---- Mfuzz soft clustering ----
# Standardize rows (genes) to mean=0, SD=1 (Mfuzz requirement)
expr_mf <- new("ExpressionSet", exprs = expr_mat)
expr_mf <- standardise(expr_mf)

# Determine optimal m (fuzzifier) using mestimate
m <- mestimate(expr_mf)
cat(sprintf("  Mfuzz fuzzifier m = %.2f\n", m))

# Cluster with c=3 (try 3-4 depending on gene set size)
n_clusters <- if (length(genes_core_recovery) > 50) 4 else 3
set.seed(42)
cl_mf <- mfuzz(expr_mf, centers = n_clusters, m = m)

# ---- Publication-quality cluster plot with ggplot2 ----
cluster_membership <- cl_mf$cluster
expr_plot <- as.data.frame(expr_mat)
expr_plot$Gene <- rownames(expr_plot)
expr_plot$Cluster <- factor(cluster_membership[rownames(expr_plot)])

expr_plot_long <- expr_plot %>%
  pivot_longer(-c(Gene, Cluster), names_to = "Timepoint", values_to = "Expression") %>%
  mutate(Timepoint = factor(Timepoint, levels = timepoints))

# Summarise: mean ± SD per cluster per timepoint
cluster_summary <- expr_plot_long %>%
  group_by(Cluster, Timepoint) %>%
  summarise(
    Mean = mean(Expression, na.rm = TRUE),
    SD = sd(Expression, na.rm = TRUE),
    N = n(),
    .groups = "drop"
  )

cluster_colors <- brewer.pal(max(3, n_clusters), "Set1")[seq_len(n_clusters)]
names(cluster_colors) <- levels(cluster_summary$Cluster)

p2d <- ggplot(cluster_summary, aes(x = Timepoint, y = Mean,
                                    color = Cluster, group = Cluster)) +
  geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD, fill = Cluster),
              alpha = 0.15, color = NA) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  scale_color_manual(values = cluster_colors) +
  scale_fill_manual(values = cluster_colors) +
  labs(title = "Mfuzz Temporal Clusters (Restoration Core)",
       x = "Timepoint",
       y = "Standardised Expression",
       color = "Cluster", fill = "Cluster") +
  theme_pub +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("output/Fig2D_mfuzz_temporal_clusters.pdf", p2d, width = 8, height = 5.5)
cat("  Saved: output/Fig2D_mfuzz_temporal_clusters.pdf\n")


# ================================================================
# Fig 2E — GO enrichment (clusterProfiler compareCluster)
# ================================================================
cat("Generating Fig 2E: GO BP enrichment\n")

# Prepare gene list: Restoration Core vs Disease Drivers
go_gene_list <- list(
  "Restoration_Core" = genes_core_recovery,
  "Disease_Drivers"  = genes_core_disease
)

# Remove empty gene sets and apply minimum size filter
go_gene_list <- Filter(function(x) length(x) >= 5, go_gene_list)

if (length(go_gene_list) >= 2) {
  # Convert gene symbols to Entrez IDs
  gene_df <- bitr(
    unique(unlist(go_gene_list)),
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db,
    drop = TRUE
  )

  # Map back to lists
  go_gene_entrez <- lapply(go_gene_list, function(genes) {
    gene_df$ENTREZID[gene_df$SYMBOL %in% genes]
  })
  go_gene_entrez <- Filter(function(x) length(x) >= 5, go_gene_entrez)

  if (length(go_gene_entrez) >= 2) {
    # Compare cluster
    cgo <- compareCluster(
      go_gene_entrez,
      fun = "enrichGO",
      OrgDb = org.Hs.eg.db,
      ont = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.2,
      readable = TRUE
    )

    if (!is.null(cgo) && nrow(cgo@compareClusterResult) > 0) {
      p2e <- dotplot(cgo, showCategory = 10,
                     title = "GO Biological Process Enrichment") +
        theme_pub +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))

      ggsave("output/Fig2E_GO_enrichment.pdf", p2e, width = 10, height = 7)
      cat("  Saved: output/Fig2E_GO_enrichment.pdf\n")
    } else {
      cat("  WARNING: No GO terms enriched. Skipping Fig 2E.\n")
    }
  } else {
    cat("  WARNING: <2 gene groups with >=5 Entrez IDs. Skipping Fig 2E.\n")
  }
} else {
  cat("  WARNING: <2 gene groups with >=5 genes. Skipping Fig 2E.\n")
}


# ================================================================
# Fig 2F — VIPER TF activity heatmap
# ================================================================
cat("Generating Fig 2F: VIPER TF activity heatmap\n")

# ---- Run VIPER via decoupleR with DoRothEA regulons ----
# Load DoRothEA regulons (human, confidence A+B+C)
regulons <- dorothea_hs %>%
  filter(confidence %in% c("A", "B", "C"))

# Run VIPER on all samples
# Use the full expression matrix (genes in rows matching regulons)
tf_activities <- decouple(
  mat = exprs_symbol,
  network = regulons,
  .source = "tf",
  .target = "target",
  .mor = "mor",
  methods = "viper",
  consensus_scores = TRUE
)

# Reshape to matrix: TFs x samples
tf_mat <- tf_activities %>%
  filter(statistic == "viper") %>%
  select(source, condition, score) %>%
  pivot_wider(names_from = condition, values_from = score) %>%
  column_to_rownames("source") %>%
  as.matrix()

# ---- Module scores: mean expression of core gene sets ----
# Restoration Core module score (per sample)
module_core <- colMeans(exprs_symbol[intersect(genes_core_recovery, rownames(exprs_symbol)), , drop = FALSE], na.rm = TRUE)
# Disease Drivers module score (per sample)
module_disease <- colMeans(exprs_symbol[intersect(genes_core_disease, rownames(exprs_symbol)), , drop = FALSE], na.rm = TRUE)

# ---- Correlation: TF activity vs module scores ----
cor_mat <- matrix(NA, nrow = nrow(tf_mat), ncol = 2,
                  dimnames = list(rownames(tf_mat), c("Restoration_Core", "Disease_Drivers")))
for (tf in rownames(tf_mat)) {
  cor_mat[tf, "Restoration_Core"] <- cor(tf_mat[tf, ], module_core, method = "spearman", use = "complete.obs")
  cor_mat[tf, "Disease_Drivers"]  <- cor(tf_mat[tf, ], module_disease, method = "spearman", use = "complete.obs")
}

# Filter TFs with at least one |rho| > 0.3
tf_select <- rownames(cor_mat)[rowSums(abs(cor_mat) > 0.3, na.rm = TRUE) > 0]

# For heatmap: Z-score normalize TF activity per row (TF)
tf_mat_plot <- tf_mat[tf_select, , drop = FALSE]
tf_mat_plot <- tf_mat_plot - rowMeans(tf_mat_plot, na.rm = TRUE)
if (ncol(tf_mat_plot) > 1) {
  tf_mat_plot <- tf_mat_plot / matrix(
    apply(tf_mat_plot, 1, sd, na.rm = TRUE),
    nrow = nrow(tf_mat_plot), ncol = ncol(tf_mat_plot), byrow = FALSE
  )
}

# Annotation: group membership
annotation_col <- data.frame(
  Group = group,
  Restoration_Core = scale(module_core)[, 1],
  Disease_Drivers  = scale(module_disease)[, 1],
  row.names = colnames(exprs_symbol)
)
# Keep only columns that are in the heatmap
annotation_col <- annotation_col[colnames(tf_mat_plot), , drop = FALSE]

# Annotation colors
ann_colors <- list(
  Group = c(
    "LTB"    = "#F7ABE8",
    "TB_0M"  = "#FC4E07",
    "TB_0.5M" = "#E7B800",
    "TB_2M"  = "#AD07E3",
    "TB_6M"  = "#57A4FD",
    "TB_12M" = "#1B9E77"
  )
)

# Cap extreme values for display
tf_mat_plot[tf_mat_plot > 3] <- 3
tf_mat_plot[tf_mat_plot < -3] <- -3

# Add correlation annotation rows
annotation_row <- data.frame(
  Cor_Core    = cor_mat[tf_select, "Restoration_Core"],
  Cor_Disease = cor_mat[tf_select, "Disease_Drivers"],
  row.names = tf_select
)

# Broken-stick: top TFs by absolute correlation
n_show <- min(50, nrow(tf_mat_plot))
tf_keep <- names(sort(rowSums(abs(cor_mat[tf_select, ])), decreasing = TRUE))[seq_len(n_show)]
tf_mat_plot <- tf_mat_plot[tf_keep, ]
annotation_row <- annotation_row[tf_keep, ]

pdf("output/Fig2F_VIPER_TF_heatmap.pdf", width = 10, height = max(6, n_show * 0.25))
pheatmap(tf_mat_plot,
         annotation_col = annotation_col,
         annotation_row = annotation_row,
         annotation_colors = ann_colors,
         scale = "none",
         color = colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(100),
         breaks = seq(-3, 3, length.out = 101),
         border_color = NA,
         fontsize = 8,
         fontsize_row = 7,
         fontsize_col = 6,
         clustering_method = "ward.D2",
         treeheight_row = 20,
         treeheight_col = 20,
         main = "VIPER TF Activity (DoRothEA A+B+C)")
dev.off()
cat("  Saved: output/Fig2F_VIPER_TF_heatmap.pdf\n")


# ================================================================
# Fig 2G — Fisher TF enrichment bubble plot
# ================================================================
cat("Generating Fig 2G: Fisher TF enrichment bubble plot\n")

# Background: all expressed genes (with DoRothEA targets info)
bg_genes <- background_genes

# DoRothEA regulons: TF -> target genes (all confidence)
tf_targets <- regulons %>%
  filter(tf %in% rownames(tf_mat)) %>%
  split(.$tf) %>%
  lapply(function(df) unique(df$target))

# Keep only TFs with >= 5 targets in the background
tf_targets <- Filter(function(genes) {
  sum(genes %in% bg_genes) >= 5
}, tf_targets)

# Fisher's exact test for each TF against restoration core genes
fisher_results <- do.call(rbind, lapply(names(tf_targets), function(tf_name) {
  targets <- unique(tf_targets[[tf_name]])
  targets_in_bg <- targets[targets %in% bg_genes]
  targets_in_core <- sum(targets_in_bg %in% genes_core_recovery)

  # 2x2 contingency table
  a <- targets_in_core                          # in core & target of TF
  b <- length(targets_in_bg) - targets_in_core  # not in core, target of TF
  c <- length(genes_core_recovery) - targets_in_core  # in core, not target
  d <- length(bg_genes) - a - b - c             # neither

  mat <- matrix(c(a, b, c, d), nrow = 2,
                dimnames = list(c("Target", "NonTarget"),
                                c("InCore", "NotInCore")))
  ft <- fisher.test(mat, alternative = "greater")
  data.frame(
    TF = tf_name,
    Overlap = a,
    TotalTargets = length(targets_in_bg),
    FoldEnrichment = (a / (a + b)) / ((a + c) / length(bg_genes)),
    Pval = ft$p.value,
    stringsAsFactors = FALSE
  )
}))

fisher_results$log10P <- -log10(fisher_results$Pval)
fisher_results$Padj <- p.adjust(fisher_results$Pval, method = "BH")
fisher_results$Sig <- ifelse(fisher_results$Padj < 0.05, "FDR < 0.05", "NS")

# Top 20 TFs by -log10 P-value (with at least 1 overlapped gene)
top20 <- fisher_results %>%
  filter(Overlap >= 1) %>%
  arrange(Pval) %>%
  head(20) %>%
  mutate(TF = factor(TF, levels = rev(TF)))

p2g <- ggplot(top20, aes(x = FoldEnrichment, y = TF)) +
  geom_point(aes(size = Overlap, fill = log10P),
             shape = 21, color = "grey30", stroke = 0.5) +
  scale_size_continuous(range = c(3, 12),
                        name = "Overlap") +
  scale_fill_gradient(low = "#F5F5F5", high = "#B2182B",
                      name = expression(-log[10](P))) +
  labs(title = "TF Enrichment in Restoration Core (Fisher's Exact Test)",
       x = "Fold Enrichment", y = "") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50", linewidth = 0.5) +
  theme_pub +
  theme(panel.grid.major.x = element_line(color = "grey90", linewidth = 0.3))

ggsave("output/Fig2G_TF_fisher_enrichment.pdf", p2g, width = 8, height = 7)
cat("  Saved: output/Fig2G_TF_fisher_enrichment.pdf\n")


# ================================================================
# Summary
# ================================================================
cat("\n========================================\n")
cat("Figure 2 generation complete.\n")
cat("Output files:\n")
cat("  output/Fig2A_volcano_TB0M_vs_LTB.pdf\n")
cat("  output/Fig2B_volcano_TB6M_vs_TB0M.pdf\n")
cat("  output/Fig2C_venn_transcriptional_reversal.pdf\n")
cat("  output/Fig2D_mfuzz_temporal_clusters.pdf\n")
if (file.exists("output/Fig2E_GO_enrichment.pdf")) {
  cat("  output/Fig2E_GO_enrichment.pdf\n")
}
cat("  output/Fig2F_VIPER_TF_heatmap.pdf\n")
cat("  output/Fig2G_TF_fisher_enrichment.pdf\n")
cat("========================================\n")
