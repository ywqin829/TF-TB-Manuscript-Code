# ============================================================
# Figure 4: Systemic Validation & Cellular Crosstalk
#  4A — CIBERSORT immune cell dynamics (GSE31348)
#  4B — PAX5/TBX21/NFATC2 longitudinal expression (GSE31348)
#  4C — NicheNet ligand-receptor (CD4+ T -> B)
#  4D — scRNA pseudobulk validation (GSE192483)
# ============================================================
# Data Dependencies:
#   GSE31348: Affymetrix HG-U133 Plus 2.0 (5 timepoints)
#   GSE192483: Lung scRNA-seq (HC / FDR-LOW / FDR-HIGH)
#   NicheNet DB: ligand_target_matrix.rds, lr_network.rds,
#                weighted_networks.rds (external, see README)
# ============================================================

source("utils/theme_memory.R")

library(tidyverse)
library(ggpubr)
library(cowplot)
library(patchwork)

OUT <- "output"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# ────────────────────────────────────────────────────────────
# Fig 4A — CIBERSORT Immune Cell Dynamics (GSE31348)
# ────────────────────────────────────────────────────────────
# NOTE: CIBERSORT must be run externally on the GSE31348
# expression matrix. The expected input is a CSV with columns:
#   SampleID, Timepoint, <22 cell types>
# Running CIBERSORT requires a licence from
#   https://cibersortx.stanford.edu/
#
# If you already have CIBERSORT results, place them at:
#   data/GSE31348_cibersort_results.csv
# or adjust the path below.

cibersort_file <- "data/GSE31348_cibersort_results.csv"

if (!file.exists(cibersort_file)) {
  message("[4A] CIBERSORT results not found at: ", cibersort_file)
  message("[4A] Skipping Fig 4A. See script header for instructions.")
  p_4a <- ggplot() +
    annotate("text", x = 0.5, y = 0.5, size = 5,
             label = "Fig 4A: Run CIBERSORT externally\nSee script for details") +
    theme_void()
} else {
  # ── Load & aggregate ────────────────────────────────────
  cib <- read.csv(cibersort_file, row.names = 1)

  cib_agg <- cib %>%
    mutate(
      Timepoint = factor(Timepoint,
                         levels = c("Week 0", "Week 1", "Week 2",
                                    "Week 4", "Week 26")),
      `Total T Cells`  = `T cells CD8` + `T cells CD4 naive` +
        `T cells CD4 memory resting` + `T cells CD4 memory activated` +
        `T cells follicular helper` + `T cells regulatory (Tregs)` +
        `T cells gamma delta`,
      `Monocytes`      = `Monocytes` + `Macrophages M0` +
        `Macrophages M1` + `Macrophages M2` +
        `Dendritic cells resting` + `Dendritic cells activated`,
      `Neutrophils`    = `Neutrophils`,
      `B Cells`        = `B cells naive` + `B cells memory`,
      `Plasma Cells`   = `Plasma cells`,
      `Eosinophils`    = `Eosinophils`
    ) %>%
    select(SampleID, Timepoint,
           `Total T Cells`, `Total B Cells` = `B Cells`,
           Monocytes, Neutrophils, `Plasma Cells`, Eosinophils) %>%
    pivot_longer(-c(SampleID, Timepoint),
                 names_to = "CellType", values_to = "Fraction")

  # ── Statistics: each timepoint vs Week 0 ────────────────
  stats_4a <- cib_agg %>%
    group_by(CellType) %>%
    rstatix::wilcox_test(Fraction ~ Timepoint,
                         comparisons = list(
                           c("Week 0", "Week 1"),
                           c("Week 0", "Week 2"),
                           c("Week 0", "Week 4"),
                           c("Week 0", "Week 26")
                         )) %>%
    rstatix::adjust_pvalue(method = "fdr") %>%
    rstatix::add_significance() %>%
    rstatix::add_xy_position(x = "Timepoint", group = "CellType")

  # ── Aggregate medians for trend line ────────────────────
  medians_4a <- cib_agg %>%
    group_by(CellType, Timepoint) %>%
    summarise(MedianFrac = median(Fraction, na.rm = TRUE), .groups = "drop")

  # ── Plot ────────────────────────────────────────────────
  p_4a <- ggplot(cib_agg, aes(x = Timepoint, y = Fraction)) +
    geom_boxplot(aes(fill = Timepoint), outlier.shape = NA,
                 alpha = 0.8, width = 0.6, color = "black", linewidth = 0.5) +
    geom_jitter(width = 0.15, size = 1.8, alpha = 0.7,
                shape = 21, fill = "grey90", color = "grey30", stroke = 0.4) +
    geom_line(data = medians_4a, aes(x = Timepoint, y = MedianFrac,
                                     group = 1),
              color = "#D32F2F", linewidth = 0.8, linetype = "dashed") +
    geom_point(data = medians_4a, aes(x = Timepoint, y = MedianFrac),
               color = "#D32F2F", size = 2.5) +
    stat_pvalue_manual(stats_4a, label = "p.signif",
                       tip.length = 0.01, hide.ns = TRUE, size = 4) +
    facet_wrap(~ CellType, scales = "free_y", ncol = 3) +
    scale_fill_manual(values = col_31348) +
    labs(title = "CIBERSORT Immune Cell Dynamics (GSE31348)",
         y = "Fraction", x = "") +
    theme_pub +
    theme(legend.position = "none",
          axis.text.x  = element_text(angle = 45, hjust = 1, size = 10),
          plot.margin = margin(10, 10, 10, 10))

  ggsave(file.path(OUT, "Fig4A_CIBERSORT_dynamics.pdf"),
         p_4a, width = 12, height = 8, device = "pdf")
  message("[4A] Saved: ", file.path(OUT, "Fig4A_CIBERSORT_dynamics.pdf"))
}

# ────────────────────────────────────────────────────────────
# Fig 4B — TF / B-cell Gene Longitudinal Expression (GSE31348)
# ────────────────────────────────────────────────────────────
# NOTE: GSE31348 normalised expression matrix expected at:
#   data/GSE31348_expression.csv
# with rows = probesets (collapsed to symbols), cols = samples.
# Sample metadata expected:
#   data/GSE31348_pheno.csv
# with columns: SampleID, Timepoint.

mat_file  <- "data/GSE31348_expression.csv"
pheno_file <- "data/GSE31348_pheno.csv"

genes_4b <- c("NFATC2", "TBX21", "IFNG", "CD40LG",
              "EBF1", "PAX5", "CD19", "CD79A")

if (!file.exists(mat_file) || !file.exists(pheno_file)) {
  message("[4B] Expression data not found at: ", mat_file, " / ", pheno_file)
  message("[4B] Skipping Fig 4B.")
  p_4b <- ggplot() +
    annotate("text", x = 0.5, y = 0.5, size = 5,
             label = "Fig 4B: GSE31348 expression data required\nSee script for details") +
    theme_void()
} else {
  expr <- read.csv(mat_file, row.names = 1)
  pheno <- read.csv(pheno_file) %>%
    mutate(Timepoint = factor(Timepoint,
                              levels = c("Week 0", "Week 1", "Week 2",
                                         "Week 4", "Week 26")))

  # Subset to target genes (those present in the matrix)
  genes_present <- intersect(genes_4b, rownames(expr))
  genes_absent  <- setdiff(genes_4b, rownames(expr))
  if (length(genes_absent) > 0) {
    message("[4B] Genes not found in expression matrix: ",
            paste(genes_absent, collapse = ", "))
  }

  expr_long <- expr[genes_present, , drop = FALSE] %>%
    as.data.frame() %>%
    rownames_to_column("Gene") %>%
    pivot_longer(-Gene, names_to = "SampleID", values_to = "Expression") %>%
    left_join(pheno, by = "SampleID") %>%
    mutate(Gene = factor(Gene, levels = genes_present))

  # ── Statistics: each timepoint vs Week 0 per gene ──────
  stats_4b <- expr_long %>%
    group_by(Gene) %>%
    rstatix::wilcox_test(Expression ~ Timepoint,
                         comparisons = comp_w0) %>%
    rstatix::adjust_pvalue(method = "fdr") %>%
    rstatix::add_significance() %>%
    rstatix::add_xy_position(x = "Timepoint", group = "Gene")

  # ── Plot ────────────────────────────────────────────────
  p_4b <- ggplot(expr_long, aes(x = Timepoint, y = Expression)) +
    geom_boxplot(aes(fill = Timepoint), outlier.shape = NA,
                 alpha = 0.8, width = 0.6, color = "black", linewidth = 0.5) +
    geom_jitter(width = 0.15, size = 1.8, alpha = 0.7,
                shape = 21, fill = "grey90", color = "grey30", stroke = 0.4) +
    stat_pvalue_manual(stats_4b, label = "p.signif",
                       tip.length = 0.01, hide.ns = TRUE, size = 3.5) +
    facet_wrap(~ Gene, scales = "free_y", ncol = 4) +
    scale_fill_manual(values = col_31348) +
    labs(title = "TF / B-cell Gene Expression (GSE31348)",
         y = "Normalised Expression", x = "") +
    theme_pub +
    theme(legend.position = "none",
          axis.text.x  = element_text(angle = 45, hjust = 1, size = 10),
          plot.margin = margin(10, 10, 10, 10))

  ggsave(file.path(OUT, "Fig4B_TF_Bcell_longitudinal.pdf"),
         p_4b, width = 14, height = 8, device = "pdf")
  message("[4B] Saved: ", file.path(OUT, "Fig4B_TF_Bcell_longitudinal.pdf"))
}

# ────────────────────────────────────────────────────────────
# Fig 4C — NicheNet Ligand-Receptor (CD4+ T -> B)
# ────────────────────────────────────────────────────────────
# NOTE: Requires external NicheNet database files downloaded
# from https://github.com/saeyslab/nichenetr
#   - nichenet_db/ligand_target_matrix.rds
#   - nichenet_db/lr_network.rds
#   - nichenet_db/weighted_networks.rds
#
# Sender: CD4+ T cells       Receiver: B cells
# Ligands of interest: CD40LG, IFNG (and co-expressed)

nichenet_dir  <- "nichenet_db"
lr_matrix     <- file.path(nichenet_dir, "ligand_target_matrix.rds")
lr_network    <- file.path(nichenet_dir, "lr_network.rds")
weighted_nets <- file.path(nichenet_dir, "weighted_networks.rds")

if (!file.exists(lr_matrix)) {
  message("[4C] NicheNet database not found in ", nichenet_dir, "/")
  message("[4C] Skipping Fig 4C. See README for download instructions.")
  p_4c <- ggplot() +
    annotate("text", x = 0.5, y = 0.5, size = 5,
             label = "Fig 4C: NicheNet DB required\nDownload from github.com/saeyslab/nichenetr") +
    theme_void()
} else {
  library(nichenetr)

  # ── Load NicheNet data ──────────────────────────────────
  ligand_target_matrix  <- readRDS(lr_matrix)
  lr_network            <- readRDS(lr_network)
  weighted_networks     <- readRDS(weighted_nets)

  # ── Define sender / receiver ────────────────────────────
  # Sender: genes expressed in CD4+ T cells
  # Receiver: genes expressed in B cells
  # (Replace with actual scRNA DE results when available)

  sender_genes <- c("CD40LG", "IFNG", "TBX21", "NFATC2",
                     "IL2", "IL21", "TNF", "CSF2")
  receiver_genes <- c("PAX5", "EBF1", "CD19", "CD79A", "CD79B",
                       "MS4A1", "BACH2", "BCL6", "IRF4", "PRDM1",
                       "XBP1", "JCHAIN", "MZB1", "SDC1")

  genes_present_lr <- intersect(sender_genes, rownames(ligand_target_matrix))
  genes_present_rec <- intersect(receiver_genes, colnames(ligand_target_matrix))

  # ── Predict ligand activities ───────────────────────────
  ligand_activities <- predict_ligand_activities(
    targets = genes_present_rec,
    ligands = genes_present_lr,
    ligand_target_matrix = ligand_target_matrix
  )

  # ── Filter to top ligands ───────────────────────────────
  top_ligands <- ligand_activities %>%
    slice_max(pearson, n = 10) %>%
    pull(test_ligand)

  # ── Heatmap: ligand -> target ───────────────────────────
  vis <- ligand_target_matrix[top_ligands, genes_present_rec, drop = FALSE]

  # Save heatmap directly via pheatmap (not a ggplot object)
  pheatmap::pheatmap(
    vis,
    color = colorRampPalette(c("#4DBBD5", "white", "#E64B35"))(100),
    main = "NicheNet: CD4+ T -> B Cell Ligand-Target",
    fontsize_row = 9,
    fontsize_col = 9,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    border_color = NA,
    filename = file.path(OUT, "Fig4C_NicheNet_CD4T_to_B.pdf"),
    width = 8,
    height = 6
  )
  message("[4C] Saved: ", file.path(OUT, "Fig4C_NicheNet_CD4T_to_B.pdf"))

  # Also create a ggplot placeholder for the composite
  p_4c <- ggplot() +
    annotate("text", x = 0.5, y = 0.5, size = 5,
             label = "Fig 4C: See saved PDF\nNicheNet heatmap (non-ggplot)") +
    theme_void()
}

# ────────────────────────────────────────────────────────────
# Fig 4D — scRNA Pseudobulk Validation (GSE192483, B/Plasma)
# ────────────────────────────────────────────────────────────
# NOTE: Preprocessed GSE192483 pseudobulk data expected at:
#   data/GSE192483_pseudobulk.csv
# Columns: SampleID, Celltype, Group (HC / FDR-LOW / FDR-HIGH),
#          PAX5, IRF4, PRDM1, XBP1, ...

pb_file <- "data/GSE192483_pseudobulk.csv"

if (!file.exists(pb_file)) {
  message("[4D] Pseudobulk data not found at: ", pb_file)
  message("[4D] Skipping Fig 4D.")
  p_4d <- ggplot() +
    annotate("text", x = 0.5, y = 0.5, size = 5,
             label = "Fig 4D: GSE192483 pseudobulk data required\nSee script for details") +
    theme_void()
} else {
  pb <- read.csv(pb_file) %>%
    mutate(
      Group    = factor(Group, levels = c("HC", "FDR-LOW", "FDR-HIGH")),
      Celltype = factor(Celltype)
    )

  genes_4d <- c("PAX5", "IRF4", "PRDM1", "XBP1")
  genes_4d <- intersect(genes_4d, colnames(pb))

  pb_long <- pb %>%
    select(SampleID, Celltype, Group, all_of(genes_4d)) %>%
    pivot_longer(all_of(genes_4d), names_to = "Gene", values_to = "Expression") %>%
    mutate(Gene = factor(Gene, levels = genes_4d))

  # ── Statistics ──────────────────────────────────────────
  stats_4d <- pb_long %>%
    group_by(Celltype, Gene) %>%
    rstatix::wilcox_test(Expression ~ Group,
                         comparisons = list(
                           c("HC", "FDR-LOW"),
                           c("HC", "FDR-HIGH"),
                           c("FDR-LOW", "FDR-HIGH")
                         )) %>%
    rstatix::adjust_pvalue(method = "fdr") %>%
    rstatix::add_significance() %>%
    rstatix::add_xy_position(x = "Group", group = "Celltype")

  # ── Plot ────────────────────────────────────────────────
  p_4d <- ggplot(pb_long, aes(x = Group, y = Expression)) +
    stat_summary(geom = "bar", aes(fill = Group),
                 fun = "mean", width = 0.65, color = "black",
                 linewidth = 0.5, alpha = 0.85) +
    geom_jitter(width = 0.12, size = 2, alpha = 0.7,
                shape = 21, fill = "grey85", color = "grey30", stroke = 0.4) +
    stat_compare_means(comparisons = list(
                         c("HC", "FDR-LOW"),
                         c("HC", "FDR-HIGH")
                       ),
                       label = "p.signif", hide.ns = TRUE,
                       tip.length = 0.01, size = 3.5) +
    facet_grid(Celltype ~ Gene, scales = "free_y") +
    scale_fill_manual(values = c("HC"        = "#4DBBD5",
                                  "FDR-LOW"   = "#F5A623",
                                  "FDR-HIGH"  = "#E64B35")) +
    labs(title = "scRNA Pseudobulk: B & Plasma Cells (GSE192483)",
         y = "Pseudobulk Expression", x = "") +
    theme_pub +
    theme(legend.position = "none",
          axis.text.x  = element_text(angle = 45, hjust = 1, size = 10),
          strip.text   = element_text(face = "bold", size = 11),
          plot.margin  = margin(10, 10, 10, 10))

  ggsave(file.path(OUT, "Fig4D_pseudobulk_GSE192483.pdf"),
         p_4d, width = 10, height = 7, device = "pdf")
  message("[4D] Saved: ", file.path(OUT, "Fig4D_pseudobulk_GSE192483.pdf"))
}

# ────────────────────────────────────────────────────────────
# Composite Figure 4
# ────────────────────────────────────────────────────────────
# Arrange all four panels into a single figure.

p_composite <- plot_grid(
  p_4a, p_4b, p_4c, p_4d,
  ncol = 2,
  labels = c("A", "B", "C", "D"),
  label_size = 16,
  label_fontface = "bold"
)

ggsave(file.path(OUT, "Figure4_composite.pdf"),
       p_composite, width = 20, height = 16, device = "pdf")
message("Figure 4 composite saved: ", file.path(OUT, "Figure4_composite.pdf"))

# ────────────────────────────────────────────────────────────
# Session Info
# ────────────────────────────────────────────────────────────
sink(file.path(OUT, "session_info_Fig4.txt"))
cat("=== Figure 4: Systemic Validation & Cellular Crosstalk ===\n\n")
cat("Date: ", format(Sys.Date(), "%Y-%m-%d"), "\n\n")
sessionInfo()
sink()

message("=== Figure 4 complete ===")
message("Output directory: ", normalizePath(OUT))
