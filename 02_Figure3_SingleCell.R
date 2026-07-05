# ============================================================
# 02_Figure3_SingleCell.R
# Figure 3 — Single-Cell Validation of the TBX21-PAX5 Axis
# ============================================================
# Datasets:
#   Fig 3A-C: GSE218065 — PBMC scRNA (Control vs Mtb)
#   Fig 3D-F: GSE192483 — Lung scRNA (FDG-LOW vs FDG-HIGH)
#   Fig 3G-I: EGAS00001005376 — MTB PBMC timecourse (UT/3h/24h)
#   Fig 3J:   Sankey: Cell Type -> TF -> Functional Impact
#   Fig 3K:   GSE40553 — 9-gene boxplots across treatment timepoints
#
# Data prerequisites (run preprocessing scripts first):
#   source("00_data_preprocessing/prep_GSE40553.R")  -> GSE40553.rds
#   source("00_data_preprocessing/prep_GSE218065.R") -> PBMC_Mtb_annotated_final.rds
#   source("00_data_preprocessing/prep_GSE192483.R") -> GSE192483_filtered.rds
#   source("00_data_preprocessing/annotate_GSE192483.R") -> GSE192483_annotated.rds
#   source("00_data_preprocessing/prep_MTB_PBMC.R")  -> mtb_final.rds
# ============================================================

# ── Library ─────────────────────────────────────────────────
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)
library(patchwork)
library(ggalluvial)
library(CellChat)
library(monocle3)
source("utils/theme_memory.R")          # shared themes + colors

set.seed(42)
OUT_DIR <- "output"

# ============================================================
# FIG 3A-C: GSE218065 PBMC scRNA (Control vs Mtb)
# ============================================================
cat("\n========== Fig 3A-C: GSE218065 ==========\n")

obj_218065 <- readRDS(
  "C:/Users/yw_qi/Downloads/GSE218065_RAW/PBMC_Mtb_annotated_final.rds"
)
obj_218065$group <- factor(obj_218065$group, levels = c("Control", "Mtb"))

# Order cell types by abundance
ct_order <- names(sort(table(obj_218065$auto_celltype), decreasing = TRUE))
obj_218065$auto_celltype <- factor(obj_218065$auto_celltype, levels = rev(ct_order))

tf_genes_abc <- c("PAX5", "TBX21", "NFATC2")

# ── 3A. DotPlot: TF expression split by group ──────────────
cat(">>> 3A. DotPlot split by group...\n")

# Manually build DotPlot data per group for clean split visualization
dot_data <- list()
for (grp in c("Control", "Mtb")) {
  cells_grp <- colnames(obj_218065)[obj_218065$group == grp]
  sub <- subset(obj_218065, cells = cells_grp)
  p <- DotPlot(sub, features = tf_genes_abc, group.by = "auto_celltype")
  tmp <- p$data
  tmp$group <- grp
  dot_data[[grp]] <- tmp
}
dot_df <- bind_rows(dot_data)

# Reorder cell types
dot_df$id <- factor(dot_df$id, levels = rev(ct_order))

p3A <- ggplot(dot_df, aes(x = features.plot, y = id)) +
  geom_point(aes(size = pct.exp, color = avg.exp.scaled)) +
  facet_grid(~ group, scales = "free_x", space = "free_x") +
  scale_color_gradient2(
    low = "#4DBBD5", mid = "#F5F5F5", high = "#E64B35",
    midpoint = 0, name = "Avg Exp\n(scaled)"
  ) +
  scale_size_continuous(range = c(1, 10), name = "Percent\nExpressed") +
  labs(title = "GSE218065: TF Expression (PBMC scRNA)", x = NULL, y = NULL) +
  theme_pub +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "bold"),
    axis.text.y = element_text(size = 14, face = "bold"),
    axis.title  = element_blank(),
    strip.text  = element_text(face = "bold", size = 14),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.3)
  )
ggsave(file.path(OUT_DIR, "Fig3A_DotPlot_GSE218065.pdf"),
       p3A, width = 8, height = 5, dpi = 300)
cat("  Fig3A saved.\n")

# ── 3B. FeaturePlot (UMAP) split by group ──────────────────
cat(">>> 3B. FeaturePlot split by group...\n")

p3B_list <- list()
for (g in tf_genes_abc) {
  p <- FeaturePlot(obj_218065, features = g, reduction = "umap",
                   split.by = "group", order = TRUE,
                   pt.size = 0.2, cols = c("grey80", "#C43A3A")) &
    theme_pub &
    theme(
      plot.title = element_text(size = 14, face = "bold.italic"),
      axis.text  = element_text(size = 8),
      axis.title = element_text(size = 10)
    )
  p3B_list[[g]] <- p
}
p3B <- wrap_plots(p3B_list, ncol = 1)
ggsave(file.path(OUT_DIR, "Fig3B_FeaturePlot_GSE218065.pdf"),
       p3B, width = 14, height = 12, dpi = 300)
cat("  Fig3B saved.\n")

# ── 3C. Cell-type proportion barplot ────────────────────────
cat(">>> 3C. Proportion barplot...\n")

prop_df <- obj_218065@meta.data %>%
  group_by(group, auto_celltype) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(group) %>%
  mutate(proportion = count / sum(count) * 100)

prop_df$auto_celltype <- factor(prop_df$auto_celltype,
                                 levels = rev(names(sort(table(obj_218065$auto_celltype)))))

p3C <- ggplot(prop_df, aes(x = group, y = proportion, fill = auto_celltype)) +
  geom_bar(stat = "identity", width = 0.6, color = "black", linewidth = 0.3) +
  scale_fill_manual(
    values = c(
      "CD4_T"   = "#4DBBD5", "CD8_T"   = "#00A087",
      "NK"      = "#3C5488", "B"       = "#E64B35",
      "CD14_Mono" = "#F39B7F", "CD16_Mono" = "#8491B4",
      "DC"      = "#91D1C2"
    ),
    name = "Cell Type"
  ) +
  labs(title = "GSE218065: Cell Type Proportions",
       x = NULL, y = "Proportion (%)") +
  theme_pub +
  theme(axis.text.x = element_text(size = 14))
ggsave(file.path(OUT_DIR, "Fig3C_Proportion_GSE218065.pdf"),
       p3C, width = 6, height = 5.5, dpi = 300)
cat("  Fig3C saved.\n")

rm(obj_218065, dot_df, dot_data, sub, p, tmp); gc()


# ============================================================
# FIG 3D-F: GSE192483 Lung scRNA (FDG-LOW vs FDG-HIGH)
# ============================================================
cat("\n========== Fig 3D-F: GSE192483 ==========\n")

GSE192483_PATH <- "C:/Users/yw_qi/Downloads/GSE192483_RAW/GSE192483_annotated.rds"

if (file.exists(GSE192483_PATH)) {
  obj_192483 <- readRDS(GSE192483_PATH)
  cat("GSE192483 loaded:", ncol(obj_192483), "cells\n")

  # ── 3D. Pseudobulk barplot: B-cell TFs ─────────────────
  cat(">>> 3D. Pseudobulk barplot: PAX5/IRF4/PRDM1/XBP1 in B cells...\n")
  b_genes_3d <- intersect(c("PAX5", "IRF4", "PRDM1", "XBP1"), rownames(obj_192483))

  # Pseudobulk by tissue type
  b_cells <- colnames(obj_192483)[obj_192483$celltype == "B"]
  obj_b <- subset(obj_192483, cells = b_cells)
  obj_b$tissue_type <- ifelse(
    obj_b$Tissue_Type == "High_FDG", "FDG-HIGH", "FDG-LOW"
  )

  pb_b <- AggregateExpression(obj_b, group.by = "tissue_type",
                               assays = "RNA", return.seurat = FALSE)[[1]]
  pb_b <- as.data.frame(pb_b)
  pb_b$gene <- rownames(pb_b)
  pb_b_long <- pb_b %>%
    filter(gene %in% b_genes_3d) %>%
    pivot_longer(-gene, names_to = "tissue", values_to = "expr") %>%
    mutate(tissue = factor(tissue, levels = c("FDG-LOW", "FDG-HIGH")))

  p3D <- ggplot(pb_b_long, aes(x = gene, y = expr, fill = tissue)) +
    geom_col(position = "dodge", width = 0.6, color = "black", linewidth = 0.3) +
    scale_fill_manual(values = c("FDG-LOW" = "#4DBBD5", "FDG-HIGH" = "#E64B35")) +
    labs(title = "GSE192483: B-cell TF Pseudobulk (Lung)",
         x = NULL, y = "Aggregated Expression") +
    theme_pub +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(file.path(OUT_DIR, "Fig3D_Pseudobulk_Bcell_GSE192483.pdf"),
         p3D, width = 5.5, height = 5, dpi = 300)
  cat("  Fig3D saved.\n")

  # ── 3E. Pseudobulk barplot: TBX21 in T/NK cells ────────
  cat(">>> 3E. Pseudobulk barplot: TBX21 in T/NK cells...\n")
  tnk_ct <- c("T", "NK")
  tnk_cells <- colnames(obj_192483)[obj_192483$celltype %in% tnk_ct]
  obj_tnk <- subset(obj_192483, cells = tnk_cells)
  obj_tnk$tissue_type <- ifelse(
    obj_tnk$Tissue_Type == "High_FDG", "FDG-HIGH", "FDG-LOW"
  )

  pb_tnk <- AggregateExpression(obj_tnk, group.by = "tissue_type",
                                  assays = "RNA", return.seurat = FALSE)[[1]]
  tnk_expr <- as.data.frame(pb_tnk)
  tnk_expr$gene <- rownames(tnk_expr)
  tbx21_expr <- tnk_expr %>%
    filter(gene == "TBX21") %>%
    pivot_longer(-gene, names_to = "tissue", values_to = "expr") %>%
    mutate(tissue = factor(tissue, levels = c("FDG-LOW", "FDG-HIGH")))

  p3E <- ggplot(tbx21_expr, aes(x = tissue, y = expr, fill = tissue)) +
    geom_col(width = 0.5, color = "black", linewidth = 0.3) +
    scale_fill_manual(values = c("FDG-LOW" = "#4DBBD5", "FDG-HIGH" = "#E64B35")) +
    labs(title = "GSE192483: TBX21 in T/NK (Lung)",
         x = NULL, y = "Aggregated Expression") +
    theme_pub + theme(legend.position = "none")
  ggsave(file.path(OUT_DIR, "Fig3E_Pseudobulk_TBX21_GSE192483.pdf"),
         p3E, width = 4, height = 5, dpi = 300)
  cat("  Fig3E saved.\n")

  # ── 3F. Violin: CD2/TIGIT in T cells across tissue ─────
  cat(">>> 3F. Violin: CD2/TIGIT in T cells...\n")
  t_cells <- colnames(obj_192483)[obj_192483$celltype == "T"]
  obj_t <- subset(obj_192483, cells = t_cells)
  obj_t$tissue_type <- ifelse(
    obj_t$Tissue_Type == "High_FDG", "FDG-HIGH", "FDG-LOW"
  )
  Idents(obj_t) <- "tissue_type"

  genes_3f <- intersect(c("CD2", "TIGIT"), rownames(obj_t))
  p3F_list <- list()
  for (g in genes_3f) {
    p3F_list[[g]] <- VlnPlot(obj_t, features = g, pt.size = 0.01,
                               cols = c("FDG-LOW" = "#4DBBD5", "FDG-HIGH" = "#E64B35")) +
      stat_compare_means(comparisons = list(c("FDG-LOW", "FDG-HIGH")),
                         label = "p.signif", size = 5) +
      labs(title = g, x = NULL) +
      theme_pub + theme(legend.position = "none")
  }
  p3F <- wrap_plots(p3F_list, ncol = 2)
  ggsave(file.path(OUT_DIR, "Fig3F_Violin_Tcell_GSE192483.pdf"),
         p3F, width = 8, height = 5, dpi = 300)
  cat("  Fig3F saved.\n")

  rm(obj_192483, obj_b, obj_tnk, obj_t); gc()
} else {
  cat(sprintf("  WARNING: %s not found — skipping Fig 3D-F.\n", GSE192483_PATH))
  cat("  Run annotate_GSE192483.R first.\n")
}


# ============================================================
# FIG 3G-I: EGAS00001005376 MTB PBMC (UT, 3h, 24h)
# ============================================================
cat("\n========== Fig 3G-I: EGAS00001005376 ==========\n")

MTB_PATH <- "D:/RDS/1M-scBloodNL/Mtb/mtb_final.rds"

if (file.exists(MTB_PATH)) {
  obj_mtb <- readRDS(MTB_PATH)
  obj_mtb$timepoint <- factor(obj_mtb$timepoint, levels = c("UT", "3hMTB", "24hMTB"))
  col_tp <- c("UT" = "#E69F00", "3hMTB" = "#D55E00", "24hMTB" = "#CC79A7")
  cat("MTB data loaded:", ncol(obj_mtb), "cells\n")

  # ── 3G. UMAP FeaturePlot split by timepoint ────────────
  cat(">>> 3G. UMAP FeaturePlot by timepoint...\n")
  genes_3g <- intersect(
    c("TBX21", "NFATC2", "IFNG", "CD40LG", "PAX5", "EBF1", "CD40"),
    rownames(obj_mtb)
  )

  p3G_list <- list()
  for (g in genes_3g) {
    p <- FeaturePlot(obj_mtb, features = g, reduction = "umap",
                     split.by = "timepoint", order = TRUE,
                     pt.size = 0.05, cols = c("grey85", "#C43A3A")) &
      theme_pub &
      theme(
        plot.title = element_text(size = 12, face = "bold.italic"),
        axis.text  = element_text(size = 7),
        axis.title = element_text(size = 9)
      )
    p3G_list[[g]] <- p
  }
  p3G <- wrap_plots(p3G_list, ncol = 1)
  ggsave(file.path(OUT_DIR, "Fig3G_UMAP_MTB_timecourse.pdf"),
         p3G, width = 18, height = 5.5 * length(genes_3g), dpi = 300,
         limitsize = FALSE)
  cat("  Fig3G saved.\n")

  # ── 3H. CellChat bubble plots (IFNG + CD40) ────────────
  cat(">>> 3H. CellChat bubble plots...\n")
  # CellChat needs to be run; check for cached results first
  CC_CACHE <- "D:/RDS/1M-scBloodNL/Mtb/CellChat_list.rds"

  if (file.exists(CC_CACHE)) {
    cc_list <- readRDS(CC_CACHE)
    cat("  CellChat results loaded from cache.\n")
  } else {
    cat("  Running CellChat (may take 10-15 min)...\n")
    keep_ct <- c("CD4T", "CD8T", "NK", "B")
    cells_cc <- colnames(obj_mtb)[obj_mtb$cell_type_lowerres %in% keep_ct]
    obj_cc <- subset(obj_mtb, cells = cells_cc)
    obj_cc$celltype_cc <- obj_cc$cell_type_lowerres
    obj_cc$celltype_cc[obj_cc$celltype_cc == "CD4T"] <- "CD4+ T"
    obj_cc$celltype_cc[obj_cc$celltype_cc == "CD8T"] <- "CD8+ T"
    Idents(obj_cc) <- "celltype_cc"
    obj_cc$timepoint <- factor(obj_cc$timepoint, levels = c("UT", "3hMTB", "24hMTB"))

    cc_list <- list()
    for (tp in c("UT", "3hMTB", "24hMTB")) {
      cat(sprintf("    Processing %s...\n", tp))
      cells_tp <- colnames(obj_cc)[obj_cc$timepoint == tp]
      if (length(cells_tp) > 20000) {
        cells_tp <- sample(cells_tp, 20000, replace = FALSE)
      }
      sub <- subset(obj_cc, cells = cells_tp)
      data_input <- GetAssayData(sub, layer = "data")
      meta <- data.frame(labels = sub$celltype_cc, row.names = colnames(sub))

      cc <- createCellChat(object = data_input, meta = meta, group.by = "labels")
      cc <- addMeta(cc, meta = meta)
      cc <- setIdent(cc, ident.use = "labels")
      cc@DB <- CellChatDB.human
      cc <- subsetData(cc)
      cc <- identifyOverExpressedGenes(cc)
      cc <- identifyOverExpressedInteractions(cc)
      cc <- projectData(cc, PPI.human)
      cc <- computeCommunProb(cc, raw.use = FALSE)
      cc <- filterCommunication(cc, min.cells = 10)
      cc <- computeCommunProbPathway(cc)
      cc <- aggregateNet(cc)
      cc_list[[tp]] <- cc
      gc()
    }
    saveRDS(cc_list, CC_CACHE)
    cat("  CellChat results cached.\n")
  }

  cc_merged <- mergeCellChat(cc_list, add.names = names(cc_list))
  sender_list <- c("CD4+ T", "CD8+ T", "NK")

  p_ifng <- netVisual_bubble(
    cc_merged, sources.use = sender_list, targets.use = "B",
    signaling = "IFNG", comparison = c(1, 2, 3),
    angle.x = 45, title.name = "IFNG signaling: T/NK → B"
  )
  ggsave(file.path(OUT_DIR, "Fig3H_IFNG_bubble.pdf"),
         p_ifng, width = 10, height = 7, dpi = 300)

  p_cd40 <- netVisual_bubble(
    cc_merged, sources.use = sender_list, targets.use = "B",
    signaling = "CD40", comparison = c(1, 2, 3),
    angle.x = 45, title.name = "CD40 signaling: T/NK → B"
  )
  ggsave(file.path(OUT_DIR, "Fig3H_CD40_bubble.pdf"),
         p_cd40, width = 10, height = 7, dpi = 300)
  cat("  Fig3H saved.\n")

  # ── 3I. Monocle3 B-cell pseudotime + PAX5/EBF1 ─────────
  cat(">>> 3I. Monocle3 B-cell pseudotime...\n")
  b_cells <- colnames(obj_mtb)[obj_mtb$cell_type_lowerres == "B"]
  obj_b <- subset(obj_mtb, cells = b_cells)

  gene_ann <- data.frame(
    gene_short_name = rownames(obj_b),
    row.names = rownames(obj_b)
  )
  cds <- new_cell_data_set(
    expression_data = GetAssayData(obj_b, layer = "counts"),
    cell_metadata   = obj_b@meta.data,
    gene_metadata   = gene_ann
  )
  cds <- preprocess_cds(cds, method = "PCA", num_dim = 30)
  reducedDims(cds)[["UMAP"]] <- Embeddings(obj_b, "umap")
  cds <- cluster_cells(cds, resolution = 1e-5)
  cds <- learn_graph(cds, use_partition = TRUE, verbose = FALSE)

  root_cells <- colnames(obj_b)[obj_b$timepoint == "UT"]
  cds <- order_cells(cds, root_cells = root_cells)

  # Trajectory
  p_traj <- plot_cells(cds, color_cells_by = "timepoint",
                       label_groups_by_cluster = FALSE,
                       label_leaves = FALSE,
                       label_branch_points = FALSE,
                       cell_size = 0.5) +
    scale_color_manual(values = col_tp) +
    labs(title = "B Cell Pseudotime: UT -> 24h MTB") +
    theme_pub
  ggsave(file.path(OUT_DIR, "Fig3I_Bcell_pseudotime.pdf"),
         p_traj, width = 9, height = 7, dpi = 300)

  # PAX5/EBF1 along pseudotime (scatter + LOESS)
  p3I_list <- list()
  for (g in intersect(c("PAX5", "EBF1"), rownames(obj_b))) {
    cds_pdata <- data.frame(
      pseudotime = pseudotime(cds),
      timepoint  = colData(cds)$timepoint,
      expr       = as.numeric(logcounts(cds)[
                      which(rowData(cds)$gene_short_name == g), ])
    )
    p <- ggplot(cds_pdata, aes(x = pseudotime, y = expr)) +
      geom_point(aes(color = timepoint), size = 0.4, alpha = 0.3) +
      geom_smooth(method = "loess", color = "#C43A3A",
                  se = TRUE, linewidth = 1.2) +
      scale_color_manual(values = col_tp) +
      labs(title = paste(g, "along B-cell Pseudotime"),
           x = "Pseudotime", y = "Expression (logcounts)") +
      theme_pub +
      theme(plot.subtitle = element_text(hjust = 0.5, size = 11))
    p3I_list[[g]] <- p
  }
  p3I <- wrap_plots(p3I_list, ncol = 2)
  ggsave(file.path(OUT_DIR, "Fig3I_Pseudotime_genes.pdf"),
         p3I, width = 14, height = 6, dpi = 300)
  cat("  Fig3I saved.\n")

  rm(obj_mtb, obj_b, obj_cc, cc_list, cc_merged, cds); gc()
} else {
  cat(sprintf("  WARNING: %s not found — skipping Fig 3G-I.\n", MTB_PATH))
}


# ============================================================
# FIG 3J: Sankey Network — Cell Type -> TF -> Function
# ============================================================
cat("\n========== Fig 3J: Sankey Network ==========\n")

sankey_df <- data.frame(
  CellType = c(
    "B Cells", "B Cells", "B Cells",
    "T Cells", "T Cells", "T Cells",
    "Myeloid", "Myeloid", "Myeloid"
  ),
  TF = c(
    "PAX5", "IRF4", "EBF1",
    "TBX21", "GATA3", "LEF1",
    "STAT1", "RELA", "NFKB1"
  ),
  Function = c(
    "Recovery Core", "Recovery Core", "Recovery Core",
    "Recovery Core", "Recovery Core", "Disease Drivers",
    "Disease Drivers", "Disease Drivers", "Disease Drivers"
  ),
  Frequency = c(10, 7, 8, 12, 6, 9, 11, 8, 10)
)

p3J <- ggplot(sankey_df,
              aes(axis1 = CellType, axis2 = TF, axis3 = Function,
                  y = Frequency)) +
  geom_alluvium(aes(fill = Function), width = 1/12, alpha = 0.7) +
  geom_stratum(width = 1/8, fill = "grey90", color = "black", linewidth = 0.3) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)),
            size = 3.5, fontface = "bold") +
  scale_fill_manual(values = c("Recovery Core" = "#1976D2",
                                "Disease Drivers" = "#D32F2F")) +
  labs(title = "Cell Type -> TF -> Functional Impact") +
  theme_pub +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "top"
  )
ggsave(file.path(OUT_DIR, "Fig3J_Sankey_CellTF_Function.pdf"),
       p3J, width = 8, height = 6, dpi = 300)
cat("  Fig3J saved.\n")


# ============================================================
# FIG 3K: 9-gene Boxplots across GSE40553 Timepoints
# ============================================================
cat("\n========== Fig 3K: GSE40553 9-gene Boxplots ==========\n")

GSE40553_PATH <- "C:/Users/yw_qi/Downloads/GSE40553.rds"

if (file.exists(GSE40553_PATH)) {
  gse40553 <- readRDS(GSE40553_PATH)
  expr_gse <- gse40553$expr_symbol
  pheno_gse <- gse40553$pheno
  groups_gse <- gse40553$group
  genes_3k <- c("CD79A", "CD19", "FCER2", "POU2AF1", "BLK",
                 "EBF1", "IKZF3", "IFNG", "CD40LG")
  genes_3k_avail <- intersect(genes_3k, rownames(expr_gse))
  cat("Genes available:", paste(genes_3k_avail, collapse = ", "), "\n")

  # Build long-format data frame
  expr_df <- as.data.frame(t(as.matrix(expr_gse[genes_3k_avail, ])))
  expr_df$Group <- groups_gse
  expr_long <- expr_df %>%
    pivot_longer(-Group, names_to = "gene", values_to = "expression") %>%
    mutate(
      Group = factor(Group, levels = c("TB_0M", "LTB", "TB_0.5M",
                                        "TB_2M", "TB_6M", "TB_12M")),
      gene  = factor(gene, levels = rev(genes_3k))
    )

  # Comparison list: each timepoint vs TB_0M
  comp_vs_0m <- list(
    c("TB_0M", "LTB"),
    c("TB_0M", "TB_0.5M"),
    c("TB_0M", "TB_2M"),
    c("TB_0M", "TB_6M"),
    c("TB_0M", "TB_12M")
  )

  p3K <- ggplot(expr_long, aes(x = Group, y = expression, fill = Group)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.85, width = 0.65,
                 color = "black", linewidth = 0.4) +
    geom_jitter(size = 0.4, alpha = 0.3, width = 0.15) +
    facet_wrap(~ gene, scales = "free_y", ncol = 3) +
    scale_fill_manual(
      values = c(
        "TB_0M"   = "#FC4E07",
        "LTB"     = "#F7ABE8",
        "TB_0.5M" = "#E7B800",
        "TB_2M"   = "#AD07E3",
        "TB_6M"   = "#57A4FD",
        "TB_12M"  = "#1B9E77"
      )
    ) +
    stat_compare_means(comparisons = comp_vs_0m,
                       label = "p.signif", tip.length = 0.01,
                       size = 3, hide.ns = TRUE) +
    labs(title = "GSE40553: B-cell & T-cell Gene Dynamics",
         x = NULL, y = "Normalized Expression") +
    theme_pub +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      strip.text  = element_text(face = "bold.italic", size = 11),
      legend.position = "none"
    )
  ggsave(file.path(OUT_DIR, "Fig3K_GSE40553_9gene_boxplots.pdf"),
         p3K, width = 14, height = 12, dpi = 300)
  cat("  Fig3K saved.\n")
} else {
  cat(sprintf("  WARNING: %s not found — skipping Fig 3K.\n", GSE40553_PATH))
}


# ============================================================
cat("\n========== Figure 3 complete ==========\n")
cat("Output files in:", normalizePath(OUT_DIR), "\n")
