# ============================================================
# annotate_GSE192483.R — Lung scRNA Cell Type Annotation
# ============================================================
# Input:  GSE192483_filtered.rds
# Output: GSE192483_annotated.rds (with celltype metadata column)
# ============================================================

library(Seurat)
library(ggplot2)
library(dplyr)

sce <- readRDS("GSE192483_filtered.rds")
sel.clust <- "RNA_snn_res.0.5"
Idents(sce) <- sel.clust

# ── 1. Marker-based annotation ─────────────────────────────
marker_sets <- list(
  Myeloid = c('AIF1', 'LYZ', 'CD68'),
  DC = c("CD74", "CD83", "THBD", "ITGAX", "CLEC4C", "CLEC9A"),
  Epithelial = c('EPCAM', 'KRT18', 'KRT8', 'CDH1', 'KRT19', 'KRT7', 'SFTPC', 'FOXJ1'),
  Endothelial = c('PECAM1', 'CD34'),
  Fibroblast = c('COL1A1', 'PDGFRA', 'COL1A2'),
  T = c('CD3D', 'CD3E', 'CD2', 'CD4', 'CD3G', 'CD8A', 'IL7R', 'CD8B'),
  B = c('MS4A1', 'CD19', 'CD79A', 'MZB1', 'CD79B'),
  Plasma = c("MZB1", "SDC1"),
  NK = c("KLRD1", "GNLY", "KLRF1", "AREG", "XCL2", "HSPA6"),
  Mast = c('CPA3', 'MS4A2', 'KIT', 'TPSAB1', 'TPSB2'),
  Proliferative = c('MKI67', 'CCND1', 'TOP2A', 'CDC20'),
  Neutrophil = c('CSF3R', 'CXCL8', 'FCGR3A')
)

# Get DotPlot expression data for scoring
marker_expr <- DotPlot(sce, features = unique(unlist(marker_sets)), assay = 'RNA')$data

cluster_ids <- sort(unique(marker_expr$id))
celltype <- data.frame(ClusterID = cluster_ids, celltype = NA)

for (cl in cluster_ids) {
  subset_data <- marker_expr %>% filter(id == cl)
  scores <- sapply(marker_sets, function(markers) {
    mean(subset_data$avg.exp.scaled[subset_data$features.plot %in% markers], na.rm = TRUE)
  })
  celltype[celltype$ClusterID == cl, "celltype"] <- names(which.max(scores))
}

# ── 2. Assign annotations ─────────────────────────────────
sce$celltype <- "Unassigned"
for (i in 1:nrow(celltype)) {
  sce@meta.data[sce@meta.data[[sel.clust]] == celltype$ClusterID[i], 'celltype'] <- celltype$celltype[i]
}

cat("Cell type distribution:\n")
print(table(sce$celltype))

# ── 3. Save ─────────────────────────────────────────────────
saveRDS(sce, "GSE192483_annotated.rds")
cat("Done: GSE192483_annotated.rds saved.\n")
