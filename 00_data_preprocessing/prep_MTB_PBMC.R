# ============================================================
# prep_MTB_PBMC.R — EGAS00001005376 MTB Timecourse scRNA
# ============================================================
# Input:  1M-scBloodNL master RDS (sc_master_unprocessed.rds)
# Output: mtb_final.rds (Seurat object: UT, 3hMTB, 24hMTB)
# ============================================================

library(Seurat)
library(ggplot2)
library(dplyr)

# ── 1. Load master object ──────────────────────────────────
sc_master <- readRDS("1M-scBloodNL/sc_master_unprocessed.rds")

# ── 2. Extract MTB timepoints ──────────────────────────────
sc_mtb <- subset(sc_master, subset = timepoint %in% c("UT", "3hMTB", "24hMTB"))
rm(sc_master); gc()

# ── 3. Join layers (Seurat v5) ─────────────────────────────
sc_mtb <- JoinLayers(sc_mtb)

# ── 4. Standard pipeline ──────────────────────────────────
sc_mtb <- NormalizeData(sc_mtb) %>%
  FindVariableFeatures(nfeatures = 2000) %>%
  ScaleData() %>%
  RunPCA(verbose = FALSE) %>%
  RunUMAP(dims = 1:15, verbose = FALSE)

sc_mtb$timepoint <- factor(sc_mtb$timepoint, levels = c("UT", "3hMTB", "24hMTB"))

# ── 5. Save ─────────────────────────────────────────────────
saveRDS(sc_mtb, "mtb_final.rds")
cat(sprintf("Done: mtb_final.rds saved.\n"))
cat(sprintf("  Cells: %d | Genes: %d | Donors: %d\n",
            ncol(sc_mtb), nrow(sc_mtb), length(unique(sc_mtb$assignment))))
