# ============================================================
# prep_GSE218065.R — PBMC scRNA (Control vs Mtb)
# ============================================================
# Input:  GSE218065_RAW/ (10x Genomics: GSM6734046_SampleA, GSM6734047_SampleB)
# Output: GSE218065_processed.rds (Seurat object)
# ============================================================
# Note: Requires scDblFinder. Set correct path to RAW data dir.

library(Seurat)
library(Matrix)
library(dplyr)
library(ggplot2)
library(harmony)
library(SingleCellExperiment)
library(scDblFinder)

# ── 1. Sample processing function ──────────────────────────
process_sample <- function(file_prefix, project_name) {
  data <- ReadMtx(
    mtx = paste0(file_prefix, "_matrix.mtx.gz"),
    features = paste0(file_prefix, "_features.tsv.gz"),
    cells = paste0(file_prefix, "_barcodes.tsv.gz"),
    feature.column = 2
  )
  seu <- CreateSeuratObject(counts = data, project = project_name,
                             min.cells = 3, min.features = 200)
  seu$group <- project_name
  seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
  seu[["percent.HBB"]] <- PercentageFeatureSet(seu, pattern = "^HB[AB]")
  seu[["percent.PPBP"]] <- PercentageFeatureSet(seu, features = c("PPBP", "PF4"))

  seu <- subset(seu, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 &
                  percent.mt < 15 & percent.HBB < 1 & percent.PPBP < 1)

  sce <- as.SingleCellExperiment(seu)
  set.seed(123)
  sce <- scDblFinder(sce, verbose = FALSE)
  seu$scDblFinder.class <- sce$scDblFinder.class
  seu <- subset(seu, subset = scDblFinder.class == "singlet")
  return(seu)
}

# ── 2. Process each sample ─────────────────────────────────
# Update the path to your GSE218065_RAW directory
raw_dir <- "GSE218065_RAW"

seu_A <- process_sample(file.path(raw_dir, "GSM6734046_SampleA"), "Control")
seu_B <- process_sample(file.path(raw_dir, "GSM6734047_SampleB"), "Mtb")

# ── 3. Merge and process ───────────────────────────────────
seu <- merge(seu_A, seu_B)
seu <- JoinLayers(seu)

seu <- NormalizeData(seu) %>%
  FindVariableFeatures(nfeatures = 2000) %>%
  ScaleData() %>%
  RunPCA(verbose = FALSE) %>%
  RunHarmony(group.by.vars = "group") %>%
  RunUMAP(reduction = "harmony", dims = 1:15, verbose = FALSE) %>%
  FindNeighbors(reduction = "harmony", dims = 1:15) %>%
  FindClusters(resolution = 0.5)

# ── 4. Cell type annotation ────────────────────────────────
# Semi-automatic annotation using canonical markers
seu <- RunUMAP(seu, reduction = "harmony", dims = 1:15, seed.use = 42)

# Save raw processed before annotation
saveRDS(seu, "GSE218065_unannotated.rds")

cat("Done: GSE218065_unannotated.rds saved.\n")
cat(sprintf("  Cells: %d | Genes: %d\n", ncol(seu), nrow(seu)))
