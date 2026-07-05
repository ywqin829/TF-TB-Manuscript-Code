# ============================================================
# prep_GSE192483.R — Lung scRNA (FDG-LOW vs FDG-HIGH)
# ============================================================
# Input:  GSE192483_RAW/ (10x Genomics, multiple samples)
# Output: GSE192483_processed.rds (Seurat object with cell types)
# ============================================================
# Note: Requires downloading and organizing CEL files into
# per-sample subdirectories first (see 192483.R part 1).

library(Seurat)
library(dplyr)
library(harmony)
library(stringr)
library(GEOquery)

# ── 1. Read all samples ────────────────────────────────────
# Assumes GSE192483_RAW/ contains subdirectories per sample
# with barcodes.tsv.gz, features.tsv.gz, matrix.mtx.gz
dir <- "GSE192483_RAW"
sample_dirs <- list.dirs(path = dir, full.names = TRUE, recursive = FALSE)

sceList <- lapply(sample_dirs, function(folder) {
  sample_id <- basename(folder)
  cat("Reading:", sample_id, "\n")
  mat <- Read10X(folder)
  CreateSeuratObject(counts = mat, project = sample_id,
                     min.cells = 5, min.features = 300)
})
names(sceList) <- basename(sample_dirs)

# ── 2. Merge ────────────────────────────────────────────────
sce.all <- merge(x = sceList[[1]], y = sceList[-1],
                 add.cell.ids = names(sceList))
sce.all <- JoinLayers(sce.all)

# ── 3. Load metadata from GEO ──────────────────────────────
gse <- getGEO("GSE192483", GSEMatrix = TRUE)[[1]]
sample_meta <- pData(gse)
meta_data <- sample_meta %>%
  select(geo_accession, title, characteristics_ch1) %>%
  mutate(
    Group = str_extract(title, "H|L"),
    Patient = str_extract(title, "^SP\\d+"),
    Tissue_Type = case_when(
      characteristics_ch1 == "tissue: lung tissue with 18F-FDG avidity" ~ "High_FDG",
      characteristics_ch1 == "tissue: lung tissue" ~ "Low_FDG",
      TRUE ~ "Unknown"
    )
  )

sce.all$Sample_geo_accession <- str_extract(sce.all$orig.ident, "^GSM\\d+")
merged_meta <- merge(sce.all@meta.data, meta_data,
                     by = "Sample_geo_accession", all.x = TRUE)
rownames(merged_meta) <- rownames(sce.all@meta.data)
sce.all@meta.data <- merged_meta

# ── 4. QC ──────────────────────────────────────────────────
sce.all <- PercentageFeatureSet(sce.all, features = grep("^MT-", rownames(sce.all), value = TRUE),
                                 col.name = "percent_mito")
sce.all <- PercentageFeatureSet(sce.all, features = grep("^RP[SL]", rownames(sce.all), value = TRUE),
                                 col.name = "percent_ribo")
sce.all <- PercentageFeatureSet(sce.all, features = grep("^HB[A-Z]", rownames(sce.all), value = TRUE),
                                 col.name = "percent_hb")

sce.all.filt <- subset(sce.all,
                       subset = nFeature_RNA > 200 & nCount_RNA > 500 &
                         percent_mito < 25 & percent_ribo > 3 & percent_hb < 1)

# ── 5. Normalization, PCA, Harmony, UMAP, Clustering ──────
sce.all.filt <- NormalizeData(sce.all.filt) %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA(verbose = FALSE) %>%
  RunHarmony(group.by.vars = "orig.ident") %>%
  RunUMAP(reduction = "harmony", dims = 1:10, verbose = FALSE) %>%
  RunTSNE(reduction = "harmony", dims = 1:10, perplexity = 30) %>%
  FindNeighbors(reduction = "harmony", dims = 1:10) %>%
  FindClusters(resolution = 0.5)

saveRDS(sce.all.filt, "GSE192483_filtered.rds")
cat("Done: GSE192483_filtered.rds saved.\n")
