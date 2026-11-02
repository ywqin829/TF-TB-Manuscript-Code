source(file.path("scripts", "00_config.R"))

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
})

assert_file <- function(path) {
  if (!file.exists(path)) stop("Missing input: ", path, call. = FALSE)
}

purrr::walk(unlist(paths[c("gse40553", "gse31348", "gse89403")]), assert_file)

gse40553 <- readRDS(paths$gse40553)
gse31348 <- readRDS(paths$gse31348)
gse89403 <- readRDS(paths$gse89403)

parse_gse40553_subject <- function(title) {
  bracketed <- grepl("\\[[0-9]+/", title)
  out <- title
  out[bracketed] <- sub(".*\\[([0-9]+)/.*", "\\1", title[bracketed])
  out[!bracketed] <- sub("-.*$", "", title[!bracketed])
  trimws(out)
}

gse40553$pheno$SubjectID <- parse_gse40553_subject(gse40553$pheno$title)
gse40553$pheno$Group <- factor(
  gse40553$pheno$Group,
  levels = c("LTB", "TB_0M", "TB_0.5M", "TB_2M", "TB_6M", "TB_12M")
)

stopifnot(
  identical(colnames(gse40553$expr_symbol), rownames(gse40553$pheno)),
  identical(colnames(gse31348$expression), rownames(gse31348$metadata)),
  identical(colnames(gse89403$expression), rownames(gse89403$metadata))
)

count_longitudinal <- function(meta, subject, group, time) {
  meta %>%
    transmute(SubjectID = .data[[subject]],
              Group = .data[[group]],
              TimePoint = .data[[time]]) %>%
    group_by(Group, TimePoint) %>%
    summarise(
      samples = n(),
      subjects = n_distinct(SubjectID),
      .groups = "drop"
    )
}

counts_40553 <- count_longitudinal(gse40553$pheno, "SubjectID", "Group", "Group") %>%
  mutate(Dataset = "GSE40553", TimePoint = as.character(Group)) %>%
  select(Dataset, Group, TimePoint, samples, subjects)

counts_31348 <- count_longitudinal(
  gse31348$metadata, "SubjectID", "TimePoint", "TimePoint"
) %>%
  mutate(Dataset = "GSE31348", Group = "Pulmonary TB") %>%
  select(Dataset, Group, TimePoint, samples, subjects)

counts_89403 <- count_longitudinal(
  gse89403$metadata, "SubjectID", "Group", "TimePoint"
) %>%
  mutate(Dataset = "GSE89403") %>%
  select(Dataset, Group, TimePoint, samples, subjects)

sample_counts <- bind_rows(counts_40553, counts_31348, counts_89403) %>%
  arrange(Dataset, Group, TimePoint)

dataset_manifest <- tibble::tribble(
  ~Dataset, ~Role, ~Species, ~Tissue, ~Samples, ~Subjects, ~BiologicalUnit,
  ~PermittedInference, ~Caveat,
  "GSE40553", "Discovery", "Homo sapiens", "Whole blood", 176L, 74L, "Patient",
  "Longitudinal transcriptional reversal and TF-regulon activity",
  "Repeated samples require donor-blocked inference.",
  "GSE31348", "Treatment-program review and model development", "Homo sapiens",
  "Whole blood", 135L, 27L, "Patient",
  "Treatment-direction review, retrospective TBRI development and one exploratory TF-screen condition",
  "Cannot serve as independent validation after feature selection.",
  "GSE89403", "Outcome association and exploratory TF screening", "Homo sapiens",
  "Whole blood", 312L, 82L, "Patient",
  "Association of a cohort-standardized score and TF regulons with treatment outcome",
  "Not independent validation for the TF screen; poor-outcome group is small (8 patients).",
  "GSE192483", "Contextual single-cell evidence", "Homo sapiens", "Tuberculous lung",
  11L, 6L, "Patient",
  "Descriptive cell-lineage localization and within-patient high/low lesion comparison",
  "No independent healthy-control cohort inside GSE192483.",
  "GSE326212", "In vivo disease-state single-cell context", "Homo sapiens", "BAL",
  77L, 64L, "Patient/sample",
  "Cross-sectional donor-level pseudobulk comparison of active TB and IGRA-defined contacts",
  "No treatment follow-up; project is partly confounded with disease group.",
  "GSE218065", "Descriptive stimulation example", "Homo sapiens", "PBMC",
  2L, 1L, "Donor/library",
  "Descriptive display of marker localization after Mtb stimulation",
  "One human control and one stimulated library; not a validation cohort.",
  "METANet", "Orthogonal regulatory prior", "Homo sapiens", "GTEx whole blood and lung",
  NA_integer_, NA_integer_, "Tissue network",
  "TF-target candidate prioritization and network concordance",
  "Context-informed prediction from GTEx and aggregated ChIP-seq labels; not disease-specific causality."
)

readr::write_tsv(sample_counts, out_path("results", "sample_counts_by_dataset.tsv"))
readr::write_tsv(dataset_manifest, out_path("results", "dataset_manifest.tsv"))

saveRDS(
  list(gse40553 = gse40553, gse31348 = gse31348, gse89403 = gse89403),
  out_path("results", "analysis_inputs.rds")
)

audit_lines <- c(
  paste0("GSE40553 subjects parsed from sample titles: ", n_distinct(gse40553$pheno$SubjectID)),
  paste0("GSE31348 subjects: ", n_distinct(gse31348$metadata$SubjectID)),
  paste0("GSE89403 subjects: ", n_distinct(gse89403$metadata$SubjectID)),
  "All expression matrices and metadata rows are aligned.",
  "GSE31348 is classified as model development, not independent validation.",
  "GSE218065 is classified as a one-donor descriptive example."
)
writeLines(audit_lines, out_path("logs", "data_audit.txt"))
