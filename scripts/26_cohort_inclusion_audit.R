source(file.path("scripts", "00_config.R"))

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
})

inputs <- readRDS(out_path("results", "analysis_inputs.rds"))
geo_dir <- out_path("data", "geo_official")
`%||%` <- function(x, y) if (is.null(x)) y else x

flow_rows <- list()
time_rows <- list()
batch_rows <- list()

read_series_metadata <- function(path) {
  con <- gzfile(path, "rt")
  on.exit(close(con), add = TRUE)
  lines <- readLines(con, warn = FALSE)
  table_start <- match("!series_matrix_table_begin", lines)
  if (is.na(table_start)) stop("series matrix table start not found: ", path)
  sample_lines <- grep("^!Sample_", lines[seq_len(table_start - 1)], value = TRUE)
  parsed <- lapply(sample_lines, function(line) {
    parts <- strsplit(sub("^!Sample_", "", line), "\t", fixed = TRUE)[[1]]
    list(
      key = parts[1],
      values = gsub('^"|"$', "", parts[-1])
    )
  })
  columns <- list()
  counters <- list()
  for (item in parsed) {
    key <- item$key
    counters[[key]] <- (counters[[key]] %||% 0L) + 1L
    column <- if (counters[[key]] == 1L) {
      key
    } else {
      paste0(key, ".", counters[[key]] - 1L)
    }
    columns[[column]] <- item$values
  }
  result <- as.data.frame(columns, check.names = FALSE, stringsAsFactors = FALSE)
  rownames(result) <- result$geo_accession
  result
}

## GSE40553 -------------------------------------------------------------
gse40553_path <- file.path(geo_dir, "GSE40553_series_matrix.txt.gz")
if (file.exists(gse40553_path)) {
  p <- read_series_metadata(gse40553_path)
  current_40553 <- inputs$gse40553$pheno
  illness <- p$characteristics_ch1
  status <- sub("^status:\\s*", "", p$characteristics_ch1.4)
  p$SubjectID <- ifelse(
    grepl("^[0-9]+-[0-9]+$", p$title),
    sub("^([0-9]+)-.*$", "\\1", p$title),
    sub(".*\\[([0-9]+)/.*$", "\\1", p$title)
  )
  group <- case_when(
    grepl("LTB", illness) & grepl("^untreated latent TB$", status) ~ "LTB",
    grepl("PTB", illness) & grepl("^active TB pre-treatment$", status) ~ "TB_0M",
    grepl("PTB", illness) & grepl("^active TB 2 weeks post treatment$", status) ~ "TB_0.5M",
    grepl("PTB", illness) & grepl("^active TB 2 months post treatment$", status) ~ "TB_2M",
    grepl("PTB", illness) & grepl("^active TB 6 months post treatment$", status) ~ "TB_6M",
    grepl("PTB", illness) & grepl("^active TB 12 months post treatment$", status) ~ "TB_12M",
    TRUE ~ NA_character_
  )
  p$AssignedGroup <- group
  included <- p[!is.na(p$AssignedGroup), , drop = FALSE]
  flow_rows[[length(flow_rows) + 1]] <- tibble(
    Dataset = "GSE40553",
    Step = c(
      "Official GEO samples",
      "Samples with assignable disease/time group",
      "Current analysis samples",
      "Current unique subjects"
    ),
    Samples = c(nrow(p), nrow(included), nrow(current_40553),
                n_distinct(current_40553$SubjectID)),
    Subjects = c(n_distinct(p$SubjectID),
                 n_distinct(included$SubjectID),
                 n_distinct(current_40553$SubjectID),
                 n_distinct(current_40553$SubjectID)),
    Detail = c(
      "All deposited GEO samples",
      "Excluded samples without a recognized LTB/PTB treatment-time label",
      paste0("Current accession overlap: ",
             length(intersect(current_40553$geo_accession, rownames(p)))),
      "Patient IDs parsed from GEO sample metadata; UK contributed seven baseline TB samples"
    )
  )
  time_rows[[length(time_rows) + 1]] <- included %>%
    count(AssignedGroup, name = "Samples") %>%
    mutate(Dataset = "GSE40553", .before = 1)
  batch_rows[[length(batch_rows) + 1]] <- included %>%
    count(
      AssignedGroup,
      Batch = sub(
        "^geographical location:\\s*",
        "",
        `characteristics_ch1.3`
      ),
      name = "Samples"
    ) %>%
    mutate(Dataset = "GSE40553", .before = 1)
}

## GSE31348 -------------------------------------------------------------
gse31348_path <- file.path(geo_dir, "GSE31348_series_matrix.txt.gz")
if (file.exists(gse31348_path)) {
  p <- read_series_metadata(gse31348_path)
  current_31348 <- inputs$gse31348$metadata
  p$PatientID <- sub("^patient_id:\\s*", "", p$characteristics_ch1.1)
  flow_rows[[length(flow_rows) + 1]] <- tibble(
    Dataset = "GSE31348",
    Step = c("Official GEO samples", "Current analysis samples",
             "Current unique subjects"),
    Samples = c(nrow(p), nrow(current_31348), nrow(current_31348)),
    Subjects = c(n_distinct(p$PatientID),
                 n_distinct(current_31348$SubjectID),
                 n_distinct(current_31348$SubjectID)),
    Detail = c(
      "All deposited longitudinal samples",
      paste0("Current accession overlap: ",
             length(intersect(current_31348$geo_accession, rownames(p)))),
      "All deposited participants were retained"
    )
  )
  time_rows[[length(time_rows) + 1]] <- current_31348 %>%
    count(TimePoint, name = "Samples") %>%
    mutate(Dataset = "GSE31348", .before = 1)
}

## GSE89403 -------------------------------------------------------------
gse89403_path <- file.path(geo_dir, "GSE89403_series_matrix.txt.gz")
if (file.exists(gse89403_path)) {
  p <- read_series_metadata(gse89403_path)
  current_89403 <- inputs$gse89403$metadata
  p$SubjectID <- sub("^subject:\\s*", "", p$characteristics_ch1.2)
  p$TimeRaw <- sub("^time:\\s*", "", p$characteristics_ch1.5)
  p$OutcomeRaw <- sub(
    "^treatmentresult:\\s*",
    "",
    p$characteristics_ch1.4
  )
  p$MatchID <- paste0("S", p$SubjectID, "_", p$TimeRaw)
  allowed_time <- p$TimeRaw %in% c("DX", "day_7", "week_4", "week_24")
  allowed_outcome <- p$OutcomeRaw %in% c("Definite Cure", "Not Cured")
  filtered <- p[allowed_time & allowed_outcome, , drop = FALSE]
  deduplicated <- filtered %>% distinct(MatchID, .keep_all = TRUE)
  current_ids <- current_89403$SampleID
  if (is.null(current_ids)) current_ids <- rownames(current_89403)
  sample_map <- p %>%
    mutate(
      AnalysisTime = allowed_time,
      EligibleOutcome = allowed_outcome,
      InFinalSubjectTimeSet = MatchID %in% current_ids,
      SelectedExpressionGSM = geo_accession %in% current_89403$geo_accession,
      InclusionReason = case_when(
        !AnalysisTime ~ "excluded: outside four analysis times",
        !EligibleOutcome ~ "excluded: outcome not definite cure/not cured",
        SelectedExpressionGSM ~ "GSM selected in processed expression metadata",
        InFinalSubjectTimeSet ~ "eligible paired GSM; expression processing unresolved",
        TRUE ~ "eligible GEO metadata; no corresponding final record"
      )
    ) %>%
    transmute(GSM=geo_accession, SubjectID, TimePoint=TimeRaw, Outcome=OutcomeRaw,
              SubjectTimeID=MatchID, AnalysisTime, EligibleOutcome, InFinalSubjectTimeSet,
              SelectedExpressionGSM,
              InclusionReason)
  write_tsv(sample_map, out_path("results", "gse89403_gsm_subject_time_map.tsv"))
  flow_rows[[length(flow_rows) + 1]] <- tibble(
    Dataset = "GSE89403",
    Step = c(
      "Official GEO samples",
      "After restriction to four analysis timepoints",
      "After restriction to definite cure or not cured",
      "Unique eligible subject-time IDs in GEO metadata",
      "Current analysis samples",
      "Current unique subjects"
    ),
    Samples = c(nrow(p), sum(allowed_time), sum(allowed_time & allowed_outcome),
                nrow(deduplicated), nrow(current_89403), nrow(current_89403)),
    Subjects = c(
      n_distinct(p$SubjectID),
      n_distinct(p$SubjectID[allowed_time]),
      n_distinct(p$SubjectID[allowed_time & allowed_outcome]),
      n_distinct(deduplicated$SubjectID),
      n_distinct(current_89403$SubjectID),
      n_distinct(current_89403$SubjectID)
    ),
    Detail = c(
      "All technical replicates, controls and cases; metadata contain 142 subject IDs",
      "DX, day_7, week_4 and week_24",
      "Definite Cure and Not Cured only",
      "Unique metadata IDs; expression-level technical-replicate processing not verified",
      paste0("Current accession overlap: ",
             length(intersect(current_89403$geo_accession, rownames(p)))),
      "Outcome-linked participants retained"
    )
  )
  time_rows[[length(time_rows) + 1]] <- current_89403 %>%
    count(TimePoint, Group, name = "Samples") %>%
    unite("TimeGroup", Group, TimePoint, sep = " | ") %>%
    mutate(Dataset = "GSE89403", .before = 1)
}

flow <- bind_rows(flow_rows)
time_counts <- bind_rows(time_rows)
batch_counts <- bind_rows(batch_rows)

write_tsv(flow, out_path("results", "cohort_inclusion_flow.tsv"))
write_tsv(time_counts, out_path("results", "cohort_timepoint_counts.tsv"))
write_tsv(batch_counts, out_path("results", "cohort_batch_group_counts.tsv"))

writeLines(c(
  paste0("GSE40553 official matrix present: ", file.exists(gse40553_path)),
  paste0("GSE31348 official matrix present: ", file.exists(gse31348_path)),
  paste0("GSE89403 official matrix present: ", file.exists(gse89403_path)),
  "GSE89403 analysis matrix has one subject-time expression column; source processing of technical replicates has not been verified.",
  "GSE40553 samples without a recognizable disease/treatment-time label are excluded by the documented preprocessing rule."
), out_path("logs", "cohort_inclusion_audit_summary.txt"))
