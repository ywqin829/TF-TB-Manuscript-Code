source(file.path("scripts", "00_config.R"))

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
})

timepoints <- read_tsv(
  out_path("results", "cohort_timepoint_counts.tsv"),
  show_col_types = FALSE
)

format_counts <- function(dataset, group_column, time_column) {
  x <- timepoints %>% filter(Dataset == dataset)
  if (dataset == "GSE40553") {
    order_levels <- c("LTB", "TB_0M", "TB_0.5M", "TB_2M", "TB_6M", "TB_12M")
    x <- x %>% arrange(match(AssignedGroup, order_levels))
    return(paste0(x$AssignedGroup, "=", x$Samples, collapse = "; "))
  }
  if (dataset == "GSE89403") {
    x <- x %>%
      separate(
        TimeGroup,
        into = c("Outcome", "TimePoint"),
        sep = " \\| ",
        remove = FALSE
      ) %>%
      mutate(
        Outcome = recode(
          Outcome,
          Successful_Treatment = "Success",
          Poor_Outcome = "Poor"
        ),
        TimePoint = recode(
          TimePoint,
          `0_Baseline` = "Baseline",
          `1_Week1` = "Week 1",
          `2_Week4` = "Week 4",
          `3_Week24` = "Week 24"
        )
      ) %>%
      arrange(match(TimePoint, c("Baseline", "Week 1", "Week 4", "Week 24")),
              desc(Outcome == "Success"))
    return(paste0(x$TimePoint, " ", x$Outcome, "=", x$Samples, collapse = "; "))
  }
  if (!is.null(group_column)) {
    x <- x %>% mutate(Label = .data[[group_column]])
  } else {
    x <- x %>% mutate(Label = .data[[time_column]])
  }
  paste0(x$Label, "=", x$Samples, collapse = "; ")
}

cohort_table <- tribble(
  ~Dataset, ~Role, ~Source, ~PMID, ~DOI, ~Platform, ~Country,
  ~Current_participants, ~Current_samples, ~Timepoint_sample_counts,
  ~Treatment, ~HIV_status, ~Drug_susceptibility,
  ~Outcome_definition, ~Reporting_note,

  "GSE40553",
  "Discovery and population sensitivity analysis",
  "Bloom et al. 2012, PLoS One",
  "23056259",
  "10.1371/journal.pone.0046191",
  "GPL10558 Illumina HumanHT-12 V4.0",
  "South Africa and United Kingdom",
  74L,
  176L,
  format_counts("GSE40553", "AssignedGroup", NULL),
  "Standard multidrug anti-tuberculosis therapy in treated TB participants",
  "Negative in all participants in the source study",
  "South African active cases were drug-sensitive, culture-confirmed pulmonary TB; UK cases were active TB and drug-susceptibility was not harmonized in GEO",
  "Longitudinal response; all treated participants were clinically cured through 12-month follow-up",
  "The current analysis retained 28 South African and 7 UK active-TB baseline samples plus 38 South African LTB samples; one South African and one UK baseline sample without an analyzable disease/time label were excluded",

  "GSE31348",
  "Treatment-program review, TBRI development, and one TF-screen condition",
  "Cliff et al. 2013, J Infect Dis",
  "22872737",
  "10.1093/infdis/jis499",
  "GPL570 Affymetrix HG-U133 Plus 2.0",
  "South Africa",
  27L,
  135L,
  format_counts("GSE31348", NULL, "TimePoint"),
  "2 months of HRZE followed by 4 months of HR",
  "Negative in all patients in the source study",
  "First-episode smear-positive pulmonary TB infected with isoniazid- and rifampicin-sensitive bacilli",
  "All patients were cured after 6 months and remained disease-free for 2 years in the source study",
  "The deposited GEO series contains 27 arrays at each of five time points and all 135 were retained",

  "GSE89403",
  "Outcome association and one exploratory TF-screen condition",
  "Thompson et al. 2017, Tuberculosis (Edinb)",
  "29050771",
  "10.1016/j.tube.2017.08.004",
  "GPL11154 Illumina HiSeq 2000",
  "South Africa",
  82L,
  312L,
  format_counts("GSE89403", NULL, "TimeGroup"),
  "Standard six-month first-line anti-tuberculosis treatment",
  "HIV-uninfected adults in the source cohort",
  "Newly diagnosed pulmonary TB treated with standard first-line therapy; complete drug-susceptibility data were not harmonized in the extracted metadata",
  "Definite cure required sustained sputum culture negativity through month 6; failure was a positive month-6 culture; probable cure, possible cure, and unevaluable outcomes were excluded",
  "Current analysis contains 74 participants with successful treatment and 8 not-cured participants after restriction and one sample per subject-time"
)

write_tsv(
  cohort_table,
  out_path("results", "cohort_characteristics.tsv")
)

markdown_lines <- c(
  "# Table 1. Cohort Characteristics",
  "",
  paste0("| ", paste(names(cohort_table), collapse = " | "), " |"),
  paste0("| ", paste(rep("---", ncol(cohort_table)), collapse = " | "), " |"),
  apply(cohort_table, 1, function(row) {
    paste0("| ", paste(gsub("\\|", "\\\\|", row), collapse = " | "), " |")
  }),
  "",
  "Sources were checked against the GEO series records, the Bloom 2012 full text,",
  "the Cliff 2013 publication record, and the Thompson 2017 full text.",
  "",
  "The Cliff 2013 HIV status was not recoverable from the accessible GEO metadata",
  "and is therefore not asserted as negative in this table."
)
writeLines(
  markdown_lines,
  out_path("docs", "Table1_Cohort_Characteristics.md")
)

writeLines(c(
  paste0("Datasets: ", nrow(cohort_table)),
  paste0("Total current samples across datasets: ", sum(cohort_table$Current_samples)),
  "HIV status is explicitly supported for GSE40553 and GSE89403.",
  "GSE31348 HIV status is marked unavailable rather than inferred."
), out_path("logs", "cohort_characteristics_summary.txt"))
