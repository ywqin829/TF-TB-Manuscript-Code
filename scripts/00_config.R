options(stringsAsFactors = FALSE)
options(timeout = 600)
set.seed(20261028)

if (.Platform$OS.type == "windows" && !isTRUE(l10n_info()$`UTF-8`)) {
  for (locale in c(
    "Chinese (Simplified)_China.65001",
    "English_United States.utf8",
    "en_US.UTF-8"
  )) {
    if (nzchar(suppressWarnings(Sys.setlocale("LC_CTYPE", locale)))) {
      break
    }
  }
}

project_root <- Sys.getenv("TF_TB_PROJECT_ROOT", unset = "")
if (!nzchar(project_root)) {
  project_root <- normalizePath(file.path(getwd()), winslash = "/", mustWork = FALSE)
  if (!nzchar(project_root)) project_root <- getwd()
}

dirs <- c("scripts", "results", "figures", "docs", "tables", "logs", "supplement")
for (d in dirs) dir.create(file.path(project_root, d), showWarnings = FALSE, recursive = TRUE)

paths <- list(
  gse40553 = Sys.getenv("TF_TB_GSE40553_RDS", unset = ""),
  gse31348 = Sys.getenv("TF_TB_GSE31348_RDS", unset = ""),
  gse89403 = Sys.getenv("TF_TB_GSE89403_RDS", unset = ""),
  original_manuscript = Sys.getenv("TF_TB_MANUSCRIPT_DOCX", unset = ""),
  original_code = Sys.getenv("TF_TB_LEGACY_CODE", unset = ""),
  legacy_figures = Sys.getenv("TF_TB_LEGACY_FIGURES", unset = "")
)

out_path <- function(...) file.path(project_root, ...)

palette_group <- c(
  "Successful_Treatment" = "#1976D2",
  "Poor_Outcome" = "#D32F2F",
  "LTB" = "#4D9221",
  "TB_0M" = "#C51B7D",
  "TB_0.5M" = "#F1A340",
  "TB_2M" = "#8073AC",
  "TB_6M" = "#3182BD",
  "TB_12M" = "#00A6A6"
)

theme_tf <- function(base_size = 12) {
  base_size <- base_size + 1
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      text = ggplot2::element_text(family = "sans", colour = "black"),
      axis.text = ggplot2::element_text(
        colour = "black", size = base_size, face = "bold"
      ),
      axis.title = ggplot2::element_text(
        face = "bold", colour = "black", size = base_size + 2
      ),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(
        face = "bold", colour = "black", size = base_size + 1
      ),
      legend.title = ggplot2::element_text(face = "bold", size = base_size + 1),
      legend.text = ggplot2::element_text(size = base_size),
      plot.title = ggplot2::element_text(
        face = "bold", hjust = 0.5, size = base_size + 4
      ),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, size = base_size + 1)
    )
}

save_plot_pair <- function(plot, stem, width, height) {
  ggplot2::ggsave(out_path("figures", paste0(stem, ".pdf")), plot,
                  width = width, height = height, device = grDevices::cairo_pdf,
                  bg = "white")
}

write_session <- function(name = "session_info.txt") {
  utils::capture.output(sessioninfo::session_info(),
                        file = out_path("logs", name))
}

message("Project root: ", project_root)

