# Run from the project root: Rscript scripts/figure_sources/FigureS5_InVivo_Candidate_TF_Context.R
# Canonical implementation: scripts/20_gse326212_invivo_figure.R.
source(file.path("scripts", "20_gse326212_invivo_figure.R"))

file.copy(out_path("figures", "Figure7_InVivo_Candidate_TF_Context.pdf"), out_path("figures", "FigureS5_InVivo_Candidate_TF_Context.pdf"), overwrite = TRUE)
