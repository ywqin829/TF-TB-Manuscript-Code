#!/usr/bin/env python3
"""Export independently runnable R source for each integrated project figure.

The existing scripts/06_figures.R remains the canonical layout source. This
script slices its complete plotting blocks into self-contained, reviewable R
files. Run it again after changing that source to keep the exports in sync.
"""
from pathlib import Path

root = Path(__file__).resolve().parents[1]
source = (root / "scripts/06_figures.R").read_text(encoding="utf-8")
out = root / "scripts/figure_sources"
out.mkdir(exist_ok=True)

def section(start, end):
    assert source.count(start) == 1, start
    assert source.count(end) == 1, end
    return source[source.index(start):source.index(end)]

prelude = section('source(file.path("scripts", "00_config.R"))', 'p2a <- volcano_plot')
render_helper = section('render_pdf_panel <- function', 'lung_tf <- read_tsv')
blocks = {
    "Figure1_StudyDesign": section('figure1_boxes <- tribble(', 'render_pdf_panel <- function'),
    "Figure2_Discovery": section('p2a <- volcano_plot', 'figure1_boxes <- tribble('),
    "Figure3_SingleCellContext": section('lung_tf <- read_tsv', 'trajectory_panel <- function'),
    "Figure4_OutcomeModels": section('trajectory_panel <- function', 'score_auc <- read_tsv'),
    "Figure5_TBRI_OutcomeAssociation": section('score_auc <- read_tsv', '# Figure 6:'),
    "Figure6_CompositionConfounding": section('# Figure 6:', '# Figure 7:'),
    "FigureS2_CellTypeSpecificMtbResponse": section('# Figure 7:', '# Supplementary figure: METANet'),
    "FigureS1_METANet": section('# Supplementary figure: METANet', '# Supplementary figure: descriptive'),
}

for name, block in blocks.items():
    helper = render_helper if name == "Figure3_SingleCellContext" else ""
    header = ("# Run from the project root: Rscript scripts/figure_sources/" + name + ".R\n"
              "# Generated from scripts/06_figures.R by scripts/39_export_integrated_figure_sources.py.\n"
              "# Uses existing analysis tables and source-context panels; no upstream analysis is rerun.\n\n")
    (out / f"{name}.R").write_text(header + prelude + "\n" + helper + "\n" + block, encoding="utf-8")

# The main Figure 7 already has its own fully independent R implementation.
(out / "Figure7_InVivo_Candidate_TF_Context.R").write_text(
    '# Run from the project root: Rscript scripts/figure_sources/Figure7_InVivo_Candidate_TF_Context.R\n'
    '# Canonical implementation: scripts/20_gse326212_invivo_figure.R.\n'
    'source(file.path("scripts", "20_gse326212_invivo_figure.R"))\n', encoding="utf-8")

for name, script in [("FigureS3_TBRI_Bcell_TimeAdjusted", "37_tbri_bcell_time_adjusted.R"),
                     ("FigureS4_GSE31348_CIBERSORT_LM22_Bcell_Trajectory", "40_editorial_core_validation.R")]:
    (out / f"{name}.R").write_text(
        f'# Run from the project root: Rscript scripts/figure_sources/{name}.R\n'
        f'# Canonical analysis and figure implementation: scripts/{script}.\n'
        '# This reruns the sensitivity model from retained inputs before drawing the figure.\n'
        f'source(file.path("scripts", "{script}"))\n', encoding="utf-8")

(out / "README.md").write_text(
    "# Integrated figure R sources\n\n"
    "Run each `.R` file from the project root. Files for Figures 1-6 and S1-S2 "
    "are generated from `scripts/06_figures.R` by "
    "`scripts/39_export_integrated_figure_sources.py`; rerun the exporter after "
    "editing the canonical plot script. Figure 7 and S3-S4 call their separate "
    "canonical R analysis scripts. The source-context panels in Figure 3 "
    "are composed in R from existing source PDF panels. No raw-data "
    "preprocessing is rerun for Figures 1-7 and S1-S2.\n\n"
    "From the project root, run a figure with "
    "`Rscript scripts/figure_sources/Figure1_StudyDesign.R` (substitute the "
    "desired script). These R entry points regenerate the files in `figures/`; "
    "Figures S3 and S4 rerun their retained-input sensitivity models. "
    "Figure 7 requires the prepared "
    "GSE326212 pseudobulk inputs. The source-context panel in Figure 3 "
    "requires its existing source PDFs.\n\n"
    "The full workflow is `powershell -ExecutionPolicy Bypass -File "
    "scripts/run_all.ps1`. It regenerates the figure sources, runs the "
    "sensitivity analyses, exports supplementary tables, and finalizes the "
    "XLSX workbook. Inspect local input requirements in the project README.\n",
    encoding="utf-8")
print(f"Exported {len(list(out.glob('*.R')))} R scripts to {out}")

# Final editorial numbering; preserve the canonical legacy plotting code.
import shutil
numbering = {'Figure3_SingleCellContext': 'Figure4_SingleCellContext', 'Figure4_OutcomeModels': 'Figure5_OutcomeModels', 'Figure5_TBRI_OutcomeAssociation': 'Figure6_TBRI_OutcomeAssociation', 'Figure6_CompositionConfounding': 'Figure3_CompositionConfounding', 'Figure7_InVivo_Candidate_TF_Context': 'FigureS5_InVivo_Candidate_TF_Context'}
for old, new in numbering.items():
    source_file = out / (old + '.R')
    text = source_file.read_text(encoding='utf-8').replace(old, new)
    if old.startswith('Figure7_'):
        text += '\nfile.copy(out_path("figures", "' + old + '.pdf"), out_path("figures", "' + new + '.pdf"), overwrite = TRUE)\n'
    (out / (new + '.R')).write_text(text, encoding='utf-8')
    legacy = out / 'archive' / 'pre_structural_20261031'
    legacy.mkdir(parents=True, exist_ok=True)
    shutil.move(str(source_file), str(legacy / source_file.name))
(out / 'README.md').write_text(
    '# Final integrated figure R sources\n\n'
    'Run from the project root. Final main figures are 1-6; supplements are S1-S5. '
    'The R entry files reproduce the corresponding PDF names. Main 1-6 and S1-S2 '
    'reuse prepared results. S3-S4 rerun retained-input models; S5 uses prepared BAL inputs. '
    'Legacy source filenames under archive are not final figure entries.\n\n'
    'Numbering: old 6 -> 3; old 3 -> 4; old 4 -> 5; old 5 -> 6; old 7 -> S5. '
    'Run scripts/39_export_integrated_figure_sources.py after changing scripts/06_figures.R.\n', encoding='utf-8')
