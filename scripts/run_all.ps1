$ErrorActionPreference = "Stop"
$ProjectRoot = Split-Path -Parent $PSScriptRoot
Set-Location -LiteralPath $ProjectRoot

Rscript scripts/01_audit_data.R
Rscript scripts/02_discovery_gse40553.R
Rscript scripts/02b_module_trajectory_models.R
Rscript scripts/02c_donor_correlation_tests.R
Rscript scripts/03_development_gse31348.R
Rscript scripts/04_outcome_gse89403.R
python scripts/05_fetch_metanet.py
Rscript scripts/24_quantiseq_scale_audit.R
Rscript scripts/12_composition_adjustment.R
Rscript scripts/15_unbiased_tf_screen.R
Rscript scripts/16_candidate_composition_robustness.R
Rscript scripts/13_metanet_candidate_validation.R
Rscript scripts/14_pseudobulk_mtb_pbmc.R
Rscript scripts/25_tbri_standardization_audit.R
Rscript scripts/26_cohort_inclusion_audit.R
Rscript scripts/30_gse89403_hg19_hg38_sensitivity.R
Rscript scripts/31_gse40553_south_africa_sensitivity.R
Rscript scripts/32_tbri_bcell_correlation_components.R
Rscript scripts/34_composition_framework_limma.R
Rscript scripts/35_cohort_characteristics_table.R
Rscript scripts/36_global_reversal_correlation.R
Rscript scripts/37_tbri_bcell_time_adjusted.R
Rscript scripts/38_gse31348_cibersort_lm22_sensitivity.R
Rscript scripts/06_figures.R
python scripts/39_export_integrated_figure_sources.py

$GseSource = "E:\Singlecell\单细胞测序分析\192483\3-Celltype"
$GseLink = Join-Path $env:TEMP "tf_tb_gse192483_3celltype"
if (-not (Test-Path -LiteralPath $GseLink)) {
    New-Item -ItemType Junction -Path $GseLink -Target $GseSource | Out-Null
}
$env:GSE192483_RDS = Join-Path $GseLink "sce_celltype.rds"
try {
    Rscript scripts/17_gse192483_candidate_localization.R
}
finally {
    if (Test-Path -LiteralPath $GseLink) {
        Remove-Item -LiteralPath $GseLink -Force -ErrorAction SilentlyContinue
    }
}

python scripts/18_gse326212_candidate_tf_pseudobulk.py
Rscript scripts/19_gse326212_candidate_tf_models.R
python scripts/27_gse326212_stat4_raw_detection.py
Rscript scripts/28_gse326212_stat4_raw_threshold_models.R
Rscript scripts/23_gse326212_stat4_robustness.R
Rscript scripts/20_gse326212_invivo_figure.R
Rscript scripts/07_export_tables.R
python scripts/40_finalize_supplement_workbook.py
python scripts/08_build_manuscript.py
python scripts/21_sync_markdown.py
python scripts/09_validate_outputs.py
python scripts/22_build_submission_package.py
python scripts/10_pre_submission_gate.py
powershell -ExecutionPolicy Bypass -File scripts/11_write_delivery_manifest.ps1
