#!/usr/bin/env python
"""Build a comprehensive Markdown report of all current result domains."""

from pathlib import Path
import math

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
DOCS = ROOT / "docs"
SUPPLEMENT = ROOT / "supplement"
OUT = DOCS / "Complete_Results_Report.md"


def read_tsv(name: str) -> pd.DataFrame:
    return pd.read_csv(RESULTS / name, sep="\t")


def fmt(value) -> str:
    if pd.isna(value):
        return "NA"
    if isinstance(value, bool):
        return "yes" if value else "no"
    if isinstance(value, (int, float)):
        value = float(value)
        if value == 0:
            return "0"
        if abs(value) < 0.001 or abs(value) >= 10000:
            return f"{value:.3e}"
        return f"{value:.4g}"
    return str(value).replace("|", "\\|").replace("\n", " ")


def md_table(
    data: pd.DataFrame,
    columns: list[str],
    headers: list[str] | None = None,
    max_rows: int | None = None,
) -> str:
    view = data.loc[:, columns].copy()
    if max_rows is not None:
        view = view.head(max_rows)
    headers = headers or columns
    lines = [
        "| " + " | ".join(headers) + " |",
        "|" + "|".join(["---"] * len(columns)) + "|",
    ]
    for _, row in view.iterrows():
        lines.append("| " + " | ".join(fmt(row[col]) for col in columns) + " |")
    return "\n".join(lines)


def file_link(relative: str, label: str | None = None) -> str:
    path = ROOT / relative
    label = label or path.name
    return f"[{label}](</{path.as_posix()}>)"


def section_cohorts() -> list[str]:
    flow = read_tsv("cohort_inclusion_flow.tsv")
    time = read_tsv("cohort_timepoint_counts.tsv")
    batch = read_tsv("cohort_batch_group_counts.tsv")
    subject = read_tsv("gse326212_subject_unit_audit.tsv")
    characteristics = read_tsv("cohort_characteristics.tsv")
    sa_overlap = read_tsv("gse40553_sa_vs_original_module_overlap.tsv")
    sa_models = read_tsv("gse40553_sa_core_gse31348_models.tsv")
    sa_tf = read_tsv("gse40553_south_africa_candidate_tf_contrasts.tsv")
    sa_go = read_tsv("gse40553_south_africa_restoration_go.tsv")
    return [
        "## 1. 队列纳入、时间点与技术重复",
        "",
        "### 2.0 主文队列表",
        "",
        md_table(
            characteristics,
            [
                "Dataset", "Source", "Platform", "Country",
                "Current_participants", "Current_samples", "HIV_status",
                "Drug_susceptibility", "Outcome_definition",
            ],
            [
                "数据集", "来源", "平台", "国家", "当前患者数",
                "当前样本数", "HIV 状态", "药敏状态", "结局定义",
            ],
        ),
        "",
        "### 2.1 纳入流程",
        "",
        md_table(
            flow,
            ["Dataset", "Step", "Samples", "Subjects", "Detail"],
            ["数据集", "步骤", "样本", "受试者/ID", "说明"],
        ),
        "",
        "GSE40553 的 28 个排除样本没有可匹配的 LTB/PTB 治疗时间标签。"
        "GSE89403 的 914 条 GSM 中包含对照和技术重复；最终分析使用 312 个 subject-time 样本。",
        "",
        "### 2.2 时间点",
        "",
        md_table(
            time,
            ["Dataset", "AssignedGroup", "TimePoint", "TimeGroup", "Samples"],
            ["数据集", "组别", "时间点", "时间/结局", "样本"],
        ),
        "",
        "### 2.3 GSE40553 地区批次",
        "",
        md_table(
            batch,
            ["Dataset", "AssignedGroup", "Batch", "Samples"],
            ["数据集", "组别", "地区批次", "样本"],
        ),
        "",
        "UK 样本仅提供 7 个基线活动性 TB 样本；其余可用的纵向时间点来自 South Africa，"
        "因此地区与研究时相不能完全分离，正文已明确披露。",
        "",
        "### 2.4 排除 UK 基线样本的 South Africa 内部敏感性分析",
        "",
        "该分析只使用 28 个 South Africa 活动性 TB 基线样本和 38 个 South Africa LTB 样本，"
        "在排除 7 个 UK 基线样本后重新拟合疾病与治疗对比。",
        "",
        md_table(
            sa_overlap,
            ["Module", "Original_genes", "SA_only_genes", "Overlap", "Percent_original_retained", "Jaccard"],
            ["模块", "原基因数", "South Africa 基因数", "交集", "原模块保留率", "Jaccard"],
        ),
        "",
        md_table(
            sa_models,
            ["module", "contrast", "estimate", "SE", "p.value", "overall_time_p"],
            ["模块", "对比", "估计值", "SE", "P", "总体时间 P"],
        ),
        "",
        "South Africa 重建的恢复核心在 GSE31348 治疗过程中仍显著升高；"
        "原 171 个恢复核心中 161 个被保留（94.2%），原 175 个疾病驱动基因中 161 个被保留（92.0%）。",
        "",
        "South-Africa-only 候选 TF 对比：",
        "",
        md_table(
            sa_tf,
            [
                "TF", "disease_effect_SA", "disease_FDR_SA",
                "treatment_effect_SA", "treatment_FDR_SA",
            ],
            ["TF", "疾病效应", "疾病 FDR", "治疗效应", "治疗 FDR"],
        ),
        "",
        "South Africa 恢复核心的前列 GO 条目：",
        "",
        md_table(
            sa_go.head(10),
            ["Description", "GeneRatio", "FoldEnrichment", "pvalue", "p.adjust", "Count"],
            ["GO 条目", "基因比例", "富集倍数", "P", "FDR", "基因数"],
        ),
        "",
        "### 2.5 GSE326212 独立受试者审计",
        "",
        md_table(
            subject,
            ["group_norm", "samples", "unique_subjects", "max_samples_per_subject"],
            ["组别", "样本", "唯一受试者", "每人最大样本数"],
        ),
        "",
        "主分析中没有同一受试者贡献多个样本，因此患者级 pseudobulk 的统计单位与受试者一致。",
        "",
    ]


def section_bulk() -> list[str]:
    module40553 = read_tsv("gse40553_module_mixed_models.tsv")
    module31348 = read_tsv("gse31348_module_mixed_models.tsv")
    model89403 = read_tsv("gse89403_mixed_model_anova.tsv")
    model89403 = model89403[
        model89403["term"].isin(["Group", "Group:TimePoint"])
        & model89403["feature"].isin([
            "PAX5", "TBX21", "NFATC2", "PAX5_target_module",
            "TBX21_target_module", "B_cell_marker_module",
            "TNK_marker_module",
        ])
    ]
    restoration = (RESULTS / "restoration_core_genes.txt").read_text(encoding="utf-8").splitlines()
    disease = (RESULTS / "disease_driver_genes.txt").read_text(encoding="utf-8").splitlines()
    global_reversal = read_tsv("gse40553_global_reversal_correlation.tsv")
    return [
        "## 2. Bulk 血液纵向转录组结果",
        "",
        "### 3.1 恢复模块与疾病驱动模块",
        "",
        f"- GSE40553 定义 {len(restoration)} 个 restoration-core 基因。",
        f"- GSE40553 定义 {len(disease)} 个 disease-driver 基因。",
        "- 恢复核心的入选方向为：活动性 TB 中较低，治疗后升高。",
        "- 疾病驱动模块的入选方向相反：活动性 TB 中较高，治疗后下降。",
        "",
        md_table(
            global_reversal,
            [
                "Comparison", "N_genes", "Pearson_r", "Pearson_P",
                "Spearman_rho", "Spearman_P",
            ],
            ["比较", "基因数", "Pearson r", "Pearson P", "Spearman rho", "Spearman P"],
        ),
        "",
        "全基因组疾病效应与治疗效应呈负相关，支持双向逆转；"
        "这与 duplicateCorrelation 的 0.075 共识相关参数不是同一个指标。",
        "",
        "### 3.2 GSE40553 模块轨迹",
        "",
        md_table(
            module40553,
            ["feature", "contrast", "estimate", "SE", "p.value", "overall_time_p"],
            ["模块", "对比", "估计值", "SE", "P", "总体时间 P"],
        ),
        "",
        "GSE40553 的恢复核心从基线到 6 月持续升高，疾病驱动模块持续下降；"
        "12 月仍与基线不同，提示部分恢复在治疗结束后仍可持续。",
        "",
        "### 3.3 GSE31348 方向复现",
        "",
        md_table(
            module31348,
            ["feature", "contrast", "estimate", "SE", "p.value", "FDR"],
            ["模块", "对比", "估计值", "SE", "P", "FDR"],
        ),
        "",
        "GSE31348 是模型开发队列，不能作为独立验证。其恢复核心在 Week 1 已上升，"
        "疾病驱动模块同步下降，并在 Week 26 达到最大差异。",
        "",
        "### 3.4 GSE89403 结局关联模型",
        "",
        md_table(
            model89403,
            ["feature", "adjustment", "term", "F value", "Pr(>F)", "FDR"],
            ["特征", "校正", "项", "F", "P", "FDR"],
        ),
        "",
        "B-cell marker 的 Group×Time 原始 P=0.036，但 FDR=0.14。"
        "因此正文只把它解释为治疗过程中的可重复恢复，不把结局关联称为已验证的独立标志。",
        "",
    ]


def section_tf() -> list[str]:
    screen = read_tsv("unbiased_tf_cross_cohort_screen.tsv")
    top_screen = screen[screen["evidence_count"] >= 3].sort_values(
        ["evidence_count", "GSE40553_treatment_FDR"], ascending=[False, True]
    )
    robust = read_tsv("candidate_composition_robustness.tsv")
    top_robust = robust[robust["evidence_count"] >= 3].sort_values(
        ["adjusted_evidence_count", "evidence_count"], ascending=[False, False]
    )
    comp = read_tsv("tf_activity_composition_models.tsv")
    comp = comp[comp["TF"].isin(["MAF", "TBX21", "PAX5", "NFATC2", "STAT4"])]
    scale = read_tsv("quantiseq_scale_correction_comparison.tsv")
    scale = scale[scale["dataset"].eq("GSE89403")]
    hg = read_tsv("gse89403_hg19_hg38_quantiseq_comparison.tsv")
    gene_comp_summary = read_tsv("gene_composition_attenuation_summary.tsv")
    gene_comp_cross = read_tsv("gene_composition_attenuation_cross_cohort.tsv")
    gene_comp_floor = read_tsv(
        "gene_composition_attenuation_cross_cohort_by_floor.tsv"
    )
    return [
        "## 3. 系统性 TF 筛查与组成校正",
        "",
        "### 4.1 271 个 TF 的跨队列筛查",
        "",
        "筛选使用四个方向标准：GSE40553 治疗后升高、GSE31348 Week 26 升高、"
        "GSE89403 poor outcome 中较低、paired donor 与 restoration module 正相关。"
        "每个标准要求 FDR<0.05。",
        "",
        md_table(
            top_screen,
            [
                "TF", "GSE40553_treatment_logFC", "GSE40553_treatment_FDR",
                "GSE31348_week26_logFC", "GSE31348_week26_FDR",
                "GSE89403_poor_vs_success_logFC", "GSE89403_outcome_FDR",
                "donor_rho_restoration", "donor_FDR_restoration", "evidence_count",
            ],
            [
                "TF", "GSE40553 logFC", "FDR", "GSE31348 logFC", "FDR",
                "GSE89403 logFC", "FDR", "供者 rho", "FDR", "命中标准数",
            ],
        ),
        "",
        "MAF 是唯一满足四项标准的 TF；TBX21 满足三项，PAX5 和 NFATC2 各满足两项。"
        "这是一个候选筛选结果，不等同于独立验证。",
        "",
        "### 4.2 组成校正后的候选稳健性",
        "",
        md_table(
            top_robust,
            [
                "TF", "evidence_count", "adjusted_evidence_count",
                "GSE40553_adjusted_FDR", "GSE31348_adjusted_FDR",
                "GSE89403_adjusted_FDR",
            ],
            ["TF", "原始命中数", "校正后命中数", "GSE40553 FDR", "GSE31348 FDR", "GSE89403 FDR"],
        ),
        "",
        "MAF 在 GSE89403 中保留结局关联（校正后 FDR=0.046），"
        "但在 GSE40553/GSE31348 的组成校正后丢失治疗方向支持。"
        "没有 TF 在三个队列中同时保留校正后支持。",
        "",
        "### 4.3 原候选 TF 的组成解释",
        "",
        md_table(
            comp,
            ["dataset", "TF", "delta_R2", "LRT_p", "FDR"],
            ["数据集", "TF", "组成增量 R²", "似然比 P", "FDR"],
        ),
        "",
        "TBX21 和 PAX5 的模型在三个队列中均因加入组成 PCs 而显著改善；"
        "这支持“表观 TF 活性恢复至少部分来自细胞组成变化”的主结论。",
        "",
        "### 4.4 逐基因组成衰减的探索性检验",
        "",
        md_table(
            gene_comp_summary,
            [
                "dataset", "module", "n_genes", "median_attenuation",
                "prop_attenuation_ge_25pct", "prop_attenuation_ge_50pct",
                "prop_sign_reversed", "prop_composition_FDR_05",
            ],
            [
                "数据集", "模块", "基因数", "中位衰减", "衰减≥25% 比例",
                "衰减≥50% 比例", "方向反转比例", "组成联合 FDR<0.05 比例",
            ],
        ),
        "",
        "跨队列逐基因衰减一致性：",
        "",
        md_table(
            gene_comp_cross,
            ["module", "n_genes", "spearman_attenuation"],
            ["模块", "基因数", "跨队列 Spearman"],
        ),
        "",
        "设置原始治疗效应下限后的跨队列一致性：",
        "",
        md_table(
            gene_comp_floor,
            ["module", "effect_floor", "n_genes", "spearman_attenuation"],
            ["模块", "效应下限", "基因数", "跨队列 Spearman"],
        ),
        "",
        "逐基因衰减比例在队列间不稳定，不能解释为稳定的 cell-intrinsic fraction；"
        "因此该分析只作为探索性审计，不升级为新的主结论。",
        "",
        "### 3.4 quanTIseq 尺度及 hg19/hg38 注释敏感性",
        "",
        "旧版曾把 log2 表达直接传入 quanTIseq，依赖函数内部的自动反 log 逻辑；"
        "GSE89403 的官方 raw counts 由 hg19 比对和计数产生，因此主分析必须使用 hg19 基因长度。",
        "",
        md_table(
            scale,
            ["cell_type", "spearman_old_vs_nonlog", "median_abs_difference"],
            ["细胞类型", "旧图/新版 Spearman", "中位绝对差"],
        ),
        "",
        "hg19 主分析与 hg38 敏感性分析的细胞类型一致性如下：",
        "",
        md_table(
            hg,
            ["cell_type", "spearman_hg19_vs_hg38", "median_abs_fraction_difference"],
            ["细胞类型", "hg19/hg38 Spearman", "中位绝对分数差"],
        ),
        "",
        "B cell、Neutrophil、NK 和 CD8 T 的 hg19/hg38 一致性很高；"
        "Monocyte 和 Myeloid DC 的注释较不稳定。因此 hg19 是唯一主分析，"
        "hg38 仅用于显示注释不一致可能造成的影响。",
        "",
    ]


def section_score() -> list[str]:
    score = read_tsv("gse89403_score_external_auc.tsv")
    incremental = read_tsv("tbri_incremental_value.tsv")
    standardization = read_tsv("gse89403_tbri_standardization_auc.tsv")
    corr_components = read_tsv("gse89403_tbri_bcell_correlation_components.tsv")
    marker_check = read_tsv("tbri_independent_bcell_marker_correlation.tsv")
    signature_overlap = read_tsv("tbri_quantiseq_signature_overlap.tsv")
    return [
        "## 4. TBRI 与结局关联",
        "",
        "### 5.1 外部队列的 TBRI AUC",
        "",
        md_table(
            score,
            ["TimePoint", "n_success", "n_poor", "AUC", "lower_95", "upper_95", "brier_score"],
            ["时间点", "成功", "不良结局", "AUC", "95% 下限", "95% 上限", "Brier"],
        ),
        "",
        "Week 1 AUC 约 0.73；Week 4 的置信区间包含 0.5；Week 24 AUC 较高但属于接近疗程结束的时点，"
        "不能解释为早期预测。",
        "",
        "### 5.2 TBRI 是否超过 B-cell fraction",
        "",
        md_table(
            incremental,
            ["TimePoint", "n", "n_poor", "AUC_Bcell_only", "AUC_Bcell_plus_TBRI", "delta_AUC", "DeLong_p", "LRT_p"],
            ["时间点", "n", "不良结局", "仅 B-cell AUC", "B-cell+TBRI", "ΔAUC", "DeLong P", "LRT P"],
        ),
        "",
        "加入 TBRI 后 ΔAUC 为 0.002–0.037，所有 DeLong 检验均不显著。"
        "因此 TBRI 没有显示超过 B-cell fraction 的明确增量价值。",
        "",
        "### 5.3 标准化敏感性",
        "",
        md_table(
            standardization,
            ["task", "predictor", "TimePoint", "n_success", "n_poor", "AUC", "lower_95", "upper_95"],
            ["分析", "评分", "时间点", "成功", "不良结局", "AUC", "95% 下限", "95% 上限"],
        ),
        "",
        "队列内标准化和 GSE31348 冻结均值/SD 的 AUC 接近，但评分仍是回顾性、跨平台标准尚未解决；"
        "患者内基线变化评分的 AUC 较低，未形成早期监测证据。",
        "",
        "### 5.4 TBRI 与 B-cell fraction 的患者间及患者内相关",
        "",
        md_table(
            corr_components,
            ["Component", "Method", "N_samples", "N_participants", "Estimate", "P_value"],
            ["成分", "方法", "样本数", "患者数", "相关系数", "P"],
        ),
        "",
        "312 个重复时间点样本的合并相关系数仅作描述；82 名患者之间的相关和患者内重复测量相关"
        "分别计算后仍均为强相关。该结果说明 TBRI 同时反映患者间稳定差异和治疗期间患者内变化，"
        "但由于 CD79A 等评分基因本身也是 B-cell 特征，不能把相关解释为独立验证。",
        "",
        "进一步使用不含 CD79A/EBF1 的 MS4A1/CD19/CD79B 标记模块复核：",
        "",
        md_table(
            marker_check,
            ["Component", "Method", "N_samples", "N_participants", "Estimate", "P_value"],
            ["成分", "方法", "样本数", "患者数", "相关系数", "P"],
        ),
        "",
        "TBRI 基因是否存在于 quanTIseq 固定 153-gene signature：",
        "",
        md_table(
            signature_overlap,
            ["Gene", "In_quanTIseq_fixed_signature", "B_cell_weight"],
            ["基因", "是否在固定签名", "B-cell 权重"],
        ),
        "",
        "即使排除 TBRI 的 CD79A/EBF1，TBRI 仍与独立 B-cell marker 强相关。"
        "但只有 EBF1 与固定 quanTIseq 签名直接重叠，因此不能简单说相关完全由共享基因制造；"
        "更准确的结论是 TBRI 为 B-cell-associated proxy，而非独立调控轴。",
        "",
    ]


def section_single_cell() -> list[str]:
    loc = read_tsv("gse192483_candidate_tf_localization.tsv")
    top_loc = loc.sort_values(["TF", "median_mean_expression"], ascending=[True, False]).groupby(
        "TF", as_index=False
    ).head(1)
    models = read_tsv("gse326212_candidate_tf_models.tsv")
    top_models = models[models["comparison"].eq("PTB_vs_IGRAneg")].sort_values("p_value").head(20)
    robustness = read_tsv("gse326212_stat4_model_robustness.tsv")
    threshold = read_tsv("gse326212_stat4_raw_threshold_models.tsv")
    lopo = read_tsv("gse326212_stat4_leave_one_project_out.tsv")
    myeloid = read_tsv("gse326212_stat4_myeloid_composite.tsv")
    myc = read_tsv("mtb_pbmc_celltype_time_models.tsv")
    key_pb = myc[(myc["gene"].isin(["TBX21", "PAX5", "EBF1", "CD79A", "IFNG", "CD40LG"]))]
    return [
        "## 5. 体内与体外单细胞结果",
        "",
        "### 6.1 GSE192483 肺组织定位",
        "",
        md_table(
            top_loc,
            ["TF", "CellType", "n_patients", "median_mean_expression", "median_percent_expressing", "top_in_vivo_compartment"],
            ["TF", "最高表达细胞类型", "患者数", "中位表达", "中位阳性率", "体内主要区间"],
        ),
        "",
        "GSE192483 仅用于细胞类型定位；没有健康对照或治疗时间序列，细胞级 P 值不用于治疗结论。",
        "",
        "### 6.2 GSE326212 患者级 BAL 疾病态比较",
        "",
        md_table(
            top_models,
            ["TF", "celltype_coarse", "estimate", "std_error", "p_value", "FDR"],
            ["TF", "细胞类型", "估计值", "SE", "P", "FDR"],
        ),
        "",
        "最稳健的跨组信号是 Mo-derived AM 中的 STAT4。它在 PTB 中高于 IGRA− 和 IGRA+；"
        "其余候选没有同样一致的跨组支持。",
        "",
        "### 6.3 STAT4 模型稳健性",
        "",
        md_table(
            robustness,
            ["analysis", "outcome", "estimate", "std_error", "lower_95", "upper_95", "p_value"],
            ["模型", "结局", "估计值", "SE", "95% 下限", "95% 上限", "P"],
        ),
        "",
        "STAT4 的 pseudobulk、项目固定效应、项目随机效应和细胞数校正模型均保持同向；"
        "原始 count≥1 和平均 count 也支持该方向，但 count≥2 衰减。",
        "",
        "### 6.4 STAT4 原始 counts 阈值敏感性",
        "",
        md_table(
            threshold,
            ["outcome", "estimate", "std_error", "lower_95", "upper_95", "p_value"],
            ["指标", "估计值", "SE", "95% 下限", "95% 上限", "P"],
        ),
        "",
        "### 6.5 STAT4 留一项目分析",
        "",
        f"共 {len(lopo)} 个项目逐一被排除，全部估计为正；"
        f"{sum(lopo['estimate'] > 0)}/{len(lopo)} 为正，"
        f"{sum(lopo['p_value'] < 0.05)}/{len(lopo)} 名义 P<0.05。",
        "",
        "### 6.6 跨髓系亚群复核",
        "",
        md_table(
            myeloid,
            ["analysis", "outcome", "n_IGRAneg_samples", "n_PTB_samples", "estimate", "std_error", "p_value"],
            ["模型", "结局", "IGRA− 样本", "PTB 样本", "估计值", "SE", "P"],
        ),
        "",
        "### 6.7 体外 Mtb 刺激：Supplementary Figure S3",
        "",
        md_table(
            key_pb,
            ["contrast", "gene", "cell_type", "estimate", "p.value", "FDR", "n_observations", "n_donors"],
            ["对比", "基因", "细胞类型", "估计值", "P", "FDR", "观察", "供者"],
        ),
        "",
        "这部分只作为急性刺激情境：TBX21 在刺激后的 NK/CD8 T 中下降，"
        "PAX5 在部分 B-cell 时间点下降；它不能替代或否定体内治疗恢复。",
        "",
    ]


def section_metanet() -> list[str]:
    cand = read_tsv("metanet_candidate_cross_cohort.tsv")
    overlap = read_tsv("metanet_module_overlap.tsv")
    direct = read_tsv("metanet_direct_edges.tsv")
    return [
        "## 6. METANet 结果",
        "",
        "### 7.1 候选上游 TF",
        "",
        md_table(
            cand.head(12),
            ["TF", "METANet_targets", "best_METANet_rank", "best_METANet_score", "FDR_GSE40553", "FDR_GSE31348", "FDR_GSE89403"],
            ["TF", "预测靶点", "最佳排名", "分数", "GSE40553 FDR", "GSE31348 FDR", "GSE89403 FDR"],
        ),
        "",
        "MAF→PAX5 是当前最相关的候选先验，但它不是疾病特异性验证。",
        "",
        "### 7.2 模块重叠",
        "",
        md_table(
            overlap,
            ["network", "TF", "module", "overlap", "fold_enrichment", "p_value", "FDR"],
            ["网络", "TF", "模块", "重叠", "富集倍数", "P", "FDR"],
        ),
        "",
        "组织特异网络的模块富集没有通过 FDR；因此不能把 METANet 写成完整 TBX21–PAX5 轴的验证。",
        "",
        "### 7.3 直接靶边",
        "",
        md_table(
            direct.head(12),
            ["network", "TF", "target", "edge_present", "rank_within_TF", "outgoing_edges"],
            ["网络", "TF", "靶点", "边", "TF 内排名", "出边数"],
        ),
        "",
        "较高排名的 TBX21–IFNG、PAX5–CD79A 只支持候选先验；未达到可用于因果叙述的独立验证。",
        "",
    ]


def section_inventory() -> list[str]:
    content = [
        "## 7. 指标解释与不可支持的主张",
        "",
        "### 可以支持的表述",
        "",
        "- B-cell-associated transcriptional recovery 在治疗过程中可重复。",
        "- 推断的 TF 活性变化与细胞组成变化有大量重叠。",
        "- 没有证据支持一个跨队列、组成独立的 TF master regulator。",
        "- TBRI 有回顾性结局关联，但没有超过 B-cell fraction 的明确增量价值。",
        "- STAT4 阳性 Mo-derived AM 比例在活动性 TB 中较高，并在批次和分层敏感性分析中保持同向。",
        "",
        "### 不能支持的表述",
        "",
        "- TBX21 或 PAX5 是上游因果主轴。",
        "- METANet 已经验证该网络。",
        "- TBRI 是 deployable early-monitoring test。",
        "- 体外 Mtb 刺激可以否定体内治疗恢复。",
        "- GSE326212 是治疗前后纵向单细胞队列。",
        "- GSE312606 已完成患者配对治疗分析。",
        "",
        "## 8. 文件与表索引",
        "",
        f"- {file_link('docs/Complete_Results_Report.md', '本报告')}",
        f"- {file_link('docs/markdown/Manuscript_Inflammation_Research.md', '主稿 Markdown')}",
        f"- {file_link('docs/markdown/Abstract.md', '摘要')}",
        f"- {file_link('docs/markdown/Cover_Letter.md', 'Cover letter')}",
        f"- {file_link('docs/External_Review_Response_20261029.md', '第二轮审稿意见回复')}",
        f"- {file_link('docs/Expert_Discussion_Publication_Assessment.md', '发表可行性专家讨论')}",
        "",
    ]
    return content


def supplement_inventory() -> list[str]:
    workbook = SUPPLEMENT / "Supplementary_Tables.xlsx"
    sheet_purpose = {
        "S1_Dataset_Manifest": "数据集角色、样本基数、允许推断和限制",
        "S2_Sample_Counts": "各队列时间点样本数",
        "S3_Restoration_Core": "恢复核心基因列表",
        "S4_Disease_Driver": "疾病驱动基因列表",
        "S5_TF_Activity": "GSE40553 的 271 个 TF activity 结果",
        "S6_TF_Target_Enrichment": "TF target set 与模块的 Fisher 富集",
        "S7_GSE40553_Donor_Correlations": "供者配对 TF/模块相关",
        "S8_GSE31348_Module_Models": "GSE31348 模块纵向模型",
        "S9_GSE31348_TF_Models": "GSE31348 关键 TF 纵向模型",
        "S10_GSE89403_Model_ANOVA": "GSE89403 主模型和细胞标记校正",
        "S11_GSE89403_Group_Contrasts": "GSE89403 各时间点组间对比",
        "S12_GSE89403_Score_AUC": "TBRI 的结局 AUC 和区间",
        "S13_METANet_Direct_Edges": "METANet 候选直接边",
        "S14_METANet_Module_Overlap": "METANet 模块重叠富集",
        "S15_METANet_Upstream_Candidates": "METANet 候选上游 TF",
        "S16_TF_Composition_Models": "TF 组成校正混合模型",
        "S17_TBRI_Incremental_Value": "TBRI 相对 B-cell fraction 的增量 AUC",
        "S18_CellType_Pseudobulk_Models": "急性刺激的 donor-aware cell-type 模型",
        "S19_CellType_Percent_Positive": "急性刺激的细胞类型阳性率",
        "S20_Systematic_TF_Screen": "271 个 TF 四方向系统筛查",
        "S21_TF_Composition_Robustness": "候选 TF 的组成校正稳健性",
        "S22_GSE326212_TF_Models": "GSE326212 患者级候选 TF 模型",
        "S23_GSE192483_TF_Localization": "肺组织候选 TF 细胞类型定位",
        "S24_GSE326212_STAT4_Robustness": "STAT4 项目/深度/细胞数敏感性",
        "S25_GSE326212_STAT4_LOPO": "STAT4 留一项目分析",
        "S26_Cohort_Inclusion_Flow": "三队列样本纳入流程",
        "S27_Cohort_Timepoint_Counts": "各时间点样本计数",
        "S28_GSE40553_Batch_Groups": "GSE40553 地区/批次 x 组别",
        "S29_QuanTIseq_Scale_Audit": "旧尺度与显式非对数/TPM 比较",
        "S30_TBRI_Standardization_AUC": "TBRI 标准化敏感性 AUC",
        "S31_GSE326212_Subject_Unit": "GSE326212 受试者单位审计",
        "S32_GSE326212_STAT4_Thresholds": "STAT4 raw-count ≥1/≥2 阈值敏感性",
        "S33_GSE89403_HG19_HG38": "GSE89403 hg19/hg38 基因长度注释敏感性",
        "S34_GSE40553_SA_Overlap": "排除 UK 基线样本后原模块与 South Africa 模块的交集",
        "S35_GSE31348_SA_Core_Models": "South Africa 重建核心在 GSE31348 的纵向反应",
        "S36_GSE40553_SA_Restoration_GO": "South Africa 恢复核心的 GO 富集",
        "S37_GSE40553_SA_Driver_GO": "South Africa 疾病驱动模块的 GO 富集",
        "S38_TBRI_BCell_Corr": "TBRI 与 B-cell fraction 的合并、患者间及患者内相关",
        "S39_Gene_Comp_Attenuation": "逐基因组成衰减模型",
        "S40_Gene_Comp_Summary": "逐基因组成衰减模块汇总",
        "S41_Gene_Comp_CrossCohort": "不同原始效应阈值下的跨队列衰减一致性",
        "S42_Gene_Comp_Floor": "原始治疗效应下限敏感性",
        "S43_TBRI_BCell_Independent": "排除 CD79A/EBF1 的 B-cell marker 复核",
        "S44_QuanTIseq_Overlap": "TBRI 基因与 quanTIseq 固定签名重叠",
        "S45_Global_Reversal": "GSE40553 全基因组双向逆转相关",
        "S46_QuanTIseq_Diag": "quanTIseq 分数和值与取值诊断",
        "S47_SA_Candidate_TF": "South-Africa-only 候选 TF 治疗与疾病对比",
    }
    xls = pd.ExcelFile(workbook)
    rows = []
    for sheet in xls.sheet_names:
        frame = pd.read_excel(xls, sheet_name=sheet)
        rows.append({
            "Sheet": sheet,
            "Rows": len(frame),
            "Columns": len(frame.columns),
            "Purpose": sheet_purpose.get(sheet, ""),
        })
    rows = pd.DataFrame(rows)
    return [
        "## 9. 补充结果表清单",
        "",
        md_table(
            rows,
            ["Sheet", "Rows", "Columns", "Purpose"],
            ["工作表", "数据行", "列数", "用途"],
        ),
        "",
        "完整行级结果位于 "
        f"{file_link('supplement/Supplementary_Tables.xlsx', 'Supplementary_Tables.xlsx')}。"
        "本报告对每个结果域解释关键统计和解释边界，补充表保留全部行级结果。",
        "",
    ]


def main() -> None:
    lines = [
        "# TB 治疗转录组研究：完整结果说明",
        "",
        "日期：2026-10-29",
        "",
        "本文档汇总当前项目可读取的全部主要结果，包括主稿 Results、补充表、"
        "组成校正、TBRI、体内 BAL、体外刺激、METANet 和数据审计。"
        "数字由当前结果文件生成；行级完整结果见 Supplementary_Tables.xlsx。",
        "",
        "## 0. 结果分层",
        "",
        "| 层级 | 内容 | 当前判断 |",
        "|---|---|---|",
        "| 直接观察 | 治疗轨迹、结局组差异、细胞类型表达 | 可描述 |",
        "| 统计关联 | 模块、TF activity、TBRI、STAT4 | 可报告，但不能自动升级为因果 |",
        "| 计算推断 | DoRothEA/VIPER、METANet、DoRothEA targets | 只能作为候选或先验 |",
        "| 体内情境 | GSE326212 BAL、GSE192483 lung | 横断面情境，不是治疗验证 |",
        "| 实验机制 | 无扰动/rescue | 尚不可声称因果机制 |",
        "",
    ]
    for fn in [
        section_cohorts,
        section_bulk,
        section_tf,
        section_score,
        section_single_cell,
        section_metanet,
        section_inventory,
        supplement_inventory,
    ]:
        lines.extend(fn())

    OUT.write_text("\n".join(lines).rstrip() + "\n", encoding="utf-8")
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
