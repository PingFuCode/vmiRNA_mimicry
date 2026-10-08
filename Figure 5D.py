import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import warnings

from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_curve, auc, accuracy_score
from sklearn.preprocessing import StandardScaler
from scipy.stats import fisher_exact

warnings.filterwarnings("ignore")

# ==================== 0. 全局绘图参数 ====================
plt.rcParams["font.family"] = "Arial"
plt.rcParams["font.size"] = 10
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42


# ==================== 1. 读取数据 ====================
file_path = r"D:\0work\0wholetransctiptome\2sRNAminic\15TCGA\EBVmiRNAdata\ebv_patients_with_virus_simulated_mirnas.csv"
df = pd.read_csv(file_path)

print("===== Basic information =====")
print(f"Raw input samples: {len(df)}")
print("Columns:")
print(df.columns.tolist())


# 保留完整数据，仅用于确定31个 EBV-like hmiRNA 特征列
all_samples_df = df.copy()


# ==================== 2. 筛选原发肿瘤并重建 stage_group ====================
# Stage I / II 及其带字母亚分期 -> Early
# Stage III / IV 及其带字母亚分期 -> Late
required_cols = {"sample_id", "stage"}
missing_cols = required_cols - set(df.columns)
if missing_cols:
    raise ValueError(f"Required columns were not found: {sorted(missing_cols)}")

# TCGA sample type 01A 为原发肿瘤；排除两个 11A 正常组织样本
before_sample_filter = len(df)
df = df[
    df["sample_id"].astype(str).str.contains("-01A", regex=False, na=False)
].copy()

stage_series = df["stage"].astype(str).str.strip()
early_mask = stage_series.str.contains(
    r"^Stage (?:I|II)[ABC]?$",
    regex=True,
    na=False
)
late_mask = stage_series.str.contains(
    r"^Stage (?:III|IV)[ABC]?$",
    regex=True,
    na=False
)

# 不沿用输入文件中的旧 stage_group，始终根据原始病理分期重建
df["stage_group"] = np.where(
    early_mask, "Early",
    np.where(late_mask, "Late", np.nan)
)

# 删除没有明确 Early/Late 分组的样本
before_drop = len(df)
df = df[df["stage_group"].isin(["Early", "Late"])].copy()
after_drop = len(df)

print("\n===== Stage group information =====")
print(f"Primary tumor samples retained: {len(df)} / {before_sample_filter}")
print(f"Samples before stage filtering: {before_drop}")
print(f"Samples after stage filtering: {after_drop}")
print(df["stage_group"].value_counts())


# ==================== 3. 提取并筛选 miRNA 列 ====================
# 默认前4列不是 miRNA
candidate_cols = all_samples_df.columns[4:].tolist()

# 去掉明显不应作为特征的列
exclude_cols = {"stage_group", "stage"}
candidate_cols = [c for c in candidate_cols if c not in exclude_cols]

# 只保留数值列
candidate_cols = [
    c for c in candidate_cols
    if pd.api.types.is_numeric_dtype(all_samples_df[c])
]

# 只保留至少一个样本中表达值 > 0 的 miRNA
expr_df = all_samples_df.loc[:, candidate_cols]
mirna_cols = expr_df.columns[(expr_df > 0).sum(axis=0) >= 1].tolist()

# 被删除的 miRNA
removed_mirnas = sorted(set(candidate_cols) - set(mirna_cols))

print("\n===== miRNA filtering =====")
print(f"Candidate numeric miRNAs before filtering: {len(candidate_cols)}")
print(f"Retained miRNAs (>0 in at least one sample): {len(mirna_cols)}")
print(f"Removed miRNAs (all values <= 0): {len(removed_mirnas)}")

if len(mirna_cols) == 0:
    raise ValueError("No miRNAs remained after filtering. Please check the input data.")


# ==================== 4. 按给定分组做 Fisher 精确检验 ====================
contingency_table = pd.DataFrame(
    [[2, 11], [9, 5]],
    index=["High expression", "Low expression"],
    columns=["Early", "Late"]
)

stage_counts = df["stage_group"].value_counts().reindex(
    ["Early", "Late"], fill_value=0
)
if not np.array_equal(contingency_table.sum(axis=0).to_numpy(), stage_counts.to_numpy()):
    raise ValueError(
        "The stage totals in the supplied expression table do not match "
        f"the analysis cohort: table={contingency_table.sum(axis=0).to_dict()}, "
        f"cohort={stage_counts.to_dict()}"
    )

# 原表方向（列为 Early, Late）对应 Early odds：OR = 0.101；
# 为便于论文解释，同时报告 High vs Low 的 Late-stage odds ratio（其倒数）。
odds_ratio_early, fisher_p = fisher_exact(
    contingency_table.to_numpy(),
    alternative="two-sided"
)
odds_ratio_late = 1.0 / odds_ratio_early
log_or_se = np.sqrt(np.sum(1.0 / contingency_table.to_numpy()))
odds_ratio_late_ci = np.exp(
    np.log(odds_ratio_late) + np.array([-1, 1]) * 1.96 * log_or_se
)

print("\n===== Expression group vs stage =====")
print(contingency_table)
print(f"Two-sided Fisher's exact P value: {fisher_p:.12g}")
print(f"Odds ratio for Early stage (High vs Low): {odds_ratio_early:.6f}")
print(f"Odds ratio for Late stage (High vs Low): {odds_ratio_late:.6f}")
print(
    "95% CI for Late-stage odds ratio: "
    f"{odds_ratio_late_ci[0]:.6f}-{odds_ratio_late_ci[1]:.6f}"
)


# ==================== 5. 构建特征矩阵和标签 ====================
X = df[mirna_cols].values
y = (df["stage_group"] == "Late").astype(int).values

print("\n===== Modeling input =====")
print(f"X shape: {X.shape}")
print(f"Positive samples (Late): {np.sum(y == 1)}")
print(f"Negative samples (Early): {np.sum(y == 0)}")


# ==================== 6. 整体模型：100次分层随机 8:2 划分 ====================
n_repeats = 100
mean_fpr = np.linspace(0, 1, 200)

tprs = []
aucs = []
accs = []

fig_mean, ax_mean = plt.subplots(figsize=(6, 5))

for i in range(n_repeats):
    X_train, X_test, y_train, y_test = train_test_split(
        X,
        y,
        test_size=0.2,
        random_state=42 + i,
        stratify=y
    )

    # 只在训练集上拟合标准化器
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    X_test_scaled = scaler.transform(X_test)

    rf = RandomForestClassifier(
        n_estimators=100,
        random_state=42 + i,
        n_jobs=-1
    )

    rf.fit(X_train_scaled, y_train)

    y_prob = rf.predict_proba(X_test_scaled)[:, 1]
    y_pred = rf.predict(X_test_scaled)

    fpr, tpr, _ = roc_curve(y_test, y_prob)
    this_auc = auc(fpr, tpr)
    this_acc = accuracy_score(y_test, y_pred)

    aucs.append(this_auc)
    accs.append(this_acc)

    interp_tpr = np.interp(mean_fpr, fpr, tpr)
    interp_tpr[0] = 0.0
    tprs.append(interp_tpr)

# 平均 ROC
mean_tpr = np.mean(tprs, axis=0)
mean_tpr[-1] = 1.0
mean_auc = np.mean(aucs)
std_auc = np.std(aucs)

# 只画一根平均 ROC 曲线
ax_mean.plot(
    mean_fpr,
    mean_tpr,
    lw=2.5,
    label=f"Mean ROC (AUC = {mean_auc:.3f} ± {std_auc:.3f})"
)

# 随机分类器基线，灰色虚线
ax_mean.plot(
    [0, 1], [0, 1],
    linestyle="--",
    color="gray",
    lw=1.5,
    label="Random Classifier"
)

ax_mean.set_xlim([0.0, 1.0])
ax_mean.set_ylim([0.0, 1.05])
ax_mean.set_xlabel("False Positive Rate (1 - Specificity)", fontsize=12, fontweight="bold")
ax_mean.set_ylabel("True Positive Rate (Sensitivity)", fontsize=12, fontweight="bold")
ax_mean.set_title("Mean ROC Curve from 100 Repeated Stratified 8:2 Splits", fontsize=13, fontweight="bold")
ax_mean.legend(loc="lower right", fontsize=8, frameon=True)
ax_mean.grid(False)

plt.tight_layout()

out_pdf_mean = r"D:\0work\0wholetransctiptome\2sRNAminic\15TCGA\EBVmiRNAdata\Random_Forest_100x_StratifiedSplit_MeanROC_8to2.pdf"
fig_mean.savefig(out_pdf_mean, dpi=300, bbox_inches="tight", format="pdf")

print("\n===== 100 repeated stratified 8:2 splits =====")
print(f"Mean ROC figure saved to: {out_pdf_mean}")
print(f"Mean AUC: {mean_auc:.4f} ± {std_auc:.4f}")
print(f"Mean accuracy: {np.mean(accs):.4f} ± {np.std(accs):.4f}")
print(f"Min AUC: {np.min(aucs):.4f}")
print(f"Max AUC: {np.max(aucs):.4f}")


# ==================== 7. 每个 miRNA：100次分层随机 8:2 划分，计算单变量 AUC ====================
n_repeats_mirna = 100
mean_fpr_single = np.linspace(0, 1, 200)

mirna_auc_dict = {}
mirna_tpr_dict = {}
mirna_all_auc_dict = {}

for mirna in mirna_cols:
    idx = mirna_cols.index(mirna)

    auc_list = []
    tpr_list = []

    for i in range(n_repeats_mirna):
        X_train, X_test, y_train, y_test = train_test_split(
            X,
            y,
            test_size=0.2,
            random_state=1000 + i,
            stratify=y
        )

        # 单独取该 miRNA 在测试集中的表达值
        single_scores = X_test[:, idx]

        fpr_single, tpr_single, _ = roc_curve(y_test, single_scores)
        auc_single = auc(fpr_single, tpr_single)

        # 若方向相反，则取反，统一按 AUC >= 0.5 计算
        if auc_single < 0.5:
            single_scores = -single_scores
            fpr_single, tpr_single, _ = roc_curve(y_test, single_scores)
            auc_single = auc(fpr_single, tpr_single)

        auc_list.append(auc_single)

        interp_tpr = np.interp(mean_fpr_single, fpr_single, tpr_single)
        interp_tpr[0] = 0.0
        tpr_list.append(interp_tpr)

    mirna_auc_dict[mirna] = {
        "mean_auc": np.mean(auc_list),
        "std_auc": np.std(auc_list),
        "median_auc": np.median(auc_list),
        "min_auc": np.min(auc_list),
        "max_auc": np.max(auc_list),
        "auc_gt_0.7_count": np.sum(np.array(auc_list) > 0.7),
        "auc_gt_0.8_count": np.sum(np.array(auc_list) > 0.8)
    }

    mirna_tpr_dict[mirna] = np.array(tpr_list)
    mirna_all_auc_dict[mirna] = auc_list


# ==================== 8. 导出每个 miRNA 的汇总表 ====================
mirna_auc_df = pd.DataFrame([
    {
        "miRNA": mirna,
        "mean_auc": mirna_auc_dict[mirna]["mean_auc"],
        "std_auc": mirna_auc_dict[mirna]["std_auc"],
        "median_auc": mirna_auc_dict[mirna]["median_auc"],
        "min_auc": mirna_auc_dict[mirna]["min_auc"],
        "max_auc": mirna_auc_dict[mirna]["max_auc"],
        "auc_gt_0.7_count": mirna_auc_dict[mirna]["auc_gt_0.7_count"],
        "auc_gt_0.8_count": mirna_auc_dict[mirna]["auc_gt_0.8_count"]
    }
    for mirna in mirna_cols
]).sort_values("mean_auc", ascending=False).reset_index(drop=True)

for col in ["mean_auc", "std_auc", "median_auc", "min_auc", "max_auc"]:
    mirna_auc_df[col] = mirna_auc_df[col].round(4)

print("\n===== Top miRNAs ranked by mean AUC across 100 repeated stratified 8:2 splits =====")
print(mirna_auc_df.head(10))

out_auc_csv = r"D:\0work\0wholetransctiptome\2sRNAminic\15TCGA\EBVmiRNAdata\All_miRNA_AUC_100repeats_8to2.csv"
out_auc_xlsx = r"D:\0work\0wholetransctiptome\2sRNAminic\15TCGA\EBVmiRNAdata\All_miRNA_AUC_100repeats_8to2.xlsx"

mirna_auc_df.to_csv(out_auc_csv, index=False, encoding="utf-8-sig")
mirna_auc_df.to_excel(out_auc_xlsx, index=False)

print(f"All miRNA AUC summary table saved to: {out_auc_csv}")
print(f"All miRNA AUC summary table saved to: {out_auc_xlsx}")


# ==================== 9. 导出每个 miRNA 每一次重复的 AUC 结果 ====================
all_auc_records = []

for mirna in mirna_cols:
    row = {"miRNA": mirna}
    auc_list = mirna_all_auc_dict[mirna]

    for j, auc_value in enumerate(auc_list, start=1):
        row[f"AUC_repeat_{j}"] = round(auc_value, 4)

    all_auc_records.append(row)

mirna_auc_all_df = pd.DataFrame(all_auc_records)

out_auc_all_xlsx = r"D:\0work\0wholetransctiptome\2sRNAminic\15TCGA\EBVmiRNAdata\All_miRNA_AUC_each_repeat_100times_8to2.xlsx"
mirna_auc_all_df.to_excel(out_auc_all_xlsx, index=False)

print(f"All per-repeat miRNA AUC table saved to: {out_auc_all_xlsx}")


# ==================== 10. 选取 mean AUC 最高的 top 3 miRNA ====================
top3_mirnas = mirna_auc_df.head(3)["miRNA"].tolist()

print("\n===== Top 3 miRNAs by mean AUC =====")
print(top3_mirnas)


# ==================== 11. 绘制 top 3 miRNA 的平均 ROC 曲线 ====================
fig_top3, ax_top3 = plt.subplots(figsize=(6, 5))

# Set1 调色板前三种颜色
set1_colors = plt.get_cmap("Set1").colors[:3]

for i, mirna in enumerate(top3_mirnas):
    tpr_array = mirna_tpr_dict[mirna]
    mean_tpr_single = np.mean(tpr_array, axis=0)
    mean_tpr_single[-1] = 1.0

    mean_auc_single = mirna_auc_dict[mirna]["mean_auc"]
    std_auc_single = mirna_auc_dict[mirna]["std_auc"]

    ax_top3.plot(
        mean_fpr_single,
        mean_tpr_single,
        lw=2,
        color=set1_colors[i],
        label=f"{mirna} (AUC = {mean_auc_single:.3f} ± {std_auc_single:.3f})"
    )

ax_top3.plot(
    [0, 1], [0, 1],
    linestyle="--",
    color="gray",
    lw=1.5,
    label="Random Classifier"
)

ax_top3.set_xlim([0.0, 1.0])
ax_top3.set_ylim([0.0, 1.05])
ax_top3.set_xlabel("False Positive Rate (1 - Specificity)", fontsize=12, fontweight="bold")
ax_top3.set_ylabel("True Positive Rate (Sensitivity)", fontsize=12, fontweight="bold")
ax_top3.set_title("Mean ROC Curves of Top 3 miRNAs Across 100 Repeated Splits", fontsize=13, fontweight="bold")
ax_top3.legend(loc="lower right", fontsize=8, frameon=True)
ax_top3.grid(False)

plt.tight_layout()

out_pdf_top3 = r"D:\0work\0wholetransctiptome\2sRNAminic\15TCGA\EBVmiRNAdata\Top3_miRNA_MeanROC_100repeats_8to2.pdf"
fig_top3.savefig(out_pdf_top3, dpi=300, bbox_inches="tight", format="pdf")

print(f"Top 3 miRNA mean ROC figure saved to: {out_pdf_top3}")

plt.show()
