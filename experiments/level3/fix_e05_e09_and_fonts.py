#!/usr/bin/env python3
"""Re-run E05 (with data-driven feasible M grid for Three-Child) & E09 (Korea feasible M grid),
and re-render all 6 figures with verified CJK font registration."""
from __future__ import annotations
import json
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager

# Register Noto Sans CJK SC directly so matplotlib finds Chinese glyphs cleanly
FONT_PATH = "/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc"
font_manager.fontManager.addfont(FONT_PATH)
cjk_font = font_manager.FontProperties(fname=FONT_PATH).get_name()
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = [cjk_font, "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False
plt.rcParams["figure.dpi"] = 180
plt.rcParams["axes.spines.top"] = False
plt.rcParams["axes.spines.right"] = False

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from cipd.causality.causal_bounds import solve_honest_smoothness_bounds

TABLES = ROOT / "results/tables"
FIGS = ROOT / "results/figures"
LOGS = ROOT / "results/logs"

df_cn = pd.read_csv(ROOT / "data/raw/china_annual_vital_1990_2025.csv")
df_kr = pd.read_csv(ROOT / "data/external/korea_vital_2014_2025.csv")

# 1. Re-run E05 with policy-specific feasible M grids (explicitly documenting M_min_pre)
sub = df_cn[df_cn["year"] >= 2005].copy().reset_index(drop=True)
yrs = sub["year"].to_numpy(dtype=int)
b_orig = sub["births_wan"].to_numpy(dtype=float)
b_rev = sub["revised_births_wan_post7th"].to_numpy(dtype=float)
y_lo = np.where(np.isnan(b_rev), b_orig * 0.99, np.minimum(b_orig, b_rev))
y_hi = np.where(np.isnan(b_rev), b_orig * 1.01, np.maximum(b_orig, b_rev))

policy_configs = [
    ("Selective Two-Child (2014-2015)", 2014, 2015, "C2", [20.0, 40.0, 60.0, 100.0]),
    ("Universal Two-Child (2016-2019)", 2016, 2019, "C2", [20.0, 40.0, 60.0, 100.0]),
    ("Three-Child Policy (2021-2024)",  2021, 2024, "C3", [160.0, 200.0, 250.0, 320.0]),
]
rows_e05 = []
for pname, start_y, end_y, grade, m_list in policy_configs:
    mask = (yrs >= start_y - 6) & (yrs <= end_y)
    pre_d2_max = float(np.max(np.abs(np.diff(b_orig[(yrs >= start_y - 6) & (yrs < start_y)], n=2))))
    for M in m_list:
        res = solve_honest_smoothness_bounds(
            yrs[mask], y_lo[mask], y_hi[mask],
            policy_start_year=start_y, eval_end_year=end_y,
            M=M, policy_name=pname, evidence_grade=grade,
        )
        rows_e05.append({
            "policy": pname,
            "pre_max_d2_wan": round(pre_d2_max, 1),
            "M_wan_per_yr2": M,
            "cum_effect_lower_wan": round(res.cum_effect_lower, 1),
            "cum_effect_upper_wan": round(res.cum_effect_upper, 1),
            "avg_annual_lower_wan": round(res.avg_annual_lower, 1),
            "avg_annual_upper_wan": round(res.avg_annual_upper, 1),
            "sign_robust_positive": res.sign_robust_positive,
            "evidence_grade": grade,
        })
df_e05 = pd.DataFrame(rows_e05)
df_e05.to_csv(TABLES / "E05_china_policy_bounds.csv", index=False)

# 2. Re-run E09 Korea with feasible M grid (pre-2024 max |d2| = 35.2k)
yrs_kr = df_kr["year"].to_numpy(dtype=int)
b_kr = df_kr["births_thousands"].to_numpy(dtype=float)
kr_pre_d2 = float(np.max(np.abs(np.diff(b_kr[yrs_kr < 2024], n=2))))
rows_e09 = []
for M in [35.0, 45.0, 60.0, 80.0]:
    res = solve_honest_smoothness_bounds(
        yrs_kr, b_kr * 0.995, b_kr * 1.005,
        policy_start_year=2024, eval_end_year=2025, M=M,
        policy_name="Korea_2024_2025_Rebound", evidence_grade="C2",
    )
    rows_e09.append({
        "country": "South Korea (Real KOSIS)",
        "pre_max_d2_k": round(kr_pre_d2, 2),
        "M_k_per_yr2": M,
        "cum_effect_lower_k": round(res.cum_effect_lower, 2),
        "cum_effect_upper_k": round(res.cum_effect_upper, 2),
        "sign_robust_positive": res.sign_robust_positive,
    })
df_e09 = pd.DataFrame(rows_e09)
df_e09.to_csv(TABLES / "E09_cross_country_korea_validation.csv", index=False)

# Load remaining tables for figures
df_e02_real = pd.read_csv(TABLES / "E02_real_census2010_masking.csv")
df_e02_dyn = pd.read_csv(TABLES / "E02_dynamic_flow_masking.csv")
df_e03 = pd.read_csv(TABLES / "E03_synthetic_causal_recovery.csv")
df_e06 = pd.read_csv(TABLES / "E06_shapley_decomposition.csv")
df_e07 = pd.read_csv(TABLES / "E07_voi_observation_selection.csv")
df_pol = pd.read_csv(TABLES / "E08_robust_policy_metrics.csv")
df_proj = pd.read_csv(TABLES / "E08_china_longterm_scenarios_2025_2070.csv")

# Re-render all 6 figures cleanly with zero missing glyphs
# Fig 1
fig, ax = plt.subplots(figsize=(9.5, 4.8))
ax.plot(df_cn["year"], df_cn["births_wan"], color="#1B4F72", lw=2.3, marker="o", ms=3.5, label="出生人口（同期统计公报，万人）")
ax.plot(df_cn["year"], df_cn["deaths_wan"], color="#922B21", lw=2.0, ls="--", label="死亡人口（万人）")
rev_mask = df_cn["revised_births_wan_post7th"].notna()
ax.scatter(df_cn.loc[rev_mask, "year"], df_cn.loc[rev_mask, "revised_births_wan_post7th"], color="#D35400", s=45, zorder=5, label="七普事后修订出生数（万人）")
for _, r in df_cn[rev_mask].iterrows():
    ax.plot([r["year"], r["year"]], [r["births_wan"], r["revised_births_wan_post7th"]], color="#D35400", lw=1.4, alpha=0.7)
for yr, txt in [(2014, "单独二孩\n(2014)"), (2016, "全面二孩\n(2016)"), (2021, "三孩政策\n(2021)"), (2025, "全国育儿补贴\n(2025)")]:
    ax.axvline(yr, color="#7F8C8D", ls=":", lw=1.2)
    ax.text(yr + 0.15, 2340, txt, fontsize=8.5, color="#2C3E50", va="top")
ax.set_title("图1  1990—2025年中国出生与死亡人口演变、政策节点及普查修订区间", fontsize=11.5, fontweight="bold")
ax.set_xlabel("年份"); ax.set_ylabel("人口数（万人）")
ax.set_ylim(550, 2500); ax.grid(alpha=0.25); ax.legend(frameon=False, loc="lower left")
fig.tight_layout(); fig.savefig(FIGS / "Fig01_china_population_policy_timeline.png"); plt.close(fig)

# Fig 2
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.5))
y_pos = np.arange(len(df_e02_real))
ax1.hlines(y_pos, df_e02_real["cipd_lower"] / 1e4, df_e02_real["cipd_upper"] / 1e4, color="#1F618D", lw=5, alpha=0.85, label="CIPD 识别区间")
ax1.axvline(df_e02_real["true_2nd_births_age25_39"].iloc[0] / 1e4, color="#C0392B", ls="--", lw=1.8, label="六普真实值")
ax1.scatter(df_e02_real["ipf_point_est"] / 1e4, y_pos + 0.15, color="#D68910", marker="D", s=38, label="IPF 点估计")
ax1.scatter(df_e02_real["maxent_point_est"] / 1e4, y_pos - 0.15, color="#1E8449", marker="^", s=42, label="MaxEnt 点估计")
ax1.set_yticks(y_pos); ax1.set_yticklabels(["Mask-A\n(仅总出生)", "Mask-B\n(仅年龄边际)", "Mask-C\n(仅孩次边际)", "Mask-F\n(年龄+孩次边际)"], fontsize=8.5)
ax1.set_xlabel("25—39岁二孩出生样本数（万人）"); ax1.set_title("(a) 2010六普真实35×3矩阵遮蔽识别边界", fontsize=10, fontweight="bold")
ax1.legend(frameon=False, fontsize=8, loc="lower right"); ax1.grid(axis="x", alpha=0.25)

ax2.fill_between(df_e02_dyn["mask_rate_pct"], df_e02_dyn["cipd_lower"], df_e02_dyn["cipd_upper"], color="#5DADE2", alpha=0.35, label="CIPD 动态流量识别区间")
ax2.plot(df_e02_dyn["mask_rate_pct"], df_e02_dyn["true_cum_2nd_births"], color="#C0392B", ls="--", lw=2, label="真实累计二孩出生")
ax2.set_xlabel("孩次观测缺失比例（%）"); ax2.set_ylabel("6期累计二孩出生数")
ax2.set_title("(b) 动态流量网络随缺失率变化的识别宽度", fontsize=10, fontweight="bold")
ax2.legend(frameon=False, fontsize=8.5); ax2.grid(alpha=0.25)
fig.tight_layout(); fig.savefig(FIGS / "Fig02_E02_masking_identification.png"); plt.close(fig)

# Fig 3
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.5, 4.6))
x = np.arange(len(df_e03))
ax1.bar(x - 0.18, df_e03["naive_wrong_sign_rate"] * 100, width=0.35, color="#C0392B", alpha=0.85, label="传统线性外推错误方向率")
ax1.bar(x + 0.18, df_e03["cipd_wrong_sign_rate"] * 100, width=0.35, color="#1F618D", alpha=0.9, label="CIPD 敏感性区间错误方向率")
ax1.set_xticks(x); ax1.set_xticklabels([f"DGP-{i+1}" for i in range(len(df_e03))], rotation=35, fontsize=8.5)
ax1.set_ylabel("错误方向判断率（%）"); ax1.set_title("(a) 10类合成机制下错误政策方向判断率对比", fontsize=10, fontweight="bold")
ax1.legend(frameon=False, fontsize=8.5); ax1.grid(axis="y", alpha=0.25)

ax2.bar(x, df_e03["cipd_coverage"] * 100, color="#2E86C1", alpha=0.85)
ax2.axhline(95.0, color="#1E8449", ls="--", lw=1.5, label="95% 目标参考线")
ax2.set_xticks(x); ax2.set_xticklabels([f"DGP-{i+1}" for i in range(len(df_e03))], rotation=35, fontsize=8.5)
ax2.set_ylim(0, 108); ax2.set_ylabel("CIPD 真实因果效应覆盖率（%）")
ax2.set_title("(b) 200次蒙特卡洛重复下的真实效应覆盖率", fontsize=10, fontweight="bold")
ax2.legend(frameon=False, fontsize=8.5); ax2.grid(axis="y", alpha=0.25)
fig.tight_layout(); fig.savefig(FIGS / "Fig03_E03_synthetic_causal_recovery.png"); plt.close(fig)

# Fig 4
fig, axes = plt.subplots(1, 3, figsize=(12.2, 4.3), sharey=False)
for ax, pname in zip(axes, df_e05["policy"].unique()):
    sub = df_e05[df_e05["policy"] == pname]
    ax.fill_between(sub["M_wan_per_yr2"], sub["avg_annual_lower_wan"], sub["avg_annual_upper_wan"], color="#2E86C1", alpha=0.30)
    ax.plot(sub["M_wan_per_yr2"], sub["avg_annual_lower_wan"], color="#1B4F72", lw=1.8, marker="o", ms=4)
    ax.plot(sub["M_wan_per_yr2"], sub["avg_annual_upper_wan"], color="#1B4F72", lw=1.8, marker="o", ms=4)
    ax.axhline(0.0, color="#C0392B", ls="--", lw=1.3)
    ax.set_title(pname, fontsize=9.5, fontweight="bold")
    ax.set_xlabel("背景趋势曲率偏离上界 M（万人/年²）", fontsize=8.5)
    ax.grid(alpha=0.25)
axes[0].set_ylabel("年均出生人口增量边界（万人/年）", fontsize=9)
fig.suptitle("图4  中国三项生育政策在平滑趋势偏离 Δ^SD(M) 与普查修订区间下的部分识别边界", fontsize=11, fontweight="bold", y=1.02)
fig.tight_layout(); fig.savefig(FIGS / "Fig04_E05_china_policy_sensitivity.png", bbox_inches="tight"); plt.close(fig)

# Fig 5
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.5))
ax1.plot(df_e07["budget"], df_e07["random_reduction_pct"], "o--", color="#7F8C8D", lw=1.8, label="随机选择")
ax1.plot(df_e07["budget"], df_e07["cheapest_reduction_pct"], "s--", color="#D68910", lw=1.8, label="最低成本优先")
ax1.plot(df_e07["budget"], df_e07["cipd_greedy_reduction_pct"], "^-", color="#1F618D", lw=2.3, label="CIPD 双轨迹 VOI 贪心")
ax1.plot(df_e07["budget"], df_e07["exhaustive_opt_reduction_pct"], "x:", color="#C0392B", lw=1.8, ms=7, label="穷举全局最优")
ax1.set_xlabel("观测获取成本预算"); ax1.set_ylabel("目标识别宽度缩减比例（%）")
ax1.set_title("(a) E07 目标导向数据价值选择与穷举最优对比", fontsize=10, fontweight="bold")
ax1.legend(frameon=False, fontsize=8.5); ax1.grid(alpha=0.25)

colors = ["#2E86C1" if safe else "#C0392B" for safe in df_pol["baseline_safe"]]
ax2.scatter(df_pol["annual_fiscal_cost_bn"], df_pol["worst_rel_loss_vs_baseline"], c=colors, s=80, zorder=4)
ax2.axhline(0.0, color="#1E8449", ls="--", lw=1.5, label="基线保护线 R(u) = 0")
for _, r in df_pol.iterrows():
    ax2.annotate(r["policy"].split(" ")[0], (r["annual_fiscal_cost_bn"] + 3, r["worst_rel_loss_vs_baseline"] + 0.008), fontsize=8)
ax2.set_xlabel("年化财政支出估算（亿元）"); ax2.set_ylabel("相对现行2025基线的最坏情形损失 R(u)")
ax2.set_title("(b) E08 基线保护型鲁棒政策前沿（蓝=满足基线保护）", fontsize=10, fontweight="bold")
ax2.legend(frameon=False, fontsize=8.5); ax2.grid(alpha=0.25)
fig.tight_layout(); fig.savefig(FIGS / "Fig05_E07_VOI_and_E08_robust_policy.png"); plt.close(fig)

# Fig 6
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.5))
dec_row = df_e06[df_e06["period"] == "2010 -> 2025"].iloc[0]
vals = [dec_row["quantum"], dec_row["tempo"], dec_row["structure"], dec_row["delta_B"]]
labels = ["生育数量\n(Quantum)", "生育时机\n(Tempo)", "育龄人口结构\n(Structure)", "总出生变动\n(2010→2025)"]
bar_c = ["#2874A6", "#D68910", "#922B21", "#2C3E50"]
ax1.bar(labels, vals, color=bar_c, alpha=0.88)
ax1.axhline(0, color="black", lw=0.8)
for i, v in enumerate(vals):
    ax1.text(i, v - 50 if v < 0 else v + 15, f"{v:+.0f}", ha="center", fontsize=8.5, fontweight="bold")
ax1.set_ylabel("出生人口变动贡献（万人）"); ax1.set_title("(a) E06 2010—2025年中国出生下滑的Shapley机制分解", fontsize=10, fontweight="bold")
ax1.grid(axis="y", alpha=0.25)

for scen, ls, col in [
    ("Low (Continued Decline)", "--", "#922B21"),
    ("Baseline (Current 2025 Pol)", "-", "#1F618D"),
    ("Robust Support (U2/U4)", "-.", "#1E8449"),
]:
    sub = df_proj[df_proj["scenario"] == scen]
    ax2.plot(sub["year"], sub["pop_wan"] / 1e4, ls=ls, color=col, lw=2.2, marker="o", ms=4, label=scen)
ax2.set_xlabel("年份"); ax2.set_ylabel("总人口（亿人）")
ax2.set_title("(b) E08 2025—2070年中国人口长期情景敏感性", fontsize=10, fontweight="bold")
ax2.legend(frameon=False, fontsize=8.5); ax2.grid(alpha=0.25)
fig.tight_layout(); fig.savefig(FIGS / "Fig06_E06_decomposition_and_E08_projections.png"); plt.close(fig)

print("=== UPDATED E05 ===")
print(df_e05.to_string(index=False))
print("\n=== UPDATED E09 ===")
print(df_e09.to_string(index=False))
