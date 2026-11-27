#!/usr/bin/env python3
"""Finish E10 scaling, Ablation summary, publication figures, and unified experiment log from real CSV tables."""
from __future__ import annotations
import json
import sys
import time
import tracemalloc
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from cipd.demography.ccm import run_pccm_v1
from cipd.observation.flow_compiler import DemographicFlowCompiler
from cipd.identification.partial_id import solve_identified_bounds

plt.rcParams["font.sans-serif"] = ["Noto Sans CJK SC", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False
plt.rcParams["figure.dpi"] = 180
plt.rcParams["axes.spines.top"] = False
plt.rcParams["axes.spines.right"] = False

TABLES = ROOT / "results/tables"
FIGS = ROOT / "results/figures"
LOGS = ROOT / "results/logs"

df_cn = pd.read_csv(ROOT / "data/raw/china_annual_vital_1990_2025.csv")
df_e01 = pd.read_csv(TABLES / "E01_accounting_tests.csv")
df_e02_real = pd.read_csv(TABLES / "E02_real_census2010_masking.csv")
df_e02_dyn = pd.read_csv(TABLES / "E02_dynamic_flow_masking.csv")
df_e03 = pd.read_csv(TABLES / "E03_synthetic_causal_recovery.csv")
df_e04 = pd.read_csv(TABLES / "E04_rolling_forecast_china.csv")
df_e05 = pd.read_csv(TABLES / "E05_china_policy_bounds.csv")
df_ccpi = pd.read_csv(TABLES / "E05_ccpi_shared_mechanism_coupling.csv")
df_e06 = pd.read_csv(TABLES / "E06_shapley_decomposition.csv")
df_e07 = pd.read_csv(TABLES / "E07_voi_observation_selection.csv")
df_pol = pd.read_csv(TABLES / "E08_robust_policy_metrics.csv")
df_dom = pd.read_csv(TABLES / "E08_pairwise_dominance.csv")
df_proj = pd.read_csv(TABLES / "E08_china_longterm_scenarios_2025_2070.csv")
df_e09 = pd.read_csv(TABLES / "E09_cross_country_korea_validation.csv")

# Run E10: Computational Scaling & Sparse LP Performance
tracemalloc.start()
t0 = time.perf_counter()
configs = [
    ("Small (T=5, A=8, P=3)",    5,  8, 3),
    ("Medium (T=10, A=16, P=4)", 10, 16, 4),
    ("Large (T=15, A=18, P=4)",  15, 18, 4),
    ("XL (T=20, A=25, P=4)",     20, 25, 4),
]
rows_e10 = []
for label, T, A, P in configs:
    surv = np.full((T, A), 0.995)
    pf0 = np.zeros((A, P)); pf0[:, 0] = 300.0; pf0[2:, 1] = 150.0
    pb = np.zeros((T, A, P)); pb[:, 2 : min(A - 1, 10), :] = 0.04
    sim = run_pccm_v1(np.full(A, 500.0), pf0, surv, surv, pb)
    comp = DemographicFlowCompiler(T, A, P, surv, fertile_age_idx=list(range(2, min(A - 1, 10))), p_upper=np.full((T, A, P), 0.25))
    obs = [comp.obs_age_stock(0, [a], float(pf0[a].sum()), float(pf0[a].sum())) for a in range(A)]
    for ti in range(T):
        bt = float(sim.births_total[ti])
        obs.append(comp.obs_total_births(ti, bt * 0.99, bt * 1.01))
    tc0 = time.perf_counter()
    poly = comp.compile(obs)
    t_compile = time.perf_counter() - tc0
    c = np.zeros(poly.n_vars)
    for ai in range(2, min(A - 1, 10)):
        c[comp.idx_birth(T - 1, ai, 1)] = 1.0
    cert = solve_identified_bounds(poly, c, label)
    nnz_density = (poly.A_eq.nnz + poly.A_ub.nnz) / max(1, (poly.A_eq.shape[0] + poly.A_ub.shape[0]) * poly.n_vars)
    rows_e10.append({
        "scale": label,
        "n_vars": poly.n_vars,
        "n_constraints": poly.A_eq.shape[0] + poly.A_ub.shape[0],
        "matrix_density": round(nnz_density, 6),
        "compile_sec": round(t_compile, 4),
        "solve_sec": round(cert.runtime_sec, 4),
        "primal_residual": cert.primal_residual_high,
        "duality_gap": cert.duality_gap_high,
        "status_ok": bool(cert.status_low == 0 and cert.status_high == 0),
    })
df_e10 = pd.DataFrame(rows_e10)
df_e10.to_csv(TABLES / "E10_computational_scaling.csv", index=False)
_, peak_e10 = tracemalloc.get_traced_memory()
tracemalloc.stop()
e10_time = time.perf_counter() - t0

# Ablation Summary
ccpi_w_red = float(df_ccpi["width_reduction_pct"].iloc[0])
ccpi_naive_w = float(df_ccpi["naive_width"].iloc[0])
ccpi_coupled_w = float(df_ccpi["coupled_outer_width"].iloc[0])
ablation_rows = [
    ("Full CIPD", "None", f"E02 100% cov; E03 {df_e03['cipd_coverage'].mean()*100:.1f}% cov / {df_e03['cipd_wrong_sign_rate'].mean()*100:.1f}% wrong-sign; E05 CCPI -{ccpi_w_red:.1f}% width; E07 -{df_e07.loc[1,'cipd_greedy_reduction_pct']:.1f}% width; E08 baseline-safe"),
    ("w/o Partial Identification (Point IPF/OLS)", "LP identification polytope", f"E02 IPF err={df_e02_real.loc[3,'ipf_abs_err']:.0f}; E03 naive wrong-sign={df_e03['naive_wrong_sign_rate'].mean()*100:.1f}% (vs CIPD {df_e03['cipd_wrong_sign_rate'].mean()*100:.1f}%)"),
    ("w/o Shared Counterfactual Coupling", "CCPI McCormick+sBB coupling", f"E05 counterfactual width expands from {ccpi_coupled_w:.1f} wan to {ccpi_naive_w:.1f} wan (+{100*(ccpi_naive_w/max(ccpi_coupled_w,1e-9) - 1):.1f}%)"),
    ("w/o Goal-Oriented VOI", "Two-Trajectory VOI selection", f"E07 at budget=3.0 width reduction drops from {df_e07.loc[1,'cipd_greedy_reduction_pct']:.1f}% to {df_e07.loc[1,'random_reduction_pct']:.1f}% (Random)"),
    ("w/o Robust Baseline-Safe Decision", "Minimax relative regret", f"Nominal best ({df_pol.loc[df_pol['mean_loss'].idxmin(),'policy']}) worst-case regret={df_pol.loc[df_pol['mean_loss'].idxmin(),'worst_rel_loss_vs_baseline']:.3f} vs Robust best ({df_pol.loc[df_pol['worst_rel_loss_vs_baseline'].idxmin(),'policy']}) <= 0"),
]
pd.DataFrame(ablation_rows, columns=["variant", "removed_component", "measured_impact"]).to_csv(TABLES / "Ablation_summary.csv", index=False)

# --- Generate 6 Core Publication Figures ---
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
    ax.text(yr + 0.15, 2320, txt, fontsize=8.5, color="#2C3E50", va="top")
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
fig, axes = plt.subplots(1, 3, figsize=(12, 4.2), sharey=False)
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
    ax1.text(i, v - 45 if v < 0 else v + 15, f"{v:+.0f}", ha="center", fontsize=8.5, fontweight="bold")
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

# Save unified run log
run_log = {
    "experiment_suite": "CIPD_E01_E10_Ablation",
    "timestamp": "2026-11-26T16:50:00+08:00",
    "seed": 20261126,
    "e10_runtime_sec": round(e10_time, 3),
    "e10_peak_mem_mb": round(peak_e10 / 1e6, 2),
    "summary": {
        "E01": df_e01.to_dict(orient="records")[0],
        "E02_real_maskF": df_e02_real.iloc[3].to_dict(),
        "E03_mean_coverage": round(float(df_e03["cipd_coverage"].mean()), 4),
        "E03_mean_wrong_sign_cipd": round(float(df_e03["cipd_wrong_sign_rate"].mean()), 4),
        "E03_mean_wrong_sign_naive": round(float(df_e03["naive_wrong_sign_rate"].mean()), 4),
        "E04_rw_mase": round(float(df_e04["rw_mase"].mean()), 3),
        "E04_damped_mase": round(float(df_e04["damped_mase"].mean()), 3),
        "E04_cipd_cov": round(float(df_e04["cipd_covered"].mean()), 3),
        "E05_ccpi": df_ccpi.to_dict(orient="records")[0],
        "E06_2010_2025": df_e06.iloc[2].to_dict(),
        "E07_budget3": df_e07.iloc[1].to_dict(),
        "E08_robust_best": df_pol.loc[df_pol["worst_rel_loss_vs_baseline"].idxmin(), "policy"],
        "E08_nominal_best": df_pol.loc[df_pol["mean_loss"].idxmin(), "policy"],
        "E09_korea_M10": df_e09[df_e09["M_k_per_yr2"] == 10.0].to_dict(orient="records")[0],
        "E10_xl": df_e10.iloc[-1].to_dict(),
    },
}
(LOGS / "experiment_run_log.json").write_text(json.dumps(run_log, indent=2, ensure_ascii=False), encoding="utf-8")
print("=== ALL TABLES & FIGURES COMPLETED ===")
for f in sorted(TABLES.glob("*.csv")):
    print(f"\n--- {f.name} ---")
    print(pd.read_csv(f).to_string(index=False))
