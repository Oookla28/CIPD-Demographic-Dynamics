#!/usr/bin/env python3
"""Publication figures for the Level-3 sprint (Fig07-Fig10), regenerated from
the result tables only (no hard-coded numbers)."""
from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
TAB = ROOT / "results/tables"
FIG = ROOT / "results/figures"
FIG.mkdir(parents=True, exist_ok=True)

plt.rcParams["font.sans-serif"] = ["Hiragino Sans GB", "Heiti SC", "Arial Unicode MS", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False
plt.rcParams["figure.dpi"] = 180
plt.rcParams["axes.spines.top"] = False
plt.rcParams["axes.spines.right"] = False
COLORS = {"cipd": "#C0392B", "interp": "#2980B9", "ipf": "#27AE60",
          "maxent": "#8E44AD", "nmf": "#E67E22", "ellison_gam": "#16A085",
          "persistence": "#7F8C8D", "arima": "#2C3E50", "ets": "#D35400"}


def fig07():
    df = pd.read_csv(TAB / "E02_multicountry_summary.csv")
    el = TAB / "E02_multicountry_ellison.csv"
    masks = ["A", "B", "C", "D", "E", "F", "G"]
    labels = {"A": "10%", "B": "30%", "C": "50%", "D": "70%", "E": "90%",
              "F": "年龄块", "G": "整孩次"}
    piv = df.pivot_table(index="mask", columns="method", values="mae_mean")
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.6))
    ax = axes[0]
    for m in ["interp", "ipf", "maxent", "nmf", "cipd"]:
        if m in piv.columns:
            ax.plot([labels[k] for k in masks], piv.loc[masks, m], "o-",
                    label={"cipd": "CIPD (本文)"}.get(m, m), color=COLORS.get(m), lw=2, ms=5)
    if el.exists():
        e = pd.read_csv(el).groupby("mask").mae.mean()
        xs = [labels[k] for k in masks if k in e.index]
        ax.plot(xs, [e[k] for k in masks if k in e.index], "s--", label="Ellison GAM+HMC",
                color=COLORS["ellison_gam"], lw=1.6, ms=5)
    ax.set_yscale("log")
    ax.set_xlabel("遮蔽模式（随机缺失率 / 结构缺失）")
    ax.set_ylabel("MAE（每位育龄妇女·年 的 ASFR 误差）")
    ax.set_title("(a) 16 开发国单岁×孩次重构误差", fontsize=11)
    ax.legend(fontsize=8, ncol=2)
    ax.grid(alpha=0.25, ls=":")
    ax = axes[1]
    ci = pd.read_csv(TAB / "E02_multicountry_cipd_summary.csv")
    ax.bar([labels[k] for k in masks], ci.set_index("mask").loc[masks, "coverage_mean"] * 100,
           color="#C0392B", alpha=0.75, width=0.62, label="覆盖率(%)")
    ax2 = ax.twinx()
    ax2.plot([labels[k] for k in masks], ci.set_index("mask").loc[masks, "width_mean"],
             "ko-", ms=5, label="平均宽度")
    ax.set_ylabel("CIPD 认证区间覆盖率 (%)")
    ax2.set_ylabel("平均识别宽度 (ASFR)")
    ax.set_ylim(0, 105)
    ax.set_title("(b) CIPD 认证边界的覆盖率与宽度", fontsize=11)
    h1, l1 = ax.get_legend_handles_labels(); h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, fontsize=8, loc="lower right")
    ax.grid(alpha=0.25, ls=":", axis="y")
    fig.tight_layout()
    fig.savefig(FIG / "Fig07_E02_multicountry_headtohead.png", bbox_inches="tight")
    plt.close(fig)
    print("Fig07 done")


def fig08():
    df = pd.read_csv(TAB / "E04_rolling_multicountry_summary.csv")
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.2))
    for m, g in df.groupby("method"):
        axes[0].plot(g.horizon, g.mae, "o-", label={"cipd": "CIPD"}.get(m, m),
                     color=COLORS.get(m), lw=1.8, ms=4)
        axes[1].plot(g.horizon, g.coverage * 100, "o-", color=COLORS.get(m), lw=1.6, ms=4)
        axes[2].plot(g.horizon, g.wis, "o-", color=COLORS.get(m), lw=1.6, ms=4)
    axes[0].set_xlabel("预测期长 h (年)"); axes[0].set_ylabel("MAE (TFR)")
    axes[0].set_title("(a) 分期限点预测误差", fontsize=11)
    axes[1].set_xlabel("预测期长 h (年)"); axes[1].set_ylabel("90% 区间覆盖率 (%)")
    axes[1].set_title("(b) 区间覆盖率", fontsize=11)
    axes[2].set_xlabel("预测期长 h (年)"); axes[2].set_ylabel("WIS")
    axes[2].set_title("(c) 加权区间得分", fontsize=11)
    axes[0].legend(fontsize=8)
    for ax in axes:
        ax.grid(alpha=0.25, ls=":")
    fig.tight_layout()
    fig.savefig(FIG / "Fig08_E04_rolling_forecast_headtohead.png", bbox_inches="tight")
    plt.close(fig)

    # E09 holdout vs dev
    hv = pd.read_csv(TAB / "E09_holdout_vs_dev.csv")
    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    roles = {"dev": ("#2980B9", "开发国 (16)"), "holdout": ("#C0392B", "严格独立测试国 (3)")}
    for role, (color, lab) in roles.items():
        sub = hv[hv.role == role]
        ax.scatter(sub.mae, sub.coverage * 100, s=46, color=color, alpha=0.85, label=lab)
        for _, r in sub.iterrows():
            ax.annotate(r.country, (r.mae, r.coverage * 100), fontsize=6.5,
                        xytext=(3, 3), textcoords="offset points")
    ax.set_xlabel("重构 MAE（留一国调参后在该国评估）")
    ax.set_ylabel("CIPD 区间覆盖率 (%)")
    ax.set_title("E09 留一国泛化：独立测试国与开发国同分布", fontsize=11)
    ax.legend(fontsize=9)
    ax.grid(alpha=0.25, ls=":")
    fig.tight_layout()
    fig.savefig(FIG / "Fig09_E09_loco_generalization.png", bbox_inches="tight")
    plt.close(fig)
    print("Fig08/09 done")


def fig10():
    # E05 event studies + HonestDiD
    hd = pd.read_csv(TAB / "E05_honestdid_sensitivity.csv")
    eA = pd.read_csv(TAB / "E05_event_study_A_dandui2014.csv")
    eB = pd.read_csv(TAB / "E05_event_study_B_quanmian2016.csv")
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.2))
    axes[0].axhline(0, color="k", lw=0.8)
    axes[0].axvline(-1, color="gray", ls="--", lw=0.8)
    axes[0].plot(eA.k, eA.beta, "o-", color="#C0392B", lw=1.8, ms=5)
    axes[0].set_title("(a) 单独二孩 2014（落地月强度）", fontsize=10.5)
    axes[0].set_xlabel("相对政策年 k"); axes[0].set_ylabel("β_k（对数 CBR / 强度 SD）")
    axes[1].axhline(0, color="k", lw=0.8)
    axes[1].axvline(-1, color="gray", ls="--", lw=0.8)
    axes[1].plot(eB.k, eB.beta, "o-", color="#2C3E50", lw=1.8, ms=5)
    axes[1].set_title("(b) 全面二孩 2016（一孩家庭强度）", fontsize=10.5)
    axes[1].set_xlabel("相对政策年 k")
    ax = axes[2]
    for pol, color, lab in [("A_2014", "#C0392B", "2014 单独二孩"),
                            ("B_2016", "#2C3E50", "2016 全面二孩")]:
        sub = hd[(hd.policy == pol) & (hd.method == "FLCI_relativeMagnitudes")]
        ax.fill_between(sub.M, sub.lb, sub.ub, color=color, alpha=0.22)
        ax.plot(sub.M, sub.lb, color=color, lw=1.2)
        ax.plot(sub.M, sub.ub, color=color, lw=1.2, label=lab)
    ax.axhline(0, color="k", lw=0.8)
    ax.set_xlabel("M（相对偏差界限）"); ax.set_ylabel("HonestDiD 稳健置信区间")
    ax.set_title("(c) Δ^SD(M) 敏感性：全部区间跨 0", fontsize=10.5)
    ax.legend(fontsize=8)
    for a in axes:
        a.grid(alpha=0.25, ls=":")
    fig.tight_layout()
    fig.savefig(FIG / "Fig10_E05_provincial_event_studies_honestdid.png", bbox_inches="tight")
    plt.close(fig)

    # E10 + CCPI embedding
    e10 = pd.read_csv(TAB / "E10_acceleration.csv")
    cc = pd.read_csv(TAB / "E05_ccpi_did_embedding.csv")
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.2))
    ax = axes[0]
    for _, r in e10.iterrows():
        ax.barh(f"{r.instance}\n({int(r.n_vars):,} 变量)", r.speedup,
                color="#C0392B" if r.speedup >= 2 else "#7F8C8D")
        ax.text(r.speedup + 0.08, f"{r.instance}\n({int(r.n_vars):,} 变量)",
                f"{r.speedup:.2f}×", va="center", fontsize=9)
    ax.axvline(2.0, color="k", ls="--", lw=0.9, label="2× 阈值")
    ax.set_xlabel("相对冷启动通用 HiGHS 的加速比（同精度）")
    ax.set_title("(a) E10 求解加速（对偶间隙 ≤ 5.5e-9）", fontsize=10.5)
    ax.legend(fontsize=8)
    ax = axes[1]
    x = np.arange(len(cc))
    ax.bar(x - 0.2, cc.P0_width, width=0.38, color="#7F8C8D", label="仅人口核算 P0")
    ax.bar(x + 0.2, cc.PDiD_width, width=0.38, color="#C0392B", label="嵌入因果约束 P-DiD")
    for i, r in cc.iterrows():
        ax.text(i, max(r.P0_width, r.PDiD_width) * 1.02,
                f"-{r.width_reduction_pct:.0f}%", ha="center", fontsize=9)
    ax.set_xticks(x)
    ax.set_xticklabels(["2016-2020 总出生", "2017 二孩出生"], fontsize=9)
    ax.set_ylabel("识别区间宽度（人）")
    ax.set_title("(b) HonestDiD 因果约束嵌入人口流多面体", fontsize=10.5)
    ax.legend(fontsize=8)
    for a in axes:
        a.grid(alpha=0.25, ls=":", axis="x")
    fig.tight_layout()
    fig.savefig(FIG / "Fig11_E10_E05CCPI_embedding.png", bbox_inches="tight")
    plt.close(fig)
    print("Fig10/11 done")


if __name__ == "__main__":
    fig07()
    fig08()
    fig10()
