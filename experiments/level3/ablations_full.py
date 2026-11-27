#!/usr/bin/env python3
"""Full 9-item ablation suite (Level-3 upgrade).

A1-A5 recompute the multi-country reconstruction on HFC development countries
with one structural component removed at a time; A6-A9 aggregate the previously
measured counterfactual / VOI / robust-decision ablations from the Level-2
pipeline tables into a single dossier.

Variants (mask C = 50% random, seed 101, all 16 dev countries):
  Full CIPD                 all constraints
  w/o Population Accounting drop the cross-year persistence constraint
  w/o Cohort Dynamics       drop the age-smoothness constraint
  w/o Parity Structure      reconstruct the aggregate age profile only and split
                            by observed parity shares (parity-free model)
  w/o Hierarchical Calibration  pooled (global) smoothness deltas instead of
                            country-specific quantiles
  w/o Partial Identification    point MAP only (no certified bounds; the
                            "coverage" degenerates to 0 by construction)
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import linprog

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "experiments/level3"))
from e02_multicountry_masking import (load_panel, pick_eval_year, make_mask,
                                      MASKS, AGES, m_interp, metrics)

TAB = ROOT / "results/tables"
LOG = ROOT / "results/logs"


def cipd_variant(mats, t_star, miss, variant="full", pooled_deltas=None, q=0.90):
    Xstar = mats[t_star]
    a, p = Xstar.shape
    obs = np.where(miss, np.nan, Xstar)
    known = ~miss
    # deltas
    if pooled_deltas is not None:
        d_age, d_year = pooled_deltas
    else:
        d_age, d_year = [], []
        for j in range(p):
            da, dy = [], []
            for y, M in mats.items():
                da += list(np.abs(np.diff(M[:, j])))
                if y - 1 in mats:
                    dy += list(np.abs(M[:, j] - mats[y - 1][:, j]))
            d_age.append(float(np.quantile(da, q)) if da else 0.01)
            d_year.append(float(np.quantile(dy, q)) if dy else 0.01)
        d_age, d_year = np.array(d_age), np.array(d_year)
    ub = np.array([3.0 * max(1e-6, np.nanmax([M[:, j].max() for M in mats.values()]))
                   for j in range(p)])
    if variant == "no_parity":
        # aggregate model: one column = total; split by observed parity shares
        tot_obs = np.where(np.all(miss, axis=1), np.nan, Xstar.sum(axis=1))
        miss_tot = np.all(miss, axis=1)
        # age interpolation of the total, then parity shares from observed cells
        tot_f = tot_obs.copy()
        idx = np.where(~miss_tot)[0]
        if len(idx) >= 2:
            tot_f[miss_tot] = np.interp(np.arange(a)[miss_tot], idx, tot_obs[idx])
        shares = np.tile(np.array([0.45, 0.42, 0.13]), (a, 1))
        for aa in range(a):
            ok = ~miss[aa]
            if ok.sum() >= 2:
                sh = np.zeros(p)
                sh[ok] = Xstar[aa][ok] / max(Xstar[aa][ok].sum(), 1e-9)
                sh[~ok] = sh[ok].mean()
                shares[aa] = sh
            elif ok.sum() == 1:
                shares[aa] = Xstar[aa] / max(Xstar[aa].sum(), 1e-9)
        X = np.where(miss, tot_f[:, None] * shares, Xstar)
        return X, obs.copy(), obs.copy(), {"variant": variant}
    mi = np.argwhere(miss)
    n = len(mi)
    index = {tuple(rc): k for k, rc in enumerate(mi)}
    A_ub, b_ub = [], []
    if variant != "no_cohort":  # age smoothness
        for (aa, jj) in mi:
            for da_ in (-1, 1):
                a2 = aa + da_
                if 0 <= a2 < a:
                    row = np.zeros(n); row[index[(aa, jj)]] = 1.0
                    rhs = d_age[jj]
                    if (a2, jj) in index:
                        row[index[(a2, jj)]] = -1.0
                    else:
                        rhs += Xstar[a2, jj]
                    A_ub.append(row); b_ub.append(rhs)
                    A_ub.append(-row); b_ub.append(rhs)
    prev = mats.get(t_star - 1)
    if variant != "no_accounting" and prev is not None:  # cross-year accounting
        for (aa, jj) in mi:
            row = np.zeros(n); row[index[(aa, jj)]] = 1.0
            A_ub.append(row); b_ub.append(prev[aa, jj] + d_year[jj])
            A_ub.append(-row); b_ub.append(-(prev[aa, jj] - d_year[jj]))
    for (aa, jj) in mi:
        row = np.zeros(n); row[index[(aa, jj)]] = 1.0
        A_ub.append(row); b_ub.append(ub[jj])
    A_ub = np.array(A_ub); b_ub = np.array(b_ub)
    lo = np.zeros(n); hi = np.zeros(n)
    if variant == "no_partial_id":
        # MAP only, no certified bounds: point estimate = midpoint of the box
        for k, (aa, jj) in enumerate(mi):
            lo[k] = max(0.0, Xstar[aa, jj] * 0.0)
            hi[k] = ub[jj]
    else:
        for k in range(n):
            c = np.zeros(n); c[k] = 1.0
            r1 = linprog(c, A_ub=A_ub, b_ub=b_ub, bounds=(0, None), method="highs")
            r2 = linprog(-c, A_ub=A_ub, b_ub=b_ub, bounds=(0, None), method="highs")
            lo[k] = float(r1.fun) if r1.status == 0 else 0.0
            hi[k] = float(-r2.fun) if r2.status == 0 else ub[mi[k][1]]
    # MAP = L1 projection toward temporal persistence
    z_prior = np.array([prev[aa, jj] if prev is not None else 0.0 for aa, jj in mi])
    c = np.concatenate([np.zeros(n), np.ones(n)])
    I = np.eye(n)
    A_map = np.vstack([np.hstack([I, -I]), np.hstack([-I, -I]),
                       np.hstack([I, np.zeros((n, n))]), np.hstack([-I, np.zeros((n, n))])])
    b_map = np.concatenate([z_prior, -z_prior, hi, -lo])
    rmap = linprog(c, A_ub=A_map, b_ub=b_map, bounds=(0, None), method="highs")
    X = obs.copy()
    if rmap.status == 0:
        for k, (aa, jj) in enumerate(mi):
            X[aa, jj] = rmap.x[k]
    else:
        for k, (aa, jj) in enumerate(mi):
            X[aa, jj] = 0.5 * (lo[k] + hi[k])
    full_lo = obs.copy(); full_hi = obs.copy()
    for k, (aa, jj) in enumerate(mi):
        full_lo[aa, jj] = lo[k]; full_hi[aa, jj] = hi[k]
    return X, full_lo, full_hi, {"variant": variant}


def main():
    t0 = time.perf_counter()
    panels = load_panel()
    mask_idx = [m for m in MASKS if m[0] == "C"][0]
    kind, level = mask_idx
    # pooled deltas (global quantiles across all countries) for A5
    da_all, dy_all = [], []
    for c, mats in panels.items():
        for j in range(3):
            for y, M in mats.items():
                da_all += list(np.abs(np.diff(M[:, j])))
                if y - 1 in mats:
                    dy_all += list(np.abs(M[:, j] - mats[y - 1][:, j]))
    pooled = (np.array([np.quantile(da_all, 0.9)] * 3), np.array([np.quantile(dy_all, 0.9)] * 3))

    rows = []
    for variant, label in [
        ("full", "Full CIPD"),
        ("no_accounting", "w/o Population Accounting"),
        ("no_cohort", "w/o Cohort Dynamics"),
        ("no_parity", "w/o Parity Structure"),
        ("pooled", "w/o Hierarchical Calibration"),
        ("no_partial_id", "w/o Partial Identification"),
    ]:
        maes, covs, widths = [], [], []
        for c, mats in panels.items():
            if len(mats) < 9:
                continue
            t_star = pick_eval_year(mats)
            true = mats[t_star]
            rng = np.random.RandomState(101 + 1000 * ord(kind))
            miss = make_mask(kind, level, true.shape, rng)
            kw = dict(pooled_deltas=pooled) if variant == "pooled" else {}
            X, lo, hi, _ = cipd_variant(mats, t_star, miss,
                                        variant=("full" if variant == "pooled" else variant), **kw)
            mt = metrics(true, X, miss)
            maes.append(mt["mae"])
            if variant == "no_partial_id":
                covs.append(0.0); widths.append(0.0)
            else:
                covs.append(float(((true >= lo - 1e-9) & (true <= hi + 1e-9))[miss].mean()))
                widths.append(float(np.mean((hi - lo)[miss])))
        rows.append({"variant": label, "mask": kind,
                     "mae_mean": float(np.mean(maes)),
                     "coverage_mean": float(np.mean(covs)),
                     "width_mean": float(np.mean(widths)),
                     "n_countries": len(maes)})
        print(f"[{label}] mae={rows[-1]['mae_mean']:.5f} cov={rows[-1]['coverage_mean']:.3f}")

    base = rows[0]["mae_mean"]
    # aggregate the Level-2 structural ablations
    rows += [
        {"variant": "w/o Shared Counterfactual Coupling", "mask": "",
         "mae_mean": None, "coverage_mean": None, "width_mean": None,
         "n_countries": None,
         "note": "E05 CCPI: width without coupling 908.8 wan vs with 409.6 wan (+121.9%); "
                 "causal-constrained P_DiD reduces 56.1% (total 2016-2020) / 88.0% (p2 2017)"},
        {"variant": "w/o Goal-Oriented VOI", "mask": "", "mae_mean": None,
         "coverage_mean": None, "width_mean": None, "n_countries": None,
         "note": "E07: width reduction at budget 3.0 falls 96.5% -> 10.4% (random selection)"},
        {"variant": "w/o Robust Baseline-Safe Decision", "mask": "", "mae_mean": None,
         "coverage_mean": None, "width_mean": None, "n_countries": None,
         "note": "E08: nominal-best worst-case regret -0.046 vs baseline-safe R(u) <= 0"},
    ]
    df = pd.DataFrame(rows)
    df["mae_increase_pct_vs_full"] = [
        (r["mae_mean"] / base - 1) * 100 if r["mae_mean"] is not None else None for r in rows]
    df.to_csv(TAB / "Ablation_summary.csv", index=False)
    (LOG / "ablation_runlog.json").write_text(json.dumps({
        "runtime_sec": round(time.perf_counter() - t0, 1),
        "mask": "C_50pct_seed101", "n_countries": rows[0]["n_countries"],
        "full_mae": base, "pooled_deltas": [float(pooled[0][0]), float(pooled[1][0])]}))
    print(df.round(5).to_string(index=False))


if __name__ == "__main__":
    main()
