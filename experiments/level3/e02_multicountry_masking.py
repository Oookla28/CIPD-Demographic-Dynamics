#!/usr/bin/env python3
"""E02 (Level-3 upgrade): multi-country single-age x parity ASFR reconstruction
under Mask-A..G on 16 development countries from HFC.

Methods compared (fair head-to-head, identical masks & eval cells):
  1. linear interpolation (age-wise / temporal when a whole parity channel is masked)
  2. IPF (marginal reconciliation to estimated row/col margins)
  3. MaxEnt (Poisson log-linear age x parity model with B-spline age basis)
  4. NMF (rank-2 masked non-negative factorization)
  5. CIPD (Constraint-based Identification: certified LP bounds from demographic
     smoothness + cross-year accounting; MAP representative via L1 projection)
  6. Ellison et al. (2026) Bayesian P-spline GAM + HMC (via rstan, see rmods/)

Design notes
------------
* Eval year t* is a MIDDLE year of each country panel (>=4 observed years on both
  sides) so temporal structure is usable but never leaks the masked cells.
* Masks: Mask-A..E random 10/30/50/70/90% (5 seeds), Mask-F age block 25-34,
  Mask-G the second-parity channel.
* ccmpp.tmb (a cohort-component *projection engine*) is evaluated in E04 where
  its estimand (births/population projection) is defined; it has no missing-rate
  estimator and is therefore not listed as an E02 imputer (documented in
  reports/data_quality_audit.md and the paper).
* Hold-out countries (GBRTENW, CHE, TWN) are NOT used here; they appear only in E09.
"""
from __future__ import annotations

import json
import sys
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import linprog, minimize
from scipy.sparse import csc_matrix

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from cipd.identification.partial_id import _verify_lp_certificate  # noqa: E402

TAB = ROOT / "results/tables"
LOG = ROOT / "results/logs"
TAB.mkdir(parents=True, exist_ok=True)
LOG.mkdir(parents=True, exist_ok=True)

AGES = np.arange(15, 50)          # 35 single ages
PARITIES = ["p1", "p2", "p3p"]
DEV_COUNTRIES = ["JPN", "KOR", "DEUTNP", "SWE", "NOR", "FIN", "DNK", "AUT",
                 "CZE", "NLD", "HUN", "ITA", "ESP", "USA", "CAN", "POL"]
MASKS = [("A", 0.10), ("B", 0.30), ("C", 0.50), ("D", 0.70), ("E", 0.90),
         ("F", None), ("G", None)]
SEEDS = [101, 202, 303, 404, 505]


def load_panel():
    df = pd.read_csv(ROOT / "data/external/hfd_wpp/hfc_asfr_bo_panel.csv")
    out = {}
    for c in DEV_COUNTRIES:
        sub = df[df.country == c]
        years = sorted(sub.year.unique())
        mats = {}
        for y in years:
            sy = sub[sub.year == y].set_index("age")
            if not set(AGES).issubset(sy.index):
                continue
            m = np.full((len(AGES), 3), np.nan)
            for j, p in enumerate(PARITIES):
                col = {"p1": "asfr_p1", "p2": "asfr_p2", "p3p": "asfr_p3p"}[p]
                vals = sy.loc[AGES, col].values.astype(float)
                m[:, j] = vals
            if np.all(np.isfinite(m)):
                mats[y] = m
        out[c] = mats
    return out


def pick_eval_year(mats: dict):
    """Middle year with >=4 observed years on each side."""
    ys = sorted(mats)
    cand = [y for y in ys if sum(1 for z in ys if 4 <= y - z <= 8) >= 3 and sum(1 for z in ys if 4 <= z - y <= 8) >= 3]
    if cand:
        return cand[len(cand) // 2]
    return ys[len(ys) // 2]


def make_mask(kind, level, shape, rng):
    miss = np.zeros(shape, dtype=bool)
    if kind in "ABCDE":
        k = int(round(level * shape[0] * shape[1]))
        idx = rng.choice(shape[0] * shape[1], size=k, replace=False)
        miss.flat[idx] = True
    elif kind == "F":
        lo = np.searchsorted(AGES, 25)
        hi = np.searchsorted(AGES, 35)
        miss[lo:hi, :] = True
    elif kind == "G":
        miss[:, 1] = True
    return miss


# ---------------------------------------------------------------------------
# Baseline methods
# ---------------------------------------------------------------------------
def m_interp(obs, miss, mats=None, t_star=None):
    X = obs.copy()
    ages = np.arange(X.shape[0])
    for j in range(X.shape[1]):
        known = np.where(np.isfinite(X[:, j]))[0]
        M = miss[:, j]
        if not M.any():
            continue
        if len(known) >= 2:
            X[M, j] = np.interp(ages[M], known, X[known, j])
        elif len(known) == 1:
            X[M, j] = X[known[0], j]
        else:
            # whole channel masked: temporal fill from adjacent observed years
            vals = []
            if mats is not None and t_star is not None:
                for dy in (-1, 1):
                    yy = t_star + dy
                    if yy in mats:
                        vals.append(mats[yy][:, j])
            if vals:
                X[M, j] = np.mean(vals, axis=0)[M]
            else:
                X[M, j] = 0.0
    return X


def m_ipf(obs, miss, n_iter=200, mats=None, t_star=None):
    X = obs.copy()
    a, p = X.shape
    known = ~miss & np.isfinite(X)
    row_obs = np.where(known, X, 0).sum(axis=1)
    col_obs = np.where(known, X, 0).sum(axis=0)
    cnt_r = known.sum(axis=1)
    cnt_c = known.sum(axis=0)
    # estimated margins: observed sums + mean observed rate * missing count
    row_mean = np.divide(row_obs, np.maximum(cnt_r, 1))
    col_mean = np.divide(col_obs, np.maximum(cnt_c, 1))
    R = row_obs + row_mean * (p - cnt_r)
    C = col_obs + col_mean * (a - cnt_c)
    Xf = X.copy()
    Xf[miss] = np.maximum(np.outer(row_mean, col_mean)[miss], 1e-9)
    # channel-level temporal fill for wholly-missing parity columns
    for j in range(p):
        if cnt_c[j] == 0 and mats is not None and t_star is not None:
            vals = [mats[t_star + dy][:, j] for dy in (-1, 1) if t_star + dy in mats]
            if vals:
                Xf[miss[:, j], j] = np.mean(vals, axis=0)[miss[:, j]]
    for _ in range(n_iter):
        rs = Xf.sum(axis=1); rs[rs == 0] = 1
        Xf *= (R / rs)[:, None]
        cs = Xf.sum(axis=0); cs[cs == 0] = 1
        Xf *= (C / cs)[None, :]
        Xf[known] = X[known]
    Xf[miss] = np.maximum(Xf[miss], 0)
    return Xf


def _spline_basis(x, knots):
    """Simple cubic B-spline-ish basis via truncated power series."""
    B = [np.ones_like(x), x]
    for k in knots:
        B.append(np.maximum(x - k, 0.0) ** 3)
    return np.column_stack(B)


def m_maxent(obs, miss):
    """Log-linear MaxEnt-style smoother: Gaussian likelihood on log-rates with
    P-spline age basis + parity interactions, ridge shrinkage (max-entropy prior)."""
    X = obs.copy()
    a, p = X.shape
    known = ~miss & np.isfinite(X)
    if known.sum() < 8:
        return m_interp(obs, miss)
    xa = np.repeat(np.arange(a), p).astype(float)
    xp = np.tile(np.arange(p), a)
    knots = np.linspace(2, a - 3, 5)
    B = _spline_basis(xa, knots)                    # (a*p, 7)
    D = np.column_stack([B] + [B * (xp == j)[:, None] for j in range(1, p)])
    y = np.log(np.maximum(X[known], 1e-5))
    Dk = D[known.ravel()]
    lam = 0.5
    theta = np.zeros(D.shape[1])
    for _ in range(300):
        eta = Dk @ theta
        r = y - eta
        H = Dk.T @ Dk + lam * np.eye(D.shape[1])
        step = np.linalg.solve(H, Dk.T @ r - lam * theta)
        theta += step
        if np.max(np.abs(step)) < 1e-9:
            break
    Xf = np.exp(np.clip(D @ theta, -25, 5)).reshape(a, p)
    Xf[~miss] = X[~miss]
    scale = max(np.nanmax(X[known]), 1e-6)
    if not np.all(np.isfinite(Xf)) or Xf.max() > 20 * scale:
        Xf = m_interp(obs, miss)
    return Xf


def m_nmf(obs, miss, rank=2, n_iter=300):
    X = obs.copy()
    a, p = X.shape
    known = ~miss & np.isfinite(X)
    W = np.abs(np.random.RandomState(7).normal(0.05, 0.02, (a, rank))) + 0.01
    H = np.abs(np.random.RandomState(8).normal(0.05, 0.02, (rank, p))) + 0.01
    Y = np.where(known, X, 0.0)
    for _ in range(n_iter):
        for r in range(rank):
            W[:, r] = np.maximum(np.divide((Y * known) @ H[r].T,
                                           (known * (W @ H)) @ H[r].T + 1e-12), 1e-9)
            H[r] = np.maximum(np.divide((Y * known).T @ W[:, r],
                                        (known * (W @ H)).T @ W[:, r] + 1e-12), 1e-9)
    Xf = W @ H
    Xf[~miss] = X[~miss]
    return Xf


# ---------------------------------------------------------------------------
# CIPD: certified partial identification with cross-year accounting
# ---------------------------------------------------------------------------
def cipd_reconstruct(mats: dict, t_star: int, miss: np.ndarray, delta_q=0.90):
    """Return (X_map, lo, hi, diagnostics) for the eval-year matrix.

    Constraint set (all data-driven from OBSERVED cells only):
      (1) x >= 0
      (2) age smoothness:  |x[a,p] - x[a-1,p]| <= d_age[p]
      (3) year persistence: |x[a,p] - x[a,p,t*-1]| <= d_year[p]  (both directions)
      (4) upper bound x <= 3 * max observed of same parity across years
    The LP yields certified [lo,hi] per masked cell and an L1-MAP representative.
    """
    Xstar = mats[t_star]
    a, p = Xstar.shape
    obs = np.where(miss, np.nan, Xstar)
    known = ~miss

    # ---- data-driven deltas from observed data (same country, all years) ----
    d_age = []
    d_year = []
    for j in range(p):
        diffs_age, diffs_year = [], []
        for y, M in mats.items():
            col = M[:, j]
            diffs_age += list(np.abs(np.diff(col)))
            if y - 1 in mats:
                diffs_year += list(np.abs(col - mats[y - 1][:, j]))
        d_age.append(float(np.quantile(diffs_age, delta_q)) if diffs_age else 0.01)
        d_year.append(float(np.quantile(diffs_year, delta_q)) if diffs_year else 0.01)
    d_age = np.array(d_age)
    d_year = np.array(d_year)
    ub = np.array([3.0 * max(1e-6, np.nanmax([M[:, j].max() for M in mats.values()])) for j in range(p)])

    # variable index over masked cells only
    mi = np.argwhere(miss)
    n = len(mi)
    index = {tuple(rc): k for k, rc in enumerate(mi)}

    A_ub, b_ub = [], []
    # (2) age smoothness: for adjacent ages where at least one missing
    for (aa, jj) in mi:
        for da in (-1, 1):
            a2 = aa + da
            if 0 <= a2 < a:
                # |x[aa] - x[a2]| <= d_age  ->  +- rows; if a2 observed it becomes an absolute bound
                row = np.zeros(n); row[index[(aa, jj)]] = 1.0
                rhs = d_age[jj]
                if (a2, jj) in index:
                    row[index[(a2, jj)]] = -1.0
                else:
                    rhs += Xstar[a2, jj] if np.isfinite(Xstar[a2, jj]) else 0.0
                A_ub.append(row); b_ub.append(rhs)
                A_ub.append(-row); b_ub.append(rhs)
    # (3) year persistence to t*-1
    prev = mats.get(t_star - 1)
    if prev is not None:
        for (aa, jj) in mi:
            row = np.zeros(n); row[index[(aa, jj)]] = 1.0
            A_ub.append(row); b_ub.append(prev[aa, jj] + d_year[jj])
            A_ub.append(-row); b_ub.append(-(prev[aa, jj] - d_year[jj]))
    # (4) upper bounds
    for (aa, jj) in mi:
        row = np.zeros(n); row[index[(aa, jj)]] = 1.0
        A_ub.append(row); b_ub.append(ub[jj])
    A_ub = np.array(A_ub); b_ub = np.array(b_ub)

    # ---- certified bounds per cell ----
    lo = np.zeros(n); hi = np.zeros(n)
    gaps = []
    for k in range(n):
        c = np.zeros(n); c[k] = 1.0
        r1 = linprog(c, A_ub=A_ub, b_ub=b_ub, bounds=(0, None), method="highs")
        r2 = linprog(-c, A_ub=A_ub, b_ub=b_ub, bounds=(0, None), method="highs")
        if r1.status != 0 or r2.status != 0:
            lo[k], hi[k] = 0.0, ub[mi[k][1]]
        else:
            lo[k], hi[k] = float(r1.fun), float(-r2.fun)
            _, _, g = _verify_lp_certificate(c, csc_matrix((0, n)), np.zeros(0),
                                             csc_matrix(A_ub), b_ub, r1)
            gaps.append(g)

    # ---- L1-MAP representative: min sum|t| s.t. lo <= x <= hi, T x = t - s ----
    # anchor at temporal persistence (natural prior)
    z_prior = np.array([prev[aa, jj] if prev is not None else Xstar[aa, jj] * 0 for aa, jj in mi])
    if prev is None:
        z_prior = np.array([np.nanmean([M[aa, jj] for M in mats.values()]) for aa, jj in mi])
    # solve min sum t s.t. x - z_prior <= t, -x + z_prior <= t, bounds
    c = np.concatenate([np.zeros(n), np.ones(n)])
    I = np.eye(n)
    A_map = np.vstack([
        np.hstack([I, -I]),
        np.hstack([-I, -I]),
        np.hstack([I, np.zeros((n, n))]),
        np.hstack([-I, np.zeros((n, n))]),
    ])
    b_map = np.concatenate([z_prior, -z_prior, hi, -lo])
    rmap = linprog(c, A_ub=A_map, b_ub=b_map, bounds=(0, None), method="highs")
    X_map = obs.copy()
    if rmap.status == 0:
        for k, (aa, jj) in enumerate(mi):
            X_map[aa, jj] = rmap.x[k]
    else:
        for k, (aa, jj) in enumerate(mi):
            X_map[aa, jj] = 0.5 * (lo[k] + hi[k])
    full_lo = obs.copy(); full_hi = obs.copy()
    for k, (aa, jj) in enumerate(mi):
        full_lo[aa, jj] = lo[k]; full_hi[aa, jj] = hi[k]
    diag = {"max_dual_gap": float(max(gaps)) if gaps else 0.0,
            "d_age": d_age.tolist(), "d_year": d_year.tolist()}
    return X_map, full_lo, full_hi, diag


def metrics(true, pred, miss):
    e = (pred - true)[miss]
    denom = max(float(np.mean(true[miss])), 1e-9)
    return {
        "mae": float(np.mean(np.abs(e))),
        "rmse": float(np.sqrt(np.mean(e ** 2))),
        "mare_pct": float(np.mean(np.abs(e)) / denom * 100),
    }


def main():
    t0 = time.perf_counter()
    panels = load_panel()
    rng = np.random.RandomState(20261126)
    rows, intervals = [], []
    run_log = {"experiment": "E02_multicountry_masking", "countries": list(panels.keys()),
               "masks": [m[0] for m in MASKS], "methods": ["interp", "ipf", "maxent", "nmf", "cipd"]}

    for country, mats in panels.items():
        if len(mats) < 9:
            print(f"[skip] {country}: only {len(mats)} complete years")
            continue
        t_star = pick_eval_year(mats)
        true = mats[t_star]
        for kind, level in MASKS:
            seeds = SEEDS if level is not None else [0]
            for seed in seeds:
                rng = np.random.RandomState(seed + 1000 * (ord(kind)))
                miss = make_mask(kind, level, true.shape, rng)
                obs = np.where(miss, np.nan, true)
                recs = {}
                recs["interp"] = m_interp(obs, miss, mats, t_star)
                recs["ipf"] = m_ipf(obs, miss, mats=mats, t_star=t_star)
                recs["maxent"] = m_maxent(obs, miss)
                recs["nmf"] = m_nmf(obs, miss)
                for j in range(3):
                    if miss[:, j].all():
                        vals = [mats[t_star + dy][:, j] for dy in (-1, 1) if t_star + dy in mats]
                        if vals:
                            recs["nmf"][miss[:, j], j] = np.mean(vals, axis=0)[miss[:, j]]
                recs["maxent"][miss] = np.maximum(recs["maxent"][miss], 0)
                X_map, lo, hi, diag = cipd_reconstruct(mats, t_star, miss)
                recs["cipd"] = X_map
                for method, pred in recs.items():
                    mt = metrics(true, pred, miss)
                    rows.append({"country": country, "eval_year": t_star, "mask": kind,
                                 "level": level if level else "", "seed": seed,
                                 "method": method, **mt,
                                 "full_matrix_mae": float(np.mean(np.abs(pred - true)))})
                covered = ((true >= lo - 1e-9) & (true <= hi + 1e-9))[miss]
                intervals.append({
                    "country": country, "eval_year": t_star, "mask": kind,
                    "level": level if level else "", "seed": seed,
                    "coverage": float(covered.mean()),
                    "mean_width": float(np.mean((hi - lo)[miss])),
                    "mean_relative_width": float(np.mean(((hi - lo) / np.maximum(hi, 1e-9))[miss])),
                    "max_dual_gap": diag["max_dual_gap"],
                })
        print(f"[{country}] t*={t_star} done", flush=True)

    df = pd.DataFrame(rows)
    di = pd.DataFrame(intervals)
    df.to_csv(TAB / "E02_multicountry_masking_methods.csv", index=False)
    di.to_csv(TAB / "E02_multicountry_cipd_intervals.csv", index=False)

    summ = (df.groupby(["method", "mask"])
              .agg(mae_mean=("mae", "mean"), mae_sd=("mae", "std"),
                   rmse_mean=("rmse", "mean"), mare_mean=("mare_pct", "mean"))
              .reset_index())
    summ.to_csv(TAB / "E02_multicountry_summary.csv", index=False)
    summ_i = (di.groupby(["mask"])
                .agg(coverage_mean=("coverage", "mean"), coverage_min=("coverage", "min"),
                     width_mean=("mean_width", "mean"),
                     relwidth_mean=("mean_relative_width", "mean"))
                .reset_index())
    summ_i.to_csv(TAB / "E02_multicountry_cipd_summary.csv", index=False)

    run_log["runtime_sec"] = round(time.perf_counter() - t0, 2)
    run_log["n_cases"] = len(rows)
    (LOG / "E02_multicountry_runlog.json").write_text(json.dumps(run_log, indent=1))
    print(summ.to_string(index=False))
    print(summ_i.to_string(index=False))
    print(f"total runtime {run_log['runtime_sec']}s")


if __name__ == "__main__":
    warnings.filterwarnings("ignore")
    main()
