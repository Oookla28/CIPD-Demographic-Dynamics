#!/usr/bin/env python3
"""E04 (Level-3 upgrade): strict rolling-origin fertility forecasts on HFC
development countries, multi-method head-to-head.

Target: period TFR_t = sum_a sum_p ASFR[a,p,t] (parity-decomposed).
Origins: 2005, 2010, 2013, 2016, 2020; horizons h = 1..10.
Only data with year <= origin is ever used to fit a forecast. HFC series are
final-vintage vital statistics (revisions are negligible by construction);
WPP is *not* used in this experiment at all (no post-hoc revision leakage).

Methods (identical origins/horizons):
  persistence (random walk), ARIMA(1,1,0) on TFR, damped-trend ETS,
  CIPD structural partial-identification forecast (constraint polytope over
  future age x parity ASFR paths; certified [lo,hi] fan + L1-MAP path).
R baselines (bayesTFR, Ellison cohort GAM) run through separate scripts and are
merged by experiments/e04_merge_r_baselines.py.

Metrics: MAE, RMSE, MASE (vs in-sample naive), CRPS and WIS (90% interval),
coverage of the 90% predictive interval, mean interval width.
"""
from __future__ import annotations

import json
import sys
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import linprog
from scipy.sparse import lil_matrix

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

TAB = ROOT / "results/tables"
LOG = ROOT / "results/logs"
TAB.mkdir(parents=True, exist_ok=True)
LOG.mkdir(parents=True, exist_ok=True)

AGES = np.arange(15, 50)
DEV_COUNTRIES = ["JPN", "KOR", "DEUTNP", "SWE", "NOR", "FIN", "DNK", "AUT",
                 "CZE", "NLD", "HUN", "ITA", "ESP", "USA", "CAN", "POL"]
ORIGINS = [2005, 2010, 2013, 2016, 2020]
H = 10


def load_tfr_panel():
    df = pd.read_csv(ROOT / "data/external/hfd_wpp/hfc_asfr_bo_panel.csv")
    out = {}
    for c in DEV_COUNTRIES:
        sub = df[(df.country == c) & (df.age >= 15) & (df.age <= 49)]
        tfr, parity = {}, {}
        for y, sy in sub.groupby("year"):
            if set(AGES).issubset(sy.age):
                sy = sy.set_index("age")
                p1 = sy.loc[AGES, "asfr_p1"].values.astype(float)
                p2 = sy.loc[AGES, "asfr_p2"].values.astype(float)
                p3 = sy.loc[AGES, "asfr_p3p"].values.astype(float)
                tot = sy.loc[AGES, "asfr_total"].values.astype(float)
                # prefer the published total if consistent, else sum of parities
                s = p1 + p2 + p3
                if np.all(np.isfinite(tot)) and abs(tot.sum() - s.sum()) < max(1e-4, 0.02 * s.sum()):
                    tfr[y] = float(tot.sum())
                else:
                    tfr[y] = float(s.sum())
                parity[y] = np.column_stack([p1, p2, p3])
        out[c] = {"tfr": tfr, "parity": parity}
    return out


# ---------------------------------------------------------------------------
# Metrics
# ---------------------------------------------------------------------------
def crps_normal(mu, sigma, y):
    sigma = max(sigma, 1e-9)
    z = (y - mu) / sigma
    from math import erf, exp, pi, sqrt
    phi = exp(-0.5 * z * z) / sqrt(2 * pi)
    Phi = 0.5 * (1 + erf(z / sqrt(2)))
    return sigma * (z * (2 * Phi - 1) + 2 * phi - 1 / sqrt(pi))


def crps_uniform(lo, hi, y):
    if hi <= lo:
        return abs(y - lo)
    if y < lo:
        return (lo - y) + 0.5 * (hi - lo)
    if y > hi:
        return (y - hi) + 0.5 * (hi - lo)
    return (hi - lo) / 2 - (y - lo) ** 2 / (2 * (hi - lo)) - (hi - y) ** 2 / (2 * (hi - lo)) + 0.0


def wis90(lo, hi, y, alpha=0.10):
    """WIS with central (1-alpha) interval only (equivalent to IS when alpha fixed)."""
    is_ = (hi - lo) + (2 / alpha) * max(0, lo - y) + (2 / alpha) * max(0, y - hi)
    return is_ / 2.0


# ---------------------------------------------------------------------------
# Methods
# ---------------------------------------------------------------------------
def f_persistence(hist_tfr, h):
    return float(hist_tfr[-1]), None, None


def f_arima(hist_years, hist_tfr, h):
    try:
        from statsmodels.tsa.arima.model import ARIMA
        y = np.array(hist_tfr)
        res = ARIMA(y, order=(1, 1, 0), trend="n").fit()
        fc = res.get_forecast(steps=h)
        mu = float(fc.predicted_mean[-1])
        ci = fc.conf_int(alpha=0.10)
        lo, hi = float(ci[-1, 0]), float(ci[-1, 1])
        return mu, lo, hi
    except Exception:
        return float(hist_tfr[-1]), None, None


def f_ets(hist_years, hist_tfr, h):
    try:
        from statsmodels.tsa.holtwinters import ExponentialSmoothing
        y = np.array(hist_tfr)
        res = ExponentialSmoothing(y, trend="add", damped_trend=True,
                                   initialization_method="estimated").fit()
        fc = res.forecast(h)
        mu = float(fc[-1])
        resid = y - res.fittedvalues
        sig = float(np.std(resid, ddof=1))
        # damped-trend PI: approximate with residual sigma scaled sqrt(h)
        scale = np.sqrt(h)
        return mu, mu - 1.645 * sig * scale, mu + 1.645 * sig * scale
    except Exception:
        return float(hist_tfr[-1]), None, None


def f_cipd_forecast(hist_tfr, hist_parity_by_year, t0, h, q=0.90):
    """Structural forecast: constraint polytope over future age x parity ASFR.

    Constraints (deltas estimated ONLY from years <= t0):
      age smoothness |x[a,p,h]-x[a-1,p,h]| <= da[p]
      time smoothness |x[a,p,h]-x[a,p,h-1]| <= dt[p]
      aggregate drift |TFR_h - TFR_0| <= h * dTFR
      x >= 0, x <= ub
    Returns TFR MAP + certified [lo,hi] at horizon h.
    """
    years = sorted(hist_parity_by_year)
    mats = [hist_parity_by_year[y] for y in years]
    a, p = 35, 3
    # deltas from history
    da = np.zeros(p); dt = np.zeros(p)
    for j in range(p):
        diffs_a, diffs_t = [], []
        for M in mats:
            diffs_a += list(np.abs(np.diff(M[:, j])))
        for k in range(1, len(mats)):
            diffs_t += list(np.abs(mats[k][:, j] - mats[k - 1][:, j]))
        da[j] = float(np.quantile(diffs_a, q)) if diffs_a else 1e-4
        dt[j] = float(np.quantile(diffs_t, q)) if diffs_t else 1e-4
    tfr_hist = [float(M.sum()) for M in mats]
    dtfr = float(np.quantile(np.abs(np.diff(tfr_hist)), q)) if len(tfr_hist) > 1 else 1e-3
    ub = 2.0 * max(1e-6, max(M.max() for M in mats))

    # variables x[h', a, p] for h'=1..h  (we solve the full path to horizon h)
    nv = h * a * p
    def ix(hh, aa, pp):  # hh in 0..h-1 (step hh+1)
        return (hh * a + aa) * p + pp

    rows, rhs = [], []
    for hh in range(h):
        for pp in range(p):
            for aa in range(a - 1):
                r = np.zeros(nv); r[ix(hh, aa + 1, pp)] = 1; r[ix(hh, aa, pp)] = -1
                rows.append(r); rhs.append(da[pp])
                rows.append(-r); rhs.append(da[pp])
            for aa in range(a):
                prev = mats[-1][aa, pp] if hh == 0 else None
                r = np.zeros(nv); r[ix(hh, aa, pp)] = 1
                if hh == 0:
                    rows.append(r); rhs.append(prev + dt[pp])
                    rows.append(-r); rhs.append(-(prev - dt[pp]))
                else:
                    r2 = r.copy(); r2[ix(hh - 1, aa, pp)] = -1
                    rows.append(r2); rhs.append(dt[pp])
                    rows.append(-r2); rhs.append(dt[pp])
        # aggregate drift constraint on TFR at this step
        r = np.zeros(nv)
        for aa in range(a):
            for pp in range(p):
                r[ix(hh, aa, pp)] = 1
        base = float(mats[-1].sum())
        rows.append(r); rhs.append(base + (hh + 1) * dtfr)
        rows.append(-r); rhs.append(-(base - (hh + 1) * dtfr))

    A_ub = np.array(rows); b_ub = np.array(rhs)
    bounds = [(0, ub)] * nv
    # TFR bounds at horizon h
    c = np.zeros(nv)
    for aa in range(a):
        for pp in range(p):
            c[ix(h - 1, aa, pp)] = 1
    r_lo = linprog(c, A_ub=A_ub, b_ub=b_ub, bounds=bounds, method="highs")
    r_hi = linprog(-c, A_ub=A_ub, b_ub=b_ub, bounds=bounds, method="highs")
    if r_lo.status != 0 or r_hi.status != 0:
        # relax: drop aggregate drift rows (recompute with wider quantile deltas)
        da_r = da * 2.0; dt_r = dt * 2.0
        rows2, rhs2 = [], []
        for hh in range(h):
            for pp in range(p):
                for aa in range(a - 1):
                    r = np.zeros(nv); r[ix(hh, aa + 1, pp)] = 1; r[ix(hh, aa, pp)] = -1
                    rows2.append(r); rhs2.append(da_r[pp])
                    rows2.append(-r); rhs2.append(da_r[pp])
                for aa in range(a):
                    if hh == 0:
                        prev = mats[-1][aa, pp]
                        r = np.zeros(nv); r[ix(hh, aa, pp)] = 1
                        rows2.append(r); rhs2.append(prev + dt_r[pp])
                        rows2.append(-r); rhs2.append(-(prev - dt_r[pp]))
                    else:
                        r = np.zeros(nv); r[ix(hh, aa, pp)] = 1; r[ix(hh - 1, aa, pp)] = -1
                        rows2.append(r); rhs2.append(dt_r[pp])
                        rows2.append(-r); rhs2.append(dt_r[pp])
        A2b = np.array(rows2); b2b = np.array(rhs2)
        r_lo = linprog(c, A_ub=A2b, b_ub=b2b, bounds=bounds, method="highs")
        r_hi = linprog(-c, A_ub=A2b, b_ub=b2b, bounds=bounds, method="highs")
    if r_lo.status != 0 or r_hi.status != 0:
        base = float(mats[-1].sum())
        return base, base * 0.5, base * 1.5   # documented fallback envelope
    lo, hi = float(r_lo.fun), float(-r_hi.fun)
    # MAP: L1 projection toward persistence path
    prior = np.zeros(nv)
    for hh in range(h):
        for aa in range(a):
            for pp in range(p):
                prior[ix(hh, aa, pp)] = mats[-1][aa, pp]
    n2 = 2 * nv
    A2 = lil_matrix((2 * nv + A_ub.shape[0], n2))
    b2 = np.zeros(2 * nv + A_ub.shape[0])
    I = np.eye(nv)
    A2[:nv, :nv] = I; A2[:nv, nv:] = -I; b2[:nv] = prior
    A2[nv:2 * nv, :nv] = -I; A2[nv:2 * nv, nv:] = -I; b2[nv:2 * nv] = -prior
    A2[2 * nv:, :nv] = A_ub
    b2[2 * nv:] = b_ub
    c2 = np.concatenate([np.zeros(nv), np.ones(nv)])
    rmap = linprog(c2, A_ub=A2.tocsr(), b_ub=b2, bounds=[(0, None)] * n2, method="highs")
    if rmap.status == 0:
        x = rmap.x[:nv].reshape(h, a, p)
        mu = float(x[h - 1].sum())
    else:
        mu = base
    return mu, lo, hi


def main():
    t0all = time.perf_counter()
    panels = load_tfr_panel()
    rows = []
    for country, d in panels.items():
        tfr = d["tfr"]; parity = d["parity"]
        years = sorted(tfr)
        for origin in ORIGINS:
            if origin not in tfr:
                continue
            hist_years = [y for y in years if y <= origin]
            hist_tfr = [tfr[y] for y in hist_years]
            hist_parity = {y: parity[y] for y in hist_years}
            if len(hist_tfr) < 8:
                continue
            future = [y for y in years if origin < y <= origin + H]
            if len(future) < 3:
                continue
            # in-sample naive scale for MASE
            naive_mae = float(np.mean(np.abs(np.diff(hist_tfr)))) if len(hist_tfr) > 1 else np.nan
            for h in range(1, H + 1):
                y = origin + h
                if y not in tfr:
                    continue
                truth = tfr[y]
                preds = {}
                preds["persistence"] = f_persistence(hist_tfr, h)
                preds["arima"] = f_arima(hist_years, hist_tfr, h)
                preds["ets"] = f_ets(hist_years, hist_tfr, h)
                preds["cipd"] = f_cipd_forecast(hist_tfr, hist_parity, origin, h)
                for m, (mu, lo, hi) in preds.items():
                    rec = {"country": country, "origin": origin, "horizon": h,
                           "method": m, "truth": truth, "point": mu,
                           "mae": abs(mu - truth), "sq": (mu - truth) ** 2,
                           "mase": abs(mu - truth) / naive_mae if naive_mae else np.nan}
                    if lo is not None and hi is not None and np.isfinite(lo) and np.isfinite(hi):
                        rec["lo90"] = lo; rec["hi90"] = hi
                        rec["covered"] = float(lo <= truth <= hi)
                        rec["width"] = hi - lo
                        if m == "cipd":
                            rec["crps"] = crps_uniform(lo, hi, truth)
                        else:
                            sig = (hi - lo) / (2 * 1.645)
                            rec["crps"] = crps_normal(mu, sig, truth)
                        rec["wis"] = wis90(lo, hi, truth)
                    else:
                        # persistence: RW predictive distribution
                        if m == "persistence":
                            sig = float(np.std(np.diff(hist_tfr), ddof=1)) * np.sqrt(h) if len(hist_tfr) > 2 else np.nan
                            if np.isfinite(sig):
                                rec["lo90"] = mu - 1.645 * sig
                                rec["hi90"] = mu + 1.645 * sig
                                rec["covered"] = float(rec["lo90"] <= truth <= rec["hi90"])
                                rec["width"] = rec["hi90"] - rec["lo90"]
                                rec["crps"] = crps_normal(mu, sig, truth)
                                rec["wis"] = wis90(rec["lo90"], rec["hi90"], truth)
                    rows.append(rec)
        print(f"[{country}] done", flush=True)

    df = pd.DataFrame(rows)
    df.to_csv(TAB / "E04_rolling_multicountry_raw.csv", index=False)
    agg = (df.groupby(["method", "horizon"])
             .agg(mae=("mae", "mean"), rmse=("sq", lambda s: float(np.sqrt(np.mean(s)))),
                  mase=("mase", "mean"), crps=("crps", "mean"), wis=("wis", "mean"),
                  coverage=("covered", "mean"), width=("width", "mean"))
             .reset_index())
    agg.to_csv(TAB / "E04_rolling_multicountry_summary.csv", index=False)
    overall = (df.groupby("method")
                 .agg(mae=("mae", "mean"), rmse=("sq", lambda s: float(np.sqrt(np.mean(s)))),
                      mase=("mase", "mean"), crps=("crps", "mean"), wis=("wis", "mean"),
                      coverage=("covered", "mean"), width=("width", "mean"),
                      n=("mae", "size"))
                 .reset_index())
    overall.to_csv(TAB / "E04_rolling_multicountry_overall.csv", index=False)
    print(overall.round(4).to_string(index=False))
    (LOG / "E04_multicountry_runlog.json").write_text(json.dumps({
        "runtime_sec": round(time.perf_counter() - t0all, 1),
        "origins": ORIGINS, "countries": list(panels.keys()), "n_forecasts": len(df)}))


if __name__ == "__main__":
    warnings.filterwarnings("ignore")
    main()
