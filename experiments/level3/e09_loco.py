#!/usr/bin/env python3
"""E09 (Level-3 upgrade): Leave-One-Country-Out cross-country generalization.

Protocol (strict, no tuning on hold-out countries):
  1. For every country c and mask m, compute CIPD reconstruction error over a
     grid of constraint-confidence quantiles q in QGRID (structural parameter).
  2. Development tuning: q* = argmin over q of mean MAE across the 16 dev
     countries; hold-out evaluation: apply q* unchanged to GBRTENW, CHE, TWN.
  3. Full LOCO: for each country c (all 19), q*_{-c} = argmin over the other 18,
     then evaluate country c at q*_{-c}. The hold-out countries never enter any
     tuning set; their dev-time evaluation is the strictest test.
  4. Forecast generalization: CIPD rolling TFR forecasts at default q=0.90 for
     dev and hold-out countries separately (same code as E04).
Outputs: E09_loco_grid.csv, E09_loco_summary.csv, E09_holdout_vs_dev.csv
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "experiments/level3"))

from e02_multicountry_masking import (make_mask, MASKS, AGES, cipd_reconstruct, metrics)
from e04_rolling_multicountry import f_cipd_forecast

TAB = ROOT / "results/tables"
LOG = ROOT / "results/logs"

DEV = ["JPN", "KOR", "DEUTNP", "SWE", "NOR", "FIN", "DNK", "AUT",
       "CZE", "NLD", "HUN", "ITA", "ESP", "USA", "CAN", "POL"]
HOLDOUT = ["GBRTENW", "CHE", "TWN"]
ALL = DEV + HOLDOUT
QGRID = [0.75, 0.80, 0.85, 0.90, 0.95]
SEED = 101


def load_all_panels():
    df = pd.read_csv(ROOT / "data/external/hfd_wpp/hfc_asfr_bo_panel.csv")
    out = {}
    for c in ALL:
        sub = df[df.country == c]
        mats = {}
        for y, sy in sub.groupby("year"):
            sy = sy.set_index("age")
            if not set(AGES).issubset(sy.index):
                continue
            m = sy.loc[AGES, ["asfr_p1", "asfr_p2", "asfr_p3p"]].values.astype(float)
            if np.all(np.isfinite(m)):
                mats[int(y)] = m
        out[c] = mats
    return out


def pick_eval_year(mats):
    ys = sorted(mats)
    return ys[len(ys) // 2]


def main():
    t0 = time.perf_counter()
    panels = load_all_panels()
    rows = []
    for c in ALL:
        mats = panels[c]
        if len(mats) < 5:
            print(f"[skip] {c}: {len(mats)} complete years")
            continue
        t_star = pick_eval_year(mats)
        true = mats[t_star]
        for kind, level in MASKS:
            seed = SEED if level is not None else 0
            rng = np.random.RandomState(seed + 1000 * ord(kind))
            miss = make_mask(kind, level, true.shape, rng)
            for q in QGRID:
                X_map, lo, hi, diag = cipd_reconstruct(mats, t_star, miss, delta_q=q)
                mt = metrics(true, X_map, miss)
                cov = float(((true >= lo - 1e-9) & (true <= hi + 1e-9))[miss].mean())
                rows.append({"country": c, "role": "dev" if c in DEV else "holdout",
                             "eval_year": t_star, "mask": kind, "q": q,
                             "mae": mt["mae"], "rmse": mt["rmse"],
                             "coverage": cov,
                             "mean_width": float(np.mean((hi - lo)[miss]))})
        print(f"[{c}] done", flush=True)

    grid = pd.DataFrame(rows)
    grid.to_csv(TAB / "E09_loco_grid.csv", index=False)

    # ---- dev tuning / holdout evaluation ----
    dev_mean = (grid[grid.role == "dev"].groupby("q")
                .agg(mae=("mae", "mean"), coverage=("coverage", "mean"),
                     width=("mean_width", "mean")).reset_index())
    q_star = float(dev_mean.loc[dev_mean.mae.idxmin(), "q"])
    hold = (grid[(grid.role == "holdout") & (np.isclose(grid.q, q_star))]
            .groupby("country").agg(mae=("mae", "mean"), coverage=("coverage", "mean"),
                                    width=("mean_width", "mean")).reset_index())
    dev_at_q = (grid[(grid.role == "dev") & (np.isclose(grid.q, q_star))]
                .groupby("country").agg(mae=("mae", "mean"), coverage=("coverage", "mean"),
                                        width=("mean_width", "mean")).reset_index())
    summ = pd.concat([dev_at_q.assign(role="dev"), hold.assign(role="holdout")])
    summ.to_csv(TAB / "E09_holdout_vs_dev.csv", index=False)

    # ---- full LOCO ----
    loco_rows = []
    for c in sorted(grid.country.unique()):
        rest = grid[grid.country != c]
        per_q = rest.groupby("q").mae.mean()
        q_c = float(per_q.idxmin())
        sub = grid[(grid.country == c) & (np.isclose(grid.q, q_c))]
        loco_rows.append({"test_country": c,
                           "role": "holdout" if c in HOLDOUT else "dev",
                           "q_tuned_without_country": q_c,
                           "mae": float(sub.mae.mean()),
                           "coverage": float(sub.coverage.mean()),
                           "width": float(sub.mean_width.mean())})
    loco = pd.DataFrame(loco_rows)
    loco.to_csv(TAB / "E09_loco_summary.csv", index=False)

    # ---- forecast generalization (default q=0.90), ALL countries ----
    fc_rows = []
    for c in ALL:
        mats_c = panels[c]
        if len(mats_c) < 12:
            continue
        years = sorted(mats_c)
        tfr_c = {y: float(mats_c[y].sum()) for y in years}
        for origin in [2005, 2010, 2016]:
            if origin not in tfr_c:
                continue
            hist = {y: mats_c[y] for y in years if y <= origin}
            if len(hist) < 8:
                continue
            for h in [3, 5, 8]:
                y = origin + h
                if y not in tfr_c:
                    continue
                mu, lo, hi = f_cipd_forecast(None, hist, origin, h, q=0.90)
                truth = tfr_c[y]
                fc_rows.append({"country": c, "role": "dev" if c in DEV else "holdout",
                                "origin": origin, "horizon": h, "truth": truth,
                                "point": mu, "err": abs(mu - truth),
                                "covered": float(lo <= truth <= hi),
                                "width": hi - lo})
        print(f"[forecast {c}] done", flush=True)
    fcd = pd.DataFrame(fc_rows)
    fcd.to_csv(TAB / "E09_forecast_generalization_raw.csv", index=False)
    fc_summ = (fcd.groupby(["role", "horizon"])
               .agg(n=("err", "size"), mae=("err", "mean"),
                    coverage=("covered", "mean"), width=("width", "mean"))
               .reset_index())
    fc_summ.to_csv(TAB / "E09_forecast_generalization.csv", index=False)

    print("dev-tuned q* =", q_star)
    print(summ.round(4).to_string(index=False))
    print(loco.round(4).to_string(index=False))
    print(fc_summ.round(4).to_string(index=False))
    (LOG / "E09_runlog.json").write_text(json.dumps({
        "runtime_sec": round(time.perf_counter() - t0, 1),
        "q_tuned_development": q_star,
        "qgrid": QGRID, "seed": SEED,
        "n_grid_cases": len(grid)}))


if __name__ == "__main__":
    main()
