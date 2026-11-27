#!/usr/bin/env python3
"""Merge bayesTFR E04 forecasts into the multi-method comparison table.

Reads results/tables/E04_bayestfr_raw.csv (country, origin, horizon, median,
q05, q95) and appends rows to E04_rolling_multicountry_raw.csv with identical
metric definitions, then regenerates the summary/overall tables and a compact
bayesTFR-specific summary.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
TAB = ROOT / "results/tables"
sys.path.insert(0, str(ROOT / "experiments/level3"))
from e04_rolling_multicountry import load_tfr_panel, crps_normal, wis90  # noqa


def main():
    p = TAB / "E04_bayestfr_raw.csv"
    if not p.exists() or p.stat().st_size < 50:
        print("no bayesTFR output; skipping merge")
        return
    bt = pd.read_csv(p)
    if len(bt) == 0:
        print("bayesTFR output empty; skipping")
        return
    panels = load_tfr_panel()
    # bayesTFR runs on WPP2019 TFR estimates; the truth series for scoring stays
    # the HFC observed TFR (identical role as for the other methods).
    rows = []
    for _, r in bt.iterrows():
        c = r.country
        if c not in panels:
            continue
        tfr = panels[c]["tfr"]
        y = int(r.origin + r.horizon)
        if y not in tfr:
            continue
        truth = tfr[y]
        hist = [tfr[yy] for yy in sorted(tfr) if yy <= r.origin]
        naive_mae = float(np.mean(np.abs(np.diff(hist)))) if len(hist) > 1 else np.nan
        mu, lo, hi = float(r["median"]), float(r.q05), float(r.q95)
        sig = max((hi - lo) / (2 * 1.645), 1e-9)
        rows.append({"country": c, "origin": int(r.origin), "horizon": int(r.horizon),
                     "method": "bayesTFR", "truth": truth, "point": mu,
                     "mae": abs(mu - truth), "sq": (mu - truth) ** 2,
                     "mase": abs(mu - truth) / naive_mae if naive_mae else np.nan,
                     "lo90": lo, "hi90": hi,
                     "covered": float(lo <= truth <= hi), "width": hi - lo,
                     "crps": crps_normal(mu, sig, truth),
                     "wis": wis90(lo, hi, truth)})
    if not rows:
        print("no scorable bayesTFR rows (truth years unavailable)")
        return
    add = pd.DataFrame(rows)
    raw = pd.read_csv(TAB / "E04_rolling_multicountry_raw.csv")
    raw = raw[raw.method != "bayesTFR"]
    base = load_tfr_panel()
    all_raw = pd.concat([raw, add], ignore_index=True)
    all_raw.to_csv(TAB / "E04_rolling_multicountry_raw.csv", index=False)
    agg = (all_raw.groupby(["method", "horizon"])
           .agg(mae=("mae", "mean"),
                rmse=("sq", lambda s: float(np.sqrt(np.mean(s)))),
                mase=("mase", "mean"), crps=("crps", "mean"), wis=("wis", "mean"),
                coverage=("covered", "mean"), width=("width", "mean"))
           .reset_index())
    agg.to_csv(TAB / "E04_rolling_multicountry_summary.csv", index=False)
    overall = (all_raw.groupby("method")
               .agg(mae=("mae", "mean"),
                    rmse=("sq", lambda s: float(np.sqrt(np.mean(s)))),
                    mase=("mase", "mean"), crps=("crps", "mean"),
                    wis=("wis", "mean"), coverage=("covered", "mean"),
                    width=("width", "mean"), n=("mae", "size"))
               .reset_index())
    overall.to_csv(TAB / "E04_rolling_multicountry_overall.csv", index=False)
    print(f"merged {len(add)} bayesTFR forecasts")
    print(overall.round(4).to_string(index=False))


if __name__ == "__main__":
    main()
