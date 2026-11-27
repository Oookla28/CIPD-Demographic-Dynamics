#!/usr/bin/env python3
"""Run (or merge) the Ellison et al. (2026) P-spline GAM + HMC baseline for E02."""
from __future__ import annotations
import json, subprocess, sys, time
from pathlib import Path
import numpy as np, pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src")); sys.path.insert(0, str(ROOT / "experiments/level3"))
TAB = ROOT / "results/tables"; LOG = ROOT / "results/logs"


def build_and_run(force=False):
    from e02_multicountry_masking import (load_panel, pick_eval_year, make_mask, MASKS, AGES)
    panels = load_panel()
    cases, meta = [], []
    for country, mats in panels.items():
        if len(mats) < 9: continue
        t_star = pick_eval_year(mats); true = mats[t_star]
        for kind, level in MASKS:
            seed = 101 if level is not None else 0
            rng = np.random.RandomState(seed + 1000 * ord(kind))
            miss = make_mask(kind, level, true.shape, rng)
            for j, pname in enumerate(["p1", "p2", "p3p"]):
                cid = f"{country}|{t_star}|{kind}|{seed}|{pname}"
                cases.append({"id": cid, "age": AGES.tolist(),
                              "y": np.nan_to_num(true[:, j], nan=0.0).tolist(),
                              "y_na": miss[:, j].tolist()})
                meta.append({"id": cid, "country": country, "eval_year": t_star,
                             "mask": kind, "seed": seed, "parity": pname, "miss": miss[:, j]})
    inp = LOG / "e02_ellison_input.json"; outp = LOG / "e02_ellison_output.json"
    inp.write_text(json.dumps({"cases": cases}))
    print(f"prepared {len(cases)} parity-cases")
    if force or not outp.exists():
        r = subprocess.run([str(Path.home() / ".local/bin/Rscript"),
                            str(ROOT / "experiments/rmods/e02_ellison_gam.R"),
                            str(inp), str(outp), "8", "2", "800"], capture_output=True, text=True)
        print(r.stdout[-500:])
        if r.returncode != 0:
            print(r.stderr[-2000:]); raise SystemExit("Rscript failed")
    return panels, meta, outp


def merge(panels, meta, outp):
    res = {c["id"]: c for c in json.loads(outp.read_text())["cases"]}
    groups = {}
    for m in meta:
        groups.setdefault((m["country"], m["eval_year"], m["mask"], m["seed"]), []).append(m)
    rows, skipped = [], 0
    for (country, t_star, mask, seed), ms in groups.items():
        mat = panels[country][t_star]
        pred = np.full_like(mat, np.nan); q05 = np.full_like(mat, np.nan); q95 = np.full_like(mat, np.nan)
        miss_all = np.zeros_like(mat, dtype=bool); ok = True
        for m in ms:
            j = ["p1", "p2", "p3p"].index(m["parity"]); c = res.get(m["id"])
            if c is None or c.get("method") == "insufficient_data": ok = False; break
            arr = np.array([float(v) if v not in ("NA", None) else np.nan for v in c["pred"]])
            if arr.shape[0] != 35 or not np.all(np.isfinite(arr)): ok = False; break
            pred[:, j] = arr
            q05[:, j] = np.array([float(v) if v not in ("NA", None) else np.nan for v in c["q05"]])
            q95[:, j] = np.array([float(v) if v not in ("NA", None) else np.nan for v in c["q95"]])
            miss_all[:, j] = m["miss"]
        if not ok: skipped += 1; continue
        e = np.abs(pred - mat)[miss_all]
        denom = max(float(np.mean(mat[miss_all])), 1e-9)
        rows.append({"country": country, "eval_year": t_star, "mask": mask, "level": "",
                     "seed": seed, "method": "ellison_gam", "mae": float(np.mean(e)),
                     "rmse": float(np.sqrt(np.mean(e ** 2))), "mare_pct": float(np.mean(e) / denom * 100),
                     "full_matrix_mae": float(np.mean(np.abs(pred - mat))),
                     "coverage_90": float(((mat >= q05) & (mat <= q95))[miss_all].mean()),
                     "mean_width_90": float(np.mean((q95 - q05)[miss_all]))})
    df = pd.DataFrame(rows)
    df.to_csv(TAB / "E02_multicountry_ellison.csv", index=False)
    print(f"merged {len(rows)} pairs, skipped {skipped}")
    if len(df):
        print(df.groupby("mask").agg(n=("mae", "size"), mae_mean=("mae", "mean"),
                                     coverage=("coverage_90", "mean"), width=("mean_width_90", "mean")).round(5).to_string())


if __name__ == "__main__":
    t0 = time.perf_counter()
    panels, meta, outp = build_and_run(force=False)
    merge(panels, meta, outp)
    print(f"runtime {time.perf_counter()-t0:.1f}s")
