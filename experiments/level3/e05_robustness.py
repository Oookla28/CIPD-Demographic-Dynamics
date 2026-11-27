#!/usr/bin/env python3
"""E05 robustness: re-estimate the two continuous-intensity event studies and
the HonestDiD bounds after dropping Xinjiang (whose 单独二孩 date is only
'approximate'), plus a placebo-lead check. Outputs E05_robustness_excl_xinjiang.csv.
"""
import json, subprocess, sys
from pathlib import Path
import numpy as np, pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "experiments/level3"))
from e05_provincial_did import build_panel, event_study  # noqa

TAB = ROOT / "results/tables"
LOG = ROOT / "results/logs"

df = build_panel()
rows = []
fits = {}
for excl in ([], ["新疆"]):
    sub = df[~df.province.isin(excl)].copy()
    tag = "all31" if not excl else "excl_xinjiang"
    esA = event_study(sub, "D_z", 2014, -6, 5, ref=-1)
    esB = event_study(sub, "E_z", 2016, -6, 8, ref=-1)
    fits[tag] = {"A": esA, "B": esB}
    for nm, es in [("A_2014", esA), ("B_2016", esB)]:
        pre = np.array(es["beta"][: es["n_pre"]])
        post = np.array(es["beta"][es["n_pre"]:])
        rows.append({"sample": tag, "policy": nm,
                     "pre_max_abs": float(np.max(np.abs(pre))),
                     "post_mean": float(np.mean(post)),
                     "post_last": float(post[-1]) if len(post) else np.nan,
                     "n_obs": es["n_obs"]})
pd.DataFrame(rows).to_csv(TAB / "E05_robustness_excl_xinjiang.csv", index=False)

# HonestDiD for both samples (reuse the R script with the exported JSONs)
for tag, f in fits.items():
    (LOG / f"e05_event_study_{tag}_A.json").write_text(json.dumps(f["A"]))
    (LOG / f"e05_event_study_{tag}_B.json").write_text(json.dumps(f["B"]))
out = []
for tag in fits:
    for pol, fname in [("A", f"e05_event_study_{tag}_A.json"), ("B", f"e05_event_study_{tag}_B.json")]:
        src = LOG / fname
        # create a temp R input by pointing the R script at a custom file list
        r = subprocess.run([str(Path.home() / ".local/bin/Rscript"),
                            str(ROOT / "experiments/rmods/e05_honestdid_single.R"),
                            str(src), str(LOG / f"hd_{tag}_{pol}.json")],
                           capture_output=True, text=True)
        if r.returncode == 0 and (LOG / f"hd_{tag}_{pol}.json").exists():
            h = json.loads((LOG / f"hd_{tag}_{pol}.json").read_text())
            for rec in h["sensitivity"]:
                out.append({"sample": tag, "policy": pol, **rec})
pd.DataFrame(out).to_csv(TAB / "E05_robustness_honestdid.csv", index=False)
print(pd.DataFrame(rows).round(4).to_string(index=False))
if len(out):
    hd = pd.DataFrame(out)
    piv = hd[hd.method.str.contains("relative")].pivot_table(
        index=["sample", "policy"], columns="M", values=["lb", "ub"])
    print(piv.round(4).to_string())
