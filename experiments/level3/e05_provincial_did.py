#!/usr/bin/env python3
"""E05 (Level-3 upgrade): provincial quasi-experimental identification of the
three Chinese fertility-policy reforms with real provincial panels.

Data
----
* 31-province annual panel 2004-2024 (CBR/CDR/NGR/population/urbanisation)
  from China Statistical Yearbook tables (see build_provincial_panel.py).
* Provincial policy implementation dates (provincial_policy_timeline.csv):
  单独二孩 2014 implementation months vary by province (Jan 17 - Sep 17).
* 2010 census long-form table 6-5: provincial women 15-64 by children ever
  born (pre-treatment parity stock) -> one-child-family intensity E_i.

Designs
-------
A. 单独二孩 (2014): continuous-timing DiD event study.
   y_it = log(CBR_it); D_i = months exposed in 2014 / 12;
   y = a_i + g_t + sum_k b_k (D_i x 1[t=2014+k]) + e,  k=-6..+5 (k=-1 ref).
   Identification: within-2014 adoption timing is set by provincial
   legislature calendars; pre-trend coefficients test this.
B. 全面二孩 (2016): continuous-intensity DiD event study with the PRE-policy
   (2010 census) one-child share E_i = ceb1/(ceb_total - ceb0). k=-6..+8.
C. 三孩 (2021) + 2025 national subsidy: no cross-province variation in
   treatment -> honest C3 evidence: national before/after + placebo breaks.

HonestDiD: Delta^SD(M) robust CIs for M in {0,0.5,1,1.5,2,3} on both event
studies (R package HonestDiD via experiments/rmods/e05_honestdid.R).

`did` package: bcallaway11/did att_gt on the 2014 dose-contrast design
(above-median-timing provinces vs below-median) -- documented as a dose
contrast, since adoption is universal within 2014.
"""
from __future__ import annotations

import json
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
TAB = ROOT / "results/tables"
LOG = ROOT / "results/logs"

REGIONS = ["北京", "天津", "河北", "山西", "内蒙古", "辽宁", "吉林", "黑龙江",
           "上海", "江苏", "浙江", "安徽", "福建", "江西", "山东", "河南",
           "湖北", "湖南", "广东", "广西", "海南", "重庆", "四川", "贵州",
           "云南", "西藏", "陕西", "甘肃", "青海", "宁夏", "新疆"]


def load_inputs():
    panel = pd.read_csv(ROOT / "data/raw/china_provincial/china_provincial_panel_2004_2024.csv")
    panel = panel[panel.province != "全国"].copy()
    timeline = pd.read_csv(ROOT / "data/raw/china_provincial/provincial_policy_timeline.csv")
    dandui = timeline[timeline.policy == "DANDUI"].copy()
    dandui["eff"] = pd.to_datetime(dandui.effective_date)
    dandui["months_2014"] = ((pd.Timestamp("2014-12-31") - dandui.eff) / pd.Timedelta(days=30.44)).round(1)
    dandui["D"] = (dandui.months_2014 / 12).clip(0, 1)
    ceb = pd.read_csv(ROOT / "data/raw/china_provincial/census2010_table6_5_provincial_ceb.csv")
    ceb = ceb[ceb.region_type == "total"].copy()
    tot = ceb[["ceb0", "ceb1", "ceb2", "ceb3", "ceb4", "ceb5p"]].sum(axis=1)
    ceb["parous"] = tot - ceb.ceb0
    ceb["E_onechild"] = ceb.ceb1 / ceb.parous
    ceb["E_two_plus"] = (ceb.ceb2 + ceb.ceb3 + ceb.ceb4 + ceb.ceb5p) / ceb.parous
    return panel, dandui[["province", "effective_date", "months_2014", "D", "date_quality"]], \
        ceb[["province", "E_onechild", "E_two_plus", "parous"]]


def build_panel():
    panel, dandui, ceb = load_inputs()
    df = panel.merge(dandui, on="province", how="left").merge(ceb, on="province", how="left")
    df = df.dropna(subset=["cbr_permille", "pop_wan"])
    df["log_cbr"] = np.log(df.cbr_permille)
    df["births_wan"] = df.cbr_permille * df.pop_wan / 1000.0
    df["log_births"] = np.log(df.births_wan.clip(lower=1e-6))
    df["year"] = df.year.astype(int)
    # standardized intensities
    df["D_z"] = (df.D - df.D.mean()) / df.D.std()
    df["E_z"] = (df.E_onechild - df.E_onechild.mean()) / df.E_onechild.std()
    return df


def event_study(df, treat_col, policy_year, kmin, kmax, ref=-1, outcome="log_cbr"):
    """PyFixest TWFE event study with continuous treatment intensity."""
    import pyfixest as pf
    d = df.copy()
    d = d[(d.year >= policy_year + kmin) & (d.year <= policy_year + kmax)]
    # 2020 is unavailable in the provincial panel (census year); documented.
    def cname(k):
        return f"{treat_col}_k{'m' if k < 0 else 'p'}{abs(k)}"
    keep_cols = []
    ks_kept = []
    for k in range(kmin, kmax + 1):
        if k == ref:
            continue
        if not (d.year == policy_year + k).any():
            continue  # e.g. 2020: no provincial yearbook rates (census year)
        col = cname(k)
        d[col] = d[treat_col] * (d.year == policy_year + k).astype(float)
        keep_cols.append(col)
        ks_kept.append(k)
    fml = f"{outcome} ~ " + " + ".join(keep_cols) + " | province + year"
    m = pf.feols(fml, data=d, vcov={"CRV1": "province"})
    coefs = m.coef().to_dict()
    _v = getattr(m, "_vcov", None)
    if _v is None:
        try:
            _v = m.vcov("cluster")
        except Exception:
            _v = m.vcov()
    vcov = np.asarray(_v, dtype=float)
    names = list(m.coef().index)
    # order: pre (kmin..ref-1), post (ref+1..kmax)
    order = [cname(k) for k in ks_kept]
    idx = [names.index(o) for o in order]
    beta = np.array([coefs[o] for o in order])
    V = vcov[np.ix_(idx, idx)]
    pre = [k for k in ks_kept if k < ref]
    post = [k for k in ks_kept if k > ref]
    return {
        "k_pre": pre, "k_post": post,
        "beta": beta.tolist(), "vcov": V.tolist(),
        "n_pre": len(pre), "n_post": len(post),
        "n_obs": int(m._N), "formula": fml,
        "r2": float(m._r2) if hasattr(m, "_r2") else None,
        "coef_raw": coefs, "names": order,
    }


def national_break_test():
    """C3: national before/after for 三孩 (2021) and 2025 subsidy with placebo breaks."""
    cn = pd.read_csv(ROOT / "data/raw/china_annual_vital_1990_2025.csv")
    cn = cn[cn.year >= 2005].copy()
    cn["log_births"] = np.log(cn.births_wan)
    rows = []
    for brk, label in [(2021, "T3_2021"), (2025, "SUBSIDY_2025"),
                       (2014, "placebo_2014"), (2018, "placebo_2018")]:
        pre = cn[(cn.year >= brk - 5) & (cn.year < brk)]
        post = cn[(cn.year >= brk) & (cn.year < brk + 5)]
        if len(pre) < 3 or len(post) < 2:
            continue
        lvl_pre = pre.births_wan.mean(); lvl_post = post.births_wan.mean()
        # simple interrupted time-series slope change on log births
        x = np.arange(len(pre) + len(post))
        y = np.r_[pre.log_births.values, post.log_births.values]
        D = (x >= len(pre)).astype(float)
        X = np.column_stack([np.ones_like(x), x, D, x * D])
        beta, *_ = np.linalg.lstsq(X, y, rcond=None)
        resid = y - X @ beta
        s2 = resid @ resid / max(1, len(y) - 4)
        se = np.sqrt(np.diag(s2 * np.linalg.pinv(X.T @ X)))
        rows.append({"break_year": brk, "label": label,
                     "mean_pre_wan": float(lvl_pre), "mean_post_wan": float(lvl_post),
                     "pct_change": float((lvl_post / lvl_pre - 1) * 100),
                     "its_level_jump": float(beta[2]), "its_level_se": float(se[2]),
                     "its_slope_change": float(beta[3]), "its_slope_se": float(se[3]),
                     "evidence_grade": "C3" if "placebo" not in label else "placebo"})
    return rows


def run_did_r(df):
    """Dose-contrast CS DiD via bcallaway11/did (see R script)."""
    d = df[["province", "year", "log_births", "log_cbr", "D", "E_onechild"]].copy()
    d["high_timing"] = (d.D > d.D.median()).astype(int)
    d["high_E"] = (d.E_onechild > d.E_onechild.median()).astype(int)
    d.to_csv(LOG / "e05_did_input.csv", index=False)
    r = subprocess.run([str(Path.home() / ".local/bin/Rscript"),
                        str(ROOT / "experiments/rmods/e05_did_honest.R"),
                        str(LOG / "e05_did_input.csv"), str(LOG / "e05_did_output.json"),
                        str(LOG / "e05_honestevent.json")],
                       capture_output=True, text=True)
    print(r.stdout[-1500:])
    if r.returncode != 0:
        print(r.stderr[-2000:])
    return


def main():
    t0 = time.perf_counter()
    df = build_panel()
    df.to_csv(LOG / "e05_provincial_analysis_panel.csv", index=False)

    esA = event_study(df, "D_z", 2014, -6, 5, ref=-1)
    esB = event_study(df, "E_z", 2016, -6, 8, ref=-1)
    (LOG / "e05_event_study_A_2014.json").write_text(json.dumps(esA))
    (LOG / "e05_event_study_B_2016.json").write_text(json.dumps(esB))

    # coefficient tables
    for nm, es in [("A_dandui2014", esA), ("B_quanmian2016", esB)]:
        rows = []
        for k, b in zip(es["k_pre"] + es["k_post"], es["beta"]):
            rows.append({"policy": nm, "k": k, "beta": b})
        pd.DataFrame(rows).to_csv(TAB / f"E05_event_study_{nm}.csv", index=False)

    nat = national_break_test()
    pd.DataFrame(nat).to_csv(TAB / "E05_national_break_C3.csv", index=False)

    run_did_r(df)

    # merge HonestDiD output if present
    hp = LOG / "e05_honestevent_results.json"
    if hp.exists():
        h = json.loads(hp.read_text())
        pd.DataFrame(h["sensitivity"]).to_csv(TAB / "E05_honestdid_sensitivity.csv", index=False)

    print(pd.DataFrame(nat).round(4).to_string(index=False))
    print(f"runtime {time.perf_counter()-t0:.1f}s")
    (LOG / "E05_runlog.json").write_text(json.dumps({
        "runtime_sec": round(time.perf_counter() - t0, 1),
        "n_provinces": int(df.province.nunique()),
        "years": [int(df.year.min()), int(df.year.max())],
        "note_2020": "2020 missing from provincial panel (census year, no yearbook rate table); excluded."}))


if __name__ == "__main__":
    main()
