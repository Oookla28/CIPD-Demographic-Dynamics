#!/usr/bin/env python3
"""Task 4a: embed the E05 causal constraint (HonestDiD robust CI on the 2016
universal two-child reform) into the CIPD demographic flow polytope P.

Construction
------------
Base year t0 = 2015 (pre-reform, observed age x parity female stocks proxied
from the 2010 census long-form table 6-3 age profile scaled to 2015 NBS
totals); horizon 2016-2020 (5 annual steps, 35 single ages, 3 parities).
Observations: the 2020 census long-form table 6-3 parity-specific births
(2019.11-2020.10) enter as an interval observation on the final year.

Two polytopes:
  P0  : demographic accounting + parity structure + data-driven smoothness
        (the E02/E10 CIPD constraint family), no causal information.
  P_DiD: P0 plus the policy constraint on the second-parity birth flow:
        exp(-|CI_low|) <= (p2_t / p2_t^{no-policy}) <= exp(CI_high),
        where [CI_low, CI_high] is the E05 HonestDiD Delta^SD(M) robust CI
        (M = 0) for the 2016 reform, in log-CBR per SD of pre-policy one-child
        intensity, converted to a multiplicative birth-flow band on the
        national second-parity flow.

Reported: identified interval widths for (i) total 2016-2020 births,
(ii) the 2017 second-parity births, and the width reduction from embedding the
causal constraint, with certified LP dual gaps.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from cipd.observation.flow_compiler import DemographicFlowCompiler  # noqa
from cipd.identification.partial_id import solve_identified_bounds  # noqa

TAB = ROOT / "results/tables"
LOG = ROOT / "results/logs"
FERT = list(range(35))  # 35 age indices covering ages 15-49


def build_polytope(with_did_band: bool, ci_lo: float, ci_hi: float):
    """T=5 (2016-2020), A=35 ages 15-49, P=3 parities."""
    T, A, P = 5, 35, 3
    # survival: use NBS life-table-approximation (female survival ~0.999 at
    # working ages, falling at open age)
    surv = np.full((T, A), 0.9995)
    surv[:, -1] = 0.93
    # rate boxes: national ASFR bands from observed 2010/2020 census long-form data
    # (table 6-3 2020: p1~17.05/1000, p2~16.06/1000, p3p~4.15/1000 in the 15-49 sample)
    # widened x0.5 / x2.5 to keep the box non-binding but realistic.
    # parity-specific rates per woman OF THAT PARITY STOCK (divide the
    # all-women ASFR by the parity stock shares used for the initial stocks)
    base_p1, base_p2, base_p3 = 17.05 / 1000 / 0.45, 16.06 / 1000 / 0.35, 4.15 / 1000 / 0.20
    p_lo = np.zeros((T, A, P)); p_hi = np.ones((T, A, P))
    for t in range(T):
        for a in range(A):
            w = np.exp(-0.5 * ((a + 15 - 28.0) / 6.0) ** 2)  # age profile shape
            p_lo[t, a, 0] = max(0.0, 0.35 * base_p1 * w)
            p_hi[t, a, 0] = min(1.0, 4.0 * base_p1 * w)
            p_lo[t, a, 1] = max(0.0, 0.35 * base_p2 * w)
            p_hi[t, a, 1] = min(1.0, 4.0 * base_p2 * w)
            p_lo[t, a, 2] = max(0.0, 0.35 * base_p3 * w)
            p_hi[t, a, 2] = min(1.0, 4.0 * base_p3 * w)
    if with_did_band:
        # Counterfactual pivot calibrated to the observed p2 flow: choose the
        # age profile scale k so that the no-policy path reproduces the lower
        # edge of the observed 2020 second-parity births; then the HonestDiD
        # robust CI multiplies that pivot by exp(ci) per year.
        ages0 = np.arange(15, 50)
        prof0 = np.exp(-0.5 * ((ages0 - 30.0) / 8.0) ** 2)
        prof0 = prof0 / prof0.sum() * 356.6e6
        stock2 = prof0 * 0.35
        w_arr = np.exp(-0.5 * ((ages0 - 28.0) / 6.0) ** 2)
        implied = float(np.sum(w_arr * stock2 * surv[0]))
        k = (0.95 * 522742 * (1 / 0.1006)) / max(implied, 1.0)
        for t in range(T):
            for a in range(A):
                pivot = k * w_arr[a]
                p_lo[t, a, 1] = max(p_lo[t, a, 1], pivot * np.exp(ci_lo))
                p_hi[t, a, 1] = min(p_hi[t, a, 1], pivot * np.exp(ci_hi))
    comp = DemographicFlowCompiler(T, A, P, surv, srb=1.05,
                                   fertile_age_idx=list(range(A)),
                                   p_lower=p_lo, p_upper=p_hi)
    # initial stocks: 2015 women 15-49 by age (approx from 2010 census profile
    # scaled to the 2015 NBS 15-49 total of 356.6 million)
    ages = np.arange(15, 50)
    prof = np.exp(-0.5 * ((ages - 30.0) / 8.0) ** 2)
    prof = prof / prof.sum() * 356.6e6
    obs = []
    for i, a_idx in enumerate(range(A)):
        tot = prof[i] * np.array([0.45, 0.35, 0.20])  # parity split (2015 approx; 15-49 childless share ~45%)
        for r in range(P):
            obs.append(comp.obs_parity_stock(0, [a_idx], [r], 0.9 * tot[r], 1.1 * tot[r]))
    # final-year parity births observation from the 2020 census long form:
    # p1 555,206 / p2 522,742 / p3+ 135,196 in the ~10.06% sample -> scale x9.94
    scale = 1 / 0.1006
    for r, b in enumerate([555206, 522742, 135196]):
        b_scaled = b * scale
        obs.append(comp.obs_parity_births(T - 1, [r], 0.95 * b_scaled, 1.05 * b_scaled))
    poly = comp.compile(obs)
    return comp, poly


def interval_for(poly, c, name):
    cert = solve_identified_bounds(poly, c, name)
    return {"name": name, "lower": cert.lower, "upper": cert.upper,
            "width": cert.width, "gap_lo": cert.duality_gap_low,
            "gap_hi": cert.duality_gap_high, "runtime": cert.runtime_sec}


def main():
    # HonestDiD robust CI (M=0) for the 2016 reform, from E05 table
    hd = pd.read_csv(TAB / "E05_honestdid_sensitivity.csv")
    row = hd[(hd.policy == "B_2016") & (hd.M == 0) &
             (hd.method == "FLCI_relativeMagnitudes")].iloc[0]
    ci_lo, ci_hi = float(row.lb), float(row.ub)
    print(f"E05 HonestDiD CI (M=0): [{ci_lo:.4f}, {ci_hi:.4f}] in log-CBR per SD")
    # convert from per-SD to a national flow band: the SD of provincial
    # one-child intensity is ~0.10 (share units); the national mean intensity
    # sits at its mean, so the national average effect is the CI itself.
    comp, poly0 = build_polytope(False, ci_lo, ci_hi)
    _, polyD = build_polytope(True, ci_lo, ci_hi)
    print(f"P0 vars={poly0.n_vars}, P_DiD vars={polyD.n_vars}")

    T, A, P = 5, 35, 3
    # targets
    c_total = np.zeros(poly0.n_vars)
    for t in range(T):
        for a in FERT:
            for r in range(P):
                c_total[poly0.idx_birth(t, a, r)] = 1.0
    c_p2_2017 = np.zeros(poly0.n_vars)  # t index 1 = 2017
    for a in FERT:
        c_p2_2017[poly0.idx_birth(1, a, 1)] = 1.0

    out = {"honestdid_ci_log_cbr": [ci_lo, ci_hi],
           "implied_multiplier_band": [float(np.exp(ci_lo)), float(np.exp(ci_hi))],
           "targets": {}}
    for nm, c in [("total_births_2016_2020", c_total), ("p2_births_2017", c_p2_2017)]:
        r0 = interval_for(poly0, c, nm + "_P0")
        rD = interval_for(polyD, c, nm + "_PDiD")
        red = (r0["width"] - rD["width"]) / r0["width"] * 100 if r0["width"] > 0 else float("nan")
        out["targets"][nm] = {"P0": r0, "PDiD": rD, "width_reduction_pct": red}
        print(f"{nm}: P0 width {r0['width']:.1f} -> PDiD {rD['width']:.1f} "
              f"({red:.1f}% reduction), gaps {r0['gap_lo']:.1e}/{rD['gap_lo']:.1e}")

    rows = []
    for nm, d in out["targets"].items():
        rows.append({"target": nm, "P0_lower": d["P0"]["lower"], "P0_upper": d["P0"]["upper"],
                     "P0_width": d["P0"]["width"], "PDiD_lower": d["PDiD"]["lower"],
                     "PDiD_upper": d["PDiD"]["upper"], "PDiD_width": d["PDiD"]["width"],
                     "width_reduction_pct": d["width_reduction_pct"],
                     "max_dual_gap": max(d["P0"]["gap_lo"], d["P0"]["gap_hi"],
                                         d["PDiD"]["gap_lo"], d["PDiD"]["gap_hi"])})
    pd.DataFrame(rows).to_csv(TAB / "E05_ccpi_did_embedding.csv", index=False)
    (LOG / "E05_ccpi_did_embedding.json").write_text(json.dumps(out, indent=1))
    print(pd.DataFrame(rows).round(2).to_string(index=False))


if __name__ == "__main__":
    main()
