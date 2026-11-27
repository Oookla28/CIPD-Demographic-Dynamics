"""V5: Dynamic Policy Causal Identification & HonestDiD-Style Smoothness Bounds (DCPI)."""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from scipy.optimize import linprog


@dataclass
class CausalSensitivityResult:
    policy_name: str
    M: float
    post_years: list[int]
    cum_effect_lower: float
    cum_effect_upper: float
    avg_annual_lower: float
    avg_annual_upper: float
    sign_robust_positive: bool
    evidence_grade: str


def frechet_only_child_eligibility(p_husband: float, p_wife: float) -> dict[str, float]:
    """Exact Fréchet bounds for P(H U W) (at least one spouse is an only child) without independence assumption."""
    p_both_lo = max(0.0, p_husband + p_wife - 1.0)
    p_both_hi = min(p_husband, p_wife)
    p_either_lo = max(p_husband, p_wife)
    p_either_hi = min(1.0, p_husband + p_wife)
    # Strictly Selective Two-Child newly eligible = P(H xor W) = P(H U W) - P(H cap W)
    # Since P(H xor W) = p_H + p_W - 2 P(H cap W):
    p_xor_lo = p_husband + p_wife - 2.0 * p_both_hi
    p_xor_hi = p_husband + p_wife - 2.0 * p_both_lo
    return {
        "p_both_lower": p_both_lo,
        "p_both_upper": p_both_hi,
        "p_either_lower": p_either_lo,
        "p_either_upper": p_either_hi,
        "p_selective_newly_eligible_lower": max(0.0, p_xor_lo),
        "p_selective_newly_eligible_upper": min(1.0, p_xor_hi),
    }


def solve_honest_smoothness_bounds(
    years: np.ndarray,
    y_lo: np.ndarray,
    y_hi: np.ndarray,
    policy_start_year: int,
    eval_end_year: int,
    M: float,
    policy_name: str = "Policy",
    evidence_grade: str = "C2",
) -> CausalSensitivityResult:
    """Exact Rambachan-Roth Delta^SD(M) LP partial identification for dynamic counterfactual births.

    Model:
      B_t = B_t^0 + delta_t,  where y_lo[t] <= B_t <= y_hi[t] (incorporating NBS revision intervals!)
      delta_t = 0 for all t < policy_start_year (no anticipation)
      | (B_{t+1}^0 - B_t^0) - (B_t^0 - B_{t-1}^0) | <= M  for all interior t (second-difference smoothness)
      B_t^0 >= 0.

    Objective: min / max  sum_{t in [policy_start_year, eval_end_year]} delta_t.
    """
    years = np.asarray(years, dtype=int)
    y_lo = np.asarray(y_lo, dtype=float)
    y_hi = np.asarray(y_hi, dtype=float)
    N = len(years)

    # Variables: [B0_0..B0_{N-1}, delta_0..delta_{N-1}] (2N vars)
    # delta_t = 0 for t < policy_start_year
    bounds = []
    for t in range(N):
        bounds.append((0.0, None))  # B0_t >= 0
    for t in range(N):
        if years[t] < policy_start_year:
            bounds.append((0.0, 0.0))
        else:
            bounds.append((None, None))

    A_ub = []
    b_ub = []

    # Observation intervals: y_lo[t] <= B0_t + delta_t <= y_hi[t]
    for t in range(N):
        r1 = np.zeros(2 * N); r1[t] = 1.0; r1[N + t] = 1.0
        A_ub.append(r1); b_ub.append(y_hi[t])
        r2 = np.zeros(2 * N); r2[t] = -1.0; r2[N + t] = -1.0
        A_ub.append(r2); b_ub.append(-y_lo[t])

    # Second-difference smoothness on counterfactual B0: |B0_{t+1} - 2*B0_t + B0_{t-1}| <= M
    for t in range(1, N - 1):
        r1 = np.zeros(2 * N)
        r1[t - 1] = 1.0; r1[t] = -2.0; r1[t + 1] = 1.0
        A_ub.append(r1); b_ub.append(M)
        A_ub.append(-r1); b_ub.append(M)

    A_ub = np.array(A_ub)
    b_ub = np.array(b_ub)

    post_mask = (years >= policy_start_year) & (years <= eval_end_year)
    post_years = [int(y) for y in years[post_mask]]
    c = np.zeros(2 * N)
    c[N + np.where(post_mask)[0]] = 1.0

    res_lo = linprog(c=c, A_ub=A_ub, b_ub=b_ub, bounds=bounds, method="highs")
    res_hi = linprog(c=-c, A_ub=A_ub, b_ub=b_ub, bounds=bounds, method="highs")

    if res_lo.status != 0 or res_hi.status != 0:
        return CausalSensitivityResult(
            policy_name=policy_name, M=M, post_years=post_years,
            cum_effect_lower=float("nan"), cum_effect_upper=float("nan"),
            avg_annual_lower=float("nan"), avg_annual_upper=float("nan"),
            sign_robust_positive=False, evidence_grade=evidence_grade,
        )

    lo = float(res_lo.fun)
    hi = float(-res_hi.fun)
    k = max(1, len(post_years))
    return CausalSensitivityResult(
        policy_name=policy_name,
        M=M,
        post_years=post_years,
        cum_effect_lower=lo,
        cum_effect_upper=hi,
        avg_annual_lower=lo / k,
        avg_annual_upper=hi / k,
        sign_robust_positive=(lo > 0.0),
        evidence_grade=evidence_grade,
    )
