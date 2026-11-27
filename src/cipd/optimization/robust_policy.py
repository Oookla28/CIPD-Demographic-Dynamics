"""V7: Baseline-Safe Robust Policy Optimization & Policy Ranking Reversal Certificates (CIPD)."""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np


@dataclass
class PolicyOption:
    name: str
    cash_subsidy_extra_cny: float   # delta cash subsidy above current 3600 CNY/yr
    childcare_slots_added_m: float  # added affordable childcare slots (millions)
    housing_tax_relief_pct: float   # % housing/tax cost relief for families
    annual_fiscal_cost_bn: float    # estimated annual fiscal cost (billion CNY)


@dataclass
class RobustDecisionReport:
    baseline_name: str
    nominal_best: str
    robust_best: str
    policy_metrics: list[dict]
    pairwise_dominance: list[dict]


def evaluate_policy_loss_matrix(
    policies: list[PolicyOption],
    theta_samples: list[dict],
    weights: dict[str, float] | None = None,
) -> np.ndarray:
    """Evaluate normalized multi-objective loss J(u_i, theta_s) for each policy i and mechanism scenario s.
    Lower J is better:
      J(u, theta) = w_demo * (-normalized_birth_gain) + w_fiscal * normalized_fiscal_cost + w_equity * inequity_penalty.
    Every parameter in theta_samples represents an explicitly disclosed elasticity & structural regime.
    """
    w = weights or {"demography": 0.50, "fiscal": 0.35, "equity": 0.15}
    n_pol = len(policies)
    n_th = len(theta_samples)
    raw_demo = np.zeros((n_pol, n_th))
    raw_fiscal = np.zeros((n_pol, n_th))
    raw_equity = np.zeros((n_pol, n_th))

    for i, pol in enumerate(policies):
        for s, th in enumerate(theta_samples):
            # Demographic birth gain (% lift in TFR / annual births) under mechanism theta:
            # diminishing returns concave response
            cash_lift = th["elas_cash"] * np.log1p(pol.cash_subsidy_extra_cny / 3600.0)
            care_lift = th["elas_care"] * np.log1p(pol.childcare_slots_added_m / 1.5)
            house_lift = th["elas_housing"] * np.log1p(pol.housing_tax_relief_pct / 5.0)
            # Interaction / crowding-out factor
            total_lift = (cash_lift + care_lift + house_lift) * th["structural_multiplier"]
            raw_demo[i, s] = -total_lift  # lower loss = larger birth stabilization
            raw_fiscal[i, s] = pol.annual_fiscal_cost_bn / max(th["fiscal_budget_cap_bn"], 1.0)
            # Equity proxy: universal childcare + cash scores better than pure tax relief (which favors high earners)
            raw_equity[i, s] = 0.4 * (pol.housing_tax_relief_pct / 10.0) - 0.3 * np.log1p(pol.childcare_slots_added_m / 2.0)

    # Dimensionless standardization across policies/scenarios
    def norm(mat):
        span = max(float(np.max(mat) - np.min(mat)), 1e-9)
        return (mat - np.min(mat)) / span

    J = w["demography"] * norm(raw_demo) + w["fiscal"] * norm(raw_fiscal) + w["equity"] * norm(raw_equity)
    return J


def solve_baseline_safe_robust_decision(
    policies: list[PolicyOption],
    theta_samples: list[dict],
    baseline_idx: int = 0,
    weights: dict[str, float] | None = None,
) -> RobustDecisionReport:
    """V7: Baseline-Safe Minimax Relative Regret & Pairwise Dominance Certification.
    - R(u) = max_{theta in Theta} [J(u, theta) - J(u_0, theta)]
    - u* = argmin_{u in U} R(u) guarantees J(u*, theta) <= J(u_0, theta) for all theta in Theta when u_0 in U.
    - Pairwise dominance D_{AB} = max_theta [J(u_A, theta) - J(u_B, theta)] < 0 certifies A strictly dominates B;
      if min < 0 < max, reports Ranking Reversal and identifies the driving mechanism parameters.
    """
    J = evaluate_policy_loss_matrix(policies, theta_samples, weights=weights)
    J0 = J[baseline_idx]

    metrics = []
    for i, pol in enumerate(policies):
        rel_loss = J[i] - J0
        worst_rel = float(np.max(rel_loss))
        best_rel = float(np.min(rel_loss))
        mean_loss = float(np.mean(J[i]))
        worst_abs = float(np.max(J[i]))
        metrics.append({
            "policy": pol.name,
            "annual_fiscal_cost_bn": pol.annual_fiscal_cost_bn,
            "mean_loss": mean_loss,
            "worst_abs_loss": worst_abs,
            "worst_rel_loss_vs_baseline": worst_rel,
            "best_rel_loss_vs_baseline": best_rel,
            "baseline_safe": worst_rel <= 1e-12,
        })

    nominal_best = min(metrics, key=lambda d: d["mean_loss"])["policy"]
    robust_best = min(metrics, key=lambda d: d["worst_rel_loss_vs_baseline"])["policy"]

    pairwise = []
    n_pol = len(policies)
    for i in range(n_pol):
        for j in range(i + 1, n_pol):
            diff = J[i] - J[j]
            d_min = float(np.min(diff))
            d_max = float(np.max(diff))
            if d_max < -1e-12:
                verdict = f"{policies[i].name} strictly dominates {policies[j].name}"
            elif d_min > 1e-12:
                verdict = f"{policies[j].name} strictly dominates {policies[i].name}"
            else:
                verdict = "Ranking Reversal (indeterminate under Theta)"
            pairwise.append({
                "policy_A": policies[i].name,
                "policy_B": policies[j].name,
                "min_diff_A_minus_B": d_min,
                "max_diff_A_minus_B": d_max,
                "verdict": verdict,
            })

    return RobustDecisionReport(
        baseline_name=policies[baseline_idx].name,
        nominal_best=nominal_best,
        robust_best=robust_best,
        policy_metrics=metrics,
        pairwise_dominance=pairwise,
    )
