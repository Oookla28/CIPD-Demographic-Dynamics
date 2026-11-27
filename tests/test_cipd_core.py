"""Core unit tests for CIPD V0-V7 mathematical invariants, primal-dual certificates, and accounting."""
import sys
import unittest
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from cipd.demography.ccm import (
    run_ccm_v0,
    run_pccm_v1,
    adjoint_policy_gradient,
    decompose_births_direct_structure,
    shapley_quantum_tempo_structure,
)
from cipd.observation.flow_compiler import DemographicFlowCompiler
from cipd.identification.partial_id import solve_identified_bounds, diagnose_infeasibility
from cipd.counterfactual.coupler import solve_shared_mechanism_coupler
from cipd.causality.causal_bounds import frechet_only_child_eligibility, solve_honest_smoothness_bounds
from cipd.information_value.voi import solve_two_trajectory_separation, select_observations_greedy_voi, select_observations_exhaustive_optimal
from cipd.optimization.robust_policy import PolicyOption, solve_baseline_safe_robust_decision


class TestCIPDCore(unittest.TestCase):
    def test_v0_accounting_exact(self):
        T, A = 10, 15
        rng = np.random.default_rng(42)
        pm0 = rng.uniform(50, 100, size=A)
        pf0 = rng.uniform(50, 100, size=A)
        sm = np.full((T, A), 0.99); sm[:, -1] = 0.85
        sf = np.full((T, A), 0.992); sf[:, -1] = 0.88
        asfr = np.zeros((T, A)); asfr[:, 3:9] = 0.08
        in_m = rng.uniform(0, 2, size=(T, A))
        out_m = rng.uniform(0, 2, size=(T, A))
        res = run_ccm_v0(pm0, pf0, sm, sf, asfr, in_mig_m=in_m, out_mig_m=out_m)
        self.assertLess(res.max_rel_residual, 1e-12)
        self.assertTrue(np.all(res.pop_male >= 0))
        self.assertTrue(np.all(res.pop_female >= 0))

    def test_v1_pccm_parity_conservation_and_adjoint(self):
        T, A, P = 6, 10, 4
        pm0 = np.full(A, 100.0)
        pf0 = np.zeros((A, P)); pf0[:, 0] = 60.0; pf0[:, 1] = 30.0; pf0[:, 2] = 8.0; pf0[:, 3] = 2.0
        sm = np.full((T, A), 0.99)
        sf = np.full((T, A), 0.995)
        pb = np.zeros((T, A, P)); pb[:, 2:7, :] = np.array([0.08, 0.05, 0.02, 0.01])
        res = run_pccm_v1(pm0, pf0, sm, sf, pb)
        self.assertLess(res.max_rel_residual, 1e-12)
        # Verify terminal parity r=3 still generates births (not an absolute ceiling)
        self.assertGreater(res.births_by_parity[0, 3], 0.0)

        # Verify reverse-mode adjoint gradient against central finite differences
        u = np.full((T, P), 0.005)
        J0, grad_adj = adjoint_policy_gradient(pf0, sf, pb, u)
        eps = 1e-6
        grad_fd = np.zeros_like(u)
        for t in range(T):
            for r in range(P):
                up = u.copy(); up[t, r] += eps
                um = u.copy(); um[t, r] -= eps
                Jp, _ = adjoint_policy_gradient(pf0, sf, pb, up)
                Jm, _ = adjoint_policy_gradient(pf0, sf, pb, um)
                grad_fd[t, r] = (Jp - Jm) / (2 * eps)
        self.assertLess(np.max(np.abs(grad_adj - grad_fd)), 1e-6)

    def test_decompositions_exact(self):
        f0 = np.array([0.05, 0.08, 0.03]); x0 = np.array([100.0, 120.0, 90.0])
        f1 = np.array([0.04, 0.09, 0.02]); x1 = np.array([90.0, 110.0, 100.0])
        d = decompose_births_direct_structure(f0, x0, f1, x1)
        self.assertLess(d["identity_residual"], 1e-12)

        g0 = np.array([0.3, 0.5, 0.2]); g1 = np.array([0.2, 0.55, 0.25])
        s = shapley_quantum_tempo_structure(1.5, g0, x0, 1.1, g1, x1)
        self.assertLess(s["identity_residual"], 1e-12)

    def test_v2_v3_polytope_primal_dual_certificate_and_monotonicity(self):
        # Observational equivalence toy example from Section 8.2: 100 women (50 p0, 50 p1), 10 births observed
        T, A, P = 2, 4, 3
        surv = np.ones((T, A))
        comp = DemographicFlowCompiler(T, A, P, surv, srb=1.0, fertile_age_idx=[1, 2])
        obs_base = [
            comp.obs_parity_stock(0, [1], [0], 50.0, 50.0, "init_p0"),
            comp.obs_parity_stock(0, [1], [1], 50.0, 50.0, "init_p1"),
            comp.obs_parity_stock(0, [1], [2], 0.0, 0.0, "init_p2"),
            comp.obs_age_stock(0, [0, 2, 3], 0.0, 0.0, "other_ages_zero"),
            comp.obs_total_births(0, 10.0, 10.0, "total_b0"),
        ]
        poly1 = comp.compile(obs_base)
        # Target: 2nd-order births at t=0
        c = np.zeros(poly1.n_vars)
        c[comp.idx_birth(0, 1, 1)] = 1.0
        cert1 = solve_identified_bounds(poly1, c, "second_births_t0")
        self.assertAlmostEqual(cert1.lower, 0.0, places=6)
        self.assertAlmostEqual(cert1.upper, 10.0, places=6)
        self.assertLess(cert1.duality_gap_low, 1e-8)
        self.assertLess(cert1.duality_gap_high, 1e-8)
        self.assertLess(cert1.primal_residual_high, 1e-8)

        # Monotonicity (Proposition 3): add consistent observation b1 >= 4
        obs_more = obs_base + [comp.obs_parity_births(0, [0], 4.0, 10.0, "b1_ge_4")]
        cert2 = solve_identified_bounds(comp.compile(obs_more), c, "second_births_t0")
        self.assertGreaterEqual(cert2.lower, cert1.lower - 1e-8)
        self.assertLessEqual(cert2.upper, cert1.upper + 1e-8)
        self.assertAlmostEqual(cert2.upper, 6.0, places=6)

        # Infeasibility diagnostic
        obs_bad = obs_base + [comp.obs_parity_births(0, [0], 15.0, 15.0, "impossible_b1")]
        diag = diagnose_infeasibility(comp.compile(obs_bad))
        self.assertFalse(diag["feasible"])
        self.assertGreater(diag["total_slack"], 4.0)

    def test_v4_shared_mechanism_coupler(self):
        res = solve_shared_mechanism_coupler(
            total_women=100.0,
            obs_births_y0=(9.8, 10.2),
            p0_bounds=(0.05, 0.20),
            p1_bounds=(0.02, 0.15),
            x0_share_bounds=(0.35, 0.65),
            delta_p1_policy=0.05,
            horizon=4,
            tol=0.05,
            max_nodes=250,
        )
        self.assertLess(res.coupled_outer_width, res.naive_width)
        self.assertLessEqual(res.coupled_outer_lower, res.coupled_inner_lower + 1e-6)
        self.assertGreaterEqual(res.coupled_outer_upper, res.coupled_inner_upper - 1e-6)

    def test_v5_causal_bounds_and_frechet(self):
        fr = frechet_only_child_eligibility(0.4, 0.4)
        self.assertAlmostEqual(fr["p_either_lower"], 0.4)
        self.assertAlmostEqual(fr["p_either_upper"], 0.8)
        years = np.arange(2010, 2020)
        y = np.array([1600, 1600, 1600, 1600, 1600, 1600, 1750, 1700, 1600, 1600], dtype=float)
        res0 = solve_honest_smoothness_bounds(years, y, y, 2016, 2017, M=0.0)
        res1 = solve_honest_smoothness_bounds(years, y, y, 2016, 2017, M=50.0)
        self.assertAlmostEqual(res0.cum_effect_lower, 250.0, places=5)
        self.assertAlmostEqual(res0.cum_effect_upper, 250.0, places=5)
        self.assertLess(res1.cum_effect_lower, res0.cum_effect_lower)
        self.assertGreater(res1.cum_effect_upper, res0.cum_effect_upper)

    def test_v6_two_trajectory_voi(self):
        T, A, P = 2, 4, 3
        surv = np.ones((T, A))
        comp = DemographicFlowCompiler(T, A, P, surv, srb=1.0, fertile_age_idx=[1, 2])
        base_obs = [
            comp.obs_age_stock(0, [1], 100.0, 100.0, "w_a1"),
            comp.obs_age_stock(0, [0, 2, 3], 0.0, 0.0, "w_other"),
            comp.obs_total_births(0, 10.0, 10.0, "b_tot"),
        ]
        c = np.zeros(comp.n_vars); c[comp.idx_birth(0, 1, 1)] = 1.0
        cands = [
            comp.obs_parity_births(0, [0], 4.5, 5.5, "cand_b1", cost=1.0),
            comp.obs_parity_births(0, [1], 4.5, 5.5, "cand_b2", cost=1.0),
            comp.obs_parity_stock(0, [1], [0], 45.0, 55.0, "cand_x0", cost=1.0),
        ]
        g = select_observations_greedy_voi(comp, base_obs, cands, c, budget=1.0, use_two_trajectory=True)
        opt = select_observations_exhaustive_optimal(comp, base_obs, cands, c, budget=1.0, use_two_trajectory=True)
        self.assertAlmostEqual(g["final_width"], opt["final_width"], places=5)
        self.assertIn("cand_b2", g["selected_names"])

    def test_v7_baseline_safe_robust_decision(self):
        policies = [
            PolicyOption("U0_Current", 0.0, 0.0, 0.0, 90.0),
            PolicyOption("U1_CashHeavy", 3600.0, 0.5, 0.0, 180.0),
            PolicyOption("U2_BalancedCare", 1800.0, 3.0, 2.0, 150.0),
        ]
        thetas = [
            {"elas_cash": 0.02, "elas_care": 0.08, "elas_housing": 0.03, "structural_multiplier": 0.9, "fiscal_budget_cap_bn": 200.0},
            {"elas_cash": 0.06, "elas_care": 0.04, "elas_housing": 0.02, "structural_multiplier": 1.1, "fiscal_budget_cap_bn": 200.0},
            {"elas_cash": 0.01, "elas_care": 0.02, "elas_housing": 0.01, "structural_multiplier": 0.6, "fiscal_budget_cap_bn": 160.0},
        ]
        rep = solve_baseline_safe_robust_decision(policies, thetas, baseline_idx=0)
        best_row = next(r for r in rep.policy_metrics if r["policy"] == rep.robust_best)
        self.assertLessEqual(best_row["worst_rel_loss_vs_baseline"], 1e-12)


if __name__ == "__main__":
    unittest.main(verbosity=2)
