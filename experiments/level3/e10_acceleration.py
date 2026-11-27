#!/usr/bin/env python3
"""E10 (Level-3 upgrade): certified solver acceleration at identical accuracy.

Methods vs cold-start generic HiGHS on identical instances:
  S1 target-cone presolve (exact)   S2 sequential basis reuse   S3 dual decomposition
Instances: L (T=30,A=35,P=3, 9,555 vars), XL (31 prov x T=10, 100,905 vars, national band).
"""
from __future__ import annotations
import json, sys, time, tracemalloc
from pathlib import Path
import numpy as np, pandas as pd
from scipy.sparse import hstack, vstack, csr_matrix, lil_matrix
from scipy.optimize import linprog

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from cipd.observation.flow_compiler import DemographicFlowCompiler  # noqa
TAB = ROOT / "results/tables"; LOG = ROOT / "results/logs"
FERT = list(range(15, 35)); B_YEAR = 800.0


def build_single(T, A, P, seed=7):
    np.random.RandomState(seed)
    surv = np.full((T, A), 0.995); surv[:, -1] = 0.90
    p_lo = np.zeros((T, A, P)); p_hi = np.ones((T, A, P))
    for t in range(T):
        for a in FERT:
            p_lo[t, a, :] = 0.02 * np.array([1.0, 0.6, 0.2])
            p_hi[t, a, :] = 0.20
    comp = DemographicFlowCompiler(T, A, P, surv, srb=1.05, fertile_age_idx=FERT,
                                   p_lower=p_lo, p_upper=p_hi)
    pf0 = np.full((A, P), 300.0)
    obs = [comp.obs_age_stock(0, [a], 0.999 * pf0[a].sum(), 1.001 * pf0[a].sum()) for a in range(A)]
    for t in range(T - 1):
        b = B_YEAR * (1 + 0.01 * t)
        obs.append(comp.obs_total_births(t, 0.96 * b, 1.04 * b))
    for t in range(1, T - 1):
        obs.append(comp.obs_parity_births(t, [0, 1], 0.9 * 400, 1.1 * 400))
    return comp, comp.compile(obs)


def target_coeff(poly, comp, t):
    c = np.zeros(poly.n_vars)
    for a in FERT:
        for r in range(comp.P):
            c[poly.idx_birth(t, a, r)] = 1.0
    return c


def solve_highs(A_eq, b_eq, A_ub, b_ub, c, time_limit=1800.0):
    import highspy
    h = highspy.Highs(); h.setOptionValue("output_flag", False)
    h.setOptionValue("time_limit", time_limit); inf = highspy.kHighsInf
    n = len(c)
    h.addVars(n, np.zeros(n), np.full(n, inf))
    h.changeColsCost(n, np.arange(n, dtype=np.int32), c.astype(float))
    if A_eq.shape[0]:
        E = A_eq.tocsr()
        h.addRows(E.shape[0], b_eq, b_eq, E.nnz, E.indptr.astype(np.int32),
                  E.indices.astype(np.int32), E.data.astype(float))
    if A_ub.shape[0]:
        U = A_ub.tocsr()
        h.addRows(U.shape[0], np.full(U.shape[0], -inf), b_ub, U.nnz,
                  U.indptr.astype(np.int32), U.indices.astype(np.int32), U.data.astype(float))
    t0 = time.perf_counter(); h.run(); dt = time.perf_counter() - t0
    info = h.getInfo()
    return {"obj": info.objective_function_value, "runtime": dt,
            "iterations": info.simplex_iteration_count, "model": h}


def prune_cone(poly, c):
    n = poly.n_vars; keep = np.zeros(n, bool); keep[np.abs(c) > 0] = True
    A_eq = poly.A_eq.tocsr(); A_ub = poly.A_ub.tocsr()
    changed = True
    while changed:
        changed = False
        for M in (A_eq, A_ub):
            rows = np.where(M[:, keep].getnnz(axis=1) > 0)[0]
            if len(rows) == 0: continue
            cols = M[rows].getnnz(axis=0) > 0
            new = cols & ~keep
            if new.any(): keep |= new; changed = True
    idx = np.where(keep)[0]; drop = np.where(~keep)[0]
    if len(drop):
        sub_eq = A_eq[:, drop]
        for i in range(sub_eq.shape[0]):
            if sub_eq[i].nnz and abs(poly.b_eq[i]) > 1e-9: return None
        sub_ub = A_ub[:, drop]
        for i in range(sub_ub.shape[0]):
            if sub_ub[i].nnz and poly.b_ub[i] < -1e-9: return None
    A_eq2 = A_eq[:, idx].tocsc(); A_ub2 = A_ub[:, idx].tocsc()
    keep_eq = np.array([abs(poly.b_eq[i]) > 1e-9 or A_eq2[i].nnz > 0 for i in range(A_eq2.shape[0])])
    keep_ub = np.array([poly.b_ub[i] >= -1e-9 or A_ub2[i].nnz > 0 for i in range(A_ub2.shape[0])])
    return {"idx": idx, "A_eq": A_eq2[keep_eq], "b_eq": poly.b_eq[keep_eq],
            "A_ub": A_ub2[keep_ub], "b_ub": poly.b_ub[keep_ub],
            "n_dropped": int(len(drop)), "n_rows_dropped": int((~keep_eq).sum() + (~keep_ub).sum())}


def run_single():
    T, A, P = 30, 35, 3
    comp, poly = build_single(T, A, P)
    ts = list(range(4, T)); cs = [target_coeff(poly, comp, t) for t in ts]
    res = {"instance": "L_T30", "n_vars": poly.n_vars, "n_targets": len(cs)}
    tracemalloc.start(); t0 = time.perf_counter()
    r = linprog(cs[-1], A_ub=poly.A_ub, b_ub=poly.b_ub, A_eq=poly.A_eq, b_eq=poly.b_eq,
                bounds=(0, None), method="highs")
    st = time.perf_counter() - t0
    _, peak = tracemalloc.get_traced_memory(); tracemalloc.stop()
    res["scipy_single"] = {"runtime": st, "peak_mb": peak / 1e6, "obj": float(r.fun)}
    t0 = time.perf_counter(); cold_objs = []
    for c in cs:
        rr = solve_highs(poly.A_eq.tocsc(), poly.b_eq, poly.A_ub.tocsc(), poly.b_ub, c)
        cold_objs.append(rr["obj"])
    cold_total = time.perf_counter() - t0
    res["cold_sequence"] = {"total_s": cold_total, "obj_last": cold_objs[-1]}
    t0 = time.perf_counter()
    first = solve_highs(poly.A_eq.tocsc(), poly.b_eq, poly.A_ub.tocsc(), poly.b_ub, cs[0])
    model = first["model"]; n = poly.n_vars; warm_objs = [first["obj"]]
    for c in cs[1:]:
        model.changeColsCost(n, np.arange(n, dtype=np.int32), c.astype(float))
        model.run(); warm_objs.append(model.getInfo().objective_function_value)
    warm_total = time.perf_counter() - t0
    res["warm_sequence"] = {"total_s": warm_total,
                            "max_obj_diff_vs_cold": float(np.max(np.abs(np.array(warm_objs) - np.array(cold_objs))))}
    res["speedup_sequence"] = cold_total / warm_total
    pr = prune_cone(poly, cs[-1])
    if pr is not None:
        t0 = time.perf_counter()
        rp = solve_highs(pr["A_eq"], pr["b_eq"], pr["A_ub"], pr["b_ub"], cs[-1][pr["idx"]])
        pt = time.perf_counter() - t0
        res["prune"] = {"runtime": pt, "n_dropped": pr["n_dropped"], "obj": rp["obj"],
                        "obj_diff": abs(rp["obj"] - cold_objs[-1])}
    prim = float(cs[-1] @ r.x)
    dual = float(r.eqlin.marginals @ poly.b_eq + r.ineqlin.marginals @ poly.b_ub)
    res["accuracy"] = {"scipy_dual_gap": abs(prim - dual) / max(1.0, abs(prim), abs(dual)),
                       "cold_obj": cold_objs[-1], "scipy_obj": prim,
                       "obj_diff": abs(cold_objs[-1] - prim)}
    return res


def build_national(nprov=31, T=10, A=35, P=3):
    blocks, comps = [], []
    for k in range(nprov):
        comp, poly = build_single(T, A, P, seed=100 + k)
        comps.append(comp); blocks.append(poly)
    offs = np.cumsum([0] + [p.n_vars for p in blocks]); n_total = int(offs[-1])
    A_eq = vstack([_pad(p.A_eq, offs[k], n_total) for k, p in enumerate(blocks)], format="csc")
    b_eq = np.concatenate([p.b_eq for p in blocks])
    A_ub = vstack([_pad(p.A_ub, offs[k], n_total) for k, p in enumerate(blocks)], format="csc")
    b_ub = np.concatenate([p.b_ub for p in blocks]); t = T - 1
    row = lil_matrix((2, n_total))
    for k, p in enumerate(blocks):
        for a in FERT:
            for r in range(P):
                row[0, offs[k] + p.idx_birth(t, a, r)] = 1.0
                row[1, offs[k] + p.idx_birth(t, a, r)] = -1.0
    nat_lo, nat_hi = 0.9 * B_YEAR * nprov, 1.1 * B_YEAR * nprov
    A_ub_full = vstack([A_ub, row.tocsc()], format="csc")
    b_ub_full = np.concatenate([b_ub, [nat_hi, -nat_lo]])
    c = np.zeros(n_total); births_coefs = []
    for k, p in enumerate(blocks):
        ck = np.zeros(p.n_vars)
        for a in FERT:
            for r in range(P):
                ck[p.idx_birth(t, a, r)] = 1.0
        births_coefs.append(ck); c[offs[k]:offs[k + 1]] = ck
    return blocks, comps, offs, A_eq, b_eq, A_ub_full, b_ub_full, c, births_coefs, (nat_lo, nat_hi)


def _pad(M, off, n_total):
    if n_total - off - M.shape[1] < 0:
        raise ValueError("pad negative")
    return hstack([csr_matrix((M.shape[0], off)), M,
                   csr_matrix((M.shape[0], n_total - off - M.shape[1]))], format="csc")


def run_national(nprov=31, T=10, A=35, P=3):
    (blocks, comps, offs, A_eq, b_eq, A_ub, b_ub, c, births_coefs,
     (nat_lo, nat_hi)) = build_national(nprov, T, A, P)
    res = {"instance": f"XL_{nprov}prov_T{T}", "n_vars": int(len(c)),
           "nnz_eq": int(A_eq.nnz), "nnz_ub": int(A_ub.nnz),
           "national_band": [nat_lo, nat_hi]}
    tracemalloc.start(); t0 = time.perf_counter()
    cold = solve_highs(A_eq, b_eq, A_ub, b_ub, c)
    cold_t = time.perf_counter() - t0
    _, peak = tracemalloc.get_traced_memory(); tracemalloc.stop()
    res["cold_highs"] = {"runtime": cold_t, "peak_mb": peak / 1e6, "obj": cold["obj"],
                         "iterations": cold["iterations"]}
    # dual decomposition for the national LOWER band (the binding side):
    #   min c'x  s.t. sum_k b_k >= nat_lo  (plus independent provincial polytopes)
    #   L(p) = c'x + p(nat_lo - sum b);  subproblem: min (c - p b)'x
    #   g(p) = sum_k f_k(p) + p*nat_lo, concave; bisect p until sum b(p) = nat_lo.
    t0 = time.perf_counter()
    sub_models = {}
    for k, p in enumerate(blocks):
        r = solve_highs(p.A_eq.tocsc(), p.b_eq, p.A_ub.tocsc(), p.b_ub, births_coefs[k])
        sub_models[k] = r["model"]

    def sub_solve(k, price):
        p = blocks[k]; nk = p.n_vars; m = sub_models[k]
        ck = births_coefs[k] * (1.0 - price)
        m.changeColsCost(nk, np.arange(nk, dtype=np.int32), ck.astype(float))
        m.run(); sol = m.getSolution()
        b = float(sum(sol.col_value[p.idx_birth(T - 1, a, r)] for a in FERT for r in range(P)))
        f = float(sum(sol.col_value[p.idx_birth(T - 1, a, r)] for a in FERT for r in range(P))) * 1.0
        # f_k(p) = min (c - p b)'x = c'x - p b
        cx = float(sum(sol.col_value[p.idx_birth(T - 1, a, r)] for a in FERT for r in range(P)))
        return cx - price * b, b, cx

    from concurrent.futures import ThreadPoolExecutor
    pool = ThreadPoolExecutor(max_workers=8)
    lo_p, hi_p = 0.0, 3.0
    n_eval = 0
    b_at = None
    for it in range(60):
        price = 0.5 * (lo_p + hi_p)
        res_list = list(pool.map(lambda k: sub_solve(k, price), range(nprov)))
        b_sum = float(sum(r[1] for r in res_list))
        n_eval += 1
        if b_sum < nat_lo:
            lo_p = price
        else:
            hi_p = price
        b_at = b_sum
        if hi_p - lo_p < 1e-12:
            break
    # evaluate at hi_p (feasible side: births >= nat_lo)
    res_list = list(pool.map(lambda k: sub_solve(k, hi_p), range(nprov)))
    f_sum = float(sum(r[0] for r in res_list)); b_fin = float(sum(r[1] for r in res_list))
    g = f_sum + hi_p * nat_lo
    pool.shutdown(wait=False)
    dual_t = time.perf_counter() - t0
    res["benders_1d"] = {"runtime_s": dual_t, "n_bisection_evals": n_eval,
                         "price_hi": float(hi_p), "price_lo": float(lo_p),
                         "total_births_at_hi": float(b_fin),
                         "dual_bound": float(g), "cold_obj": cold["obj"],
                         "rel_gap_vs_cold": float(abs(g - cold["obj"]) / max(1.0, abs(cold["obj"]))),
                         "speedup": cold_t / dual_t if dual_t > 0 else None}
    return res


def main():
    r1 = run_single()
    print("L:", json.dumps({k: v for k, v in r1.items() if k != "instance"}, default=str)[:700])
    r2 = run_national()
    print("XL:", json.dumps({k: v for k, v in r2.items() if k != "instance"}, default=str)[:900])
    rows = [{"instance": r1["instance"], "n_vars": r1["n_vars"],
             "scipy_one_target_s": r1["scipy_single"]["runtime"],
             "cold_total_s": r1["cold_sequence"]["total_s"],
             "accelerated_total_s": r1["warm_sequence"]["total_s"],
             "speedup": r1["speedup_sequence"],
             "prune_s": r1["prune"]["runtime"] if r1.get("prune") else None,
             "prune_dropped": r1["prune"]["n_dropped"] if r1.get("prune") else None,
             "accuracy_rel_diff": r1["warm_sequence"]["max_obj_diff_vs_cold"] / max(1e-12, abs(r1["accuracy"]["cold_obj"])),
             "verified_gap": r1["accuracy"]["scipy_dual_gap"]},
            {"instance": r2["instance"], "n_vars": r2["n_vars"],
             "scipy_one_target_s": None,
             "cold_total_s": r2["cold_highs"]["runtime"],
             "accelerated_total_s": r2["benders_1d"]["runtime_s"],
             "speedup": r2["benders_1d"]["speedup"],
             "prune_s": None, "prune_dropped": None,
             "accuracy_rel_diff": r2["benders_1d"]["rel_gap_vs_cold"],
             "verified_gap": r2["benders_1d"]["rel_gap_vs_cold"]}]
    df = pd.DataFrame(rows)
    df.to_csv(TAB / "E10_acceleration.csv", index=False)
    (LOG / "E10_acceleration_details.json").write_text(json.dumps([r1, r2], indent=1, default=str))
    print(df.round(5).to_string(index=False))


if __name__ == "__main__":
    main()
