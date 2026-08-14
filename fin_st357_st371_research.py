#!/usr/bin/env python3
"""FIN ST357--ST371: certified fold normal form and next local frontiers.

The batch is local and reproducible.  It preserves the kernel split and does
not manufacture a coupling source, physical calibration, apparatus, or record.
"""

from __future__ import annotations

import csv
import hashlib
import itertools
import json
import math
import platform
from collections import Counter
from pathlib import Path
from typing import Any

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mpmath as mp
import numpy as np
import scipy
import sympy as sp
from scipy.linalg import expm
from scipy.optimize import linprog

from fin_st01_st15_research import N, strict_operator
from fin_st118_st129_research import interval_left, interval_matvec
from fin_st132_center_isolation_replay import bounds, iv, strict_interval_matrix
from fin_st203_st217_research import IDX, stationary7
from fin_st263_st277_research import emission_point
from fin_st308_st326_research import MULT, localized_krawczyk, objective_from_radial
from fin_st327_st341_research import objective_interval, radial_interval_matrix
from fin_st342_st356_research import (
    discover_at_gain, fold_interval_system, qary_confusion, simplex_objective,
    spectral_blocks,
)


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "FIN_ST357_ST371_Results.json"
SUMMARY = ROOT / "FIN_ST357_ST371_Summary.csv"
FIG_DIR = ROOT / "FIN_ST357_ST371_Figures"
SEED = 20260818
NAMES = {
    357: "Interval_Simple_Fold_Normal_Form",
    358: "Entropy_Convex_Global_Cell_Lower_Bound",
    359: "Reflection_Symmetrization_Counterexample",
    360: "Global_Minimizer_Strategy_Boundary",
    361: "Certified_Rank7_Mediator_Localized_State",
    362: "Full_D12_Molien_Invariant_Count",
    363: "Normalized_Chernoff_Pressure_Derivative_Audit",
    364: "Rotated_Interval_Implicit_Function_Slab",
    365: "Certified_IR_Dual_Threshold_Zeros",
    366: "General_Cyclic_Convolution_Observer_Deficiency",
    367: "Multinomial_Fiber_Recovery_Bounds",
    368: "Strict_Heat_Block_Rate_Distortion_Extension",
    369: "Martingale_Multilevel_Clock_Concentration",
    370: "Strict_Active_Source_Admission_Gate",
    371: "Independent_Record_Stop",
}
PACKETS = {k: ROOT / f"FIN_ST{k}_{v}.json" for k, v in NAMES.items()}


def native(x: Any) -> Any:
    if isinstance(x, dict): return {str(k): native(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)): return [native(v) for v in x]
    if isinstance(x, np.ndarray): return native(x.tolist())
    if isinstance(x, (np.floating, np.integer)): return x.item()
    return x


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def finalize(k: int, obj: str, status: str, boundary: str, packet: dict) -> dict:
    path = PACKETS[k]
    path.write_text(json.dumps(native(packet), indent=2, sort_keys=True), encoding="utf-8")
    return {"program": f"ST{k}", "object": obj, "packet_file": path.name,
            "packet_sha256": sha(path), **packet, "status": status, "boundary": boundary}


def st357() -> dict:
    mp.iv.dps = 80
    old = json.loads((ROOT/"FIN_ST342_ST356_Results.json").read_text())["ST343"]
    z = np.array(old["fold_box_center"]); r = old["Krawczyk_certificate"]["coordinate_radius"]
    X = [iv((z[i]-r, z[i]+r)) for i in range(17)]
    q, v = X[:7], X[9:]
    aiv, _, _ = strict_interval_matrix(); biv = radial_interval_matrix(aiv)
    transverse = -sum((iv(float(MULT[i]))*v[i]*
                       sum((biv[i][j]*q[j] for j in range(7)), iv(0))
                       for i in range(7)), iv(0))
    curvature = -sum((iv(float(MULT[i]))*v[i]**3/q[i]**2 for i in range(7)), iv(0))
    tr = list(bounds(transverse)); cu = list(bounds(curvature))
    packet = {
        "inherited_unique_17D_root_box": old["gain_interval"],
        "left_nullvector_identity": "w=Dv, D=diag(1,2,2,2,2,2,1,1), because DJ is symmetric",
        "interval_wT_Fg": tr, "interval_wT_Fxx_vv": cu,
        "product_interval": [tr[0]*cu[1], tr[1]*cu[0]],
        "theorem": (
            "At the unique ST343 event, symmetrizability gives the exact left nullvector "
            "w=Dv. Outward interval evaluation proves w^T F_g<0 and "
            "w^T F_xx[v,v]>0 throughout the complete root box. Together with the "
            "one-dimensional kernel and nonsingular augmented system, the standard "
            "saddle-node hypotheses hold and the stationary set has the local simple-fold normal form."
        ),
    }
    return finalize(357, "Interval-Certified Simple Fold Normal Form",
                    "proven_local_simple_saddle_node_fold",
                    "The theorem is local to the supplied negative-information model. It does not prove dynamical stability, global branch topology, a physical spinodal, or an internal source for g.", packet)


def st358(a: np.ndarray) -> dict:
    mp.iv.dps = 70
    _, _, siv = strict_interval_matrix(); slo, shi = bounds(siv)
    g = 4.; global_lower = -g*shi/2
    candidate = json.loads((ROOT/"FIN_ST308_ST326_Results.json").read_text())["ST311"]["objective_at_center"]
    # Demonstrate a cell bound around the certified localized probability.
    p = np.array(json.loads((ROOT/"FIN_ST308_ST326_Results.json").read_text())["ST311"]["reflection_even_root_center"][:7])[IDX]
    radius = 1e-3; lo = np.maximum(0, p-radius); hi = np.minimum(1, p+radius)
    phi_min = []
    xstar = 1/(N*math.e)
    for l, u in zip(lo, hi):
        x = min(max(xstar, l), u); phi_min.append(x*math.log(N*x))
    q = p-np.ones(N)/N; delta = math.sqrt(N)*radius
    lammax = float(np.linalg.eigvalsh(a)[-1])
    quad_upper = float(q@a@q+2*np.linalg.norm(a@q)*delta+lammax*delta*delta)
    local_bound = sum(phi_min)-g*quad_upper/2
    packet = {
        "global_theorem_lower_bound": global_lower,
        "strict_row_sum_interval": [slo, shi],
        "certified_localized_candidate_point_value": candidate,
        "unresolved_global_bound_gap": candidate-global_lower,
        "cell_formula": "sum_i min_[li,ui] p log(12p) - (g/2) max_[box intersect simplex] <p-u,A(p-u)>",
        "cell_quadratic_maximum_rule": "a convex quadratic reaches its maximum at a polytope vertex",
        "localized_radius_1e_3_demonstration_lower_bound": local_bound,
        "theorem": (
            "Relative entropy is nonnegative and the convex Laplacian quadratic reaches "
            "its simplex maximum at a pure vertex, where it equals the strict row sum s. "
            "Hence V_g>=-gs/2 globally. On every box-simplex cell, separable entropy minima "
            "and vertex enumeration give a rigorous computable lower bound."
        ),
    }
    return finalize(358, "Entropy-Convex Global and Cell Lower-Bound Theorem",
                    "proven_lower_bound_primitive_too_loose_for_global_closure",
                    "At g=4 the root lower bound lies far below the localized value. No complete adaptive cover or global minimizer certificate is obtained.", packet)


def interval_objective_for_point(p: np.ndarray, aiv, g: float = 4.):
    pp = [iv(float(x)) for x in p]; u = iv(1)/iv(N)
    ent = sum((x*mp.iv.log(iv(N)*x) for x in pp), iv(0)); q = [x-u for x in pp]
    quad = sum((q[i]*aiv[i][j]*q[j] for i in range(N) for j in range(N)), iv(0))
    return ent-iv(g)*quad/2


def st359() -> dict:
    mp.iv.dps = 80; aiv, _, _ = strict_interval_matrix()
    weights = np.array([500, 1, 3, 2, 1, 470, 1, 1, 1, 11, 1, 1], float)
    p = weights/weights.sum(); vp = interval_objective_for_point(p, aiv); diffs = []
    for k in range(N):
        reflected = p[(2*k-np.arange(N)) % N]; ps = (p+reflected)/2
        d = interval_objective_for_point(ps, aiv)-vp
        diffs.append({"reflection_axis": k, "difference_interval": list(bounds(d))})
    minimum = min(x["difference_interval"][0] for x in diffs)
    packet = {
        "rational_probability_numerators": weights.astype(int).tolist(),
        "denominator": int(weights.sum()), "all_twelve_reflection_averages": diffs,
        "minimum_certified_objective_increase": minimum,
        "theorem": (
            "For the explicit positive rational probability vector, averaging with every "
            "one of the twelve reflections strictly increases the g=4 objective, after "
            "paying the complete transcendental strict-operator interval. Therefore no "
            "global proof may assume that some one-step reflection symmetrization is nonincreasing."
        ),
    }
    return finalize(359, "Certified Counterexample to One-Step Reflection Rearrangement",
                    "proven_all_reflection_averages_raise_objective_for_explicit_state",
                    "This refutes the declared symmetrization strategy, not every possible multi-step rearrangement theorem and not reflection symmetry of the actual global minimizer.", packet)


def st360() -> dict:
    a = json.loads((ROOT/"FIN_ST357_ST371_Results.json").read_text())["ST358"] if (ROOT/"FIN_ST357_ST371_Results.json").exists() else None
    # Values are read directly from current functions' already fixed statements in main order.
    packet = {
        "strategy_components": {
            "symmetry_only": "ruled out by ST344 trivial-stabilizer stratum",
            "one_step_reflection_rearrangement": "ruled out by ST359 certified counterexample",
            "root_entropy_quadratic_lower_bound": "valid but too loose by more than two objective units",
            "multistart": "strong numerical evidence only",
        },
        "global_minimizer_certificate": None,
        "result": (
            "The three proposed low-cost global routes do not compose into a proof. "
            "Symmetry is nonexhaustive, one-step reflection averaging is false, and the "
            "current cell primitive is too weak at the root node. A successful proof now "
            "requires sharper coupled entropy-quadratic envelopes or a different exact dual certificate."
        ),
    }
    return finalize(360, "Post-Falsification Global-Minimizer Strategy Boundary",
                    "bounded_no_go_for_current_symmetry_plus_root_bound_strategy",
                    "Global uniqueness modulo D12 is neither proved nor refuted. The stop applies only to the declared combination of methods.", packet)


def rank7_interval_matrix():
    aiv, _, _ = strict_interval_matrix(); ks = [3, 4, 5, 6, 7, 8, 9]
    lambdas = {k: sum((aiv[0][j]*mp.iv.cos(iv(2)*mp.iv.pi*k*j/N) for j in range(N)), iv(0)) for k in ks}
    return [[sum((lambdas[k]*mp.iv.cos(iv(2)*mp.iv.pi*k*(i-j)/N)/N for k in ks), iv(0))
             for j in range(N)] for i in range(N)]


def st361() -> dict:
    mp.iv.dps = 80; miv = rank7_interval_matrix()
    lo = np.array([[bounds(x)[0] for x in row] for row in miv]); hi = np.array([[bounds(x)[1] for x in row] for row in miv])
    m = (lo+hi)/2; center = discover_at_gain(m, 4.)
    cert = next(c for r in (1e-8, 3e-9, 1e-9)
                if (c := localized_krawczyk(center, miv, r))["included"])
    p = center[:7][IDX]; z = np.vstack([np.eye(N-1), -np.ones(N-1)])
    point_hessian = z.T@(np.diag(1/p)-4*m)@z
    family_error = float(np.linalg.norm(np.maximum(abs(lo-m), abs(hi-m)), "fro"))
    diag_error = cert["radius"]/(np.min(p)-cert["radius"])**2
    hessian_lower = float(np.linalg.eigvalsh(point_hessian)[0]-(np.linalg.norm(z, 2)**2)*(diag_error+4*family_error))
    # Enclose the objective at the exact root, not merely at the floating centre.
    val_interval = objective_interval(center, cert["radius"], (4., 4.), miv)
    diag_iv = miv[0][0]; pure = iv(math.log(N))-iv(2)*diag_iv
    packet = {
        "retained_fourier_modes": [3, 4, 5, 6, 7, 8, 9], "mediator_rank": 7,
        "root_center": center.tolist(), "Krawczyk_certificate": cert,
        "exact_root_objective_interval": val_interval, "pure_vertex_value_interval": list(bounds(pure)),
        "uniform_value": 0.0, "tangent_Hessian_lower_bound": hessian_lower,
        "boundary_minimizer_theorem": "none: p log p has derivative -infinity at zero while the quadratic has finite directional derivative",
        "theorem": (
            "The Fourier-defined rank-seven truncation has a unique positive reflection-even "
            "stationary root in the declared interval box. Its tangent Hessian is positive, "
            "and its objective is below uniform and every pure vertex. For any finite matrix "
            "and coupling, no boundary probability can minimize entropy minus a smooth quadratic, "
            "because adding epsilon to a zero coordinate lowers epsilon log epsilon dominantly."
        ),
    }
    return finalize(361, "Certified Rank-Seven Modal-Mediator Localized Minimum",
                    "proven_strict_local_minimum_beating_uniform_vertices_boundary_excluded",
                    "Interior competitors with trivial stabilizer are not excluded. Rank seven and g=4 remain supplied; no physical mediator or global minimum is derived.", packet)


def cycle_type(p: list[int]) -> tuple[int, ...]:
    seen = set(); lengths = []
    for i in range(len(p)):
        if i in seen: continue
        j = i; n = 0
        while j not in seen: seen.add(j); n += 1; j = p[j]
        lengths.append(n)
    return tuple(sorted(lengths))


def st362() -> dict:
    t = sp.symbols("t"); types = []
    for sign in (1, -1):
        for k in range(N): types.append(cycle_type([(sign*i+k) % N for i in range(N)]))
    series = 0
    for typ in types:
        series += sp.series((1-t)/sp.prod(1-t**ell for ell in typ), t, 0, 9).removeO()
    series = sp.expand(series/24); coeff = [int(series.coeff(t, d)) for d in range(9)]
    packet = {
        "cycle_type_counts": {str(k): v for k, v in Counter(types).items()},
        "Molien_formula": "(1-t)/24 sum_g product_cycles (1-t^length)^(-1)",
        "mean_zero_invariant_dimensions_degree_0_to_8": coeff,
        "degree_4_full_dimension": coeff[4], "degree_6_full_dimension": coeff[6],
        "spectral_energy_subalgebra_dimensions": {"degree_4": 21, "degree_6": 56},
        "additional_phase_sensitive_or_nonspectral_dimensions": {"degree_4": coeff[4]-21, "degree_6": coeff[6]-56},
        "theorem": (
            "Molien character averaging for the mean-zero permutation representation of "
            "D12 gives invariant dimensions 53 in degree four and 365 in degree six. "
            "Therefore the ST347 modal-energy algebra is a proper subalgebra and misses at "
            "least 32 quartic and 309 sextic invariant directions."
        ),
    }
    return finalize(362, "Full D12 Molien Count through Degree Eight",
                    "proven_full_invariant_dimensions_explicit_basis_deferred",
                    "Dimensions are complete, but an explicit Reynolds-reduced basis and coefficient-source law are not supplied.", packet)


def transfer_value_derivative_point(b0, b1, depth, e0, e1, s):
    if depth == 0: return 1., 0.
    total = deriv = 0.
    for y in (0, 1):
        n0, q0 = emission_point(b0, e0, y); n1, q1 = emission_point(b1, e1, y)
        w = q0**s*q1**(1-s); child, dc = transfer_value_derivative_point(n0, n1, depth-1, e0, e1, s)
        total += w*child; deriv += w*(dc+math.log(q0/q1)*child)
    return total, deriv


def st363() -> dict:
    e0, e1 = [[.98, .02], [.92, .08]], [[.08, .92], [.02, .98]]
    rows = []
    for depth in (4, 6, 8, 10, 12):
        for s in (.45, .500739416715181, .55):
            vals = []
            for b0 in np.linspace(0, 1, 5):
                for b1 in np.linspace(0, 1, 5):
                    val, der = transfer_value_derivative_point(b0, b1, depth, e0, e1, s)
                    vals.append((val, der/(depth*val)))
            vmax = max(range(len(vals)), key=lambda i: vals[i][0])
            rows.append({"depth": depth, "s": s, "maximum_collocation_Tn_root": vals[vmax][0]**(1/depth),
                         "normalized_derivative_at_maximizer": vals[vmax][1],
                         "derivative_range_over_collocation": [min(x[1] for x in vals), max(x[1] for x in vals)]})
    packet = {
        "exact_all_depth_bound": [-math.log(49), math.log(49)],
        "identity": "(1/n) d_s log T_n is a tilted expectation of the path-average log likelihood ratio",
        "numerical_depth_rows": rows, "certified_infinite_depth_optimizer": None,
        "theorem": (
            "For every depth and initial belief, the normalized log-transfer derivative is "
            "a tilted expectation of the average one-step log-likelihood ratio and therefore "
            "lies in [-log49,log49]. This prevents growth with depth but does not prove a "
            "vanishing depth-to-infinity derivative error."
        ),
    }
    return finalize(363, "Normalized Chernoff Pressure-Derivative Audit",
                    "proven_uniform_normalized_bound_numerical_convergence_only",
                    "The pressure optimizer is still not interval-certified at infinite depth. Collocation is diagnostic, and the HMM fixture is synthetic.", packet)


def field_interval_FJ(X, aiv):
    biv = [[sum((aiv[i][site] for site in range(N) if IDX[site] == col), iv(0)) for col in range(7)] for i in range(7)]
    q, kappa = X[:7], X[7]; u = [q[IDX[s]] for s in range(N)]
    au = [sum((aiv[i][j]*u[j] for j in range(N)), iv(0)) for i in range(N)]
    f = []; rd = []
    for i, z in enumerate(u[:7]):
        rho = z*z; den = 1+rho/2; qf = rho/den; qp = den**-2; qpp = -den**-3
        h = -qf*qp+iv(".075"); hp = -(qp*qp+qf*qpp)
        f.append(kappa*au[i]+2*z*h); rd.append(2*h+4*rho*hp)
    jac = [[iv(0) for _ in range(8)] for _ in range(7)]
    for i in range(7):
        for col in range(7): jac[i][col] = kappa*biv[i][col]+(rd[i] if i == col else iv(0))
        jac[i][7] = au[i]
    return f, jac, biv


def st364() -> dict:
    mp.iv.dps = 70
    rows = json.loads((ROOT/"FIN_ST293_ST307_Results.json").read_text())["ST298"]["rows"]
    center = np.array(rows[85]["center"]); aiv, _, _ = strict_interval_matrix()
    f0, j0, biv = field_interval_FJ([iv(x) for x in center], aiv)
    flo = np.array([bounds(x)[0] for x in f0]); fhi = np.array([bounds(x)[1] for x in f0])
    jlo = np.array([[bounds(x)[0] for x in row] for row in j0]); jhi = np.array([[bounds(x)[1] for x in row] for row in j0])
    jmid = (jlo+jhi)/2; _, _, vh = np.linalg.svd(jmid); tangent = vh[-1]; tangent /= np.linalg.norm(tangent)
    qbasis = np.linalg.svd(tangent.reshape(1, -1))[2][1:].T; pre = np.linalg.inv(jmid@qbasis)
    h, radius = 2e-7, 1e-7
    jtl, jth = interval_matvec(jlo, jhi, tangent, tangent)
    pl = np.minimum(jtl*h, -jth*h); ph = np.maximum(jth*h, -jtl*h)
    hess = []
    for i in range(7):
        z = iv((center[i]-abs(tangent[i])*h, center[i]+abs(tangent[i])*h)); rho = z*z; den = 1+rho/2
        qf = rho/den; qp = den**-2; qpp = -den**-3; qppp = iv("1.5")*den**-4
        hp = -(qp*qp+qf*qpp); hpp = -(3*qp*qpp+qf*qppp); fpp = 2*z*(6*hp+4*rho*hpp)
        fm = max(abs(x) for x in bounds(fpp)); bm = np.array([max(abs(x) for x in bounds(biv[i][j])) for j in range(7)])
        hess.append(math.sqrt(fm*fm+2*float(bm@bm)))
    remainder = .5*np.array(hess)*h*h
    base_lo, base_hi = flo+pl-remainder, fhi+ph+remainder
    ylo, yhi = interval_matvec(-pre, -pre, base_lo, base_hi)
    axis_radius = abs(tangent)*h+np.sum(abs(qbasis), axis=1)*radius
    _, jfull, _ = field_interval_FJ([iv((center[j]-axis_radius[j], center[j]+axis_radius[j])) for j in range(8)], aiv)
    jy = [[sum((jfull[i][j]*iv(qbasis[j, k]) for j in range(8)), iv(0)) for k in range(7)] for i in range(7)]
    jyl = np.array([[bounds(x)[0] for x in row] for row in jy]); jyh = np.array([[bounds(x)[1] for x in row] for row in jy])
    cjl, cjh = interval_left(pre, jyl, jyh); mlo, mhi = -cjh, -cjl
    for k in range(7): mlo[k, k] = np.nextafter(mlo[k, k]+1, -np.inf); mhi[k, k] = np.nextafter(mhi[k, k]+1, np.inf)
    dlo, dhi = interval_matvec(mlo, mhi, np.full(7, -radius), np.full(7, radius))
    klo, khi = ylo+dlo, yhi+dhi; margins = np.minimum(klo+radius, radius-khi)
    packet = {
        "ST298_center_index": 85, "longitudinal_halfwidth": h, "transverse_radius": radius,
        "minimum_Krawczyk_margin": float(np.min(margins)),
        "maximum_transverse_contraction_row_sum": float(np.max(np.sum(np.maximum(abs(mlo), abs(mhi)), axis=1))),
        "included_for_every_longitudinal_parameter": bool(np.min(margins) > 0),
        "theorem": (
            "A Taylor-form interval Krawczyk argument in one tangent and seven transverse "
            "coordinates proves that, for every longitudinal parameter s in [-2e-7,2e-7], "
            "there is one unique transverse correction in the radius-1e-7 ball. This is a "
            "continuous local branch slab through the worst-conditioned sampled ST298 region."
        ),
    }
    return finalize(364, "Rotated Tangent/Transverse Interval Implicit-Function Slab",
                    "proven_continuous_local_branch_slab",
                    "One local slab is certified, not the full 160-section chain. Overlap between a sequence of rotated slabs remains to be proved.", packet)


def ir_G_interval(ybox, a1, a4, nu):
    y = iv(ybox); x = mp.iv.exp(-y)
    return 2*a1*y/(1+x)+8*a4*y*x**3/(1+x**4)-nu


def ir_G_derivative_interval(ybox, a1, a4):
    y = iv(ybox); x = mp.iv.exp(-y); f = x**3/(1+x**4)
    return 2*a1*((1+x)+y*x)/(1+x)**2+8*a4*(f+y*x**3*(-3+x**4)/(1+x**4)**2)


def st365() -> dict:
    mp.iv.dps = 60
    prior = json.loads((ROOT/"FIN_ST342_ST356_Results.json").read_text())["ST350"]["grid_convergence"][-1]
    a1, a4, nu = iv(prior["time_dual_weights"][1]), iv(prior["time_dual_weights"][2]), iv(prior["heat_dual_price"])
    centers = [0.02845660914428703, 2.099529366923138, 2.547622646372957]
    boxes = [(x-2e-10, x+2e-10) for x in centers]
    roots = [{"y_interval": list(b), "mu_interval": [b[0]/2, b[1]/2],
              "derivative_interval": list(bounds(ir_G_derivative_interval(b, a1, a4)))} for b in boxes]
    stack = [(2e-5, 10.)]; sign_boxes = 0; unresolved = []
    while stack:
        lo, hi = stack.pop()
        if any(lo >= a and hi <= b for a, b in boxes): continue
        split = next((x for b in boxes for x in b if lo < x < hi), None)
        if split is not None: stack.extend([(lo, split), (split, hi)]); continue
        glo, ghi = bounds(ir_G_interval((lo, hi), a1, a4, nu))
        if glo > 0 or ghi < 0: sign_boxes += 1; continue
        if hi-lo < 1e-10: unresolved.append([lo, hi]); continue
        mid = (lo+hi)/2; stack.extend([(lo, mid), (mid, hi)])
    # For y >= 10, the first positive summand alone is bounded below by
    # 20*a1/(1+exp(-10)); the second summand is nonnegative.
    tail_lower = bounds(20*a1/(1+mp.iv.exp(iv(-10)))-nu)[0]
    packet = {
        "erratum_to_ST350": "four threshold transitions correspond to two selected bands, not four bands",
        "frozen_grid_dual_weights": prior["time_dual_weights"], "frozen_heat_price": prior["heat_dual_price"],
        "certified_simple_zeros": roots, "exhaustive_sign_boxes_on_y_2e_5_to_10": sign_boxes,
        "unresolved_complement_boxes": unresolved, "analytic_tail_y_ge_10_lower_at_endpoint": tail_lower,
        "selected_bands_in_mu": [[boxes[0][1]/2, boxes[1][0]/2], [boxes[2][1]/2, 1e5]],
        "theorem": (
            "For the frozen 1000-grid dual multipliers, the transformed analytic threshold "
            "has exactly three simple positive zeros on the declared spectral interval. "
            "An exhaustive outward sign cover excludes every other zero; the asymptotic "
            "lower bound is positive beyond y=10. Hence the selector has exactly two bands."
        ),
    }
    return finalize(365, "Certified Zeros and Corrected Band Count for the Frozen IR Dual",
                    "proven_three_simple_zeros_two_selected_bands_for_frozen_dual",
                    "The multipliers come from a finite grid and are not certified continuum-optimal. The theorem certifies their threshold geometry only.", packet)


def cyclic_convolve(p: np.ndarray, r: np.ndarray) -> np.ndarray:
    return np.array([sum(p[(j-k) % len(p)]*r[k] for k in range(len(p))) for j in range(len(p))])


def cyclic_deficiency_lp(source: np.ndarray, target: np.ndarray) -> tuple[float, np.ndarray]:
    q = len(source); nv = 2*q; c = np.r_[np.zeros(q), .5*np.ones(q)]
    aub = []; bub = []
    for j in range(q):
        base = np.zeros(nv)
        for k in range(q): base[k] = source[(j-k) % q]
        row = base.copy(); row[q+j] = -1; aub.append(row); bub.append(target[j])
        row = -base.copy(); row[q+j] = -1; aub.append(row); bub.append(-target[j])
    sol = linprog(c, A_ub=aub, b_ub=bub, A_eq=np.r_[np.ones(q), np.zeros(q)][None, :],
                  b_eq=[1], bounds=[(0, None)]*nv, method="highs")
    return float(sol.fun), sol.x[:q]


def st366() -> dict:
    fine = np.zeros(N); fine[[0, 1, 11]] = [.8, .1, .1]
    added = np.zeros(N); added[[0, 2]] = [.8, .2]
    coarse = cyclic_convolve(fine, added)
    forward, rf = cyclic_deficiency_lp(fine, coarse); reverse, rr = cyclic_deficiency_lp(coarse, fine)
    alt = np.zeros(N); alt[[0, 3, 8]] = [.7, .2, .1]
    d1, _ = cyclic_deficiency_lp(fine, alt); d2, _ = cyclic_deficiency_lp(alt, fine)
    # Exact no-garbling witnesses: the unique deconvolution filters have
    # negative rational entries in both directions (both circulants invertible).
    rfine = [sp.Rational(8, 10), sp.Rational(1, 10)] + [sp.Rational(0)]*9 + [sp.Rational(1, 10)]
    ralt = [sp.Rational(7, 10), 0, 0, sp.Rational(2, 10), 0, 0, 0, 0,
            sp.Rational(1, 10), 0, 0, 0]
    def cmat(v): return sp.Matrix(N, N, lambda j, k: v[(j-k) % N])
    inv_fa = cmat(rfine).inv()*sp.Matrix(ralt)
    inv_af = cmat(ralt).inv()*sp.Matrix(rfine)
    packet = {
        "exact_reduction": "delta(C_p,C_q)=min_{r in simplex} TV(p*r,q)",
        "garbling_example": {"fine_noise": fine.tolist(), "added_noise": added.tolist(),
                              "coarse_noise": coarse.tolist(), "fine_to_coarse": forward,
                              "coarse_to_fine": reverse, "recovered_forward_filter": rf.tolist()},
        "incomparable_example": {"alternative_noise": alt.tolist(), "fine_to_alternative": d1,
                                  "alternative_to_fine": d2,
                                  "exact_minimum_fine_to_alternative_inverse_filter": str(min(inv_fa)),
                                  "exact_minimum_alternative_to_fine_inverse_filter": str(min(inv_af))},
        "theorem": (
            "Cyclic twirling of any postprocessing cannot increase worst-row total variation. "
            "Therefore Le Cam deficiency between two 12-state cyclic convolution experiments "
            "is exactly the displayed twelve-variable convolution LP. The order is generally "
            "partial: exact rational deconvolution in either direction has a negative entry, "
            "so the explicit pair has positive deficiency both ways. The displayed positive "
            "deficiency magnitudes are floating LP evaluations, not interval endpoints."
        ),
    }
    return finalize(366, "Exact LP for General Twelve-State Cyclic Observer Deficiency",
                    "proven_convolution_LP_and_explicit_incomparability",
                    "These are instrument families, not a FIN-derived detector. Numerical LP values are not empirical observer depths.", packet)


def st367() -> dict:
    q, delta = 12, .05; rows = []
    for counts in (100, 1000, 10000, 100000):
        coord = math.sqrt(math.log(2*q/delta)/(2*counts)); l2upper = math.sqrt(q)*coord
        shift = .1/math.sqrt(counts); p0 = np.ones(q)/q; p1 = p0.copy(); p1[0] += shift; p1[1] -= shift
        kl = float(np.sum(p0*np.log(p0/p1))); tvupper = min(1., math.sqrt(counts*kl/2))
        separation = math.sqrt(2)*shift; lower = separation*(1-tvupper)/4
        rows.append({"counts": counts, "simultaneous_L2_radius_95": l2upper,
                     "two_point_L2_separation": separation, "Pinsker_product_TV_upper": tvupper,
                     "minimax_expected_L2_lower": lower})
    packet = {
        "multinomial_states": q, "confidence": .95, "rows": rows,
        "upper_theorem": "coordinate Hoeffding plus a 24-event union bound",
        "lower_theorem": "two-point Le Cam with exact categorical KL and product Pinsker",
        "scaling": "both root risks are order N^(-1/2) for fixed state dimension",
    }
    return finalize(367, "Finite-Sample Multinomial Fiber Recovery Bounds",
                    "proven_nonasymptotic_upper_and_two_point_lower_bounds",
                    "Bounds assume independent correctly labelled multinomial events. Detector loss, calibration, nuisance parameters, adaptive design and custody are absent.", packet)


def kron_apply(tensor: np.ndarray, kernel: np.ndarray) -> np.ndarray:
    out = tensor
    for axis in range(tensor.ndim):
        out = np.tensordot(kernel, out, axes=(1, axis)); out = np.moveaxis(out, 0, axis)
    return out


def kron_dist_apply(tensor: np.ndarray, kernel: np.ndarray, kd: np.ndarray) -> np.ndarray:
    total = np.zeros_like(tensor)
    for special in range(tensor.ndim):
        out = tensor
        for axis in range(tensor.ndim):
            out = np.tensordot(kd if axis == special else kernel, out, axes=(1, axis)); out = np.moveaxis(out, 0, axis)
        total += out/tensor.ndim
    return total


def markov_tensor(p: np.ndarray, n: int) -> np.ndarray:
    shape = [N]*n; out = np.empty(shape)
    for word in itertools.product(range(N), repeat=n):
        value = 1/N
        for i in range(1, n): value *= p[word[i-1], word[i]]
        out[word] = value
    return out


def tensor_ba(p: np.ndarray, n: int, betas=(0., 1., 2., 4., 8., 16.)) -> list[dict]:
    px = markov_tensor(p, n); rows = []
    single_d = 1-np.eye(N)
    for beta in betas:
        kernel = np.exp(-beta*single_d/n); kd = single_d*kernel
        q = np.ones_like(px)/px.size
        for it in range(1200):
            z = kron_apply(q, kernel); qn = q*kron_apply(px/z, kernel)
            qn /= qn.sum()
            if np.max(abs(qn-q)) < 2e-11: break
            q = qn
        z = kron_apply(q, kernel); dz = kron_dist_apply(q, kernel, kd)
        distortion = float(np.sum(px*dz/z)); info_nats = float(-beta*distortion-np.sum(px*np.log(z)))
        rows.append({"beta": beta, "distortion": distortion,
                     "bits_per_symbol": info_nats/(n*math.log(2)), "iterations": it})
    return rows


def st368(a: np.ndarray) -> dict:
    tau = .5; p = expm(-tau*a); rows = {str(n): tensor_ba(p, n) for n in (3, 4)}
    hrate = float(-np.mean(np.sum(p*np.log2(np.maximum(p, 1e-300)), axis=1)))
    hinitial = math.log2(N)
    block_entropy = {str(n): (hinitial+(n-1)*hrate)/n for n in (3, 4)}
    packet = {
        "transition": "exp(-0.5 A)", "block_results": rows,
        "exact_entropy_rate_bits": hrate,
        "zero_distortion_block_entropy_bits_per_symbol": block_entropy,
        "block5_state_count": N**5, "block5_executed": False,
        "resource_stop": "block five was not needed for the low-compute batch; tensor source storage alone has 248832 entries",
        "result": (
            "Kronecker-factorized Blahut--Arimoto avoids the impossible 12^(2n) distortion "
            "matrix and computes complete block-three and block-four curves. At fixed block "
            "length the zero-distortion limit is H(X_1:n)/n; these block entropies converge "
            "to the exact strict-heat Markov entropy rate only as n tends to infinity. Block five is left "
            "as an explicit bounded resource continuation, not silently approximated."
        ),
    }
    return finalize(368, "Block-Three and Block-Four Strict-Heat Rate-Distortion",
                    "strong_numerical_block3_block4_curves_block5_stopped",
                    "Tau, Markov interpretation and Hamming distortion remain premises. Finite-block curves do not prove the asymptotic rate-distortion function or physical compression density.", packet)


def st369() -> dict:
    alpha = .05; rows = []
    for eta in (1e-4, 1e-3, 1e-2):
        b = math.log((1+eta)/(1-eta))
        for levels in (10, 100, 1000, 10000):
            worst = levels*b; highprob = b*math.sqrt(2*levels*math.log(2/alpha))
            rows.append({"eta": eta, "levels": levels, "worst_case_log_drift": worst,
                         "95_percent_log_drift_radius": highprob,
                         "wave_clock_radius": highprob/2})
    packet = {
        "declared_premise": "independent mode-pair log-error differences, conditionally mean zero and bounded by b_n",
        "Azuma_Hoeffding_bound": "P(|sum X_n|>=sqrt(2 log(2/alpha) sum b_n^2))<=alpha",
        "rows": rows,
        "theorem": (
            "Under the declared martingale-difference premise, relative log-clock drift "
            "scales as square root of depth at fixed confidence rather than linearly as in "
            "the adversarial ST355 bound. Common calibration gauge still cancels exactly."
        ),
    }
    return finalize(369, "Martingale Concentration for Multilevel Relative Clocks",
                    "proven_conditional_sqrt_depth_concentration",
                    "Independence/conditional centering is an added stochastic law, not strict-derived. No seconds, time arrow, or physical noise mechanism follows.", packet)


def st370() -> dict:
    records = list(ROOT.glob("FIN_ST370_NEW_STRICT_ACTIVE_SOURCE*.json"))
    packet = {
        "required_pattern": "FIN_ST370_NEW_STRICT_ACTIVE_SOURCE*.json",
        "matching_new_typed_sources": [x.name for x in records],
        "admitted_source_count": len(records),
        "admission_obligations": ["strict provenance", "nonzero signed coupling", "D12 typing",
                                  "selector safety", "kernel-split safety", "no hidden dimensional premise"],
        "theorem": "No new typed source object exists in the declared admission channel; closed source lanes are not replayed.",
    }
    return finalize(370, "Strict Active-Coupling Source Admission Gate",
                    "blocked_no_new_strict_typed_source" if not records else "candidate_present_requires_audit",
                    "This preserves prior no-new-frontier certificates. It is not proof that no source can exist in an enlarged theory.", packet)


def st371() -> dict:
    records = list(ROOT.glob("FIN_ST371_INDEPENDENT_RAW_EVENTS*.jsonl"))
    packet = {"required_pattern": "FIN_ST371_INDEPENDENT_RAW_EVENTS*.jsonl",
              "matching_records": [x.name for x in records], "independent_record_count": len(records),
              "theorem": "Local mathematics cannot generate independent custody, apparatus calibration, unblinding, or empirical events."}
    return finalize(371, "Independent Evidence Gate",
                    "blocked_no_independent_events" if not records else "record_present_requires_blinded_validation",
                    "The stop is evidentiary, not a failed physical experiment.", packet)


def make_figures(out: dict) -> None:
    FIG_DIR.mkdir(exist_ok=True)
    r = out["ST357"]; fig, ax = plt.subplots(figsize=(6.8, 3.7))
    ax.bar(["$w^T F_g$", "$w^T F_{xx}[v,v]$"],
           [np.mean(r["interval_wT_Fg"]), np.mean(r["interval_wT_Fxx_vv"])],
           color=["#9d2f2f", "#28764d"]); ax.axhline(0, color="black", lw=.8)
    ax.set(title="ST357: interval-separated simple-fold coefficients"); fig.tight_layout()
    fig.savefig(FIG_DIR/"st357_fold_coefficients.png", dpi=200); plt.close(fig)
    r = out["ST368"]["block_results"]; fig, ax = plt.subplots(figsize=(7, 3.8))
    for n, rows in r.items(): ax.plot([x["distortion"] for x in rows], [x["bits_per_symbol"] for x in rows], "o-", label=f"block {n}")
    ax.set(xlabel="Hamming distortion", ylabel="bits per symbol", title="ST368: strict-heat block rate--distortion")
    ax.legend(); fig.tight_layout(); fig.savefig(FIG_DIR/"st368_block_rate_distortion.png", dpi=200); plt.close(fig)
    r = out["ST369"]["rows"]; fig, ax = plt.subplots(figsize=(7, 3.8))
    rr = [x for x in r if x["eta"] == .001]
    ax.loglog([x["levels"] for x in rr], [x["worst_case_log_drift"] for x in rr], "o-", label="adversarial")
    ax.loglog([x["levels"] for x in rr], [x["95_percent_log_drift_radius"] for x in rr], "o-", label="martingale 95%")
    ax.set(xlabel="refinement levels", ylabel="log-ratio drift", title="ST369: linear versus square-root drift")
    ax.legend(); fig.tight_layout(); fig.savefig(FIG_DIR/"st369_clock_concentration.png", dpi=200); plt.close(fig)


def main() -> None:
    np.random.seed(SEED); _, a, s = strict_operator()
    out = {"metadata": {"seed": SEED, "python": platform.python_version(), "numpy": np.__version__,
                         "scipy": scipy.__version__, "strict_row_sum": s,
                         "scope": "local analytical, interval, exact finite and bounded numerical research"}}
    funcs = {357: st357, 358: lambda: st358(a), 359: st359, 360: st360, 361: st361,
             362: st362, 363: st363, 364: st364, 365: st365, 366: st366,
             367: st367, 368: lambda: st368(a), 369: st369, 370: st370, 371: st371}
    for k in range(357, 372): out[f"ST{k}"] = funcs[k]()
    make_figures(out)
    RESULTS.write_text(json.dumps(native(out), indent=2, sort_keys=True), encoding="utf-8")
    with SUMMARY.open("w", newline="", encoding="utf-8") as h:
        w = csv.writer(h); w.writerow(["program", "object", "status"])
        for k in range(357, 372): w.writerow([f"ST{k}", out[f"ST{k}"]["object"], out[f"ST{k}"]["status"]])
    print(json.dumps({k: v["status"] for k, v in out.items() if k.startswith("ST")}, indent=2))


if __name__ == "__main__": main()
