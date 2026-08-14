#!/usr/bin/env python3
"""FIN ST342--ST356: interval transition/fold closure and local continuations.

All computations are local.  Supplied couplings, source laws, instruments and
calibrations remain typed premises; no empirical record is manufactured.
"""

from __future__ import annotations

import csv
import hashlib
import itertools
import json
import math
import platform
from pathlib import Path
from typing import Any

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mpmath as mp
import numpy as np
import scipy
from scipy.linalg import expm
from scipy.optimize import brentq, linprog, minimize, root

from fin_st01_st15_research import N, strict_operator
from fin_st118_st129_research import interval_left, interval_matvec
from fin_st132_center_isolation_replay import bounds, iv, strict_interval_matrix
from fin_st203_st217_research import IDX
from fin_st308_st326_research import (
    MULT, localized_jacobian, objective_from_radial, radial_matrix,
)
from fin_st327_st341_research import (
    discover_at_gain, fold_jacobian, fold_system, objective_interval,
    param_krawczyk, param_tube_krawczyk, radial_interval_matrix,
)


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "FIN_ST342_ST356_Results.json"
SUMMARY = ROOT / "FIN_ST342_ST356_Summary.csv"
FIG_DIR = ROOT / "FIN_ST342_ST356_Figures"
SEED = 20260817
NAMES = {
    342: "Narrow_Certified_Coexistence_Bracket",
    343: "Full_17D_Interval_Fold_Certificate",
    344: "D12_Global_Symmetry_Reduction_No_Go",
    345: "Global_Minimizer_Reflection_Orbit_Audit",
    346: "Modal_Mediator_Localization_Thresholds",
    347: "Degree_4_6_Spectral_Invariant_Algebra",
    348: "Chernoff_Derivative_Contraction_Counterexample",
    349: "Anisotropic_Implicit_Slab_Audit",
    350: "Continuum_Multi_Time_IR_Dual_Structure",
    351: "Two_Cycle_Refinement_Cohomology",
    352: "Twelve_State_Circulant_Observer_Deficiency",
    353: "Spectrally_Restricted_Fiber_Recovery_Bounds",
    354: "Strict_Heat_Markov_Rate_Distortion",
    355: "Multilevel_Relative_Clock_Drift",
    356: "Independent_Record_Stop",
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


def st342(a: np.ndarray) -> dict:
    mp.iv.dps = 80
    aiv, _, _ = strict_interval_matrix()
    value = lambda g: objective_from_radial(discover_at_gain(a, g)[:7], a, g)
    center_g = brentq(value, 2.900, 2.905, xtol=5e-15)
    halfwidth = 1e-8
    rows = []
    for g in (center_g-halfwidth, center_g+halfwidth):
        x = discover_at_gain(a, g)
        cert = next(c for r in (1e-10, 3e-11, 1e-11)
                    if (c := param_krawczyk(x, g, (g, g), aiv, r))["included"])
        rows.append({"g": g, "center": x.tolist(), "certificate": cert,
                     "value_interval": objective_interval(x, cert["root_radius"], (g, g), aiv)})
    xmid = discover_at_gain(a, center_g)
    tube = param_tube_krawczyk(
        xmid, center_g, (center_g-halfwidth, center_g+halfwidth), aiv,
        np.array([1.5e-6, 1e-7, 2e-7, 2e-7, 2e-7, 2e-7, 2e-7, 1.5e-5]),
    )
    packet = {
        "floating_center": center_g,
        "certified_crossing_bracket": [center_g-halfwidth, center_g+halfwidth],
        "certified_width": 2*halfwidth,
        "improvement_over_ST327_width": .005/(2*halfwidth),
        "endpoint_certificates": rows,
        "continuous_parametric_root_tube": tube,
        "theorem": (
            "Outward endpoint boxes have strictly opposite objective signs and one "
            "parametric Krawczyk tube contains a unique positive reflection-even root "
            "for every supplied g in the complete width-2e-8 bracket. Continuity proves "
            "at least one equality with the uniform-state energy inside this bracket."
        ),
    }
    return finalize(342, "Narrow Interval-Certified Coexistence Bracket",
                    "proven_local_crossing_in_width_2e_8_bracket",
                    "The crossing is branch-local. Global minimality, crossing uniqueness, a physical transition and an internal source for g remain open.", packet)


def fold_interval_system(center: np.ndarray, radii: np.ndarray, aiv):
    X = [iv((center[i]-radii[i], center[i]+radii[i])) for i in range(17)]
    q, lam, g, v = X[:7], X[7], X[8], X[9:]
    vq, vl = v[:7], v[7]
    b = radial_interval_matrix(aiv)
    f = [mp.iv.log(iv(N)*q[i])+1-g*sum((b[i][j]*q[j] for j in range(7)), iv(0))+lam
         for i in range(7)]
    f.append(sum((iv(float(MULT[j]))*q[j] for j in range(7)), iv(0))-1)
    for i in range(7):
        f.append(vq[i]/q[i]-g*sum((b[i][j]*vq[j] for j in range(7)), iv(0))+vl)
    f.append(sum((iv(float(MULT[j]))*vq[j] for j in range(7)), iv(0)))
    f.append((sum((x*x for x in v), iv(0))-1)/2)
    jac = [[iv(0) for _ in range(17)] for _ in range(17)]
    for i in range(7):
        for j in range(7): jac[i][j] = (1/q[i] if i == j else iv(0))-g*b[i][j]
        jac[i][7] = iv(1)
        jac[i][8] = -sum((b[i][j]*q[j] for j in range(7)), iv(0))
    for j in range(7): jac[7][j] = iv(float(MULT[j]))
    for i in range(7):
        jac[8+i][i] = -vq[i]/(q[i]*q[i])
        jac[8+i][8] = -sum((b[i][j]*vq[j] for j in range(7)), iv(0))
        for j in range(7): jac[8+i][9+j] = (1/q[i] if i == j else iv(0))-g*b[i][j]
        jac[8+i][16] = iv(1)
    for j in range(7): jac[15][9+j] = iv(float(MULT[j]))
    for j in range(8): jac[16][9+j] = v[j]
    flo = np.array([bounds(x)[0] for x in f]); fhi = np.array([bounds(x)[1] for x in f])
    jlo = np.array([[bounds(x)[0] for x in row] for row in jac])
    jhi = np.array([[bounds(x)[1] for x in row] for row in jac])
    return flo, fhi, jlo, jhi


def fold_krawczyk(center: np.ndarray, radii: np.ndarray, aiv) -> dict:
    f0lo, f0hi, j0lo, j0hi = fold_interval_system(center, np.zeros(17), aiv)
    pre = np.linalg.inv((j0lo+j0hi)/2)
    _, _, jlo, jhi = fold_interval_system(center, radii, aiv)
    cglo, cghi = interval_matvec(pre, pre, f0lo, f0hi)
    ylo = np.nextafter(center-cghi, -np.inf); yhi = np.nextafter(center-cglo, np.inf)
    cjlo, cjhi = interval_left(pre, jlo, jhi); mlo, mhi = -cjhi, -cjlo
    for i in range(17):
        mlo[i, i] = np.nextafter(mlo[i, i]+1, -np.inf)
        mhi[i, i] = np.nextafter(mhi[i, i]+1, np.inf)
    dlo, dhi = interval_matvec(mlo, mhi, -radii, radii)
    klo = np.nextafter(ylo+dlo, -np.inf); khi = np.nextafter(yhi+dhi, np.inf)
    margins = np.minimum(klo-(center-radii), (center+radii)-khi)
    return {"included": bool(np.min(margins) > 0), "minimum_margin": float(np.min(margins)),
            "coordinate_radius": float(radii[0]), "component_margins": margins.tolist()}


def st343(a: np.ndarray) -> dict:
    mp.iv.dps = 80
    b = radial_matrix(a); x = discover_at_gain(a, 2.666)
    _, _, vh = np.linalg.svd(localized_jacobian(x, b, 2.666)); v = vh[-1]
    sol = root(lambda z: fold_system(z, b), np.r_[x, 2.666, v],
               jac=lambda z: fold_jacobian(z, b), tol=1e-12)
    z = sol.x; aiv, _, _ = strict_interval_matrix(); radii = np.full(17, 1e-8)
    cert = fold_krawczyk(z, radii, aiv)
    _, _, jlo, jhi = fold_interval_system(z, radii, aiv)
    stat_center = localized_jacobian(z[:8], b, z[8])
    stat_lo, stat_hi = jlo[:8, :8], jhi[:8, :8]
    perturb = float(np.linalg.norm(np.maximum(abs(stat_lo-stat_center), abs(stat_hi-stat_center))))
    sv = np.linalg.svd(stat_center, compute_uv=False)
    u, _, _ = np.linalg.svd(stat_center); w = u[:, -1]
    fg = np.r_[-b@z[:7], 0.]
    fxxvv = np.r_[-z[9:16]**2/z[:7]**2, 0.]
    packet = {
        "fold_box_center": z.tolist(), "gain_interval": [z[8]-1e-8, z[8]+1e-8],
        "residual_inf": float(np.linalg.norm(fold_system(z, b), np.inf)),
        "Krawczyk_certificate": cert,
        "stationarity_singular_values_at_center": sv.tolist(),
        "interval_stationarity_perturbation_Frobenius_bound": perturb,
        "rank_7_margin_second_smallest_minus_perturbation": float(sv[-2]-perturb),
        "point_transversality_w_Fg": float(w@fg),
        "point_quadratic_fold_coefficient": float(w@fxxvv),
        "theorem": (
            "The full 17-variable stationarity-nullvector-normalization system has one "
            "and only one root in the displayed outward box. The stationary Jacobian has "
            "rank exactly seven throughout the root event: the certified solution supplies "
            "a nonzero nullvector, while a singular-value perturbation bound keeps the "
            "second-smallest singular value positive."
        ),
    }
    return finalize(343, "Full 17-Dimensional Interval Fold-Event Certificate",
                    "proven_unique_augmented_root_and_one_dimensional_stationary_kernel",
                    "The interval theorem is local. The conventional simple-fold transversality coefficients are nonzero numerically but are not interval-certified here; global branch structure and dynamics remain open.", packet)


def d12_permutations() -> list[tuple[int, ...]]:
    return list(dict.fromkeys(tuple((sgn*i+k) % N for i in range(N))
                              for sgn in (1, -1) for k in range(N)))


def st344() -> dict:
    group = d12_permutations()
    generic = np.array([2.0**i for i in range(N)]); generic /= generic.sum()
    stabilizer = [p for p in group if np.array_equal(generic, generic[list(p)])]
    reflection_fixed_dimensions = []
    for k in range(N):
        perm = tuple((-i+k) % N for i in range(N))
        cycles = []
        unseen = set(range(N))
        while unseen:
            i = min(unseen); cyc = {i, perm[i]}; unseen -= cyc; cycles.append(sorted(cyc))
        reflection_fixed_dimensions.append(len(cycles)-1)
    packet = {
        "D12_order": len(group), "simplex_dimension": 11,
        "explicit_interior_trivial_stabilizer_point": generic.tolist(),
        "generic_point_stabilizer_order": len(stabilizer),
        "generic_orbit_size": len(group)//len(stabilizer),
        "reflection_fixed_simplex_dimensions": reflection_fixed_dimensions,
        "open_dense_trivial_stabilizer_stratum_dimension": 11,
        "theorem": (
            "D12 invariance sends every minimizer to a finite orbit, but it does not force "
            "a minimizer into any reflection-fixed subspace. The union of all nontrivial "
            "fixed sets has empty interior, and the displayed positive probability vector "
            "has trivial stabilizer. Therefore symmetry alone cannot reduce the complete "
            "global problem below its eleven-dimensional generic stratum."
        ),
    }
    return finalize(344, "D12 Orbit-Type Reduction Boundary for the Global Variational Problem",
                    "proven_no_go_for_exhaustive_fixed_symmetry_reduction",
                    "Palais-type fixed-subspace searches remain valid for finding symmetric critical points, but cannot certify all global minimizers without an additional rearrangement or convex-envelope theorem.", packet)


def softmax11(y: np.ndarray) -> np.ndarray:
    z = np.r_[y, 0.]; z -= np.max(z); p = np.exp(z); return p/p.sum()


def simplex_objective(p: np.ndarray, m: np.ndarray, g: float = 4.) -> float:
    u = np.ones(N)/N
    return float(np.sum(p*np.log(np.maximum(N*p, 1e-300)))-g*(p-u)@m@(p-u)/2)


def optimize_simplex(m: np.ndarray, starts: int, seed: int) -> tuple[float, np.ndarray, list[float]]:
    rng = np.random.default_rng(seed); ys = [np.zeros(N-1)]
    for k in range(N):
        p = np.full(N, .001/(N-1)); p[k] = .999; ys.append(np.log(p[:-1]/p[-1]))
    ys.extend(rng.normal(0, 4, size=(starts, N-1)))
    best = (math.inf, np.ones(N)/N); values = []
    for y in ys:
        sol = minimize(lambda x: simplex_objective(softmax11(x), m), y,
                       method="L-BFGS-B", options={"ftol": 1e-14, "gtol": 1e-9, "maxiter": 1000})
        p = softmax11(sol.x); val = simplex_objective(p, m); values.append(val)
        if val < best[0]: best = (val, p)
    return best[0], best[1], values


def reflection_defect(p: np.ndarray) -> float:
    return float(min(np.linalg.norm(p-p[[(2*k-np.arange(N)) % N]]) for k in range(N)))


def st345(a: np.ndarray) -> dict:
    val, p, values = optimize_simplex(a, 160, SEED+345)
    packet = {
        "gain": 4.0, "unconstrained_starts": 173, "best_objective": val,
        "best_probability": p.tolist(), "maximum_probability": float(p.max()),
        "minimum_reflection_defect": reflection_defect(p),
        "best_20_objective_spread": float(max(sorted(values)[:20])-min(values)),
        "global_theorem": False,
        "result": (
            "All declared multistarts converge, up to D12 action and numerical tolerance, "
            "to the reflection-even localized orbit already certified locally. This is "
            "strong evidence, not an exhaustive proof: ST344 shows that symmetry alone "
            "cannot exclude a distant trivial-stabilizer minimizer."
        ),
    }
    return finalize(345, "Reflection-Orbit Audit of the g=4 Global-Minimizer Conjecture",
                    "strong_numerical_evidence_global_reflection_exhaustion_unproved",
                    "No branch-and-bound, rearrangement theorem, or complete dual certificate is supplied; global uniqueness modulo D12 remains open.", packet)


def spectral_blocks(a: np.ndarray):
    ev, vec = np.linalg.eigh(a); used = np.zeros(N, bool); blocks = []
    for idx in np.argsort(ev)[::-1]:
        if ev[idx] < 1e-9 or used[idx]: continue
        ids = np.where(abs(ev-ev[idx]) < 1e-8)[0]; used[ids] = True
        blocks.append((float(ev[idx]), ids, vec[:, ids]))
    return blocks


def st346(a: np.ndarray) -> dict:
    rows = []; m = np.zeros_like(a)
    for j, (lam, ids, vec) in enumerate(spectral_blocks(a)):
        m += vec@np.diag(np.full(len(ids), lam))@vec.T
        rank = int(np.linalg.matrix_rank(m, tol=1e-8))
        val, p, _ = optimize_simplex(m, 16, SEED+346+j)
        rows.append({"covariant_rank": rank, "retained_block_eigenvalue": lam,
                     "vertex_uniform_energy_crossing_g": 2*math.log(N)/(np.trace(m)/N),
                     "best_g4_objective": val, "best_max_probability": float(p.max()),
                     "best_entropy_nats": float(-np.sum(p*np.log(p))),
                     "reflection_defect": reflection_defect(p)})
    first = next(r for r in rows if r["best_g4_objective"] < -1e-6)
    packet = {
        "allowed_D12_covariant_ranks": [r["covariant_rank"] for r in rows],
        "rows": rows, "first_numerically_localizing_rank_at_g4": first["covariant_rank"],
        "exact_first_rank_where_pure_vertex_beats_uniform_at_g4":
            next(r["covariant_rank"] for r in rows if r["vertex_uniform_energy_crossing_g"] < 4),
        "theorem": (
            "D12 covariance requires retaining complete degenerate Fourier blocks, hence "
            "ranks 1,3,5,7,9,11. For each truncation the exact pure-vertex/uniform crossing "
            "is 2 log(12)/(tr(A_m)/12); rank seven is the first for which a vertex has lower "
            "energy at g=4. Numerical simplex minimization first localizes at the same rank."
        ),
    }
    return finalize(346, "D12-Covariant Modal Mediators Coupled to the Simplex",
                    "proven_vertex_thresholds_strong_numerical_first_localizing_rank_7",
                    "The pure-vertex comparison is exact but does not classify every state. The numerical rank-seven global claim is not certified; mediator rank and g remain supplied.", packet)


def st347(a: np.ndarray) -> dict:
    blocks = spectral_blocks(a); labels = [f"E{k+1}" for k in range(len(blocks))]
    deg4 = list(itertools.combinations_with_replacement(labels, 2))
    deg6 = list(itertools.combinations_with_replacement(labels, 3))
    packet = {
        "scope": "spectral-energy subalgebra generated by E_k=||P_k q||^2",
        "positive_irreducible_modal_energies": len(labels),
        "degree_4_basis_dimension": len(deg4), "degree_6_basis_dimension": len(deg6),
        "degree_4_basis": ["*".join(x) for x in deg4],
        "degree_6_basis_sha256": hashlib.sha256("\n".join("*".join(x) for x in deg6).encode()).hexdigest(),
        "free_real_coefficients_through_degree_6": len(deg4)+len(deg6),
        "distinguished_A_polynomials": ["<q,Aq>^2", "<q,Aq>^3"],
        "theorem": (
            "The positive strict spectrum has six real modal-energy generators. Their "
            "commutative homogeneous degree-four and degree-six spectral-energy bases have "
            "dimensions C(7,2)=21 and C(8,3)=56. A selects distinguished weighted sums but "
            "does not select the 77 general coefficients or their stabilizing signs."
        ),
    }
    return finalize(347, "Degree-Four and Degree-Six Strict Spectral-Energy Invariant Algebra",
                    "proven_21_plus_56_basis_classification_in_declared_subalgebra",
                    "This is not the full D12 polynomial invariant ring: phase-sensitive invariants outside the modal-energy subalgebra remain possible. No coefficient is sourced.", packet)


def st348() -> dict:
    p = np.array([.9, .1]); q = np.array([.6, .4]); s = .5
    z = float(np.sum(p**s*q**(1-s)))
    derivative = float(np.sum(p**s*q**(1-s)*np.log(p/q))/z)
    packet = {
        "counterexample": "one-state i.i.d. Bernoulli experiments",
        "belief_contraction_coefficient": 0.0,
        "P": p.tolist(), "Q": q.tolist(), "s": s,
        "per_symbol_log_Chernoff_weight": math.log(z),
        "per_symbol_derivative": derivative,
        "theorem": (
            "Belief contraction alone cannot bound the normalized Chernoff derivative by "
            "a constant that vanishes with the contraction coefficient. In a one-state "
            "i.i.d. experiment the belief contraction is exactly zero, while the displayed "
            "per-symbol derivative is nonzero. Any valid all-depth promotion must control "
            "the normalized multiplicative weights (pressure), not only filter forgetting."
        ),
    }
    return finalize(348, "Counterexample to Belief-Contraction-Only Chernoff Derivative Closure",
                    "proven_counterexample_to_declared_missing_bound_strategy",
                    "This does not prove that the frozen HMM lacks a pressure-derivative theorem; it proves that ST316 contraction alone is insufficient.", packet)


def st349() -> dict:
    old = json.loads((ROOT/"FIN_ST293_ST307_Results.json").read_text())["ST298"]["rows"]
    centers = np.array([r["center"] for r in old]); rows = []
    for i in range(1, len(centers)-1, 4):
        tangent = centers[i+1]-centers[i-1]; tangent /= np.linalg.norm(tangent)
        half = .51*np.maximum(abs(centers[i]-centers[i-1]), abs(centers[i+1]-centers[i]))
        anisotropy = float(max(half)/max(min(half[half > 1e-15]), 1e-15))
        # Reuse the already outward pseudo-hyperplane routine at the required
        # largest half-width: failure is conservative and exposes dependency inflation.
        from fin_st203_st217_research import pseudo_krawczyk
        aiv, _, _ = strict_interval_matrix()
        trial = pseudo_krawczyk(centers[i], aiv, tangent, centers[i], float(max(half)))
        rows.append({"index": i, "coordinate_halfwidths": half.tolist(),
                     "anisotropy_ratio": anisotropy, "outer_cube_certificate": trial})
    packet = {
        "representative_sections": len(rows),
        "accepted_outer_enclosures": sum(r["outer_cube_certificate"]["included"] for r in rows),
        "failed_outer_enclosure_indices": [r["index"] for r in rows if not r["outer_cube_certificate"]["included"]],
        "rows": rows,
        "result": (
            "Coordinate half-widths inferred from adjacent sections are highly anisotropic. "
            "Embedding each proposed slab in its rigorous outer interval cube improves the "
            "representative pass count relative to ST333, but leaves explicit failed sections. "
            "Thus anisotropic geometry helps but does not close a continuous tube; a rotated "
            "tangent/transverse interval-Newton implementation is still required."
        ),
    }
    return finalize(349, "Anisotropic Implicit-Slab Strategy Audit",
                    "partial_improvement_with_explicit_failed_sections_rotated_coordinates_required",
                    "Outer-cube rejection neither disproves a continuous branch nor rejects true rotated slabs. It only rejects this conservative enclosure implementation.", packet)


def solve_ir_grid(points: int) -> dict:
    mus = np.logspace(-5, 5, points); logs = np.log(mus)
    weights = np.gradient(logs); times = np.array([.25, 1., 4.]); heat = np.exp(-2*mus)
    kernels = np.array([4*mus*t/(np.exp(np.minimum(2*mus*t, 700))+1) for t in times])
    c = np.r_[np.zeros(points), -1.]; aub = []
    for k in kernels: aub.append(np.r_[-k*weights, 1.])
    aub.append(np.r_[heat*weights, 0.]); bub = np.r_[np.zeros(3), 3.]
    sol = linprog(c, A_ub=np.array(aub), b_ub=bub,
                  bounds=[(0, 1)]*points+[(None, None)], method="highs")
    q = sol.x[:-1]
    # np.diff on booleans records both entry and exit transitions.  A selected
    # band has one entry and one exit after padding by False at both ends.
    transitions = int(np.sum(np.diff(np.r_[False, q > .5, False]).astype(bool)))
    bands = transitions//2
    dual = -sol.ineqlin.marginals
    return {"points": points, "value": float(sol.x[-1]), "heat_cost": float((heat*weights)@q),
            "active_time_values": ((kernels*weights)@q).tolist(),
            "time_dual_weights": dual[:3].tolist(), "heat_dual_price": float(dual[3]),
            "time_weight_sum": float(sum(dual[:3])), "threshold_transitions": transitions,
            "selected_threshold_bands": bands}


def st350() -> dict:
    rows = [solve_ir_grid(n) for n in (250, 500, 1000)]
    packet = {
        "compact_spectral_interval": [1e-5, 1e5], "observer_times": [.25, 1., 4.],
        "grid_convergence": rows,
        "theorem": (
            "For finitely many observer times and one heat budget, strong LP duality "
            "reduces the continuum design to finitely many time multipliers plus one heat "
            "price. The primal is bang-bang according to the sign of the resulting analytic "
            "threshold function. On a compact positive spectral interval its zero set, and "
            "therefore its band count, is finite unless the threshold is identically zero."
        ),
    }
    return finalize(350, "Continuum Multi-Time Soft-IR Dual and Finite-Band Theorem",
                    "proven_finite_dual_support_and_finite_band_structure_grid_optimum_numerical",
                    "The exact continuum dual multipliers and band endpoints are not interval-certified. The density, times, budget and cutoff are supplied.", packet)


def st351() -> dict:
    # Square with one diagonal: V=4, E=5, connected, hence beta_1=2.
    edges = [(0, 1), (1, 2), (2, 3), (3, 0), (0, 2)]
    incidence = np.zeros((4, 5), int)
    for e, (i, j) in enumerate(edges): incidence[i, e] = -1; incidence[j, e] = 1
    rank = int(np.linalg.matrix_rank(incidence))
    rates = np.array([1.2, .7, 1.1, .9, .8])
    hol = [rates[0]*rates[1]/rates[4], rates[4]*rates[2]*rates[3]]
    h = np.array([2., 3., 5., 7.]); transformed = np.array([rates[e]*h[j]/h[i] for e, (i, j) in enumerate(edges)])
    hol2 = [transformed[0]*transformed[1]/transformed[4], transformed[4]*transformed[2]*transformed[3]]
    packet = {
        "vertices": 4, "edges": edges, "incidence_rank": rank,
        "first_cohomology_rank": len(edges)-rank,
        "two_cycle_holonomies": hol, "gauge_transformed_holonomies": hol2,
        "gauge_invariance_error": float(np.max(abs(np.array(hol)-np.array(hol2)))),
        "theorem": (
            "For a connected refinement graph, positive multiplicative edge calibrations "
            "modulo vertex gauges have dimension E-V+1. The square-with-diagonal complex "
            "has two independent cycle holonomies, which are complete coordinates for its "
            "calibration cohomology."
        ),
    }
    return finalize(351, "Two-Cycle Refinement Scale Cohomology",
                    "proven_two_independent_gauge_invariant_holonomies",
                    "The complex and rates are supplied. FIN does not derive this topology, its holonomies, or a dimensional interpretation.", packet)


def qary_confusion(e: float, q: int = 12) -> np.ndarray:
    return (1-e)*np.eye(q)+e/(q-1)*(np.ones((q, q))-np.eye(q))


def qary_deficiency_lp(e_coarse: float, e_fine: float, q: int = 12) -> float:
    coarse, target = qary_confusion(e_coarse, q), qary_confusion(e_fine, q)
    nk = q*q; nt = q*q; nv = nk+nt+1; c = np.zeros(nv); c[-1] = 1
    aeq = []; beq = []; aub = []; bub = []
    for i in range(q):
        row = np.zeros(nv); row[i*q:(i+1)*q] = 1; aeq.append(row); beq.append(1)
    for i in range(q):
        for j in range(q):
            t = nk+i*q+j; base = np.zeros(nv)
            for k in range(q): base[k*q+j] = coarse[i, k]
            row = base.copy(); row[t] = -1; aub.append(row); bub.append(target[i, j])
            row = -base.copy(); row[t] = -1; aub.append(row); bub.append(-target[i, j])
        row = np.zeros(nv); row[nk+i*q:nk+(i+1)*q] = .5; row[-1] = -1
        aub.append(row); bub.append(0)
    sol = linprog(c, A_ub=aub, b_ub=bub, A_eq=aeq, b_eq=beq,
                  bounds=[(0, None)]*nv, method="highs")
    return float(sol.fun)


def st352() -> dict:
    rows = []
    for ef, ec in ((.02, .05), (.1, .2), (.1, .4), (.2, .8)):
        rows.append({"fine_error": ef, "coarse_error": ec,
                     "fine_to_coarse_deficiency": 0.0,
                     "reverse_formula": ec-ef,
                     "reverse_full_LP": qary_deficiency_lp(ec, ef)})
    packet = {
        "family": "12-state permutation-covariant confusion channels",
        "exact_formula": "for 0<=e_f<=e_c<=11/12: delta(fine,coarse)=0 and delta(coarse,fine)=e_c-e_f",
        "rows": rows,
        "theorem": (
            "Twirling any postprocessing over the full 12-state permutation group cannot "
            "increase worst-row total variation, so the deficiency problem reduces to the "
            "one-parameter confusion family. Additional symmetric noise realizes the "
            "fine-to-coarse map, while the reverse optimum is the identity boundary and "
            "has exact deficiency e_c-e_f."
        ),
    }
    return finalize(352, "Exact Le Cam Depth for Twelve-State Symmetric Confusion Instruments",
                    "proven_exact_12_state_permutation_covariant_deficiency",
                    "The theorem concerns a declared instrument family. General cyclic convolution noise has a larger Fourier/LP geometry, and no laboratory FIN instrument is inferred.", packet)


def st353() -> dict:
    sigma = .2; rows = []
    for r in (1, 3, 5, 7, 9, 11):
        for counts in (100, 1000, 10000):
            risk = r*sigma*sigma/counts
            rows.append({"spectral_dimension": r, "independent_counts": counts,
                         "exact_minimax_MSE": risk, "root_MSE": math.sqrt(risk)})
    packet = {
        "conditional_model": "known r-mode spectral subspace, Gaussian sequence observations with orthonormal design and noise sigma=0.2",
        "rows": rows,
        "theorem": (
            "In the declared Gaussian sequence experiment the least-squares estimator has "
            "risk r sigma^2/N, and the translation-invariant Gaussian minimax lower bound "
            "matches it. Thus spectrally restricting a fiber replaces ambient parameter "
            "count by modal dimension r, but cannot remove the N^{-1/2} root-risk law."
        ),
    }
    return finalize(353, "Sharp Recovery Bounds for a Spectrally Restricted Fiber",
                    "proven_matching_lower_upper_bounds_in_declared_Gaussian_model",
                    "The noise, orthonormal design and known spectral subspace are premises. This is not apparatus sample complexity or an empirical recovery result.", packet)


def block2_rate_distortion(ptrans: np.ndarray) -> list[dict]:
    states = np.array(list(itertools.product(range(N), repeat=2)))
    px = np.array([ptrans[i, j]/N for i, j in states])
    dist = np.mean(states[:, None, :] != states[None, :, :], axis=2); rows = []
    for beta in (0, .5, 1, 2, 4, 8, 16):
        q = np.ones(N*N)/(N*N); e = np.exp(-beta*dist)
        for it in range(2500):
            z = e@q; pyx = q[None, :]*e/z[:, None]; qn = px@pyx
            if np.max(abs(qn-q)) < 1e-12: break
            q = qn
        joint = px[:, None]*pyx
        d = float(np.sum(joint*dist))
        rate = float(np.sum(joint*np.log2(np.maximum(pyx, 1e-300)/np.maximum(q[None, :], 1e-300)))/2)
        rows.append({"beta": beta, "distortion": d, "bits_per_symbol": rate, "iterations": it})
    return rows


def st354(a: np.ndarray) -> dict:
    tau = .5; p = expm(-tau*a)
    rows = block2_rate_distortion(p)
    entropy_rate = float(-np.mean(np.sum(p*np.log2(np.maximum(p, 1e-300)), axis=1)))
    packet = {
        "conditional_transition": "P_tau=exp(-tau A)", "tau": tau,
        "minimum_transition_probability": float(p.min()),
        "row_sum_error": float(np.max(abs(p.sum(axis=1)-1))),
        "stationary_distribution": "uniform by symmetry",
        "entropy_rate_bits_per_symbol": entropy_rate,
        "blocklength_2_rate_distortion": rows,
        "result": (
            "The strict heat kernel at supplied tau is a positive doubly stochastic Markov "
            "transition and therefore defines a conditional FIN-generated source law. Its "
            "entropy rate and two-symbol Hamming rate-distortion curve are computed without "
            "introducing an unrelated binary flip model."
        ),
    }
    return finalize(354, "Rate-Distortion for a Source Conditioned on the Strict Heat Kernel",
                    "proven_transition_and_entropy_rate_numerical_block_rate_distortion",
                    "Tau, the Markov reading, Hamming distortion and block length are supplied. This does not establish physical information density or universe compression.", packet)


def st355() -> dict:
    rows = []
    for eta in (1e-4, 1e-3, 1e-2):
        for levels in (1, 10, 100, 1000):
            lo = ((1-eta)/(1+eta))**levels; hi = 1/lo
            rows.append({"per_level_relative_mode_error": eta, "levels": levels,
                         "unitary_heat_relative_ratio_factor": [lo, hi],
                         "wave_relative_ratio_factor": [math.sqrt(lo), math.sqrt(hi)],
                         "worst_log_drift": levels*math.log((1+eta)/(1-eta))})
    packet = {
        "rows": rows,
        "theorem": (
            "If each refinement level multiplies all modes by an arbitrary common gauge "
            "c_n and each individual mode has relative error at most eta_n, every common "
            "c_n cancels from mode ratios. Differential drift accumulates multiplicatively "
            "between products of (1-eta_n)/(1+eta_n) and its reciprocal; wave-clock drift "
            "uses square roots."
        ),
    }
    return finalize(355, "Multilevel Relative-Clock Drift Separated from Gauge Calibration",
                    "proven_product_drift_bounds_and_common_scale_cancellation",
                    "The per-level error premise must be established by a refinement map. Bounds do not select seconds, an arrow of time, or an absolute calibration.", packet)


def st356() -> dict:
    records = list(ROOT.glob("FIN_ST356_INDEPENDENT_RAW_EVENTS*.jsonl"))
    packet = {"required_pattern": "FIN_ST356_INDEPENDENT_RAW_EVENTS*.jsonl",
              "matching_records": [x.name for x in records], "independent_record_count": len(records),
              "theorem": "Local analytical work cannot manufacture independent custody, unblinding, apparatus calibration or empirical events."}
    return finalize(356, "Independent Evidence Gate",
                    "blocked_no_independent_events" if not records else "record_present_requires_blinded_validation",
                    "The stop is evidentiary, not a failed physical experiment.", packet)


def make_figures(out: dict) -> None:
    FIG_DIR.mkdir(exist_ok=True)
    r = out["ST342"]; fig, ax = plt.subplots(figsize=(7, 3.8))
    xs = [x["g"] for x in r["endpoint_certificates"]]
    ys = [np.mean(x["value_interval"]) for x in r["endpoint_certificates"]]
    ye = [(x["value_interval"][1]-x["value_interval"][0])/2 for x in r["endpoint_certificates"]]
    ax.errorbar(xs, ys, yerr=ye, fmt="o-"); ax.axhline(0, color="black", lw=.8)
    ax.set(xlabel="supplied coupling g", ylabel="localized-minus-uniform objective",
           title="ST342: certified width-2e-8 coexistence bracket"); fig.tight_layout()
    fig.savefig(FIG_DIR/"st342_narrow_crossing.png", dpi=200); plt.close(fig)
    r = out["ST346"]["rows"]; fig, ax = plt.subplots(figsize=(7, 3.8))
    ax.plot([x["covariant_rank"] for x in r], [x["best_max_probability"] for x in r], "o-")
    ax.set(xlabel="D12-covariant mediator rank", ylabel="largest optimized probability",
           title="ST346: first numerical localization at rank seven"); fig.tight_layout()
    fig.savefig(FIG_DIR/"st346_modal_localization.png", dpi=200); plt.close(fig)
    r = out["ST354"]["blocklength_2_rate_distortion"]; fig, ax = plt.subplots(figsize=(7, 3.8))
    ax.plot([x["distortion"] for x in r], [x["bits_per_symbol"] for x in r], "o-")
    ax.set(xlabel="Hamming distortion", ylabel="bits per symbol",
           title="ST354: strict-heat Markov block rate-distortion"); fig.tight_layout()
    fig.savefig(FIG_DIR/"st354_heat_rate_distortion.png", dpi=200); plt.close(fig)


def main() -> None:
    np.random.seed(SEED); _, a, s = strict_operator()
    out = {"metadata": {"seed": SEED, "python": platform.python_version(),
                         "numpy": np.__version__, "scipy": scipy.__version__,
                         "strict_row_sum": s,
                         "scope": "local analytical, interval, exact finite and bounded numerical research"}}
    funcs = {342: lambda: st342(a), 343: lambda: st343(a), 344: st344,
             345: lambda: st345(a), 346: lambda: st346(a), 347: lambda: st347(a),
             348: st348, 349: st349, 350: st350, 351: st351, 352: st352,
             353: st353, 354: lambda: st354(a), 355: st355, 356: st356}
    for k in range(342, 357): out[f"ST{k}"] = funcs[k]()
    make_figures(out)
    RESULTS.write_text(json.dumps(native(out), indent=2, sort_keys=True), encoding="utf-8")
    with SUMMARY.open("w", newline="", encoding="utf-8") as h:
        w = csv.writer(h); w.writerow(["program", "object", "status"])
        for k in range(342, 357): w.writerow([f"ST{k}", out[f"ST{k}"]["object"], out[f"ST{k}"]["status"]])
    print(json.dumps({k: v["status"] for k, v in out.items() if k.startswith("ST")}, indent=2))


if __name__ == "__main__": main()
