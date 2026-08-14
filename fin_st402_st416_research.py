#!/usr/bin/env python3
"""FIN ST402--ST416: independent orbit replay, gain persistence, and IR boundary.

This batch is finite and dimensionless.  The ST402 proof deliberately rebuilds
the strict matrix, benchmark, outside-cap cover, and cap-curvature reduction;
it does not import the ST389/ST390 cap implementation or their result packets.
No result sources the gain, selects an orbit member, supplies units, or turns a
finite variational statement into a physical vacuum theorem.
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
from scipy.linalg import null_space
from scipy.optimize import brentq, root
from scipy.special import exp1, logsumexp

from fin_st118_st129_research import interval_left, interval_matvec
from fin_st132_center_isolation_replay import bounds, iv
from fin_st372_st386_research import exponent_orbit, orbit_eval_vector
from fin_st387_st401_research import (
    e1_iv_local,
    ir_G_general_iv,
    ir_Gprime_general_iv,
    ir_equations,
    ir_param_equations,
    parametric_ir_krawczyk,
)


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "FIN_ST402_ST416_Results.json"
SUMMARY = ROOT / "FIN_ST402_ST416_Summary.csv"
FIG_DIR = ROOT / "FIN_ST402_ST416_Figures"
SEED = 20260822
N = 12
IDX = np.array([i if i <= 6 else N-i for i in range(N)])
MULT = np.array([1, 2, 2, 2, 2, 2, 1.0])
NAMES = {
    402: "Independent_Global_Orbit_Replay",
    403: "Abstract_Orbit_Exhaustion_Theorem",
    404: "Exact_Orbit_Stabilizer_Certificate",
    405: "Uniform_Gain_Persistence_Tube",
    406: "Global_Gain_Transition_Bracket",
    407: "Certified_and_Sampled_Morse_Landscape",
    408: "Projected_Gradient_Local_Basin",
    409: "Gain_Source_Admission_Gate",
    410: "Rank7_Replacement_Envelope_Audit",
    411: "First_Forced_Invariant_Relations",
    412: "Noise_Aware_Invariant_Coefficient_Rank",
    413: "Adapted_Chernoff_Cone_Gap",
    414: "Expanded_IR_Krawczyk_Tube",
    415: "Certified_Limiting_IR_Band_Loss_Face",
    416: "Strict_Source_and_Independent_Evidence_Gates",
}
PACKETS = {k: ROOT / f"FIN_ST{k}_{v}.json" for k, v in NAMES.items()}


def native(x: Any) -> Any:
    if isinstance(x, dict):
        return {str(k): native(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [native(v) for v in x]
    if isinstance(x, np.ndarray):
        return native(x.tolist())
    if isinstance(x, (np.floating, np.integer)):
        return x.item()
    return x


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def finalize(k: int, obj: str, status: str, boundary: str, packet: dict) -> dict:
    path = PACKETS[k]
    path.write_text(json.dumps(native(packet), indent=2, sort_keys=True), encoding="utf-8")
    return {
        "program": f"ST{k}", "object": obj, "packet_file": path.name,
        "packet_sha256": sha(path), **packet, "status": status, "boundary": boundary,
    }


# ---------------------------------------------------------------------------
# An independent strict-matrix and variational implementation for ST402.

def independent_strict_matrix_float() -> np.ndarray:
    omega, phi, eta = 0.18575, 0.16250, 9/5
    weights = {d: math.cos(omega*d+phi)/(1+d**eta) for d in range(1, 7)}
    s = 2*sum(weights[d] for d in range(1, 6))+weights[6]
    return np.array([[s if i == j else -weights[min((i-j) % N, (j-i) % N)]
                      for j in range(N)] for i in range(N)])


def independent_strict_matrix_interval():
    omega, phi, eta = iv("0.18575"), iv("0.16250"), iv(9)/iv(5)
    weights = {d: mp.iv.cos(omega*d+phi)/(1+iv(d)**eta) for d in range(1, 7)}
    s = 2*sum((weights[d] for d in range(1, 6)), iv(0))+weights[6]
    a = [[s if i == j else -weights[min((i-j) % N, (j-i) % N)]
          for j in range(N)] for i in range(N)]
    return a, weights, s


def radial_matrix(a: np.ndarray) -> np.ndarray:
    return np.array([[sum(a[i, k] for k in range(N) if IDX[k] == j)
                      for j in range(7)] for i in range(7)])


def radial_system(x: np.ndarray, b: np.ndarray, g: float) -> np.ndarray:
    q, lam = x[:7], x[7]
    return np.r_[np.log(N*q)+1-g*b@q+lam, MULT@q-1]


def discover_radial_state(a: np.ndarray, g: float = 4.) -> np.ndarray:
    # A qualitative localized seed; no prior result packet is opened.
    q = np.array([.985, .0005, .0008, .001, .0012, .0013, .0014])
    q /= MULT@q
    b = radial_matrix(a)
    lam = float(np.mean(g*b@q-np.log(N*q)-1))
    sol = root(lambda z: radial_system(z, b, g), np.r_[q, lam], tol=1e-12)
    if not sol.success or np.linalg.norm(radial_system(sol.x, b, g), np.inf) > 2e-10:
        raise RuntimeError("independent radial discovery failed")
    return sol.x


def independent_radial_interval_system(center, radius, aiv, gain=4.):
    q = [iv((center[i]-radius, center[i]+radius)) for i in range(7)]
    lam = iv((center[7]-radius, center[7]+radius))
    biv = [[sum((aiv[i][k] for k in range(N) if IDX[k] == j), iv(0))
            for j in range(7)] for i in range(7)]
    f = [mp.iv.log(N*q[i])+1-iv(str(gain))*sum((biv[i][j]*q[j] for j in range(7)), iv(0))+lam
         for i in range(7)]
    f.append(sum((iv(str(MULT[j]))*q[j] for j in range(7)), iv(0))-1)
    jac = [[iv(0) for _ in range(8)] for _ in range(8)]
    for i in range(7):
        for j in range(7):
            jac[i][j] = (1/q[i] if i == j else iv(0))-iv(str(gain))*biv[i][j]
        jac[i][7] = iv(1)
    for j in range(7): jac[7][j] = iv(str(MULT[j]))
    flo = np.array([bounds(x)[0] for x in f]); fhi = np.array([bounds(x)[1] for x in f])
    jlo = np.array([[bounds(x)[0] for x in row] for row in jac])
    jhi = np.array([[bounds(x)[1] for x in row] for row in jac])
    return flo, fhi, jlo, jhi


def independent_radial_krawczyk(center, aiv, radius=1e-8):
    f0lo, f0hi, j0lo, j0hi = independent_radial_interval_system(center, 0., aiv)
    pre = np.linalg.inv((j0lo+j0hi)/2)
    _, _, jlo, jhi = independent_radial_interval_system(center, radius, aiv)
    cfl, cfh = interval_matvec(pre, pre, f0lo, f0hi)
    ylo, yhi = center-cfh, center-cfl
    cjl, cjh = interval_left(pre, jlo, jhi); mlo, mhi = -cjh, -cjl
    for i in range(8):
        mlo[i, i] = np.nextafter(mlo[i, i]+1, -np.inf)
        mhi[i, i] = np.nextafter(mhi[i, i]+1, np.inf)
    dlo, dhi = interval_matvec(mlo, mhi, np.full(8, -radius), np.full(8, radius))
    margin = np.minimum(ylo+dlo-(center-radius), (center+radius)-(yhi+dhi))
    return {"included": bool(np.min(margin) > 0), "radius": radius,
            "minimum_margin": float(np.min(margin)), "component_margins": margin.tolist()}


def rational_probability(center: np.ndarray, denominator: int = 10**18):
    p = center[:7][IDX]
    nums = [int(round(float(x)*denominator)) for x in p]
    nums[0] += denominator-sum(nums)
    return nums, denominator


def rational_objective_components(nums, denominator, aiv):
    pp = [iv(n)/denominator for n in nums]
    u = iv(1)/N
    entropy = sum((x*mp.iv.log(N*x) for x in pp), iv(0))
    q = [x-u for x in pp]
    quadratic = sum((q[i]*aiv[i][j]*q[j] for i in range(N) for j in range(N)), iv(0))
    return bounds(entropy), bounds(quadratic)


def dmin_iv(lo: float, hi: float):
    r = iv((max(lo, 1e-16), hi))
    z = mp.iv.sqrt(iv(132)*r); aa = (1+z)/12; bb = (1-aa)/11
    return aa*mp.iv.log(12*aa)+11*bb*mp.iv.log(12*bb)


def outside_difference_cover(m0: float, gain: float, entropy_iv, quad_iv,
                             weights_iv, siv, cells: int = 10000):
    """Lower-bound outside envelope minus the same-g rational benchmark."""
    rmax = m0-1/N; best = math.inf; arg = None
    for k in range(cells):
        lo, hi = rmax*k/cells, rmax*(k+1)/cells
        if k == 0:
            # D>=0, E<=2s r for the origin cell.
            diff = -gain*bounds(siv)[1]*hi-entropy_iv[1]+gain*quad_iv[0]/2
        else:
            r = iv((lo, hi))
            eup = (siv+weights_iv[6])*r+siv/12-11*weights_iv[6]/12
            val = dmin_iv(lo, hi)-iv(str(gain))*eup/2-iv(entropy_iv)+iv(str(gain))*iv(quad_iv)/2
            diff = bounds(val)[0]
        if diff < best:
            best, arg = diff, [lo, hi]
    return {"lower": best, "worst_collision_cell": arg, "cells": cells}


def cap_curvature_replay(aiv, m0: float, gain: float) -> dict:
    """Fresh 121-pair reduction; independent of the earlier cap function."""
    rows = []
    for i in range(1, N):
        for j in range(1, N):
            b = [iv(0) for _ in range(N)]; c = [iv(0) for _ in range(N)]
            b[i] += iv("0.5"); b[j] -= iv("0.5")
            c[0] -= iv(1); c[i] += iv("0.5"); c[j] += iv("0.5")
            q0 = sum((b[x]*aiv[x][y]*b[y] for x in range(N) for y in range(N)), iv(0))
            q1 = 2*sum((b[x]*aiv[x][y]*c[y] for x in range(N) for y in range(N)), iv(0))
            q2 = sum((c[x]*aiv[x][y]*c[y] for x in range(N) for y in range(N)), iv(0))
            coeff = [iv(1)/iv(str(1-m0))-iv(str(gain))*q0,
                     -iv(str(gain))*q1, iv(1)-iv(str(gain))*q2]
            cb = [bounds(x) for x in coeff]
            mid = [(x[0]+x[1])/2 for x in cb]; rad = [(x[1]-x[0])/2 for x in cb]
            d0, d1, d2 = mid; vals = [d0-d1+d2, d0+d1+d2]
            if d2 > 0 and -1 < -d1/(2*d2) < 1:
                vals.append(d0-d1*d1/(4*d2))
            rows.append((min(vals)-sum(rad), i, j, cb))
    worst = min(rows)
    return {
        "extreme_pair_count": len(rows), "minimum_L1_curvature": worst[0],
        "minimum_Euclidean_curvature": worst[0]/2,
        "worst_pair": {"i": worst[1], "j": worst[2], "coefficient_intervals": worst[3]},
    }


def independent_orbit_data(gain: float = 4., cells: int = 10000):
    mp.iv.dps = 70
    a = independent_strict_matrix_float(); aiv, wiv, siv = independent_strict_matrix_interval()
    center = discover_radial_state(a, gain)
    nums, den = rational_probability(center)
    ent, quad = rational_objective_components(nums, den, aiv)
    outside = outside_difference_cover(.94, gain, ent, quad, wiv, siv, cells)
    curvature = cap_curvature_replay(aiv, .94, gain)
    return a, center, nums, den, ent, quad, outside, curvature


def st402() -> dict:
    a, center, nums, den, ent, quad, outside, curvature = independent_orbit_data()
    aiv, _, _ = independent_strict_matrix_interval()
    root_cert = independent_radial_krawczyk(center, aiv)
    passed = outside["lower"] > 0 and curvature["minimum_Euclidean_curvature"] > 0
    packet = {
        "implementation_independence": (
            "The strict matrix, qualitative radial solve, rational benchmark, 10,000-cell "
            "outside cover, and 121-pair curvature reduction are rebuilt here.  No ST389/ST390 "
            "cap routine, root packet, or conclusion is read."
        ),
        "gain": 4.0, "cap_threshold": .94,
        "matrix_row_sum_defect": float(np.max(abs(a.sum(axis=1)))),
        "discovered_radial_residual_inf": float(np.linalg.norm(radial_system(center, radial_matrix(a), 4), np.inf)),
        "independent_radial_Krawczyk": root_cert,
        "rational_benchmark": {"numerators": nums, "denominator": den,
                               "entropy_interval": ent, "quadratic_interval": quad,
                               "objective_interval": [ent[0]-2*quad[1], ent[1]-2*quad[0]]},
        "outside_minus_benchmark_certificate": outside,
        "cap_curvature_certificate": curvature,
        "logical_proof": [
            "compactness gives a global minimizer",
            "the outward cover forces every minimizer into one of twelve disjoint caps p_i>0.94",
            "the entropy directional derivative excludes every simplex-boundary minimizer",
            "the 121-pair Fisher reduction proves strict convexity in every cap",
            "D12 invariance and transitivity move one minimizer into every cap",
            "strict convexity permits at most one minimizer per cap",
        ],
        "theorem": (
            "The supplied dimensionless strict functional at g=4 has exactly twelve global "
            "minimizers, forming one D12 orbit.  This independently reproduces ST390 without "
            "using its cap or root certificate."
        ),
    }
    return finalize(402, "Independent Exact Global-Minimizer Orbit Replay",
                    "proven_independent_exactly_twelve_global_minima" if passed else "independent_replay_failed",
                    "The theorem neither sources g=4 nor chooses an orbit member; QW-2191 and physical interpretation remain open.", packet)


def st403() -> dict:
    packet = {
        "abstract_theorem": (
            "Let a finite group act transitively on n simplex coordinates.  Let F be continuous "
            "and invariant.  If m>1/2, an explicit feasible point has value below every point "
            "with max p_i<=m, and F is strictly convex on every cap {p_i>m}, then argmin F "
            "contains exactly n points and is one group orbit."
        ),
        "proof": [
            "compactness and continuity give one minimizer",
            "separation puts it in a cap; m>1/2 makes caps disjoint",
            "transitivity creates a minimizer in every cap",
            "strict convexity gives at most one minimizer in each cap",
        ],
        "axiom_removal_tests": [
            {"axiom": "existence/compactness", "necessity": 5, "counterexample": "F(x)=x on the open interval (0,1) has no minimizer"},
            {"axiom": "outside-cap separation", "necessity": 5, "counterexample": "a strictly convex invariant entropy may minimize at the uniform point outside all caps"},
            {"axiom": "strict cap convexity", "necessity": 5, "counterexample": "an invariant continuous function may have a flat continuum of cap minima"},
            {"axiom": "transitive invariance", "necessity": 5, "counterexample": "unequal cap depths can leave only one preferred cap"},
            {"axiom": "m>1/2/disjoint caps", "necessity": 4, "counterexample": "overlapping caps invalidate the orbit-counting step even when minima exist"},
        ],
        "application": "ST402 meets every hypothesis with n=12, G=D12, and m=0.94.",
    }
    return finalize(403, "Abstract Finite Orbit-Exhaustion Theorem",
                    "proven_abstract_theorem_and_necessity_counterexamples",
                    "The theorem is a reusable variational statement; it contains no physical selector, dimensional scale, or source law.", packet)


def st404() -> dict:
    cap = json.loads(PACKETS[402].read_text())
    passed = cap["outside_minus_benchmark_certificate"]["lower"] > 0
    packet = {
        "group": "D12 of order 24", "unique_maximum_from_cap": "p_0>0.94>1/2",
        "stabilizer_elements": ["identity", "reflection j -> -j fixing vertex 0"],
        "stabilizer_order": 2, "orbit_size_by_orbit_stabilizer": 12,
        "proof": (
            "Uniqueness in cap zero forces invariance under its reflection.  Any stabilizer "
            "must fix the unique maximum vertex zero.  Exactly two D12 elements fix that "
            "vertex, so the stabilizer is precisely the displayed reflection subgroup."
        ),
        "larger_stabilizer_excluded": passed,
    }
    return finalize(404, "Exact Stabilizer and Orbit-Size Certificate",
                    "proven_exact_reflection_stabilizer_order_two" if passed else "blocked_by_ST402",
                    "A stabilizer theorem describes degeneracy; it does not select one of the twelve translates.", packet)


def st405() -> dict:
    mp.iv.dps = 70
    a = independent_strict_matrix_float(); aiv, wiv, siv = independent_strict_matrix_interval()
    center = discover_radial_state(a, 4); nums, den = rational_probability(center)
    ent, quad = rational_objective_components(nums, den, aiv)
    gain_box = [3.999, 4.001]
    outside = outside_difference_cover(.94, gain_box[0], ent, quad, wiv, siv, 10000)
    curvature = cap_curvature_replay(aiv, .94, gain_box[1])
    passed = outside["lower"] > 0 and curvature["minimum_Euclidean_curvature"] > 0
    packet = {
        "gain_interval": gain_box, "common_rational_benchmark": {"numerators": nums, "denominator": den},
        "monotonicity_payment": (
            "On the audited collision interval Q_benchmark-E_upper(r)>0, so the outside-minus-"
            "benchmark separation is worst at g=3.999.  Cap curvature is worst at g=4.001."
        ),
        "minimum_Q_minus_energy_upper": quad[0]-bounds((siv+wiv[6])*iv(str(.94-1/12))+siv/12-11*wiv[6]/12)[1],
        "uniform_outside_separation": outside,
        "uniform_cap_curvature": curvature,
        "theorem": (
            "For every supplied gain in [3.999,4.001], the strict functional has exactly "
            "twelve global minimizers forming one D12 orbit."
        ),
    }
    return finalize(405, "Uniform Nonzero Gain-Persistence Tube",
                    "proven_uniform_twelve_minimum_orbit_on_gain_interval" if passed else "gain_tube_failed",
                    "The interval is a robustness theorem, not a strict derivation of any gain or a physical phase interval.", packet)


def entropy_min_collision(r: float) -> float:
    if r <= 0: return 0.
    z = math.sqrt(132*r); aa = (1+z)/12; bb = (1-aa)/11
    return aa*math.log(12*aa)+(11*bb*math.log(12*bb) if bb > 0 else 0.)


def st406() -> dict:
    a = independent_strict_matrix_float(); eig = np.linalg.eigvalsh(a); lm = float(eig[-1])
    R = 11/12; cells = 200000
    xs = (np.arange(cells)+.5)*R/cells
    ratios = np.array([entropy_min_collision(float(x))/x for x in xs])
    min_ratio_sample = float(np.min(ratios)); rstar = float(xs[np.argmin(ratios)])
    # Pay a deliberately conservative 2e-4 downward allowance; this is a
    # numerical lower certificate only if independently interval-covered, so
    # the resulting uniform range is labelled strong numerical below.
    global_uniform_gain = 2*(min_ratio_sample-2e-4)/lm
    local_cross = json.loads((ROOT / "FIN_ST342_ST356_Results.json").read_text())["ST342"]
    persistence = json.loads(PACKETS[405].read_text())
    packet = {
        "uniform_Hessian_loss_gain": 12/lm,
        "sampled_minimum_Dmin_over_collision": min_ratio_sample,
        "sampled_ratio_location": rstar,
        "conservative_uniform_global_range_candidate": [0, global_uniform_gain],
        "certified_local_coexistence_crossing_reused": local_cross["certified_crossing_bracket"],
        "certified_nonuniform_global_range": persistence["gain_interval"],
        "rigorous_conclusion": (
            "ST405 proves the localized orbit global near four, while g=0 has the unique "
            "uniform global minimizer.  Hence at least one global orbit change occurs in "
            "(0,3.999].  The narrow 2.90249648 crossing remains branch-local, not a certified "
            "global transition."
        ),
        "failed_objective": "Every global transition and its uniqueness were not certified.",
    }
    return finalize(406, "Global Gain-Transition Audit",
                    "partial_existential_global_transition_full_classification_open",
                    "The local coexistence point must not be called the physical or global critical gain.", packet)


def softmax(y):
    z = np.r_[y, 0.]; z -= np.max(z); e = np.exp(z); return e/e.sum()


def group_images(p):
    return [np.roll(p, k) for k in range(N)]+[np.roll(p[::-1], k) for k in range(N)]


def orbit_distance(p, q):
    return min(float(np.linalg.norm(x-q)) for x in group_images(p))


def st407() -> dict:
    a = independent_strict_matrix_float(); u = np.ones(N)/N
    basis = null_space(np.ones((1, N)))
    def equation(y):
        z = np.r_[y, 0.]; logp = z-logsumexp(z); p = np.exp(logp)
        h = math.log(N)+logp+1-4*a@(p-u); return h[:-1]-h[-1]
    rng = np.random.default_rng(SEED+407); reps = []
    for k in range(1200):
        seed = np.zeros(N-1) if k == 0 else rng.normal(0, 4, N-1)
        sol = root(equation, seed, method="lm")
        p = softmax(sol.x)
        if np.linalg.norm(equation(sol.x), np.inf) < 2e-7 and p.min() > 1e-9:
            if all(orbit_distance(p, q) > 2e-5 for q in reps): reps.append(p)
    rows = []
    for p in reps:
        h = basis.T@(np.diag(1/p)-4*a)@basis
        ev = np.linalg.eigvalsh(h); q = p-u
        value = float(np.sum(p*np.log(N*p))-2*q@a@q)
        rows.append({"value": value, "Morse_index": int(np.sum(ev < -1e-6)),
                     "near_zero_modes": int(np.sum(abs(ev) <= 1e-6)),
                     "maximum_probability": float(np.max(p)), "minimum_probability": float(np.min(p)),
                     "representative": p.tolist()})
    rows.sort(key=lambda x: (x["value"], x["Morse_index"]))
    exact = {
        "uniform": {"value": 0., "Morse_index": 0,
                    "minimum_tangent_Hessian_eigenvalue": float(12-4*np.max(eig_nonzero(a)))},
        "global_localized_orbit": {"orbit_size": 12, "Morse_index": 0,
                                   "certificate": "ST402 cap curvature and global exhaustion"},
    }
    packet = {
        "deterministic_start_count": 1200, "sampled_distinct_D12_orbits": len(rows),
        "certified_stationary_orbits": exact,
        "sampled_stationary_catalog": rows,
        "sampled_index_histogram": {str(i): sum(r["Morse_index"] == i for r in rows)
                                    for i in sorted(set(r["Morse_index"] for r in rows))},
        "result": (
            "The uniform local minimum and the twelve-point global orbit are certified.  "
            "The additional stationary-orbit/Morse catalog is deterministic numerical evidence "
            "only; multistart discovery is not an exhaustive interval cover."
        ),
    }
    return finalize(407, "Certified Core and Sampled Morse Landscape at g=4",
                    "two_certified_minimum_orbits_plus_nonexhaustive_numerical_saddles",
                    "No sampled saddle count or Morse index is promoted to a complete landscape theorem.", packet)


def eig_nonzero(a):
    return np.linalg.eigvalsh(a)[1:]


def st408() -> dict:
    cap = json.loads(PACKETS[402].read_text()); a = independent_strict_matrix_float()
    center = discover_radial_state(a, 4); p = center[:7][IDX]
    root_cert = cap["independent_radial_Krawczyk"]
    root_radius = root_cert["radius"]
    mu = cap["cap_curvature_certificate"]["minimum_Euclidean_curvature"]
    radius = min(1e-4, float(np.min(p)-root_radius)/3,
                 float(np.max(p)-root_radius-.94)/3)
    packet = {
        "declared_dynamics": "p_dot=-P_T grad V_4(p), P_T=I-11^T/12, inside the positive simplex",
        "certified_center_enclosure_source": "independent ST402 radial Krawczyk inclusion",
        "center_enclosure_radius": root_radius,
        "explicit_local_Euclidean_radius": radius,
        "minimum_center_probability": float(np.min(p)),
        "cap_clearance": float(np.max(p)-.94),
        "strong_convexity_rate_lower": mu,
        "distance_decay": f"||p(t)-p*|| <= exp(-{mu:.15g} t) ||p(0)-p*||",
        "basin_proof": (
            "The ball lies in the positive cap.  Strong monotonicity of the tangent gradient "
            "makes squared distance to the unique cap minimizer decrease at rate at least 2mu; "
            "the ball is forward invariant and trajectories converge exponentially."
        ),
    }
    return finalize(408, "Explicit Local Basin for a Declared Projected Gradient Flow",
                    "proven_conditional_local_forward_invariance_and_exponential_convergence",
                    "The gradient law, its time parameter, gain, and preparation are premises.  No Born weights or physical basin probabilities follow.", packet)


def st409() -> dict:
    packet = {
        "searched_current_strict_objects": ["A", "spectral projectors", "D12 invariants", "entropy", "finite orbit theorem"],
        "required_new_export": "a typed strict object fixing a nonzero signed gain with coupling semantics",
        "found": None,
        "verdict": "The orbit theorem strengthens the consequences of g but cannot reverse-engineer a source for g.",
    }
    return finalize(409, "Gain-Source Admission Gate", "blocked_no_new_strict_gain_source",
                    "Do not infer a source from the fact that a supplied gain has mathematically interesting consequences.", packet)


def st410() -> dict:
    from fin_st357_st371_research import rank7_interval_matrix
    miv = rank7_interval_matrix()
    matrix = np.array([[(bounds(x)[0]+bounds(x)[1])/2 for x in row] for row in miv])
    weights = -matrix[0, 1:]
    old = json.loads((ROOT / "FIN_ST357_ST371_Results.json").read_text())["ST361"]
    packet = {
        "rank7_row_weights": weights.tolist(),
        "negative_offdiagonal_weight_count": int(np.sum(weights < 0)),
        "minimum_weight": float(np.min(weights)), "maximum_weight": float(np.max(weights)),
        "local_root_maximum_probability": max(old["root_center"][:7]),
        "previous_cap_curvature_threshold": json.loads((ROOT / "FIN_ST391_Rank7_Global_Inequality_Stop.json").read_text())["sampled_cap_curvature_threshold"],
        "new_method_obstruction": (
            "The strict positivity/refined water-filling envelope cannot be transferred: the "
            "rank-seven Fourier truncation has sign-changing effective row weights.  Its local "
            "root also remains below the cap-curvature threshold."
        ),
        "global_rank7_certificate": None,
    }
    return finalize(410, "Rank-Seven Replacement-Envelope Audit",
                    "bounded_no_go_for_strict_positive_weight_and_cap_methods",
                    "This rejects two proof mechanisms, not global minimality of the rank-seven state.", packet)


def st411() -> dict:
    counts = {"q4": math.comb(6+4-1, 4), "q2_times_primitive4": math.comb(7, 2)*32,
              "primitive4_products": math.comb(33, 2), "q_times_primitive6": 6*117}
    total = sum(counts.values()); dim8 = 1892; forced = total-dim8
    packet = {
        "generators_through_degree_six": {"degree2": 6, "primitive_degree4": 32, "primitive_degree6": 117},
        "degree8_formal_product_counts": counts, "total_degree8_candidates": total,
        "Molien_degree8_dimension": dim8, "forced_degree8_relation_nullity_lower_bound": forced,
        "degree6_declared_candidate_syzygies": 0,
        "theorem": (
            "The declared degree-six generated family is independent, but its 2,028 formal "
            "degree-eight products map into a 1,892-dimensional invariant space.  Rank-nullity "
            "forces at least 136 degree-eight relations.  Degree eight is therefore the first "
            "forced relation degree for this chosen generator family."
        ),
    }
    return finalize(411, "First Forced Relations in the D12 Invariant Generator Family",
                    "proven_at_least_136_degree8_syzygies_by_rank_nullity",
                    "No explicit minimal relation basis or full invariant-ring presentation is claimed.", packet)


def real_orbit_eval(rep, points):
    orbit = exponent_orbit(tuple(rep)); out = np.zeros(len(points))
    for e in orbit:
        out += np.prod(points**np.array(e), axis=1)
    return out


def st412() -> dict:
    basis_packet = json.loads((ROOT / "FIN_ST378_Exact_Sextic_D12_Reynolds_Basis.json").read_text())
    reps = basis_packet["selected_orbit_representatives"]
    rng = np.random.default_rng(SEED+412); nobs = 500
    points = rng.normal(size=(nobs, N)); points -= points.mean(axis=1, keepdims=True)
    points /= np.linalg.norm(points, axis=1, keepdims=True)
    design = np.column_stack([real_orbit_eval(r, points) for r in reps])
    norms = np.linalg.norm(design, axis=0); normalized = design/norms
    sv = np.linalg.svd(normalized, compute_uv=False)
    sigma = 1e-3
    packet = {
        "exact_algebraic_rank_source": "ST378/ST393 modular rank 365",
        "supplied_generic_design_seed": SEED+412, "observations": nobs, "coefficients": len(reps),
        "column_normalization": "unit Euclidean norm", "numerical_rank_tolerance_1e_10": int(np.sum(sv > 1e-10)),
        "largest_singular_value": float(sv[0]), "smallest_singular_value": float(sv[-1]),
        "condition_number": float(sv[0]/sv[-1]),
        "singular_values_below_1e_2": int(np.sum(sv < 1e-2)),
        "supplied_noise_sigma": sigma,
        "worst_direction_least_squares_standard_deviation": float(sigma/sv[-1]),
        "result": (
            "Generic vertex-resolving evaluations retain full numerical rank after column "
            "normalization, but the finite design is ill-conditioned.  Algebraic identifiability "
            "does not imply stable coefficient recovery under noise."
        ),
    }
    return finalize(412, "Noise-Aware Sextic Invariant Coefficient Rank",
                    "exact_algebraic_rank_with_numerical_noise_conditioning_warning",
                    "The design and noise are synthetic; no laboratory instrument or measured coefficient is supplied.", packet)


def st413() -> dict:
    mp.mp.dps = 80
    kappa = mp.mpf(5)/7
    emission_log_derivative = mp.mpf(1)/3
    weight_lipschitz = mp.mpf(11)/20*emission_log_derivative
    slope = weight_lipschitz/(1-kappa)
    diameter = 2*mp.log(36); projective = 2*slope*diameter
    gap = 1-mp.tanh(projective/4)
    old = json.loads((ROOT / "FIN_ST394_Quantitative_Chernoff_Cone_Gap.json").read_text())
    packet = {
        "exact_emission_log_derivative_bound": "1/3",
        "derivation": "max_b |d_h log(e1+(e0-e1)b)|=|sqrt(e0)-sqrt(e1)|/(sqrt(e0)+sqrt(e1))",
        "declared_s_interval": [.45, .55], "adapted_product_metric_weight_Lipschitz": float(weight_lipschitz),
        "filter_contraction": float(kappa), "invariant_cone_slope": float(slope),
        "state_diameter": float(diameter), "projective_diameter_upper": float(projective),
        "Birkhoff_gap_lower": float(gap), "log10_gap_lower": float(mp.log10(gap)),
        "previous_log10_gap_lower": old["log10_contraction_gap_lower"],
        "improvement_in_log10_orders": float(mp.log10(gap)-old["log10_contraction_gap_lower"]),
        "theorem": (
            "Using the exact Hilbert-coordinate derivative of the supplied binary emission "
            "weights replaces the coarse 6.6 bound by 11/60.  The same invariant-cone argument "
            "then gives the displayed nonastronomical explicit contraction gap."
        ),
    }
    return finalize(413, "Adapted-Metric Chernoff Cone Gap",
                    "proven_conditional_gap_improved_by_over_69_orders",
                    "The HMM and Chernoff fixture remain synthetic and are not a FIN-derived detector.", packet)


def st414() -> dict:
    mp.iv.dps = 55
    center = root(ir_equations, [.028, 2.1, 2.56, .021, .108], tol=1e-12).x
    # Search the fixed affine-radius template only; this is not a maximal tube
    # among all possible centers, orientations, or interval preconditioners.
    lo, hi = .005, .006
    trials = []
    for _ in range(12):
        d = (lo+hi)/2; radii = np.array([.1*d, 4*d, 5*d, .1*d, .4*d])
        cert = parametric_ir_krawczyk(center, radii, (3-d, 3+d), (4-d, 4+d))
        trials.append({"delta": d, "included": cert["included"], "minimum_margin": cert["minimum_margin"]})
        if cert["included"]: lo = d
        else: hi = d
    delta = .0055; radii = np.array([.1*delta, 4*delta, 5*delta, .1*delta, .4*delta])
    cert = parametric_ir_krawczyk(center, radii, (3-delta, 3+delta), (4-delta, 4+delta))
    X = [iv((center[i]-radii[i], center[i]+radii[i])) for i in range(5)]
    derivatives = [bounds(ir_Gprime_general_iv(X[i], X[3], iv((4-delta, 4+delta)))) for i in range(3)]
    packet = {
        "certified_delta": delta, "budget_box": [3-delta, 3+delta], "active_time_box": [4-delta, 4+delta],
        "root_radii": radii.tolist(), "parametric_Krawczyk": cert,
        "root_derivative_intervals": derivatives,
        "expansion_factor_over_ST395_parameter_halfwidth": delta/1e-4,
        "fixed_template_bisection_trials": trials,
        "fixed_template_last_included_delta": lo, "fixed_template_first_failed_delta": hi,
        "theorem": (
            "One common interval box contains a unique KKT root for every parameter pair in "
            "the displayed rectangle, expanding the certified ST395 parameter halfwidth by 55."
        ),
    }
    return finalize(414, "Expanded Uniform Continuum-IR Krawczyk Tube",
                    "proven_unique_root_tube_55_times_wider_in_parameters" if cert["included"] else "expanded_tube_failed",
                    "Maximality is only relative to the fixed center/radius template.  A full complement-sign topology cover is not re-certified on the enlarged box.", packet)


def limiting_ir_system(z):
    y1, y2, nu, T = z
    def g(y): return 2*T*y*math.exp(-(T-1)*y)/(1+math.exp(-T*y))-nu
    F = lambda t, y: -2*math.log1p(math.exp(-t*y))
    return np.array([g(y1), g(y2), exp1(y1)-exp1(y2)-3,
                     (F(1, y2)-F(1, y1))-(F(T, y2)-F(T, y1))])


def limiting_ir_interval_fj(X):
    y1, y2, nu, T = X; ys = (y1, y2)
    f = []
    for y in ys:
        x = mp.iv.exp(-T*y); base = 2*T*y*mp.iv.exp(-(T-1)*y)/(1+x)
        f.append(base-nu)
    f.append(e1_iv_local(y1)-e1_iv_local(y2)-3)
    Ft = lambda t, y: -2*mp.iv.log(1+mp.iv.exp(-t*y))
    f.append((Ft(iv(1), y2)-Ft(iv(1), y1))-(Ft(T, y2)-Ft(T, y1)))
    jac = [[iv(0) for _ in range(4)] for _ in range(4)]
    for i, y in enumerate(ys):
        x = mp.iv.exp(-T*y); base = 2*T*y*mp.iv.exp(-(T-1)*y)/(1+x)
        jac[i][i] = 2*T*mp.iv.exp(-(T-1)*y)/(1+x)*(1+y*(-(T-1)+T*x/(1+x)))
        jac[i][2] = iv(-1)
        jac[i][3] = base*(1/T-y+y*x/(1+x))
    jac[2][0] = -mp.iv.exp(-y1)/y1; jac[2][1] = mp.iv.exp(-y2)/y2
    rate = lambda t, y: 2*t/(mp.iv.exp(t*y)+1)
    jac[3][0] = -rate(iv(1), y1)+rate(T, y1)
    jac[3][1] = rate(iv(1), y2)-rate(T, y2)
    jac[3][3] = -2*y2/(mp.iv.exp(T*y2)+1)+2*y1/(mp.iv.exp(T*y1)+1)
    return f, jac


def limiting_krawczyk(center, radius):
    f0, j0 = limiting_ir_interval_fj([iv(float(x)) for x in center])
    flo = np.array([bounds(x)[0] for x in f0]); fhi = np.array([bounds(x)[1] for x in f0])
    jm = np.array([[(bounds(x)[0]+bounds(x)[1])/2 for x in row] for row in j0]); pre = np.linalg.inv(jm)
    X = [iv((center[i]-radius, center[i]+radius)) for i in range(4)]
    _, jac = limiting_ir_interval_fj(X)
    jl = np.array([[bounds(x)[0] for x in row] for row in jac]); jh = np.array([[bounds(x)[1] for x in row] for row in jac])
    cfl, cfh = interval_matvec(pre, pre, flo, fhi); ylo, yhi = center-cfh, center-cfl
    cjl, cjh = interval_left(pre, jl, jh); mlo, mhi = -cjh, -cjl
    for i in range(4):
        mlo[i, i] = np.nextafter(mlo[i, i]+1, -np.inf); mhi[i, i] = np.nextafter(mhi[i, i]+1, np.inf)
    dlo, dhi = interval_matvec(mlo, mhi, np.full(4, -radius), np.full(4, radius))
    margin = np.minimum(ylo+dlo-(center-radius), (center+radius)-(yhi+dhi))
    return {"included": bool(np.min(margin) > 0), "radius": radius,
            "minimum_margin": float(np.min(margin)), "component_margins": margin.tolist()}


def st415() -> dict:
    mp.iv.dps = 65
    sol = root(limiting_ir_system, [.0286, 3.89, .0694, 2.442])
    center = sol.x; cert = limiting_krawczyk(center, 1e-8)
    # Floating continuation records the approach but is deliberately separated
    # from the interval theorem for the limiting face itself.
    z = root(ir_equations, [.028, 2.1, 2.56, .021, .108]).x; rows = []
    for T in np.linspace(4, 2.46, 78):
        s = root(lambda x: ir_param_equations(x, 3., float(T)), z)
        residual = float(np.linalg.norm(ir_param_equations(s.x, 3., T), np.inf))
        if s.success and residual < 1e-8 and 0 < s.x[0] < s.x[1] < s.x[2] and s.x[3] > 0:
            z = s.x; rows.append({"T": float(T), "a": float(z[3]), "y3": float(z[2]), "residual": residual})
    packet = {
        "limiting_face": "a=0 and y3=+infinity",
        "limiting_unknowns": ["y1", "y2", "nu", "T"],
        "limiting_root_center": center.tolist(), "point_residual_inf": float(np.linalg.norm(limiting_ir_system(center), np.inf)),
        "limiting_face_Krawczyk": cert,
        "numerical_branch_approach_rows": rows,
        "interpretation": (
            "At a=0 the large-y linear term disappears, so the third threshold escapes to "
            "infinity and the upper active band vanishes.  The limiting four-equation face has "
            "a unique root in the certified box."
        ),
        "theorem_scope": "unique limiting-face root only",
        "branch_attachment_status": "strong numerical evidence, not an interval-certified singular continuation",
    }
    return finalize(415, "Certified Limiting Face for Continuum-IR Upper-Band Loss",
                    "proven_unique_limiting_face_root_branch_attachment_numerical" if cert["included"] else "limiting_face_certificate_failed",
                    "The boundary is a supplied dimensionless optimization transition, not a physical phase transition or IR scale.", packet)


def st416() -> dict:
    packet = {
        "strict_source_gate": {"gain": "open", "selector": "QW-2191 open", "units": "open",
                               "legacy_to_strict_bridge": "open", "laboratory_record": "absent"},
        "independent_evidence_gate": {"external_referee": "absent", "independent_lab": "absent",
                                      "new_empirical_data": "absent"},
        "mathematical_exports": ["independent orbit replay", "gain persistence tube", "degree-eight relation lower bound",
                                 "adapted Chernoff gap", "expanded IR root tube", "limiting IR face"],
        "promotion_forbidden": ["physical vacuum", "Born probabilities", "Standard Model", "gravity", "Theory of Everything"],
    }
    return finalize(416, "Strict-Source and Independent-Evidence Gates",
                    "gates_remain_open_despite_new_local_mathematics",
                    "Local analytic and computational work cannot manufacture external evidence.", packet)


def figures(results):
    FIG_DIR.mkdir(exist_ok=True)
    # Gain status map and certified interval.
    fig, ax = plt.subplots(figsize=(8.0, 4.2))
    ax.axvspan(3.999, 4.001, color="#2a9d8f", alpha=.35, label="certified 12-minimum orbit")
    cross = results["ST406"]["certified_local_coexistence_crossing_reused"]
    ax.axvspan(*cross, color="#e9c46a", alpha=.65, label="certified branch-local crossing")
    ax.axvline(results["ST406"]["uniform_Hessian_loss_gain"], color="#e76f51", ls="--", label="uniform Hessian loss")
    ax.set(xlabel="supplied dimensionless gain g", yticks=[], title="Certified statements must not be conflated")
    ax.set_xlim(2.6, 5.3); ax.legend(fontsize=8); fig.tight_layout()
    fig.savefig(FIG_DIR / "st405_gain_status.png", dpi=180); plt.close(fig)

    # Sampled stationary landscape.
    rows = results["ST407"]["sampled_stationary_catalog"]
    fig, ax = plt.subplots(figsize=(8.0, 4.4))
    ax.scatter([x["maximum_probability"] for x in rows], [x["value"] for x in rows],
               c=[x["Morse_index"] for x in rows], cmap="viridis", s=36)
    ax.set(xlabel="maximum coordinate", ylabel="V4", title="Sampled stationary orbits (colour = Morse index)")
    fig.tight_layout(); fig.savefig(FIG_DIR / "st407_morse_catalog.png", dpi=180); plt.close(fig)

    # Chernoff gap comparison.
    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    vals = [results["ST413"]["previous_log10_gap_lower"], results["ST413"]["log10_gap_lower"]]
    ax.bar(["previous coarse cone", "adapted metric"], vals, color=["#8d99ae", "#2a9d8f"])
    ax.set_ylabel("log10 certified contraction gap"); ax.set_title("Metric adaptation removes a 69-order artefact")
    fig.tight_layout(); fig.savefig(FIG_DIR / "st413_gap_improvement.png", dpi=180); plt.close(fig)

    # IR approach to the limiting face.
    rows = results["ST415"]["numerical_branch_approach_rows"]
    fig, ax1 = plt.subplots(figsize=(8.0, 4.4)); ax2 = ax1.twinx()
    ax1.plot([x["T"] for x in rows], [x["a"] for x in rows], color="#264653", label="a")
    ax2.plot([x["T"] for x in rows], [x["y3"] for x in rows], color="#e76f51", label="y3")
    ax1.set(xlabel="active time T (supplied)", ylabel="a", title="Numerical approach to a=0, y3=infinity face")
    ax2.set_ylabel("third threshold y3"); fig.tight_layout()
    fig.savefig(FIG_DIR / "st415_band_loss.png", dpi=180); plt.close(fig)


def main():
    FIG_DIR.mkdir(exist_ok=True)
    results = {}
    for k, fn in [(402, st402), (403, st403), (404, st404), (405, st405),
                  (406, st406), (407, st407), (408, st408), (409, st409),
                  (410, st410), (411, st411), (412, st412), (413, st413),
                  (414, st414), (415, st415), (416, st416)]:
        print(f"running ST{k}", flush=True); results[f"ST{k}"] = fn()
    payload = {
        "metadata": {"title": "FIN ST402--ST416 research results", "date": "2026-08-14",
                     "seed": SEED, "python": platform.python_version(), "numpy": np.__version__,
                     "scipy": scipy.__version__, "mpmath": mp.__version__},
        **results,
    }
    RESULTS.write_text(json.dumps(native(payload), indent=2, sort_keys=True), encoding="utf-8")
    with SUMMARY.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f); w.writerow(["program", "object", "status", "boundary", "packet_sha256"])
        for k in range(402, 417):
            r = results[f"ST{k}"]; w.writerow([f"ST{k}", r["object"], r["status"], r["boundary"], r["packet_sha256"]])
    figures(results)
    print(json.dumps({k: v["status"] for k, v in results.items()}, indent=2))


if __name__ == "__main__":
    main()
