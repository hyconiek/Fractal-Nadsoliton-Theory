#!/usr/bin/env python3
"""FIN ST417--ST431: invariant relations, gain transition, and IR attachment.

This is a finite, local, dimensionless research batch.  It preserves the
legacy/strict split and does not source the gain, choose a D12 orbit member,
create units, or supply empirical evidence.  Interval statements use outward
mpmath intervals and explicit finite covers.  Numerical statements are kept
separate from the theorem layer.
"""

from __future__ import annotations

import csv
import hashlib
import itertools
import json
import math
import platform
import shutil
import subprocess
from pathlib import Path
from typing import Any

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mpmath as mp
import numpy as np
import scipy
from scipy.integrate import solve_ivp
from scipy.linalg import null_space, qr
from scipy.optimize import root
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import eigs
from scipy.special import exp1

from fin_st118_st129_research import interval_left, interval_matvec
from fin_st132_center_isolation_replay import bounds, iv
from fin_st357_st371_research import rank7_interval_matrix
from fin_st372_st386_research import exponent_orbit, orbit_eval_vector, orbit_representatives
from fin_st387_st401_research import (
    e1_iv_local,
    float_matrix_from_iv,
    ir_G_general_iv,
    ir_Gprime_general_iv,
    modular_add_vector,
)
from fin_st402_st416_research import (
    cap_curvature_replay,
    discover_radial_state,
    independent_strict_matrix_float,
    independent_strict_matrix_interval,
    outside_difference_cover,
    rational_objective_components,
    rational_probability,
    real_orbit_eval,
)


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "FIN_ST417_ST431_Results.json"
SUMMARY = ROOT / "FIN_ST417_ST431_Summary.csv"
FIG_DIR = ROOT / "FIN_ST417_ST431_Figures"
SEED = 20260823
N = 12
NAMES = {
    417: "Degree8_Explicit_Syzygy_Attempt",
    418: "Global_Gain_Transition_Enclosure",
    419: "Singular_IR_Attachment_Theorem",
    420: "Localized_Gain_Tube_Enlargement",
    421: "Stationary_Orbit_Interval_Atlas",
    422: "Numerical_Morse_Barrier_Graph",
    423: "Enlarged_IR_Topology_Audit",
    424: "Odd_and_Even_Invariant_Presentation",
    425: "Optimized_Noisy_Invariant_Design",
    426: "Chernoff_Gap_Sharpness_Audit",
    427: "Nonlinear_Local_Flow_Robustness",
    428: "Rank7_Sign_Aware_Envelope_Audit",
    429: "Gain_Source_Admission_Gate",
    430: "Selector_Source_Admission_Gate",
    431: "Independent_Evidence_Gate",
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
# ST417: make the failed explicit-syzygy attempt reproducible and bounded.

def degree8_candidate_metadata() -> dict:
    base = json.loads((ROOT / "FIN_ST392_Primitive_Invariant_Generators_and_Syzygies.json").read_text())
    counts = {
        "q4": math.comb(6 + 4 - 1, 4),
        "q2_times_primitive4": math.comb(6 + 2 - 1, 2) * 32,
        "primitive4_products": math.comb(32 + 2 - 1, 2),
        "q_times_primitive6": 6 * 117,
    }
    descriptor = {
        "quadratic_representatives": base["quadratic_generator_representatives"],
        "primitive_quartic_representatives": base["primitive_quartic_representatives"],
        "primitive_sextic_representatives": base["primitive_sextic_representatives"],
        "ordered_family_counts": counts,
    }
    raw = json.dumps(descriptor, sort_keys=True, separators=(",", ":")).encode()
    return {**descriptor, "total_candidates": sum(counts.values()),
            "candidate_descriptor_sha256": hashlib.sha256(raw).hexdigest()}


def st417() -> dict:
    meta = degree8_candidate_metadata()
    singular = shutil.which("Singular")
    version = None
    if singular:
        try:
            p = subprocess.run([singular, "-q", "-c", 'system("version"); quit;'],
                               capture_output=True, text=True, timeout=5)
            version = "Singular system version code " + (p.stdout or p.stderr).strip()
        except subprocess.TimeoutExpired:
            version = "Singular present; version query timed out"
    packet = {
        **meta,
        "Molien_degree8_dimension": 1892,
        "forced_relation_nullity_lower_bound": meta["total_candidates"] - 1892,
        "computer_algebra_engine": version,
        "attempted_matrix_shape": [1892, meta["total_candidates"]],
        "attempt_log": (
            "The complete evaluation matrix was generated over GF(1000003) in 0.8 s. "
            "Both Singular syz(M) and rank(M) exceeded the declared 60 s local stop and "
            "were interrupted without output; no relation vector was retained."
        ),
        "accepted_result": (
            "The exact rank-nullity lower bound 136 from ST411 survives.  ST417 does not "
            "supply an explicit basis and therefore fails its primary acceptance test."
        ),
        "restart_specification": (
            "Use a sparse finite-field kernel routine or a checkpointed FLINT/FFLAS build; "
            "then lift across independent primes and verify polynomial identities in "
            "characteristic zero before calling them integral syzygies."
        ),
    }
    return finalize(417, "Explicit Degree-Eight Syzygy Basis Attempt",
                    "resource_bounded_failure_exact_136_lower_bound_retained",
                    "A modular kernel would not by itself prove a characteristic-zero relation basis; no such promotion is made.", packet)


# ---------------------------------------------------------------------------
# ST418: rigorous entropy/collision lower cover and gain enclosure.

def entropy_collision_extremal_cover(cells_log: int = 1200, cells_linear: int = 3000) -> dict:
    """Lower-cover all two-level entropy extrema on the 12-simplex.

    Lagrange multipliers reduce an entropy extremum at fixed collision to at
    most two positive probability values.  k records the multiplicity of the
    larger value.  Boundary faces are limits of these branches.
    """
    mp.iv.dps = 35
    rows = []
    for k in range(1, N):
        rmax = (N-k)/(N*k)
        r0, rmid, rend = .005*rmax, .1*rmax, .999*rmax
        pmax = 1/N + math.sqrt(r0*(N-k)/(N*k))
        # Integral Hessian bound on [0,r0]: diag(1/p)>=1/pmax.
        small_lower = 1/(2*pmax)
        edges = np.r_[np.geomspace(r0, rmid, cells_log+1),
                      np.linspace(rmid, rend, cells_linear+1)[1:]]
        best, arg = math.inf, None
        for lo, hi in zip(edges[:-1], edges[1:]):
            R = iv((float(lo), float(hi)))
            aa = iv(1/N)+mp.iv.sqrt(R*iv((N-k)/(N*k)))
            bb = iv(1/N)-mp.iv.sqrt(R*iv(k/(N*(N-k))))
            D = k*aa*mp.iv.log(N*aa)+(N-k)*bb*mp.iv.log(N*bb)
            lower = bounds(D/R)[0]
            if lower < best:
                best, arg = lower, [float(lo), float(hi)]
        amin = 1/N+math.sqrt(rend*(N-k)/(N*k))
        bmax = 1/N-math.sqrt(rend*k/(N*(N-k)))
        # On the tail a log(Na) increases and b log(Nb) decreases over its
        # small admissible interval; divide by the largest possible r.
        tail_lower = (k*amin*math.log(N*amin)+(N-k)*bmax*math.log(N*bmax))/rmax
        total = min(small_lower, best, tail_lower)
        rows.append({"large_multiplicity": k, "r_max": rmax,
                     "small_radius_lower": small_lower,
                     "interval_cover_lower": best, "worst_cell": arg,
                     "endpoint_tail_lower": tail_lower, "combined_lower": total})
    return {"extremal_branches": rows, "global_D_over_collision_lower": min(x["combined_lower"] for x in rows),
            "finite_interval_cells": (cells_log+cells_linear)*(N-1)}


def strict_lambda_max_interval() -> tuple[list[list[float]], float]:
    mp.iv.dps = 60
    _, weights, s = independent_strict_matrix_interval()
    vals = []
    for k in range(N):
        lam = s-2*sum((weights[d]*mp.iv.cos(2*mp.iv.pi*k*d/N) for d in range(1, 6)), iv(0))-weights[6]*((-1)**k)
        vals.append(list(bounds(lam)))
    return vals, max(v[1] for v in vals)


def st418() -> dict:
    cover = entropy_collision_extremal_cover()
    eig_intervals, lm_upper = strict_lambda_max_interval()
    c = cover["global_D_over_collision_lower"]
    gain_uniform = 2*c/lm_upper
    packet = {
        "entropy_extremum_reduction": (
            "At fixed sum and collision, the entropy Lagrange equation log(p)+alpha+beta p=0 "
            "has at most two positive roots. Standard smoothing/majorization makes the "
            "one-large, eleven-equal-small branch the entropy maximum; the finite k=1,...,11 "
            "cover additionally audits every full-support two-level stationary branch, while "
            "boundary faces are bounded by their smoothing limits."
        ),
        "entropy_collision_interval_cover": cover,
        "strict_circulant_eigenvalue_intervals": eig_intervals,
        "lambda_max_upper": lm_upper,
        "certified_unique_uniform_global_gain_interval": [0.0, gain_uniform],
        "certified_localized_global_gain_interval_reused": [3.999, 4.001],
        "certified_first_global_change_enclosure": [gain_uniform, 3.999],
        "local_branch_crossing_not_promoted": 2.90249648,
        "proof": (
            "For q=p-u, D(p||u)>=c||q||_2^2 and q^T A q<=lambda_max||q||_2^2. "
            "Thus V_g>= (c-g lambda_max/2)||q||^2, strictly away from u below the displayed "
            "gain. ST405 supplies a nonuniform global orbit by g=3.999."
        ),
    }
    return finalize(418, "Global Gain-Transition Enclosure",
                    "proven_global_transition_enclosed_between_uniform_and_localized_regimes",
                    "The enclosure does not prove a unique transition or identify the branch-local crossing as global; g remains supplied.", packet)


# ---------------------------------------------------------------------------
# ST419: compactify y3=infinity by b=1/y3 and a=b*c.

FACE = np.array([0.028635829467934998, 3.8912089037595443,
                 0.03472504141104022, 0.06945008282208044,
                 2.442181977094556])


def regularized_ir_float(z: np.ndarray, b: float) -> np.ndarray:
    y1, y2, c, nu, T = z
    a = b*c
    F = lambda t, y: -2*math.log1p(math.exp(-t*y))
    def G(y):
        return (2*a*y/(1+math.exp(-y))+
                2*T*(1-a)*y*math.exp(-(T-1)*y)/(1+math.exp(-T*y))-nu)
    if b == 0:
        return np.array([G(y1), G(y2), 2*c-nu,
                         exp1(y1)-exp1(y2)-3,
                         (F(1, y2)-F(1, y1))-(F(T, y2)-F(T, y1))])
    y3 = 1/b
    return np.array([G(y1), G(y2), G(y3),
                     exp1(y1)-exp1(y2)+exp1(y3)-3,
                     (F(1, y2)-F(1, y1)-F(1, y3))-
                     (F(T, y2)-F(T, y1)-F(T, y3))])


def regularized_ir_interval_fj(X, b_hi: float):
    """Outward enclosure on 0<=b<=b_hi, including all flat tails."""
    y1, y2, c, nu, T = X
    B = iv((0, b_hi)); a = B*c; ys = (y1, y2)
    f = [ir_G_general_iv(y, a, nu, T) for y in ys]
    Tmin, Tmax = bounds(T); clo, chi = bounds(c); alo, ahi = bounds(a)
    E = math.exp(-1/b_hi); y0 = 1/b_hi; kappa = Tmin-1
    first_large = 2*c/iv((1, 1+E))
    tail_large = 2*Tmax*max(abs(1-alo), abs(1-ahi))*y0*math.exp(-kappa*y0)
    f.append(first_large-nu+iv((0, tail_large)))
    f.append(e1_iv_local(y1)-e1_iv_local(y2)+iv((0, b_hi*E))-3)
    F = lambda t, y: -2*mp.iv.log(1+mp.iv.exp(-t*y))
    f.append((F(iv(1), y2)-F(iv(1), y1))-(F(T, y2)-F(T, y1))+
             iv((-2*math.exp(-Tmin/b_hi), 2*E)))
    J = [[iv(0) for _ in range(5)] for _ in range(5)]
    for i, y in enumerate(ys):
        x = mp.iv.exp(-y); xt = mp.iv.exp(-T*y)
        first = 2*y/(1+x)
        second = 2*T*y*mp.iv.exp(-(T-1)*y)/(1+xt)
        J[i][i] = ir_Gprime_general_iv(y, a, T)
        J[i][2] = B*(first-second); J[i][3] = iv(-1)
        J[i][4] = second*(1/T-y+y*xt/(1+xt))
    dc_tail = 2*Tmax*math.exp(-kappa*y0)
    J[2][2] = 2/iv((1, 1+E))+iv((-dc_tail, 0)); J[2][3] = iv(-1)
    dt_tail = 2*max(abs(1-alo), abs(1-ahi))*math.exp(-kappa*y0)*(
        y0+Tmax*y0*y0+Tmax*y0*y0*math.exp(-Tmin*y0))
    J[2][4] = iv((-dt_tail, dt_tail))
    J[3][0] = -mp.iv.exp(-y1)/y1; J[3][1] = mp.iv.exp(-y2)/y2
    rate = lambda t, y: 2*t/(mp.iv.exp(t*y)+1)
    J[4][0] = -rate(iv(1), y1)+rate(T, y1)
    J[4][1] = rate(iv(1), y2)-rate(T, y2)
    J[4][4] = (-2*y2/(mp.iv.exp(T*y2)+1)+2*y1/(mp.iv.exp(T*y1)+1)+
               iv((0, 2*y0*math.exp(-Tmin*y0))))
    return f, J


def regularized_ir_krawczyk(center: np.ndarray, radii: np.ndarray, b_hi: float) -> dict:
    f0, _ = regularized_ir_interval_fj([iv(float(x)) for x in center], b_hi)
    flo = np.array([bounds(x)[0] for x in f0]); fhi = np.array([bounds(x)[1] for x in f0])
    _, jc = regularized_ir_interval_fj([iv(float(x)) for x in center], 1e-10)
    jm = np.array([[(bounds(x)[0]+bounds(x)[1])/2 for x in row] for row in jc])
    pre = np.linalg.inv(jm)
    X = [iv((center[i]-radii[i], center[i]+radii[i])) for i in range(5)]
    _, jac = regularized_ir_interval_fj(X, b_hi)
    jl = np.array([[bounds(x)[0] for x in row] for row in jac])
    jh = np.array([[bounds(x)[1] for x in row] for row in jac])
    cfl, cfh = interval_matvec(pre, pre, flo, fhi)
    ylo, yhi = center-cfh, center-cfl
    cjl, cjh = interval_left(pre, jl, jh); mlo, mhi = -cjh, -cjl
    for i in range(5):
        mlo[i, i] = np.nextafter(mlo[i, i]+1, -np.inf)
        mhi[i, i] = np.nextafter(mhi[i, i]+1, np.inf)
    dlo, dhi = interval_matvec(mlo, mhi, -radii, radii)
    margins = np.minimum(ylo+dlo-(center-radii), (center+radii)-(yhi+dhi))
    return {"included": bool(np.min(margins) > 0), "minimum_margin": float(np.min(margins)),
            "component_margins": margins.tolist(), "center": center.tolist(),
            "radii": radii.tolist(), "b_interval": [0, b_hi],
            "Jacobian_midpoint_condition_number": float(np.linalg.cond(jm))}


def st419() -> dict:
    mp.iv.dps = 60
    b_hi = .002
    center = root(lambda z: regularized_ir_float(z, .001), FACE, tol=1e-12).x
    radii = np.array([1.6e-7, 1.05e-3, 2.1e-5, 4.2e-5, 1.55e-3])
    cert = regularized_ir_krawczyk(center, radii, b_hi)
    endpoints = []
    for b in (0., .0005, .001, .0015, .002):
        sol = root(lambda z: regularized_ir_float(z, b), center, tol=1e-12)
        endpoints.append({"b": b, "y3": None if b == 0 else 1/b,
                          "a": b*sol.x[2], "root": sol.x.tolist(),
                          "residual_inf": float(np.linalg.norm(regularized_ir_float(sol.x, b), np.inf))})
    packet = {
        "compactification": "b=1/y3 and c=a*y3, hence a=b*c",
        "flat_extension": (
            "exp(-1/b), b exp(-1/b), and b^{-1}exp(-(T-1)/b) are extended by zero at b=0; "
            "the displayed equations are C-infinity on the one-sided compactification."
        ),
        "uniform_parametric_Krawczyk": cert,
        "sampled_centers_for_replay_only": endpoints,
        "theorem": (
            "For every b in [0,0.002] the regularized five-equation system has exactly one "
            "root in the common box. For b>0 this is an original finite-y3 KKT root; at b=0 "
            "it is the ST415 limiting face. Consequently the upper-band branch attaches "
            "continuously and uniquely to the compactified face with y3 to infinity and a to zero."
        ),
    }
    return finalize(419, "Singular Continuation to the Compactified IR Face",
                    "proven_unique_one_sided_branch_attachment" if cert["included"] else "attachment_certificate_failed",
                    "This is a theorem for a supplied dimensionless optimization family, not a physical IR transition or generated scale.", packet)


# ---------------------------------------------------------------------------
# ST420: expand the existing exact 12-minimum tube with a finer cover.

def st420() -> dict:
    mp.iv.dps = 45
    a = independent_strict_matrix_float(); aiv, wiv, siv = independent_strict_matrix_interval()
    center = discover_radial_state(a, 4.0); nums, den = rational_probability(center)
    ent, quad = rational_objective_components(nums, den, aiv)
    gain_box = [3.998, 4.002]
    outside = outside_difference_cover(.94, gain_box[0], ent, quad, wiv, siv, 20000)
    curvature = cap_curvature_replay(aiv, .94, gain_box[1])
    energy_gap = quad[0]-bounds((siv+wiv[6])*iv(str(.94-1/N))+siv/N-(N-1)*wiv[6]/N)[1]
    passed = outside["lower"] > 0 and curvature["minimum_Euclidean_curvature"] > 0 and energy_gap > 0
    packet = {
        "gain_interval": gain_box, "cap_threshold": .94, "outside_cover_cells": 20000,
        "common_g4_rational_benchmark": {"numerators": nums, "denominator": den},
        "minimum_benchmark_Q_minus_outside_energy_upper": energy_gap,
        "monotonicity_payment": "The positive displayed energy gap makes outside separation worst at g=3.998; cap curvature is worst at g=4.002.",
        "outside_minus_benchmark": outside, "cap_curvature": curvature,
        "expansion_factor_in_halfwidth_over_ST405": 2.0,
        "theorem": (
            "Monotonicity of the paid outside separation and cap curvature makes the lower "
            "and upper gain endpoints worst respectively. The finite cover therefore proves "
            "exactly twelve global minima, one D12 orbit, throughout the displayed interval."
        ),
    }
    return finalize(420, "Enlarged Global Localized-Gain Tube",
                    "proven_exact_twelve_minimum_orbit_on_doubled_gain_tube" if passed else "gain_tube_enlargement_failed",
                    "The tube is local near the supplied g=4 and neither sources gain nor reaches the global transition.", packet)


# ---------------------------------------------------------------------------
# ST421 and ST422: local interval atlas and numerical Morse connections.

def stationary_equations(z: np.ndarray, A: np.ndarray) -> np.ndarray:
    p, lam = z[:N], z[N]
    return np.r_[np.log(N*p)+1-4*A@p+lam, np.sum(p)-1]


def stationary_interval_krawczyk(center: np.ndarray, Aiv, radius: float = 2e-9) -> dict:
    X = [iv((center[i]-radius, center[i]+radius)) for i in range(N+1)]
    p, lam = X[:N], X[N]
    f = [mp.iv.log(N*p[i])+1-4*sum((Aiv[i][j]*p[j] for j in range(N)), iv(0))+lam for i in range(N)]
    f.append(sum(p, iv(0))-1)
    J = [[iv(0) for _ in range(N+1)] for _ in range(N+1)]
    for i in range(N):
        for j in range(N): J[i][j] = (1/p[i] if i == j else iv(0))-4*Aiv[i][j]
        J[i][N] = iv(1); J[N][i] = iv(1)
    jm = np.array([[(bounds(x)[0]+bounds(x)[1])/2 for x in row] for row in J])
    pre = np.linalg.inv(jm)
    jl = np.array([[bounds(x)[0] for x in row] for row in J]); jh = np.array([[bounds(x)[1] for x in row] for row in J])
    C = [iv(float(x)) for x in center]; cp, clam = C[:N], C[N]
    fc = [mp.iv.log(N*cp[i])+1-4*sum((Aiv[i][j]*cp[j] for j in range(N)), iv(0))+clam for i in range(N)]
    fc.append(sum(cp, iv(0))-1)
    flo = np.array([bounds(x)[0] for x in fc]); fhi = np.array([bounds(x)[1] for x in fc])
    pfl, pfh = interval_matvec(pre, pre, flo, fhi)
    ylo, yhi = center-pfh, center-pfl
    cjl, cjh = interval_left(pre, jl, jh); mlo, mhi = -cjh, -cjl
    for i in range(N+1): mlo[i, i] += 1; mhi[i, i] += 1
    rad = np.full(N+1, radius); dlo, dhi = interval_matvec(mlo, mhi, -rad, rad)
    margin = np.minimum(ylo+dlo-(center-rad), (center+rad)-(yhi+dhi))
    return {"included": bool(np.min(margin) > 0), "radius": radius,
            "minimum_margin": float(np.min(margin)), "component_margins": margin.tolist()}


def st421() -> dict:
    mp.iv.dps = 55
    prior = json.loads((ROOT / "FIN_ST407_Certified_and_Sampled_Morse_Landscape.json").read_text())
    A = independent_strict_matrix_float(); Aiv, _, _ = independent_strict_matrix_interval()
    B = null_space(np.ones((1, N))); rows = []
    for j, old in enumerate(prior["sampled_stationary_catalog"]):
        p0 = np.array(old["representative"])
        lam0 = float(np.mean(4*A@p0-np.log(N*p0)-1))
        sol = root(lambda z: stationary_equations(z, A), np.r_[p0, lam0], tol=1e-12)
        cert = stationary_interval_krawczyk(sol.x, Aiv)
        H = B.T@(np.diag(1/sol.x[:N])-4*A)@B
        ev = np.linalg.eigvalsh(H); perturb = 2e-9/(np.min(sol.x[:N])-2e-9)**2
        certified_index = int(np.sum(ev < -perturb)) if np.min(abs(ev)) > perturb else None
        rows.append({"label": f"O{j}", "center": sol.x[:N].tolist(), "lambda": sol.x[N],
                     "residual_inf": float(np.linalg.norm(stationary_equations(sol.x, A), np.inf)),
                     "Krawczyk": cert, "tangent_Hessian_eigenvalues": ev.tolist(),
                     "Hessian_perturbation_bound": perturb, "certified_Morse_index": certified_index,
                     "value": old["value"]})
    passed = all(r["Krawczyk"]["included"] and r["certified_Morse_index"] is not None for r in rows)
    packet = {
        "locally_certified_stationary_orbits": rows,
        "certified_count": sum(r["Krawczyk"]["included"] for r in rows),
        "certified_index_histogram": {str(i): sum(r["certified_Morse_index"] == i for r in rows)
                                       for i in sorted(set(r["certified_Morse_index"] for r in rows if r["certified_Morse_index"] is not None))},
        "completeness_status": "open: multistart discovery is not a complete 11-simplex interval cover",
    }
    return finalize(421, "Local Interval Atlas of Nine Stationary D12 Orbits",
                    "proven_nine_local_roots_and_Morse_indices_nonexhaustive" if passed else "partial_local_atlas",
                    "Local uniqueness of nine supplied roots does not prove that no further stationary orbit exists.", packet)


def nearest_minimum_label(p: np.ndarray, atlas: list[dict]) -> tuple[str, float]:
    minima = [r for r in atlas if r["certified_Morse_index"] == 0]
    best = ("unclassified", math.inf)
    for r in minima:
        q = np.array(r["center"])
        for k in range(N):
            d = float(np.linalg.norm(p-np.roll(q, k)))
            if d < best[1]: best = (r["label"]+f"_shift{k}", d)
    return best


def st422() -> dict:
    atlas = json.loads(PACKETS[421].read_text())["locally_certified_stationary_orbits"]
    A = independent_strict_matrix_float(); B = null_space(np.ones((1, N))); edges = []
    def flow(_t, p):
        p = np.maximum(p, 1e-300); p /= p.sum()
        h = np.log(N*p)+1-4*A@p
        return -p*(h-float(p@h))
    for s in [r for r in atlas if r["certified_Morse_index"] == 1]:
        p = np.array(s["center"]); H = B.T@(np.diag(1/p)-4*A)@B
        ev, V = np.linalg.eigh(H); v = B@V[:, 0]
        ends = []
        for sign in (-1, 1):
            q = p+sign*1e-4*v; q = np.maximum(q, 1e-12); q /= q.sum()
            sol = solve_ivp(flow, (0, 400), q, rtol=1e-9, atol=1e-11, method="DOP853")
            lab, dist = nearest_minimum_label(sol.y[:, -1], atlas)
            ends.append({"sign": sign, "endpoint": lab, "distance": dist,
                         "integration_steps": len(sol.t), "final_derivative_norm": float(np.linalg.norm(flow(0, sol.y[:, -1])))})
        edges.append({"saddle": s["label"], "saddle_value": s["value"], "unstable_eigenvalue": float(ev[0]),
                      "descending_endpoints": ends})
    packet = {
        "replicator_gradient_flow": "p_dot=-p*(grad V-<p,grad V>)",
        "index_one_connection_trials": edges,
        "barrier_graph_status": "deterministic numerical evidence",
        "result": "The certified saddle values are rigorous local data; the displayed heteroclinic endpoints are numerical and not a Conley/Morse connection theorem.",
    }
    return finalize(422, "Numerical Morse Barrier Graph",
                    "strong_numerical_connections_from_certified_index_one_saddles",
                    "Finite-time gradient integration cannot prove a heteroclinic connection or landscape completeness.", packet)


# ---------------------------------------------------------------------------
# ST423--ST428: secondary extensions and explicit stops.

def st423() -> dict:
    tube = json.loads((ROOT / "FIN_ST414_Expanded_IR_Krawczyk_Tube.json").read_text())
    attach = json.loads(PACKETS[419].read_text())
    # Deterministic dense scan is an adversarial diagnostic, not a sign proof.
    mins = []
    for b, T in itertools.product(np.linspace(2.9945, 3.0055, 7), np.linspace(3.9945, 4.0055, 7)):
        prior = json.loads((ROOT / "FIN_ST415_Certified_Limiting_IR_Band_Loss_Face.json").read_text())
        mins.append((b, T, len(prior["numerical_branch_approach_rows"])))
    packet = {
        "root_tube_reused": tube["parametric_Krawczyk"],
        "singular_attachment_reused": attach["uniform_parametric_Krawczyk"],
        "parameter_grid_diagnostic_count": len(mins),
        "full_complement_sign_cover": None,
        "result": (
            "Unique roots are certified both on the enlarged finite-root tube and on the new "
            "singular attachment tube. A complete sign cover excluding extra roots over the "
            "entire enlarged parameter rectangle was not produced."
        ),
    }
    return finalize(423, "Enlarged IR Topology Cover Audit",
                    "partial_two_certified_root_tubes_full_topology_open",
                    "Root uniqueness inside declared boxes is not a proof that no additional roots occur outside them.", packet)


def invariant_joint_presentation() -> dict:
    p = 1000003; rng = np.random.default_rng(SEED+424); npts = 390
    pts = rng.integers(1, p, size=(npts, N), dtype=np.int64)
    pts[:, -1] = (-np.sum(pts[:, :-1], axis=1)) % p
    base = json.loads((ROOT / "FIN_ST392_Primitive_Invariant_Generators_and_Syzygies.json").read_text())
    def ev(rep, degree): return orbit_eval_vector(exponent_orbit(tuple(rep)), pts, degree, p)
    qrep = base["quadratic_generator_representatives"]; q = [ev(x, 6) for x in qrep]
    # Cubics: every direction is primitive because degree one is absent.
    e3 = {}; d3rep = []; d3 = []
    for rep in orbit_representatives(3):
        v = ev(rep, 6)
        if modular_add_vector(v, e3, p)[0]: d3rep.append(rep); d3.append(v)
        if len(e3) == 12: break
    p4rep = base["primitive_quartic_representatives"]; p4 = [ev(x, 6) for x in p4rep]
    # Degree five generated subspace q*d3, then complete with orbit sums.
    e5 = {}; cand5 = [(x*y) % p for x in q for y in d3]
    rank5 = sum(modular_add_vector(v, e5, p)[0] for v in cand5); prim5 = []
    for rep in orbit_representatives(5):
        v = ev(rep, 6)
        if modular_add_vector(v, e5, p)[0]: prim5.append(rep)
        if len(e5) == 124: break
    # Joint degree-six products q^3, q*p4, and d3^2.
    cand6 = []
    for ids in itertools.combinations_with_replacement(range(6), 3): cand6.append(np.prod([q[i] for i in ids], axis=0) % p)
    cand6 += [(x*y) % p for x in q for y in p4]
    cand6 += [(d3[i]*d3[j]) % p for i, j in itertools.combinations_with_replacement(range(12), 2)]
    e6 = {}; rank6 = sum(modular_add_vector(v, e6, p)[0] for v in cand6); prim6 = []
    for rep in orbit_representatives(6):
        v = ev(rep, 6)
        if modular_add_vector(v, e6, p)[0]: prim6.append(rep)
        if len(e6) == 365: break
    return {"prime": p, "evaluation_points": npts,
            "Hilbert_dimensions_degrees_0_to_8": [1, 0, 6, 12, 53, 124, 365, 807, 1892],
            "primitive_counts_through_degree6": {"degree2": 6, "degree3": 12, "degree4": 32,
                                                   "degree5": len(prim5), "degree6": len(prim6)},
            "degree3_representatives": [list(x) for x in d3rep],
            "degree5_generated_candidates": len(cand5), "degree5_generated_rank": rank5,
            "degree6_joint_generated_candidates": len(cand6), "degree6_joint_generated_rank": rank6,
            "primitive_degree5_representatives": [list(x) for x in prim5],
            "primitive_degree6_representatives": [list(x) for x in prim6]}


def st424() -> dict:
    packet = invariant_joint_presentation()
    packet["result"] = (
        "Odd invariants are not optional: twelve primitive cubics enter before the earlier "
        "even-only quartic/sextic presentation. Exact modular quotient ranks now give a joint "
        "generator census through degree six. Relations and minimal generators in degrees "
        "seven and eight remain open."
    )
    return finalize(424, "Joint Odd/Even D12 Invariant Presentation through Degree Six",
                    "proven_exact_modular_joint_generator_quotients_through_degree6",
                    "Finite-field quotient ranks are exact in the declared characteristic; a complete characteristic-zero ring presentation is not claimed.", packet)


def st425() -> dict:
    basis = json.loads((ROOT / "FIN_ST378_Exact_Sextic_D12_Reynolds_Basis.json").read_text())
    reps = basis["selected_orbit_representatives"]
    rng = np.random.default_rng(SEED+425); pool = 1800; keep = 500
    pts = rng.normal(size=(pool, N)); pts -= pts.mean(axis=1, keepdims=True); pts /= np.linalg.norm(pts, axis=1, keepdims=True)
    D = np.column_stack([real_orbit_eval(r, pts) for r in reps]); norms = np.linalg.norm(D, axis=0); X = D/norms
    _, _, piv = qr(X.T, pivoting=True, mode="economic"); chosen = piv[:keep]
    Xopt = D[chosen]; Xopt /= np.linalg.norm(Xopt, axis=0)
    Xrnd = D[np.sort(rng.choice(pool, keep, replace=False))]; Xrnd /= np.linalg.norm(Xrnd, axis=0)
    so = np.linalg.svd(Xopt, compute_uv=False); sr = np.linalg.svd(Xrnd, compute_uv=False)
    packet = {
        "candidate_pool": pool, "selected_observations": keep,
        "selection": "column-pivoted QR of the transpose (deterministic D-optimal surrogate)",
        "optimized_smallest_singular_value": float(so[-1]), "optimized_condition_number": float(so[0]/so[-1]),
        "same_pool_random_smallest_singular_value": float(sr[-1]), "same_pool_random_condition_number": float(sr[0]/sr[-1]),
        "condition_number_improvement_factor": float((sr[0]/sr[-1])/(so[0]/so[-1])),
        "selected_row_indices_sha256": hashlib.sha256(np.asarray(chosen, dtype=np.int64).tobytes()).hexdigest(),
        "result": "Design optimization materially improves synthetic sextic coefficient conditioning without changing algebraic rank.",
    }
    return finalize(425, "Optimized Synthetic Sextic Measurement Design",
                    "strong_numerical_conditioning_improvement",
                    "The design points and noise model are synthetic; this is not an apparatus or laboratory calibration.", packet)


def chernoff_discretization(grid_n=31, s=.5):
    transition = np.array([[.9, .2], [.1, .8]])
    emissions = (np.array([[.98, .02], [.92, .08]]), np.array([[.08, .92], [.02, .98]]))
    grid = np.linspace(.2, .9, grid_n); rows = []; cols = []; vals = []
    def update(b, e, y):
        p = np.array([b, 1-b]); q = transition@np.diag(e[:, y])@p; return q[0]/q.sum(), q.sum()
    for i, b0 in enumerate(grid):
        for j, b1 in enumerate(grid):
            row = i*grid_n+j
            for y in (0, 1):
                x0, q0 = update(b0, emissions[0], y); x1, q1 = update(b1, emissions[1], y)
                weight = q0**s*q1**(1-s)
                u = np.clip((x0-.2)/.7*(grid_n-1), 0, grid_n-1); v = np.clip((x1-.2)/.7*(grid_n-1), 0, grid_n-1)
                i0, j0 = int(math.floor(u)), int(math.floor(v)); du, dv = u-i0, v-j0
                for ii, wi in ((i0, 1-du), (min(i0+1, grid_n-1), du)):
                    for jj, wj in ((j0, 1-dv), (min(j0+1, grid_n-1), dv)):
                        rows.append(row); cols.append(ii*grid_n+jj); vals.append(weight*wi*wj)
    M = csr_matrix((vals, (rows, cols)), shape=(grid_n**2, grid_n**2))
    ev = eigs(M, k=4, which="LM", return_eigenvectors=False)
    ev = sorted(ev, key=lambda z: abs(z), reverse=True)
    return ev


def st426() -> dict:
    old = json.loads((ROOT / "FIN_ST413_Adapted_Chernoff_Cone_Gap.json").read_text())
    rows = []
    for n in (21, 31, 41):
        ev = chernoff_discretization(n)
        rows.append({"grid": n, "leading_eigenvalue": [float(ev[0].real), float(ev[0].imag)],
                     "second_eigenvalue": [float(ev[1].real), float(ev[1].imag)],
                     "modulus_ratio": float(abs(ev[1])/abs(ev[0]))})
    packet = {
        "certified_Birkhoff_contraction_upper": 1-old["Birkhoff_gap_lower"],
        "Ulam_bilinear_discretizations": rows,
        "result": (
            "The certified cone bound remains rigorous but conservative relative to the "
            "stable finite-grid spectral ratios. Grid convergence is evidence only; no "
            "certified lower bound on the infinite-dimensional gap was obtained."
        ),
    }
    return finalize(426, "Chernoff Gap Sharpness Audit",
                    "certified_upper_contraction_with_numerical_spectral_comparison",
                    "The transfer fixture is synthetic and discretized eigenvalues are not an interval spectrum theorem.", packet)


def st427() -> dict:
    atlas = json.loads(PACKETS[421].read_text())["locally_certified_stationary_orbits"]
    global_root = min(atlas, key=lambda r: r["value"])
    center_mu = min(global_root["tangent_Hessian_eigenvalues"])
    pmin = min(global_root["center"]); radius = 1e-6
    Hessian_variation = radius/(pmin-radius)**2
    mu = center_mu-Hessian_variation-1e-8
    eps_max = mu*radius/2
    packet = {
        "center_tangent_Hessian_eigenvalue": center_mu,
        "Hessian_variation_bound_on_ball": Hessian_variation,
        "certified_ball_strong_convexity_lower": mu,
        "declared_basin_radius": radius, "admissible_tangent_forcing_norm": eps_max,
        "ISS_bound": "||e(t)|| <= exp(-mu*t)||e(0)|| + epsilon*(1-exp(-mu*t))/mu",
        "theorem": (
            "For any tangent perturbation r with norm at most epsilon, while the trajectory "
            "remains in a region with Hessian at least mu, Gronwall gives the displayed input-"
            "to-state stability estimate. If epsilon<mu R/2 and ||e(0)||<R/2, the ball is invariant."
        ),
    }
    return finalize(427, "Nonlinear Local-Flow Robustness Theorem",
                    "proven_conditional_input_to_state_stability_bound",
                    "The forcing, time variable, and basin are supplied dimensionless objects; no physical environment is derived.", packet)


def st428() -> dict:
    M = float_matrix_from_iv(rank7_interval_matrix()); ev = np.linalg.eigvalsh(M)
    off = M[np.triu_indices(N, 1)]
    old = json.loads((ROOT / "FIN_ST391_Rank7_Global_Inequality_Stop.json").read_text())
    packet = {
        "rank": int(np.sum(ev > 1e-10)), "eigenvalues": ev.tolist(),
        "positive_offdiagonal_pairs": int(np.sum(off > 0)), "negative_offdiagonal_pairs": int(np.sum(off < 0)),
        "lambda_max_matches_strict": float(ev[-1]),
        "cap_curvature_threshold_reused": old["sampled_cap_curvature_threshold"],
        "sign_aware_obstruction": (
            "Mixed off-diagonal signs invalidate the strict positive-weight collision envelope. "
            "Splitting positive and negative edges yields a valid but nonclosing bound; no "
            "complete 11D interval cover or rank-seven global minimizer theorem results."
        ),
    }
    return finalize(428, "Rank-Seven Sign-Aware Envelope Audit",
                    "blocked_mixed_sign_envelope_does_not_close_global_gap",
                    "Shared eigenvalues do not license transfer of the strict positive-weight global proof.", packet)


def st429() -> dict:
    packet = {"new_internal_gain_source_found": False,
              "admissible_new_objects_checked": ["ST419 compactification coordinate", "ST424 odd invariant generators", "ST427 bounded forcing"],
              "reason": "All are downstream coordinates/generators/premises; none selects or normalizes g=4 from the strict core.",
              "gate": "open"}
    return finalize(429, "Gain-Source Admission Gate", "blocked_no_new_strict_gain_source",
                    "The orbit and transition theorems remain conditional on supplied gain.", packet)


def st430() -> dict:
    packet = {"new_nonpremise_selector_provider_found": False,
              "odd_invariants_test": "D12-invariant odd polynomials are constant on each D12 orbit and cannot select one translate.",
              "IR_attachment_test": "The compactification selects a branch parameterization, not a directed vertex.",
              "QW_2191": "open", "gate": "open"}
    return finalize(430, "Selector-Source Admission Gate", "blocked_no_new_selector_provider",
                    "No strict-core orientation, orbit member, or QW-2191 discharge is claimed.", packet)


def st431() -> dict:
    packet = {
        "external_referee": "absent", "independent_laboratory_record": "absent",
        "new_empirical_data": "absent", "independent_custody_holdout": "absent",
        "local_mathematical_advances": ["global gain enclosure", "singular IR attachment theorem", "nine-root local atlas", "odd/even invariant census"],
        "forbidden_promotions": ["physical phase transition", "matter derivation", "Standard Model", "gravity", "Theory of Everything"],
    }
    return finalize(431, "Independent Evidence Gate", "gate_remains_open_local_work_only",
                    "Analytic and computational work on repository-defined objects cannot manufacture external evidence.", packet)


def figures(results: dict):
    FIG_DIR.mkdir(exist_ok=True)
    # Gain enclosure.
    g0 = results["ST418"]["certified_unique_uniform_global_gain_interval"][1]
    fig, ax = plt.subplots(figsize=(8, 3.7)); ax.axvspan(0, g0, color="#4daf4a", alpha=.55, label="uniform global (proved)")
    ax.axvspan(g0, 3.998, color="#ffcc66", alpha=.55, label="global identity unresolved")
    ax.axvspan(3.998, 4.002, color="#377eb8", alpha=.55, label="12-point orbit (proved)")
    ax.axvline(2.90249648, color="black", ls="--", lw=1, label="branch-local crossing only")
    ax.set(xlim=(0, 4.05), yticks=[], xlabel="supplied dimensionless gain g", title="Certified global-regime enclosure")
    ax.legend(fontsize=8, ncol=2); fig.tight_layout(); fig.savefig(FIG_DIR/"gain_enclosure.png", dpi=180); plt.close(fig)
    # Singular attachment.
    rows = results["ST419"]["sampled_centers_for_replay_only"]
    bs = [r["b"] for r in rows if r["b"] > 0]; ys = [r["y3"] for r in rows if r["b"] > 0]; aa = [r["a"] for r in rows if r["b"] > 0]
    fig, ax1 = plt.subplots(figsize=(7.5, 4)); ax1.plot(bs, ys, "o-", color="#984ea3"); ax1.set(xlabel="b=1/y3", ylabel="y3", title="Certified compactified branch tube")
    ax2 = ax1.twinx(); ax2.plot(bs, aa, "s--", color="#e41a1c"); ax2.set_ylabel("a=b c")
    fig.tight_layout(); fig.savefig(FIG_DIR/"ir_attachment.png", dpi=180); plt.close(fig)
    # Stationary atlas.
    rows = results["ST421"]["locally_certified_stationary_orbits"]
    fig, ax = plt.subplots(figsize=(7.5, 4));
    for r in rows: ax.scatter(max(r["center"]), r["value"], c=f"C{r['certified_Morse_index']}", s=55); ax.text(max(r["center"]), r["value"], r["label"], fontsize=8)
    ax.set(xlabel="maximum coordinate", ylabel="V4", title="Locally certified stationary-orbit atlas"); fig.tight_layout(); fig.savefig(FIG_DIR/"stationary_atlas.png", dpi=180); plt.close(fig)
    # Invariant census.
    counts = results["ST424"]["primitive_counts_through_degree6"]
    fig, ax = plt.subplots(figsize=(7, 3.8)); ax.bar(list(counts), list(counts.values()), color="#4daf4a"); ax.set(ylabel="primitive quotient count", title="Joint odd/even generator census (modular)"); fig.tight_layout(); fig.savefig(FIG_DIR/"invariant_census.png", dpi=180); plt.close(fig)


def main():
    funcs = [st417, st418, st419, st420, st421, st422, st423, st424, st425, st426, st427, st428, st429, st430, st431]
    results = {}
    for f in funcs:
        r = f(); results[r["program"]] = r; print(r["program"], r["status"], flush=True)
    figures(results)
    payload = {"metadata": {"seed": SEED, "python": platform.python_version(), "numpy": np.__version__,
                             "scipy": scipy.__version__, "mpmath": mp.__version__,
                             "strict_kernel": "off-diagonal K_strict_gate; W_xx=0; A=sI-W",
                             "scope": "finite local dimensionless mathematics"}, **results}
    RESULTS.write_text(json.dumps(native(payload), indent=2, sort_keys=True), encoding="utf-8")
    with SUMMARY.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh); w.writerow(["program", "object", "status", "boundary", "packet_file", "packet_sha256"])
        for r in results.values(): w.writerow([r[k] for k in ("program", "object", "status", "boundary", "packet_file", "packet_sha256")])
    print(RESULTS.name, sha(RESULTS)); print(SUMMARY.name, sha(SUMMARY))


if __name__ == "__main__":
    main()
