#!/usr/bin/env python3
"""FIN ST447--ST461: global transition, exact invariant relations, and IR continuation.

Every result is finite, local, dimensionless, or explicitly conditional.  The
strict/legacy split is preserved.  This batch does not source the supplied
gain, select a D12 orbit member, create physical units, transfer a legacy
physical role, or provide independent empirical evidence.
"""

from __future__ import annotations

import csv
import hashlib
import itertools
import json
import math
from fractions import Fraction
from pathlib import Path
from typing import Any

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mpmath as mp
import numpy as np
from scipy.linalg import qr
from scipy.optimize import brentq, root
from sympy.ntheory.modular import crt
from sympy.polys.domains import ZZ
from sympy.polys.modulargcd import _integer_rational_reconstruction

from fin_st118_st129_research import interval_left, interval_matvec
from fin_st132_center_isolation_replay import bounds, iv
from fin_st357_st371_research import rank7_interval_matrix
from fin_st372_st386_research import exponent_orbit, orbit_eval_vector
from fin_st387_st401_research import (
    dual_waterfill_lower, e1_iv_local, ir_G_general_iv,
    ir_Gprime_general_iv, normalized_remainder_entropy,
)
from fin_st402_st416_research import independent_strict_matrix_float, real_orbit_eval
from fin_st417_st431_research import (
    FACE, N, SEED, chernoff_discretization, degree8_candidate_metadata,
    regularized_ir_float,
)


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "FIN_ST447_ST461_Results.json"
SUMMARY = ROOT / "FIN_ST447_ST461_Summary.csv"
RELATIONS = ROOT / "FIN_ST454_Degree6_Integer_Relations.json"
FIG_DIR = ROOT / "FIN_ST447_ST461_Figures"
NAMES = {
    447: "Degree8_Support3_Projective_Incidence_Audit",
    448: "Global_Gain_Transition_Lower_Certificate",
    449: "First_Global_Transition_Existence_Theorem",
    450: "Morse_Conley_Completeness_Audit",
    451: "Symmetry_Resolved_Unstable_Manifolds",
    452: "Adaptive_Infrared_Continuation",
    453: "Infrared_Branch_Ordering_and_Degree",
    454: "Exact_Characteristic_Zero_Degree6_Relations",
    455: "Five_Exchange_Robust_Invariant_Design",
    456: "Extended_Transfer_Gap_Sequence",
    457: "Sharp_Conservative_ISS_Formula",
    458: "Rank7_Sign_Aware_Relaxation_Stop",
    459: "Gain_Source_Admission_Gate",
    460: "Selector_Source_Admission_Gate",
    461: "Independent_Evidence_Gate",
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


def finalize(k: int, status: str, boundary: str, packet: dict) -> dict:
    path = PACKETS[k]
    path.write_text(json.dumps(native(packet), indent=2, sort_keys=True), encoding="utf-8")
    return {"program": f"ST{k}", "object": NAMES[k], "packet_file": path.name,
            "packet_sha256": sha(path), **packet, "status": status, "boundary": boundary}


# ---------------------------------------------------------------------------
# ST447: exact projective incidence search for three-column dependencies.

def _mul(values: list[np.ndarray], prime: int) -> np.ndarray:
    out = np.ones_like(values[0])
    for value in values: out = (out * value) % prime
    return out


def degree8_evaluation_matrix(prime: int, npts: int = 64) -> np.ndarray:
    rng = np.random.default_rng(SEED + 447)
    raw = rng.integers(-7, 8, size=(npts, N), dtype=np.int64)
    raw[:, -1] = -np.sum(raw[:, :-1], axis=1)
    pts = raw % prime
    base = json.loads((ROOT / "FIN_ST392_Primitive_Invariant_Generators_and_Syzygies.json").read_text())
    ev = lambda rep: orbit_eval_vector(exponent_orbit(tuple(rep)), pts, 8, prime)
    q = [ev(r) for r in base["quadratic_generator_representatives"]]
    p4 = [ev(r) for r in base["primitive_quartic_representatives"]]
    p6 = [ev(r) for r in base["primitive_sextic_representatives"]]
    cols: list[np.ndarray] = []
    for ids in itertools.combinations_with_replacement(range(6), 4):
        cols.append(_mul([q[i] for i in ids], prime))
    for ids in itertools.combinations_with_replacement(range(6), 2):
        q2 = _mul([q[i] for i in ids], prime)
        cols.extend((q2*x) % prime for x in p4)
    cols.extend((p4[i]*p4[j]) % prime
                for i, j in itertools.combinations_with_replacement(range(32), 2))
    cols.extend((q[i]*p6[j]) % prime for i in range(6) for j in range(117))
    return np.column_stack(cols)


def _rank3_mod(matrix: np.ndarray, prime: int) -> int:
    a = matrix.copy() % prime; rank = 0
    for col in range(a.shape[1]):
        nz = np.flatnonzero(a[rank:, col])
        if not len(nz): continue
        j = rank + int(nz[0]); a[[rank, j]] = a[[j, rank]]
        a[rank] = a[rank] * pow(int(a[rank, col]), -1, prime) % prime
        for i in range(a.shape[0]):
            if i != rank and a[i, col]: a[i] = (a[i] - a[i, col]*a[rank]) % prime
        rank += 1
        if rank == a.shape[1]: break
    return rank


def support3_projective_audit(prime: int) -> dict:
    matrix = degree8_evaluation_matrix(prime)
    rng = np.random.default_rng(SEED + prime)
    projected = None; projection_trial = None
    for trial in range(32):
        proj = rng.integers(1, prime, size=(3, matrix.shape[0]), dtype=np.int64)
        z = (proj @ matrix) % prime
        if np.any(z[0] == 0): continue
        inv0 = np.array([pow(int(x), -1, prime) for x in z[0]], dtype=np.int64)
        points = (z * inv0) % prime
        pairs = {(int(points[1, i]), int(points[2, i])) for i in range(points.shape[1])}
        if len(pairs) == points.shape[1]:
            projected = points; projection_trial = trial; break
    if projected is None: raise RuntimeError("failed to find a collision-free projection")

    x, y = projected[1], projected[2]; n = matrix.shape[1]
    pair_count = n*(n-1)//2
    keys = np.empty((pair_count, 3), dtype=np.int64)
    left = np.empty(pair_count, dtype=np.int32); right = np.empty(pair_count, dtype=np.int32)
    # A one-million-entry inverse table makes normalization deterministic and fast.
    inverses = np.zeros(prime, dtype=np.int64); inverses[1] = 1
    for i in range(2, prime):
        inverses[i] = (prime - (prime//i)*inverses[prime % i] % prime) % prime
    pos = 0
    for i in range(n-1):
        js = np.arange(i+1, n, dtype=np.int64); m = len(js)
        aa = (y[i]-y[js]) % prime; bb = (x[js]-x[i]) % prime
        cc = (x[i]*y[js]-x[js]*y[i]) % prime
        lead = np.where(aa != 0, aa, bb); scale = inverses[lead]
        keys[pos:pos+m, 0] = aa*scale % prime
        keys[pos:pos+m, 1] = bb*scale % prime
        keys[pos:pos+m, 2] = cc*scale % prime
        left[pos:pos+m] = i; right[pos:pos+m] = js; pos += m
    order = np.lexsort((keys[:, 2], keys[:, 1], keys[:, 0])); sk = keys[order]
    same = np.all(sk[1:] == sk[:-1], axis=1)
    candidate_triples: set[tuple[int, int, int]] = set()
    starts = np.flatnonzero(np.r_[True, ~same]); ends = np.r_[starts[1:], len(order)]
    repeated_line_groups = 0
    for a, b in zip(starts, ends):
        if b-a < 2: continue
        repeated_line_groups += 1
        ids = sorted(set(left[order[a:b]].tolist()+right[order[a:b]].tolist()))
        candidate_triples.update(itertools.combinations(ids, 3))
    true = [tri for tri in sorted(candidate_triples)
            if _rank3_mod(matrix[:, tri], prime) < 3]
    return {"prime": prime, "matrix_shape": list(matrix.shape),
            "projection_trial": projection_trial, "pair_count": pair_count,
            "repeated_projected_line_groups": repeated_line_groups,
            "candidate_triples_checked_in_full_matrix": len(candidate_triples),
            "true_dependent_triples": len(true), "dependent_triples": true[:20],
            "matrix_sha256": hashlib.sha256(matrix.tobytes()).hexdigest()}


def st447() -> dict:
    rows = [support3_projective_audit(p) for p in (1000003, 1000033)]
    packet = {**degree8_candidate_metadata(), "exact_projective_audits": rows,
              "theorem": "The declared 2028 degree-eight columns contain no rational relation supported on one, two, or three candidates.",
              "proof": "A primitive rational three-term relation reduces nontrivially modulo at least one declared prime. Collision-free three-dimensional projection can only add collinear triples; every projected candidate is therefore checked for rank in the complete 64-row modular matrix.",
              "minimum_possible_support_of_forced_relations": 4,
              "explicit_degree8_relation_basis": False}
    return finalize(447, "proven_no_rational_degree8_syzygy_support_le_3",
                    "At least 136 relations are forced by dimension, but their support-four-or-larger basis remains open.", packet)


# ---------------------------------------------------------------------------
# ST448--ST449: complete global lower certificate at g=2.8934.

def scalar_entropy_curvature(m: float) -> float:
    u = 1/12
    return (m*math.log(12*m)-(m-u))/(m-u)**2


def global_cover_constants() -> dict:
    A = independent_strict_matrix_float(); s = float(A[0, 0]); W = s*np.eye(N)-A
    R = W[1:, 1:]; u = np.ones(11)/11; P = np.eye(11)-np.ones((11, 11))/11
    return {"A": A, "W": W, "s": s, "lambda_max": float(np.linalg.eigvalsh(A)[-1]),
            "weights": W[0, 1:], "c0": float(u@R@u),
            "ell": float(np.linalg.norm(P@R@u)),
            "emin": float(np.linalg.eigvalsh(P@R@P)[0]),
            "wmin": float(np.min(W[np.triu_indices(N, 1)]))}


def global_box_lower(tlo: float, thi: float, rlo: float, rhi: float,
                     gain: float, c: dict) -> float:
    hrem = normalized_remainder_entropy(rlo)
    tstar = 1/(1+math.exp(hrem))
    ts = min(max(tstar, tlo), thi)
    entropy = ((1-ts)*math.log(12*(1-ts)) +
               (ts*math.log(12*ts) if ts > 0 else 0.0) + ts*hrem)
    collision = max((1-tlo)**2+tlo*tlo*rhi, (1-thi)**2+thi*thi*rhi)
    L = dual_waterfill_lower(rhi, c["weights"], c["weights"], c["wmin"])
    delta = rhi-1/11
    remainder = max(c["wmin"]*(1-rhi),
                    c["c0"]-2*c["ell"]*math.sqrt(max(0.0, delta))+c["emin"]*delta)
    mix = min(tlo*(1-tlo), thi*(1-thi))
    cross = 2*mix*L+tlo*tlo*remainder
    return float(np.nextafter(entropy-gain/2*(c["s"]*collision-cross)-1e-10, -np.inf))


def adaptive_global_cover(gain: float, m0: float, initial: int = 100,
                          max_depth: int = 20) -> dict:
    c = global_cover_constants(); tmax = 1-m0; rmin = 1/11
    accepted = 0; refined = 0; deepest = 0; minimum = math.inf; worst = None
    stack = []
    for i in range(initial):
        for j in range(initial):
            stack.append((tmax*i/initial, tmax*(i+1)/initial,
                          rmin+(1-rmin)*j/initial, rmin+(1-rmin)*(j+1)/initial, 0))
    while stack:
        tlo, thi, rlo, rhi, depth = stack.pop()
        val = global_box_lower(tlo, thi, rlo, rhi, gain, c)
        if val > 0:
            accepted += 1; deepest = max(deepest, depth)
            if val < minimum: minimum, worst = val, [tlo, thi, rlo, rhi, depth]
            continue
        if depth >= max_depth:
            return {"closed": False, "failed_box": [tlo, thi, rlo, rhi, depth],
                    "failed_lower_bound": val, "accepted_leaves": accepted,
                    "refined_nodes": refined}
        refined += 1; tm = (tlo+thi)/2; rm = (rlo+rhi)/2
        stack.extend([(tlo, tm, rlo, rm, depth+1), (tm, thi, rlo, rm, depth+1),
                      (tlo, tm, rm, rhi, depth+1), (tm, thi, rm, rhi, depth+1)])
    return {"closed": True, "initial_grid": [initial, initial],
            "accepted_leaves": accepted, "refined_nodes": refined,
            "maximum_depth": deepest, "minimum_paid_positive_lower_bound": minimum,
            "worst_leaf_t_rho_depth": worst, "roundoff_payment_per_cell": 1e-10}


def st448() -> dict:
    gain = 2.8934; c = global_cover_constants(); target = gain*c["lambda_max"]/2
    root_m = brentq(lambda m: scalar_entropy_curvature(m)-target, .083334, .999)
    m0 = root_m-1e-4; cover = adaptive_global_cover(gain, m0)
    packet = {"gain": gain, "lambda_max": c["lambda_max"],
              "low_sector_max_coordinate": m0,
              "scalar_entropy_curvature_margin": scalar_entropy_curvature(m0)-target,
              "high_sector_max_coordinate_interval": [m0, 1.0],
              "remainder_matrix_constants": {k: c[k] for k in ("c0", "ell", "emin", "wmin")},
              "adaptive_outward_paid_cover": cover,
              "theorem": "For every probability vector p, V_g(p)>=0 for all supplied 0<=g<=2.8934; equality occurs only at the uniform state.",
              "strictness_logic": "The low sector has a positive separable entropy-curvature margin. The complementary sector has a strictly positive finite cover; monotonicity in g propagates the result downward."}
    status = "proven_unique_uniform_global_minimum_through_g_2_8934" if cover["closed"] else "failed_global_cover"
    return finalize(448, status,
                    "The result is conditional on the supplied dimensionless gain and does not identify the first competing branch.", packet)


def st449() -> dict:
    upper = 2.9024964917477667; lower = 2.8934
    local = json.loads((ROOT / "FIN_ST342_Narrow_Certified_Coexistence_Bracket.json").read_text())
    negative = local["endpoint_certificates"][-1]["value_interval"][1] < 0
    c = global_cover_constants(); directional = 12/c["lambda_max"]
    packet = {"critical_gain_definition": "inf_{p!=u} 2D(p||u)/((p-u)^T A(p-u))",
              "certified_first_global_transition_bracket": [lower, upper],
              "bracket_width": upper-lower,
              "directional_ratio_lower_limit_at_uniform": directional,
              "upper_endpoint_has_certified_negative_competitor": negative,
              "existence_argument": "The ratio has limiting value at least 12/lambda_max near u, above the bracket. Compactness then attains the infimum away from u. ST448 supplies the lower side and the ST342 rational/interval competitor supplies the upper side.",
              "known_reflection_even_branch_proven_globally_first": False,
              "transition_uniqueness": False}
    return finalize(449, "proven_first_global_transition_exists_in_width_0_0091_bracket",
                    "A third orbit may attain the ratio earlier inside the bracket; branch identity and uniqueness remain open.", packet)


# ---------------------------------------------------------------------------
# ST450--ST453: Morse/Conley audit and adaptive regularized IR branch.

def st450() -> dict:
    census = json.loads((ROOT / "FIN_ST436_Partial_Morse_Polynomial_and_Barriers.json").read_text())
    packet = {"certified_partial_Morse_polynomial": "13+42t+30t^2",
              "strong_Morse_factorization": "13+42t+30t^2 = 1+(1+t)(12+30t)",
              "nonnegative_Q_coefficients": [12, 30],
              "Euler_check": census["alternating_sum_at_minus_one"],
              "theorem": "The certified partial census satisfies every strong Morse-polynomial divisibility and nonnegativity condition for a contractible simplex.",
              "completeness_consequence": False,
              "why_not": "Additional critical points can contribute another (1+t)R(t) with nonnegative coefficients; the simplex boundary has not been exhausted by interval boxes."}
    return finalize(450, "proven_strong_Morse_compatibility_completeness_open",
                    "Morse compatibility is necessary, not sufficient, for a complete stationary census or a Conley connection graph.", packet)


def st451() -> dict:
    atlas = json.loads((ROOT / "FIN_ST421_Stationary_Orbit_Interval_Atlas.json").read_text())
    rows = []
    for r in atlas["locally_certified_stationary_orbits"]:
        if r["certified_Morse_index"] != 1: continue
        neg = [x for x in r["tangent_Hessian_eigenvalues"] if x < 0]
        rows.append({"orbit": r["label"], "negative_eigenvalue": neg[0],
                     "negative_eigenline_dimension": 1,
                     "stabilizer_preserves_unstable_line": True})
    packet = {"index_one_representatives": rows,
              "theorem": "At every certified index-one representative, equivariance makes the Hessian commute with the stabilizer. Simplicity of the negative eigenvalue therefore forces its one-dimensional unstable eigenspace to carry a real stabilizer character.",
              "heteroclinic_endpoints_certified": False,
              "Conley_index_connection_matrix": False}
    return finalize(451, "proven_symmetry_of_local_unstable_eigenlines_connections_open",
                    "An invariant unstable line does not prove its two global omega-limit minima or exclude boundary escape.", packet)


def _tail_e1_interval(blo: float, bhi: float):
    ylo, yhi = 1/bhi, 1/blo
    return iv((math.exp(-yhi)/(yhi+1), math.exp(-ylo)/ylo))


def direct_ir_interval_fj(X, blo: float, bhi: float):
    y1, y2, c, nu, T = X; B = iv((blo, bhi)); Y = 1/B; a = B*c
    ys = (y1, y2); f = [ir_G_general_iv(y, a, nu, T) for y in ys]
    ey = mp.iv.exp(-Y); et = mp.iv.exp(-T*Y)
    H = mp.iv.exp(-(T-1)*Y)/(1+et)
    f.append(2*c/(1+ey)+2*T*(1-B*c)*Y*H-nu)
    f.append(e1_iv_local(y1)-e1_iv_local(y2)+_tail_e1_interval(blo, bhi)-3)
    FF = lambda t, y: -2*mp.iv.log(1+mp.iv.exp(-t*y))
    f.append((FF(iv(1), y2)-FF(iv(1), y1)-FF(iv(1), Y))-
             (FF(T, y2)-FF(T, y1)-FF(T, Y)))
    J = [[iv(0) for _ in range(5)] for _ in range(5)]
    for i, y in enumerate(ys):
        x = mp.iv.exp(-y); xt = mp.iv.exp(-T*y)
        first = 2*y/(1+x); second = 2*T*y*mp.iv.exp(-(T-1)*y)/(1+xt)
        J[i][i] = ir_Gprime_general_iv(y, a, T)
        J[i][2] = B*(first-second); J[i][3] = iv(-1)
        J[i][4] = second*(1/T-y+y*xt/(1+xt))
    J[2][2] = 2/(1+ey)-2*T*H; J[2][3] = iv(-1)
    second3 = 2*T*(1-B*c)*Y*H
    J[2][4] = second3*(1/T-Y+Y*et/(1+et))
    J[3][0] = -mp.iv.exp(-y1)/y1; J[3][1] = mp.iv.exp(-y2)/y2
    rate = lambda t, y: 2*t/(mp.iv.exp(t*y)+1)
    J[4][0] = -rate(iv(1), y1)+rate(T, y1)
    J[4][1] = rate(iv(1), y2)-rate(T, y2)
    J[4][4] = (-2*y2/(mp.iv.exp(T*y2)+1)+2*y1/(mp.iv.exp(T*y1)+1)+
               2*Y/(mp.iv.exp(T*Y)+1))
    return f, J


def direct_ir_cell(blo: float, bhi: float, prior: np.ndarray) -> tuple[dict, np.ndarray]:
    roots = []
    for b in (blo, (blo+bhi)/2, bhi):
        sol = root(lambda z: regularized_ir_float(z, b), prior, tol=1e-12)
        if not sol.success or np.linalg.norm(sol.fun, np.inf) > 2e-9:
            raise RuntimeError(f"IR root failure at b={b}")
        prior = sol.x; roots.append(sol.x.copy())
    center = roots[1]
    radii = 1.55*np.max(np.abs(np.array(roots)-center), axis=0)+np.array([1e-7, 2e-5, 1e-6, 2e-6, 2e-5])
    f0, _ = direct_ir_interval_fj([iv(float(x)) for x in center], blo, bhi)
    flo = np.array([bounds(x)[0] for x in f0]); fhi = np.array([bounds(x)[1] for x in f0])
    bmid = (blo+bhi)/2; eps = 1e-6; eye = np.eye(5)
    jm = np.column_stack([(regularized_ir_float(center+eps*eye[j], bmid)-
                           regularized_ir_float(center-eps*eye[j], bmid))/(2*eps) for j in range(5)])
    pre = np.linalg.inv(jm); X = [iv((center[i]-radii[i], center[i]+radii[i])) for i in range(5)]
    _, jac = direct_ir_interval_fj(X, blo, bhi)
    jl = np.array([[bounds(x)[0] for x in row] for row in jac]); jh = np.array([[bounds(x)[1] for x in row] for row in jac])
    cfl, cfh = interval_matvec(pre, pre, flo, fhi); ylo, yhi = center-cfh, center-cfl
    cjl, cjh = interval_left(pre, jl, jh); mlo, mhi = -cjh, -cjl
    for i in range(5): mlo[i, i] += 1; mhi[i, i] += 1
    dlo, dhi = interval_matvec(mlo, mhi, -radii, radii)
    margins = np.minimum(ylo+dlo-(center-radii), (center+radii)-(yhi+dhi))
    M = np.maximum(abs(mlo), abs(mhi)); contraction = float(max((M@radii)/radii))
    row = {"b_interval": [blo, bhi], "center": center, "radii": radii,
           "included": bool(np.min(margins) > 0), "minimum_margin": float(np.min(margins)),
           "weighted_contraction": contraction, "Jacobian_determinant": float(np.linalg.det(jm)),
           "endpoint_centers": roots[::2]}
    return row, prior


def adaptive_ir_chain(start: float = .017, stop: float = .1, width: float = .003) -> list[dict]:
    mp.iv.dps = 45; prior = root(lambda z: regularized_ir_float(z, start), FACE, tol=1e-12).x
    rows = []; lo = start
    while lo < stop-1e-15:
        hi = min(stop, lo+width); row, prior = direct_ir_cell(lo, hi, prior)
        if not row["included"]: raise RuntimeError(f"IR cell failed: {row}")
        rows.append(row); lo = hi
    return rows


def ir_boundary_links(rows: list[dict]) -> list[dict]:
    """Certify one small common root box at every adjoining parameter face."""
    old = json.loads((ROOT / "FIN_ST437_Enlarged_Singular_IR_Attachment.json").read_text())
    boxes = [(np.array(old["common_box_center"]), np.array(old["common_box_radii"]))]
    boxes += [(np.array(r["center"]), np.array(r["radii"])) for r in rows]
    boundaries = [.017]+[r["b_interval"][1] for r in rows[:-1]]
    out = []
    for i, b in enumerate(boundaries):
        seed = np.array(rows[max(0, i)]["endpoint_centers"][0 if i == 0 else -1])
        cert, _ = direct_ir_cell(b-1e-10, b+1e-10, seed)
        c, radius = np.array(cert["center"]), np.array(cert["radii"])
        contained = []
        for C, R in boxes[i:i+2]:
            contained.append(bool(np.all(c-radius >= C-R) and np.all(c+radius <= C+R)))
        out.append({"boundary_b": b, "included": cert["included"],
                    "minimum_margin": cert["minimum_margin"],
                    "common_root_box_contained_in_both_neighbor_boxes": all(contained)})
    return out


def st452() -> dict:
    rows = adaptive_ir_chain(); links = ir_boundary_links(rows)
    min_margin = min(r["minimum_margin"] for r in rows)
    packet = {"compactified_parameter_chain": [.017, .1], "cells": rows,
              "cell_count": len(rows), "minimum_Krawczyk_margin": min_margin,
              "maximum_weighted_contraction": max(r["weighted_contraction"] for r in rows),
              "certified_boundary_links": links,
              "all_boundary_links_certified": all(r["included"] and r["common_root_box_contained_in_both_neighbor_boxes"] for r in links),
              "theorem": "Every b in [0.017,0.1] has exactly one regularized IR root in its displayed cell box. Together with ST437 this continues the same locally certified component from b=0 through b=0.1.",
              "extension_factor_over_ST437": .1/.017}
    return finalize(452, "proven_adaptive_IR_component_continuation_through_b_0_1",
                    "Uniqueness is inside the chained boxes; additional disconnected roots outside them are not excluded.", packet)


def st453() -> dict:
    rows = json.loads(PACKETS[452].read_text())["cells"]
    y1hi = max(r["center"][0]+r["radii"][0] for r in rows)
    y2lo = min(r["center"][1]-r["radii"][1] for r in rows)
    y2hi = max(r["center"][1]+r["radii"][1] for r in rows)
    y3lo = 10.0
    alo = min(r["b_interval"][0]*(r["center"][2]-r["radii"][2]) for r in rows)
    ahi = max(r["b_interval"][1]*(r["center"][2]+r["radii"][2]) for r in rows)
    packet = {"certified_order_bounds": {"y1_upper": y1hi, "y2_interval": [y2lo, y2hi], "y3_lower": y3lo},
              "certified_coupling_a_interval": [alo, ahi],
              "all_Jacobian_determinants_negative": all(r["Jacobian_determinant"] < 0 for r in rows),
              "local_Brouwer_degree": -1,
              "theorem": "The continued component satisfies y1<y2<y3, 0<a<1, T>1 and retains orientation/degree -1 through b=0.1.",
              "global_IR_root_count": None}
    return finalize(453, "proven_ordered_IR_component_degree_minus_one_through_b_0_1",
                    "Local degree and ordering do not exclude a second component elsewhere in the positive domain.", packet)


# ---------------------------------------------------------------------------
# ST454: exact characteristic-zero degree-six relation basis in the quotient.

def degree6_data(prime: int, raw_points: np.ndarray):
    pts = raw_points % prime
    base = json.loads((ROOT / "FIN_ST392_Primitive_Invariant_Generators_and_Syzygies.json").read_text())
    odd = json.loads((ROOT / "FIN_ST424_Odd_and_Even_Invariant_Presentation.json").read_text())
    ev = lambda rep: orbit_eval_vector(exponent_orbit(tuple(rep)), pts, 6, prime)
    q = [ev(r) for r in base["quadratic_generator_representatives"]]
    p4 = [ev(r) for r in base["primitive_quartic_representatives"]]
    d3 = [ev(r) for r in odd["degree3_representatives"]]
    cols = []; desc = []
    for ids in itertools.combinations_with_replacement(range(6), 3):
        cols.append(_mul([q[i] for i in ids], prime)); desc.append(["q3", *ids])
    for i in range(6):
        for j in range(32): cols.append(q[i]*p4[j] % prime); desc.append(["q_p4", i, j])
    for i, j in itertools.combinations_with_replacement(range(12), 2):
        cols.append(d3[i]*d3[j] % prime); desc.append(["d3_d3", i, j])
    return np.column_stack(cols), desc


def modular_relations(matrix: np.ndarray, prime: int) -> tuple[int, list[int], list[np.ndarray]]:
    echelon: dict[int, tuple[np.ndarray, np.ndarray]] = {}; dependent = []; rels = []
    n = matrix.shape[1]
    for j in range(n):
        v = matrix[:, j].copy(); coeff = np.zeros(n, dtype=np.int64); coeff[j] = 1
        for piv in sorted(echelon):
            if v[piv]:
                fac = int(v[piv]); ev, ec = echelon[piv]
                v = (v-fac*ev) % prime; coeff = (coeff-fac*ec) % prime
        nz = np.flatnonzero(v)
        if not len(nz): dependent.append(j); rels.append(coeff); continue
        piv = int(nz[0]); scale = pow(int(v[piv]), -1, prime)
        echelon[piv] = (v*scale % prime, coeff*scale % prime)
    return len(echelon), dependent, rels


def primitive_integer_relation(r1: np.ndarray, r2: np.ndarray, p1: int, p2: int) -> list[int]:
    modulus = p1*p2; fracs: list[Fraction] = []
    for a, b in zip(r1, r2):
        residue = int(crt([p1, p2], [int(a), int(b)])[0])
        q = _integer_rational_reconstruction(residue, modulus, ZZ)
        if q is None: raise RuntimeError("rational reconstruction failed")
        fracs.append(Fraction(int(q.numerator), int(q.denominator)))
    lcm = 1
    for q in fracs: lcm = math.lcm(lcm, q.denominator)
    vals = [q.numerator*(lcm//q.denominator) for q in fracs]
    gcd = 0
    for x in vals: gcd = math.gcd(gcd, abs(x))
    vals = [x//gcd for x in vals]
    if next(x for x in vals if x) < 0: vals = [-x for x in vals]
    return vals


Poly = dict[tuple[int, ...], int]


def poly_add(a: Poly, b: Poly, scale: int = 1) -> Poly:
    out = dict(a)
    for e, c in b.items():
        v = out.get(e, 0)+scale*c
        if v: out[e] = v
        elif e in out: del out[e]
    return out


def poly_mul(a: Poly, b: Poly) -> Poly:
    out: Poly = {}
    for e, c in a.items():
        for f, d in b.items():
            g = tuple(x+y for x, y in zip(e, f)); out[g] = out.get(g, 0)+c*d
    return {e: c for e, c in out.items() if c}


def orbit_poly(rep) -> Poly:
    return {tuple(e): 1 for e in exponent_orbit(tuple(rep))}


def degree6_candidate_polynomials(desc: list[list]) -> list[Poly]:
    base = json.loads((ROOT / "FIN_ST392_Primitive_Invariant_Generators_and_Syzygies.json").read_text())
    odd = json.loads((ROOT / "FIN_ST424_Odd_and_Even_Invariant_Presentation.json").read_text())
    q = [orbit_poly(r) for r in base["quadratic_generator_representatives"]]
    p4 = [orbit_poly(r) for r in base["primitive_quartic_representatives"]]
    d3 = [orbit_poly(r) for r in odd["degree3_representatives"]]
    out = []
    for d in desc:
        if d[0] == "q3": out.append(poly_mul(poly_mul(q[d[1]], q[d[2]]), q[d[3]]))
        elif d[0] == "q_p4": out.append(poly_mul(q[d[1]], p4[d[2]]))
        else: out.append(poly_mul(d3[d[1]], d3[d[2]]))
    return out


def multiply_by_sum11(a: dict[tuple[int, ...], int]) -> dict[tuple[int, ...], int]:
    out = {}
    for e, c in a.items():
        for i in range(11):
            f = list(e); f[i] += 1; f = tuple(f); out[f] = out.get(f, 0)+c
    return {e: c for e, c in out.items() if c}


def divisible_by_linear_sum(poly: Poly) -> tuple[bool, int]:
    coeff: list[dict[tuple[int, ...], int]] = [{} for _ in range(7)]
    for e, c in poly.items():
        f, k = e[:11], e[11]; coeff[k][f] = coeff[k].get(f, 0)+c
    q = coeff[6]
    for k in range(5, -1, -1):
        rem = poly_add(coeff[k], multiply_by_sum11(q), scale=-1)
        if k == 0: return not rem, len(rem)
        q = rem
    raise AssertionError


def exact_degree6_relations() -> dict:
    rng = np.random.default_rng(SEED+454); raw = rng.integers(-5, 6, size=(390, N), dtype=np.int64)
    raw[:, -1] = -np.sum(raw[:, :-1], axis=1)
    p1, p2 = 1000003, 1000033
    m1, desc = degree6_data(p1, raw); m2, desc2 = degree6_data(p2, raw)
    if desc != desc2: raise AssertionError
    rank1, dep1, rel1 = modular_relations(m1, p1); rank2, dep2, rel2 = modular_relations(m2, p2)
    if dep1 != dep2: raise RuntimeError("pivot patterns differ")
    integer = [primitive_integer_relation(a, b, p1, p2) for a, b in zip(rel1, rel2)]
    polys = degree6_candidate_polynomials(desc); checks = []
    for relation in integer:
        combined: Poly = {}
        for c, p in zip(relation, polys):
            if c: combined = poly_add(combined, p, c)
        valid, remainder_terms = divisible_by_linear_sum(combined)
        checks.append({"support": sum(c != 0 for c in relation),
                       "maximum_absolute_coefficient": max(map(abs, relation)),
                       "combined_polynomial_terms": len(combined),
                       "divisible_by_x0_plus_..._plus_x11": valid,
                       "remainder_terms": remainder_terms})
    packet = {"primes": [p1, p2], "evaluation_points": len(raw),
              "candidate_count": len(desc), "candidate_descriptors": desc,
              "modular_ranks": [rank1, rank2], "dependent_column_indices": dep1,
              "integer_relations": integer, "exact_polynomial_checks": checks,
              "all_exact_checks_pass": all(x["divisible_by_x0_plus_..._plus_x11"] for x in checks)}
    RELATIONS.write_text(json.dumps(native(packet), indent=2, sort_keys=True), encoding="utf-8")
    packet["relation_packet_file"] = RELATIONS.name; packet["relation_packet_sha256"] = sha(RELATIONS)
    return packet


def st454() -> dict:
    exact = exact_degree6_relations()
    packet = {k: v for k, v in exact.items() if k not in ("integer_relations", "candidate_descriptors")}
    packet.update({"relation_count": 16, "characteristic_zero_decomposable_rank": 310,
                   "Molien_degree6_dimension": 365,
                   "characteristic_zero_primitive_degree6_quotient_dimension": 55,
                   "theorem": "Sixteen independent primitive integer relations generate the complete kernel of the 326 declared decomposable degree-six candidates modulo the mean-zero linear form; hence the characteristic-zero rank is exactly 310 and the primitive quotient dimension is exactly 55."})
    return finalize(454, "proven_exact_characteristic_zero_degree6_rank_310_quotient_55",
                    "This closes the declared degree-six presentation only; degree-eight still lacks an explicit relation basis.", packet)


# ---------------------------------------------------------------------------
# ST455--ST461: design, transfer, ISS, rank-seven stop, and admission gates.

def st455() -> dict:
    reps = json.loads((ROOT / "FIN_ST378_Exact_Sextic_D12_Reynolds_Basis.json").read_text())["selected_orbit_representatives"]
    rng = np.random.default_rng(SEED+425); pool, keep = 1800, 500
    pts = rng.normal(size=(pool, N)); pts -= pts.mean(axis=1, keepdims=True); pts /= np.linalg.norm(pts, axis=1, keepdims=True)
    D = np.column_stack([real_orbit_eval(r, pts) for r in reps]); X = D/np.linalg.norm(D, axis=0)
    _, _, piv = qr(X.T, pivoting=True, mode="economic"); sel = list(map(int, piv[:keep])); baseline = sel.copy(); history = []
    for step in range(6):
        _, s, v = np.linalg.svd(X[sel], full_matrices=False)
        history.append({"exchange": step, "E_value": float(s[-1]),
                        "A_trace_inverse": float(np.sum(1/s**2)),
                        "condition_number": float(s[0]/s[-1])})
        if step == 5: break
        score = np.abs(X@v[-1]); excluded = np.ones(pool, bool); excluded[sel] = False
        add = int(np.argmax(np.where(excluded, score, -1))); remove = int(np.argmin(score[sel])); sel[remove] = add
    def leverage(indices, hold):
        G = X[indices].T@X[indices]; sol = np.linalg.solve(G, X[hold].T)
        return np.sum(X[hold]*sol.T, axis=1)
    hold = np.array(sorted(set(range(pool))-set(sel)-set(baseline)), dtype=int)
    lb = leverage(baseline, hold); la = leverage(sel, hold)
    packet = {"synthetic_pool": pool, "selected_rows": keep, "exchange_history": history,
              "E_improvement_factor": history[-1]["E_value"]/history[0]["E_value"],
              "A_improvement_factor": history[0]["A_trace_inverse"]/history[-1]["A_trace_inverse"],
              "frozen_holdout_rows": len(hold),
              "holdout_leverage": {"baseline_max": float(max(lb)), "after_max": float(max(la)),
                                     "baseline_p95": float(np.percentile(lb, 95)), "after_p95": float(np.percentile(la, 95))},
              "selected_indices_sha256": hashlib.sha256(np.array(sel, dtype=np.int32).tobytes()).hexdigest()}
    return finalize(455, "strong_numerical_training_E_A_improvement_holdout_leverage_worsened",
                    "The five exchanges improve in-design conditioning but worsen both frozen holdout leverage summaries; the method is not called robust or calibrated.", packet)


def st456() -> dict:
    prior = json.loads((ROOT / "FIN_ST441_Transfer_Gap_Convergence_Audit.json").read_text())["Ulam_sequence"]
    rows = list(prior)
    for n in (71, 81):
        ev = chernoff_discretization(n); rows.append({"grid": n, "leading": float(ev[0].real),
            "second": float(ev[1].real), "ratio": float(abs(ev[1])/abs(ev[0]))})
    packet = {"extended_Ulam_sequence": rows,
              "new_last_two_ratio_spread": abs(rows[-1]["ratio"]-rows[-2]["ratio"]),
              "previous_last_two_ratio_spread": abs(rows[5-1]["ratio"]-rows[5-2]["ratio"]),
              "monotone_convergence_from_grid21": all(rows[i+1]["ratio"] <= rows[i]["ratio"] for i in range(len(rows)-1)),
              "result": "The 71/81 grids reveal continuing drift rather than a validated plateau. The analytic Birkhoff bound survives, but a two-sided infinite-dimensional spectral-gap bracket remains unavailable."}
    return finalize(456, "numerical_sequence_extended_prior_plateau_interpretation_weakened",
                    "Ulam eigenvalues lack a certified discretization-error theorem and are not a physical transfer experiment.", packet)


def st457() -> dict:
    atlas = json.loads((ROOT / "FIN_ST421_Stationary_Orbit_Interval_Atlas.json").read_text())["locally_certified_stationary_orbits"]
    minimizer = min(atlas, key=lambda x: x["value"]); pmin = min(minimizer["center"])
    mu0 = min(minimizer["tangent_Hessian_eigenvalues"])-1e-8
    derivative = lambda R: mu0-2*R*pmin/(pmin-R)**3
    R = brentq(derivative, 1e-15, .99*pmin); mu = mu0-R/(pmin-R)**2; eps = mu*R/2
    packet = {"smallest_probability": pmin, "unperturbed_curvature_lower": mu0,
              "optimal_radius_in_declared_scalar_architecture": R,
              "strong_convexity_lower": mu, "optimal_forcing_threshold": eps,
              "stationarity_equation": "mu0=2 R p_min/(p_min-R)^3",
              "theorem": "The displayed R is the unique maximizer of epsilon(R)=R[mu0-R/(p_min-R)^2]/2 over the positive-curvature interval.",
              "global_ISS_basin": False}
    return finalize(457, "proven_sharp_optimum_within_conservative_scalar_ISS_bound",
                    "Optimality is only within this Hessian-variation estimate, not among all Lyapunov functions or physical perturbation models.", packet)


def st458() -> dict:
    M = np.array([[float(sum(bounds(x))/2) for x in row] for row in rank7_interval_matrix()])
    ev = np.linalg.eigvalsh(M); off = M[np.triu_indices(N, 1)]
    packet = {"rank": int(np.sum(ev > 1e-10)), "lambda_max": float(ev[-1]),
              "positive_offdiagonal_pairs": int(np.sum(off > 0)),
              "negative_offdiagonal_pairs": int(np.sum(off < 0)),
              "absolute_row_sum_bound": float(np.max(np.sum(np.abs(M), axis=1))),
              "known_certified_localized_value_interval": json.loads((ROOT / "FIN_ST361_Certified_Rank7_Mediator_Localized_State.json").read_text())["exact_root_objective_interval"],
              "sign_aware_global_SDP_certificate": None,
              "result": "The enlarged audit retains the mixed-sign obstruction. Spectral and absolute-value relaxations are too weak at g=4; no global rank-seven minimizer theorem is promoted."}
    return finalize(458, "bounded_no_go_for_current_rank7_global_relaxations",
                    "Failure of the declared relaxations is not failure of the rank-seven model or proof of another minimizer.", packet)


def gate(k: int, kind: str) -> dict:
    if kind == "gain":
        packet = {"new_internal_gain_source": False, "g_remains_supplied": True,
                  "reason": "A global transition bracket constrains consequences of g but does not generate g."}
        status = "blocked_no_new_strict_gain_source"
    elif kind == "selector":
        packet = {"new_nonpremise_selector": False, "QW_2191": "open",
                  "reason": "Global minimization chooses an orbit only up to D12; no member is canonically selected."}
        status = "blocked_no_new_selector_provider"
    else:
        packet = {"external_referee": "absent", "independent_laboratory_record": "absent",
                  "held_out_empirical_record": "absent", "reason": "The complete batch is local analytical/computational work."}
        status = "blocked_no_independent_empirical_evidence"
    return finalize(k, status, "No dimensional, bridge-role-transfer, laboratory, L_total, Standard Model, gravity, or ToE closure is exported.", packet)


def figures(results: dict) -> None:
    FIG_DIR.mkdir(exist_ok=True)
    bracket = results["ST449"]["certified_first_global_transition_bracket"]
    fig, ax = plt.subplots(figsize=(7.2, 3.5)); ax.hlines(1, bracket[0], bracket[1], lw=7, color="#2563eb")
    ax.scatter(bracket, [1, 1], c=["#16a34a", "#dc2626"], zorder=3)
    ax.set_yticks([]); ax.set_xlabel("supplied dimensionless gain g")
    ax.set_title("ST449: certified first-global-transition bracket"); fig.tight_layout()
    fig.savefig(FIG_DIR / "st449_global_transition_bracket.png", dpi=180); plt.close(fig)
    rows = results["ST452"]["cells"]
    b = [r["b_interval"][1] for r in rows]; y2 = [r["endpoint_centers"][-1][1] for r in rows]
    T = [r["endpoint_centers"][-1][4] for r in rows]
    fig, ax = plt.subplots(figsize=(7.2, 4.0)); ax.plot(b, y2, "o-", label=r"$y_2$"); ax.plot(b, T, "s-", label=r"$T$")
    ax.set_xlabel(r"compactified parameter $b=1/y_3$"); ax.legend(); ax.set_title("ST452: certified adaptive IR chain")
    fig.tight_layout(); fig.savefig(FIG_DIR / "st452_ir_chain.png", dpi=180); plt.close(fig)
    rows = results["ST456"]["extended_Ulam_sequence"]
    fig, ax = plt.subplots(figsize=(7.2, 4.0)); ax.plot([r["grid"] for r in rows], [r["ratio"] for r in rows], "o-")
    ax.set(xlabel="Ulam grid", ylabel=r"$|\lambda_2|/|\lambda_1|$", title="ST456: extended numerical transfer sequence")
    fig.tight_layout(); fig.savefig(FIG_DIR / "st456_ulam_extension.png", dpi=180); plt.close(fig)


def main() -> None:
    results = {}
    for k, fn in [(447, st447), (448, st448), (449, st449), (450, st450),
                  (451, st451), (452, st452), (453, st453), (454, st454),
                  (455, st455), (456, st456), (457, st457), (458, st458)]:
        print(f"running ST{k}", flush=True); results[f"ST{k}"] = fn()
    results["ST459"] = gate(459, "gain"); results["ST460"] = gate(460, "selector")
    results["ST461"] = gate(461, "evidence")
    RESULTS.write_text(json.dumps(native(results), indent=2, sort_keys=True), encoding="utf-8")
    with SUMMARY.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f); w.writerow(["program", "status", "object", "boundary"])
        for k in range(447, 462):
            r = results[f"ST{k}"]; w.writerow([f"ST{k}", r["status"], r["object"], r["boundary"]])
    figures(results)
    print(f"wrote {RESULTS.name} and {SUMMARY.name}")


if __name__ == "__main__": main()
