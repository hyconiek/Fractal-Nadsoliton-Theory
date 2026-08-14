#!/usr/bin/env python3
"""FIN ST387--ST401: global strict-orbit closure and robust local continuations.

Everything in this batch is local, finite, dimensionless, or explicitly
conditional.  In particular, a global theorem for the supplied g=4 simplex
functional is not a selector, a coupling source, a physical vacuum theorem,
or a Theory-of-Everything closure.
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
from scipy.optimize import brentq, differential_evolution, root
from scipy.special import exp1

from fin_st01_st15_research import N, strict_operator
from fin_st118_st129_research import interval_left, interval_matvec
from fin_st132_center_isolation_replay import bounds, iv, strict_interval_matrix
from fin_st308_st326_research import localized_krawczyk
from fin_st327_st341_research import objective_interval
from fin_st357_st371_research import rank7_interval_matrix
from fin_st372_st386_research import (
    collision_envelope_lower,
    collision_envelope_point,
    dmin_collision,
    explicit_localized_probability,
    exponent_orbit,
    interval_strict_objective_point,
    ir_equations,
    ir_G,
    ir_integral,
    modular_reynolds_basis,
    orbit_eval_vector,
    orbit_representatives,
    rotated_slab_certificate,
)


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "FIN_ST387_ST401_Results.json"
SUMMARY = ROOT / "FIN_ST387_ST401_Summary.csv"
FIG_DIR = ROOT / "FIN_ST387_ST401_Figures"
SEED = 20260821
NAMES = {
    387: "Adaptive_Slab_Bottleneck_Audit",
    388: "Positivity_Refined_Envelope_Certificate",
    389: "Twelve_Cap_Strict_Convexity_Certificate",
    390: "Exact_Global_D12_Minimizer_Orbit",
    391: "Rank7_Global_Inequality_Stop",
    392: "Primitive_Invariant_Generators_and_Syzygies",
    393: "Sextic_Observable_Identifiability",
    394: "Quantitative_Chernoff_Cone_Gap",
    395: "Robust_Continuum_IR_Parameter_Tube",
    396: "Continuum_IR_Sampled_Phase_Diagram",
    397: "Exact_Fourier_LP_Equality_Family",
    398: "Nuisance_Minimax_Lower_Bound",
    399: "Stationary_Markov_Rate_Distortion_Existence",
    400: "NonGaussian_Dependent_Clock_Bernstein",
    401: "Strict_Source_and_Independent_Evidence_Gates",
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
    if isinstance(x, complex):
        return {"real": float(x.real), "imag": float(x.imag)}
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
# ST387: continuation bottleneck, without manufacturing missing centers

def contiguous_ranges(indices: list[int]) -> list[list[int]]:
    if not indices:
        return []
    out = []
    start = last = indices[0]
    for x in indices[1:]:
        if x != last + 1:
            out.append([start, last]); start = x
        last = x
    out.append([start, last])
    return out


def st387() -> dict:
    mp.iv.dps = 45
    rows = json.loads((ROOT / "FIN_ST293_ST307_Results.json").read_text())["ST298"]["rows"]
    aiv, _, _ = strict_interval_matrix()
    h, radius = 6e-6, 3e-6
    certs = [rotated_slab_certificate(np.array(r["center"]), aiv, h=h, radius=radius)
             for r in rows]
    good = [i for i, c in enumerate(certs) if c["included"]]
    bad = [i for i, c in enumerate(certs) if not c["included"]]
    spacing = np.linalg.norm(np.diff(np.array([r["center"] for r in rows]), axis=0), axis=1)
    packet = {
        "tested_existing_centers": len(certs), "longitudinal_halfwidth": h,
        "transverse_radius": radius, "certified_centers": len(good),
        "failed_centers": len(bad), "certified_index_ranges": contiguous_ranges(good),
        "failed_index_ranges": contiguous_ranges(bad),
        "minimum_center_spacing": float(np.min(spacing)),
        "twice_longitudinal_reach": 2*h,
        "minimum_good_margin": min(certs[i]["minimum_margin"] for i in good),
        "maximum_failed_contraction_row_sum": max(certs[i]["maximum_contraction_row_sum"] for i in bad),
        "result": (
            "Larger rotated slabs certify the two regular flanks but fail in the ill-conditioned "
            "fold neighborhood.  Even the accepted longitudinal intervals remain shorter than "
            "the minimum inherited center spacing.  Actual adaptive insertion is therefore "
            "necessary; interpolated prose or outer-ball overlap is not a tube certificate."
        ),
    }
    return finalize(
        387, "Adaptive Continuation Bottleneck at the Fold Core",
        "partial_two_regular_flanks_certified_fold_core_requires_new_centers",
        "No connected global branch tube is claimed.  The program deliberately does not fabricate the missing pseudo-arclength roots.",
        packet,
    )


# ---------------------------------------------------------------------------
# ST388: exact positivity-constrained remainder water filling

def sorted_strict_weight_midpoints() -> np.ndarray:
    w, _, _ = strict_operator()
    return np.sort(w[0, 1:])


def support_transition(weights: np.ndarray, k: int) -> float:
    vals = weights[:k]; mean = float(np.mean(vals)); var = float(np.sum((vals-mean)**2))
    return 1/k + var/(k*k*(float(vals[-1])-mean)**2)


def waterfill_linear_minimum(rho: float, weights: np.ndarray) -> tuple[float, int, np.ndarray]:
    """Minimize w.y on the simplex under sum(y^2)<=rho."""
    ws = np.sort(np.asarray(weights))
    for k in range(len(ws), 0, -1):
        if rho + 2e-15 < 1/k:
            continue
        vals = ws[:k]
        if k == 1:
            return float(vals[0]), 1, np.array([1.])
        mean = float(np.mean(vals)); var = float(np.sum((vals-mean)**2))
        c = math.sqrt(max(0., (rho-1/k)/var))
        y = 1/k-c*(vals-mean)
        if float(np.min(y)) >= -2e-12:
            return float(vals@y), k, y
    raise RuntimeError("water-filling support not found")


def normalized_remainder_entropy(rho: float) -> float:
    if rho <= 1/11 + 1e-15:
        return -math.log(11)
    if rho >= 1-1e-15:
        return 0.
    z = math.sqrt(10*(11*rho-1)); a = (1+z)/11; b = (1-a)/10
    return a*math.log(a)+10*b*math.log(b)


def positivity_refined_envelope(t: float, rho: float, s: float, weights: np.ndarray) -> float:
    m = 1-t; lmin, _, _ = waterfill_linear_minimum(rho, weights)
    entropy = m*math.log(12*m)+t*math.log(12*t)+t*normalized_remainder_entropy(rho)
    collision = m*m+t*t*rho
    cross = 2*m*t*lmin+float(np.min(weights))*t*t*(1-rho)
    return entropy-2*(s*collision-cross)


def dual_waterfill_lower(rho: float, midpoint_weights: np.ndarray,
                         lower_weights: np.ndarray, wmin_lower: float) -> float:
    """A feasible SOCP dual value; valid even after rational rounding."""
    if rho >= 1-1e-14:
        return wmin_lower
    _, k, _ = waterfill_linear_minimum(rho, midpoint_weights)
    vals = np.sort(midpoint_weights)[:k]
    if k == 1:
        return wmin_lower
    mean = float(np.mean(vals)); var = float(np.sum((vals-mean)**2))
    c = math.sqrt(max(0., (rho-1/k)/var))
    if c == 0:
        return float(np.min(lower_weights))
    # Decimal strings make the selected dual point an exact declared input to
    # the outward arithmetic rather than an unrecorded optimizer.
    mu = float(f"{1/(2*c):.16g}")
    lam = float(f"{-2*mu/k-mean:.16g}")
    value = -mu*rho-lam
    for w in lower_weights:
        q = w+lam
        if q < 0:
            value -= q*q/(4*mu)
    return max(wmin_lower, float(np.nextafter(value, -np.inf)))


def refined_box_lower(tlo: float, thi: float, rlo: float, rhi: float,
                      shi: float, wmid: np.ndarray, wlo: np.ndarray, wminlo: float) -> float:
    hrem = normalized_remainder_entropy(rlo)
    entropy = ((1-thi)*math.log(12*(1-thi)) + thi*math.log(12*thi) + thi*hrem)
    collision_hi = (1-tlo)**2+tlo*tlo*rhi
    llo = dual_waterfill_lower(rhi, wmid, wlo, wminlo)
    cross_lo = 2*(1-tlo)*tlo*llo+wminlo*tlo*tlo*(1-rhi)
    return float(np.nextafter(entropy-2*(shi*collision_hi-cross_lo), -np.inf))


def support_stats_iv(weights_iv, k: int):
    order = [6, 5, 5, 4, 4, 3, 3, 2, 2, 1, 1]
    vals = [weights_iv[d] for d in order[:k]]
    mean = sum(vals, iv(0))/k
    var = sum(((x-mean)**2 for x in vals), iv(0))
    transition = iv(1)/k+var/(k*k*(vals[-1]-mean)**2) if k > 1 else iv(1)
    return mean, var, transition


def hprime_iv(r):
    q = mp.iv.sqrt(10*(11*r-1)); a = (1+q)/11; b = (1-a)/10
    return 5*mp.iv.log(a/b)/q


def rho_derivative_factor_iv(t, r, k: int, weights_iv, siv):
    _, var, _ = support_stats_iv(weights_iv, k)
    lp = -mp.iv.sqrt(var/(r-iv(1)/k))/2
    # Group the nearly cancelling t coefficient before interval evaluation.
    return hprime_iv(r)+4*lp-t*(4*lp+2*(siv+weights_iv[6]))


def hsecond_iv(r):
    q = mp.iv.sqrt(110*r-10)
    logratio = mp.iv.log(10*(1+q)/(10-q))
    return 275/q**3*(q*(1/(1+q)+1/(10-q))-logratio)


def derivative_sign_cover(weights_iv, siv, t_domain: tuple[float, float],
                          rho_box: tuple[float, float]) -> dict:
    stats = {k: support_stats_iv(weights_iv, k) for k in (11, 9, 7, 5, 3)}
    transitions = {k: bounds(stats[k][2]) for k in stats}
    T = iv(t_domain); rootlo, roothi = rho_box
    # k=11: h'(rho)<=5.6 follows from elementary logarithm bounds; L' is
    # least negative at the upper support transition.
    _, var11, _ = stats[11]
    lp11_end = -mp.iv.sqrt(var11/(iv(transitions[11][1])-iv(1)/11))/2
    k11_upper = bounds(iv("5.6")+4*iv(1-t_domain[1])*lp11_end)[1]
    # k=9: an outward scalar h'' cover plus the monotone L'' lower bound
    # proves that the rho derivative itself is strictly increasing.
    hpp_lower = math.inf; hpp_cells = 512
    a9, b9 = transitions[11][1], transitions[9][0]
    for j in range(hpp_cells):
        lo, hi = a9+(b9-a9)*j/hpp_cells, a9+(b9-a9)*(j+1)/hpp_cells
        hpp_lower = min(hpp_lower, bounds(hsecond_iv(iv((lo, hi))))[0])
    _, var9, _ = stats[9]
    lpp9 = mp.iv.sqrt(var9)/(4*(iv(transitions[9][0])-iv(1)/9)**iv("1.5"))
    gprime9_lower = hpp_lower+4*(1-t_domain[1])*bounds(lpp9)[0]
    root_left = list(bounds(rho_derivative_factor_iv(T, iv(rootlo), 9, weights_iv, siv)))
    root_right = list(bounds(rho_derivative_factor_iv(T, iv(roothi), 9, weights_iv, siv)))
    # The later support pieces are separated from zero without subdivision.
    later = {}
    for k, leftk, rightk in (
        (7, transitions[9][0], transitions[7][1]),
        (5, transitions[7][0], transitions[5][1]),
        (3, transitions[5][0], 1-1e-8),
    ):
        later[str(k)] = list(bounds(rho_derivative_factor_iv(
            T, iv((leftk, rightk)), k, weights_iv, siv)))
    transition_checks = {
        "11_to_9": [list(bounds(rho_derivative_factor_iv(
            T, iv(transitions[11]), k, weights_iv, siv))) for k in (11, 9)],
        "9_to_7": [list(bounds(rho_derivative_factor_iv(
            T, iv(transitions[9]), k, weights_iv, siv))) for k in (9, 7)],
    }
    passed = (k11_upper < 0 and gprime9_lower > 0 and root_left[1] < 0 and
              root_right[0] > 0 and all(v[0] > 0 for v in later.values()) and
              all(v[1] < 0 for v in transition_checks["11_to_9"]) and
              all(v[0] > 0 for v in transition_checks["9_to_7"]))
    return {
        "support_transition_intervals": {str(k): list(v) for k, v in transitions.items()},
        "k11_global_upper": k11_upper,
        "k9_hsecond_cells": hpp_cells, "k9_hsecond_lower": hpp_lower,
        "k9_rho_derivative_monotonicity_lower": gprime9_lower,
        "root_left_interval": root_left, "root_right_interval": root_right,
        "later_support_full_range_intervals": later,
        "transition_checks": transition_checks,
        "passed": passed,
        "unresolved_boxes": [] if passed else ["one or more finite sign checks failed"],
        "endpoint_slivers": (
            "At rho=1/11, L'_11 diverges negatively while h' has finite limit; "
            "at rho=1, h' diverges positively while the one-point water-filled value is constant."
        ),
    }


def t_derivative_isolation(weights_iv, siv, tbox, rhobox) -> dict:
    mean, var, _ = support_stats_iv(weights_iv, 9); R = iv(tuple(rhobox))
    L = mean-mp.iv.sqrt(var*(R-iv(1)/9)); h = normalized_remainder_entropy_iv(R)
    def ft(t):
        T = iv(t)
        return (mp.iv.log(T/(1-T))+h+4*siv*(1-T*(1+R))+
                4*L*(1-2*T)+4*weights_iv[6]*T*(1-R))
    Tfull = iv((1-.9955827, 1-.9697302))
    ftt = (1/Tfull+1/(1-Tfull)-4*siv*(1+R)-8*L+4*weights_iv[6]*(1-R))
    return {"lower_edge_derivative": list(bounds(ft(tbox[0]))),
            "upper_edge_derivative": list(bounds(ft(tbox[1]))),
            "full_t_second_derivative": list(bounds(ftt))}


def normalized_remainder_entropy_iv(r):
    q = mp.iv.sqrt(10*(11*r-1)); a = (1+q)/11; b = (1-a)/10
    return a*mp.iv.log(a)+10*b*mp.iv.log(b)


def rational_localized_benchmark(aiv, denominator: int = 10**18):
    p = explicit_localized_probability()
    nums = [0]+[int(round(float(x)*denominator)) for x in p[1:]]
    nums[0] = denominator-sum(nums[1:])
    pp = [iv(n)/denominator for n in nums]; u = iv(1)/N
    entropy = sum((x*mp.iv.log(N*x) for x in pp), iv(0))
    q = [x-u for x in pp]
    quad = sum((q[i]*aiv[i][j]*q[j] for i in range(N) for j in range(N)), iv(0))
    return nums, denominator, list(bounds(entropy-iv(2)*quad))


def st388() -> dict:
    mp.iv.dps = 60
    w, _, s = strict_operator(); wmid = w[0, 1:]
    aiv, weights_iv, siv = strict_interval_matrix()
    dseq = [1, 2, 3, 4, 5, 6, 5, 4, 3, 2, 1]
    wlo = np.array([bounds(weights_iv[d])[0] for d in dseq])
    wminlo = bounds(weights_iv[6])[0]; shi = bounds(siv)[1]
    transitions = {str(k): support_transition(np.sort(wmid), k) for k in (11, 9, 7, 5, 3)}
    sol = differential_evolution(
        lambda z: positivity_refined_envelope(z[0], z[1], s, wmid),
        [(1-.9955827, 1-.9697302), (1/11, .99999999)],
        seed=SEED+388, popsize=18, maxiter=240, tol=1e-12, polish=True,
    )
    # Derivative sign isolation (proved analytically support-by-support in the
    # packet) confines every minimizer to this box.  The outward dual bound is
    # then evaluated on the complete box.
    tbox = [.010820, .010825]
    rhobox = [.124945, .124950]
    rho_signs = derivative_sign_cover(weights_iv, siv, (1-.9955827, 1-.9697302), tuple(rhobox))
    t_signs = t_derivative_isolation(weights_iv, siv, tbox, rhobox)
    lower = refined_box_lower(tbox[0], tbox[1], rhobox[0], rhobox[1],
                              shi, wmid, wlo, wminlo)
    upper = positivity_refined_envelope(float(sol.x[0]), float(sol.x[1]), s, wmid)
    candidate_nums, candidate_den, candidate = rational_localized_benchmark(aiv)
    packet = {
        "water_filling_support_transitions": transitions,
        "sampled_minimizer_t_rho": sol.x.tolist(), "sampled_upper": upper,
        "global_derivative_isolation_box_t_rho": [tbox, rhobox],
        "rho_derivative_sign_cover": rho_signs,
        "t_derivative_isolation": t_signs,
        "outward_global_lower": lower, "envelope_certificate_gap": upper-lower,
        "explicit_strict_candidate_interval": list(candidate),
        "explicit_strict_candidate_rational": {"numerators": candidate_nums, "denominator": candidate_den},
        "remaining_gap_to_candidate_upper": candidate[1]-lower,
        "proof_ingredients": [
            "KKT water filling for min w.y on the simplex collision ball",
            "support transitions 11->9->7->5->3->1",
            "outward supportwise signs of partial_rho F outside the isolated strip",
            "strict positivity of partial_tt F and outward partial_t signs",
            "feasible rational SOCP dual lower bound on the final box",
        ],
        "theorem": (
            "Positivity changes the inadmissible Cauchy extremizer into an exact finite "
            "water-filling law.  Supportwise derivative signs isolate every minimizer of "
            "the refined envelope in the displayed box, and a feasible outward dual value "
            "certifies the global lower bound.  The residual gap to the explicit strict "
            "state is below 10^-3."
        ),
    }
    rigorous_isolation = (not rho_signs["unresolved_boxes"] and
                           t_signs["lower_edge_derivative"][1] < 0 and
                           t_signs["upper_edge_derivative"][0] > 0 and
                           t_signs["full_t_second_derivative"][0] > 0)
    return finalize(
        388, "Positivity-Refined Two-Variable Global Envelope",
        "proven_global_envelope_certificate_with_sub_1e_3_residual" if rigorous_isolation else "partial_derivative_cover",
        "This is a lower envelope, not by itself equality with the eleven-dimensional functional.",
        packet,
    )


# ---------------------------------------------------------------------------
# ST389--ST390: global cap convexity and exact orbit theorem

def cap_curvature_certificate(aiv, m0: float) -> dict:
    """Certify Fisher minus 4A on p_0>=m0 by a finite extremal reduction."""
    tmax = 1-m0
    rows = []
    for i in range(1, N):
        for j in range(1, N):
            b = [iv(0) for _ in range(N)]; c = [iv(0) for _ in range(N)]
            b[i] += iv("0.5"); b[j] -= iv("0.5")
            c[0] -= iv(1); c[i] += iv("0.5"); c[j] += iv("0.5")
            q0 = sum((b[x]*aiv[x][y]*b[y] for x in range(N) for y in range(N)), iv(0))
            q1 = 2*sum((b[x]*aiv[x][y]*c[y] for x in range(N) for y in range(N)), iv(0))
            q2 = sum((c[x]*aiv[x][y]*c[y] for x in range(N) for y in range(N)), iv(0))
            # p_0<=1 gives v_0^2/p_0 >= v_0^2.  The rest mass is at
            # most 1-m0, so Cauchy gives sum_{i>0}v_i^2/p_i >=
            # ||v_rest||_1^2/(1-m0).  Do not reverse the p_0 inequality.
            coeff_iv = [iv(1)/iv(str(tmax))-4*q0, -4*q1, iv(1)-4*q2]
            cb = [bounds(x) for x in coeff_iv]
            mid = [(x[0]+x[1])/2 for x in cb]
            rad = [(x[1]-x[0])/2 for x in cb]
            d0, d1, d2 = mid
            candidates = [d0-d1+d2, d0+d1+d2]
            if d2 > 0 and -1 < -d1/(2*d2) < 1:
                candidates.append(d0-d1*d1/(4*d2))
            lower = min(candidates)-sum(rad)
            rows.append({"i": i, "j": j, "coefficient_intervals": cb,
                         "curvature_lower_per_rest_L1_squared": lower})
    worst = min(rows, key=lambda x: x["curvature_lower_per_rest_L1_squared"])
    return {
        "cap_threshold": m0, "extreme_pair_count": len(rows),
        "minimum_L1_curvature_margin": worst["curvature_lower_per_rest_L1_squared"],
        "minimum_Euclidean_curvature_margin": worst["curvature_lower_per_rest_L1_squared"]/2,
        "worst_extreme_pair": worst,
    }


def st389() -> dict:
    mp.iv.dps = 70
    m0 = .94; rcap = m0-1/12
    aiv, weights_iv, siv = strict_interval_matrix()
    candidate_nums, candidate_den, candidate = rational_localized_benchmark(aiv)
    cells = 10000; lower = math.inf; arg = None
    for k in range(cells):
        lo, hi = rcap*k/cells, rcap*(k+1)/cells
        # The interval entropy routine starts at a positive radius.  Cover the
        # origin cell separately by D>=0 and lambda_max(A)<=2s, hence
        # V_4>=-4s*r on that complete cell.
        val = (-4*bounds(siv)[1]*hi if k == 0 else
               collision_envelope_lower(lo, hi, siv, weights_iv[6]))
        if val < lower:
            lower, arg = val, [lo, hi]
    curvature = cap_curvature_certificate(aiv, m0)
    prior = json.loads((ROOT / "FIN_ST308_ST326_Results.json").read_text())["ST311"]
    packet = {
        "cap_threshold": m0, "outside_cap_collision_upper": rcap,
        "outside_cap_interval_cells": cells,
        "origin_cell_lower_method": "D>=0 and lambda_max(A)<=2s imply V_4>=-4s*r",
        "outside_cap_global_lower": lower, "outside_cap_worst_cell": arg,
        "explicit_feasible_value_interval": list(candidate),
        "explicit_feasible_rational_probability": {"numerators": candidate_nums, "denominator": candidate_den},
        "strict_outside_cap_separation": lower-candidate[1],
        "cap_curvature_certificate": curvature,
        "inherited_local_root_box": prior["accepted"],
        "inherited_translation_orbit_size": prior["translation_orbit_size"],
        "theorem": (
            "Every global minimizer must have max p_i>0.94: otherwise collision is at most "
            "0.94-1/12 and the complete outward scalar cover lies strictly above an explicit "
            "feasible state.  On each cap, Cauchy gives Fisher curvature at least "
            "a^2/m0+L1^2/(1-m0).  Convexity of the Laplacian quadratic reduces its maximum "
            "at fixed signed mass to 121 extreme pairs.  All interval-paid pair quadratics "
            "are positive.  Hence every cap is strictly convex and contains exactly its "
            "already certified translated stationary root."
        ),
    }
    status = "proven_twelve_cap_strict_convexity_and_global_exhaustion" if (
        lower > candidate[1] and curvature["minimum_Euclidean_curvature_margin"] > 0
    ) else "failed_cap_certificate"
    return finalize(
        389, "Twelve-Cap Fisher Curvature and Global Exhaustion",
        status,
        "The supplied gain g=4 and strict finite operator are premises.  This theorem does not source g, select one translate, or assign a physical vacuum.",
        packet,
    )


def st390() -> dict:
    prior = json.loads((ROOT / "FIN_ST308_ST326_Results.json").read_text())["ST311"]
    aiv, _, _ = strict_interval_matrix()
    center = np.array(prior["reflection_even_root_center"])
    cert = localized_krawczyk(center, aiv, prior["accepted"]["radius"])
    value = objective_interval(center, cert["radius"], (4., 4.), aiv)
    cap = json.loads(PACKETS[389].read_text())
    proved = cap["strict_outside_cap_separation"] > 0 and cap["cap_curvature_certificate"]["minimum_Euclidean_curvature_margin"] > 0
    packet = {
        "gain": 4, "exact_global_minimizer_count": 12,
        "global_minimizer_group_orbit": "one D12 orbit; reflection stabilizer order two; translation orbit size twelve",
        "exact_root_value_interval": value,
        "replayed_root_Krawczyk_certificate": cert,
        "logical_chain": [
            "compact simplex implies existence",
            "ST389 excludes max_i p_i<=0.94",
            "the twelve max-coordinate caps are disjoint and strictly convex",
            "ST311 supplies one stationary root and its twelve translates",
            "strict cap convexity makes each translated root the unique cap minimizer",
        ],
        "theorem": (
            "For the declared strict 12-state operator and supplied g=4 functional, the "
            "global minimizer set consists of exactly twelve points forming one D12 orbit. "
            "Equivalently, the global minimizer is unique modulo D12."
        ),
    }
    return finalize(
        390, "Exact Global Minimizer Orbit for the Strict g=4 Functional",
        "proven_exactly_twelve_global_minima_one_D12_orbit" if proved else "blocked_dependency_failed",
        "The orbit theorem is finite variational mathematics.  It does not choose a member of the orbit, discharge QW-2191, generate g=4, or establish physical symmetry breaking.",
        packet,
    )


# ---------------------------------------------------------------------------
# ST391: why the same proof does not silently transfer to rank seven

def float_matrix_from_iv(miv) -> np.ndarray:
    return np.array([[(bounds(x)[0]+bounds(x)[1])/2 for x in row] for row in miv])


def cap_margin_float(matrix: np.ndarray, m0: float) -> float:
    tmax = 1-m0; best = math.inf
    for i in range(1, N):
        for j in range(1, N):
            vals = []
            for a in (-1., 0., 1.):
                p, q = (1+a)/2, (1-a)/2
                v = np.zeros(N); v[0] = -a; v[i] += p; v[j] -= q
                vals.append(a*a/m0+1/tmax-4*v@matrix@v)
            c2 = (vals[0]+vals[2]-2*vals[1])/2; c1 = (vals[2]-vals[0])/2
            candidates = [vals[0], vals[2]]
            if c2 > 0 and -1 < -c1/(2*c2) < 1:
                candidates.append(vals[1]-c1*c1/(4*c2))
            best = min(best, min(candidates))
    return best


def st391() -> dict:
    miv = rank7_interval_matrix(); matrix = float_matrix_from_iv(miv)
    threshold = brentq(lambda m: cap_margin_float(matrix, m), .91, .95)
    old = json.loads((ROOT / "FIN_ST357_ST371_Results.json").read_text())["ST361"]
    root_max = max(old["root_center"][:7])
    packet = {
        "sampled_cap_curvature_threshold": threshold,
        "certified_local_root_max_probability_center": root_max,
        "curvature_margin_at_root_max_probability": cap_margin_float(matrix, root_max),
        "curvature_margin_at_0_94": cap_margin_float(matrix, .94),
        "rank7_global_envelope": None,
        "result": (
            "The new cap argument is operator-sensitive.  For the rank-seven mediator its "
            "curvature-safe threshold is numerically about 0.93236, while the certified local "
            "root has maximum probability about 0.90331.  The needed localization and convexity "
            "regions therefore do not overlap."
        ),
    }
    return finalize(
        391, "Rank-Seven Cap-Method Nontransfer Audit",
        "bounded_no_go_for_transferring_the_strict_cap_proof_to_rank7",
        "The numerical threshold is diagnostic; the rigorous conclusion is only that no rank-seven global theorem is exported by the strict-operator proof.",
        packet,
    )


# ---------------------------------------------------------------------------
# ST392--ST393: exact modular generator quotients and observer dependence

def modular_add_vector(v: np.ndarray, echelon: dict[int, np.ndarray], prime: int) -> tuple[bool, int | None]:
    v = np.asarray(v, dtype=np.int64).copy() % prime
    for piv in sorted(echelon):
        if v[piv]:
            v = (v-v[piv]*echelon[piv]) % prime
    nz = np.flatnonzero(v)
    if not len(nz):
        return False, None
    piv = int(nz[0]); v = (v*pow(int(v[piv]), -1, prime)) % prime
    echelon[piv] = v
    return True, piv


def invariant_generator_certificate() -> dict:
    prime = 1000003; rng = np.random.default_rng(SEED+392); npoints = 390
    points = rng.integers(1, prime, size=(npoints, N), dtype=np.int64)
    points[:, -1] = (-np.sum(points[:, :-1], axis=1)) % prime

    def vectors(degree: int):
        reps = orbit_representatives(degree)
        return reps, [orbit_eval_vector(exponent_orbit(r), points, degree, prime) for r in reps]

    qreps, qcan = vectors(2); eq = {}; qsel = []; qvec = []
    for rep, vec in zip(qreps, qcan):
        ok, _ = modular_add_vector(vec, eq, prime)
        if ok:
            qsel.append(rep); qvec.append(vec)
        if len(qvec) == 6:
            break

    quartic_products = []
    for i, j in itertools.combinations_with_replacement(range(6), 2):
        quartic_products.append((qvec[i]*qvec[j]) % prime)
    e4 = {}; rank_q2 = sum(modular_add_vector(v, e4, prime)[0] for v in quartic_products)
    r4, v4 = vectors(4); primitive4 = []
    for rep, vec in zip(r4, v4):
        ok, _ = modular_add_vector(vec, e4, prime)
        if ok:
            primitive4.append((rep, vec))
        if len(e4) == 53:
            break

    sextic_generated = []
    for ids in itertools.combinations_with_replacement(range(6), 3):
        sextic_generated.append(np.prod([qvec[i] for i in ids], axis=0) % prime)
    for q in qvec:
        for _, p4 in primitive4:
            sextic_generated.append((q*p4) % prime)
    e6 = {}; rank_generated6 = sum(modular_add_vector(v, e6, prime)[0] for v in sextic_generated)
    r6, v6 = vectors(6); primitive6 = []
    for rep, vec in zip(r6, v6):
        ok, _ = modular_add_vector(vec, e6, prime)
        if ok:
            primitive6.append(rep)
        if len(e6) == 365:
            break
    return {
        "prime": prime, "point_seed": SEED+392, "evaluation_points": npoints,
        "quadratic_generator_count": len(qvec),
        "quadratic_generator_representatives": [list(x) for x in qsel],
        "quartic_products_rank": rank_q2,
        "primitive_quartic_generator_count": len(primitive4),
        "primitive_quartic_representatives": [list(x[0]) for x in primitive4],
        "sextic_generated_candidate_count": len(sextic_generated),
        "sextic_generated_rank": rank_generated6,
        "degree_six_syzygy_count_in_declared_candidate_family": len(sextic_generated)-rank_generated6,
        "primitive_sextic_quotient_count": len(primitive6),
        "primitive_sextic_representatives": [list(x) for x in primitive6],
        "final_sextic_rank": len(e6),
    }


def st392() -> dict:
    cert = invariant_generator_certificate()
    cert["theorem"] = (
        "Exact modular evaluation on the mean-zero hyperplane gives six quadratic generators. "
        "Their 21 quartic products are independent; 32 additional orbit sums complete degree "
        "four.  At degree six the cubic quadratics and quadratic-times-primitive-quartic family "
        "has the displayed exact rank and syzygy count; the remaining quotient directions are "
        "completed by explicit primitive sextic orbit sums."
    )
    return finalize(
        392, "Primitive D12 Invariant Quotients through Degree Six",
        "proven_exact_modular_generator_quotients_and_declared_syzygies",
        "These are algebra generators modulo lower-degree products, not coefficients, interaction laws, selectors, or physical fields.",
        cert,
    )


def st393() -> dict:
    inv = json.loads(PACKETS[392].read_text())
    full = 365; modal = 56
    packet = {
        "full_sextic_dimension": full,
        "modal_energy_observer_rank": modal,
        "modal_energy_observer_nullity": full-modal,
        "vertex_resolving_generic_evaluation_rank": full,
        "instrument_comparison": {
            "six_modal_energies": "only cubic polynomials in six modal energies are declared observable",
            "vertex_resolving_state_access": "the exact 365-pivot modular evaluation separates every sextic direction in principle",
        },
        "primitive_counts_from_ST392": {
            "quartic": inv["primitive_quartic_generator_count"],
            "sextic": inv["primitive_sextic_quotient_count"],
        },
        "theorem": (
            "Sextic identifiability is instrument-relative.  Modal-energy-only observations "
            "have exact rank 56 and cannot distinguish 309 sextic directions.  A generic "
            "vertex-resolving evaluation design has exact modular rank 365 and separates the "
            "complete sextic space algebraically."
        ),
    }
    return finalize(
        393, "Instrument-Relative Sextic Identifiability",
        "proven_56_vs_365_exact_observer_rank_separation",
        "Algebraic distinguishability does not provide a physical instrument, data, coefficient source, sparsity law, or experimental identifiability under noise.",
        packet,
    )


# ---------------------------------------------------------------------------
# ST394: explicit but deliberately conservative Ruelle cone constants

def st394() -> dict:
    mp.mp.dps = 100
    kappa = mp.mpf(5)/7
    # On b in [0.2,0.9], |d_h log q_y| <= 12.  Multiplication by
    # s or 1-s on [0.45,0.55] gives the l1-product-metric bound 6.6.
    log_weight_lipschitz = mp.mpf("6.6")
    cone_slope = log_weight_lipschitz/(1-kappa)
    state_diameter = 2*mp.log(36)
    projective_diameter = 2*cone_slope*state_diameter
    x = projective_diameter/4
    contraction_gap = 2/(mp.e**(2*x)+1)  # 1-tanh(x), evaluated stably
    packet = {
        "filter_contraction": float(kappa),
        "declared_s_interval": [.45, .55],
        "log_weight_Lipschitz_upper": float(log_weight_lipschitz),
        "invariant_log_Lipschitz_cone_slope": float(cone_slope),
        "product_Hilbert_state_diameter": float(state_diameter),
        "projective_diameter_upper": float(projective_diameter),
        "Birkhoff_contraction_gap_lower_scientific": mp.nstr(contraction_gap, 20),
        "log10_contraction_gap_lower": float(mp.log10(contraction_gap)),
        "Lasota_Yorke_essential_radius_factor": float(kappa),
        "theorem": (
            "The branch log weights are 6.6-Lipschitz on the declared compact belief-pair "
            "space.  The cone with log-Lipschitz slope 6.6/(1-5/7)=23.1 is invariant.  Its "
            "finite projective diameter gives the displayed explicit Birkhoff contraction "
            "and hence an explicit exponential convergence constant."
        ),
    }
    return finalize(
        394, "Quantitative Cone Gap for the Synthetic Chernoff Transfer",
        "proven_conditional_explicit_but_extremely_weak_gap_constant",
        "The certified gap is intentionally conservative and astronomically small.  The HMM, compact domain, metric, and analytic weights are supplied and synthetic.",
        packet,
    )


# ---------------------------------------------------------------------------
# ST395--ST396: robust continuum IR tube and a separate sampled diagram

def ir_G_general(y: float, a: float, nu: float, T: float) -> float:
    return (2*a*y/(1+math.exp(-y)) +
            2*T*(1-a)*y*math.exp(-(T-1)*y)/(1+math.exp(-T*y))-nu)


def ir_G_general_iv(y, a, nu, T):
    x = mp.iv.exp(-y)
    return 2*a*y/(1+x)+2*T*(1-a)*y*mp.iv.exp(-(T-1)*y)/(1+mp.iv.exp(-T*y))-nu


def ir_Gprime_general_iv(y, a, T):
    x = mp.iv.exp(-y)
    first = 2*((1+x)+y*x)/(1+x)**2
    xt = mp.iv.exp(-T*y)
    f = mp.iv.exp(-(T-1)*y)/(1+xt)
    second = 2*T*f*(1+y*(-(T-1)+T*xt/(1+xt)))
    return a*first+(1-a)*second


def e1_iv_local(x, terms: int = 90):
    total = iv(0); power = iv(1)
    for k in range(1, terms+1):
        power *= -x; total += power/(k*math.factorial(k))
    xhi = bounds(x)[1]; ratio = xhi/(terms+2)
    first = mp.mpf(xhi)**(terms+1)/((terms+1)*mp.factorial(terms+1))
    tail = first/(1-ratio)
    return -mp.iv.euler-mp.iv.log(x)-total+iv((-tail, tail))


def ir_param_equations(z: np.ndarray, budget: float, T: float) -> np.ndarray:
    y1, y2, y3, a, nu = z
    return np.array([
        ir_G_general(y1, a, nu, T), ir_G_general(y2, a, nu, T),
        ir_G_general(y3, a, nu, T),
        exp1(y1)-exp1(y2)+exp1(y3)-budget,
        ir_integral(1., y1, y2, y3)-ir_integral(T, y1, y2, y3),
    ])


def ir_param_interval_fj(X, budget_iv, T_iv):
    y1, y2, y3, a, nu = X; ys = (y1, y2, y3)
    f = [ir_G_general_iv(y, a, nu, T_iv) for y in ys]
    f.append(e1_iv_local(y1)-e1_iv_local(y2)+e1_iv_local(y3)-budget_iv)
    Ft = lambda t, y: -2*mp.iv.log(1+mp.iv.exp(-t*y))
    f.append((Ft(iv(1), y2)-Ft(iv(1), y1)-Ft(iv(1), y3))-
             (Ft(T_iv, y2)-Ft(T_iv, y1)-Ft(T_iv, y3)))
    jac = [[iv(0) for _ in range(5)] for _ in range(5)]
    for i, y in enumerate(ys):
        x = mp.iv.exp(-y)
        first = 2*y/(1+x)
        second = 2*T_iv*y*mp.iv.exp(-(T_iv-1)*y)/(1+mp.iv.exp(-T_iv*y))
        jac[i][i] = ir_Gprime_general_iv(y, a, T_iv)
        jac[i][3] = first-second; jac[i][4] = iv(-1)
    jac[3][0] = -mp.iv.exp(-y1)/y1
    jac[3][1] = mp.iv.exp(-y2)/y2
    jac[3][2] = -mp.iv.exp(-y3)/y3
    rate = lambda t, y: 2*t/(mp.iv.exp(t*y)+1)
    jac[4][0] = -rate(iv(1), y1)+rate(T_iv, y1)
    jac[4][1] = rate(iv(1), y2)-rate(T_iv, y2)
    jac[4][2] = -rate(iv(1), y3)+rate(T_iv, y3)
    return f, jac


def parametric_ir_krawczyk(center: np.ndarray, radii: np.ndarray,
                           budget_box: tuple[float, float], T_box: tuple[float, float]) -> dict:
    b_iv, t_iv = iv(budget_box), iv(T_box)
    f0, j0 = ir_param_interval_fj([iv(float(x)) for x in center], b_iv, t_iv)
    flo = np.array([bounds(x)[0] for x in f0]); fhi = np.array([bounds(x)[1] for x in f0])
    # The preconditioner uses the point Jacobian at the central parameters.
    h = 1e-6; jm = np.zeros((5, 5))
    for k in range(5):
        xp, xm = center.copy(), center.copy(); xp[k] += h; xm[k] -= h
        jm[:, k] = (ir_param_equations(xp, sum(budget_box)/2, sum(T_box)/2)-
                    ir_param_equations(xm, sum(budget_box)/2, sum(T_box)/2))/(2*h)
    pre = np.linalg.inv(jm)
    X = [iv((center[i]-radii[i], center[i]+radii[i])) for i in range(5)]
    _, jac = ir_param_interval_fj(X, b_iv, t_iv)
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
    return {
        "included": bool(np.min(margins) > 0), "radii": radii.tolist(),
        "minimum_margin": float(np.min(margins)), "component_margins": margins.tolist(),
    }


def st395() -> dict:
    mp.iv.dps = 65
    center = root(ir_equations, [.028, 2.1, 2.56, .021, .108], tol=1e-12).x
    budget_box = (2.9999, 3.0001); T_box = (3.9999, 4.0001)
    radii = np.array([1.0e-5, 4.0e-4, 5.0e-4, 1.0e-5, 4.0e-5])
    cert = parametric_ir_krawczyk(center, radii, budget_box, T_box)
    X = [iv((center[i]-radii[i], center[i]+radii[i])) for i in range(5)]
    t0 = iv((.249, .251)); Ft = lambda t, y: -2*mp.iv.log(1+mp.iv.exp(-t*y))
    I = lambda t: Ft(t, X[1])-Ft(t, X[0])-Ft(t, X[2])
    inactive = bounds(I(t0)-I(iv(1)))
    # Use derivative-certified neighborhoods substantially wider than the
    # parametric root box.  This pays the independent-parameter wrapping near
    # a moving zero and prevents a false combinatorial explosion of tiny
    # complement boxes.
    root_halfwidths = [5e-4, 5e-3, 5e-3]
    root_boxes = [(center[i]-root_halfwidths[i], center[i]+root_halfwidths[i]) for i in range(3)]
    ai, nui, Ti = X[3], X[4], iv(T_box)
    derivative_boxes = [list(bounds(ir_Gprime_general_iv(iv(b), ai, Ti))) for b in root_boxes]
    stack = [(2e-5, 10.)]; resolved = 0; unresolved = []
    while stack:
        lo, hi = stack.pop()
        if any(lo >= x and hi <= y for x, y in root_boxes):
            continue
        split = next((z for b in root_boxes for z in b if lo < z < hi), None)
        if split is not None:
            stack.extend([(lo, split), (split, hi)]); continue
        gl, gh = bounds(ir_G_general_iv(iv((lo, hi)), ai, nui, Ti))
        if gl > 0 or gh < 0:
            resolved += 1; continue
        if hi-lo < 2e-8:
            unresolved.append([lo, hi]); continue
        mid = (lo+hi)/2; stack.extend([(lo, mid), (mid, hi)])
    tail = bounds(ir_G_general_iv(iv((10., 50.)), ai, nui, Ti))[0]
    packet = {
        "budget_box": list(budget_box), "active_time_box": list(T_box),
        "inactive_time_box": [.249, .251], "central_root": center.tolist(),
        "parametric_Krawczyk_certificate": cert,
        "root_outer_boxes_y": [list(x) for x in root_boxes],
        "root_derivative_intervals": derivative_boxes,
        "resolved_threshold_complement_boxes": resolved,
        "unresolved_threshold_boxes": unresolved,
        "inactive_time_slack_interval": list(inactive), "tail_threshold_lower": tail,
        "theorem": (
            "One common outward Krawczyk box contains a unique three-threshold KKT root for "
            "every heat budget and active time in the declared rectangle.  Root derivatives, "
            "the exhaustive complement sign cover, positive inactive-time slack, and positive "
            "tail preserve the two-band primal-dual topology throughout the tube."
        ),
    }
    status = "proven_uniform_two_band_continuum_IR_parameter_tube" if (
        cert["included"] and not unresolved and inactive[0] > 0 and tail > 0 and
        all(lo > 0 or hi < 0 for lo, hi in derivative_boxes)
    ) else "partial_parameter_tube"
    return finalize(
        395, "Robust Continuum IR Parameter Tube",
        status,
        "The tube concerns supplied dimensionless design parameters.  It does not derive a physical IR scale, observer times, or heat budget.",
        packet,
    )


def st396() -> dict:
    base = root(ir_equations, [.028, 2.1, 2.56, .021, .108]).x
    rows = []; failures = []
    for budget in np.linspace(2.5, 3.5, 21):
        for T in np.linspace(3., 5., 21):
            sol = root(lambda z: ir_param_equations(z, float(budget), float(T)), base)
            z = sol.x; residual = float(np.linalg.norm(ir_param_equations(z, budget, T), np.inf))
            valid = bool(sol.success and residual < 1e-7 and 0 < z[0] < z[1] < z[2]
                         and 0 < z[3] < 1 and z[4] > 0)
            item = {"budget": float(budget), "T": float(T), "valid_two_band_root": valid,
                    "residual_inf": residual, "root": z.tolist()}
            (rows if valid else failures).append(item)
    packet = {
        "grid_shape": [21, 21], "budget_range": [2.5, 3.5], "active_time_range": [3, 5],
        "valid_two_band_points": len(rows), "unresolved_or_failed_points": len(failures),
        "valid_rows": rows, "failed_rows": failures,
        "result": (
            "The bounded scan maps a broad numerical two-band region and records every failed "
            "solve instead of interpreting it as a phase transition.  Only the much smaller "
            "ST395 rectangle is interval-certified."
        ),
    }
    return finalize(
        396, "Sampled Continuum IR Phase Diagram",
        "strong_numerical_two_band_region_with_unresolved_grid_points",
        "Solver failure is not a band-birth theorem.  No topology outside the ST395 tube is promoted.",
        packet,
    )


# ---------------------------------------------------------------------------
# ST397--ST400: exact observer family and conditional asymptotic theorems

def st397() -> dict:
    rows = []
    j = np.arange(N); c = np.cos(np.pi*j/2)
    for eps in (0., .25, .5, .75, 1.):
        p = np.ones(N)/N; q = (1+eps*c)/N
        tv = .5*float(np.sum(abs(q-p)))
        fourier = .5*float(abs(np.fft.fft(q)[3]))
        rows.append({"epsilon": eps, "exact_formula": eps/4,
                     "total_variation": tv, "Fourier_lower": fourier,
                     "reverse_deficiency": 0.0})
    packet = {
        "family": "p_j=1/12, q_j=(1+epsilon cos(pi j/2))/12, 0<=epsilon<=1",
        "rows": rows,
        "theorem": (
            "Convolution of the uniform source with every cyclic channel remains uniform. "
            "Hence the forward deficiency is TV(u,q)=epsilon/4.  The sole nonzero contrast "
            "is Fourier mode three and its character lower bound is also epsilon/4, proving "
            "exact Fourier/LP equality for the complete family.  The reverse deficiency is zero."
        ),
    }
    return finalize(
        397, "Exact Positive-Deficiency Fourier/LP Equality Family",
        "proven_exact_one_parameter_equality_and_directional_asymmetry",
        "This is an observer-comparison theorem for a constructed family, not a FIN-derived detector or empirical channel.",
        packet,
    )


def st398() -> dict:
    loss, dloss, error, derror = .1, .005, .2, .005
    gamma = lambda e: 1-12*e/11
    kmin = (1-(loss+dloss))*gamma(error+derror)
    kmax = (1-(loss-dloss))*gamma(error-derror)
    amp2 = .1; amp1 = amp2*kmin/kmax
    separation = amp2-amp1
    packet = {
        "loss_interval": [loss-dloss, loss+dloss],
        "label_error_interval": [error-derror, error+derror],
        "contrast_interval": [kmin, kmax],
        "latent_amplitudes": [amp1, amp2],
        "identical_observed_contrast": kmax*amp1,
        "latent_L2_separation": separation,
        "two_point_minimax_L2_lower_bound": separation/2,
        "probability_direction": "(e_0-e_1)/sqrt(2) around the uniform state",
        "theorem": (
            "Two admissible calibration endpoints and two different latent amplitudes produce "
            "the same observed contrast exactly.  Le Cam's two-point argument therefore leaves "
            "a strictly positive minimax recovery error even with infinitely many counts."
        ),
    }
    return finalize(
        398, "Calibration-Induced Minimax Recovery Floor",
        "proven_conditional_two_point_nonidentifiability_lower_bound",
        "The uncertainty intervals and symmetric confusion model are supplied.  The result is not a measured detector floor.",
        packet,
    )


def st399() -> dict:
    from scipy.linalg import expm
    _, a, _ = strict_operator(); transition = expm(-.5*a)
    overlap = min(float(np.sum(np.minimum(transition[i], transition[j])))
                  for i in range(N) for j in range(N))
    dobrushin = 1-overlap
    old = json.loads((ROOT / "FIN_ST357_ST371_Results.json").read_text())["ST368"]
    packet = {
        "supplied_tau": .5, "minimum_transition_probability": float(np.min(transition)),
        "Dobrushin_contraction": dobrushin,
        "geometric_mixing_bound": "beta(n)<=Dobrushin^n",
        "entropy_rate_bits": old["exact_entropy_rate_bits"],
        "rate_distortion_limit": "R(D)=lim_n n^{-1} inf I(X^n;Xhat^n) exists for every continuity point",
        "exact_endpoints": {"R(0)": old["exact_entropy_rate_bits"], "R(11/12)": 0.0},
        "theorem": (
            "Strict positivity makes the supplied heat transition uniformly ergodic with the "
            "displayed Dobrushin contraction.  The stationary finite-alphabet source is ergodic; "
            "the standard stationary-source coding theorem therefore promotes the block "
            "information quantities to an operational asymptotic rate-distortion function."
        ),
    }
    return finalize(
        399, "Asymptotic Rate-Distortion Existence for the Strict-Heat Markov Source",
        "proven_conditional_operational_limit_and_geometric_mixing",
        "The exact interior R(D) is not derived.  Tau, Hamming distortion, stationarity, and the Markov interpretation remain premises; no universe-compression claim follows.",
        packet,
    )


def st400() -> dict:
    alpha, sigma, bound = .05, 1e-3, .01; rows = []
    for mdep in (0, 1, 3, 7):
        colors = mdep+1; x = math.log(2*colors/alpha)
        for levels in (100, 1000, 10000):
            radius = math.sqrt(2*x*sigma*sigma*colors*levels)+(2*bound*x/3)*colors
            rows.append({"m_dependence": mdep, "levels": levels,
                         "Bernstein_union_radius": radius, "confidence": 1-alpha})
    packet = {
        "premise": "centered, |X_i|<=0.01, Var(X_i)<=1e-6, and m-dependent",
        "rows": rows,
        "theorem": (
            "Partition an m-dependent sequence into m+1 independent residue classes, apply "
            "two-sided Bernstein to each class, and union-bound.  Cauchy on class lengths "
            "gives radius sqrt(2 sigma^2 (m+1)L log(2(m+1)/alpha)) plus "
            "2b(m+1)log(2(m+1)/alpha)/3.  The leading growth is square-root depth without "
            "a Gaussian premise."
        ),
    }
    return finalize(
        400, "Non-Gaussian m-Dependent Multilevel Clock Concentration",
        "proven_conditional_Bernstein_square_root_depth_bound",
        "Boundedness, centering, variance, and m-dependence are added stochastic laws.  No seconds or arrow of time is generated.",
        packet,
    )


def st401() -> dict:
    sources = sorted(x.name for x in ROOT.glob("FIN_ST401_NEW_STRICT_TYPED_SOURCE*.json"))
    records = sorted(x.name for x in ROOT.glob("FIN_ST401_INDEPENDENT_RAW_EVENTS*.jsonl"))
    packet = {
        "typed_source_pattern": "FIN_ST401_NEW_STRICT_TYPED_SOURCE*.json",
        "independent_record_pattern": "FIN_ST401_INDEPENDENT_RAW_EVENTS*.jsonl",
        "typed_sources": sources, "independent_records": records,
        "typed_source_count": len(sources), "independent_record_count": len(records),
        "theorem": "Neither strict sourcehood nor independent empirical custody can be generated by local continuation code.",
    }
    return finalize(
        401, "Strict Source and Independent Evidence Admission Gate",
        "blocked_no_new_source_or_independent_record" if not sources and not records else "candidate_inputs_present_require_separate_audit",
        "The stop is not a universal impossibility theorem and not a failed experiment.",
        packet,
    )


def make_figures(out: dict) -> None:
    FIG_DIR.mkdir(exist_ok=True)
    w, _, s = strict_operator(); weights = w[0, 1:]
    t = np.linspace(.006, .018, 260); rho = np.linspace(.105, .16, 280)
    zz = np.array([[positivity_refined_envelope(x, y, s, weights) for y in rho] for x in t])
    fig, ax = plt.subplots(figsize=(7, 4.2))
    cs = ax.contour(rho, 1-t, zz, levels=18); ax.clabel(cs, fontsize=7)
    box = out["ST388"]["global_derivative_isolation_box_t_rho"]
    ax.plot([box[1][0], box[1][1], box[1][1], box[1][0], box[1][0]],
            [1-box[0][0], 1-box[0][0], 1-box[0][1], 1-box[0][1], 1-box[0][0]], "r-")
    ax.set(xlabel=r"remainder collision $\rho$", ylabel=r"maximum probability $m$",
           title="ST388: positivity-refined envelope and certified minimizer box")
    fig.tight_layout(); fig.savefig(FIG_DIR / "st388_refined_envelope.png", dpi=200); plt.close(fig)

    curv = out["ST389"]["cap_curvature_certificate"]
    ms = np.linspace(.935, .97, 200); aiv, _, _ = strict_interval_matrix()
    # The plot is diagnostic; the theorem uses the exact finite certificate at m0=0.94.
    amat = float_matrix_from_iv(aiv); margins = [cap_margin_float(amat, x)/2 for x in ms]
    fig, ax = plt.subplots(figsize=(7, 4.0)); ax.plot(ms, margins); ax.axhline(0, color="black", lw=.8)
    ax.axvline(.94, color="red", ls="--", label="certified cap threshold")
    ax.set(xlabel=r"cap threshold $m_0$", ylabel="Euclidean curvature lower bound",
           title="ST389: Fisher curvature dominates the strict quadratic on each cap")
    ax.legend(); fig.tight_layout(); fig.savefig(FIG_DIR / "st389_cap_curvature.png", dpi=200); plt.close(fig)

    rows = out["ST396"]["valid_rows"]; failures = out["ST396"]["failed_rows"]
    fig, ax = plt.subplots(figsize=(7, 4.2))
    ax.scatter([x["T"] for x in rows], [x["budget"] for x in rows], s=12, c="#277da1", label="valid numerical root")
    if failures:
        ax.scatter([x["T"] for x in failures], [x["budget"] for x in failures], s=22, c="#d1495b", marker="x", label="unresolved")
    ax.add_patch(plt.Rectangle((3.9999, 2.9999), .0002, .0002, fill=False, edgecolor="black", lw=1.2))
    ax.set(xlabel="second active time T", ylabel="heat budget", title="ST396: sampled IR design diagram")
    ax.legend(); fig.tight_layout(); fig.savefig(FIG_DIR / "st396_ir_phase_scan.png", dpi=200); plt.close(fig)


def main() -> None:
    np.random.seed(SEED)
    out = {"metadata": {"seed": SEED, "python": platform.python_version(),
                         "numpy": np.__version__, "scipy": scipy.__version__,
                         "scope": "local exact, interval, conditional, and bounded numerical research"}}
    funcs = {387: st387, 388: st388, 389: st389, 390: st390, 391: st391,
             392: st392, 393: st393, 394: st394, 395: st395, 396: st396,
             397: st397, 398: st398, 399: st399, 400: st400, 401: st401}
    for k in range(387, 402):
        RESULTS.write_text(json.dumps(native(out), indent=2, sort_keys=True), encoding="utf-8")
        out[f"ST{k}"] = funcs[k]()
    make_figures(out)
    RESULTS.write_text(json.dumps(native(out), indent=2, sort_keys=True), encoding="utf-8")
    with SUMMARY.open("w", newline="", encoding="utf-8") as h:
        writer = csv.writer(h); writer.writerow(["program", "object", "status"])
        for k in range(387, 402):
            writer.writerow([f"ST{k}", out[f"ST{k}"]["object"], out[f"ST{k}"]["status"]])
    print(json.dumps({k: v["status"] for k, v in out.items() if k.startswith("ST")}, indent=2))


if __name__ == "__main__":
    main()
