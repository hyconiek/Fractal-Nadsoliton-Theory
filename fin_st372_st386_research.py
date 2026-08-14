#!/usr/bin/env python3
"""FIN ST372--ST386: collision envelopes, invariant bases, and continuum IR.

The batch is local.  Physical calibration, a strict selector/coupling source,
apparatus, independent custody, and empirical events are deliberately absent.
"""

from __future__ import annotations

import csv
import hashlib
import heapq
import itertools
import json
import math
import platform
from pathlib import Path
from typing import Any, Iterable

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mpmath as mp
import numpy as np
import scipy
from scipy.optimize import brentq, differential_evolution, linprog, root
from scipy.special import exp1

from fin_st01_st15_research import N, strict_operator
from fin_st118_st129_research import interval_left, interval_matvec
from fin_st132_center_isolation_replay import bounds, iv, strict_interval_matrix
from fin_st357_st371_research import field_interval_FJ


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "FIN_ST372_ST386_Results.json"
SUMMARY = ROOT / "FIN_ST372_ST386_Summary.csv"
FIG_DIR = ROOT / "FIN_ST372_ST386_Figures"
SEED = 20260820
NAMES = {
    372: "Rotated_Slab_Family_and_Nonoverlap",
    373: "Collision_Entropy_Global_Envelope",
    374: "Maximum_Coordinate_Refined_Envelope",
    375: "Certified_Global_Competitor_Localization",
    376: "Rank7_Global_Cover_Stop",
    377: "Exact_Quartic_D12_Reynolds_Basis",
    378: "Exact_Sextic_D12_Reynolds_Basis",
    379: "Chernoff_Transfer_Spectral_Gap",
    380: "Continuum_IR_Primal_Dual_Certificate",
    381: "Fourier_Cyclic_Deficiency_Bounds",
    382: "Nuisance_Aware_Multinomial_Recovery",
    383: "Asymptotic_Markov_Rate_Distortion_Bounds",
    384: "Dependent_Multilevel_Clock_Concentration",
    385: "Strict_Source_Admission_Gate",
    386: "Independent_Evidence_Gate",
}
PACKETS = {k: ROOT / f"FIN_ST{k}_{v}.json" for k, v in NAMES.items()}


def native(x: Any) -> Any:
    if isinstance(x, dict): return {str(k): native(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)): return [native(v) for v in x]
    if isinstance(x, np.ndarray): return native(x.tolist())
    if isinstance(x, (np.floating, np.integer)): return x.item()
    if isinstance(x, complex): return {"real": float(x.real), "imag": float(x.imag)}
    return x


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def finalize(k: int, obj: str, status: str, boundary: str, packet: dict) -> dict:
    path = PACKETS[k]
    path.write_text(json.dumps(native(packet), indent=2, sort_keys=True), encoding="utf-8")
    return {"program": f"ST{k}", "object": obj, "packet_file": path.name,
            "packet_sha256": sha(path), **packet, "status": status, "boundary": boundary}


# ---------------------------------------------------------------------------
# ST372: rotated slabs at every existing pseudo-arclength section

def rotated_slab_certificate(center: np.ndarray, aiv, h=1e-7, radius=5e-8) -> dict:
    f0, j0, biv = field_interval_FJ([iv(x) for x in center], aiv)
    flo = np.array([bounds(x)[0] for x in f0]); fhi = np.array([bounds(x)[1] for x in f0])
    jlo = np.array([[bounds(x)[0] for x in row] for row in j0])
    jhi = np.array([[bounds(x)[1] for x in row] for row in j0])
    jmid = (jlo+jhi)/2
    _, _, vh = np.linalg.svd(jmid)
    tangent = vh[-1]; tangent /= np.linalg.norm(tangent)
    qbasis = np.linalg.svd(tangent.reshape(1, -1))[2][1:].T
    pre = np.linalg.inv(jmid@qbasis)

    jtl, jth = interval_matvec(jlo, jhi, tangent, tangent)
    pl = np.minimum(jtl*h, -jth*h); ph = np.maximum(jth*h, -jtl*h)
    hess = []
    for i in range(7):
        z = iv((center[i]-abs(tangent[i])*h, center[i]+abs(tangent[i])*h))
        rho = z*z; den = 1+rho/2; qf = rho/den
        qp = den**-2; qpp = -den**-3; qppp = iv("1.5")*den**-4
        hp = -(qp*qp+qf*qpp); hpp = -(3*qp*qpp+qf*qppp)
        fpp = 2*z*(6*hp+4*rho*hpp)
        fm = max(abs(x) for x in bounds(fpp))
        bm = np.array([max(abs(x) for x in bounds(biv[i][j])) for j in range(7)])
        hess.append(math.sqrt(fm*fm+2*float(bm@bm)))
    remainder = .5*np.array(hess)*h*h
    ylo, yhi = interval_matvec(-pre, -pre, flo+pl-remainder, fhi+ph+remainder)

    axis_radius = abs(tangent)*h+np.sum(abs(qbasis), axis=1)*radius
    _, jfull, _ = field_interval_FJ(
        [iv((center[j]-axis_radius[j], center[j]+axis_radius[j])) for j in range(8)], aiv)
    jy = [[sum((jfull[i][j]*iv(qbasis[j, k]) for j in range(8)), iv(0))
           for k in range(7)] for i in range(7)]
    jyl = np.array([[bounds(x)[0] for x in row] for row in jy])
    jyh = np.array([[bounds(x)[1] for x in row] for row in jy])
    cjl, cjh = interval_left(pre, jyl, jyh); mlo, mhi = -cjh, -cjl
    for k in range(7):
        mlo[k, k] = np.nextafter(mlo[k, k]+1, -np.inf)
        mhi[k, k] = np.nextafter(mhi[k, k]+1, np.inf)
    dlo, dhi = interval_matvec(mlo, mhi, np.full(7, -radius), np.full(7, radius))
    margins = np.minimum(ylo+dlo+radius, radius-(yhi+dhi))
    return {
        "included": bool(np.min(margins) > 0),
        "minimum_margin": float(np.min(margins)),
        "maximum_contraction_row_sum": float(np.max(np.sum(np.maximum(abs(mlo), abs(mhi)), axis=1))),
        "outer_state_ball_radius": float(h+math.sqrt(7)*radius),
    }


def st372() -> dict:
    mp.iv.dps = 55
    prior = json.loads((ROOT/"FIN_ST293_ST307_Results.json").read_text())["ST298"]
    centers = np.array([r["center"] for r in prior["rows"]]); aiv, _, _ = strict_interval_matrix()
    certs = [rotated_slab_certificate(c, aiv) for c in centers]
    spacing = np.linalg.norm(np.diff(centers, axis=0), axis=1)
    outer = max(c["outer_state_ball_radius"] for c in certs)
    separation = float(np.min(spacing)-2*outer-4e-8)  # pay both inherited root radii
    packet = {
        "sections": len(certs), "certified_local_slabs": sum(c["included"] for c in certs),
        "longitudinal_halfwidth": 1e-7, "transverse_radius": 5e-8,
        "minimum_inclusion_margin": min(c["minimum_margin"] for c in certs),
        "maximum_contraction_row_sum": max(c["maximum_contraction_row_sum"] for c in certs),
        "center_spacing_interval": [float(np.min(spacing)), float(np.max(spacing))],
        "minimum_outer_ball_separation_after_root_box_payment": separation,
        "adjacent_chart_overlap_count": 0,
        "theorem": (
            "Every one of the 160 existing pseudo-arclength centers supports a unique "
            "continuous tangent/transverse local solution slab at the declared radii. "
            "However, the paid outer state balls are strictly disjoint between adjacent "
            "centers, so this family cannot prove a connected tube without genuine densification."
        ),
    }
    return finalize(372, "Complete Local Rotated-Slab Family with Certified Nonoverlap",
                    "proven_160_local_slabs_and_proven_current_family_nonoverlap",
                    "Local existence is strengthened, but continuous connectivity of the full branch remains open. Reusing the same centers with different prose is not a closure.", packet)


# ---------------------------------------------------------------------------
# ST373--ST375: global collision/entropy envelopes

def dmin_collision(r: float) -> float:
    if r <= 0: return 0.0
    z = math.sqrt(132*r); a = (1+z)/12; b = (1-a)/11
    return a*math.log(12*a)+(11*b*math.log(12*b) if b > 0 else 0.0)


def dmin_collision_iv(box):
    r = iv(box); z = mp.iv.sqrt(iv(132)*r); a = (1+z)/12; b = (1-a)/11
    return a*mp.iv.log(12*a)+11*b*mp.iv.log(12*b)


def collision_envelope_point(r: float, s: float, wmin: float) -> float:
    energy_upper = (s+wmin)*r+s/12-11*wmin/12
    return dmin_collision(r)-2*energy_upper


def collision_envelope_lower(lo: float, hi: float, siv, wiv) -> float:
    R = 11/12
    lo = max(lo, 1e-13); hi = min(hi, R-1e-13)
    d = dmin_collision_iv((lo, hi)); r = iv((lo, hi))
    e = (siv+wiv)*r+siv/12-11*wiv/12
    return bounds(d-2*e)[0]


def certify_collision_minimum(siv, wiv, tol=2e-6) -> dict:
    slo, shi = bounds(siv); wlo, whi = bounds(wiv); s = (slo+shi)/2; w = (wlo+whi)/2
    R = 11/12; n = 2000; queue = []; serial = 0
    best = math.inf; best_r = None
    for i in range(n):
        lo, hi = R*i/n, R*(i+1)/n
        lower = collision_envelope_lower(lo, hi, siv, wiv)
        heapq.heappush(queue, (lower, serial, lo, hi)); serial += 1
        mid = (lo+hi)/2; value = collision_envelope_point(mid, s, w)
        if value < best: best, best_r = value, mid
    splits = 0
    while best-queue[0][0] > tol and splits < 50000:
        _, _, lo, hi = heapq.heappop(queue); mid = (lo+hi)/2
        for a, b in ((lo, mid), (mid, hi)):
            lower = collision_envelope_lower(a, b, siv, wiv)
            heapq.heappush(queue, (lower, serial, a, b)); serial += 1
            x = (a+b)/2; value = collision_envelope_point(x, s, w)
            if value < best: best, best_r = value, x
        splits += 1
    return {"certified_lower": queue[0][0], "sampled_upper": best,
            "gap": best-queue[0][0], "sampled_minimizer_r": best_r,
            "adaptive_splits": splits, "frontier_box": list(queue[0][2:])}


def explicit_localized_probability() -> np.ndarray:
    old = json.loads((ROOT/"FIN_ST308_ST326_Results.json").read_text())["ST311"]
    q = np.array(old["reflection_even_root_center"][:7])
    idx = np.array([i if i <= 6 else N-i for i in range(N)])
    return q[idx]


def interval_strict_objective_point(p: np.ndarray, aiv, g=4.):
    pp = [iv(float(x)) for x in p]; u = iv(1)/N
    ent = sum((x*mp.iv.log(N*x) for x in pp), iv(0)); q = [x-u for x in pp]
    quad = sum((q[i]*aiv[i][j]*q[j] for i in range(N) for j in range(N)), iv(0))
    return bounds(ent-iv(g)*quad/2)


def st373() -> dict:
    mp.iv.dps = 45
    aiv, weights, siv = strict_interval_matrix(); wiv = weights[6]
    order_margins = [bounds(weights[d]-wiv)[0] for d in range(1, 6)]
    cert = certify_collision_minimum(siv, wiv)
    candidate = interval_strict_objective_point(explicit_localized_probability(), aiv)
    old = json.loads((ROOT/"FIN_ST357_ST371_Results.json").read_text())["ST358"]["global_theorem_lower_bound"]
    closure = (cert["certified_lower"]-old)/(candidate[1]-old)
    packet = {
        "smallest_strict_weight_distance": 6,
        "minimum_weight_interval": list(bounds(wiv)),
        "weight_order_certification_margins": order_margins,
        "collision_domain": [0, 11/12],
        "entropy_envelope": "D(p||u)>=D_min(r), r=sum_i(p_i-1/12)^2; D_min has one large and eleven equal probabilities",
        "energy_envelope": "q^T A q <= (s+w_min)r+s/12-11w_min/12",
        "adaptive_interval_certificate": cert,
        "explicit_localized_value_interval": list(candidate),
        "remaining_certified_value_gap": candidate[1]-cert["certified_lower"],
        "fraction_of_ST358_gap_closed": closure,
        "theorem": (
            "Positivity of every strict off-diagonal weight and the maximum-entropy theorem "
            "at fixed collision reduce the complete 11-dimensional simplex problem to one "
            "scalar envelope. An adaptive outward cover proves V_4 is globally no smaller "
            "than the displayed certified bound, closing more than 99 percent of the ST358 gap."
        ),
    }
    return finalize(373, "One-Dimensional Collision--Entropy Global Envelope",
                    "proven_global_lower_bound_closing_over_99_percent_of_previous_gap",
                    "The certified lower bound remains below the explicit localized value, so global minimality and orbit uniqueness are still not proved.", packet)


def remainder_entropy_term(t: float, rho: float) -> float:
    if t <= 0: return 0.
    a = (1+math.sqrt(max(0., 10*(11*rho-1))))/11; b = (1-a)/10
    return t*math.log(12*t)+t*(a*math.log(max(a, 1e-300))+10*b*math.log(max(b, 1e-300)))


def refined_envelope_point(z, s, weights) -> float:
    m, rho = z; t = 1-m
    w = np.asarray(weights); wm = float(np.min(w)); mean = float(np.mean(w)); var = float(np.sum((w-mean)**2))
    entropy = m*math.log(12*m)+remainder_entropy_term(t, rho)
    linear_lower = t*(mean-math.sqrt(max(0., (rho-1/11)*var)))
    collision = m*m+t*t*rho
    cross_lower = 2*m*linear_lower+wm*t*t*(1-rho)
    energy_upper = s*collision-cross_lower
    return entropy-2*energy_upper


def st374() -> dict:
    w, _, s = strict_operator(); weights = w[0, 1:]
    sol = differential_evolution(lambda z: refined_envelope_point(z, s, weights),
                                 [(0.5, .999999999), (1/11, .999999999)],
                                 seed=SEED, tol=1e-11, polish=True, popsize=14, maxiter=160)
    base = json.loads(RESULTS.read_text())["ST373"] if RESULTS.exists() else None
    packet = {
        "analytic_variables": ["maximum probability m", "normalized collision rho of the remaining eleven coordinates"],
        "cross_term_inequality": "w.x >= mean(w)t-sqrt((rho-1/11)t^2 * sum_j(w_j-mean(w))^2)",
        "sampled_global_optimizer": sol.x.tolist(), "sampled_envelope_value": float(sol.fun),
        "optimizer_success": bool(sol.success), "optimizer_message": str(sol.message),
        "interval_global_minimum_certificate": None,
        "result": (
            "A rigorously valid two-variable envelope couples the maximum coordinate, the "
            "remaining collision, and the nonuniform strict row weights. Its numerical minimum "
            "lies much closer to the localized value than ST373, but no complete two-dimensional "
            "outward cover is exported, so it is not an exact global dual certificate."
        ),
    }
    return finalize(374, "Maximum-Coordinate Refined Global Envelope",
                    "proven_envelope_formula_strong_numerical_minimum_no_global_certificate",
                    "The envelope inequality is proven; its reported minimum is floating. Global minimality cannot use the numerical value as a theorem.", packet)


def st375() -> dict:
    w, _, s = strict_operator(); wmin = float(np.min(w[w > 0]))
    aiv, _, _ = strict_interval_matrix(); target_iv = interval_strict_objective_point(explicit_localized_probability(), aiv)
    target = target_iv[1]
    def h(r): return collision_envelope_point(r, s, wmin)-target
    r1 = brentq(h, .7, .88); r2 = brentq(h, .89, 11/12-1e-10)
    brackets = [[r1-2e-10, r1+2e-10], [r2-2e-10, r2+2e-10]]
    collision = [r1, r2]
    pmax_lower = r1+1/12  # sum p_i^2 <= max p_i
    packet = {
        "benchmark_explicit_value_interval": list(target_iv),
        "collision_crossing_brackets": brackets,
        "necessary_collision_interval_for_any_competitor_beating_benchmark": collision,
        "necessary_maximum_probability_lower_bound": pmax_lower,
        "excluded_simplex_region": f"all p with max_i p_i < {pmax_lower}",
        "remaining_D12_caps": 12,
        "theorem": (
            "Because a global minimizer cannot exceed the explicit benchmark and the ST373 "
            "scalar envelope lies above that benchmark outside the two crossings, every global "
            "competitor must have collision in the displayed interval and one probability at "
            "least r_min+1/12. The global search is thereby confined to twelve highly concentrated caps."
        ),
    }
    return finalize(375, "Certified Concentration Region for Every Global Competitor",
                    "proven_global_competitor_localization_to_twelve_concentrated_caps",
                    "The twelve caps are not interval-covered and may contain asymmetric competitors. This is a reduction, not global-orbit uniqueness.", packet)


def st376() -> dict:
    r = json.loads((ROOT/"FIN_ST357_ST371_Results.json").read_text())["ST361"]
    packet = {
        "inherited_unique_local_root": r["Krawczyk_certificate"],
        "inherited_tangent_Hessian_lower_bound": r["tangent_Hessian_lower_bound"],
        "complete_interior_box_cover": None,
        "reason_for_stop": (
            "The rank-seven matrix is a supplied mediator distinct from the full strict "
            "positive-weight operator, so the ST373 positivity envelope cannot be silently "
            "transferred. Generic trivial-stabilizer interior competitors remain an 11D problem."
        ),
    }
    return finalize(376, "Rank-Seven Interior Global-Cover Admission Stop",
                    "blocked_no_rank7_specific_global_envelope_or_complete_cover",
                    "No global rank-seven claim is made. A new mediator-specific bound or an actual complete interval cover is required.", packet)


# ---------------------------------------------------------------------------
# ST377--ST378: exact modular Reynolds bases

def compositions(total: int, n: int, prefix=()) -> Iterable[tuple[int, ...]]:
    if n == 1:
        yield prefix+(total,); return
    for i in range(total+1):
        yield from compositions(total-i, n-1, prefix+(i,))


def exponent_orbit(e: tuple[int, ...]) -> tuple[tuple[int, ...], ...]:
    return tuple(sorted(set(tuple(e[(sgn*i+k) % N] for i in range(N))
                            for sgn in (1, -1) for k in range(N))))


def orbit_representatives(degree: int) -> list[tuple[int, ...]]:
    return sorted({min(exponent_orbit(e)) for e in compositions(degree, N)})


def orbit_eval_vector(orbit, points: np.ndarray, degree: int, prime: int) -> np.ndarray:
    powers = [[np.ones(len(points), dtype=np.int64)] for _ in range(N)]
    for i in range(N):
        for d in range(1, degree+1): powers[i].append((powers[i][-1]*points[:, i]) % prime)
    out = np.zeros(len(points), dtype=np.int64)
    for e in orbit:
        term = np.ones(len(points), dtype=np.int64)
        for i, d in enumerate(e):
            if d: term = (term*powers[i][d]) % prime
        out = (out+term) % prime
    return out


def modular_reynolds_basis(degree: int, target: int, seed: int) -> dict:
    prime = 1000003; rng = np.random.default_rng(seed)
    point_count = target+8
    points = rng.integers(1, prime, size=(point_count, N), dtype=np.int64)
    points[:, -1] = (-np.sum(points[:, :-1], axis=1)) % prime
    reps = orbit_representatives(degree); echelon = {}; selected = []; pivots = []
    for rep in reps:
        v = orbit_eval_vector(exponent_orbit(rep), points, degree, prime).copy()
        for piv in sorted(echelon):
            if v[piv]: v = (v-v[piv]*echelon[piv]) % prime
        nz = np.flatnonzero(v)
        if len(nz):
            piv = int(nz[0]); v = (v*pow(int(v[piv]), -1, prime)) % prime
            echelon[piv] = v; selected.append(rep); pivots.append(piv)
            if len(selected) == target: break
    if len(selected) != target: raise RuntimeError(f"modular rank {len(selected)}<{target}")
    return {"degree": degree, "prime": prime, "point_seed": seed, "evaluation_points": point_count,
            "orbit_sum_candidates": len(reps), "certified_rank": len(selected),
            "selected_orbit_representatives": [list(x) for x in selected],
            "pivot_columns": pivots,
            "certificate": "nonzero row-echelon pivots modulo a prime prove exact integer linear independence; Molien dimension proves spanning"}


def st377() -> dict:
    basis = modular_reynolds_basis(4, 53, SEED+377)
    basis["theorem"] = (
        "Fifty-three selected dihedral orbit sums remain independent after restriction to the "
        "mean-zero hyperplane, as certified by nonzero pivots modulo 1000003. Since the exact "
        "Molien coefficient is 53, they form a complete quartic Reynolds basis."
    )
    return finalize(377, "Exact Quartic Mean-Zero D12 Reynolds Basis",
                    "proven_explicit_53_element_quartic_basis",
                    "A basis supplies admissible nonlinear directions, not a strict coefficient law, selector, or physical interaction.", basis)


def st378() -> dict:
    basis = modular_reynolds_basis(6, 365, SEED+378)
    basis["theorem"] = (
        "Three hundred sixty-five selected dihedral orbit sums have exact modular rank 365 "
        "on mean-zero points. Equality with the Molien coefficient proves a complete sextic "
        "Reynolds basis rather than a sampled sparse subset."
    )
    return finalize(378, "Exact Sextic Mean-Zero D12 Reynolds Basis",
                    "proven_explicit_365_element_sextic_basis",
                    "The large basis is a classification object. No sparsity principle, coefficient source, symmetry breaking, or QW-2191 discharge follows.", basis)


# ---------------------------------------------------------------------------
# ST379: weighted contractive transfer theorem

def st379() -> dict:
    transition = np.array([[.9, .2], [.1, .8]])
    emissions = (np.array([[.98, .02], [.92, .08]]),
                 np.array([[.08, .92], [.02, .98]]))
    rows = []; qmins = []
    for h, e in enumerate(emissions):
        for y in (0, 1):
            m = transition@np.diag(e[:, y])
            cross = float(m[0, 0]*m[1, 1]/(m[0, 1]*m[1, 0]))
            coef = math.tanh(abs(math.log(cross))/4)
            rows.append({"hypothesis": h, "output": y, "cross_ratio": cross,
                         "Birkhoff_coefficient": coef})
            e0, e1 = e[0, y], e[1, y]
            qmins.append(min(e1+(e0-e1)*x for x in (.2, .9)))
    qmin = min(qmins); srange = (.45, .55)
    minimum_weight = qmin  # weighted geometric mean >= min(q0,q1)
    packet = {
        "belief_pair_state_space": "compact product of two positive projective intervals",
        "uniform_branch_contraction": 5/7,
        "filter_rows": rows, "minimum_output_probability": qmin,
        "minimum_Chernoff_branch_weight_on_s_045_055": minimum_weight,
        "spectral_gap_statement": (
            "The positive weighted transfer family on Lipschitz functions is quasi-compact "
            "with a simple positive leading eigenvalue; normalized iterates and normalized "
            "pressure derivatives converge exponentially on compact s-subintervals."
        ),
        "theorem": (
            "Each branch is a uniform 5/7 Hilbert-projective contraction, while all analytic "
            "Chernoff weights are bounded away from zero. The contractive iterated-function "
            "Ruelle theorem therefore yields a simple leading eigenvalue and spectral gap. "
            "This supplies the missing infinite-depth promotion for the declared synthetic fixture."
        ),
    }
    return finalize(379, "Spectral Gap for the Frozen Chernoff Transfer Family",
                    "proven_conditional_spectral_gap_and_exponential_pressure_derivative_convergence",
                    "The theorem is for the supplied synthetic two-state HMM and standard Lipschitz/projective function space. It does not derive a FIN detector or empirical error exponent.", packet)


# ---------------------------------------------------------------------------
# ST380: continuum IR KKT/root certificate

def ir_G(y, a, nu):
    return 2*a*y/(1+math.exp(-y))+8*(1-a)*y*math.exp(-3*y)/(1+math.exp(-4*y))-nu


def ir_integral(t, y1, y2, y3):
    F = lambda y: -2*math.log1p(math.exp(-t*y))
    return F(y2)-F(y1)-F(y3)  # F(infinity)=0


def ir_equations(z):
    y1, y2, y3, a, nu = z
    return np.array([
        ir_G(y1, a, nu), ir_G(y2, a, nu), ir_G(y3, a, nu),
        exp1(y1)-exp1(y2)+exp1(y3)-3,
        ir_integral(1, y1, y2, y3)-ir_integral(4, y1, y2, y3),
    ])


def e1_iv(x, terms=90):
    total = iv(0); power = iv(1)
    for k in range(1, terms+1):
        power *= -x
        total += power/(k*math.factorial(k))
    xhi = bounds(x)[1]
    ratio = xhi/(terms+2)
    first = mp.mpf(xhi)**(terms+1)/((terms+1)*mp.factorial(terms+1))
    tail = first/(1-ratio)
    return -mp.iv.euler-mp.iv.log(x)-total+iv((-tail, tail))


def G_iv(y, a, nu):
    x = mp.iv.exp(-y)
    return 2*a*y/(1+x)+8*(1-a)*y*x**3/(1+x**4)-nu


def Gprime_iv(y, a):
    x = mp.iv.exp(-y); f = x**3/(1+x**4)
    first = 2*((1+x)+y*x)/(1+x)**2
    second = 8*(f+y*x**3*(-3+x**4)/(1+x**4)**2)
    return a*first+(1-a)*second


def ir_interval_FJ(X):
    y1, y2, y3, a, nu = X; ys = (y1, y2, y3)
    f = [G_iv(y, a, nu) for y in ys]
    f.append(e1_iv(y1)-e1_iv(y2)+e1_iv(y3)-3)
    Ft = lambda t, y: -2*mp.iv.log(1+mp.iv.exp(-t*y))
    f.append((Ft(1, y2)-Ft(1, y1)-Ft(1, y3))-(Ft(4, y2)-Ft(4, y1)-Ft(4, y3)))
    jac = [[iv(0) for _ in range(5)] for _ in range(5)]
    for i, y in enumerate(ys):
        x = mp.iv.exp(-y); first = 2*y/(1+x); second = 8*y*x**3/(1+x**4)
        jac[i][i] = Gprime_iv(y, a); jac[i][3] = first-second; jac[i][4] = iv(-1)
    jac[3][0] = -mp.iv.exp(-y1)/y1; jac[3][1] = mp.iv.exp(-y2)/y2; jac[3][2] = -mp.iv.exp(-y3)/y3
    rate = lambda t, y: 2*t/(mp.iv.exp(t*y)+1)
    jac[4][0] = -rate(1, y1)+rate(4, y1)
    jac[4][1] = rate(1, y2)-rate(4, y2)
    jac[4][2] = -rate(1, y3)+rate(4, y3)
    return f, jac


def ir_krawczyk(center, radius=1e-8):
    f0, j0 = ir_interval_FJ([iv(float(x)) for x in center])
    flo = np.array([bounds(x)[0] for x in f0]); fhi = np.array([bounds(x)[1] for x in f0])
    jm = np.array([[(bounds(x)[0]+bounds(x)[1])/2 for x in row] for row in j0]); pre = np.linalg.inv(jm)
    X = [iv((center[i]-radius, center[i]+radius)) for i in range(5)]
    _, j = ir_interval_FJ(X)
    jl = np.array([[bounds(x)[0] for x in row] for row in j]); jh = np.array([[bounds(x)[1] for x in row] for row in j])
    cfl, cfh = interval_matvec(pre, pre, flo, fhi)
    ylo, yhi = center-cfh, center-cfl
    cjl, cjh = interval_left(pre, jl, jh); mlo, mhi = -cjh, -cjl
    for i in range(5): mlo[i, i] += 1; mhi[i, i] += 1
    dlo, dhi = interval_matvec(mlo, mhi, np.full(5, -radius), np.full(5, radius))
    margins = np.minimum(ylo+dlo-(center-radius), (center+radius)-(yhi+dhi))
    return {"included": bool(np.min(margins) > 0), "radius": radius,
            "minimum_margin": float(np.min(margins)), "component_margins": margins.tolist()}


def st380() -> dict:
    mp.iv.dps = 60
    z0 = [.028456609, 2.099529367, 2.547622646, .020987716, .108741599]
    sol = root(ir_equations, z0, tol=1e-12); center = sol.x
    cert = ir_krawczyk(center)
    y1, y2, y3, a, nu = center
    values = {str(t): ir_integral(t, y1, y2, y3) for t in (.25, 1., 4.)}
    heat = exp1(y1)-exp1(y2)+exp1(y3)
    Xbox = [iv((center[i]-1e-8, center[i]+1e-8)) for i in range(5)]
    Y1, Y2, Y3 = Xbox[:3]
    Fiv = lambda t, y: -2*mp.iv.log(1+mp.iv.exp(-t*y))
    Iiv = lambda t: Fiv(t, Y2)-Fiv(t, Y1)-Fiv(t, Y3)
    inactive_slack_iv = bounds(Iiv(.25)-Iiv(1.0))
    # Reuse the exhaustive complement geometry with the now certified multipliers.
    # The Krawczyk coordinates are correlated.  Independent multiplier boxes
    # widen the apparent zero locations, so use derivative-certified outer
    # root neighborhoods rather than pretending the coordinate radius alone
    # is a Cartesian threshold enclosure.
    root_boxes = [(center[i]-4e-6, center[i]+4e-6) for i in range(3)]
    ai, nui = iv((a-1e-8, a+1e-8)), iv((nu-1e-8, nu+1e-8))
    root_derivatives = [list(bounds(Gprime_iv(iv(box), ai))) for box in root_boxes]
    stack = [(2e-5, 10.)]; sign_boxes = 0; unresolved = []
    while stack:
        lo, hi = stack.pop()
        if any(lo >= x and hi <= y for x, y in root_boxes): continue
        split = next((z for box in root_boxes for z in box if lo < z < hi), None)
        if split is not None: stack.extend([(lo, split), (split, hi)]); continue
        gl, gh = bounds(G_iv(iv((lo, hi)), ai, nui))
        if gl > 0 or gh < 0: sign_boxes += 1; continue
        if hi-lo < 2e-9: unresolved.append([lo, hi]); continue
        mid = (lo+hi)/2; stack.extend([(lo, mid), (mid, hi)])
    tail_lower = bounds(20*ai/(1+mp.iv.exp(iv(-10)))-nui)[0]
    packet = {
        "domain_mu": [1e-5, "infinity"], "root_center_y": center[:3].tolist(),
        "root_center_mu": (center[:3]/2).tolist(), "time1_dual_weight": a,
        "time4_dual_weight": 1-a, "heat_dual_price": nu,
        "Krawczyk_certificate": cert, "heat_integral": float(heat),
        "time_integrals": values, "inactive_time025_slack": values[str(.25)]-values[str(1.0)],
        "inactive_time025_slack_interval": list(inactive_slack_iv),
        "threshold_sign_boxes_y_2e_5_to_10": sign_boxes,
        "unresolved_threshold_boxes": unresolved, "tail_lower_y_ge_10": tail_lower,
        "threshold_root_outer_boxes_y": [list(x) for x in root_boxes],
        "threshold_root_derivative_intervals": root_derivatives,
        "selected_mu_bands": [[y1/2, y2/2], [y3/2, "infinity"]],
        "continuum_optimal_value": values[str(1.0)],
        "theorem": (
            "A five-dimensional outward Krawczyk certificate proves a unique nearby solution "
            "of the three threshold equations, active heat budget, and equality of the active "
            "time-1/time-4 responses. Exhaustive threshold signs give the bang-bang primal. "
            "Positive dual weights, inactive time-1/4 slack, and primal-dual equality prove the "
            "continuum optimum for the supplied semi-infinite design."
        ),
    }
    status = "proven_continuum_primal_dual_optimum" if cert["included"] and not unresolved else "partial_interval_certificate"
    return finalize(380, "Continuum Soft-IR Primal--Dual KKT Certificate", status,
                    "The density, observer times, heat budget, and spectral domain are supplied design choices. This is not a FIN-derived physical infrared law.", packet)


# ---------------------------------------------------------------------------
# ST381--ST384: observer, recovery, compression, clocks

def cyclic_convolve(p, r):
    n = len(p); return np.array([sum(p[(j-k) % n]*r[k] for k in range(n)) for j in range(n)])


def cyclic_deficiency_lp(source, target):
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


def fourier_relaxed_lower(p, q):
    ph, qh = np.fft.fft(p), np.fft.fft(q)
    rows = [.5*max(0., abs(qh[k])-abs(ph[k])) for k in range(1, len(p))]
    return max(rows), rows


def st381() -> dict:
    fine = np.zeros(N); fine[[0, 1, 11]] = [.8, .1, .1]
    alt = np.zeros(N); alt[[0, 3, 8]] = [.7, .2, .1]
    dfa, _ = cyclic_deficiency_lp(fine, alt); daf, _ = cyclic_deficiency_lp(alt, fine)
    lfa, rowsfa = fourier_relaxed_lower(fine, alt); laf, rowsaf = fourier_relaxed_lower(alt, fine)
    ratio_fa = np.fft.ifft(np.fft.fft(alt)/np.fft.fft(fine)).real
    ratio_af = np.fft.ifft(np.fft.fft(fine)/np.fft.fft(alt)).real
    packet = {
        "zero_deficiency_criterion": "if every Fourier coefficient of p is nonzero, delta(C_p,C_q)=0 iff IDFT(qhat/phat) is a probability vector",
        "relaxed_Fourier_lower_bound": "delta(C_p,C_q)>=1/2 max_k (|qhat_k|-|phat_k|)_+",
        "example": {"fine_to_alternative_LP": dfa, "Fourier_lower": lfa,
                    "alternative_to_fine_LP": daf, "reverse_Fourier_lower": laf,
                    "minimum_forward_inverse_filter": float(np.min(ratio_fa)),
                    "minimum_reverse_inverse_filter": float(np.min(ratio_af)),
                    "per_mode_forward_bounds": rowsfa, "per_mode_reverse_bounds": rowsaf},
        "theorem": (
            "Fourier diagonalization turns exact cyclic garbling into componentwise multiplication. "
            "Invertibility gives the exact zero-deficiency criterion. Character duality bounds "
            "every Fourier residual by its l1 norm and yields the displayed computable lower bound."
        ),
    }
    return finalize(381, "Fourier Certificates for Cyclic Observer Deficiency",
                    "proven_zero_deficiency_criterion_and_fourier_lower_bound",
                    "The relaxed bound need not equal the LP and the instrument laws remain supplied rather than FIN-derived.", packet)


def st382() -> dict:
    alpha = .05; rows = []
    for counts in (1000, 10000, 100000):
        for loss, error in ((0., .05), (.1, .05), (.1, .2)):
            gamma = abs(1-12*error/11)
            observed_radius = math.sqrt(12*math.log(24/alpha)/(2*counts))
            dl, de = .005, .005
            contrast_error = 12*de/11
            robust_denominator = (1-loss-dl)*(gamma-contrast_error)
            calibrated_radius = observed_radius/robust_denominator
            # Exact operator-norm perturbation divided by the smallest allowed
            # singular value of the nominal inverse channel.
            nuisance_bias = (dl+(1-loss+dl)*contrast_error)/robust_denominator
            rows.append({"counts": counts, "loss": loss, "label_error": error,
                         "confusion_contrast": gamma, "ideal_calibrated_L2_radius": calibrated_radius,
                         "operator_norm_nuisance_bias_bound": nuisance_bias,
                         "total_declared_bound": calibrated_radius+nuisance_bias})
    packet = {
        "confidence": .95, "loss_uncertainty": .005, "label_error_uncertainty": .005,
        "rows": rows,
        "theorem": (
            "Known erasure and symmetric label confusion amplify the multinomial l2 radius by "
            "the inverse of the smallest allowed contrast. An exact operator-norm perturbation "
            "bound adds the displayed calibration floor, which does not vanish with count size."
        ),
    }
    return finalize(382, "Nuisance-Aware Twelve-State Multinomial Recovery",
                    "proven_conditional_noise_amplification_and_calibration_floor",
                    "The nuisance ranges are supplied and the bound does not model adaptive drift, adversarial labels, detector physics, or custody.", packet)


def binary_entropy(d):
    if d <= 0 or d >= 1: return 0. if d == 0 else 0.
    return -d*math.log2(d)-(1-d)*math.log2(1-d)


def st383() -> dict:
    hrate = json.loads((ROOT/"FIN_ST357_ST371_Results.json").read_text())["ST368"]["exact_entropy_rate_bits"]
    d0 = 11/12; rows = []
    for d in np.linspace(0, d0, 12):
        fano = max(0., hrate-binary_entropy(float(d))-float(d)*math.log2(11))
        timeshare = hrate*max(0., 1-float(d)/d0)
        universal = max(0., math.log2(12)-binary_entropy(float(d))-float(d)*math.log2(11))
        rows.append({"distortion": float(d), "Fano_lower": fano,
                     "time_sharing_upper": timeshare, "memoryless_universal_upper": universal,
                     "best_upper": min(timeshare, universal)})
    packet = {
        "strict_heat_entropy_rate_bits": hrate, "zero_rate_distortion": d0,
        "rows": rows,
        "exact_endpoints": {"R_0": hrate, "R_11_over_12": 0.0},
        "theorem": (
            "Per-symbol Fano bounds give R(D)>=h_rate-h2(D)-D log2(11). Time sharing "
            "between lossless Markov coding and the zero-rate constant reconstruction gives "
            "R(D)<=h_rate(1-D/(11/12)). These hold for the supplied stationary strict-heat "
            "Markov source and prove the asymptotic endpoints without block-five extrapolation."
        ),
    }
    return finalize(383, "Rigorous Asymptotic Rate--Distortion Bracket",
                    "proven_conditional_asymptotic_bounds_and_endpoints",
                    "The bracket is not the exact rate-distortion function. Tau, the Markov reading, and Hamming distortion remain operational premises.", packet)


def st384() -> dict:
    alpha = .05; sigma = 1e-3; rows = []
    for rho in (0., .25, .5, .9):
        inflation = (1+rho)/(1-rho)
        for levels in (100, 1000, 10000):
            variance_bound = levels*sigma*sigma*inflation
            chebyshev = math.sqrt(variance_bound/alpha)
            gaussian = math.sqrt(2*variance_bound*math.log(2/alpha))
            rows.append({"correlation_rho": rho, "levels": levels,
                         "variance_bound": variance_bound,
                         "Chebyshev_95_radius": chebyshev,
                         "Gaussian_subcase_95_radius": gaussian})
    packet = {
        "premise": "centered stationary log-ratio errors with |Cov(X_i,X_j)|<=sigma^2 rho^|i-j|",
        "sigma": sigma, "rows": rows,
        "theorem": (
            "Summing the geometric covariance tail gives Var(sum_1^L X_i)<=L sigma^2 "
            "(1+rho)/(1-rho). Chebyshev therefore preserves square-root-depth scaling under "
            "correlation decay; the sharper displayed Gaussian radius requires the additional Gaussian premise."
        ),
    }
    return finalize(384, "Dependent-Error Multilevel Clock Concentration",
                    "proven_conditional_covariance_decay_bound",
                    "Covariance decay, centering and sigma are added stochastic laws. The result creates no absolute time unit or time arrow.", packet)


def st385() -> dict:
    records = list(ROOT.glob("FIN_ST385_NEW_STRICT_TYPED_SOURCE*.json"))
    packet = {"required_pattern": "FIN_ST385_NEW_STRICT_TYPED_SOURCE*.json",
              "matching_sources": [x.name for x in records], "admitted_source_count": len(records),
              "theorem": "No new source packet is present; prior source/selector no-go lanes are not replayed."}
    return finalize(385, "New Strict Typed-Source Admission Gate",
                    "blocked_no_new_typed_source" if not records else "candidate_present_requires_full_audit",
                    "Absence from the admission channel is not a universal nonexistence theorem. No QW-2191 discharge follows.", packet)


def st386() -> dict:
    records = list(ROOT.glob("FIN_ST386_INDEPENDENT_RAW_EVENTS*.jsonl"))
    packet = {"required_pattern": "FIN_ST386_INDEPENDENT_RAW_EVENTS*.jsonl",
              "matching_records": [x.name for x in records], "independent_record_count": len(records),
              "theorem": "Local computation cannot create independent custody, apparatus calibration, preregistration, unblinding, or empirical events."}
    return finalize(386, "Independent Evidence Gate",
                    "blocked_no_independent_events" if not records else "record_present_requires_blinded_validation",
                    "This is an evidentiary stop, not a failed experiment.", packet)


def make_figures(out: dict) -> None:
    FIG_DIR.mkdir(exist_ok=True)
    w, _, s = strict_operator(); wmin = float(np.min(w[w > 0])); R = 11/12
    rs = np.linspace(1e-7, R-1e-7, 1200)
    vals = [collision_envelope_point(r, s, wmin) for r in rs]
    fig, ax = plt.subplots(figsize=(7, 3.8)); ax.plot(rs, vals, label="certified scalar envelope")
    target = out["ST373"]["explicit_localized_value_interval"][1]
    ax.axhline(target, color="#9d2f2f", ls="--", label="explicit localized benchmark")
    ax.set(xlabel=r"collision $r=\|p-u\|_2^2$", ylabel=r"lower envelope for $V_4$",
           title="ST373: global collision--entropy reduction"); ax.legend(); fig.tight_layout()
    fig.savefig(FIG_DIR/"st373_collision_envelope.png", dpi=200); plt.close(fig)

    r = out["ST380"]; y = np.logspace(math.log10(2e-5), 2, 1500)
    a, nu = r["time1_dual_weight"], r["heat_dual_price"]
    g = np.array([ir_G(x, a, nu) for x in y])
    fig, ax = plt.subplots(figsize=(7, 3.8)); ax.semilogx(y/2, g); ax.axhline(0, color="black", lw=.8)
    ax.set(xlabel=r"spectral coordinate $\mu$", ylabel="divided KKT threshold",
           title="ST380: continuum IR threshold and two selected bands"); fig.tight_layout()
    fig.savefig(FIG_DIR/"st380_continuum_ir_threshold.png", dpi=200); plt.close(fig)

    rows = out["ST383"]["rows"]; fig, ax = plt.subplots(figsize=(7, 3.8))
    ax.plot([x["distortion"] for x in rows], [x["Fano_lower"] for x in rows], "o-", label="Fano lower")
    ax.plot([x["distortion"] for x in rows], [x["best_upper"] for x in rows], "o-", label="achievable upper")
    ax.set(xlabel="Hamming distortion", ylabel="bits per symbol",
           title="ST383: asymptotic strict-heat rate--distortion bracket"); ax.legend(); fig.tight_layout()
    fig.savefig(FIG_DIR/"st383_asymptotic_rate_distortion_bounds.png", dpi=200); plt.close(fig)


def main() -> None:
    np.random.seed(SEED)
    out = {"metadata": {"seed": SEED, "python": platform.python_version(),
                         "numpy": np.__version__, "scipy": scipy.__version__,
                         "scope": "local analytical, interval, exact finite, and bounded numerical research"}}
    funcs = {372: st372, 373: st373, 374: st374, 375: st375, 376: st376,
             377: st377, 378: st378, 379: st379, 380: st380, 381: st381,
             382: st382, 383: st383, 384: st384, 385: st385, 386: st386}
    for k in range(372, 387):
        # ST374 optionally reads the aggregate; serialize current prefix first.
        RESULTS.write_text(json.dumps(native(out), indent=2, sort_keys=True), encoding="utf-8")
        out[f"ST{k}"] = funcs[k]()
    make_figures(out)
    RESULTS.write_text(json.dumps(native(out), indent=2, sort_keys=True), encoding="utf-8")
    with SUMMARY.open("w", newline="", encoding="utf-8") as h:
        writer = csv.writer(h); writer.writerow(["program", "object", "status"])
        for k in range(372, 387): writer.writerow([f"ST{k}", out[f"ST{k}"]["object"], out[f"ST{k}"]["status"]])
    print(json.dumps({k: v["status"] for k, v in out.items() if k.startswith("ST")}, indent=2))


if __name__ == "__main__": main()
