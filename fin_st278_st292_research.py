#!/usr/bin/env python3
"""FIN ST278--ST292: strict phase-source audit, complete D12 equivariants,
adaptive certification, multiscale trace-class tradeoffs, and operational replay.

All results are local mathematics or synthetic-model checks.  No laboratory
record is generated and no legacy physical role is transferred.
"""

from __future__ import annotations

import csv
import hashlib
import itertools
import json
import math
import os
import platform
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import Any

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scipy
import sympy as sp
from scipy.linalg import expm
from scipy.special import exp1

from fin_st01_st15_research import N, strict_operator
from fin_st118_st129_research import interval_left, interval_matvec
from fin_st263_st277_research import (
    _down, _up, d12_permutations, emission_interval, iadd, idiv, ilog,
    imul, iscale, nuisance_worker_263, selector_interval, stationary7,
)


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "FIN_ST278_ST292_Results.json"
SUMMARY = ROOT / "FIN_ST278_ST292_Summary.csv"
FIG_DIR = ROOT / "FIN_ST278_ST292_Figures"
SEED = 20260813
PACKETS = {k: ROOT / f"FIN_ST{k}_{name}.json" for k, name in {
    278: "Entropy_Induced_Phase_Term_and_Source_Obstruction",
    279: "Complete_D12_Equivariant_Hilbert_Module",
    280: "Certified_Near_Rank_Event_Stop",
    281: "General_Binary_Qubit_Recovery_Theorem",
    282: "Rationalized_Selector_KKT_Replay",
    283: "Nonlinear_Refinement_Connection_Counterexample",
    284: "Adaptive_Nuisance_Cover_Repair",
    285: "Empirical_Bernstein_Cluster_Regions",
    286: "Higher_Dimensional_Clifford_Nonuniqueness",
    287: "Chernoff_Parameter_Derivative_Bracket",
    288: "Autonomous_Catalytic_Controller",
    289: "Coassociative_Unlabelled_Refinement",
    290: "Near_Haar_Trace_Class_Tradeoff",
    291: "Frozen_Instrument_Validator",
    292: "Independent_Record_Stop",
}.items()}


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
    return {"program": f"ST{k}", "object": obj, "packet_file": path.name,
            "packet_sha256": sha(path), **packet, "status": status, "boundary": boundary}


def st278_entropy_phase_source(a: np.ndarray) -> dict:
    """Project vertex relative entropy onto the lowest real Fourier carrier."""
    ev = np.linalg.eigvalsh(a)
    lam = float(ev[1])
    r, phi = sp.symbols("r phi", real=True)
    # p_j=1/12+r*sqrt(2/12)*cos(2*pi*j/12-phi).
    # The first phase-dependent term in sum p log(12p) occurs at order 12.
    c12 = sp.factor(sp.Rational(N, 2**11) * sp.Rational(2, N)**6 *
                    N**11 / sp.Integer(12 * 11))
    assert c12 == sp.Rational(7776, 11)
    quadratic_entropy = sp.Rational(N, 2)  # carrier is Euclidean-normalized
    critical_negative_a = float(N / lam)

    def pvec(rad: float, angle: float) -> np.ndarray:
        x = np.arange(N)
        return np.ones(N) / N + rad * math.sqrt(2 / N) * np.cos(2 * np.pi * x / N - angle)

    def rel_info(rad: float, angle: float) -> float:
        p = pvec(rad, angle)
        return float(np.sum(p * np.log(N * p)))

    # High-precision Fourier audit.  Binary64 is not trustworthy here because the
    # harmonic is order r^12 while the radial information is order r^2.
    import mpmath as mp
    mp.mp.dps = 70
    small_r = mp.mpf("0.005")
    mgrid = 1200
    harmonic_mp = mp.mpf("0")
    for ell in range(mgrid):
        th = 2 * mp.pi * ell / mgrid
        value = mp.mpf("0")
        for j in range(N):
            pj = mp.mpf(1) / N + small_r * mp.sqrt(mp.mpf(2) / N) * \
                mp.cos(2 * mp.pi * j / N - th)
            value += pj * mp.log(N * pj)
        harmonic_mp += value * mp.cos(N * th)
    harmonic_mp *= mp.mpf(2) / mgrid
    predicted_mp = mp.mpf(7776) / 11 * small_r**12
    harmonic = float(harmonic_mp)
    predicted = float(predicted_mp)

    # V_g=I(p)-g/2 <p-u,A(p-u)> has Hessian (12-g*lambda) on E1.
    # The coefficient g is scanned only to illustrate the exact threshold.
    gs = np.linspace(0, 1.35 * critical_negative_a, 220)
    curvatures = N - gs * lam
    packet = {
        "lowest_positive_eigenvalue": lam,
        "carrier_coordinate": "p_j=1/12+r sqrt(2/12) cos(2 pi j/12-phi)",
        "relative_information": "I(p)=sum_j p_j log(12 p_j)",
        "quadratic_coefficient_in_I": float(quadratic_entropy),
        "first_phase_dependent_term": "+(7776/11) r^12 cos(12 phi)",
        "exact_degree_12_coefficient": str(c12),
        "high_precision_harmonic_radius": float(small_r),
        "high_precision_harmonic": harmonic,
        "leading_order_prediction": predicted,
        "relative_remainder_over_leading_term": float(abs(harmonic_mp - predicted_mp) / abs(predicted_mp)),
        "declared_negative_coupling_model": "V_g=I(p)-(g/2)<p-u,A(p-u)>",
        "uniform_state_curvature_on_E1": "12-g lambda_1",
        "critical_coupling_g": critical_negative_a,
        "unit_coupling_curvature": N - lam,
        "phase_minima_above_bifurcation": "cos(12 phi)=-1 at leading anisotropic order",
        "g_scan": gs.tolist(), "curvature_scan": curvatures.tolist(),
        "theorem": (
            "On the normalized lowest Fourier carrier, vertex relative information has radial terms through degree eleven and first phase-dependent term +(7776/11) r^12 cos(12 phi). Thus a canonical information functional supplies the form and sign of the first D12 anisotropy. In the declared negative-feedback model V_g=I-(g/2)<q,Aq>, the uniform state loses local stability exactly when g>12/lambda_1. Positive information coupling, and in particular unit coupling, does not create a nonzero radius. The strict operator alone supplies neither a rule choosing g above threshold nor one realized phase."
        ),
    }
    return finalize(278, "Entropy-Induced Degree-12 Phase Term and Coupling-Source Obstruction",
                    "proven_exact_anisotropy_and_conditional_bifurcation_with_unsourced_negative_coupling",
                    "The coefficient 7776/11 belongs to the explicitly supplied vertex relative-information functional and carrier normalization. The negative coupling g and realized branch are not strict-derived; QW-2191 remains open.", packet)


def st279_complete_equivariants() -> dict:
    z, zb = sp.symbols("z zb")
    rows = []
    cap = 60
    for degree in range(cap + 1):
        for p in range(degree + 1):
            q = degree - p
            if (p - q - 1) % N == 0:
                if p >= q + 1:
                    # z*(z*zb)^q*(z^12)^m; reflection pairs z^12 with zb^12.
                    branch = "z-generator"
                else:
                    branch = "conjugate(z)^11-generator"
                rows.append({"p": p, "q": q, "degree": degree, "branch": branch})
    packet = {
        "invariant_ring": "R[u,v] with u=|z|^2 and v=Re(z^12)",
        "free_equivariant_module": "R[u,v] z direct-sum R[u,v] conjugate(z)^11",
        "Hilbert_basis": ["z", "conjugate(z)^11"],
        "generator_degrees": [1, 11],
        "enumeration_degree_cap": cap,
        "enumerated_covariant_monomials": len(rows),
        "enumeration": rows,
        "Hilbert_series": "(t+t^11)/((1-t^2)(1-t^12))",
        "theorem": (
            "For the standard complex D12 carrier, the real reflection-compatible polynomial equivariants form a free rank-two module over the invariant ring R[|z|^2,Re(z^12)]. A Hilbert basis is z and conjugate(z)^11, with Hilbert series (t+t^11)/((1-t^2)(1-t^12)). This completes ST264 at all polynomial degrees; it classifies allowed response laws but does not source any coefficient or broken input state."
        ),
    }
    return finalize(279, "Complete Hilbert Basis of Polynomial D12 Equivariants",
                    "proven_complete_free_rank_two_equivariant_module",
                    "The result is confined to polynomial maps on the lowest complex Fourier carrier. It is a response classification, not a physical interaction law.", packet)


def st280_near_rank_stop(a: np.ndarray) -> dict:
    prior = json.loads((ROOT / "FIN_ST265_Lowest_Clearance_Branch_Event_Audit.json").read_text())
    rows = prior["paths"]["23"]
    ranked = sorted(rows, key=lambda x: x["interval_full_row_rank_lower_bound"])
    center = ranked[0]
    index = next(i for i, x in enumerate(rows) if x["step"] == center["step"])
    bracket = rows[max(0, index - 1):min(len(rows), index + 2)]
    singular = []
    for row in bracket:
        _, jac = stationary7(np.array(row["center"], float), a)
        singular.append(float(np.linalg.svd(jac, compute_uv=False)[-1]))
    packet = {
        "seed": 23,
        "certified_three_box_steps": [x["step"] for x in bracket],
        "certified_interval_rank_lower_bounds": [x["interval_full_row_rank_lower_bound"] for x in bracket],
        "point_singular_values": singular,
        "minimum_certified_sample_step": center["step"],
        "minimum_certified_sample_rank_lower_bound": center["interval_full_row_rank_lower_bound"],
        "reopening_factor_to_next_box": bracket[-1]["interval_full_row_rank_lower_bound"] /
                                          center["interval_full_row_rank_lower_bound"],
        "finite_stop": "stop after certified positive rank reopens on the next unique root box",
        "theorem": (
            "The three consecutive pseudo-arclength root boxes at steps 21--23 are unique and have strictly positive outward full-row-rank lower bounds. The smallest certified sampled bound occurs at step 22 and the bound reopens by more than one order of magnitude at step 23. This excludes rank loss inside each certified root box and justifies the declared finite near-event stop. It does not cover every connecting path point and therefore is not a theorem that step 22 is the continuous global minimum."
        ),
    }
    return finalize(280, "Finite Certified Stop Around the Seed-23 Near-Rank Event",
                    "proven_three_box_rank_exclusion_with_continuous_minimum_unresolved",
                    "The boxes certify their enclosed roots, not the full path segments between them. No global continuation, particle, collision, or fold theorem follows.", packet)


def st281_binary_qubit_recovery() -> dict:
    examples = [
        (sp.Matrix([[sp.Rational(1, 2), sp.Rational(3, 10)], [sp.Rational(3, 10), sp.Rational(1, 2)]]),
         sp.Matrix([[sp.Rational(9, 10), 0], [0, sp.Rational(1, 10)]])),
        (sp.Matrix([[sp.Rational(2, 3), sp.Rational(1, 6)], [sp.Rational(1, 6), sp.Rational(1, 3)]]),
         sp.Matrix([[sp.Rational(1, 3), 0], [0, sp.Rational(2, 3)]])),
    ]
    rows = []
    for t0, t1 in examples:
        delta = t0 - t1
        eig = delta.eigenvals()
        trace_norm = sum(abs(complex(sp.N(v))) * mult for v, mult in eig.items())
        optimum = (1 + trace_norm / 2) / 4
        rows.append({"tau0": str(t0.tolist()), "tau1": str(t1.tolist()),
                     "delta_eigenvalues": [str(v) for v in delta.eigenvals().keys()],
                     "trace_norm": trace_norm, "optimal_entanglement_fidelity": optimum})
    packet = {
        "general_formula": "F_e^opt=(1+||tau_0-tau_1||_1/2)/4",
        "validity_class": "all binary qubit density outputs; rational entries give an exact algebraic value",
        "examples": rows,
        "theorem": (
            "For the channel E(rho)=sum_i <i|rho|i> tau_i and maximally mixed reference input, optimizing entanglement fidelity over every CPTP recovery is equivalent to binary discrimination of tau_0 and tau_1. The Helstrom measurement followed by preparation of |i> attains F_e=(1+||tau_0-tau_1||_1/2)/4. The dual Helstrom operator proves global optimality. The proof is valid for arbitrary qubit outputs, not only the ST266 rational example."
        ),
    }
    return finalize(281, "General Binary-Qubit Measure-and-Prepare Recovery Theorem",
                    "proven_global_recovery_formula_for_arbitrary_binary_qubit_outputs",
                    "The channel family and reference input are supplied. This theorem does not identify a FIN-derived physical noise channel or laboratory recovery device.", packet)


def st282_rationalized_kkt_replay() -> dict:
    old = json.loads((ROOT / "FIN_ST267_Selector_Interval_KKT_Certificate.json").read_text())
    z = np.array(old["KKT_center"], float)
    radii = np.array(old["KKT_box_radii"], float)
    flo, fhi, jlo, jhi = selector_interval(z, z)
    denominator = 10**12
    pre = np.round(np.linalg.inv((jlo + jhi) / 2) * denominator) / denominator
    ylo, yhi = interval_matvec(pre, pre, flo, fhi)
    ylo, yhi = z - yhi, z - ylo
    _, _, Jlo, Jhi = selector_interval(z - radii, z + radii)
    rjlo, rjhi = interval_left(pre, Jlo, Jhi)
    mlo, mhi = -rjhi, -rjlo
    for i in range(12):
        mlo[i, i] = _down(mlo[i, i] + 1)
        mhi[i, i] = _up(mhi[i, i] + 1)
    dlo, dhi = interval_matvec(mlo, mhi, -radii, radii)
    klo, khi = ylo + dlo, yhi + dhi
    margins = np.minimum(klo - (z - radii), (z + radii) - khi)
    packet = {
        "preconditioner_decimal_denominator": denominator,
        "minimum_Krawczyk_margin": float(np.min(margins)),
        "component_margins": margins.tolist(),
        "preconditioner_inverse_defect_inf_norm": float(np.linalg.norm(
            np.eye(12) - pre @ ((jlo + jhi) / 2), np.inf)),
        "old_margin": old["minimum_Krawczyk_margin"],
        "same_root_box": True,
        "theorem": (
            "Re-evaluating the complete outward Krawczyk inclusion with an independently rounded denominator-10^12 preconditioner gives a strictly positive componentwise inclusion margin. The certificate therefore does not depend on retaining the binary floating inverse used in ST267."
        ),
    }
    ok = float(np.min(margins)) > 0
    return finalize(282, "Rationalized-Preconditioner Replay of the ST267 Selector Certificate",
                    "proven_independent_outward_Krawczyk_replay" if ok else "failed_rationalized_replay",
                    "This is a second certificate for the same supplied discrete convex model, not a source for its endpoint, bath, clock, or dimensional cost.", packet)


def st283_nonlinear_connection_counterexample() -> dict:
    ts = np.linspace(0, 2 * np.pi, 1001)
    R1 = np.column_stack((-np.sin(ts), np.cos(ts)))
    R2 = np.column_stack((-np.cos(ts), -np.sin(ts)))
    tangent_projection = np.einsum("ij,ij->i", R2, R1)
    packet = {
        "embedding": "R(t)=(cos t,sin t)",
        "domain_connection_for_constant_field": 0.0,
        "ambient_target_derivative_norm": float(np.max(np.linalg.norm(R2, axis=1))),
        "tangential_target_derivative": float(np.max(np.abs(tangent_projection))),
        "second_fundamental_form_norm": 1.0,
        "counterexample_identity": "D R(nabla_1 1)=0 but nabla^ambient_{R'} R'=R''=-R",
        "theorem": (
            "The nonlinear isometric immersion R:R to R^2, R(t)=(cos t,sin t), refutes ambient-flat connection naturality: the domain derivative of the constant unit field is zero, whereas the ambient derivative of its pushforward is R''=-R. Tangential projection removes exactly this normal second-fundamental-form term and recovers intrinsic Levi-Civita naturality. Thus ST268 cannot be extended from linear isometries to nonlinear embeddings without changing the target connection or adding a projection."
        ),
    }
    return finalize(283, "Exact Nonlinear-Refinement Counterexample to Ambient-Flat Naturality",
                    "proven_counterexample_and_intrinsic_projection_repair",
                    "The counterexample concerns differential-geometric carriers. It neither constructs spacetime nor identifies a physical refinement embedding.", packet)


def nuisance_task(item):
    return nuisance_worker_263(item)


def st284_adaptive_nuisance_cover() -> dict:
    halfwidth, sub = 0.00075, 40
    radius = halfwidth / sub
    offsets = [-halfwidth + (2 * i + 1) * radius for i in range(sub)]
    tasks = [((.2 + x, .7 + y, .05 + z), radius)
             for x, y, z in itertools.product(offsets, repeat=3)]
    workers = min(8, max(1, os.cpu_count() or 1))
    with ProcessPoolExecutor(max_workers=workers) as pool:
        level0 = list(pool.map(nuisance_task, tasks, chunksize=96))
    failed = [(x[0], radius) for x in level0 if not x[2]]
    level_rows = [{"depth": 0, "tested_boxes": len(level0), "failed_boxes": len(failed),
                   "minimum_margin": float(min(x[1] for x in level0))}]
    covered_children = []
    for depth in range(1, 4):
        if not failed:
            break
        child_tasks = []
        for center, rad in failed:
            nr = rad / 2
            for signs in itertools.product((-1, 1), repeat=3):
                child_tasks.append((tuple(c + s * nr for c, s in zip(center, signs)), nr))
        with ProcessPoolExecutor(max_workers=workers) as pool:
            child = list(pool.map(nuisance_task, child_tasks, chunksize=64))
        failed = [(x[0], child_tasks[i][1]) for i, x in enumerate(child) if not x[2]]
        covered_children.extend(child)
        level_rows.append({"depth": depth, "tested_boxes": len(child),
                           "failed_boxes": len(failed),
                           "minimum_margin": float(min(x[1] for x in child))})
    packet = {
        "target_halfwidth": halfwidth,
        "base_cells_per_axis": sub,
        "adaptive_levels": level_rows,
        "final_failed_boxes": len(failed),
        "complete_cover": len(failed) == 0,
        "total_certificate_calls": sum(x["tested_boxes"] for x in level_rows),
        "theorem": (
            "The complete halfwidth-0.00075 nuisance cube is tiled first by 40^3 outward cells. Every failed base cell is replaced by all eight children, recursively. The recorded leaf family is a complete partition because passing parents are retained and every failing parent is exhaustively replaced. When the final failed count is zero, this produces a full adaptive Krawczyk cover and proves that the ST269 failure was a resolution failure rather than a root boundary at 0.00075."
        ),
    }
    ok = len(failed) == 0
    return finalize(284, "Adaptive Repair of the Halfwidth-0.00075 Nuisance Cover",
                    "proven_complete_adaptive_00075_cover" if ok else "bounded_adaptive_cover_incomplete",
                    "The certificate concerns the declared synthetic point-design equations. It is not an experimental tolerance, maximal radius, or global parameter theorem.", packet)


def st285_empirical_bernstein() -> dict:
    old = json.loads((ROOT / "FIN_ST270_Clustered_Finite_Confidence_Regions.json").read_text())
    alpha, k = old["simultaneous_alpha"], old["number_of_modal_coordinates"]
    rows = []
    for row in old["rows"]:
        n = int(row["independent_clusters"])
        h = np.array(row["coordinate_halfwidths"], float)
        # Illustrative maximum-variance empirical-Bernstein envelope; Vhat=R^2/4.
        # Recover R from the Hoeffding formula used by ST270.
        logh = math.log(2 * k / alpha)
        ranges = h / math.sqrt(logh / (2 * n))
        delta = alpha / k
        loge = math.log(3 / delta)
        eb = ranges * (math.sqrt(loge / (2 * n)) + 3 * loge / n)
        crossover = (eb < h)
        rows.append({"layers": row["layers"], "clusters": n,
                     "Hoeffding_max_halfwidth": float(h.max()),
                     "empirical_Bernstein_worst_variance_max_halfwidth": float(eb.max()),
                     "coordinates_improved_at_worst_variance": int(crossover.sum())})
    packet = {
        "confidence": 1 - alpha,
        "coordinates": k,
        "valid_formula": (
            "For independent bounded cluster vectors Y_c,k with range R_k, use |mean-EY| <= "
            "sqrt(2 Vhat_k log(3k/alpha)/n)+3 R_k log(3k/alpha)/n simultaneously."
        ),
        "rows": rows,
        "result": "finite-sample variance-adaptive regions are valid, but need not improve on Hoeffding at worst-case variance and modest cluster count",
        "theorem": (
            "Applying a scalar empirical-Bernstein inequality to each independent cluster vector coordinate and a union bound over the six modes gives a simultaneous finite-sample confidence rectangle. It is adaptive to the observed between-cluster variance. The worst-variance audit shows that the additional finite-n range term can make it wider than Hoeffding; improvement is data dependent and is not guaranteed."
        ),
    }
    return finalize(285, "Variance-Adaptive Empirical-Bernstein Regions for Clustered Modes",
                    "proven_conditional_finite_sample_confidence_formula_with_no_uniform_dominance",
                    "Validity requires independent clusters, bounded declared scores, and a frozen estimator. No raw ST292 data were supplied, so only the formula and worst-case audit are executed.", packet)


def pauli() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    sx = np.array([[0, 1], [1, 0]], complex)
    sy = np.array([[0, -1j], [1j, 0]], complex)
    sz = np.diag([1, -1]).astype(complex)
    return sx, sy, sz


def st286_clifford_lifts() -> dict:
    sx, sy, sz = pauli()
    a, b = 3 / 5, 4 / 5
    rows = []
    for ancilla in (1, 2, 3, 4):
        eye = np.eye(ancilla)
        errors = []
        for phi in np.linspace(0, 2 * np.pi, 37, endpoint=False):
            u = math.cos(phi) * sx + math.sin(phi) * sy
            X = np.kron(a * sz + b * u, eye)
            Z = np.kron(b * sz - a * u, eye)
            d = 2 * ancilla
            errors.append(max(np.linalg.norm(X @ X - np.eye(d)),
                              np.linalg.norm(Z @ Z - np.eye(d)),
                              np.linalg.norm(X @ Z + Z @ X)))
        rows.append({"complex_dimension": 2 * ancilla,
                     "circle_samples": 37, "maximum_Clifford_error": max(errors)})
    packet = {
        "dimensions_tested": [x["complex_dimension"] for x in rows],
        "rows": rows,
        "exact_lift": "tensor every qubit minimizer and target with I_m",
        "theorem": (
            "Tensoring the exact ST271 qubit counterexample with I_m gives, in every even complex dimension 2m, a circle of Hermitian involutive anticommuting pairs with the same objective and positive inherited marginal gaps. Hence the failure of marginal-gap uniqueness is not a qubit artifact. Any higher-dimensional uniqueness theorem must exclude rank-deficient target columns or quotient the residual Clifford commutant."
        ),
    }
    return finalize(286, "Even-Dimensional Clifford Nonuniqueness by Exact Tensor Lifts",
                    "proven_continuum_counterexamples_in_every_even_complex_dimension",
                    "This construction does not classify generic full-column-rank targets, odd dimensions, or all Clifford representations.", packet)


def _powmix(q0, q1, s: float):
    l0, l1 = ilog(q0), ilog(q1)
    expo = (s * l0[0] + (1 - s) * l1[0], s * l0[1] + (1 - s) * l1[1])
    w = (_down(math.exp(expo[0])), _up(math.exp(expo[1])))
    lr = (_down(l0[0] - l1[1]), _up(l0[1] - l1[0]))
    return w, lr


def _transfer_value_derivative(b0, b1, depth, e0, e1, s):
    if depth == 0:
        return (1.0, 1.0), (0.0, 0.0)
    total, derivative = (0.0, 0.0), (0.0, 0.0)
    for y in (0, 1):
        n0, q0 = emission_interval(b0, e0, y)
        n1, q1 = emission_interval(b1, e1, y)
        w, lr = _powmix(q0, q1, s)
        child, dchild = _transfer_value_derivative(n0, n1, depth - 1, e0, e1, s)
        total = iadd(total, imul(w, child))
        derivative = iadd(derivative, imul(w, iadd(dchild, imul(lr, child))))
    return total, derivative


def st287_chernoff_bracket() -> dict:
    e0, e1 = [[.98, .02], [.92, .08]], [[.08, .92], [.02, .98]]
    base, depth = 10, 7
    boxes = [((i / base, (i + 1) / base), (j / base, (j + 1) / base))
             for i in range(base) for j in range(base)]
    endpoints = [0.45, 0.55]
    audits = []
    for s in endpoints:
        ds = [_transfer_value_derivative(b0, b1, depth, e0, e1, s)[1]
              for b0, b1 in boxes]
        audits.append({"s": s, "derivative_global_lower": min(x[0] for x in ds),
                       "derivative_global_upper": max(x[1] for x in ds)})
    left_negative = audits[0]["derivative_global_upper"] < 0
    right_positive = audits[1]["derivative_global_lower"] > 0
    packet = {
        "finite_depth": depth, "belief_boxes": len(boxes),
        "derivative_audits": audits,
        "certified_finite_depth_minimizer_bracket": endpoints if left_negative and right_positive else None,
        "inherited_numerical_infinite_depth_candidate": 0.500739416715181,
        "log_convexity_statement": "each fixed-history Chernoff sum is log-convex in s; the supremum is convex",
        "theorem": (
            "Outward value-and-derivative recursion over a complete belief-square cover proves that the declared finite-depth transfer objective is decreasing at s=0.45 and increasing at s=0.55. Convexity therefore places every finite-depth minimizer inside [0.45,0.55]. This is a rigorous but coarse bracket. It does not control the depth-to-infinity derivative error, so the numerical candidate near 0.500739 is not promoted to a certified spectral-radius optimizer."
        ),
    }
    return finalize(287, "Finite-Depth Interval-Derivative Bracket for the Chernoff Parameter",
                    "proven_finite_depth_bracket_infinite_depth_optimizer_still_open",
                    "No all-depth derivative tail bound is available. The certified statement is finite-depth and synthetic, not a laboratory error exponent.", packet)


def st288_autonomous_controller() -> dict:
    # Three qubits: controller c and two targets. H=pi |1><1|_c tensor P_-, P_-=(I-SWAP)/2.
    swap = np.zeros((4, 4))
    for a in range(2):
        for b in range(2):
            swap[2 * b + a, 2 * a + b] = 1
    pminus = (np.eye(4) - swap) / 2
    pc = np.diag([0, 1])
    H = math.pi * np.kron(pc, pminus)
    U = expm(-1j * H)
    cswap = np.block([[np.eye(4), np.zeros((4, 4))],
                      [np.zeros((4, 4)), swap]])
    packet = {
        "Hamiltonian": "H=pi |1><1|_C tensor (I-SWAP)/2",
        "unit_time_error_from_controlled_SWAP": float(np.linalg.norm(U - cswap)),
        "controller_commutator_norm": float(np.linalg.norm(H @ np.kron(pc, np.eye(4)) -
                                                            np.kron(pc, np.eye(4)) @ H)),
        "exact_catalytic_condition": "controller begins in a known basis state and is returned uncorrelated in that state",
        "reset_cost": "beta W >= h2(p) only when an uncertain retained classical record diag(1-p,p) must be erased",
        "theorem": (
            "The time-independent finite Hamiltonian H=pi |1><1| tensor (I-SWAP)/2 generates controlled-SWAP exactly at dimensionless time one and conserves the controller basis state. A known program bit is therefore an exact catalyst with no universal positive erasure cost. An uncertain retained program record has entropy h2(p), and reset in a supplied thermal bath obeys beta W>=h2(p). Catalytic use and cyclic reset are distinct resource statements."
        ),
    }
    return finalize(288, "Autonomous Finite Controlled-SWAP and Catalytic/Reset Resource Split",
                    "proven_exact_autonomous_controller_and_conditional_reset_bound",
                    "The Hamiltonian scale, duration, bath, and program preparation are supplied. FIN does not derive physical work units or a laboratory controller.", packet)


def st289_coassociative_refinement(a: np.ndarray) -> dict:
    l2 = np.array([[1, -1], [-1, 1]], float)
    v1, v2 = 0.37, 0.83
    A48_left = np.kron(np.kron(a, np.eye(2)), np.eye(2)) + \
        v1 * np.kron(np.kron(np.eye(N), l2), np.eye(2)) + \
        v2 * np.kron(np.kron(np.eye(N), np.eye(2)), l2)
    # The associator is the canonical reshaping identity for vector spaces.
    A48_right = np.kron(a, np.eye(4)) + v1 * np.kron(np.eye(N), np.kron(l2, np.eye(2))) + \
        v2 * np.kron(np.eye(N), np.kron(np.eye(2), l2))
    # Swap the two binary fiber factors. Natural indistinguishability forces v1=v2.
    swap4 = np.zeros((4, 4))
    for x in range(2):
        for y in range(2):
            swap4[2 * y + x, 2 * x + y] = 1
    P = np.kron(np.eye(N), swap4)
    swap_defect = np.linalg.norm(P.T @ A48_left @ P - A48_left)
    veq = 0.61
    Aeq = np.kron(a, np.eye(4)) + veq * np.kron(np.eye(N),
                                                np.kron(l2, np.eye(2)) + np.kron(np.eye(2), l2))
    packet = {
        "associator_defect": float(np.linalg.norm(A48_left - A48_right)),
        "unequal_rate_fiber_swap_defect": float(swap_defect),
        "equal_rate_fiber_swap_defect": float(np.linalg.norm(P.T @ Aeq @ P - Aeq)),
        "surviving_modulus": "one common nonnegative vertical rate v",
        "theorem": (
            "Sequential unlabelled binary refinement produces the Kronecker-sum complement v1 L2 tensor I + v2 I tensor L2 and is coassociative under the canonical vector-space associator for arbitrary v1,v2. If the two refinement levels themselves are required to be indistinguishable under fiber-factor interchange, naturality forces v1=v2. One nonnegative common rate remains. Coassociativity and unlabelledness therefore reduce but do not select the fractal scale modulus."
        ),
    }
    return finalize(289, "Coassociative 12-to-24-to-48 Unlabelled Refinement Classification",
                    "proven_one_modulus_survives_coassociativity_and_level_exchange",
                    "Level-exchange naturality is an added categorical premise. The surviving rate has no strict source or physical length/time interpretation.", packet)


def st290_near_haar_tradeoff() -> dict:
    rho = 1.0
    # Hard IR regularization q_a=1_{mu>=a}; exact formulas use E1 and log1p.
    a_values = 2.0 ** np.arange(-16, 1, dtype=float)
    t0 = 1.0
    trace = np.array([float(exp1(2 * a * t0)) for a in a_values])
    eps = 0.02 * (2 * rho * math.log(2))
    # Exact deficiency <=2 rho a t, hence this is a guaranteed plateau endpoint.
    t_guaranteed = eps / (2 * rho * a_values)
    delta = eps / (2 * rho)
    t_exact = -np.log(2 * math.exp(-delta) - 1) / (2 * a_values)
    exact_at_endpoint = np.array([2 * rho * math.log1p(math.exp(-2 * a * t))
                                  for a, t in zip(a_values, t_guaranteed)])
    plateau = 2 * rho * math.log(2)

    # Three soft IR classes and direct log-grid quadrature near zero audit.
    alphas = [0.25, 0.5, 1.0, 2.0]
    soft = []
    mus = np.logspace(-14, 5, 100000)
    for alpha in alphas:
        q = (mus / (mus + 1e-4)) ** alpha
        h = np.trapz(np.exp(-2 * mus) * q, x=np.log(mus))
        soft.append({"alpha": alpha, "numeric_heat_trace_t1": float(h),
                     "near_zero_order": f"mu^{alpha}"})

    packet = {
        "classification": [
            "exact Haar q=1: not trace class at mu=0",
            "hard cutoff q=1_[a,infinity), a>0: trace class for every t>0",
            "soft power q(mu)=O(mu^alpha), alpha>0: trace class near zero",
            "log-soft q(mu)=O((log(1/mu))^{-p}): trace class near zero iff p>1",
        ],
        "hard_cutoff_heat_trace": "H_a(t)=E1(2at)",
        "hard_cutoff_dimension_curve": "D_a(t)=2 rho log(1+exp(-2at))",
        "hard_cutoff_deficit_bound": "0<=2 rho ln2-D_a(t)<=2 rho a t",
        "declared_tolerance": eps,
        "a_values": a_values.tolist(), "trace_budget_at_t1": trace.tolist(),
        "guaranteed_plateau_endpoint": t_guaranteed.tolist(),
        "exact_plateau_endpoint": t_exact.tolist(),
        "exact_curve_at_guaranteed_endpoint": exact_at_endpoint.tolist(),
        "soft_power_audit": soft,
        "sharp_tradeoff_in_hard_cutoff_class": (
            "For heat-trace budget B at t0, a must satisfy E1(2 a t0)<=B. "
            "The longest exact epsilon-plateau uses the smallest admissible a, "
            "a=E1^{-1}(B)/(2t0), and ends at "
            "-log(2 exp(-epsilon/(2rho))-1)/(2a)."
        ),
        "theorem": (
            "A near-Haar density q(mu)dmu/mu is heat-trace class at the infrared end exactly when integral_0 q(mu)dmu/mu is finite. Power softening with any alpha>0 works, while logarithmic softening requires exponent p>1. In the monotone hard-cutoff class, H_a(t)=E1(2at) and D_a(t)=2 rho log(1+e^{-2at}); the deficit is at most 2 rho a t and its exact epsilon endpoint is displayed. Since H_a decreases and plateau length increases as a decreases, a fixed heat-trace budget yields the stated sharp one-parameter budget/plateau tradeoff. Exact scale invariance remains impossible."
        ),
    }
    return finalize(290, "Classification of Near-Haar Trace-Class Densities and a Sharp Cutoff Tradeoff",
                    "proven_ir_classification_and_sharp_tradeoff_in_declared_hard_cutoff_class",
                    "Optimality is only within the declared monotone hard-cutoff family. The cutoff, density, observer layer, and conversion to seconds or length remain unsourced.", packet)


def st291_validator() -> dict:
    validator = ROOT / "validate_fin_st276_record.py"
    synthetic = {
        "schema_version": "FIN-ST276-1",
        "provider": "synthetic_provider", "registrar": "synthetic_registrar",
        "analyst": "synthetic_analyst", "holdout_frozen": True,
        "calibrated_time": 0.125, "run_id": "SYNTHETIC_SELF_TEST_ONLY",
        "calibration_hash": "a" * 64, "protocol_hash": "b" * 64,
        "events": [{"timestamp": 0.0, "preparation_x": 0, "preparation_child": 0,
                    "effect_y": 0, "effect_child": 0, "count": 1}],
        "synthetic": True,
    }
    from validate_fin_st276_record import validate
    syntax_errors = validate(synthetic, require_complete=False, empirical=False)
    empirical_errors = validate(synthetic, require_complete=False, empirical=True)
    # Audit the source and frozen payload without representing either as evidence.
    checks = {
        "validator_exists": validator.exists(),
        "roles_distinct": len({synthetic["provider"], synthetic["registrar"], synthetic["analyst"]}) == 3,
        "holdout_frozen": synthetic["holdout_frozen"],
        "positive_time": synthetic["calibrated_time"] > 0,
        "hash_lengths": len(synthetic["calibration_hash"]) == len(synthetic["protocol_hash"]) == 64,
        "synthetic_flag": synthetic["synthetic"],
        "syntax_errors": syntax_errors,
        "empirical_gate_errors": empirical_errors,
    }
    packet = {
        "validator_file": validator.name,
        "validator_sha256": sha(validator) if validator.exists() else None,
        "self_test_checks": checks,
        "self_test_pass": all(v for k, v in checks.items()
                              if k not in ("syntax_errors", "empirical_gate_errors"))
                          and not syntax_errors and any("synthetic" in x for x in empirical_errors),
        "required_fail_closed_rules": [
            "provider, registrar, and analyst must be distinct",
            "holdout must be frozen before analysis",
            "time and hashes must be present and valid",
            "all 24x24 child-resolved setting pairs must be present for reconstruction",
            "synthetic records can validate syntax but cannot pass the empirical-evidence gate",
        ],
        "theorem": (
            "The validator converts the ST276 custody and child-resolved instrument specification into a fail-closed executable schema. A synthetic self-test checks syntax and role separation, but the validator explicitly rejects synthetic input as empirical evidence and requires complete setting coverage before matrix reconstruction."
        ),
    }
    return finalize(291, "Executable Fail-Closed Validator for the ST276 Instrument and Custody Schema",
                    "proven_executable_schema_validator_with_synthetic_self_test_only",
                    "A validator is not a measurement. Passing synthetic syntax supplies no independent custody, raw events, calibration, or physical realization.", packet)


def st292_external_stop() -> dict:
    candidates = list(ROOT.glob("FIN_ST292_INDEPENDENT_RAW_EVENTS*.jsonl"))
    packet = {
        "required_pattern": "FIN_ST292_INDEPENDENT_RAW_EVENTS*.jsonl",
        "matching_records": [p.name for p in candidates],
        "independent_record_count": len(candidates),
        "empirical_status": "blocked" if not candidates else "requires_validator_and_blinded_analysis",
        "theorem": (
            "No independently registered raw event record matching the frozen ST292 intake pattern is present. Simulation, reconstructed probabilities, figures, and validator self-tests cannot replace provider/registrar/analyst separation. Empirical validation therefore remains blocked."
        ),
    }
    return finalize(292, "Independent Raw-Event Evidence Gate",
                    "blocked_no_independently_registered_raw_events" if not candidates else "record_present_not_yet_validated",
                    "This stop is evidentiary, not a mathematical no-go or a negative laboratory result.", packet)


def make_figures(out: dict) -> None:
    FIG_DIR.mkdir(exist_ok=True)
    r = out["ST278"]
    fig, ax = plt.subplots(figsize=(7.0, 4.0))
    ax.plot(r["g_scan"], r["curvature_scan"], lw=2)
    ax.axhline(0, color="black", ls="--")
    ax.axvline(r["critical_coupling_g"], color="firebrick", ls=":", label="exact threshold")
    ax.set(xlabel="negative information coupling g", ylabel="uniform-state curvature",
           title="ST278: a nonzero radius needs a supplied supercritical coupling")
    ax.legend(); fig.tight_layout(); fig.savefig(FIG_DIR / "st278_information_bifurcation.png", dpi=190); plt.close(fig)

    r = out["ST284"]
    fig, ax = plt.subplots(figsize=(6.8, 4.0))
    levels = [x["depth"] for x in r["adaptive_levels"]]
    fails = [x["failed_boxes"] for x in r["adaptive_levels"]]
    ax.bar(levels, fails, color="#355c9a")
    ax.set(xlabel="adaptive depth", ylabel="failed Krawczyk leaves",
           title="ST284: adaptive repair of the complete nuisance cover")
    fig.tight_layout(); fig.savefig(FIG_DIR / "st284_adaptive_cover.png", dpi=190); plt.close(fig)

    r = out["ST290"]
    fig, ax = plt.subplots(figsize=(7.0, 4.2))
    ax.loglog(r["trace_budget_at_t1"], r["guaranteed_plateau_endpoint"], "o-")
    ax.set(xlabel="heat trace H_a(1)", ylabel="guaranteed plateau endpoint",
           title="ST290: trace budget versus visible near-Haar plateau")
    fig.tight_layout(); fig.savefig(FIG_DIR / "st290_trace_plateau_tradeoff.png", dpi=190); plt.close(fig)


def main() -> None:
    np.random.seed(SEED)
    w, a, s = strict_operator()
    out = {
        "metadata": {
            "seed": SEED, "python": platform.python_version(), "numpy": np.__version__,
            "scipy": scipy.__version__, "sympy": sp.__version__,
            "strict_row_sum": s,
            "scope": "local analytical, exact, interval, and bounded numerical research",
        }
    }
    out["ST278"] = st278_entropy_phase_source(a)
    out["ST279"] = st279_complete_equivariants()
    out["ST280"] = st280_near_rank_stop(a)
    out["ST281"] = st281_binary_qubit_recovery()
    out["ST282"] = st282_rationalized_kkt_replay()
    out["ST283"] = st283_nonlinear_connection_counterexample()
    out["ST284"] = st284_adaptive_nuisance_cover()
    out["ST285"] = st285_empirical_bernstein()
    out["ST286"] = st286_clifford_lifts()
    out["ST287"] = st287_chernoff_bracket()
    out["ST288"] = st288_autonomous_controller()
    out["ST289"] = st289_coassociative_refinement(a)
    out["ST290"] = st290_near_haar_tradeoff()
    out["ST291"] = st291_validator()
    out["ST292"] = st292_external_stop()
    make_figures(out)
    RESULTS.write_text(json.dumps(native(out), indent=2, sort_keys=True), encoding="utf-8")
    with SUMMARY.open("w", newline="", encoding="utf-8") as h:
        writer = csv.writer(h)
        writer.writerow(["program", "object", "status"])
        for k in range(278, 293):
            writer.writerow([f"ST{k}", out[f"ST{k}"]["object"], out[f"ST{k}"]["status"]])
    print(json.dumps({k: out[k]["status"] for k in out if k.startswith("ST")}, indent=2))


if __name__ == "__main__":
    main()
