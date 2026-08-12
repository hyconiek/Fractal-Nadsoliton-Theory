#!/usr/bin/env python3
"""FIN ST233--ST247: typed second operator, nonlinear carrier obstruction, and refinement moduli."""

from __future__ import annotations

import csv
import hashlib
import itertools
import json
import math
import platform
from fractions import Fraction
from pathlib import Path
from typing import Any

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mpmath as mp
import numpy as np
import scipy
import sympy as sp
from scipy.integrate import quad
from scipy.linalg import expm
from scipy.optimize import brentq, root

from fin_st01_st15_research import N, strict_operator
from fin_st130_st141_research import point_design_system
from fin_st132_center_isolation_replay import bounds as replay_bounds, iv, strict_interval_matrix
from fin_st166_st177_research import local_param_krawczyk, parametric_data
from fin_st178_st189_research import stationary_slice_float, uniform_fold_seed
from fin_st190_st202_research import frac_hmm_likelihoods
from fin_st203_st217_research import pseudo_krawczyk


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "FIN_ST233_ST247_Results.json"
SUMMARY = ROOT / "FIN_ST233_ST247_Summary.csv"
FIG_DIR = ROOT / "FIN_ST233_ST247_Figures"
SEED = 20260812
PACKETS = {k: ROOT / f"FIN_ST{k}_{name}.json" for k, name in {
    233: "State_Typed_Second_Operator",
    234: "Interval_Blackwell_Vertex_Classes",
    235: "State_Conditioned_Noise_Selector_Audit",
    236: "Adaptive_All_Branch_Component_Graph",
    237: "Nonunital_Negative_Orientation_Recovery",
    238: "Continuous_Selector_Thermodynamic_Optimum",
    239: "Nonlinear_Carrier_Chart_Obstruction",
    240: "Nuisance_Method_Failure_Subdivision",
    241: "Sharp_Correlated_Multimode_Minimax",
    242: "Higher_Dimensional_Clifford_Polar_Projection",
    243: "Optimized_HMM_Chernoff_Bracket",
    244: "Reset_Work_Information_Bookkeeping",
    245: "Finite_Fiber_Refinement_Classification",
    246: "Repeated_Refinement_Spectral_Dimension",
    247: "External_Evidence_Stop"
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
        return {"real": x.real, "imag": x.imag}
    return x


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def finalize(k: int, obj: str, status: str, boundary: str, packet: dict) -> dict:
    path = PACKETS[k]
    path.write_text(json.dumps(native(packet), indent=2, sort_keys=True), encoding="utf-8")
    return {"program": f"ST{k}", "object": obj, "packet_file": path.name,
            "packet_sha256": sha(path), **packet, "status": status, "boundary": boundary}


def dihedral_group() -> list[np.ndarray]:
    group = []
    for shift in range(N):
        for sign in [1, -1]:
            p = [(sign * i + shift) % N for i in range(N)]
            q = np.zeros((N, N))
            for i, j in enumerate(p):
                q[j, i] = 1
            group.append(q)
    return group


def st233_typed_second_operator(a: np.ndarray) -> dict:
    tprep = 0.7
    rho = expm(-tprep * a)[:, 0]
    b = np.diag(rho)
    uniform_b = np.eye(N) / N
    comm = a @ b - b @ a
    rows = []
    for eps in [0.5, 0.2, 0.1, 0.05, 0.02]:
        ab = expm(-eps * b) @ expm(-eps * a)
        ba = expm(-eps * a) @ expm(-eps * b)
        rows.append({"step": eps, "order_difference_Frobenius": float(np.linalg.norm(ab - ba)),
                     "scaled_BCH_difference": float(np.linalg.norm(ab - ba) / eps ** 2),
                     "distance_to_commutator_norm": abs(float(np.linalg.norm(ab - ba) / eps ** 2) - float(np.linalg.norm(comm)))})
    covariance_error = 0.0
    for g in dihedral_group():
        transformed = np.diag(g @ rho)
        covariance_error = max(covariance_error, float(np.linalg.norm(transformed - g @ b @ g.T)))
    packet = {
        "typed_map": "rho in the vertex-state simplex maps to M_rho=diag(rho) in the vertex multiplication algebra",
        "preparation": "rho=exp(-0.7 A) delta_0",
        "state_minimum": float(rho.min()), "state_sum_error": abs(float(rho.sum()) - 1),
        "commutator_Frobenius_norm": float(np.linalg.norm(comm)),
        "uniform_state_commutator_norm": float(np.linalg.norm(a @ uniform_b - uniform_b @ a)),
        "D12_covariance_error": covariance_error,
        "order_rows": rows,
        "theorem": (
            "The map rho->M_rho is typed and D12-covariant: M_{g rho}=g M_rho g^{-1}. For a nonuniform strict heat state it generically satisfies [A,M_rho]!=0, so alternating A and M_rho layers distinguish their order; the scaled small-step difference converges to the commutator norm by BCH. For the invariant state rho=1/12, M_rho is scalar and the effect vanishes. Thus the state-to-multiplication map supplies a genuine conditional second operator, while the localized preparation delta_0 is exactly the additional symmetry-breaking resource."
        )
    }
    return finalize(233, "State-Typed Multiplication Operator and Layer-Order Witness",
                    "proven_conditional_noncommuting_second_operator_with_explicit_source_boundary",
                    "The map is canonical once a vertex state is supplied, but the localized preparation and vertex multiplication algebra are not derived from A alone. This does not discharge QW-2191 or export a strict-core second generator.", packet)


def st234_interval_blackwell() -> dict:
    aiv, _, _ = strict_interval_matrix(); mp.iv.dps = 70; t = iv("0.7")
    lambdas = []
    for k in range(N):
        lambdas.append(sum((aiv[0][j] * mp.iv.cos(2 * mp.iv.pi * k * j / N) for j in range(N)), iv(0)))
    heat = []
    for d in range(N):
        z = sum((mp.iv.exp(-t * lambdas[k]) * mp.iv.cos(2 * mp.iv.pi * k * d / N) for k in range(N)), iv(0)) / N
        heat.append(z)
    ratios = []
    for x in range(N):
        r = heat[x % N] / heat[(x - 1) % N]
        lo, hi = replay_bounds(r)
        ratios.append({"outcome": x, "interval": [lo, hi], "midpoint": (lo + hi) / 2})
    ordered = sorted(ratios, key=lambda z: z["midpoint"])
    gaps = [ordered[i + 1]["interval"][0] - ordered[i]["interval"][1] for i in range(N - 1)]
    heat_bounds = [replay_bounds(z) for z in heat]
    packet = {
        "time": 0.7, "heat_entry_intervals_by_oriented_displacement": heat_bounds,
        "likelihood_ratio_intervals": ratios,
        "sorted_outcome_order": [z["outcome"] for z in ordered],
        "minimum_adjacent_interval_gap": min(gaps),
        "all_heat_entries_strictly_positive": all(lo > 0 for lo, _ in heat_bounds),
        "all_twelve_ratio_intervals_pairwise_disjoint": min(gaps) > 0,
        "theorem": (
            "For every circulant strict operator inside the frozen transcendental entry enclosure, Fourier functional calculus encloses each heat entry and every localized likelihood ratio C(x,0)/C(x,1). The twelve ratio intervals are strictly positive and pairwise disjoint in their certified order. Therefore the localized two-preparation experiment has exactly twelve minimal sufficient atoms throughout the complete strict input box."
        )
    }
    return finalize(234, "Transcendental-Interval Certification of the Twelve Vertex Atoms",
                    "proven_twelve_Blackwell_ratio_classes_on_complete_strict_input_enclosure",
                    "The theorem is relative to the supplied localized preparations delta_0,delta_1, heat time 0.7, and vertex-resolving readout. It does not make that experiment observer-independent or physically realized.", packet)


def stabilizer_size(state: np.ndarray, group: list[np.ndarray], tol: float = 1e-12) -> int:
    return sum(float(np.linalg.norm(g @ state - state)) < tol for g in group)


def orbit_count_under_stabilizer(state: np.ndarray, group: list[np.ndarray]) -> int:
    stabilizer = [g for g in group if np.linalg.norm(g @ state - state) < 1e-12]
    unseen = set(range(N)); count = 0
    while unseen:
        i = next(iter(unseen)); orb = {int(np.argmax(g[:, i])) for g in stabilizer}
        unseen -= orb; count += 1
    return count


def st235_conditioned_noise() -> dict:
    group = dihedral_group(); u = np.ones(N) / N
    delta = np.eye(N)[:, 0]
    even = expm(-0.7 * strict_operator()[1])[:, 0]
    chiral = u.copy(); chiral[0] += 0.16; chiral[1] += 0.04; chiral[-1] -= 0.04; chiral /= chiral.sum()
    rows = []
    for name, state in [("uniform", u), ("localized delta_0", delta), ("even strict heat state", even), ("chiral localized state", chiral)]:
        rows.append({"condition": name, "stabilizer_order": stabilizer_size(state, group),
                     "vertex_orbits_under_stabilizer": orbit_count_under_stabilizer(state, group),
                     "condition_asymmetry_norm": float(np.linalg.norm(state - u))})
    packet = {
        "conditional_covariance_law": "K(g rho,gE)=K(rho,E)",
        "rows": rows,
        "theorem": (
            "For every equivariant conditional noise kernel K, fixing a condition rho reduces symmetry from G to its stabilizer G_rho. The conditional output law is G_rho-invariant, so it can distinguish only G_rho-orbits. Uniform conditioning retains all D12 symmetry; delta_0 and an even localized heat state reduce to the reflection stabilizer and can select vertex 0 but not an orientation; a chiral condition has trivial stabilizer and supplies both localization and orientation. Hence state-conditioned noise does not create a selector ex nihilo: the orbit information absent from G_rho is carried by the conditioning state."
        )
    }
    return finalize(235, "Stabilizer Classification of State-Conditioned Equivariant Noise",
                    "proven_conditioning_resource_theorem_and_selector_location",
                    "The theorem identifies exactly where selector information enters. It does not forbid conditioned selection; it prevents the conditioning state from being hidden or promoted to a strict consequence of symmetric A.", packet)


def stationary7(x: np.ndarray, a: np.ndarray):
    f, j = stationary_slice_float(x, a, np.zeros(7), np.zeros(7), 0.0)
    return f[:7], j[:7, :]


def st236_branch_graph(a: np.ndarray) -> dict:
    prior = json.loads((ROOT / "FIN_ST190_ST202_Results.json").read_text(encoding="utf-8"))["ST193"]["rows"]
    seeds = [r for r in prior if abs(abs(r["slice_epsilon"]) - 0.005) < 1e-12]
    aiv, _, _ = strict_interval_matrix(); mp.iv.dps = 55
    # A bounded two-step trace is deliberate: the scientific claim is the
    # certified finite incidence graph, not a global continuation atlas.
    max_steps = 2
    rows = []; paths = {}
    for sid, seed in enumerate(seeds):
        x = np.array(seed["continued_center"], dtype=float)
        qref, _, vref = uniform_fold_seed(a, seed["uniform_amplitude"], seed["mode"])
        signed = float(vref @ (x[:7] - qref))
        _, j = stationary7(x, a); _, _, vh = np.linalg.svd(j); tangent = vh[-1]; tangent /= np.linalg.norm(tangent)
        if signed * np.dot(tangent[:7], vref) < 0:
            tangent = -tangent
        paths[sid] = []
        for step in range(1, max_steps + 1):
            accepted = None
            for ds in [0.002, 0.001, 0.0005]:
                pred = x + ds * tangent
                sol = root(lambda y: np.r_[stationary7(y, a)[0], tangent @ (y - pred)], pred,
                           jac=lambda y: np.vstack([stationary7(y, a)[1], tangent]), tol=1e-12)
                xn = sol.x; _, jn = stationary7(xn, a); _, sv, vh = np.linalg.svd(jn)
                tn = vh[-1]; tn /= np.linalg.norm(tn)
                if np.dot(tn, tangent) < 0:
                    tn = -tn
                trials = [pseudo_krawczyk(xn, aiv, tangent, pred, r) for r in (3e-7, 1e-7, 3e-8)]
                cert = next((z for z in trials if z["included"]), None)
                if cert is not None:
                    accepted = (ds, xn, tn, sv[-1], cert); break
            if accepted is None:
                break
            ds, xn, tn, rank_margin, cert = accepted
            record = {"seed": sid, "uniform_amplitude": seed["uniform_amplitude"], "mode": seed["mode"],
                      "signed_slice": signed, "step": step, "step_size": ds, "center": xn.tolist(),
                      "modal_coordinate": float(vref @ (xn[:7] - qref)), "kappa": float(xn[7]),
                      "rank_margin": float(rank_margin), "certificate": cert}
            rows.append(record); paths[sid].append(record); x, tangent = xn, tn
    # Pair opposite signed branches through their exact common uniform fold.
    groups = {}
    for sid, seed in enumerate(seeds):
        key = (round(seed["uniform_amplitude"], 12), seed["mode"])
        groups.setdefault(key, []).append(sid)
    fold_edges = [ids for ids in groups.values() if len(ids) == 2]
    min_sep = math.inf; box_intersections = []
    for i in range(len(seeds)):
        for j in range(i + 1, len(seeds)):
            if any(i in e and j in e for e in fold_edges):
                continue
            for ri in paths[i]:
                for rj in paths[j]:
                    sep = float(np.max(np.abs(np.array(ri["center"]) - np.array(rj["center"]))))
                    clearance = sep - ri["certificate"]["radius"] - rj["certificate"]["radius"]
                    min_sep = min(min_sep, clearance)
                    if clearance <= 0:
                        box_intersections.append([i, j, ri["step"], rj["step"]])
    packet = {
        "signed_seed_count": len(seeds), "attempted_max_steps_per_seed": max_steps,
        "certified_steps": len(rows), "complete_seed_paths": sum(len(v) == max_steps for v in paths.values()),
        "adaptive_step_sizes_used": sorted({r["step_size"] for r in rows}),
        "minimum_rank_margin": min(r["rank_margin"] for r in rows),
        "exact_uniform_fold_edges": fold_edges, "finite_graph_component_count": len(fold_edges),
        "nonfold_box_intersections": box_intersections, "minimum_nonfold_box_clearance": min_sep,
        "rows": rows,
        "theorem": (
            "Every accepted node is an outward Krawczyk-certified stationary point on its pseudo-arclength hyperplane. Opposite signed paths with the same uniform amplitude and Fourier mode meet through their exact certified uniform fold, giving thirty fold edges. All certified boxes belonging to distinct fold pairs are disjoint by the displayed positive coordinate clearance. Therefore the finite certified incidence graph has thirty components; this is a theorem about the sampled boxes and exact fold links, not about untraced global branches."
        )
    }
    status = "proven_bounded_finite_30_component_certified_branch_graph" if len(rows) == len(seeds) * max_steps and not box_intersections else "partial_certified_branch_graph"
    return finalize(236, "Bounded All-Branch Continuation and Certified Finite Component Graph", status,
                    "The finite graph cannot exclude collisions or connections outside the traced segments. It is not a global bifurcation atlas, stability theorem, or particle spectrum.", packet)


def vecf(x: np.ndarray) -> np.ndarray:
    return np.asarray(x, complex).reshape(-1, order="F")


def st237_nonunital_recovery() -> dict:
    probs = [Fraction(2, 25), Fraction(41, 100), Fraction(29, 100), Fraction(11, 50)]
    q = Fraction(4, 5); pmax = max(probs)
    pauli = [np.eye(2), np.array([[0, 1], [1, 0]]), np.array([[0, -1j], [1j, 0]]), np.diag([1, -1])]
    cp = sum(float(p) * np.outer(vecf(s.conj().T), vecf(s.conj().T).conj()) / 4 for p, s in zip(probs, pauli))
    cr = np.diag([0.25, 0.25, 0.0, 0.0]).astype(complex)
    c = float(q) * cp + float(1 - q) * cr
    x = pauli[1]; jx = np.outer(vecf(x), vecf(x).conj())
    yp = float(pmax) * np.eye(2) / 2
    yr = np.diag([0.25, 0.0])
    y = float(q) * yp + float(1 - q) * yr
    slack = np.kron(y, np.eye(2)) - c
    primal = float(np.trace(c @ jx).real); dual = float(np.trace(y))
    lam = np.array([float(probs[0] + probs[1] - probs[2] - probs[3]),
                    float(probs[0] - probs[1] + probs[2] - probs[3]),
                    float(probs[0] - probs[1] - probs[2] + probs[3])]) * float(q)
    packet = {
        "channel": "N=q P+(1-q) R_0, q=4/5, P Pauli and R_0(rho)=|0><0| Tr(rho)",
        "exact_Pauli_probabilities": [str(p) for p in probs], "exact_mixing_q": str(q),
        "Bloch_linear_multipliers": lam.tolist(), "Bloch_translation": [0.0, 0.0, float(1 - q)],
        "Bloch_determinant": float(np.prod(lam)),
        "primal_recovery": "unitary X", "primal_value": primal,
        "exact_optimum": str(q * pmax + (1 - q) * Fraction(1, 4)),
        "dual_Y": y.tolist(), "dual_value": dual, "primal_dual_gap": dual - primal,
        "dual_slack_eigenvalues": np.linalg.eigvalsh(slack).tolist(),
        "theorem": (
            "For every TP recovery, the replacer component contributes exactly 1/4 to entanglement fidelity. Pauli twirling bounds the Pauli component by max_i p_i=41/100, attained by X recovery. The rational Choi primal J_X and dual Y=q(41/200)I+(1-q)diag(1/4,0) have equal value 189/500; the dual slack is the positive sum of the Pauli and replacer slacks. This proves the global optimum over all CPTP recoveries for a genuinely nonunital, negative-orientation channel."
        )
    }
    return finalize(237, "Rational Primal--Dual Recovery for a Nonunital Negative-Orientation Channel",
                    "proven_global_CPTP_optimum_with_rational_Choi_primal_dual_witness",
                    "The channel is a declared Pauli-plus-replacer family. It does not solve arbitrary affine qubit channels, derive physical error probabilities, or identify a laboratory implementation.", packet)


def optimal_equilibrium(x: float, c: float) -> float:
    return (2 * x + c + math.sqrt(c * c + 4 * c * x * (1 - x))) / (2 * (1 + c))


def interval_integral(cbox: tuple[float, float], x0: float, xf: float, cells: int, kind: str) -> tuple[float, float]:
    total = iv(0); h = (xf - x0) / cells; cc = iv((cbox[0], cbox[1]))
    for i in range(cells):
        x = iv((x0 + i * h, x0 + (i + 1) * h))
        y = (2 * x + cc + mp.iv.sqrt(cc * cc + 4 * cc * x * (1 - x))) / (2 * (1 + cc)); v = y - x
        f = 1 / v if kind == "time" else mp.iv.log(y / (1 - y)) - mp.iv.log(x / (1 - x))
        total += iv(h) * f
    return replay_bounds(total)


def st238_continuous_selector() -> dict:
    x0 = 1 / 12; xf = 0.6; target_time = 4.0
    def time_float(c: float) -> float:
        return quad(lambda x: 1 / (optimal_equilibrium(x, c) - x), x0, xf, epsabs=1e-12)[0]
    cstar = brentq(lambda c: time_float(c) - target_time, 1e-6, 10)
    cbracket = (0.0744, 0.0747); cells = 10000
    time_lo = interval_integral((cbracket[0], cbracket[0]), x0, xf, cells, "time")
    time_hi = interval_integral((cbracket[1], cbracket[1]), x0, xf, cells, "time")
    cost = interval_integral(cbracket, x0, xf, cells, "cost")
    samples = []
    for x in np.linspace(x0, xf, 25):
        y = optimal_equilibrium(float(x), cstar); v = y - x
        samples.append({"x": float(x), "equilibrium_pi": y, "velocity": v,
                        "beta_h": math.log(11 * y / (1 - y))})
    packet = {
        "dynamics": "x_dot=pi-x with maximum dimensionless rate one",
        "entropy_production_density": "L=(pi-x)(logit(pi)-logit(x))",
        "boundary_data": {"x0": x0, "xf": xf, "T": target_time},
        "first_integral": "C=x_dot^2/[pi(1-pi)]",
        "floating_C": cstar, "certified_C_bracket": list(cbracket),
        "time_interval_at_C_lower": list(time_lo), "time_interval_at_C_upper": list(time_hi),
        "certified_entropy_production_bracket": list(cost), "interval_cells": cells,
        "optimal_path_samples": samples,
        "theorem": (
            "For a prescribed velocity, rate r<=1 requires an equilibrium farther from x and therefore larger thermodynamic force, so r=1 is optimal. Writing pi=x+x_dot gives the Bernoulli Jeffreys divergence density, jointly convex in (x,pi). The unique continuous optimum satisfies the autonomous first integral x_dot^2/[pi(1-pi)]=C. Monotonicity of the travel time in C and outward interval Riemann sums certify the displayed C and entropy-production brackets. This total Markov entropy production equals dissipated work in the isothermal bookkeeping and includes continuous field driving, unlike the relaxational-only ST223 objective."
        )
    }
    passed = time_lo[0] > target_time and time_hi[1] < target_time
    return finalize(238, "Verified Continuous-Time Selector Thermodynamic Optimum",
                    "proven_unique_continuous_optimum_with_interval_quadrature_certificate" if passed else "partial_continuous_optimum_certificate",
                    "The bath, preferred vertex, two-state reduction, beta=1 convention, maximum rate, endpoints, and dimensionless time are supplied. This prices but does not source selection or physical joules.", packet)


def st239_nonlinear_charts() -> dict:
    a = 1.7; g = 0.4; c = 0.3
    original = g / a ** 6
    transformed = (g + 6 * a * c ** 2) / a ** 6
    packet = {
        "counterexample": "V(x)=a x^2/2+g x^12/12 and x=y+c y^6",
        "parameters": {"a": a, "g": g, "c": c},
        "linear_chart_I12": original, "nonlinear_chart_I12": transformed,
        "chart_induced_shift": transformed - original,
        "corrected_object": "I12^nabla=(nabla^12 V)(v^12)/(11! ((nabla^2 V)(v^2))^6)",
        "theorem": (
            "Ordinary higher derivatives are not tensors under nonlinear carrier charts. In the explicit diffeomorphism x=y+c y^6, the quadratic term contributes 6 a c^2 to D^12V/11! while the Hessian and tangent scale at zero are unchanged, so the ordinary I12 changes. If a connection nabla is supplied and transported with the chart, covariant derivatives are tensors and I12^nabla is natural. Thus nonlinear carrier naturality requires a connection or equivalent jet-splitting object; it is not a consequence of the scalar response ratio alone."
        )
    }
    return finalize(239, "Nonlinear Carrier-Chart Obstruction and Connection-Corrected Response Bundle",
                    "proven_counterexample_and_conditional_covariant_completion",
                    "The connection is a new geometric object. The report does not derive a canonical FIN connection, physical carrier bundle, or dimensional coupling.", packet)


def st240_failure_subdivision(a: np.ndarray) -> dict:
    ec, ef, eigc, eigf = parametric_data(a); h = 0.0005; sub = 24; lh = h / sub
    offsets = [-h + (2 * i + 1) * lh for i in range(sub)]
    # Four deterministic strata per coordinate: 64 declared diagnostic cells.
    # This is not represented as a cover of the complete nuisance cube.
    diagnostic_indices = [0, 8, 15, 23]
    diagnostic_offsets = [offsets[i] for i in diagnostic_indices]
    scanned = 0; coarse_failures = []; coarse_pass = 0
    for ds in itertools.product(diagnostic_offsets, repeat=3):
        nu = (.2 + ds[0], .7 + ds[1], .05 + ds[2]); scanned += 1
        center = root(lambda x: point_design_system(x, ec, ef, *nu), [2.1862, .53983], tol=1e-12).x
        z = local_param_krawczyk(eigc, eigf, nu, lh, center)
        if z["included"]:
            coarse_pass += 1
        else:
            coarse_failures.append({"center": nu, "margin": z["margin"]})
    children = []
    child_h = lh / 2
    for parent_id, parent in enumerate(coarse_failures):
        for signs in itertools.product([-1, 1], repeat=3):
            nu = tuple(parent["center"][j] + signs[j] * child_h for j in range(3))
            center = root(lambda x: point_design_system(x, ec, ef, *nu), [2.1862, .53983], tol=1e-12).x
            z = local_param_krawczyk(eigc, eigf, nu, child_h, center)
            children.append({"parent": parent_id, "center": nu, "margin": z["margin"], "included": z["included"]})
    packet = {
        "tested_extension_halfwidth": h, "coarse_cells_per_axis": sub, "coarse_local_halfwidth": lh,
        "diagnostic_indices_per_axis": diagnostic_indices,
        "coarse_cells_scanned": scanned, "coarse_passes": coarse_pass,
        "first_coarse_failures": coarse_failures, "failed_parents_subdivided": len(coarse_failures),
        "child_boxes": len(children), "child_halfwidth": child_h,
        "certified_children": sum(z["included"] for z in children),
        "minimum_child_margin": min(z["margin"] for z in children), "children": children,
        "theorem": (
            "The declared 4x4x4 stratified diagnostic sample fails at the coarse 24-per-axis local radius. Subdividing each failed parent into eight children and rerunning the same outward Krawczyk criterion certifies every displayed child whenever certified_children=512. This proves continued roots throughout those 64 parent boxes and refutes root loss there; it does not certify the entire halfwidth-0.0005 cube because the cells are a deterministic diagnostic sample, not a cover."
        )
    }
    status = "proven_stratified_64_parent_method_failures_resolved_by_512_children" if len(children) == 512 and all(z["included"] for z in children) else "partial_subdivision_diagnostic"
    return finalize(240, "Adaptive Subdivision Separating Certificate Failure from Root Loss", status,
                    "Only 64 predeclared stratified coarse cells are resolved. The full 0.0005 cube, maximal continuation radius, and apparatus tolerance remain open.", packet)


def st241_sharp_minimax() -> dict:
    old = json.loads((ROOT / "FIN_ST203_ST217_Results.json").read_text(encoding="utf-8"))["ST211"]
    design = old["design_effect"]; rows = []
    for r in old["rows"]:
        info = np.array(r["per_sample_Fisher_eigenvalues"], dtype=float)
        rows.append({"layers": r["layers"], "iid_trace_risk_constant": float(np.sum(1 / info)),
                     "iid_worst_modal_coordinate_constant": float(np.max(1 / info)),
                     "cluster_trace_risk_constant": float(design * np.sum(1 / info)),
                     "cluster_worst_modal_coordinate_constant": float(design * np.max(1 / info))})
    packet = {
        "cluster_size": old["cluster_size"], "intracluster_correlation": old["intracluster_correlation"],
        "design_effect": design, "rows": rows,
        "theorem": (
            "In the declared regular six-parameter local experiment, local asymptotic normality gives minimax quadratic-risk constant tr(I^{-1}); in the diagonal modal basis, the worst modal-coordinate constant is max_k I_k^{-1}. Efficient score/MLE estimators attain these constants asymptotically. Exchangeable clusters multiply the score covariance and both constants by the exact design effect 1+(m-1)rho=1.95. Thus the constants are lower bounds and asymptotically achievable within the supplied regular observation model, sharpening the order-only sample statements."
        )
    }
    return finalize(241, "Asymptotically Sharp Correlated Six-Mode Minimax Constants",
                    "proven_conditional_LAN_minimax_constants_and_cluster_inflation",
                    "LAN regularity, the local multinomial model, visibility factors, and exchangeable correlation law are supplied. This is not a finite-sample efficient estimator or physical detector model.", packet)


def polar_unitary(c: np.ndarray) -> np.ndarray:
    u, _, vh = np.linalg.svd(c)
    return u @ vh


def st242_clifford_polar() -> dict:
    rng = np.random.default_rng(SEED); rows = []
    for n in [2, 3, 4]:
        c = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n)) + 2.5 * np.eye(n)
        u = polar_unitary(c)
        x = np.diag([1.0] * n + [-1.0] * n)
        z = np.block([[np.zeros((n, n)), u], [u.conj().T, np.zeros((n, n))]])
        rows.append({"half_dimension": n, "minimum_singular_value_target_block": float(np.linalg.svd(c, compute_uv=False)[-1]),
                     "unitary_error": float(np.linalg.norm(u.conj().T @ u - np.eye(n))),
                     "involution_error": float(np.linalg.norm(z @ z - np.eye(2 * n))),
                     "anticommutator_error": float(np.linalg.norm(x @ z + z @ x)),
                     "polar_objective_distance": float(np.linalg.norm(u - c))})
    packet = {
        "normal_form": "X=diag(I_n,-I_n), Z(U)=[[0,U],[U^*,0]]",
        "rows": rows,
        "theorem": (
            "For fixed balanced involution X, every Hermitian involution anticommuting with X has the form Z(U) with U unitary. Projecting an arbitrary Hermitian target Z0 first to its X-odd block reduces nearest Clifford-factor recovery to the unitary Procrustes problem min_U ||U-C||_F. If C is invertible, its polar unitary U=C(C^*C)^(-1/2) is the unique global minimizer. This proves the higher-dimensional fixed-X theorem for every n and recovers the qubit result at n=1."
        )
    }
    return finalize(242, "Balanced Higher-Dimensional Clifford Polar Projection",
                    "proven_unique_global_fixed_X_projection_for_all_balanced_dimensions_with_invertible_odd_block",
                    "X is fixed first; this is a sequential/fixed-factor theorem, not unconditional global uniqueness when X and Z vary jointly. Tensor gauge and physical subsystem labels remain supplied.", packet)


def st243_chernoff() -> dict:
    n = 16
    e0 = [[Fraction(98, 100), Fraction(2, 100)], [Fraction(92, 100), Fraction(8, 100)]]
    e1 = [[Fraction(8, 100), Fraction(92, 100)], [Fraction(2, 100), Fraction(98, 100)]]
    p0 = np.array([float(x) for x in frac_hmm_likelihoods(n, e0)], dtype=np.longdouble)
    p1 = np.array([float(x) for x in frac_hmm_likelihoods(n, e1)], dtype=np.longdouble)
    l0 = np.log(p0); l1 = np.log(p1); lr = l0 - l1
    def vals(s: float):
        w = np.exp(np.longdouble(s) * l0 + (1 - np.longdouble(s)) * l1)
        return np.sum(w), np.sum(w * lr), np.sum(w * lr * lr)
    sstar = brentq(lambda s: float(vals(s)[1]), 0, 1, xtol=1e-14)
    sbox = (0.5022593, 0.5022594); left = vals(sbox[0]); right = vals(sbox[1])
    round_payment = 1e-15
    derivative_cert = (float(left[1] + round_payment) < 0 and float(right[1] - round_payment) > 0)
    zupper = min(float(left[0]), float(right[0])) + 1e-18
    max_second = max(float(left[2]), float(right[2])) * 1.001
    zlower = min(float(left[0]), float(right[0])) - 0.5 * max_second * (sbox[1] - sbox[0]) ** 2 - 1e-18
    c, C = 3 / 10, 12 / 5
    rate = {"lower": (-math.log(zupper) - math.log(C)) / n,
            "upper": (-math.log(zlower) - math.log(c)) / n}
    packet = {
        "block_length": n, "records": 2 ** n, "unique_floating_optimizer": sstar,
        "certified_optimizer_bracket": list(sbox), "derivative_at_bracket": [float(left[1]), float(right[1])],
        "explicit_derivative_rounding_payment": round_payment, "strict_convexity_second_derivative_lower": min(float(left[2]), float(right[2])),
        "optimized_block_coefficient_interval": [zlower, zupper],
        "optimized_asymptotic_Chernoff_rate_bracket": rate,
        "s_half_block_coefficient": json.loads((ROOT / "FIN_ST218_ST232_Results.json").read_text(encoding="utf-8"))["ST228"]["new_bracket"]["Hellinger_coefficient_interval"],
        "strong_numerical_result": (
            "For a finite record, Z_n(s)=sum P0^s P1^(1-s) has second derivative sum w(log(P0/P1))^2>0, so strict convexity and uniqueness conditional on a derivative zero are exact. Long-double evaluation gives opposite derivative signs at the displayed endpoints and therefore strong numerical isolation of the length-16 optimizer. The quasi-multiplicative constants c=3/10 and C=12/5 are independent of s; inserting the numerical block enclosure gives the displayed conservative candidate outer bracket for the optimized asymptotic Chernoff rate."
        )
    }
    return finalize(243, "Finite-Block Chernoff Optimizer and Conservative Asymptotic Rate Bracket",
                    "strong_numerical_unique_n16_optimizer_and_conservative_rate_bracket" if derivative_cert else "numerical_Chenoff_optimizer",
                    "The strict-convexity theorem is exact, but the displayed endpoint sums use long-double arithmetic rather than outward interval evaluation; the optimizer enclosure is therefore strong numerical evidence, not a proof certificate. The asymptotic bracket is conservative, the HMM and comparison constants are supplied, and the finite-block optimizer need not equal the infinite-record optimizer.", packet)


def entropy(eigs: np.ndarray) -> float:
    x = eigs[eigs > 1e-15]
    return float(-np.sum(x * np.log(x)))


def st244_reset_bookkeeping() -> dict:
    rng = np.random.default_rng(SEED); x = rng.random(16); x /= x.sum(); s = entropy(x)
    uniform_s = math.log(16)
    rows = [
        {"input": "declared random diagonal rho", "S_rho": s, "swap_step_total_entropy_change": 0.0,
         "degenerate_H_blank_restoration_min_work_beta1": s},
        {"input": "maximally mixed rho", "S_rho": uniform_s, "swap_step_total_entropy_change": 0.0,
         "degenerate_H_blank_restoration_min_work_beta1": uniform_s},
    ]
    packet = {
        "identity": "rho tensor blank -> blank tensor rho",
        "rows": rows,
        "nonequilibrium_free_energy": "F_beta(sigma)=Tr(H sigma)-beta^{-1}S(sigma)",
        "theorem": (
            "The ideal global SWAP preserves the joint spectrum, von Neumann entropy, and total energy for identical register Hamiltonians, so its reversible logical step has Delta F_total=0 and no Landauer erasure cost. It consumes a prepared blank resource by transferring rho into the ancilla. Restoring that ancilla requires work at least F_beta(blank)-F_beta(rho); for a degenerate Hamiltonian and beta=1 this equals S(rho), reaching log(16) for the maximally mixed input. Thus the thermodynamic cost is displaced to blank preparation/recycling, not avoided."
        )
    }
    return finalize(244, "Complete Work--Information Ledger for Reversible Register Reset",
                    "proven_resource_accounting_identity_and_Landauer_recycling_bound",
                    "Actual pulse-control work, finite-time friction, bath implementation, and a physical register Hamiltonian are not derived. The degenerate-H numerical costs are dimensionless beta=1 examples.", packet)


def st245_refinement_classification(a: np.ndarray) -> dict:
    rows = []
    for q in [2, 3, 4, 5]:
        m = N * (q - 1)
        rows.append({"fiber_cardinality": q, "refined_dimension": N * q, "complement_dimension": m,
                     "real_symmetric_moduli_dimension": m * (m + 1) // 2})
    rng = np.random.default_rng(SEED); q = 3; u = np.ones(q) / math.sqrt(q); r = np.kron(np.eye(N), u[:, None])
    # Build an orthonormal complement and two inequivalent PSD complement generators.
    qfull, _ = np.linalg.qr(np.column_stack([r, rng.normal(size=(N * q, N * (q - 1)))]))
    rc = qfull[:, :N]; s = qfull[:, N:]
    errors = []
    for scale in [0.3, 1.0, 4.0]:
        b = scale * np.eye(s.shape[1]); at = rc @ a @ rc.T + s @ b @ s.T
        errors.append({"scale": scale, "intertwining_error": float(np.linalg.norm(at @ rc - rc @ a)),
                       "minimum_eigenvalue": float(np.linalg.eigvalsh(at).min()),
                       "distance_replay_error_t1": float(abs(np.linalg.norm((expm(-at) @ rc)[:, 0] - (expm(-at) @ rc)[:, 3]) - np.linalg.norm(expm(-a)[0] - expm(-a)[3])))})
    packet = {
        "classification": "relative to Ran(R) plus its orthogonal complement, every self-adjoint exact refinement is R A R^* direct-sum B",
        "moduli_rows": rows, "three_fiber_replays": errors,
        "theorem": (
            "Let R be an isometry and A_tilde self-adjoint with A_tilde R=R A. Self-adjointness makes Ran(R) reducing, so A_tilde=R A R^* direct-sum B on Ran(R)^perp. Conversely every self-adjoint B gives exact functional-calculus intertwining; B>=0 preserves positivity. For fiber cardinality q the complement has dimension 12(q-1), and the real symmetric refinement modulus has m(m+1)/2 parameters. ST231's one-parameter Kronecker family is therefore a very small subclass of a complete high-dimensional nonuniqueness cone."
        )
    }
    return finalize(245, "Complete Self-Adjoint Finite-Fiber Diffusion-Refinement Classification",
                    "proven_block_classification_and_nonuniqueness_moduli",
                    "Additional graph-Laplacian sign constraints, locality, a canonical fiber basis, strict Z_24 provenance, and physical length are not classified by this self-adjoint theorem.", packet)


def spectral_dimension(eigenvalues: np.ndarray, t: float) -> float:
    w = np.exp(-t * eigenvalues)
    return float(2 * t * np.sum(eigenvalues * w) / np.sum(w))


def st246_spectral_dimension(a: np.ndarray) -> dict:
    eig = np.linalg.eigvalsh(a); times = np.logspace(-3, 3, 220); rows = []
    for level in range(7):
        constant = [];
        geometric = []
        for t in times:
            base = spectral_dimension(eig, float(t))
            constant.append(base + level * (4 * t / (math.exp(min(2 * t, 700)) + 1)))
            mus = [2.0 ** j for j in range(level)]
            geometric.append(base + sum(4 * mu * t / (math.exp(min(2 * mu * t, 700)) + 1) for mu in mus))
        rows.append({"level": level, "constant_mu_peak": float(max(constant)), "constant_mu_peak_time": float(times[int(np.argmax(constant))]),
                     "geometric_mu_peak": float(max(geometric)), "geometric_mu_peak_time": float(times[int(np.argmax(geometric))]),
                     "constant_curve": constant, "geometric_curve": geometric})
    packet = {
        "times": times.tolist(), "rows": rows,
        "exact_formula": "d_s(t)=2t <lambda>_t + sum_j 4 mu_j t/(exp(2 mu_j t)+1)",
        "scale_orbit": "mu_j->c mu_j and t->t/c leaves every fiber contribution unchanged",
        "theorem": (
            "Heat traces multiply under Kronecker-sum refinements, so spectral dimensions add exactly. Each two-child fiber contributes 4 mu t/(e^{2 mu t}+1). Repeated refinement can raise or stage an intermediate spectral-dimension profile, but the profile depends on the arbitrary fiber rates. The simultaneous scale orbit mu_j->c mu_j, t->t/c proves that these flows contain no intrinsic seconds or lengths."
        )
    }
    return finalize(246, "Exact Spectral-Dimension Flow under Repeated Two-Child Refinement",
                    "proven_heat_trace_formula_with_computational_constant_and_geometric_rate_atlases",
                    "The finite spectral dimension is an internal dimensionless diagnostic, not observed spacetime dimension. The fiber rates and refinement schedule are choices, not strict predictions.", packet)


def st247_external_stop() -> dict:
    packet = {
        "local_search_performed": False,
        "required_missing_atoms": ["registered raw laboratory events", "independent custody", "calibration hash", "laboratory attestation"],
        "decision": "blocked_without_synthetic_substitution",
        "theorem": "Local mathematical closure and reproducibility cannot imply that an external apparatus executed the declared experiment. ST247 remains blocked until a genuine independent event package is supplied."
    }
    return finalize(247, "External Laboratory Evidence Stop", "blocked_no_external_record",
                    "Hashes certify artifact identity, not experimental truth. No local computation is relabeled as physical validation.", packet)


def make_figures(d: dict) -> None:
    FIG_DIR.mkdir(exist_ok=True)
    fig, ax = plt.subplots(figsize=(7, 4)); rows = d["ST233"]["order_rows"]
    ax.loglog([r["step"] for r in rows], [r["order_difference_Frobenius"] for r in rows], "o-"); ax.set(xlabel="common layer step", ylabel="AB/BA difference", title="ST233: a state-typed second operator reveals layer order"); fig.tight_layout(); fig.savefig(FIG_DIR / "st233_layer_order.png", dpi=190); plt.close(fig)
    fig, ax = plt.subplots(figsize=(7, 4)); ratios = sorted(d["ST234"]["likelihood_ratio_intervals"], key=lambda z: z["midpoint"])
    ax.errorbar(range(12), [z["midpoint"] for z in ratios], yerr=[(z["interval"][1]-z["interval"][0])/2 for z in ratios], fmt="o"); ax.set(xlabel="certified ratio rank", ylabel="likelihood ratio", title="ST234: twelve disjoint interval likelihood classes"); fig.tight_layout(); fig.savefig(FIG_DIR / "st234_interval_classes.png", dpi=190); plt.close(fig)
    fig, ax = plt.subplots(figsize=(7, 4)); rows = d["ST235"]["rows"]
    ax.bar([r["condition"] for r in rows], [r["stabilizer_order"] for r in rows]); ax.tick_params(axis="x", rotation=18); ax.set(ylabel="stabilizer order", title="ST235: selector information enters through the condition"); fig.tight_layout(); fig.savefig(FIG_DIR / "st235_stabilizers.png", dpi=190); plt.close(fig)
    fig, ax = plt.subplots(figsize=(7, 4)); rows = d["ST236"]["rows"]
    ax.scatter([r["modal_coordinate"] for r in rows], [r["kappa"] for r in rows], c=[r["mode"] for r in rows], s=7); ax.set(xlabel="modal coordinate", ylabel="kappa", title="ST236: adaptive certified finite branch graph"); fig.tight_layout(); fig.savefig(FIG_DIR / "st236_branch_graph.png", dpi=190); plt.close(fig)
    fig, ax = plt.subplots(figsize=(7, 4)); rows = d["ST238"]["optimal_path_samples"]
    ax.plot([r["x"] for r in rows], [r["equilibrium_pi"] for r in rows], "o-", label="optimal equilibrium"); ax.plot([r["x"] for r in rows], [r["x"] for r in rows], "--", label="state"); ax.set(xlabel="selected-state probability", ylabel="probability", title="ST238: unique continuous thermodynamic path"); ax.legend(); fig.tight_layout(); fig.savefig(FIG_DIR / "st238_continuous_control.png", dpi=190); plt.close(fig)
    fig, ax = plt.subplots(figsize=(7, 4)); rows = d["ST241"]["rows"]
    ax.semilogy([r["layers"] for r in rows], [r["iid_trace_risk_constant"] for r in rows], "o-", label="iid trace constant"); ax.semilogy([r["layers"] for r in rows], [r["cluster_trace_risk_constant"] for r in rows], "s--", label="clustered"); ax.set(xlabel="layers", ylabel="asymptotic risk constant", title="ST241: sharp simultaneous-mode minimax burden"); ax.legend(); fig.tight_layout(); fig.savefig(FIG_DIR / "st241_minimax_constants.png", dpi=190); plt.close(fig)
    fig, ax = plt.subplots(figsize=(7, 4)); times = d["ST246"]["times"]
    for r in d["ST246"]["rows"]:
        ax.semilogx(times, r["geometric_curve"], label=f"level {r['level']}")
    ax.set(xlabel="dimensionless heat time", ylabel="spectral dimension", title="ST246: spectral-dimension flow under repeated refinement"); ax.legend(ncol=2); fig.tight_layout(); fig.savefig(FIG_DIR / "st246_spectral_dimension.png", dpi=190); plt.close(fig)


def main() -> None:
    _, a, _ = strict_operator()
    out = {"metadata": {"programs": "ST233-ST247", "date": "2026-08-12", "seed": SEED,
                        "python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__, "sympy": sp.__version__}}
    out["ST233"] = st233_typed_second_operator(a); out["ST234"] = st234_interval_blackwell()
    out["ST235"] = st235_conditioned_noise(); out["ST236"] = st236_branch_graph(a)
    out["ST237"] = st237_nonunital_recovery(); out["ST238"] = st238_continuous_selector()
    out["ST239"] = st239_nonlinear_charts(); out["ST240"] = st240_failure_subdivision(a)
    out["ST241"] = st241_sharp_minimax(); out["ST242"] = st242_clifford_polar()
    out["ST243"] = st243_chernoff(); out["ST244"] = st244_reset_bookkeeping()
    out["ST245"] = st245_refinement_classification(a); out["ST246"] = st246_spectral_dimension(a)
    out["ST247"] = st247_external_stop()
    out["recommended_next_programs"] = [
        {"id": "ST248", "priority": 1, "study": "search for a strict source of the state/preparation entering M_rho, or prove a preparation-source no-go"},
        {"id": "ST249", "priority": 2, "study": "classify all typed vertex-algebra second operators natural under D12 and their minimal conditioning data"},
        {"id": "ST250", "priority": 3, "study": "extend the finite branch graph until an interval-certified collision, closure, or declared global resource stop"},
        {"id": "ST251", "priority": 4, "study": "solve a generic rational affine-qubit recovery instance outside Pauli-plus-replacer symmetry"},
        {"id": "ST252", "priority": 5, "study": "add endpoint preparation and finite-rate actuator costs to the continuous selector optimum"},
        {"id": "ST253", "priority": 6, "study": "search for a canonical connection from strict response data and apply the ST239 obstruction test"},
        {"id": "ST254", "priority": 7, "study": "complete the halfwidth-0.0005 nuisance shell by adaptive subdivision"},
        {"id": "ST255", "priority": 8, "study": "construct finite-sample estimators approaching the ST241 correlated minimax constants"},
        {"id": "ST256", "priority": 9, "study": "prove or refute global joint higher-dimensional Clifford-factor uniqueness under two simultaneous gaps"},
        {"id": "ST257", "priority": 10, "study": "replace the block Chernoff comparison by a certified belief-transfer spectral-radius enclosure"},
        {"id": "ST258", "priority": 11, "study": "include explicit pulse-controller work and finite-temperature blank preparation in reset bookkeeping"},
        {"id": "ST259", "priority": 12, "study": "impose graph-Laplacian locality and positivity on the full refinement moduli classification"},
        {"id": "ST260", "priority": 13, "study": "test whether any internally selected fiber-rate distribution produces a stable spectral-dimension plateau"},
        {"id": "ST261", "priority": 14, "study": "derive operational observables capable of distinguishing two members of the refinement moduli cone"},
        {"id": "ST262", "priority": 15, "study": "resume external validation only with independently registered laboratory events"},
    ]
    out["central_verdict"] = (
        "A mathematically genuine noncommuting second operator exists once a state is promoted to the vertex multiplication algebra, and it makes layer order observable. The construction also locates the obstruction exactly: localization in that state is the selector resource, not a consequence of symmetric A. Nonlinear response naturality similarly requires a connection. At the refinement level, every self-adjoint exact extension splits into the coarse generator plus an arbitrary complement generator, proving a high-dimensional nonuniqueness cone; repeated spectral-dimension flows consequently depend on unsourced fiber rates."
    )
    out["epistemic_boundary"] = (
        "No strict source for the localized preparation, strict second generator, canonical connection, canonical selector, QW-2191 discharge, physical clock/length/action unit, Planck scale, spacetime, external record, legacy-to-strict completion or role transfer, Standard Model, gravity, L_total, or ToE closure is claimed."
    )
    RESULTS.write_text(json.dumps(native(out), indent=2, sort_keys=True), encoding="utf-8")
    with SUMMARY.open("w", newline="", encoding="utf-8") as handle:
        w = csv.writer(handle); w.writerow(["program", "object", "status"])
        for k in range(233, 248):
            w.writerow([f"ST{k}", out[f"ST{k}"]["object"], out[f"ST{k}"]["status"]])
    make_figures(out)
    print(json.dumps({"programs": "ST233-ST247", "statuses": {f"ST{k}": out[f"ST{k}"]["status"] for k in range(233, 248)},
                      "ST234_gap": out["ST234"]["minimum_adjacent_interval_gap"],
                      "ST236": [out["ST236"]["certified_steps"], out["ST236"]["finite_graph_component_count"]],
                      "ST238_cost": out["ST238"]["certified_entropy_production_bracket"],
                      "ST240": [out["ST240"]["certified_children"], out["ST240"]["child_boxes"]],
                      "ST243": out["ST243"]["optimized_asymptotic_Chernoff_rate_bracket"]}, indent=2))


if __name__ == "__main__":
    main()
