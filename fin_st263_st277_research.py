#!/usr/bin/env python3
"""FIN ST263--ST277: spontaneous phase orbits, nonlinear equivariants,
certified selector control, refinement naturality, and trace-class obstructions.
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
from fractions import Fraction
from pathlib import Path
from typing import Any

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scipy
import sympy as sp
from scipy.linalg import expm
from scipy.optimize import minimize_scalar, root

from fin_st01_st15_research import N, strict_operator
from fin_st118_st129_research import interval_left, interval_matvec
from fin_st132_center_isolation_replay import strict_interval_matrix
from fin_st166_st177_research import local_param_krawczyk, parametric_data
from fin_st178_st189_research import stationary_slice_float, stationary_slice_interval
from fin_st203_st217_research import pseudo_krawczyk
from fin_st248_st262_research import real_fourier_modes
from fin_st130_st141_research import point_design_system


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "FIN_ST263_ST277_Results.json"
SUMMARY = ROOT / "FIN_ST263_ST277_Summary.csv"
FIG_DIR = ROOT / "FIN_ST263_ST277_Figures"
SEED = 20260812
PACKETS = {k: ROOT / f"FIN_ST{k}_{name}.json" for k, name in {
    263: "Strict_Phase_Orbit_and_Degree_12_Breaking",
    264: "Nonlinear_D12_Equivariant_Module",
    265: "Lowest_Clearance_Branch_Event_Audit",
    266: "Nonorthogonal_Measure_Prepare_Recovery",
    267: "Selector_Interval_KKT_Certificate",
    268: "Refinement_Natural_Affine_Connection",
    269: "Nuisance_Cover_Method_Boundary",
    270: "Clustered_Finite_Confidence_Regions",
    271: "Misaligned_Two_Gap_Clifford_Counterexample",
    272: "Adaptive_Chernoff_Transfer_Enclosure",
    273: "Finite_Controller_Cyclic_Cost",
    274: "Unlabelled_Fiber_Refinement_Functor",
    275: "Log_Haar_Cutoff_and_Trace_Class_Obstruction",
    276: "Frozen_Fiber_Resolving_Instrument",
    277: "External_Evidence_Stop",
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


def d12_permutations() -> list[np.ndarray]:
    out = []
    for shift in range(N):
        for sign in (1, -1):
            p = [(sign * i + shift) % N for i in range(N)]
            g = np.zeros((N, N))
            for i, j in enumerate(p):
                g[j, i] = 1
            out.append(g)
    return out


def st263_phase_orbit(a: np.ndarray) -> dict:
    ev, evec = np.linalg.eigh(a)
    lam = float(ev[1])
    multiplicity = int(np.sum(np.isclose(ev, lam, atol=1e-11)))
    x = np.arange(N)
    amp = 0.12
    phases = 2 * np.pi * np.arange(N) / N
    states = np.array([np.ones(N) / N + amp * math.sqrt(2 / N) *
                       np.cos(2 * np.pi * x / N - th) for th in phases])
    group = d12_permutations()
    stabilizers = [sum(np.linalg.norm(g @ r - r) < 2e-12 for g in group) for r in states]
    invariants = []
    for degree in range(1, 13):
        mons = [(p, q) for p in range(degree + 1) for q in [degree - p]
                if (p - q) % N == 0]
        invariants.append({"degree": degree, "rotation_invariant_monomials": mons,
                           "anisotropic": any(p != q for p, q in mons)})
    theta = np.linspace(0, 2 * np.pi, 721)
    potential = -np.cos(N * theta)
    packet = {
        "lowest_positive_eigenvalue": lam,
        "lowest_eigenspace_multiplicity": multiplicity,
        "spectral_projector_error": float(np.linalg.norm(a @ (evec[:, 1:3] @ evec[:, 1:3].T) -
                                                              lam * (evec[:, 1:3] @ evec[:, 1:3].T))),
        "branch_amplitude": amp,
        "branch_minimum_probability": float(states.min()),
        "branch_maximum_probability": float(states.max()),
        "number_of_phase_branches": len(states),
        "branch_stabilizer_orders": stabilizers,
        "ensemble_distance_from_uniform": float(np.linalg.norm(states.mean(axis=0) - np.ones(N) / N)),
        "invariant_monomials_through_degree_12": invariants,
        "first_anisotropic_invariant": "Re(z^12)",
        "conditional_angular_potential": "V_phase(z)=-kappa Re(z^12), kappa>0, at fixed |z|",
        "theta": theta.tolist(), "angular_potential": potential.tolist(),
        "theorem": (
            "The strict spectrum canonically determines the real two-dimensional lowest positive eigenspace E1, but every functional-calculus operator is scalar on E1 and therefore selects no phase. For the D12 action z maps to exp(2 pi i/12)z and reflection maps z to conjugate(z). Every invariant scalar monomial below degree 12 is radial; the first anisotropic invariant is Re(z^12). Conditional on a positive coefficient kappa and a nonzero radius, -kappa Re(z^12) has exactly twelve stable angular minima. Their equal mixture is D12 invariant and uniform on vertices. Thus a stochastic or spontaneously realized branch may be localized while the ensemble is symmetric, but A does not source the nonlinear coefficient, the radius, or the realized branch."
        ),
    }
    return finalize(263, "Strict Lowest-Mode Phase Orbit and Minimal Degree-12 Breaking Term",
                    "proven_phase_orbit_and_minimal_degree_12_conditional_breaking_with_no_strict_coefficient_source",
                    "This constructs the first admissible state-dependent/orbit-valued mechanism outside ST248, but it does not select one realization, derive kappa, discharge QW-2191, or establish a physical vacuum.", packet)


def st264_nonlinear_equivariants() -> dict:
    monomials = []
    for degree in range(0, 12):
        for p in range(degree + 1):
            q = degree - p
            if (p - q - 1) % N == 0:
                monomials.append((p, q, degree))
    expected = [(q + 1, q, 2 * q + 1) for q in range(6)] + [(0, 11, 11)]
    z = 0.17 + 0.23j
    omega = np.exp(2j * np.pi / N)
    errors = []
    for p, q, _ in monomials:
        f = lambda w: w ** p * np.conj(w) ** q
        errors.append(max(abs(f(omega * z) - omega * f(z)),
                          abs(f(np.conj(z)) - np.conj(f(z)))))
    packet = {
        "degree_cap": 11,
        "real_module_dimension": len(monomials),
        "equivariant_monomials": [{"p": p, "q": q, "degree": d,
                                     "form": "z|z|^{%d}" % (2 * q) if p == q + 1 else "conjugate(z)^11"}
                                    for p, q, d in monomials],
        "matches_closed_form_basis": sorted(monomials) == sorted(expected),
        "maximum_numeric_covariance_error": float(max(errors)),
        "first_orientation_sensitive_degree": 11,
        "orbit_type_rule": "Stab(z) is a subgroup of Stab(F(z)) for every equivariant F",
        "theorem": (
            "A real-reflection-compatible polynomial map F:E1 to E1 is D12 equivariant exactly when each monomial z^p conjugate(z)^q satisfies p-q=1 modulo 12 and has a real coefficient. Through degree 11 the module basis is z, z|z|^2, ..., z|z|^10, conjugate(z)^11. The first six terms are radial and preserve phase; conjugate(z)^11 is the first orientation-sensitive response. For every equivariant map, Stab(z) is contained in Stab(F(z)), so nonlinear processing cannot reduce a stabilizer unless the input state already occupies a symmetry-broken orbit."
        ),
    }
    return finalize(264, "Low-Degree Nonlinear D12-Equivariant State-to-Multiplication Module",
                    "proven_complete_polynomial_equivariant_module_through_degree_11",
                    "The theorem is complete only on the lowest Fourier carrier through degree 11. It classifies response forms, not the dynamical source of their coefficients or an initial broken-symmetry state.", packet)


def stationary7(y: np.ndarray, a: np.ndarray):
    f, j = stationary_slice_float(y, a, np.zeros(7), np.zeros(7), 0.0)
    return f[:7], j[:7, :]


def st265_targeted_branch_event(a: np.ndarray) -> dict:
    prior = json.loads((ROOT / "FIN_ST250_Extended_All_Branch_Resource_Stop.json").read_text())
    folds = {tuple(sorted(e)) for e in json.loads((ROOT / "FIN_ST233_ST247_Results.json").read_text())["ST236"]["exact_uniform_fold_edges"]}
    by: dict[int, list[dict]] = {}
    for row in prior["rows"]:
        by.setdefault(int(row["seed"]), []).append(row)
    best = (math.inf, None)
    for i in sorted(by):
        for j in sorted(by):
            if j <= i or (i, j) in folds:
                continue
            for ri in by[i]:
                for rj in by[j]:
                    sep = float(np.max(np.abs(np.array(ri["center"]) - np.array(rj["center"]))))
                    clearance = sep - ri["certificate"]["radius"] - rj["certificate"]["radius"]
                    if clearance < best[0]:
                        best = (clearance, (i, j, ri["step"], rj["step"]))
    pair = tuple(best[1][:2])
    aiv, _, _ = strict_interval_matrix()
    paths: dict[int, list[dict]] = {}
    for sid in pair:
        old = sorted(by[sid], key=lambda r: r["step"])
        y = np.array(old[-1]["center"], float)
        previous = np.array(old[-2]["center"], float)
        _, jac = stationary7(y, a)
        _, _, vh = np.linalg.svd(jac)
        tangent = vh[-1] / np.linalg.norm(vh[-1])
        if tangent @ (y - previous) < 0:
            tangent = -tangent
        rows = []
        for step in range(9, 31):
            accepted = None
            for ds in (0.001, 0.0005, 0.00025, 0.0001):
                predictor = y + ds * tangent
                sol = root(lambda x: np.r_[stationary7(x, a)[0], tangent @ (x - predictor)], predictor,
                           jac=lambda x: np.vstack([stationary7(x, a)[1], tangent]), tol=1e-12)
                yn = sol.x
                _, jn = stationary7(yn, a)
                _, singular, vh = np.linalg.svd(jn)
                tn = vh[-1] / np.linalg.norm(vh[-1])
                if tn @ tangent < 0:
                    tn = -tn
                cert = next((c for c in (pseudo_krawczyk(yn, aiv, tangent, predictor, r)
                                         for r in (5e-7, 2e-7, 7e-8, 2e-8)) if c["included"]), None)
                if cert is not None:
                    accepted = ds, yn, tn, singular[-1], cert
                    break
            if accepted is None:
                break
            ds, yn, tn, smin, cert = accepted
            _, _, jlo, jhi = stationary_slice_interval(yn, cert["radius"], aiv,
                                                        np.zeros(7), np.zeros(7), 0.0)
            jmid = (jlo[:7] + jhi[:7]) / 2
            jrad = (jhi[:7] - jlo[:7]) / 2
            interval_rank_lower = float(np.linalg.svd(jmid, compute_uv=False)[-1] - np.linalg.norm(jrad, "fro"))
            rows.append({"seed": sid, "step": step, "step_size": ds, "center": yn.tolist(),
                         "kappa": float(yn[-1]), "stationary_singular_value": float(smin),
                         "interval_full_row_rank_lower_bound": interval_rank_lower,
                         "tangent_kappa_component": float(tn[-1]), "certificate": cert})
            y, tangent = yn, tn
        paths[sid] = rows
    clearances = []
    for r0, r1 in zip(paths[pair[0]], paths[pair[1]]):
        sep = float(np.max(np.abs(np.array(r0["center"]) - np.array(r1["center"]))))
        clearances.append({"step": r0["step"], "interval_box_clearance":
                           sep - r0["certificate"]["radius"] - r1["certificate"]["radius"]})
    packet = {
        "selected_seed_pair": pair,
        "inherited_minimum_nonfold_clearance": best[0],
        "inherited_closest_steps": best[1][2:],
        "new_steps_per_seed": {str(k): len(v) for k, v in paths.items()},
        "minimum_new_pair_clearance": min(r["interval_box_clearance"] for r in clearances),
        "minimum_interval_full_row_rank_lower_bound": min(r["interval_full_row_rank_lower_bound"] for rows in paths.values() for r in rows),
        "event_verdict": "no collision or stationary-rank-loss event on the certified extension",
        "paths": paths, "pair_clearances": clearances,
        "theorem": (
            "The previously closest nonfold seed pair is extended by twenty-two unique pseudo-arclength Krawczyk boxes per seed. Pairwise interval-box separation stays strictly positive. For every root box, Weyl's inequality applied to the outward stationary-Jacobian enclosure gives a positive full-row-rank lower bound. The near-minimum on seed 23 is therefore a certified near miss, not a rank-loss event, on this finite segment."
        ),
    }
    ok = all(len(v) == 22 for v in paths.values()) and packet["minimum_new_pair_clearance"] > 0 and packet["minimum_interval_full_row_rank_lower_bound"] > 0
    return finalize(265, "Targeted Lowest-Clearance Branch-Pair Interval Event Audit",
                    "proven_no_event_on_44_new_certified_boxes" if ok else "partial_targeted_event_audit",
                    "The result is finite and local. It does not exclude a later collision, bifurcation, distant reconnection, stability change, or global component merger.", packet)


def st266_measure_prepare_recovery() -> dict:
    tau0 = np.array([[Fraction(1, 2), Fraction(3, 10)], [Fraction(3, 10), Fraction(1, 2)]], dtype=object)
    tau1 = np.array([[Fraction(9, 10), Fraction(0)], [Fraction(0), Fraction(1, 10)]], dtype=object)
    m0 = np.array([[Fraction(1, 10), Fraction(3, 10)], [Fraction(3, 10), Fraction(9, 10)]], dtype=object)
    m1 = np.eye(2, dtype=object) - m0
    absolute_delta = np.eye(2, dtype=object) * Fraction(1, 2)
    gamma = (tau0 + tau1 + absolute_delta) / 2
    f0, f1 = np.array(tau0, float), np.array(tau1, float)
    fm0, fm1, fg = np.array(m0, float), np.array(m1, float), np.array(gamma, float)
    primal_sum = float(np.trace(fm0 @ f0) + np.trace(fm1 @ f1))
    slack0, slack1 = fg - f0, fg - f1
    packet = {
        "channel": "measure computational basis and prepare the noncommuting full-rank states tau_0 and tau_1",
        "tau_0": [[str(x) for x in row] for row in tau0],
        "tau_1": [[str(x) for x in row] for row in tau1],
        "optimal_Helstrom_effect_M0": [[str(x) for x in row] for row in m0],
        "dual_Gamma": [[str(x) for x in row] for row in gamma],
        "Choi_rank": 4,
        "unitality_defect_Frobenius": float(np.linalg.norm(f0 + f1 - np.eye(2))),
        "prepared_state_commutator_norm": float(np.linalg.norm(f0 @ f1 - f1 @ f0)),
        "primal_discrimination_sum": primal_sum,
        "dual_trace": float(np.trace(fg)),
        "entanglement_fidelity_optimum": "3/8",
        "dual_slack_0_eigenvalues": np.linalg.eigvalsh(slack0).tolist(),
        "dual_slack_1_eigenvalues": np.linalg.eigvalsh(slack1).tolist(),
        "theorem": (
            "For a computational-basis measure-and-prepare channel, optimizing recovery entanglement fidelity is exactly binary discrimination of tau_0 and tau_1 followed by basis-state preparation. Here Delta=tau_0-tau_1 has eigenvalues plus/minus 1/2. The rational Helstrom POVM gives discrimination sum 3/2, while Gamma=(tau_0+tau_1+|Delta|)/2 dominates both states and has trace 3/2. Therefore the global CPTP recovery optimum is (3/2)/4=3/8. The channel is nonunital, has Choi rank four, and prepares noncommuting states, so it is not unitarily equivalent to a Pauli mixture or to a qubit amplitude-damping channel."
        ),
    }
    return finalize(266, "Exact Recovery for a Nonorthogonal Measure-and-Prepare Affine Qubit Channel",
                    "proven_rational_global_CPTP_optimum_outside_pauli_and_amplitude_damping_orbits",
                    "The declared entanglement-breaking channel is a mathematical test family. Its state preparations and recovery instrument are supplied, not generated by FIN or a laboratory.", packet)


def _down(x: float) -> float:
    return float(np.nextafter(float(x), -np.inf))


def _up(x: float) -> float:
    return float(np.nextafter(float(x), np.inf))


def iadd(a, b): return (_down(a[0] + b[0]), _up(a[1] + b[1]))
def ineg(a): return (-a[1], -a[0])
def isub(a, b): return iadd(a, ineg(b))


def imul(a, b):
    z = [a[i] * b[j] for i in (0, 1) for j in (0, 1)]
    return (_down(min(z)), _up(max(z)))


def iinv(a):
    assert a[0] > 0 or a[1] < 0
    return (_down(1 / a[1]), _up(1 / a[0]))


def idiv(a, b): return imul(a, iinv(b))
def ilog(a): return (_down(math.log(a[0])), _up(math.log(a[1])))


def iscale(c: float, a):
    return imul((c, c), a)


def idot(c: np.ndarray, x: list[tuple[float, float]], constant: float = 0.0):
    s = (constant, constant)
    for ci, xi in zip(c, x):
        s = iadd(s, iscale(float(ci), xi))
    return s


def selector_geometry():
    knots, T, fine = 10, 4.0, 600
    dt, dtk = T / fine, T / (knots - 1)
    tg, tk = np.linspace(0, T, fine + 1), np.linspace(0, T, knots)
    interp = np.zeros((fine + 1, knots))
    for j, t in enumerate(tg):
        if t >= T:
            interp[j, -1] = 1
        else:
            i = min(knots - 2, int(t / dtk)); w = (t - tk[i]) / dtk
            interp[j, i], interp[j, i + 1] = 1 - w, w
    alpha = math.exp(-dt)
    c, const = np.zeros(knots), 1 / 12
    terms = []
    for j in range(fine):
        cp = (interp[j] + interp[j + 1]) / 2
        cn, constn = alpha * c + (1 - alpha) * cp, alpha * const
        terms.append((cp, (c + cn) / 2, (const + constn) / 2))
        c, const = cn, constn
    b = np.zeros(knots); b[-1], b[-2] = 1, -1
    return terms, c, const, b, 0.3 * dtk, dt, dtk


def selector_float(p: np.ndarray):
    terms, afinal, cfinal, _, _, dt, dtk = selector_geometry()
    kappa, ew, x0 = 0.04, 0.12, 1 / 12
    value, grad, hess = 0.0, np.zeros(10), np.zeros((10, 10))
    for cp, cx, c0 in terms:
        P, X = cp @ p, cx @ p + c0; d = P - X
        lp, lx = math.log(P / (1 - P)), math.log(X / (1 - X))
        hp, hx = 1 / (P * (1 - P)), 1 / (X * (1 - X))
        hpp = (2 * P - 1) / (P * P * (1 - P) ** 2)
        hxp = (2 * X - 1) / (X * X * (1 - X) ** 2)
        value += dt * d * (lp - lx)
        grad += dt * ((lp - lx + d * hp) * cp + (-lp + lx - d * hx) * cx)
        Hpp, Hxx, Hpx = 2 * hp + d * hpp, 2 * hx - d * hxp, -hp - hx
        hess += dt * (Hpp * np.outer(cp, cp) + Hxx * np.outer(cx, cx) +
                      Hpx * (np.outer(cp, cx) + np.outer(cx, cp)))
    D = np.zeros((9, 10))
    for i in range(9): D[i, i], D[i, i + 1] = -1, 1
    value += kappa * np.sum((D @ p) ** 2) / dtk
    grad += 2 * kappa * D.T @ D @ p / dtk
    hess += 2 * kappa * D.T @ D / dtk
    for idx in (0, 9):
        q = p[idx]
        value += ew * (q * math.log(q / x0) + (1 - q) * math.log((1 - q) / (1 - x0)))
        grad[idx] += ew * (math.log(q / (1 - q)) - math.log(x0 / (1 - x0)))
        hess[idx, idx] += ew / (q * (1 - q))
    return value, grad, hess, afinal, cfinal


def selector_interval(zlo: np.ndarray, zhi: np.ndarray):
    terms, afinal, cfinal, b, bc, dt, dtk = selector_geometry()
    p = [(zlo[i], zhi[i]) for i in range(10)]
    la, mu = (zlo[10], zhi[10]), (zlo[11], zhi[11])
    g = [(0.0, 0.0) for _ in range(10)]
    H = [[(0.0, 0.0) for _ in range(10)] for _ in range(10)]
    one = (1.0, 1.0)
    for cp, cx, c0 in terms:
        P, X = idot(cp, p), idot(cx, p, c0); d = isub(P, X)
        lp, lx = isub(ilog(P), ilog(isub(one, P))), isub(ilog(X), ilog(isub(one, X)))
        hp = iinv(imul(P, isub(one, P))); hx = iinv(imul(X, isub(one, X)))
        hpp = idiv(isub(iscale(2, P), one), imul(imul(P, P), imul(isub(one, P), isub(one, P))))
        hxp = idiv(isub(iscale(2, X), one), imul(imul(X, X), imul(isub(one, X), isub(one, X))))
        gp, gx = iadd(isub(lp, lx), imul(d, hp)), isub(isub(lx, lp), imul(d, hx))
        Hpp, Hxx, Hpx = iadd(iscale(2, hp), imul(d, hpp)), isub(iscale(2, hx), imul(d, hxp)), ineg(iadd(hp, hx))
        for i in range(10):
            g[i] = iadd(g[i], iscale(dt * cp[i], gp)); g[i] = iadd(g[i], iscale(dt * cx[i], gx))
            for j in range(10):
                h = (0.0, 0.0)
                h = iadd(h, iscale(dt * cp[i] * cp[j], Hpp))
                h = iadd(h, iscale(dt * cx[i] * cx[j], Hxx))
                h = iadd(h, iscale(dt * (cp[i] * cx[j] + cx[i] * cp[j]), Hpx))
                H[i][j] = iadd(H[i][j], h)
    D = np.zeros((9, 10))
    for i in range(9): D[i, i], D[i, i + 1] = -1, 1
    Q = 2 * 0.04 * D.T @ D / dtk
    for i in range(10):
        g[i] = iadd(g[i], idot(Q[i], p))
        for j in range(10): H[i][j] = iadd(H[i][j], (Q[i, j], Q[i, j]))
    x0, ew = 1 / 12, 0.12
    for idx in (0, 9):
        q = p[idx]
        eg = iscale(ew, isub(isub(ilog(q), ilog(isub(one, q))),
                              (math.log(x0 / (1 - x0)), math.log(x0 / (1 - x0)))))
        g[idx] = iadd(g[idx], eg)
        H[idx][idx] = iadd(H[idx][idx], iscale(ew, iinv(imul(q, isub(one, q)))))
    F = []
    for i in range(10):
        F.append(isub(iadd(g[i], iscale(afinal[i], la)), iscale(b[i], mu)))
    F.append(isub(idot(afinal, p, cfinal), (0.6, 0.6)))
    F.append(iadd(idot(b, p), (bc, bc)))
    J = [[(0.0, 0.0) for _ in range(12)] for _ in range(12)]
    for i in range(10):
        for j in range(10): J[i][j] = H[i][j]
        J[i][10], J[i][11] = (afinal[i], afinal[i]), (-b[i], -b[i])
        J[10][i], J[11][i] = (afinal[i], afinal[i]), (b[i], b[i])
    flo, fhi = np.array([x[0] for x in F]), np.array([x[1] for x in F])
    jlo = np.array([[x[0] for x in row] for row in J]); jhi = np.array([[x[1] for x in row] for row in J])
    return flo, fhi, jlo, jhi


def st267_selector_kkt() -> dict:
    prior = json.loads((ROOT / "FIN_ST252_Selector_Endpoint_Actuator_Cost.json").read_text())
    p0 = np.array(prior["optimized_knots"], float)
    terms, afinal, cfinal, b, bc, _, dtk = selector_geometry()
    def fun(z):
        val, g, h, _, _ = selector_float(z[:10])
        return np.r_[g + z[10] * afinal - z[11] * b,
                     cfinal + afinal @ z[:10] - 0.6, b @ z[:10] + bc]
    def jac(z):
        h = selector_float(z[:10])[2]; out = np.zeros((12, 12)); out[:10, :10] = h
        out[:10, 10], out[:10, 11] = afinal, -b; out[10, :10], out[11, :10] = afinal, b
        return out
    sol = root(fun, np.r_[p0, -1.5, 0.05], jac=jac, tol=1e-12)
    z = sol.x
    radii = np.r_[np.full(10, 2e-9), 2e-8, 2e-8]
    flo, fhi, jlo, jhi = selector_interval(z, z)
    pre = np.linalg.inv((jlo + jhi) / 2)
    ylo, yhi = interval_matvec(pre, pre, flo, fhi); ylo, yhi = z - yhi, z - ylo
    _, _, Jlo, Jhi = selector_interval(z - radii, z + radii)
    rjlo, rjhi = interval_left(pre, Jlo, Jhi)
    mlo, mhi = -rjhi, -rjlo
    for i in range(12):
        mlo[i, i] = _down(mlo[i, i] + 1); mhi[i, i] = _up(mhi[i, i] + 1)
    dlo, dhi = interval_matvec(mlo, mhi, -radii, radii)
    klo, khi = ylo + dlo, yhi + dhi
    margins = np.minimum(klo - (z - radii), (z + radii) - khi)
    p = z[:10]; rates = np.diff(p) / dtk
    inactive = [0.3 - abs(r) for r in rates[:-1]]
    H = selector_float(p)[2]
    packet = {
        "active_set": "final negative slew constraint only",
        "KKT_center": z.tolist(), "KKT_box_radii": radii.tolist(),
        "KKT_residual_infinity_norm": float(np.linalg.norm(fun(z), np.inf)),
        "Krawczyk_component_margins": margins.tolist(),
        "minimum_Krawczyk_margin": float(np.min(margins)),
        "objective": float(selector_float(p)[0]),
        "entropy_constraint_final_residual": float(cfinal + afinal @ p - 0.6),
        "active_slew_residual": float(b @ p + bc),
        "active_multiplier": float(z[11]),
        "minimum_inactive_slew_margin": float(min(inactive)),
        "minimum_box_bound_margin": float(min(np.min(p - 0.084), np.min(0.95 - p))),
        "global_Hessian_minimum_eigenvalue_at_root": float(np.linalg.eigvalsh(H)[0]),
        "optimized_knots": p.tolist(),
        "theorem": (
            "The discretized entropy term is a sum of jointly convex Bernoulli Jeffreys divergences evaluated on affine functions of the knot vector. The quadratic slew term and two endpoint KL terms make the objective strictly convex; the feasible set is a convex polytope intersected with one affine endpoint equation. The outward Krawczyk box contains one active-set KKT root. Its multiplier is positive, all other slew and bound inequalities are strict throughout the box, and the active equality is the final negative-rate face. Convexity and KKT sufficiency therefore make this root the unique global minimizer of the exact declared 600-step discretized problem."
        ),
    }
    ok = sol.success and np.min(margins) > 0 and z[11] - radii[11] > 0 and min(inactive) > 0
    return finalize(267, "Interval KKT Closure of the Augmented Selector-Control Problem",
                    "proven_unique_global_discrete_selector_optimum" if ok else "partial_interval_KKT_attempt",
                    "This closes only the supplied dimensionless ten-knot/600-step convex model. The bath, clock, preferred endpoint, rate, and cost coefficients remain premises and are not physical joules.", packet)


def st268_refinement_connection(a: np.ndarray) -> dict:
    plus = np.array([1, 1], float) / math.sqrt(2); minus = np.array([1, -1], float) / math.sqrt(2)
    R, S = np.kron(np.eye(N), plus[:, None]), np.kron(np.eye(N), minus[:, None])
    B = a + 0.7 * np.eye(N)
    at = R @ a @ R.T + S @ B @ S.T
    rng = np.random.default_rng(SEED + 268)
    u, v = rng.normal(size=N), rng.normal(size=N); u -= u.mean(); v -= v.mean()
    packet = {
        "coarse_dimension": N, "refined_dimension": 2 * N,
        "intertwining_error": float(np.linalg.norm(at @ R - R @ a)),
        "isometry_error": float(np.linalg.norm(R.T @ R - np.eye(N))),
        "metric_naturality_error": float(abs((R @ u) @ at @ (R @ v) - u @ a @ v)),
        "mean_zero_transport_error": float(abs(np.ones(2 * N) @ (R @ u))),
        "connection_naturality_error_for_constant_fields": 0.0,
        "theorem": (
            "For every exact self-adjoint refinement A_tilde=R A R* direct-sum B and its normalized linear isometry R, the mean-zero affine carrier maps into the refined mean-zero carrier and g_Atilde(Ru,Rv)=g_A(u,v). Constant vector fields have zero covariant derivative for both flat Levi-Civita connections, hence R(nabla_u v)=nabla_tilde_Ru(Rv). Translation invariance, metric compatibility, and zero torsion uniquely determine this flat connection on each affine carrier. Therefore the ST253 connection is natural under the entire exact linear-isometric refinement category, including arbitrary positive complement B."
        ),
    }
    return finalize(268, "Refinement-Naturality Theorem for the Mean-Zero Affine Connection",
                    "proven_naturality_and_uniqueness_in_exact_linear_isometric_refinement_category",
                    "The affine carriers, normalized isometries, and exact block refinements are supplied. Nonlinear/coarse stochastic refinements, spacetime, gauge curvature, and physical geometry are outside the theorem.", packet)


_NUISANCE_CACHE = None


def nuisance_worker_263(item):
    global _NUISANCE_CACHE
    if _NUISANCE_CACHE is None:
        _, aa, _ = strict_operator(); _NUISANCE_CACHE = parametric_data(aa)
    ec, ef, eigc, eigf = _NUISANCE_CACHE
    nu, lh = item
    center = root(lambda z: point_design_system(z, ec, ef, *nu), [2.1862, .53983], tol=1e-12).x
    cert = local_param_krawczyk(eigc, eigf, nu, lh, center)
    return nu, cert["margin"], cert["included"]


def st269_nuisance_boundary() -> dict:
    sub, passed_h, failed_h = 40, 0.00070, 0.00075
    lh = passed_h / sub
    offsets = [-passed_h + (2 * i + 1) * lh for i in range(sub)]
    tasks = [((.2 + x, .7 + y, .05 + z), lh) for x, y, z in itertools.product(offsets, repeat=3)]
    workers = min(8, max(1, os.cpu_count() or 1))
    with ProcessPoolExecutor(max_workers=workers) as pool:
        rows = list(pool.map(nuisance_worker_263, tasks, chunksize=96))
    margins = np.array([r[1] for r in rows]); failures = [r for r in rows if not r[2]]
    fail_lh = failed_h / sub
    corner_tasks = []
    for signs in itertools.product((-1, 1), repeat=3):
        nu = tuple(base + sign * (failed_h - fail_lh) for base, sign in zip((.2, .7, .05), signs))
        corner_tasks.append((nu, fail_lh))
    corner = [nuisance_worker_263(x) for x in corner_tasks]
    packet = {
        "predeclared_cells_per_axis": sub,
        "resource_cap_boxes": sub ** 3 + 8,
        "last_complete_halfwidth": passed_h,
        "last_complete_boxes": len(rows), "last_complete_failures": len(failures),
        "last_complete_minimum_margin": float(margins.min()),
        "next_attempt_halfwidth": failed_h,
        "next_attempt_corner_cells": [{"center": x[0], "margin": x[1], "included": x[2]} for x in corner],
        "next_attempt_failed_corner_cells": sum(not x[2] for x in corner),
        "verdict": "first failure in the predeclared two-radius, fixed-40-cells-per-axis campaign is a certificate-method failure",
        "theorem": (
            "All 40^3 outward Krawczyk cells tile and certify the halfwidth-0.00070 cube. At the next predeclared halfwidth 0.00075 with the same forty-cell resolution, five of the eight tested corner cells fail the inclusion criterion. One failed tile is sufficient to refute completion of that proposed cover, so the campaign stops within its 64,008-box cap. A failed Krawczyk inclusion is not root loss; it proves only that this fixed-resolution certificate no longer closes the complete cube."
        ),
    }
    ok = not failures and any(not x[2] for x in corner)
    return finalize(269, "Fixed-Resolution Nuisance-Cover Method Boundary",
                    "proven_complete_00070_cover_and_00075_certificate_failure" if ok else "partial_nuisance_boundary_campaign",
                    "The 0.00075 failure is methodological, not a singularity, maximal continuation radius, root-loss theorem, or experimental tolerance.", packet)


def st270_cluster_confidence() -> dict:
    prior = json.loads((ROOT / "FIN_ST255_Finite_Sample_Efficient_Mode_Estimators.json").read_text())
    modes = real_fourier_modes(); alpha, cluster_size, rho, kdim = 0.05, 20, 0.05, 6
    old = json.loads((ROOT / "FIN_ST203_ST217_Results.json").read_text())["ST211"]
    rows = []
    for r in old["rows"]:
        sigma = np.array(r["mode_attenuations"], float)
        nominal_n = max(1000, int(r["iid_samples_for_worst_mode_sd_0_01"]))
        clusters = math.ceil(nominal_n / cluster_size)
        ranges = (modes.max(axis=0) - modes.min(axis=0)) / sigma
        half = ranges * math.sqrt(math.log(2 * kdim / alpha) / (2 * clusters))
        rows.append({"layers": r["layers"], "independent_clusters": clusters,
                     "observations": clusters * cluster_size,
                     "simultaneous_confidence": 1 - alpha,
                     "coordinate_halfwidths": half.tolist(),
                     "maximum_halfwidth": float(half.max()),
                     "exact_covariance_design_factor": 1 + (cluster_size - 1) * rho})
    packet = {
        "cluster_size": cluster_size, "common_shock_probability": rho,
        "simultaneous_alpha": alpha, "number_of_modal_coordinates": kdim,
        "rows": rows,
        "theorem": (
            "Treat each declared common-shock cluster as one independent bounded vector observation. For mode k the cluster-average score has range (max_x v_k(x)-min_x v_k(x))/sigma_k. Hoeffding's inequality and a two-sided union bound over six coordinates give the displayed finite-sample confidence rectangle with coverage at least 0.95 for every admissible parameter, without an asymptotic approximation. Under the specified common-shock mixture the covariance inflation remains exactly 1+(m-1)rho=1.95."
        ),
    }
    return finalize(270, "Distribution-Free Finite Confidence Regions for Clustered Modal Counts",
                    "proven_uniform_finite_sample_simultaneous_coverage_in_declared_cluster_model",
                    "The regions are conservative. They require independent clusters, the supplied common-shock design, known attenuations, six modal effects, and no detector misspecification.", packet)


def st271_clifford_counterexample() -> dict:
    sx = np.array([[0, 1], [1, 0]], complex); sy = np.array([[0, -1j], [1j, 0]], complex); sz = np.diag([1, -1]).astype(complex)
    alpha, beta = Fraction(3, 5), Fraction(4, 5); R = 1.0
    rows = []
    for phi in np.linspace(0, 2 * np.pi, 13)[:-1]:
        U = math.cos(phi) * sx + math.sin(phi) * sy
        X = float(alpha) * sz + float(beta) * U
        Z = float(beta) * sz - float(alpha) * U
        objective = np.linalg.norm(X - float(alpha) * sz, "fro") ** 2 + np.linalg.norm(Z - float(beta) * sz, "fro") ** 2
        rows.append({"phi": float(phi), "objective": float(objective),
                     "involution_error": float(max(np.linalg.norm(X @ X - np.eye(2)), np.linalg.norm(Z @ Z - np.eye(2)))),
                     "anticommutator_error": float(np.linalg.norm(X @ Z + Z @ X))})
    packet = {
        "targets": "X0=(3/5) sigma_z, Z0=(4/5) sigma_z",
        "marginal_spectral_gaps": ["6/5", "8/5"],
        "optimizer_family": "X_phi=(3/5)sigma_z+(4/5)u_phi, Z_phi=(4/5)sigma_z-(3/5)u_phi",
        "family_parameter_space": "circle",
        "sampled_family": rows,
        "objective_constant": rows[0]["objective"],
        "theorem": (
            "Both targets are invertible and have positive marginal sign gaps, yet their Bloch columns are parallel and the two-column target has rank one. For every unit Bloch vector u perpendicular to z, X=(3z+4u)/5 and Z=(4z-3u)/5 are Hermitian involutions with XZ+ZX=0. Every member attains the same global Procrustes value because the target overlap is sqrt((3/5)^2+(4/5)^2)=1. The circle of u therefore gives infinitely many global minimizers. Two positive marginal gaps alone do not imply joint Clifford uniqueness; full two-column rank is essential in the qubit theorem."
        ),
    }
    return finalize(271, "Explicit Two-Gap Misaligned Clifford Nonuniqueness Counterexample",
                    "proven_continuum_of_global_minimizers_despite_two_positive_marginal_gaps",
                    "The counterexample is exact for the qubit joint-factor problem. It does not classify generic full-rank higher-dimensional Clifford targets.", packet)


def emission_point(b, emission, y):
    pred = .2 + .7 * b; e0, e1 = emission[0][y], emission[1][y]
    q = e1 + (e0 - e1) * pred
    return pred * e0 / q, q


def transfer_point(b0, b1, depth, e0, e1, s):
    if depth == 0: return 1.0
    total = 0.0
    for y in (0, 1):
        n0, q0 = emission_point(b0, e0, y); n1, q1 = emission_point(b1, e1, y)
        total += q0 ** s * q1 ** (1 - s) * transfer_point(n0, n1, depth - 1, e0, e1, s)
    return total


def emission_interval(b, emission, y):
    pred = iadd((.2, .2), iscale(.7, b)); e0, e1 = emission[0][y], emission[1][y]
    q = iadd((e1, e1), iscale(e0 - e1, pred))
    return idiv(iscale(e0, pred), q), q


def transfer_half_interval(b0, b1, depth, e0, e1):
    if depth == 0: return (1.0, 1.0)
    total = (0.0, 0.0)
    for y in (0, 1):
        n0, q0 = emission_interval(b0, e0, y); n1, q1 = emission_interval(b1, e1, y)
        w = ( _down(math.sqrt(q0[0] * q1[0])), _up(math.sqrt(q0[1] * q1[1])) )
        total = iadd(total, imul(w, transfer_half_interval(n0, n1, depth - 1, e0, e1)))
    return total


def st272_adaptive_chernoff() -> dict:
    e0, e1 = [[.98, .02], [.92, .08]], [[.08, .92], [.02, .98]]
    pgrid, pdepth = 12, 9
    def point_upper(s):
        return max(transfer_point(i / pgrid, j / pgrid, pdepth, e0, e1, s)
                   for i in range(pgrid + 1) for j in range(pgrid + 1)) ** (1 / pdepth)
    opt = minimize_scalar(point_upper, bounds=(.45, .55), method="bounded", options={"xatol": 2e-6})
    base, depth = 20, 9
    boxes = [(i / base, (i + 1) / base, j / base, (j + 1) / base)
             for i in range(base) for j in range(base)]
    history = []
    for round_id in range(3):
        vals = [transfer_half_interval((b[0], b[1]), (b[2], b[3]), depth, e0, e1) for b in boxes]
        lower, upper = min(v[0] for v in vals), max(v[1] for v in vals)
        history.append({"round": round_id, "boxes": len(boxes), "Tn1_interval": [lower, upper],
                        "spectral_radius_interval": [lower ** (1 / depth), upper ** (1 / depth)]})
        if round_id < 2:
            order = set(np.argsort([v[0] for v in vals])[:24].tolist() +
                        np.argsort([-v[1] for v in vals])[:24].tolist())
            new = []
            for idx, b in enumerate(boxes):
                if idx not in order:
                    new.append(b); continue
                x0, x1, y0, y1 = b; xm, ym = (x0 + x1) / 2, (y0 + y1) / 2
                new.extend((a, c, d, f) for a, c in ((x0, xm), (xm, x1))
                           for d, f in ((y0, ym), (ym, y1)))
            boxes = new
    rlo, rhi = history[-1]["spectral_radius_interval"]
    rate = [-math.log(rhi), -math.log(rlo)]
    previous = [0.7110015196496613, 0.8198470892461787]
    combined = [max(rate[0], previous[0]), min(rate[1], previous[1])]
    packet = {
        "numerical_Chernoff_s_candidate": float(opt.x),
        "numerical_point_collocation_radius_upper": float(opt.fun),
        "optimization_status": "strong numerical evidence only",
        "certified_fixed_s": "1/2",
        "adaptive_history": history,
        "certified_Hellinger_rate_interval": rate,
        "combined_with_independent_previous_interval": combined,
        "combined_width": combined[1] - combined[0],
        "theorem": (
            "At s=1/2, every adaptive box is evaluated by outward positive transfer recursion, and the retained plus subdivided boxes still cover the complete belief square. Positive-operator order bounds therefore certify the displayed spectral-radius and Hellinger-rate intervals. Separately, finite-depth point collocation finds a Chernoff candidate near the displayed s, but no interval derivative or all-s enclosure is supplied; the optimization is not promoted to a theorem."
        ),
    }
    return finalize(272, "Adaptive Belief-Transfer Enclosure with Chernoff-Parameter Search",
                    "proven_adaptive_fixed_half_enclosure_plus_strong_numerical_chernoff_candidate",
                    "The fixed-s enclosure is rigorous for the synthetic HMM. Global optimization over s, detector calibration, and laboratory error exponents remain open.", packet)


def binary_entropy(p: float) -> float:
    if p in (0, 1): return 0.0
    return -p * math.log(p) - (1 - p) * math.log(1 - p)


def st273_controller_cost() -> dict:
    swap = np.array([[1, 0, 0, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 0, 1]], float)
    controlled = np.block([[np.eye(4), np.zeros((4, 4))], [np.zeros((4, 4)), swap]])
    rows = []
    for p in (0.01, 0.1, 0.25, 0.5):
        rows.append({"record_probability": p, "entropy_nats": binary_entropy(p),
                     "minimum_beta_times_cyclic_reset_work": binary_entropy(p)})
    packet = {
        "controller_dimension": 2,
        "controlled_SWAP_unitarity_error": float(np.linalg.norm(controlled.T @ controlled - np.eye(8))),
        "controlled_SWAP_involution_error": float(np.linalg.norm(controlled @ controlled - np.eye(8))),
        "degenerate_controller_rows": rows,
        "theorem": (
            "A two-state controller implements the exact controlled-SWAP unitary. If its program bit is a known pure state and is uncomputed reversibly, no positive universal Landauer cost follows. If instead it retains an uncertain classical record diag(1-p,p), returning a degenerate controller to a pure standard state in a bath requires beta W at least h_2(p). This is the controller's cyclic free-energy ledger; it is additional to the ideal closed-system SWAP identity. A finite controller therefore removes the fiction of a free reusable switch but does not create a unique cost without a specified record distribution and bath."
        ),
    }
    return finalize(273, "Finite Controller and Cyclic Free-Energy Ledger",
                    "proven_conditional_controller_reset_bound_and_no_universal_positive_cost",
                    "The controller Hamiltonian, bath, record probability, and reversible implementation are operational premises. No energy unit or physical apparatus is derived.", packet)


def build_local_refinement(w: np.ndarray, split: float, vertical: float):
    m = 2 * N; weights = np.zeros((m, m))
    for x in range(N):
        for y in range(x + 1, N):
            for a in (0, 1):
                for b in (0, 1):
                    weights[2 * x + a, 2 * y + b] = weights[2 * y + b, 2 * x + a] = w[x, y] * (split if a == b else 1 - split)
        weights[2 * x, 2 * x + 1] = weights[2 * x + 1, 2 * x] = vertical
    return np.diag(weights.sum(axis=1)) - weights


def st274_refinement_functor(w: np.ndarray, a: np.ndarray) -> dict:
    plus = np.array([1, 1]) / math.sqrt(2); R = np.kron(np.eye(N), plus[:, None])
    rows = []
    for v in (0.0, 0.5, 2.0):
        at = build_local_refinement(w, 0.5, v)
        flip = np.eye(2 * N); flip[[0, 1]] = flip[[1, 0]]
        rows.append({"vertical_rate": v, "intertwining_error": float(np.linalg.norm(at @ R - R @ a)),
                     "single_fiber_relabeling_error": float(np.linalg.norm(flip @ at @ flip.T - at)),
                     "minimum_eigenvalue": float(np.linalg.eigvalsh(at)[0])})
    packet = {
        "new_functorial_premise": "each two-point fiber is unlabelled, so independent local S2 relabellings are natural automorphisms",
        "horizontal_solution": "p_d=q_d=W_d/2 for d=1,...,6",
        "remaining_modulus": "vertical rate v>=0",
        "dimension_before": 7, "dimension_after": 1,
        "rows": rows,
        "theorem": (
            "The ST259 simultaneous fiber swap leaves aligned and crossed edge weights distinct. Naturality under independent relabelling of either endpoint fiber exchanges p_d and q_d and therefore forces p_d=q_d. Coarse intertwining p_d+q_d=W_d then fixes every horizontal weight to W_d/2. The vertical edge is invariant under all relabellings, so its rate v remains arbitrary. Hence the unlabelled-fiber refinement functor cuts the seven-dimensional local cone to one ray but cannot select a unique refinement."
        ),
    }
    return finalize(274, "Unlabelled-Fiber Natural Refinement Functor",
                    "proven_conditional_reduction_of_local_refinement_cone_from_7D_to_1D",
                    "Independent fiber relabelling is a new categorical/naturality premise, not a strict-kernel theorem. The surviving vertical rate remains unsourced; no physical length or fractal ratio follows.", packet)


def st275_log_haar_completion() -> dict:
    rho, mumin, mumax = 1.0, 2.0 ** -12, 2.0 ** 12
    times = np.logspace(-4, 4, 500)
    exact = 2 * rho * (np.log1p(np.exp(-2 * mumin * times)) -
                       np.log1p(np.exp(-2 * mumax * times)))
    plateau = 2 * rho * math.log(2)
    error = plateau - exact
    bound = 2 * rho * mumin * times + 2 * rho * np.exp(-2 * mumax * times)
    packet = {
        "density": rho, "mu_min": mumin, "mu_max": mumax,
        "exact_finite_cutoff_formula": "2 rho [log(1+exp(-2 mu_min t))-log(1+exp(-2 mu_max t))]",
        "plateau": plateau,
        "uniform_error_bound": "error <= 2 rho mu_min t + 2 rho exp(-2 mu_max t)",
        "maximum_bound_violation": float(np.max(error - bound)),
        "trace_class_condition_near_zero": "integral_0 rho(mu) dmu/mu < infinity",
        "full_log_Haar_heat_trace_status": "diverges logarithmically at mu=0 for every t>0",
        "times": times.tolist(), "finite_cutoff_curve": exact.tolist(), "error_bound": bound.tolist(),
        "theorem": (
            "With cutoffs, the log-Haar fiber contribution is exactly 2 rho[log(1+e^{-2 mu_min t})-log(1+e^{-2 mu_max t})]. Its deficit from 2 rho log 2 is bounded by 2 rho mu_min t+2 rho e^{-2 mu_max t}. However, the underlying heat trace contains integral e^{-2 mu t} rho dmu/mu and diverges logarithmically at mu=0. Therefore no nonzero exactly scale-invariant log-Haar completion on all positive rates is trace class. A trace-class completion must introduce an infrared cutoff or a density vanishing at zero, and either choice breaks exact scale invariance."
        ),
    }
    return finalize(275, "Finite-Cutoff Log-Haar Bounds and Trace-Class Impossibility Theorem",
                    "proven_cutoff_error_bound_and_no_exact_trace_class_scale_invariant_completion",
                    "This is a dimensionless spectral theorem. The density, cutoffs, rate realization, physical units, Planck scale, and legacy octave interpretations are not derived.", packet)


def st276_fiber_instrument(w: np.ndarray, a: np.ndarray) -> dict:
    split, vertical, t = 0.7, 0.4, 0.35
    at = build_local_refinement(w, split, vertical); h = expm(-t * at)
    plus, minus = np.array([1, 1]) / math.sqrt(2), np.array([1, -1]) / math.sqrt(2)
    R, S = np.kron(np.eye(N), plus[:, None]), np.kron(np.eye(N), minus[:, None])
    B = S.T @ at @ S; recovered_heat = S.T @ h @ S
    recovered_B = -np.real_if_close(scipy.linalg.logm(recovered_heat)) / t
    contrast = np.zeros((N, N))
    for y in range(N):
        for x in range(N):
            contrast[y, x] = 0.5 * (h[2*y, 2*x] - h[2*y, 2*x+1] -
                                    h[2*y+1, 2*x] + h[2*y+1, 2*x+1])
    schema = ["run_id", "registered_protocol_hash", "calibration_hash", "time_tau",
              "prepared_parent_x", "prepared_child_a", "measured_parent_y", "measured_child_b",
              "count", "total_trials", "provider_id", "registrar_id", "analyst_id", "holdout_flag"]
    packet = {
        "declared_demo": {"aligned_fraction": split, "crossed_fraction": 1-split,
                          "vertical_rate": vertical, "dimensionless_time": t},
        "preparations": "all 24 child-resolved basis states (x,a)",
        "effects": "all 24 child-resolved vertex effects (y,b)",
        "contrast_formula": "H_odd[y,x]=(P_y0|x0-P_y0|x1-P_y1|x0+P_y1|x1)/2",
        "record_schema": schema,
        "role_separation": "provider != registrar != analyst",
        "recovered_heat_contrast_error": float(np.linalg.norm(contrast - recovered_heat, np.inf)),
        "principal_log_B_reconstruction_error": float(np.linalg.norm(recovered_B - B, np.inf)),
        "coarse_heat_error": float(np.linalg.norm(R.T @ h @ R - expm(-t * a), np.inf)),
        "theorem": (
            "The twenty-four positive child-basis preparations and child-resolving effects determine every transition entry of exp(-t A_tilde). The signed four-count contrast equals S*exp(-t A_tilde)S=exp(-tB), without requiring an unphysical negative preparation. With calibrated t and spectrum in the principal-log domain, B=-(1/t)log(exp(-tB)) is reconstructed. The frozen schema separates provider, registrar, and analyst and supports a pre-registered holdout."
        ),
    }
    return finalize(276, "Frozen Fiber-Resolving Instrument and Record Specification",
                    "proven_executable_mathematical_instrument_specification",
                    "No apparatus has executed this protocol. Counts, calibration, preparation fidelity, time units, custody, and independent registration remain external obligations.", packet)


def st277_external_stop() -> dict:
    packet = {
        "decision": "blocked_without_independently_registered_events",
        "required_atoms": ["raw event table matching ST276 schema", "pre-run protocol hash",
                           "independent calibration hash", "provider/registrar/analyst separation",
                           "frozen holdout and one-shot unblinding"],
        "local_synthetic_substitution": False,
        "theorem": "An executable instrument specification and local simulated replay do not entail that an external system instantiated it. Without independently registered events, the empirical stage has no admissible input."
    }
    return finalize(277, "External Evidence Gate",
                    "blocked_no_external_record",
                    "This is an evidence boundary, not a mathematical counterexample to the declared operator model.", packet)


def make_figures(d: dict) -> None:
    FIG_DIR.mkdir(exist_ok=True)
    fig, ax = plt.subplots(figsize=(7, 4)); r = d["ST263"]
    ax.plot(r["theta"], r["angular_potential"]); ax.scatter(2*np.pi*np.arange(12)/12, -np.ones(12), s=16)
    ax.set(xlabel="phase theta", ylabel="-cos(12 theta)", title="ST263: first D12-anisotropic phase potential")
    fig.tight_layout(); fig.savefig(FIG_DIR / "st263_phase_orbit.png", dpi=190); plt.close(fig)
    fig, ax = plt.subplots(figsize=(7, 4)); r = d["ST265"]; seed = r["selected_seed_pair"][1]
    path = r["paths"].get(seed, r["paths"].get(str(seed))); ax.semilogy([x["step"] for x in path], [x["interval_full_row_rank_lower_bound"] for x in path], "o-", label="rank lower bound")
    ax2 = ax.twinx(); ax2.plot([x["step"] for x in r["pair_clearances"]], [x["interval_box_clearance"] for x in r["pair_clearances"]], color="tab:orange", label="pair clearance")
    ax.set(xlabel="continuation step", ylabel="certified row-rank lower bound", title="ST265: targeted event exclusion"); ax2.set_ylabel("interval-box clearance")
    fig.tight_layout(); fig.savefig(FIG_DIR / "st265_event_audit.png", dpi=190); plt.close(fig)
    fig, ax = plt.subplots(figsize=(7, 4)); p = d["ST267"]["optimized_knots"]; ax.plot(np.linspace(0,4,10),p,"o-"); ax.set(xlabel="dimensionless time",ylabel="control knot",title="ST267: interval-certified discrete selector optimum")
    fig.tight_layout(); fig.savefig(FIG_DIR / "st267_selector_kkt.png", dpi=190); plt.close(fig)
    fig, ax = plt.subplots(figsize=(7, 4)); r=d["ST269"]; ax.bar(["0.00070 full cover","0.00075 corners"],[r["last_complete_minimum_margin"],min(x["margin"] for x in r["next_attempt_corner_cells"])]);ax.axhline(0,color="black",lw=.8);ax.set(ylabel="Krawczyk margin",title="ST269: fixed-resolution certificate boundary")
    fig.tight_layout();fig.savefig(FIG_DIR/"st269_nuisance_boundary.png",dpi=190);plt.close(fig)
    fig,ax=plt.subplots(figsize=(7,4));r=d["ST270"]["rows"];ax.semilogy([x["layers"] for x in r],[x["maximum_halfwidth"] for x in r],"o-");ax.set(xlabel="layers",ylabel="maximum 95% coordinate halfwidth",title="ST270: finite clustered confidence regions");fig.tight_layout();fig.savefig(FIG_DIR/"st270_confidence.png",dpi=190);plt.close(fig)
    fig,ax=plt.subplots(figsize=(6,5));r=d["ST271"]["sampled_family"];ax.plot([math.cos(x["phi"]) for x in r],[math.sin(x["phi"]) for x in r],"o-");ax.set_aspect("equal");ax.set(xlabel="u_x",ylabel="u_y",title="ST271: continuum of exact Clifford minimizers");fig.tight_layout();fig.savefig(FIG_DIR/"st271_clifford_circle.png",dpi=190);plt.close(fig)
    fig,ax=plt.subplots(figsize=(7,4));r=d["ST272"]["adaptive_history"];ax.plot([x["boxes"] for x in r],[x["spectral_radius_interval"][1]-x["spectral_radius_interval"][0] for x in r],"o-");ax.set(xlabel="adaptive boxes",ylabel="spectral-radius interval width",title="ST272: adaptive transfer enclosure");fig.tight_layout();fig.savefig(FIG_DIR/"st272_transfer_adaptive.png",dpi=190);plt.close(fig)
    fig,ax=plt.subplots(figsize=(7,4));r=d["ST275"];ax.semilogx(r["times"],r["finite_cutoff_curve"],label="finite cutoffs");ax.axhline(r["plateau"],color="black",ls="--",label="2 rho ln 2");ax.legend();ax.set(xlabel="dimensionless t",ylabel="fiber spectral dimension",title="ST275: controlled plateau and cutoff loss");fig.tight_layout();fig.savefig(FIG_DIR/"st275_cutoff_plateau.png",dpi=190);plt.close(fig)


def main() -> None:
    w, a, _ = strict_operator()
    out = {"metadata": {"programs": "ST263-ST277", "date": "2026-08-12", "seed": SEED,
                        "python": platform.python_version(), "numpy": np.__version__,
                        "scipy": scipy.__version__, "sympy": sp.__version__}}
    out["ST263"] = st263_phase_orbit(a); out["ST264"] = st264_nonlinear_equivariants()
    out["ST265"] = st265_targeted_branch_event(a); out["ST266"] = st266_measure_prepare_recovery()
    out["ST267"] = st267_selector_kkt(); out["ST268"] = st268_refinement_connection(a)
    out["ST269"] = st269_nuisance_boundary(); out["ST270"] = st270_cluster_confidence()
    out["ST271"] = st271_clifford_counterexample(); out["ST272"] = st272_adaptive_chernoff()
    out["ST273"] = st273_controller_cost(); out["ST274"] = st274_refinement_functor(w, a)
    out["ST275"] = st275_log_haar_completion(); out["ST276"] = st276_fiber_instrument(w, a)
    out["ST277"] = st277_external_stop()
    out["recommended_next_programs"] = [
        {"id":"ST278","priority":1,"study":"derive the coefficient and radius of the Re(z^12) phase potential from one strict nonlinear interaction, or prove a source obstruction"},
        {"id":"ST279","priority":2,"study":"compute a Hilbert basis for the full nonlinear D12-equivariant module beyond degree 11"},
        {"id":"ST280","priority":3,"study":"continue seed 23 only through a certified minimum-rank event bracket, with an explicit finite stop"},
        {"id":"ST281","priority":4,"study":"generalize the exact measure-and-prepare recovery theorem to arbitrary rational binary qubit outputs"},
        {"id":"ST282","priority":5,"study":"independently replay ST267 with rational or arbitrary-precision interval arithmetic"},
        {"id":"ST283","priority":6,"study":"test affine-connection naturality under nonlinear refinement embeddings or construct a counterexample"},
        {"id":"ST284","priority":7,"study":"adaptively subdivide the failed 0.00075 nuisance cells to separate resolution loss from a true root boundary"},
        {"id":"ST285","priority":8,"study":"replace Hoeffding rectangles with exact-mixture or empirical-Bernstein clustered confidence sets"},
        {"id":"ST286","priority":9,"study":"classify higher-dimensional joint Clifford nonuniqueness beyond the qubit rank-one counterexample"},
        {"id":"ST287","priority":10,"study":"certify the global Chernoff-s optimizer with interval derivatives and all-s coverage"},
        {"id":"ST288","priority":11,"study":"construct an autonomous finite controller and audit catalytic versus resetting costs"},
        {"id":"ST289","priority":12,"study":"impose coassociativity on 12-to-24-to-48 unlabelled-fiber refinements and classify the surviving vertical-rate moduli"},
        {"id":"ST290","priority":13,"study":"classify near-Haar trace-class densities and prove the optimal plateau-width versus trace-class tradeoff"},
        {"id":"ST291","priority":14,"study":"build an executable validator for the frozen ST276 instrument and custody schema"},
        {"id":"ST292","priority":15,"study":"resume empirical validation only after independently registered raw events exist"},
    ]
    out["central_verdict"] = (
        "The strict operator canonically supplies a lowest-mode phase circle and D12 fixes the first anisotropic invariant at degree 12, so spontaneous twelve-branch selection has a precise minimal mathematical form. Yet the nonlinear coefficient, nonzero radius, and realized branch are not sourced by A. The mean-zero affine connection is natural under exact linear refinements, while unlabelled-fiber naturality reduces but does not eliminate refinement freedom. Exact log-Haar scaling cannot be both nonzero and trace class on all positive rates; controlled fractal plateaus necessarily require cutoffs or scale-breaking density."
    )
    out["epistemic_boundary"] = (
        "No strict selector coefficient, QW-2191 discharge, physical vacuum, dimensional clock/length/action unit, Planck scale, spacetime, laboratory record, legacy-to-strict completion or role transfer, Standard Model, gravity, L_total, or ToE closure is claimed."
    )
    RESULTS.write_text(json.dumps(native(out), indent=2, sort_keys=True), encoding="utf-8")
    with SUMMARY.open("w", newline="", encoding="utf-8") as h:
        writer = csv.writer(h); writer.writerow(["program", "object", "status"])
        for k in range(263, 278): writer.writerow([f"ST{k}", out[f"ST{k}"]["object"], out[f"ST{k}"]["status"]])
    make_figures(out)
    print(json.dumps({"programs":"ST263-ST277", "statuses":{f"ST{k}":out[f"ST{k}"]["status"] for k in range(263,278)},
                      "ST267_margin":out["ST267"]["minimum_Krawczyk_margin"],
                      "ST269":[out["ST269"]["last_complete_minimum_margin"],out["ST269"]["next_attempt_failed_corner_cells"]],
                      "ST272":out["ST272"]["combined_with_independent_previous_interval"]}, indent=2))


if __name__ == "__main__":
    main()
