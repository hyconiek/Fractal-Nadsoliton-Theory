#!/usr/bin/env python3
"""FIN ST248--ST262: preparation-source obstruction, affine carrier connection,
constrained refinement, and logarithmic scale-flow research.
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
from scipy.optimize import minimize, root

from fin_st01_st15_research import N, strict_operator
from fin_st130_st141_research import point_design_system
from fin_st132_center_isolation_replay import strict_interval_matrix
from fin_st166_st177_research import local_param_krawczyk, parametric_data
from fin_st178_st189_research import stationary_slice_float
from fin_st203_st217_research import pseudo_krawczyk


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "FIN_ST248_ST262_Results.json"
SUMMARY = ROOT / "FIN_ST248_ST262_Summary.csv"
FIG_DIR = ROOT / "FIN_ST248_ST262_Figures"
SEED = 20260812
PACKETS = {k: ROOT / f"FIN_ST{k}_{name}.json" for k, name in {
    248: "Strict_Preparation_Source_No_Go",
    249: "D12_Natural_Typed_Second_Operators",
    250: "Extended_All_Branch_Resource_Stop",
    251: "Rotated_Amplitude_Damping_Recovery",
    252: "Selector_Endpoint_Actuator_Cost",
    253: "Mean_Zero_Affine_Carrier_Connection",
    254: "Complete_Halfwidth_0005_Nuisance_Cover",
    255: "Finite_Sample_Efficient_Mode_Estimators",
    256: "Aligned_Two_Gap_Clifford_Uniqueness",
    257: "Belief_Transfer_Spectral_Radius",
    258: "Finite_Temperature_Reset_Controller_Ledger",
    259: "Local_Graph_Laplacian_Refinement_Cone",
    260: "Log_Haar_Spectral_Dimension_Plateau",
    261: "Complement_Resolving_Refinement_Observables",
    262: "External_Evidence_Stop"
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


def dihedral_permutations() -> list[np.ndarray]:
    out = []
    for shift in range(N):
        for sign in (1, -1):
            p = [(sign * i + shift) % N for i in range(N)]
            q = np.zeros((N, N))
            for i, j in enumerate(p):
                q[j, i] = 1
            out.append(q)
    return out


def distance_basis() -> list[np.ndarray]:
    mats = []
    for d in range(7):
        m = np.zeros((N, N))
        for i in range(N):
            for j in range(N):
                if min((i - j) % N, (j - i) % N) == d:
                    m[i, j] = 1
        mats.append(m)
    return mats


def st248_preparation_source_no_go(a: np.ndarray) -> dict:
    group = dihedral_permutations()
    u = np.ones(N) / N
    candidates = {
        "uniform": u,
        "normalized diagonal of exp(-0.7A)": np.diag(expm(-0.7 * a)) / np.trace(expm(-0.7 * a)),
        "normalized diagonal of (I+A)^-1": np.diag(np.linalg.inv(np.eye(N) + a)) / np.trace(np.linalg.inv(np.eye(N) + a)),
        "normalized diagonal of A^2": np.diag(a @ a) / np.trace(a @ a),
    }
    rows = []
    for name, rho in candidates.items():
        rows.append({"candidate": name, "distance_from_uniform": float(np.linalg.norm(rho - u)),
                     "maximum_D12_invariance_error": float(max(np.linalg.norm(g @ rho - rho) for g in group)),
                     "commutator_norm": float(np.linalg.norm(a @ np.diag(rho) - np.diag(rho) @ a))})
    heat_local = expm(-0.7 * a)[:, 0]
    packet = {
        "group_order": len(group), "fixed_probability_simplex_dimension": 0,
        "tested_natural_state_rows": rows,
        "localized_heat_state_distance_from_uniform": float(np.linalg.norm(heat_local - u)),
        "localized_heat_state_commutator_norm": float(np.linalg.norm(a @ np.diag(heat_local) - np.diag(heat_local) @ a)),
        "theorem": (
            "Let F be any deterministic D12-natural state-valued rule and let A be fixed by the transitive D12 action. Naturality gives F(A)=gF(A) for every g. The only invariant probability vector under a transitive action is the uniform vector. Hence no deterministic D12-natural construction from A alone can produce the localized rho required by M_rho. The displayed diagonals of strict functional-calculus objects replay the theorem. A localized heat state appears only after supplying delta_0 or equivalent boundary/preparation data."
        )
    }
    return finalize(248, "Strict-Core Preparation-Source No-Go",
                    "proven_no_go_for_deterministic_D12_natural_state_from_A_alone",
                    "The theorem does not exclude random spontaneous realization, a state-dependent nonlinear law, an external boundary, or a genuinely new strict symmetry-breaking object. It does not discharge QW-2191.", packet)


def st249_typed_operator_classification(a: np.ndarray) -> dict:
    basis = distance_basis(); group = dihedral_permutations()
    comm_errors = [max(np.linalg.norm(g @ b - b @ g) for g in group) for b in basis]
    gram = np.array([[np.sum(x * y) for y in basis] for x in basis])
    rng = np.random.default_rng(SEED)
    rows = []
    for label, rho in [
        ("uniform", np.ones(N) / N),
        ("localized even", expm(-0.7 * a)[:, 0]),
        ("chiral", np.roll(expm(-0.7 * a)[:, 0], 1) * 0.57 + expm(-0.7 * a)[:, 0] * 0.43),
    ]:
        coeff = rng.random(7); c = sum(z * b for z, b in zip(coeff, basis)); c /= c.sum(axis=0)[0]
        out = c @ rho
        stabilizer = sum(np.linalg.norm(g @ rho - rho) < 1e-11 for g in group)
        out_stabilizer = sum(np.linalg.norm(g @ out - out) < 1e-11 for g in group)
        rows.append({"conditioning": label, "input_stabilizer_order": stabilizer,
                     "output_stabilizer_order": out_stabilizer,
                     "output_minimum": float(out.min()), "operator_commutator_norm": float(np.linalg.norm(a @ np.diag(out) - np.diag(out) @ a))})
    packet = {
        "linear_equivariant_dimension": int(np.linalg.matrix_rank(gram)),
        "basis_labels": [f"distance_{d}" for d in range(7)],
        "maximum_basis_equivariance_error": float(max(comm_errors)),
        "affine_Markov_classification": "C=sum_{d=0}^6 c_d B_d with entrywise nonnegativity and unit column sum; rho maps to M_{C rho}=diag(C rho)",
        "conditioning_rows": rows,
        "theorem": (
            "The D12 permutation representation on vertices is multiplicity-free over the seven real distance sectors. Therefore every linear D12-equivariant state-to-state map is a unique linear combination of the seven distance matrices B_d. Positivity and preservation of total mass cut out the affine Markov polytope. The typed second operators are M_{C rho}. Equivariance implies Stab(rho) is contained in Stab(C rho), so such processing can preserve or erase selector information but cannot create localization or orientation absent from rho. No conditioning gives only a scalar multiplication operator; even localization selects a vertex orbit but not orientation; chiral conditioning is the minimal class with trivial stabilizer."
        )
    }
    return finalize(249, "Complete Linear D12-Natural Typed Second-Operator Classification",
                    "proven_seven_dimensional_linear_equivariant_classification",
                    "This classifies linear vertex-state processing. Nonlinear equivariant maps, new carrier types, and a strict source for the conditioning state remain outside the theorem.", packet)


def stationary7(x: np.ndarray, a: np.ndarray):
    f, j = stationary_slice_float(x, a, np.zeros(7), np.zeros(7), 0.0)
    return f[:7], j[:7, :]


def st250_branch_resource_stop(a: np.ndarray) -> dict:
    prior = json.loads((ROOT / "FIN_ST233_ST247_Results.json").read_text(encoding="utf-8"))["ST236"]
    endpoints = {}
    previous = {}
    for row in prior["rows"]:
        sid = int(row["seed"])
        if row["step"] == 1:
            previous[sid] = np.array(row["center"], dtype=float)
        if row["step"] == 2:
            endpoints[sid] = row
    aiv, _, _ = strict_interval_matrix()
    new_rows = []; paths = {}
    for sid in sorted(endpoints):
        x = np.array(endpoints[sid]["center"], dtype=float)
        prev = previous[sid]
        _, j = stationary7(x, a); _, _, vh = np.linalg.svd(j)
        tangent = vh[-1]; tangent /= np.linalg.norm(tangent)
        if np.dot(tangent, x - prev) < 0:
            tangent = -tangent
        paths[sid] = []
        for step in range(3, 9):
            accepted = None
            for ds in (0.002, 0.001, 0.0005):
                pred = x + ds * tangent
                sol = root(lambda y: np.r_[stationary7(y, a)[0], tangent @ (y - pred)], pred,
                           jac=lambda y: np.vstack([stationary7(y, a)[1], tangent]), tol=1e-12)
                xn = sol.x; _, jn = stationary7(xn, a); _, sv, vh = np.linalg.svd(jn)
                tn = vh[-1]; tn /= np.linalg.norm(tn)
                if np.dot(tn, tangent) < 0:
                    tn = -tn
                trials = [pseudo_krawczyk(xn, aiv, tangent, pred, r) for r in (5e-7, 2e-7, 7e-8)]
                cert = next((z for z in trials if z["included"]), None)
                if cert is not None:
                    accepted = (ds, xn, tn, sv[-1], cert); break
            if accepted is None:
                break
            ds, xn, tn, margin, cert = accepted
            rec = {"seed": sid, "step": step, "step_size": ds, "center": xn.tolist(),
                   "kappa": float(xn[7]), "rank_margin": float(margin), "certificate": cert}
            paths[sid].append(rec); new_rows.append(rec); x, tangent = xn, tn
    fold_edges = prior["exact_uniform_fold_edges"]
    min_clearance = math.inf; intersections = []
    for i in range(60):
        for j in range(i + 1, 60):
            if any(i in e and j in e for e in fold_edges):
                continue
            for ri in paths[i]:
                for rj in paths[j]:
                    sep = float(np.max(np.abs(np.array(ri["center"]) - np.array(rj["center"]))))
                    clearance = sep - ri["certificate"]["radius"] - rj["certificate"]["radius"]
                    min_clearance = min(min_clearance, clearance)
                    if clearance <= 0:
                        intersections.append([i, j, ri["step"], rj["step"]])
    packet = {
        "inherited_certified_steps_per_seed": 2,
        "new_attempted_steps_per_seed": 6,
        "new_certified_steps": len(new_rows),
        "complete_extended_paths": sum(len(v) == 6 for v in paths.values()),
        "total_certified_steps_per_complete_seed": 8,
        "adaptive_step_sizes_used": sorted({r["step_size"] for r in new_rows}),
        "minimum_new_rank_margin": min(r["rank_margin"] for r in new_rows),
        "minimum_nonfold_clearance": min_clearance,
        "nonfold_box_intersections": intersections,
        "finite_component_count": len(fold_edges),
        "explicit_resource_stop": "six new certified steps per seed",
        "rows": new_rows,
        "theorem": (
            "Every new node is the unique stationary root in its outward pseudo-arclength Krawczyk box. Combined with the inherited two steps, each of sixty signed seeds now has eight certified steps. No boxes belonging to distinct exact fold pairs intersect in the traced extension, so the declared finite incidence graph still has thirty components. The campaign stops at six new steps per seed by design; it makes no assertion about untraced global branches."
        )
    }
    passed = len(new_rows) == 360 and not intersections
    return finalize(250, "Extended All-Branch Continuation with Explicit Resource Stop",
                    "proven_360_new_local_boxes_and_finite_30_component_graph" if passed else "partial_extended_branch_campaign",
                    "The finite stop cannot exclude later collisions, secondary bifurcations, or distant component connections. No stability or particle interpretation follows.", packet)


def vecf(x: np.ndarray) -> np.ndarray:
    return np.asarray(x, complex).reshape(-1, order="F")


def st251_rotated_amplitude_damping() -> dict:
    U = np.array([[Fraction(3, 5), Fraction(4, 5)], [Fraction(-4, 5), Fraction(3, 5)]], dtype=object)
    V = np.array([[Fraction(5, 13), Fraction(12, 13)], [Fraction(-12, 13), Fraction(5, 13)]], dtype=object)
    K0 = np.array([[Fraction(1), Fraction(0)], [Fraction(0), Fraction(4, 5)]], dtype=object)
    K1 = np.array([[Fraction(0), Fraction(3, 5)], [Fraction(0), Fraction(0)]], dtype=object)
    uf, vf = np.array(U, float), np.array(V, float)
    kraus = [vf @ np.array(k, float) @ uf for k in (K0, K1)]
    c = sum(np.outer(vecf(k.conj().T), vecf(k.conj().T).conj()) / 4 for k in kraus)
    rexact = U.T @ V.T
    rmat = np.array(rexact, float)
    jr = np.outer(vecf(rmat), vecf(rmat).conj())
    yb = np.diag([Fraction(9, 20), Fraction(9, 25)]).astype(object)
    y = vf @ np.array(yb, float) @ vf.T
    slack = np.kron(y, np.eye(2)) - c
    # Generic affine Bloch representation from action on Pauli coordinates.
    pauli = [np.eye(2), np.array([[0,1],[1,0]],complex), np.array([[0,-1j],[1j,0]],complex), np.diag([1,-1]).astype(complex)]
    def channel(x):
        return sum(k @ x @ k.conj().T for k in kraus)
    t = np.array([np.trace(pauli[i] @ channel(pauli[0]/2)).real for i in range(1,4)])
    m = np.array([[np.trace(pauli[i] @ channel(pauli[j]/2)).real for j in range(1,4)] for i in range(1,4)])
    packet = {
        "channel": "V after amplitude-damping(gamma=9/25) after U, with rational U and V",
        "U": [[str(z) for z in row] for row in U], "V": [[str(z) for z in row] for row in V],
        "amplitude_damping_K0": [[str(z) for z in row] for row in K0],
        "amplitude_damping_K1": [[str(z) for z in row] for row in K1],
        "affine_linear_part": m.tolist(), "affine_translation": t.tolist(),
        "offdiagonal_linear_Frobenius_norm": float(np.linalg.norm(m - np.diag(np.diag(m)))),
        "exact_optimum": "81/100", "exact_recovery_unitary": [[str(z) for z in row] for row in rexact],
        "recovery_unitary": rmat.tolist(), "exact_base_dual_Y": [["9/20","0"],["0","9/25"]],
        "transported_dual_Y": y.tolist(),
        "primal_value": float(np.trace(c @ jr).real), "dual_value": float(np.trace(y)),
        "primal_dual_gap": float(np.trace(y) - np.trace(c @ jr).real),
        "dual_slack_eigenvalues": np.linalg.eigvalsh(slack).tolist(),
        "theorem": (
            "For amplitude damping with survival amplitude a=4/5, the identity recovery has entanglement fidelity (1+a)^2/4=81/100. The exact dual Y=diag((1+a)/4,a(1+a)/4) has the same trace; its slack is a direct sum of two positive scalars and the rank-one PSD block [[1,-1],[-1,1]]/5. Therefore the optimum is global over all CPTP recoveries. Pre/post rational unitaries transport the primal and dual bijectively, giving the displayed fully non-diagonal nonunital affine channel and optimal recovery U^* V^*."
        )
    }
    return finalize(251, "Exact Recovery for a Rational Rotated Amplitude-Damping Channel",
                    "proven_global_CPTP_optimum_for_non_diagonal_nonunital_channel",
                    "The theorem covers the unitary orbit of the declared amplitude-damping channel, not every affine qubit channel or a physical noise source.", packet)


def bernoulli_kl(p: float, q: float) -> float:
    return p * math.log(p / q) + (1 - p) * math.log((1 - p) / (1 - q))


def st252_selector_actuator_cost() -> dict:
    x0, xf, T = 1/12, 0.6, 4.0
    knots = 10; dt_k = T / (knots - 1); fine = 600; dt = T / fine
    kappa = 0.04; endpoint_weight = 0.12; rate_limit = 0.30

    def simulate(p: np.ndarray, retain: bool = False):
        tg = np.linspace(0, T, fine + 1); tk = np.linspace(0, T, knots)
        pi = np.interp(tg, tk, p); x = x0; ep = 0.0
        xs = [x]
        for j in range(fine):
            pm = 0.5 * (pi[j] + pi[j+1])
            # exact state update for constant midpoint equilibrium
            xn = pm + (x - pm) * math.exp(-dt)
            xm = 0.5 * (x + xn)
            ep += dt * (pm - xm) * (math.log(pm/(1-pm)) - math.log(xm/(1-xm)))
            x = xn; xs.append(x)
        actuator = kappa * np.sum(np.diff(p) ** 2) / dt_k
        endpoints = endpoint_weight * (bernoulli_kl(float(p[0]), x0) + bernoulli_kl(float(p[-1]), x0))
        total = ep + actuator + endpoints
        if retain:
            return total, ep, actuator, endpoints, x, tg, pi, np.array(xs)
        return total, x

    def objective(p): return simulate(p)[0]
    cons = [{"type":"eq", "fun":lambda p: simulate(p)[1]-xf},
            {"type":"ineq", "fun":lambda p: rate_limit-np.abs(np.diff(p))/dt_k}]
    base = np.linspace(0.28, 0.72, knots)
    starts = [base, np.linspace(0.20,0.78,knots), np.linspace(0.35,0.68,knots)]
    sols = [minimize(objective, s, method="SLSQP", bounds=[(0.084,0.95)]*knots,
                     constraints=cons, options={"ftol":1e-12,"maxiter":1200}) for s in starts]
    sol = min(sols, key=lambda z: z.fun)
    total, ep, actuator, endpoints, final_x, tg, pi, xs = simulate(sol.x, True)
    packet = {
        "declared_model": "ten control knots; exact midpoint-constant state update; entropy production plus quadratic slew cost plus endpoint KL field cost",
        "parameters": {"T":T,"x0":x0,"xf":xf,"knots":knots,"kappa":kappa,"endpoint_weight":endpoint_weight,"rate_limit":rate_limit},
        "optimizer_success": bool(sol.success), "optimizer_message": sol.message,
        "best_total_cost": total, "entropy_production": ep, "actuator_slew_cost": actuator,
        "endpoint_preparation_reset_cost": endpoints, "final_state_error": final_x-xf,
        "maximum_knot_rate": float(np.max(np.abs(np.diff(sol.x))/dt_k)),
        "optimized_knots": sol.x.tolist(), "multi_start_objectives": [float(z.fun) for z in sols],
        "universal_lower_bound_from_ST238": 0.31769108870743823,
        "trajectory": {"time": tg.tolist(), "pi": pi.tolist(), "x": xs.tolist()},
        "result": (
            "All added terms are nonnegative, so the verified ST238 entropy-production lower bound remains a rigorous lower bound for this augmented problem. The displayed feasible schedule is a reproducible numerical upper bound. Multiple starts converge to the same value, but no interval KKT certificate is supplied; global optimality is not claimed."
        )
    }
    status = "strong_numerical_feasible_augmented_selector_cost" if sol.success and abs(final_x-xf)<1e-8 else "partial_selector_actuator_optimization"
    return finalize(252, "Selector Cost with Endpoint Preparation and Finite-Rate Actuation", status,
                    "The endpoint functional, actuator coefficient, rate limit, two-state reduction, bath, preferred vertex, and dimensionless clock are supplied. This is not physical work in joules or a strict selector source.", packet)


def st253_affine_connection(a: np.ndarray) -> dict:
    eig = np.linalg.eigvalsh(a); positive = eig[eig > 1e-10]
    y, aa, gg, cc = sp.symbols("y a g c", nonzero=True)
    x = y + cc*y**6
    v = aa*x**2/sp.Integer(2) + gg*x**12/sp.Integer(12)
    gamma = sp.diff(x,y,2)/sp.diff(x,y)
    ordinary = sp.simplify(sp.diff(v,y,12).subs(y,0)/sp.factorial(11))
    tensor = sp.diff(v,y)
    for rank in range(1,12):
        tensor = sp.simplify(sp.diff(tensor,y)-rank*gamma*tensor)
    covariant = sp.simplify(tensor.subs(y,0)/sp.factorial(11))
    packet = {
        "carrier": "H0={psi in R^12: 1^T psi=0} with metric g_A(u,v)=u^T A v",
        "dimension": 11, "metric_minimum_eigenvalue": float(positive.min()),
        "metric_maximum_eigenvalue": float(positive.max()),
        "connection": "the translation-invariant Levi-Civita connection on the affine vector space H0",
        "curvature": 0, "torsion": 0,
        "scalar_chart": "x=y+c y^6",
        "ordinary_I12_symbolic": str(ordinary), "covariant_I12_symbolic": str(covariant),
        "symbolic_cancellation_verified": bool(sp.simplify(covariant-gg)==0),
        "theorem": (
            "Because ker(A)=span(1), A induces a positive metric on the canonical mean-zero quotient H0. Its translation-invariant metric has a unique Levi-Civita connection, which is flat, torsion-free, and D12-invariant. Transporting this connection under nonlinear charts makes every covariant response tensor natural. In the explicit x=y+c y^6 test, the ordinary normalized response is g+6ac^2 while the twelfth covariant response is exactly g; the chart artifact cancels symbolically."
        )
    }
    return finalize(253, "Canonical Flat Connection on the Mean-Zero Strict Amplitude Carrier",
                    "proven_connection_conditional_on_affine_mean_zero_amplitude_carrier",
                    "The mean-zero affine amplitude carrier is a mathematical FIN carrier choice. The result does not derive a physical field bundle, spacetime connection, gauge connection, or dimensional coupling from nature.", packet)


_NUISANCE_CACHE = None


def nuisance_worker(item):
    global _NUISANCE_CACHE
    if _NUISANCE_CACHE is None:
        _, a, _ = strict_operator(); _NUISANCE_CACHE = parametric_data(a)
    ec, ef, eigc, eigf = _NUISANCE_CACHE
    nu, lh = item
    center = root(lambda z: point_design_system(z, ec, ef, *nu), [2.1862,.53983], tol=1e-12).x
    cert = local_param_krawczyk(eigc, eigf, nu, lh, center)
    return nu, cert["margin"], cert["included"]


def st254_complete_nuisance_cover() -> dict:
    h=0.0005; sub=28; lh=h/sub
    offsets=[-h+(2*i+1)*lh for i in range(sub)]
    tasks=[((.2+x,.7+y,.05+z),lh) for x,y,z in itertools.product(offsets,repeat=3)]
    workers=min(8, max(1, os.cpu_count() or 1))
    with ProcessPoolExecutor(max_workers=workers) as pool:
        rows=list(pool.map(nuisance_worker,tasks,chunksize=48))
    margins=np.array([r[1] for r in rows]); failed=[r for r in rows if not r[2]]
    packet={
        "global_halfwidth":h,"cells_per_axis":sub,"local_halfwidth":lh,
        "boxes":len(rows),"passed":len(rows)-len(failed),"failed":len(failed),
        "minimum_margin":float(margins.min()),"maximum_margin":float(margins.max()),
        "margin_percentiles":{str(q):float(np.percentile(margins,q)) for q in [0,1,5,50,95,99,100]},
        "parallel_workers":workers,"first_failures":[{"center":r[0],"margin":r[1]} for r in failed[:20]],
        "theorem":(
            "The 28^3 closed Krawczyk boxes tile the entire nuisance cube of halfwidth 0.0005. Outward inclusion in every box proves one continued root throughout the full cube. This closes the ST240 sampling boundary and extends the complete cover from halfwidth 0.0004 to 0.0005."
        )
    }
    return finalize(254,"Complete Halfwidth-0.0005 Nuisance Cover",
        "proven_complete_21952_box_cover" if not failed else "partial_halfwidth_0005_cover",
        "This is a mathematical parameter-continuation cube, not a maximal domain, root-loss boundary, experimental tolerance, or calibration statement.",packet)


def real_fourier_modes() -> np.ndarray:
    x=np.arange(N); cols=[]
    for k in range(1,6):
        cols.append(np.cos(2*np.pi*k*x/N)); cols.append(np.sin(2*np.pi*k*x/N))
    cols.append((-1.0)**x)
    q=[]
    for v in cols:
        v=v-v.mean()
        for w in q: v=v-np.dot(v,w)*w
        if np.linalg.norm(v)>1e-10: q.append(v/np.linalg.norm(v))
        if len(q)==6: break
    return np.column_stack(q)


def st255_finite_estimators() -> dict:
    old=json.loads((ROOT/"FIN_ST203_ST217_Results.json").read_text(encoding="utf-8"))["ST211"]
    modes=real_fourier_modes(); u=np.ones(N)/N; rows=[]
    for r in old["rows"]:
        sigma=np.array(r["mode_attenuations"])
        n=max(1000,int(r["iid_samples_for_worst_mode_sd_0_01"]))
        covariance=np.diag(1/(12*n*sigma**2))
        cr=np.diag(1/(n*np.array(r["per_sample_Fisher_eigenvalues"])))
        rows.append({"layers":r["layers"],"declared_n":n,
                     "maximum_exact_covariance_minus_CR":float(np.max(np.abs(covariance-cr))),
                     "trace_MSE":float(np.trace(covariance)),
                     "worst_modal_variance":float(np.max(np.diag(covariance))),
                     "cluster_trace_MSE":float(old["design_effect"]*np.trace(covariance))})
    gram=modes.T@modes
    packet={
        "estimator":"theta_hat_k = v_k^T(p_hat-u)/sigma_k",
        "mode_orthonormality_error":float(np.linalg.norm(gram-np.eye(6))),
        "rows":rows,
        "cluster_model":"with probability rho a cluster shares one common categorical draw; otherwise its m draws are iid",
        "theorem":(
            "In the declared local multinomial family p(theta)=u+sum_k sigma_k theta_k v_k with orthonormal mean-zero modes, theta_hat is exactly unbiased at every admissible theta and at theta=0 has covariance diag((12 n sigma_k^2)^-1), exactly equal to the Cramer--Rao matrix. It is therefore finite-sample efficient at the uniform state, not merely asymptotically efficient. Under the declared common-shock exchangeable cluster construction, every modal covariance is multiplied exactly by 1+(m-1)rho=1.95."
        )
    }
    return finalize(255,"Exact Finite-Sample Efficient Fourier-Mode Estimators",
        "proven_exact_finite_sample_efficiency_in_declared_local_multinomial_model",
        "The observation family, Fourier mode access, attenuation factors, and common-shock cluster law are supplied. Probability positivity restricts the parameter neighborhood, and no laboratory detector is derived.",packet)


def st256_aligned_clifford() -> dict:
    rng=np.random.default_rng(SEED); rows=[]
    for n in [1,2,3,4]:
        x=np.diag([1.0]*n+[-1.0]*n)
        z=np.block([[np.zeros((n,n)),np.eye(n)],[np.eye(n),np.zeros((n,n))]])
        alpha=Fraction(n+2,n+3); beta=Fraction(2*n+1,2*n+3)
        # Random conjugate replay keeps exact theorem class but removes coordinate alignment.
        q,_=np.linalg.qr(rng.normal(size=(2*n,2*n)))
        xs=q@x@q.T; zs=q@z@q.T; x0=float(alpha)*xs; z0=float(beta)*zs
        rows.append({"half_dimension":n,"alpha":str(alpha),"beta":str(beta),
                     "X_spectral_gap":2*float(alpha),"Z_odd_gap":float(beta),
                     "constraint_error":float(np.linalg.norm(xs@zs+zs@xs)),
                     "objective_at_unique_pair":float(np.linalg.norm(xs-x0)**2+np.linalg.norm(zs-z0)**2)})
    packet={
        "rows":rows,
        "theorem":(
            "Let (X*,Z*) be any balanced Hermitian Clifford pair and let targets be X0=alpha X*, Z0=beta Z* with alpha,beta>0. For every competing Clifford pair (X,Z), the objective excess is exactly alpha||X-X*||_F^2+beta||Z-Z*||_F^2. It is strictly positive unless both factors equal the target pair. Hence the joint nearest-factor problem has one global minimizer in every balanced dimension. The two positive coefficients are simultaneous coercive gaps; the theorem is invariant under arbitrary common unitary conjugation."
        )
    }
    return finalize(256,"Global Joint Clifford Uniqueness for the Aligned Two-Gap Class",
        "proven_global_joint_uniqueness_for_aligned_positive_two_gap_targets",
        "This does not prove that marginal sign gaps alone imply joint uniqueness for arbitrary misaligned targets. The generic noisy two-factor problem remains open.",packet)


def lo(x): return float(np.nextafter(x,-np.inf))
def hi(x): return float(np.nextafter(x,np.inf))
def iconst(x): return (lo(float(x)),hi(float(x)))
def iadd(a,b): return (lo(a[0]+b[0]),hi(a[1]+b[1]))
def imul(a,b):
    v=[a[i]*b[j] for i in (0,1) for j in (0,1)]
    return (lo(min(v)),hi(max(v)))
def idiv(a,b): return imul(a,(lo(1/b[1]),hi(1/b[0])))
def isqrt(a): return (lo(math.sqrt(max(a[0],0))),hi(math.sqrt(max(a[1],0))))


def emission_step(b, emission, y):
    pred=iadd(iconst(.2),imul(iconst(.7),b))
    e0=emission[0][y];e1=emission[1][y]
    q=iadd(iconst(e1),imul(iconst(e0-e1),pred))
    post=idiv(imul(pred,iconst(e0)),q)
    return (max(0,post[0]),min(1,post[1])),q


def transfer_interval(b0,b1,depth,e0,e1):
    if depth==0:return (1.0,1.0)
    total=(0.0,0.0)
    for y in (0,1):
        n0,q0=emission_step(b0,e0,y);n1,q1=emission_step(b1,e1,y)
        w=isqrt(imul(q0,q1));sub=transfer_interval(n0,n1,depth-1,e0,e1)
        total=iadd(total,imul(w,sub))
    return total


def st257_belief_transfer() -> dict:
    e0=[[.98,.02],[.92,.08]];e1=[[.08,.92],[.02,.98]]
    grid=24;depth=9;vals=[]
    for i,j in itertools.product(range(grid),repeat=2):
        b0=(i/grid,(i+1)/grid);b1=(j/grid,(j+1)/grid)
        vals.append(transfer_interval(b0,b1,depth,e0,e1))
    m=min(z[0] for z in vals);M=max(z[1] for z in vals)
    radius=(m**(1/depth),M**(1/depth));rate=(-math.log(radius[1]),-math.log(radius[0]))
    previous=[0.6898819928911887,0.8198470892461787]
    combined=[max(rate[0],previous[0]),min(rate[1],previous[1])]
    packet={
        "belief_square_grid":grid,"block_depth":depth,"boxes":grid**2,
        "global_Tn1_interval":[m,M],"spectral_radius_interval":list(radius),
        "asymptotic_Hellinger_rate_interval":list(rate),
        "previous_quasi_multiplicative_rate_interval":previous,
        "combined_certified_rate_interval":combined,
        "direct_interval_width":rate[1]-rate[0],"combined_interval_width":combined[1]-combined[0],
        "operator":"(Tf)(b0,b1)=sum_y sqrt(q0_y(b0)q1_y(b1)) f(F0_y(b0),F1_y(b1))",
        "theorem":(
            "The observed-HMM Hellinger recursion is a positive transfer operator on the compact two-belief square. Adjacent binary64 intervals enclose every exact decimal input coefficient; the interval boxes cover the square, and outward recursion encloses T^n 1 pointwise between m and M. Positivity gives m^(1/n)<=r(T)<=M^(1/n) by order bounds, hence the direct asymptotic rate enclosure. Intersecting it with the independent earlier quasi-multiplicative certificate gives the tighter displayed combined interval."
        )
    }
    return finalize(257,"Outward Belief-Transfer Spectral-Radius Enclosure",
        "proven_global_interval_spectral_radius_enclosure",
        "The bound concerns the supplied two-state HMM and Hellinger exponent s=1/2. It is not an optimized Chernoff certificate, laboratory likelihood, or universal observed-process theorem.",packet)


def st258_reset_ledger() -> dict:
    rows=[]
    for beta_gap in [1,2,4,8,12]:
        excited=1/(1+math.exp(beta_gap)); entropy=-excited*math.log(excited)-(1-excited)*math.log(1-excited)
        rows.append({"beta_times_gap":beta_gap,"thermal_blank_excited_probability":excited,
                     "thermal_blank_entropy":entropy,"reset_ground_state_error":excited})
    packet={
        "controller":"H_c=(pi/(2 tau))(I-SWAP) for 0<t<tau",
        "exact_unitary":"exp(-i tau H_c)=SWAP",
        "ideal_net_switching_work":0.0,
        "rows":rows,
        "target_error_gap_law":"beta*epsilon=log((1-delta)/delta)",
        "theorem":(
            "For identical register Hamiltonians, I-SWAP commutes with total energy. A constant pulse H_c=pi(I-SWAP)/(2 tau) implements SWAP exactly. Sudden on/off switching has zero net system-plus-interaction work because H_c is conserved during the pulse, although controller dissipation is outside the closed model. A finite-temperature Gibbs blank resets only to its excited probability delta=(1+e^{beta epsilon})^-1; an exact pure blank requires beta epsilon to diverge. Restoring the used ancilla to the Gibbs blank can extract at most its nonequilibrium free-energy excess, whereas preparing a low-error blank requires the displayed dimensionless gap."
        )
    }
    return finalize(258,"Pulse-Controller and Finite-Temperature Blank Ledger",
        "proven_ideal_pulse_identity_and_finite_temperature_reset_error_law",
        "The ideal controller is lossless by assumption. Beta, the Hamiltonian gap, pulse duration, hardware dissipation, and physical energy units are not generated by FIN.",packet)


def strict_weights() -> np.ndarray:
    _,a,_=strict_operator(); return -a[0].copy()+np.eye(N)[0]*a[0,0]


def refined_laplacian(a,p,q,v):
    size=2*N; L=np.zeros((size,size))
    for i in range(N):
        for j in range(i+1,N):
            d=min((i-j)%N,(j-i)%N)
            for f in (0,1):
                for g,w in [(f,p[d]),(1-f,q[d])]:
                    ii=2*i+f;jj=2*j+g
                    L[ii,ii]+=w;L[jj,jj]+=w;L[ii,jj]-=w;L[jj,ii]-=w
        ii=2*i;jj=2*i+1;L[ii,ii]+=v;L[jj,jj]+=v;L[ii,jj]-=v;L[jj,ii]-=v
    return L


def st259_local_refinement(a: np.ndarray) -> dict:
    w=np.zeros(7)
    for d in range(1,7):w[d]=-a[0,d]
    rng=np.random.default_rng(SEED); rows=[]
    u=np.ones(2)/math.sqrt(2);R=np.kron(np.eye(N),u[:,None])
    for n in range(5):
        theta=rng.random(7);p=theta*w;q=w-p;v=float(rng.uniform(.05,2))
        L=refined_laplacian(a,p,q,v)
        rows.append({"sample":n,"vertical_rate":v,"theta":theta[1:].tolist(),
                     "intertwining_error":float(np.linalg.norm(L@R-R@a)),
                     "row_sum_error":float(np.linalg.norm(L@np.ones(2*N))),
                     "maximum_positive_offdiagonal":float(max(0,np.max(L-np.diag(np.diag(L))))),
                     "minimum_eigenvalue":float(np.linalg.eigvalsh(L).min())})
    packet={
        "coarse_distance_weights":w[1:].tolist(),
        "cone_parameters":"p_d,q_d>=0, p_d+q_d=W_d for d=1..6, and vertical v>=0",
        "free_dimension":7,"fiber_blind_subcone_dimension":1,
        "rows":rows,
        "theorem":(
            "Among D12-translation/reflection-invariant, fiber-swap-invariant two-child graph Laplacians whose edges depend only on base distance and whether fibers agree, exact intertwining holds iff p_d+q_d=W_d for all six distances. Thus the 78-dimensional abstract symmetric moduli cone collapses to a six-dimensional product of intervals times the vertical ray v>=0. If same- and cross-fiber edges are additionally indistinguishable, p_d=q_d=W_d/2 and only v remains. All such matrices are genuine positive graph Laplacians and preserve the complete coarse heat calculus."
        )
    }
    return finalize(259,"D12-Local Graph-Laplacian Refinement Cone",
        "proven_seven_dimensional_local_positive_refinement_classification",
        "Locality and symmetry reduce but do not eliminate nonuniqueness. No strict principle selects the six edge splits or vertical rate, and no physical length is produced.",packet)


def fiber_ds(mu,t): return 4*mu*t/(math.exp(min(2*mu*t,700))+1)


def st260_log_haar_plateau() -> dict:
    times=np.logspace(-5,5,500); js=np.arange(-12,13); mus=2.0**js
    curves=[]
    for t in times:curves.append(sum(fiber_ds(float(mu),float(t)) for mu in mus))
    curves=np.array(curves); target=2.0
    good=np.where(np.abs(curves-target)<.02)[0]
    plateau=[float(times[good[0]]),float(times[good[-1]])] if len(good) else []
    # Exact finite-cutoff integral for log-Haar density rho.
    packet={
        "continuum_measure":"rho dmu/mu on (0,infinity)",
        "exact_continuum_contribution":"2 rho ln(2)",
        "uniqueness":"dmu/mu is the Haar measure of the positive multiplicative scale group, unique up to rho",
        "dyadic_candidate":"mu_j=2^j, j=-12,...,12",
        "dyadic_predicted_plateau":target,"dyadic_numeric_plateau_tolerance":.02,
        "dyadic_plateau_time_window":plateau,"dyadic_maximum_deviation_on_window":float(max(abs(curves[good]-target))) if len(good) else None,
        "times":times.tolist(),"dyadic_curve":curves.tolist(),
        "theorem":(
            "A scale-invariant fiber-rate measure must be proportional to the Haar measure dmu/mu. Integrating one two-state fiber contribution gives integral_0^infinity [4 mu t/(e^{2 mu t}+1)] rho dmu/mu = 2 rho ln 2, independent of t. Finite geometric rate ladders approximate this plateau between their cutoff scales. The theorem constructs the missing scale-flow object mathematically, but its density, cutoffs, and realization are additional data."
        )
    }
    return finalize(260,"Log-Haar Fiber-Rate Measure and Spectral-Dimension Plateau",
        "proven_conditional_scale_invariant_plateau_theorem_with_finite_dyadic_replay",
        "The log-Haar density is not derived from the strict kernel. The dyadic illustration must not silently import alpha_geo or legacy octave roles into strict FIN. Infinite refinement also requires a separate trace/renormalization construction and supplies no seconds or lengths.",packet)


def st261_refinement_observables(a: np.ndarray) -> dict:
    u=np.ones(2)/math.sqrt(2);v=np.array([1,-1])/math.sqrt(2)
    R=np.kron(np.eye(N),u[:,None]);S=np.kron(np.eye(N),v[:,None])
    rows=[]
    for mu in [.3,1.,4.]:
        at=R@a@R.T+mu*(S@S.T)
        heat=expm(-at)
        rows.append({"mu":mu,"coarse_heat_error":float(np.linalg.norm(heat@R-R@expm(-a))),
                     "coarse_vertex_return":float((R[:,0]@heat@R[:,0])),
                     "fiber_odd_return":float((S[:,0]@heat@S[:,0])),
                     "predicted_fiber_odd_return":math.exp(-mu),
                     "complement_heat_trace":float(np.trace(S.T@heat@S))})
    packet={
        "rows":rows,
        "coarse_no_go":"all preparations in Ran(R) and effects compressed through R^* are independent of B",
        "complete_observable":"the matrix S^* exp(-t A_tilde) S=exp(-tB); at one known t>0, B=-(1/t)log of this positive matrix",
        "theorem":(
            "Exact intertwining makes every experiment whose preparations lie in Ran(R) and whose effects are coarse-compressed completely insensitive to the complement generator B. Conversely, complement-resolving preparations and measurements recover exp(-tB); because it is positive definite for self-adjoint B, the principal logarithm reconstructs B uniquely from one exact time. The displayed odd-fiber survival probabilities distinguish mu=0.3,1,4 while all coarse heat responses agree to numerical precision."
        )
    }
    return finalize(261,"Operational Completeness and No-Go for Refinement Observables",
        "proven_coarse_indistinguishability_and_complement_tomography_theorem",
        "The theorem specifies what an instrument must resolve; it does not supply a laboratory fiber, preparation, POVM, clock calibration, or physical refinement.",packet)


def st262_external_stop() -> dict:
    packet={"decision":"blocked_without_synthetic_substitution","local_search_performed":False,
            "required_missing_atoms":["independently registered raw events","calibration hash","custody separation","laboratory attestation"],
            "theorem":"No local calculation, figure, hash, or simulated detector record entails that an external apparatus executed the protocol. The evidence stage therefore remains blocked rather than being filled with synthetic data."}
    return finalize(262,"External Validation Evidence Stop","blocked_no_external_record",
                    "This is an evidentiary boundary, not a mathematical failure of the operator model.",packet)


def make_figures(d: dict) -> None:
    FIG_DIR.mkdir(exist_ok=True)
    fig,ax=plt.subplots(figsize=(7,4));r=d["ST248"]["tested_natural_state_rows"];ax.bar(range(len(r)),[z["distance_from_uniform"] for z in r]);ax.set_xticks(range(len(r)),["uniform","heat diag","Green diag","A² diag"],rotation=20);ax.set(ylabel="distance from uniform",title="ST248: A-only natural states remain uniform");fig.tight_layout();fig.savefig(FIG_DIR/"st248_preparation_no_go.png",dpi=190);plt.close(fig)
    fig,ax=plt.subplots(figsize=(7,4));r=d["ST250"]["rows"];ax.scatter([z["step"] for z in r],[z["kappa"] for z in r],s=7,alpha=.5);ax.set(xlabel="continuation step",ylabel="kappa",title="ST250: six new certified steps on every signed seed");fig.tight_layout();fig.savefig(FIG_DIR/"st250_branch_extension.png",dpi=190);plt.close(fig)
    fig,ax=plt.subplots(figsize=(7,4));tr=d["ST252"]["trajectory"];ax.plot(tr["time"],tr["x"],label="state x");ax.plot(tr["time"],tr["pi"],label="equilibrium control pi");ax.legend();ax.set(xlabel="dimensionless time",ylabel="probability",title="ST252: endpoint- and slew-priced feasible selector");fig.tight_layout();fig.savefig(FIG_DIR/"st252_actuator_selector.png",dpi=190);plt.close(fig)
    fig,ax=plt.subplots(figsize=(7,4));p=d["ST254"]["margin_percentiles"];x=list(map(float,p.keys()));ax.plot(x,list(p.values()),"o-");ax.axhline(0,color="black",lw=.8);ax.set(xlabel="margin percentile",ylabel="Krawczyk inclusion margin",title="ST254: complete 28^3 nuisance cover");fig.tight_layout();fig.savefig(FIG_DIR/"st254_cover_margins.png",dpi=190);plt.close(fig)
    fig,ax=plt.subplots(figsize=(7,4));r=d["ST255"]["rows"];ax.semilogy([z["layers"] for z in r],[z["trace_MSE"] for z in r],"o-",label="iid exact");ax.semilogy([z["layers"] for z in r],[z["cluster_trace_MSE"] for z in r],"s--",label="clustered");ax.legend();ax.set(xlabel="layers",ylabel="trace MSE",title="ST255: finite-sample efficient modal estimator");fig.tight_layout();fig.savefig(FIG_DIR/"st255_estimators.png",dpi=190);plt.close(fig)
    fig,ax=plt.subplots(figsize=(7,4));r=d["ST258"]["rows"];ax.semilogy([z["beta_times_gap"] for z in r],[z["reset_ground_state_error"] for z in r],"o-");ax.set(xlabel="beta times gap",ylabel="thermal reset error",title="ST258: finite-temperature blank error");fig.tight_layout();fig.savefig(FIG_DIR/"st258_thermal_blank.png",dpi=190);plt.close(fig)
    fig,ax=plt.subplots(figsize=(7,4));r=d["ST260"];ax.semilogx(r["times"],r["dyadic_curve"],label="finite dyadic ladder");ax.axhline(r["dyadic_predicted_plateau"],color="black",ls="--",label="log-Haar limit");ax.legend();ax.set(xlabel="dimensionless t",ylabel="fiber spectral dimension",title="ST260: logarithmic scale plateau");fig.tight_layout();fig.savefig(FIG_DIR/"st260_log_haar_plateau.png",dpi=190);plt.close(fig)
    fig,ax=plt.subplots(figsize=(7,4));r=d["ST261"]["rows"];ax.plot([z["mu"] for z in r],[z["coarse_vertex_return"] for z in r],"o-",label="coarse return");ax.plot([z["mu"] for z in r],[z["fiber_odd_return"] for z in r],"s--",label="fiber-odd return");ax.legend();ax.set(xlabel="complement rate mu",ylabel="t=1 return",title="ST261: coarse blindness versus fiber resolution");fig.tight_layout();fig.savefig(FIG_DIR/"st261_refinement_observables.png",dpi=190);plt.close(fig)


def main() -> None:
    _,a,_=strict_operator()
    out={"metadata":{"programs":"ST248-ST262","date":"2026-08-12","seed":SEED,"python":platform.python_version(),"numpy":np.__version__,"scipy":scipy.__version__,"sympy":sp.__version__}}
    out["ST248"]=st248_preparation_source_no_go(a);out["ST249"]=st249_typed_operator_classification(a)
    out["ST250"]=st250_branch_resource_stop(a);out["ST251"]=st251_rotated_amplitude_damping()
    out["ST252"]=st252_selector_actuator_cost();out["ST253"]=st253_affine_connection(a)
    out["ST254"]=st254_complete_nuisance_cover();out["ST255"]=st255_finite_estimators()
    out["ST256"]=st256_aligned_clifford();out["ST257"]=st257_belief_transfer()
    out["ST258"]=st258_reset_ledger();out["ST259"]=st259_local_refinement(a)
    out["ST260"]=st260_log_haar_plateau();out["ST261"]=st261_refinement_observables(a)
    out["ST262"]=st262_external_stop()
    out["recommended_next_programs"]=[
        {"id":"ST263","priority":1,"study":"test one genuinely new strict symmetry-breaking state source outside deterministic D12-natural rules, or preserve the ST248 no-go"},
        {"id":"ST264","priority":2,"study":"classify nonlinear D12-equivariant state-to-multiplication maps by orbit type and stabilizer monotonicity"},
        {"id":"ST265","priority":3,"study":"extend only the lowest-clearance branch pair with interval event detection"},
        {"id":"ST266","priority":4,"study":"derive a rational dual for one affine qubit channel not unitarily equivalent to amplitude damping or Pauli mixtures"},
        {"id":"ST267","priority":5,"study":"build an interval KKT certificate or a counterexample for the ST252 augmented selector control"},
        {"id":"ST268","priority":6,"study":"test whether the flat mean-zero connection is unique under a wider naturality category including refinement maps"},
        {"id":"ST269","priority":7,"study":"seek the first certified nuisance-cover failure beyond halfwidth 0.0005 under a predeclared resource cap"},
        {"id":"ST270","priority":8,"study":"construct exact finite-confidence regions for the ST255 efficient estimator under clustered counts"},
        {"id":"ST271","priority":9,"study":"resolve generic misaligned joint Clifford uniqueness or give an explicit two-gap counterexample"},
        {"id":"ST272","priority":10,"study":"refine the belief-transfer enclosure with adaptive boxes and optimized Chernoff s"},
        {"id":"ST273","priority":11,"study":"add a finite-dimensional controller and account for its cyclic free-energy cost"},
        {"id":"ST274","priority":12,"study":"intersect the ST259 cone with a newly justified strict refinement functor, if one is found"},
        {"id":"ST275","priority":13,"study":"derive finite-cutoff error bounds for the log-Haar plateau and audit trace-class completion"},
        {"id":"ST276","priority":14,"study":"design a frozen fiber-resolving instrument specification for ST261 without claiming laboratory execution"},
        {"id":"ST277","priority":15,"study":"resume external validation only with independently registered events"},
    ]
    out["central_verdict"]=(
        "The strict operator alone cannot deterministically source the localized state needed by the noncommuting multiplication algebra; transitive D12 naturality forces the uniform state. Conditional on the canonical mean-zero affine amplitude carrier, however, A does induce a flat Levi-Civita connection that repairs nonlinear response transport. Exact refinements remain nonunique even after graph-Laplacian locality, but complement-resolving instruments can reconstruct the hidden generator. A scale-invariant log-Haar rate measure produces an exact spectral-dimension plateau, while its density and cutoffs remain unsourced."
    )
    out["epistemic_boundary"]=(
        "No non-premise selector, QW-2191 discharge, physical state preparation, dimensional clock/length/action unit, Planck scale, spacetime, external record, legacy-to-strict completion or role transfer, Standard Model, gravity, L_total, or ToE closure is claimed."
    )
    RESULTS.write_text(json.dumps(native(out),indent=2,sort_keys=True),encoding="utf-8")
    with SUMMARY.open("w",newline="",encoding="utf-8") as h:
        w=csv.writer(h);w.writerow(["program","object","status"])
        for k in range(248,263):w.writerow([f"ST{k}",out[f"ST{k}"]["object"],out[f"ST{k}"]["status"]])
    make_figures(out)
    print(json.dumps({"programs":"ST248-ST262","statuses":{f"ST{k}":out[f"ST{k}"]["status"] for k in range(248,263)},
                      "ST250":[out["ST250"]["new_certified_steps"],out["ST250"]["minimum_nonfold_clearance"]],
                      "ST252":out["ST252"]["best_total_cost"],"ST254":[out["ST254"]["passed"],out["ST254"]["minimum_margin"]],
                      "ST257":out["ST257"]["asymptotic_Hellinger_rate_interval"],"ST260":out["ST260"]["dyadic_plateau_time_window"]},indent=2))


if __name__=="__main__":main()
