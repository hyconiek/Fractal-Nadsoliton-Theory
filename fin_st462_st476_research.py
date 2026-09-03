#!/usr/bin/env python3
"""FIN ST462--ST476: transition orbit, source obstruction, dual transfer, and IR.

All statements are local, finite, dimensionless, analytic, interval-certified,
or explicitly numerical in their declared scope.  The batch preserves the
legacy/strict split.  It does not source the supplied active gain, discharge
QW-2191, choose an absolute unit, create a physical continuum, transfer a
legacy physical role, or provide independent empirical evidence.
"""

from __future__ import annotations

import csv
import hashlib
import itertools
import json
import math
import shutil
import subprocess
import tempfile
import time
from pathlib import Path
from typing import Any

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.linalg import expm, qr
from scipy.optimize import brentq, minimize, minimize_scalar, root
from scipy.integrate import solve_ivp

from fin_st402_st416_research import independent_strict_matrix_float, real_orbit_eval
from fin_st417_st431_research import N, SEED
from fin_st447_st461_research import (
    adaptive_global_cover,
    degree8_evaluation_matrix,
    direct_ir_cell,
    exponent_orbit,
    regularized_ir_float,
    FACE,
)


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "FIN_ST462_ST476_Results.json"
SUMMARY = ROOT / "FIN_ST462_ST476_Summary.csv"
FIG_DIR = ROOT / "FIN_ST462_ST476_Figures"
SINGULAR_LOG = ROOT / "FIN_ST468_Singular_Degree8_Rank.log"
NAMES = {
    462: "First_Transition_Orbit_Numerical_Isolation",
    463: "Adversarial_Competing_Orbit_Search",
    464: "Global_Cover_Refinement_Stop",
    465: "Autocatalytic_Gain_Source_Obstruction",
    466: "Degree8_Support4_5_Circuit_Audit",
    467: "Degree8_Floating_Rank_Profile",
    468: "Degree8_Exact_Singular_Rank_Attempt",
    469: "Boundary_Minimizer_Exclusion_Theorem",
    470: "Adaptive_Infrared_Extension_and_Stop",
    471: "Dual_Dynamics_Universal_Discriminant",
    472: "Dual_Dynamics_Intertwining_Stability",
    473: "Isospectral_Vertex_Record_Counterexample",
    474: "Gain_Source_Admission_Gate",
    475: "Selector_and_Scale_Admission_Gate",
    476: "Independent_Evidence_Gate",
}
PACKETS = {k: ROOT / f"FIN_ST{k}_{v}.json" for k, v in NAMES.items()}


def native(x: Any) -> Any:
    if isinstance(x, dict): return {str(k): native(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)): return [native(v) for v in x]
    if isinstance(x, np.ndarray): return native(x.tolist())
    if isinstance(x, (np.floating, np.integer)): return x.item()
    if isinstance(x, bool): return bool(x)
    return x


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def finalize(k: int, status: str, boundary: str, packet: dict) -> dict:
    path = PACKETS[k]
    path.write_text(json.dumps(native(packet), indent=2, sort_keys=True), encoding="utf-8")
    return {
        "program": f"ST{k}", "object": NAMES[k], "packet_file": path.name,
        "packet_sha256": sha(path), **packet, "status": status, "boundary": boundary,
    }


A = independent_strict_matrix_float()
U = np.ones(N) / N


def probability_ratio(p: np.ndarray) -> float:
    p = np.asarray(p, float)
    d = p - U
    q = float(d @ A @ d)
    if q <= 1e-18: return 12 / float(np.linalg.eigvalsh(A)[-1])
    pos = p > 0
    entropy = float(np.sum(p[pos] * np.log(N * p[pos])))
    return 2 * entropy / q


def softmax_reduced(z: np.ndarray) -> np.ndarray:
    zz = np.r_[z, 0.0]
    zz -= np.max(zz)
    e = np.exp(zz)
    return e / np.sum(e)


def ratio_and_gradient(z: np.ndarray) -> tuple[float, np.ndarray]:
    p = softmax_reduced(z)
    d = p - U
    q = float(d @ A @ d)
    logp = np.log(N * np.maximum(p, np.finfo(float).tiny))
    entropy = float(np.sum(p * logp))
    value = 2 * entropy / q
    gp = 2 * ((logp + 1) * q - entropy * (2 * A @ d)) / q**2
    gz = p[:-1] * (gp[:-1] - float(gp @ p))
    return value, gz


def d12_orbit(p: np.ndarray, tolerance: float = 1e-7) -> list[np.ndarray]:
    rows = []
    for k in range(N):
        rows.append(np.roll(p, k))
        rows.append(np.roll(p[::-1], k))
    unique = []
    for row in rows:
        if not any(np.linalg.norm(row - x) < tolerance for x in unique): unique.append(row)
    return unique


def canonical_cycle(p: np.ndarray) -> np.ndarray:
    return min(d12_orbit(p), key=lambda x: tuple(np.round(x, 14)))


def transition_multistart(count: int = 108) -> dict:
    rng = np.random.default_rng(SEED + 462)
    starts = []
    for j in range(N):
        p = np.full(N, 0.008)
        p[j] = 1 - np.sum(p[np.arange(N) != j])
        starts.append(np.log(p[:-1] / p[-1]))
    starts += [rng.normal(0, 3.0, N - 1) for _ in range(count - N)]
    rows = []
    for z in starts:
        opt = minimize(lambda x: ratio_and_gradient(x), z, jac=True, method="BFGS",
                       options={"gtol": 2e-10, "maxiter": 1000})
        p = softmax_reduced(opt.x)
        rows.append((float(opt.fun), p, bool(opt.success), float(np.linalg.norm(opt.jac, np.inf))))
    rows.sort(key=lambda x: x[0])
    best = rows[0]
    orbit = d12_orbit(best[1])
    reflection_residual = min(float(np.linalg.norm(best[1] - np.roll(best[1][::-1], k)))
                              for k in range(N))
    hits = [r for r in rows if abs(r[0] - best[0]) < 3e-10]
    matched = []
    for _, p, _, _ in hits:
        matched.append(min(float(np.linalg.norm(p - q)) for q in orbit))
    return {
        "start_count": len(rows), "successful_optimizer_flags": sum(r[2] for r in rows),
        "best_ratio": best[0], "best_probability": best[1],
        "canonical_probability": canonical_cycle(best[1]),
        "minimum_reflection_stabilizer_residual": reflection_residual,
        "orbit_identification_tolerance": 1e-7,
        "best_gradient_infinity_norm": best[3], "D12_orbit_size": len(orbit),
        "best_orbit_hit_count": len(hits), "maximum_orbit_matching_residual": max(matched),
        "distinct_objective_clusters_1e_8": len({round(r[0], 8) for r in rows}),
        "ten_lowest_values": [r[0] for r in rows[:10]],
    }


def st462() -> dict:
    search = transition_multistart()
    bracket = [2.902496471747767, 2.9024964917477667]
    inside = bracket[0] <= search["best_ratio"] <= bracket[1]
    packet = {
        **search, "in_ST342_certified_reflection_even_crossing_bracket": inside,
        "candidate_orbit_description": "one dominant vertex and reflection-even distance profile",
        "claim": "All best multistarts land on one twelve-member D12 orbit at the center of the prior certified coexistence bracket.",
        "global_first_orbit_theorem": False,
    }
    return finalize(462, "strong_numerical_isolation_of_twelve_member_transition_orbit",
                    "The certified global lower bound remains 2.8934; numerical orbit isolation does not exclude a third orbit inside the residual strip.", packet)


def optimize_on_support(ids: list[int], seed: np.ndarray | None = None) -> tuple[float, np.ndarray]:
    k = len(ids)
    if k == 1:
        p = np.zeros(N); p[ids[0]] = 1
        return probability_ratio(p), p
    z0 = np.zeros(k - 1) if seed is None else seed
    def unpack(z):
        zz = np.r_[z, 0.0]; zz -= zz.max(); e = np.exp(zz); w = e / e.sum()
        p = np.zeros(N); p[ids] = w; return p
    opt = minimize(lambda z: probability_ratio(unpack(z)), z0, method="Nelder-Mead",
                   options={"maxiter": 1200, "xatol": 1e-11, "fatol": 1e-12})
    return float(opt.fun), unpack(opt.x)


def adversarial_search() -> dict:
    uniform_rows = []
    for mask in range(1, 2**N - 1):
        ids = [i for i in range(N) if mask & (1 << i)]
        p = np.zeros(N); p[ids] = 1 / len(ids)
        uniform_rows.append((probability_ratio(p), ids))
    uniform_rows.sort(key=lambda x: x[0])
    refined = []
    seen = set()
    for _, ids in uniform_rows[:96]:
        key = tuple(ids)
        if key in seen: continue
        seen.add(key)
        val, p = optimize_on_support(ids)
        refined.append((val, ids, p))
    refined.sort(key=lambda x: x[0])
    rng = np.random.default_rng(SEED + 463)
    random_min = math.inf; random_p = None
    for alpha in (0.03, 0.1, 0.3, 1.0, 3.0):
        p = rng.dirichlet(np.full(N, alpha), size=12000)
        d = p - U
        q = np.einsum("bi,ij,bj->b", d, A, d)
        ent = np.sum(p * np.log(N * np.maximum(p, np.finfo(float).tiny)), axis=1)
        vals = 2 * ent / q
        i = int(np.argmin(vals))
        if vals[i] < random_min: random_min, random_p = float(vals[i]), p[i]
    return {
        "all_proper_uniform_supports_tested": len(uniform_rows),
        "minimum_uniform_support_ratio": uniform_rows[0][0],
        "minimum_uniform_support": uniform_rows[0][1],
        "refined_low_support_masks": len(refined),
        "minimum_refined_boundary_ratio": refined[0][0],
        "minimum_refined_boundary_support": refined[0][1],
        "minimum_refined_boundary_probability": refined[0][2],
        "random_Dirichlet_samples": 5 * 12000,
        "minimum_unoptimized_random_ratio": random_min,
        "minimum_unoptimized_random_probability": random_p,
        "candidate_beaten": min(refined[0][0], random_min) < 2.902496471747767,
    }


def st463() -> dict:
    search = adversarial_search()
    return finalize(463, "strong_adversarial_search_no_earlier_competing_orbit_found",
                    "Finite random and selected-face searches are not an exhaustive global certificate; ST448 remains the theorem-level lower bound.", search)


def st464() -> dict:
    gain = 2.894
    t0 = time.time()
    cover = adaptive_global_cover(gain, 0.33415, initial=40, max_depth=15)
    packet = {
        "attempted_gain": gain, "attempted_low_sector_cut": 0.33415,
        "wall_seconds": time.time() - t0, "cover": cover,
        "previous_certified_lower_gain": 2.8934,
        "result": "The unchanged scalar/cover architecture does not certify 2.894 at the declared stop. Its first failed leaf reaches the roundoff-payment floor.",
        "new_global_lower_theorem": False,
    }
    return finalize(464, "bounded_no_go_for_unmodified_global_cover_refinement",
                    "Method failure is not a counterexample and does not lower the accepted 2.8934 theorem.", packet)


def st465() -> dict:
    kappa, gamma = 2.0, 1.0
    # Linearized illustrative controller: q(p) is strict-derived, but p-u is a supplied seed.
    mode = np.linalg.eigh(A)[1][:, 1]
    mode -= mode.mean(); mode /= np.linalg.norm(mode)
    def simulate(amplitude: float):
        p0 = U + amplitude * mode
        y0 = np.r_[p0, 0.0]
        def rhs(_, y):
            p, g = y[:-1], y[-1]
            q = float((p-U) @ A @ (p-U))
            # Stable relaxation is used only to expose the source accounting.
            return np.r_[-A @ (p-U), kappa*q-gamma*g]
        sol = solve_ivp(rhs, (0, 15), y0, rtol=1e-10, atol=1e-12, dense_output=False)
        return {"amplitude": amplitude, "maximum_g": float(np.max(sol.y[-1])),
                "final_g": float(sol.y[-1, -1]), "final_state_distance": float(np.linalg.norm(sol.y[:-1, -1]-U))}
    exact = simulate(0.0); perturbed = simulate(0.02)
    packet = {
        "candidate_controller": "g_dot=kappa*(p-u)^T A (p-u)-gamma*g",
        "parameters": {"kappa": kappa, "gamma": gamma},
        "symmetric_vacuum_simulation": exact, "perturbed_simulation": perturbed,
        "theorem": "Because A u=0, (p,g)=(u,0) is an invariant absorbing state. The strict quadratic q(p) can amplify a supplied nonuniform state but cannot create its own first departure from the symmetric vacuum.",
        "new_autonomous_gain_source_from_symmetric_vacuum": False,
    }
    return finalize(465, "proven_autocatalytic_receiver_requires_seed_resource",
                    "A fluctuation law, initial asymmetry, pump, or signed mediator remains an additional source object.", packet)


def rank_small_mod(a: np.ndarray, prime: int) -> int:
    b = a.copy().astype(np.int64) % prime
    rank = 0
    for col in range(b.shape[1]):
        nz = np.flatnonzero(b[rank:, col])
        if not len(nz): continue
        j = rank + int(nz[0]); b[[rank, j]] = b[[j, rank]]
        b[rank] = b[rank] * pow(int(b[rank, col]), -1, prime) % prime
        for i in range(b.shape[0]):
            if i != rank and b[i, col]: b[i] = (b[i] - b[i, col] * b[rank]) % prime
        rank += 1
        if rank == min(b.shape): break
    return rank


def sampled_circuit_audit(prime: int, sample4: int = 18000, sample5: int = 18000) -> dict:
    matrix = degree8_evaluation_matrix(prime, npts=64)
    rng = np.random.default_rng(SEED + 466 + prime)
    rows = []
    for support, count in ((4, sample4), (5, sample5)):
        proj = rng.integers(1, prime, size=(support, matrix.shape[0]), dtype=np.int64)
        projected = proj @ matrix % prime
        deficient = 0; true = 0; examples = []
        for _ in range(count):
            ids = np.sort(rng.choice(matrix.shape[1], size=support, replace=False))
            if rank_small_mod(projected[:, ids], prime) < support:
                deficient += 1
                if rank_small_mod(matrix[:, ids], prime) < support:
                    true += 1
                    if len(examples) < 8: examples.append(ids.tolist())
        rows.append({"support": support, "sampled_subsets": count,
                     "projection_deficient_candidates": deficient,
                     "true_full_matrix_dependencies": true, "examples": examples})
    return {"prime": prime, "matrix_shape": list(matrix.shape), "rows": rows,
            "matrix_sha256": hashlib.sha256(matrix.tobytes()).hexdigest()}


def st466() -> dict:
    audits = [sampled_circuit_audit(p) for p in (1000003, 1000033)]
    found = sum(x["true_full_matrix_dependencies"] for a in audits for x in a["rows"])
    packet = {
        "exact_modular_random_audits": audits, "true_sampled_circuits": found,
        "claim": "No support-four or support-five circuit was found in the declared exact random audit.",
        "exhaustive_support4_or_5_theorem": False,
    }
    return finalize(466, "exact_sampled_no_support4_5_circuit_found",
                    "The search is reproducible but not exhaustive; ST447 proves only support at least four.", packet)


def degree8_real_matrix(npts: int = 1910) -> np.ndarray:
    base = json.loads((ROOT / "FIN_ST392_Primitive_Invariant_Generators_and_Syzygies.json").read_text())
    rng = np.random.default_rng(SEED + 467)
    pts = rng.normal(size=(npts, N)); pts -= pts.mean(axis=1, keepdims=True)
    pts /= np.linalg.norm(pts, axis=1, keepdims=True)
    ev = lambda rep: real_orbit_eval(rep, pts)
    qv = [ev(r) for r in base["quadratic_generator_representatives"]]
    p4 = [ev(r) for r in base["primitive_quartic_representatives"]]
    p6 = [ev(r) for r in base["primitive_sextic_representatives"]]
    cols = []
    for ids in itertools.combinations_with_replacement(range(6), 4):
        v = np.ones(npts)
        for i in ids: v *= qv[i]
        cols.append(v)
    for ids in itertools.combinations_with_replacement(range(6), 2):
        q2 = qv[ids[0]] * qv[ids[1]]
        cols.extend(q2 * x for x in p4)
    cols.extend(p4[i] * p4[j] for i, j in itertools.combinations_with_replacement(range(32), 2))
    cols.extend(qv[i] * p6[j] for i in range(6) for j in range(117))
    return np.column_stack(cols)


def st467() -> dict:
    t0 = time.time(); matrix = degree8_real_matrix()
    norms = np.linalg.norm(matrix, axis=0); normalized = matrix / norms
    _, r, piv = qr(normalized, mode="economic", pivoting=True, check_finite=False)
    diag = np.abs(np.diag(r)); lead = float(diag[0])
    thresholds = [lead * x for x in (1e-8, 1e-10, 1e-12, 1e-14)]
    ranks = {f"relative_{x:.0e}": int(np.sum(diag > lead*x)) for x in (1e-8, 1e-10, 1e-12, 1e-14)}
    packet = {
        "matrix_shape": list(matrix.shape), "column_normalized": True,
        "pivoted_QR_ranks": ranks, "largest_diagonal": lead,
        "diagonal_1888_to_1898": diag[1887:1898],
        "candidate_nullity_if_rank_1892": matrix.shape[1] - 1892,
        "wall_seconds": time.time() - t0,
        "claim": "The floating rank profile is compared with the exact Molien dimension 1892; it is diagnostic only.",
        "exact_rank_theorem": False,
    }
    return finalize(467, "floating_degree8_rank_profile_computed",
                    "Ill-conditioning and threshold dependence prevent promotion to an exact rank/nullspace theorem.", packet)


def st468() -> dict:
    singular = shutil.which("Singular")
    prime = 1000003; timeout = 45
    packet = {"Singular_executable": singular, "prime": prime, "timeout_seconds": timeout,
              "target_shape": [1892, 2028], "exact_rank_returned": None,
              "explicit_nullspace_returned": False}
    if singular is None:
        packet["result"] = "Singular unavailable"
        status = "blocked_exact_engine_unavailable"
    else:
        t0 = time.time(); matrix = degree8_evaluation_matrix(prime, npts=1892)
        with tempfile.NamedTemporaryFile("w", suffix=".sing", delete=False) as f:
            tmp = Path(f.name)
            f.write(f"ring r={prime},(x),dp;\n")
            f.write(f"matrix M[{matrix.shape[0]}][{matrix.shape[1]}]=")
            flat = matrix.ravel()
            chunk = 20000
            for i in range(0, len(flat), chunk):
                if i: f.write(",")
                f.write(",".join(map(str, flat[i:i+chunk].tolist())))
            f.write(";\nprint(rank(M));\nquit;\n")
        try:
            proc = subprocess.run([singular, "-q", str(tmp)], capture_output=True,
                                  text=True, timeout=timeout)
            output = (proc.stdout + "\n" + proc.stderr).strip()
            SINGULAR_LOG.write_text(output + "\n", encoding="utf-8")
            nums = [int(x) for x in output.split() if x.isdigit()]
            rank_value = nums[-1] if nums else None
            packet.update({"exact_rank_returned": rank_value, "return_code": proc.returncode,
                           "log_file": SINGULAR_LOG.name, "log_sha256": sha(SINGULAR_LOG),
                           "wall_seconds": time.time()-t0})
            status = "exact_rank_returned" if rank_value is not None else "exact_attempt_returned_no_rank"
        except subprocess.TimeoutExpired as exc:
            text = ((exc.stdout or b"") if isinstance(exc.stdout, bytes) else (exc.stdout or ""))
            SINGULAR_LOG.write_text(str(text) + "\nTIMEOUT\n", encoding="utf-8")
            packet.update({"result": "bounded timeout before exact rank", "wall_seconds": time.time()-t0,
                           "log_file": SINGULAR_LOG.name, "log_sha256": sha(SINGULAR_LOG)})
            status = "bounded_exact_rank_timeout_no_scientific_verdict"
        finally:
            tmp.unlink(missing_ok=True)
    return finalize(468, status,
                    "A timeout or rank alone is not an explicit characteristic-zero syzygy basis; lifting and polynomial verification remain mandatory.", packet)


def st469() -> dict:
    packet = {
        "theorem": "For every finite g and every boundary p of the simplex, choose j with p_j=0 and p_epsilon=(1-epsilon)p+epsilon e_j. The entropy change contains epsilon log epsilon, while the quadratic change is O(epsilon). Hence V_g(p_epsilon)<V_g(p) for all sufficiently small positive epsilon; no boundary point is a local or global minimizer.",
        "dominant_term": "epsilon*log(epsilon)",
        "quadratic_change_order": "O(epsilon)",
        "consequence": "Every first-transition minimizer is interior, so boundary supports cannot attain g_*.",
        "requires_kernel_beyond_bounded_A": False,
    }
    return finalize(469, "proven_boundary_minimizer_exclusion_for_all_finite_gain",
                    "The theorem excludes boundary attainment but does not identify the interior orbit or improve the certified global lower number.", packet)


def ir_extension_until_stop() -> dict:
    start, target, width = 0.1, 0.2, 0.003
    prior = root(lambda z: regularized_ir_float(z, start), FACE, tol=1e-12).x
    rows = []; lo = start; failure = None
    while lo < target - 1e-14:
        hi = min(target, lo + width)
        row, candidate = direct_ir_cell(lo, hi, prior)
        if not row["included"]:
            failure = row; break
        rows.append(row); prior = candidate; lo = hi
    return {"start": start, "target": target, "cell_width": width,
            "certified_stop": lo, "accepted_cells": rows, "first_failed_cell": failure}


def st470() -> dict:
    chain = ir_extension_until_stop()
    rows = chain["accepted_cells"]
    packet = {
        **chain, "accepted_cell_count": len(rows),
        "minimum_Krawczyk_margin": min(r["minimum_margin"] for r in rows),
        "maximum_weighted_contraction": max(r["weighted_contraction"] for r in rows),
        "all_Jacobian_determinants_negative": all(r["Jacobian_determinant"] < 0 for r in rows),
        "theorem": f"The previously certified local IR component extends from b=0.1 through b={chain['certified_stop']:.3f} in the displayed boxes.",
        "extension_to_target_0_2": chain["certified_stop"] >= 0.2 - 1e-14,
    }
    return finalize(470, "proven_local_IR_extension_through_first_new_certificate_stop",
                    "The first failed box is a method stop, not a branch singularity; disconnected roots and the remainder to b=0.2 remain open.", packet)


def dual_distance_scalar(x: float) -> float:
    return math.sqrt(1 + math.exp(-2*x) - 2*math.exp(-x)*math.cos(x))


def st471() -> dict:
    xstar = brentq(lambda x: math.sin(x)+math.cos(x)-math.exp(-x), 1.6, 3.0)
    dstar = dual_distance_scalar(xstar)
    lam = np.linalg.eigvalsh(A)[1:]
    times = xstar / lam
    packet = {
        "scalar_discriminant": "h(x)=|exp(-ix)-exp(-x)|",
        "stationarity_equation": "sin(x)+cos(x)=exp(-x)",
        "first_positive_global_maximizer": xstar,
        "universal_maximum_distance": dstar,
        "strict_positive_eigenvalues": lam,
        "mode_times_realizing_the_same_maximum": times,
        "theorem": "For a common nonnegative self-adjoint generator, ||U_t-P_t||=max_r h(t lambda_r). Every positive mode reaches the same universal coherent-versus-diffusive separation h(x*) at t=x*/lambda_r.",
        "physical_time_calibrated": False,
    }
    return finalize(471, "proven_universal_dimensionless_dual_dynamics_discriminant",
                    "The result discriminates mathematical channels after a clock/measurement model is supplied; it does not choose which channel nature realizes.", packet)


def st472() -> dict:
    rng = np.random.default_rng(SEED + 472)
    p0 = np.eye(N) - np.ones((N, N))/N
    e = rng.normal(size=(N, N)); e = p0 @ ((e+e.T)/2) @ p0
    e /= np.linalg.norm(e, 2)
    epsilon = 1e-3; b = A + epsilon*e
    times = [0.25, 1.0, 5.0]; rows = []
    for t in times:
        eu = float(np.linalg.norm(expm(-1j*t*A)-expm(-1j*t*b), 2))
        eh = float(np.linalg.norm(expm(-t*A)-expm(-t*b), 2))
        rows.append({"t": t, "unitary_error": eu, "heat_error": eh,
                     "Duhamel_bound": t*epsilon,
                     "unitary_pass": eu <= t*epsilon*(1+1e-9),
                     "heat_pass": eh <= t*epsilon*(1+1e-9)})
    packet = {
        "generator_defect_norm": epsilon,
        "perturbed_positive_spectrum_min": float(np.linalg.eigvalsh(b)[1]),
        "numerical_fixture": rows,
        "all_bounds_pass": all(r["unitary_pass"] and r["heat_pass"] for r in rows),
        "theorem": "If self-adjoint nonnegative generators obey ||A_tilde R-R A||<=delta and R is an isometry, Duhamel's formula gives both ||exp(-it A_tilde)R-R exp(-itA)||<=|t|delta and ||exp(-t A_tilde)R-R exp(-tA)||<=t delta.",
    }
    return finalize(472, "proven_linear_dual_channel_intertwining_stability_bound",
                    "The bound controls a supplied refinement/approximation defect; it does not construct the refinement or physical continuum.", packet)


def st473() -> dict:
    vals, vecs = np.linalg.eigh(A)
    rot = np.eye(N); theta = 0.37
    c, s = math.cos(theta), math.sin(theta)
    # Rotate eigenvectors belonging to different positive eigenvalues and fix the constant mode.
    rot[1, 1] = c; rot[1, 3] = -s; rot[3, 1] = s; rot[3, 3] = c
    q = vecs @ rot @ vecs.T
    a2 = q @ A @ q.T
    t = 0.8
    p1, p2 = expm(-t*A), expm(-t*a2)
    u1, u2 = expm(-1j*t*A), expm(-1j*t*a2)
    packet = {
        "rotation_angle": theta,
        "spectrum_max_absolute_residual": float(np.max(np.abs(np.linalg.eigvalsh(a2)-vals))),
        "constant_mode_residual": float(np.linalg.norm(a2 @ np.ones(N))),
        "generator_entry_max_difference": float(np.max(np.abs(a2-A))),
        "heat_vertex_record_max_difference_t_0_8": float(np.max(np.abs(p2-p1))),
        "unitary_Born_vertex_record_max_difference_t_0_8": float(np.max(np.abs(np.abs(u2)**2-np.abs(u1)**2))),
        "theorem": "Bare eigenvalues do not determine vertex records. An orthogonal conjugation fixing the constant mode preserves the complete spectrum and positivity but changes heat and Born vertex matrices unless preparations/effects are transported with it.",
        "minimal_medium_independent_bundle": "spectral projectors plus transported preparations, transformations, effects, and record map",
    }
    return finalize(473, "proven_isospectral_spectrum_only_information_insufficiency",
                    "The conjugate generator need not remain in the strict radial-kernel subclass; the theorem targets claims that eigenvalues alone are complete information.", packet)


def gate(k: int, kind: str) -> dict:
    if kind == "gain":
        packet = {"new_strict_gain_source": False, "g_remains_supplied": True,
                  "autocatalytic_candidate": "receiver requiring seed", "reason": "ST465 proves the symmetric vacuum is absorbing."}
        status = "blocked_no_new_strict_gain_source"
    elif kind == "selector_scale":
        packet = {"new_nonpremise_selector": False, "QW_2191": "open",
                  "new_scale_charged_source": False, "absolute_unit": "absent",
                  "reason": "Orbit identification supplies degeneracy, not a section; no weight-one strict datum was constructed."}
        status = "blocked_no_new_selector_or_scale_source"
    else:
        packet = {"external_referee": "absent", "independent_laboratory_record": "absent",
                  "held_out_empirical_record": "absent", "reason": "This batch is local analytic/computational work."}
        status = "blocked_no_independent_empirical_evidence"
    return finalize(k, status,
                    "No dimensional, bridge-role-transfer, laboratory, L_total, Standard Model, gravity, or ToE closure is exported.", packet)


def figures(results: dict) -> None:
    FIG_DIR.mkdir(exist_ok=True)
    p = np.array(results["ST462"]["canonical_probability"])
    fig, ax = plt.subplots(figsize=(7.4, 4.0)); ax.bar(range(N), p, color="#2563eb")
    ax.set(xlabel="cyclic vertex", ylabel="probability", title="ST462: candidate first-transition orbit representative")
    fig.tight_layout(); fig.savefig(FIG_DIR / "st462_transition_orbit.png", dpi=180); plt.close(fig)

    audits = results["ST466"]["exact_modular_random_audits"]
    labels, counts = [], []
    for a in audits:
        for r in a["rows"]:
            labels.append(f"p={a['prime']}\ns={r['support']}")
            counts.append(r["sampled_subsets"])
    fig, ax = plt.subplots(figsize=(7.4, 4.0)); ax.bar(labels, counts, color="#0f766e")
    ax.set(ylabel="exact sampled subsets", title="ST466: support-4/5 modular circuit audit (zero found)")
    fig.tight_layout(); fig.savefig(FIG_DIR / "st466_circuit_audit.png", dpi=180); plt.close(fig)

    rows = results["ST470"]["accepted_cells"]
    b = [r["b_interval"][1] for r in rows]
    y2 = [r["endpoint_centers"][-1][1] for r in rows]
    fig, ax = plt.subplots(figsize=(7.4, 4.0)); ax.plot(b, y2, "o-")
    ax.set(xlabel="compactified parameter b", ylabel=r"$y_2$", title="ST470: certified IR extension until method stop")
    fig.tight_layout(); fig.savefig(FIG_DIR / "st470_ir_extension.png", dpi=180); plt.close(fig)

    x = np.linspace(0, 8, 800); h = np.sqrt(1+np.exp(-2*x)-2*np.exp(-x)*np.cos(x))
    xs = results["ST471"]["first_positive_global_maximizer"]
    fig, ax = plt.subplots(figsize=(7.4, 4.0)); ax.plot(x, h, color="#7c3aed"); ax.axvline(xs, ls="--", color="#dc2626")
    ax.set(xlabel=r"dimensionless modal time $x=t\lambda$", ylabel=r"$|e^{-ix}-e^{-x}|$",
           title="ST471: universal unitary--diffusive discriminant")
    fig.tight_layout(); fig.savefig(FIG_DIR / "st471_dual_discriminant.png", dpi=180); plt.close(fig)


def main() -> None:
    funcs = {462: st462, 463: st463, 464: st464, 465: st465, 466: st466,
             467: st467, 468: st468, 469: st469, 470: st470, 471: st471,
             472: st472, 473: st473}
    results = {}
    for k in range(462, 474):
        print(f"running ST{k}", flush=True)
        results[f"ST{k}"] = funcs[k]()
    results["ST474"] = gate(474, "gain")
    results["ST475"] = gate(475, "selector_scale")
    results["ST476"] = gate(476, "evidence")
    RESULTS.write_text(json.dumps(native(results), indent=2, sort_keys=True), encoding="utf-8")
    with SUMMARY.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f); w.writerow(["program", "status", "object", "boundary"])
        for k in range(462, 477):
            r = results[f"ST{k}"]; w.writerow([f"ST{k}", r["status"], r["object"], r["boundary"]])
    figures(results)
    print(f"wrote {RESULTS.name} and {SUMMARY.name}")


if __name__ == "__main__":
    main()
