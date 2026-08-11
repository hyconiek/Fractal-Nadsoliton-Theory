#!/usr/bin/env python3
"""FIN ST82--ST93: projection-hypothesis falsification batch.

The hypothesis H_PROJ is deliberately counterfactual: observable physics is an
operational projection/quotient of a deeper FIN object.  The batch tests what
the present finite mathematics proves, what it merely permits, and which new
resources remain independent.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import platform
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scipy
from scipy.linalg import expm
from scipy.optimize import differential_evolution, minimize_scalar, root

from fin_st01_st15_research import N, random_orthogonal_fixing_uniform, strict_operator
from fin_st28_st45_research import dyadic_lift, saturation_energy_gradient_hessians
from fin_st46_st57_research import carrier_probability_table
from fin_st70_st81_research import reflection_expansion, st77_pseudo_arclength_fold


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "FIN_ST82_ST93_Results.json"
SUMMARY = ROOT / "FIN_ST82_ST93_Summary.csv"
SCHEMA = ROOT / "FIN_ST93_Projection_Dependency_Schema.json"
FIG_DIR = ROOT / "FIN_ST82_ST93_Figures"
SEED = 20260818


def native(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): native(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [native(item) for item in value]
    if isinstance(value, np.ndarray):
        return native(value.tolist())
    if isinstance(value, (np.floating, np.integer)):
        return value.item()
    if isinstance(value, complex):
        return {"real": value.real, "imag": value.imag}
    return value


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def st82_independent_certificate() -> dict:
    certificate_path = ROOT / "FIN_ST82_Independent_Rational_Certificate.json"
    if certificate_path.exists():
        certificate = json.loads(certificate_path.read_text(encoding="utf-8"))
        status = "proven_source_code_exact_rational_interval_regeneration"
        boundary = certificate["proof_boundary"]
    else:
        previous = json.loads((ROOT / "FIN_ST70_Rational_Replay_Certificate.json").read_text(encoding="utf-8"))
        certificate = {
            "fallback_replay": True,
            "all_checks_pass": previous["all_checks_pass"],
            "collision_speed_interval": previous["collision_speed_interval"],
        }
        status = "blocked_independent_regeneration_fallback_replay_only"
        boundary = "The independent exact-rational regeneration did not finish; only the accepted ST70 replay is carried over."
    return {
        "program": "ST82",
        "object": "Independent Rational-Interval Regeneration of the ST58 Memory Collision",
        "certificate": certificate,
        "certificate_sha256": digest(certificate_path) if certificate_path.exists() else None,
        "status": status,
        "boundary": boundary,
    }


def st83_fold_certification_attempt(a: np.ndarray) -> dict:
    previous = json.loads((ROOT / "FIN_ST70_ST81_Results.json").read_text(encoding="utf-8"))["ST77"]
    fold = np.r_[previous["fold_state_reduced"], previous["fold_kappa"], previous["fold_null_vector"]]
    expansion = reflection_expansion()
    selected = np.arange(7)

    def system(candidate: np.ndarray) -> np.ndarray:
        q, kappa, vector = candidate[:7], candidate[7], candidate[8:]
        state = expansion @ q
        _, gradient, hessian, _ = saturation_energy_gradient_hessians(kappa * a, state)
        matrix = hessian[np.ix_(selected, np.arange(N))] @ expansion
        return np.r_[gradient[selected], matrix @ vector, 0.5 * (vector @ vector - 1.0)]

    refined = root(system, fold, method="lm", options={"ftol": 1e-14, "xtol": 1e-14, "gtol": 1e-14, "maxiter": 20000})
    x = refined.x
    residual = float(np.linalg.norm(system(x), ord=np.inf))
    scales = [1e-4, 3e-5, 1e-5, 3e-6, 1e-6, 3e-7]
    minimum_singular_values = []
    for h in scales:
        jacobian = np.column_stack([(system(x + h * np.eye(len(x))[j]) - system(x - h * np.eye(len(x))[j])) / (2 * h) for j in range(len(x))])
        minimum_singular_values.append(float(np.linalg.svd(jacobian, compute_uv=False)[-1]))
    return {
        "program": "ST83",
        "object": "Simple-Fold Certification Attempt",
        "refined_fold_kappa": float(x[7]),
        "fold_system_residual_inf": residual,
        "finite_difference_scales": scales,
        "augmented_J_minimum_singular_values": minimum_singular_values,
        "minimum_over_scales": min(minimum_singular_values),
        "interval_certificate_obtained": False,
        "failure_reason": "No outward interval enclosure of the complete 15-variable Jacobian and its inverse was obtained in this local batch.",
        "status": "strong_numerical_simple_fold_evidence_but_interval_goal_not_met",
        "boundary": "Small residual and stable nonzero floating singular values do not constitute an interval theorem or global branch classification.",
    }


def st84_passive_memory_obstruction() -> dict:
    previous = json.loads((ROOT / "FIN_ST58_Full_Interval_Certificate.json").read_text(encoding="utf-8"))
    poles = previous["pole_intervals"]
    residues = previous["residue_intervals"]
    positive = all(float(row[0]) > 0 for row in poles + residues)
    samples = []
    for time in [0.0, 0.1, 0.5, 1.0, 3.0]:
        bounds = []
        for derivative in range(5):
            lo = sum(float(r[0]) * float(p[0]) ** derivative * math.exp(-float(p[1]) * time) for r, p in zip(residues, poles))
            bounds.append(lo)
        samples.append({"time": time, "lower_bounds_for_signed_derivatives_0_to_4": bounds})
    return {
        "program": "ST84",
        "object": "Passive Stieltjes-Memory Temporal-Gain Obstruction",
        "all_certified_poles_and_residues_positive": positive,
        "sampled_complete_monotonicity_lower_bounds": samples,
        "theorem": (
            "For M(z)=sum_j r_j/(z+p_j), p_j,r_j>0, the impulse response m(t)=sum_j r_j exp(-p_j t) "
            "is completely monotone: (-1)^n m^(n)(t)=sum_j r_j p_j^n exp(-p_j t)>0. It is a passive positive-real "
            "memory and cannot by itself create an active positive-gain temporal selector. Active gain requires a pump, "
            "negative-residue/nonequilibrium component, boundary condition, or another time-oriented premise."
        ),
        "status": "proven_passive_memory_no_go_for_intrinsic_active_gain",
        "boundary": "The theorem does not prohibit passive retardation or a time arrow supplied by an external preparation/environment.",
    }


def st85_spin_lift_classification(a: np.ndarray) -> dict:
    theta_periodic = 2 * math.pi * np.arange(N) / N
    theta_antiperiodic = 2 * math.pi * (np.arange(N) + 0.5) / N
    scalar_periodic = np.sort(np.linalg.eigvalsh(a))
    # A half-angle boundary condition is a representation choice, not a new scalar spectrum.
    reflection = (-np.arange(N)) % N
    projector_plus = np.zeros((2 * N, N), dtype=complex)
    projector_minus = np.zeros((2 * N, N), dtype=complex)
    for x in range(N):
        projector_plus[x, x] = 1 / math.sqrt(2)
        projector_plus[x + N, x] = 1 / math.sqrt(2)
        projector_minus[x, x] = 1 / math.sqrt(2)
        projector_minus[x + N, x] = -1 / math.sqrt(2)
    return {
        "program": "ST85",
        "object": "State/Boundary Classification of C24 Spin Lifts",
        "spin_structures_on_cycle": ["periodic holonomy +1", "antiperiodic holonomy -1"],
        "periodic_momenta": theta_periodic,
        "antiperiodic_momenta": theta_antiperiodic,
        "strict_scalar_eigenvalues": scalar_periodic,
        "plus_minus_isometry_residuals": [
            float(np.linalg.norm(projector_plus.conj().T @ projector_plus - np.eye(N))),
            float(np.linalg.norm(projector_minus.conj().T @ projector_minus - np.eye(N))),
            float(np.linalg.norm(projector_plus.conj().T @ projector_minus)),
        ],
        "theorem": (
            "A cycle admits two Z2 spin boundary classes, and the antiperiodic class realizes half-integer momenta. "
            "Choosing either class is extra boundary/state data: the real scalar strict operator and its dihedral stabilizer "
            "do not distinguish the two lifts or select an oriented section."
        ),
        "status": "proven_conditional_spin_lift_classification_and_strict_nonselection",
        "boundary": "Existence of the C24 lift is not a strict source theorem for physical spin, chirality, or fermions.",
    }


def st86_projection_fiber_theorem(a: np.ndarray) -> dict:
    plus = np.zeros((2 * N, N))
    minus = np.zeros((2 * N, N))
    for x in range(N):
        plus[x, x] = plus[x + N, x] = 1 / math.sqrt(2)
        minus[x, x] = 1 / math.sqrt(2)
        minus[x + N, x] = -1 / math.sqrt(2)
    rows = []
    fine_blocks = []
    for q in [0.0, 0.07, 0.2, 0.7, 1.5]:
        lifted = dyadic_lift(a, q)
        coarse = plus.T @ lifted @ plus
        mixed = plus.T @ lifted @ minus
        fine = minus.T @ lifted @ minus
        fine_blocks.append(fine)
        heat_residual = float(np.linalg.norm(plus.T @ expm(-0.37 * lifted) @ plus - expm(-0.37 * a)))
        unitary_residual = float(np.linalg.norm(plus.T @ expm(-0.23j * lifted) @ plus - expm(-0.23j * a)))
        boltzmann = expm(-0.8 * lifted)
        coarse_boltzmann = expm(-0.8 * a)
        global_gibbs = boltzmann / np.trace(boltzmann)
        compressed_gibbs = plus.T @ global_gibbs @ plus
        coarse_probability = float(np.trace(compressed_gibbs))
        conditioned_gibbs = compressed_gibbs / coarse_probability
        coarse_gibbs = coarse_boltzmann / np.trace(coarse_boltzmann)
        rows.append({
            "q": q,
            "coarse_generator_residual": float(np.linalg.norm(coarse - a)),
            "coarse_fine_mixing_residual": float(np.linalg.norm(mixed)),
            "fine_trace": float(np.trace(fine)),
            "heat_compression_residual": heat_residual,
            "unitary_compression_residual": unitary_residual,
            "unnormalized_boltzmann_compression_residual": float(np.linalg.norm(plus.T @ boltzmann @ plus - coarse_boltzmann)),
            "globally_normalized_gibbs_coarse_probability": coarse_probability,
            "coarse_conditioned_gibbs_residual": float(np.linalg.norm(conditioned_gibbs - coarse_gibbs)),
        })
    differences = [float(np.linalg.norm(fine_blocks[i] - fine_blocks[j])) for i in range(len(fine_blocks)) for j in range(i)]
    return {
        "program": "ST86",
        "object": "Exact Coarse-Projection Fiber Theorem for Dyadic FIN Lifts",
        "rows": rows,
        "maximum_coarse_generator_residual": max(row["coarse_generator_residual"] for row in rows),
        "maximum_coarse_fine_mixing_residual": max(row["coarse_fine_mixing_residual"] for row in rows),
        "maximum_sampled_fine_block_difference": max(differences),
        "theorem": (
            "Let P=(I,I)^T/sqrt(2) and F=(I,-I)^T/sqrt(2). For every dyadic lift A24(q), the intertwining identity "
            "A24(q)P=P A12 and symmetry give P* A24(q)P=A12 and P* A24(q)F=0. Hence P* f(A24(q))P=f(A12) "
            "for every spectral function f. The fine block F* A24(q)F varies with q, so infinitely many inequivalent "
            "deeper generators have exactly the same complete coarse spectral dynamics."
        ),
        "projection_hypothesis_verdict": (
            "H_PROJ is mathematically compatible with FIN and has an explicit exact model, but is nonidentifiable: "
            "coarse observations cannot infer which q, or even uniqueness of the deeper lift."
        ),
        "normalization_caveat": (
            "Unnormalized Boltzmann operators compress exactly. A globally normalized 24-state Gibbs state assigns a q-dependent total probability to the coarse sector; only the Gibbs state conditioned and renormalized inside that sector equals the 12-state Gibbs state."
        ),
        "status": "proven_exact_projection_fiber_and_nonidentifiability_theorem",
        "boundary": "This proves a finite operator quotient, not that the physical universe is its image or that q is physically real.",
    }


def st87_robust_likelihood(a: np.ndarray) -> dict:
    rng = np.random.default_rng(SEED + 87)
    q = random_orthogonal_fixing_uniform(rng)
    p = carrier_probability_table(a, np.eye(N), transported=False) / N
    qdist = carrier_probability_table(a, q, transported=False) / N
    uniform = np.full_like(p, 1.0 / p.size)
    epsilon_max = 0.12

    def chernoff(eps: np.ndarray) -> tuple[float, float]:
        pp = (1 - eps[0]) * p + eps[0] * uniform
        qq = (1 - eps[1]) * qdist + eps[1] * uniform
        objective = lambda s: float(np.sum(pp**s * qq ** (1 - s)))
        opt = minimize_scalar(objective, bounds=(0.0, 1.0), method="bounded", options={"xatol": 1e-13})
        return -math.log(max(opt.fun, 1e-300)), float(opt.x)

    adversary = differential_evolution(lambda eps: chernoff(eps)[0], [(0, epsilon_max), (0, epsilon_max)], seed=SEED + 87, tol=1e-10, polish=True)
    info, sstar = chernoff(adversary.x)
    nominal, nominal_s = chernoff(np.zeros(2))
    target = 0.01
    required_events = math.ceil(math.log(1 / target) / info)
    return {
        "program": "ST87",
        "object": "Robust Chernoff-Likelihood Nuisance Optimization",
        "contamination_model": "independent uniform-mixture fractions epsilon_P,epsilon_Q in [0,0.12]",
        "nominal_chernoff_information_per_event": nominal,
        "nominal_s": nominal_s,
        "worst_case_epsilons": adversary.x,
        "worst_case_chernoff_information_per_event": info,
        "worst_case_s": sstar,
        "asymptotic_events_for_exp_bound_0_01": required_events,
        "optimization_success": bool(adversary.success),
        "status": "strong_numerical_global_nuisance_optimization_for_declared_synthetic_model",
        "boundary": "The carriers, contamination law, calibration and observations are synthetic; the asymptotic Chernoff bound is not a finite-laboratory guarantee.",
    }


def st88_reversible_faces(a: np.ndarray) -> dict:
    eigenvalues = np.linalg.eigvalsh(a)
    groups = []
    for value in eigenvalues:
        if not groups or abs(value - groups[-1][0]) > 1e-9:
            groups.append([float(value), 1])
        else:
            groups[-1][1] += 1
    dimensions = [row[1] for row in groups]
    quotient_dimension = sum(d * d for d in dimensions)
    kernel_dimension = N * N - quotient_dimension
    return {
        "program": "ST88",
        "object": "Reversible Faces of the ST75 Intertwiner Cone",
        "spectral_block_dimensions": dimensions,
        "quotient_algebra": "C direct-sum M2^5 direct-sum C",
        "quotient_linear_dimension": quotient_dimension,
        "pinching_kernel_dimension": kernel_dimension,
        "reversible_quotient_group": "S2 x (PU(2)^5 semidirect S5)",
        "theorem": (
            "Every full-space ST75 intertwiner factors through spectral pinching, whose kernel has positive dimension, "
            "so none has a CPTP inverse on M12. On the pinched finite C*-algebra, reversible channels are exactly "
            "*-automorphisms: unitary conjugations inside blocks and permutations only among isomorphic blocks."
        ),
        "status": "proven_full_space_irreversibility_and_quotient_reversible_face_classification",
        "boundary": "This operational quotient is not evidence that physical measurement uses this pinching channel.",
    }


def st89_gibbs_preserving_separation(a: np.ndarray) -> dict:
    beta = 0.8
    energies = np.array([0.0, float(np.linalg.eigvalsh(a)[1])])
    gamma = np.exp(-beta * energies); gamma /= gamma.sum()
    delta = 0.15 * math.sqrt(gamma[0] * gamma[1])
    sigma0 = np.array([[gamma[0], delta], [delta, gamma[1]]], dtype=complex)
    sigma1 = np.array([[gamma[0], -gamma[0] * delta / gamma[1]], [-gamma[0] * delta / gamma[1], gamma[1]]], dtype=complex)
    fixed = gamma[0] * sigma0 + gamma[1] * sigma1
    covariance_witness = abs(sigma0[0, 1])
    return {
        "program": "ST89",
        "object": "Explicit Gibbs-Preserving versus Thermal-Operation Separation",
        "two_level_gamma": gamma,
        "prepared_sigma0_eigenvalues": np.linalg.eigvalsh(sigma0),
        "prepared_sigma1_eigenvalues": np.linalg.eigvalsh(sigma1),
        "gibbs_fixed_residual": float(np.linalg.norm(fixed - np.diag(gamma))),
        "time_covariance_violation_witness": covariance_witness,
        "theorem": (
            "The measure-in-energy/prepare channel Phi(rho)=rho00 sigma0+rho11 sigma1 is CPTP and fixes gamma by construction. "
            "It maps the stationary input |0><0| to sigma0 with nonzero energy coherence, so it is not time-translation "
            "covariant. Energy-conserving thermal operations are covariant; therefore Gibbs preservation is strictly weaker."
        ),
        "status": "proven_explicit_GP_not_thermal_operation_counterexample",
        "boundary": "The two-level scale and bath interpretation are added conversion resources, not derived FIN physics.",
    }


def st90_disconnected_branch_search(a: np.ndarray) -> dict:
    expansion = reflection_expansion()
    selected = np.arange(7)
    rng = np.random.default_rng(SEED + 90)
    rows = []
    for kappa in [0.02, 0.05, 0.1, 0.5, 1.0]:
        solutions = []
        for _ in range(80):
            guess = rng.normal(scale=2.0, size=7)
            def equation(q: np.ndarray) -> np.ndarray:
                state = expansion @ q
                return saturation_energy_gradient_hessians(kappa * a, state)[1][selected]
            solved = root(equation, guess, method="lm")
            if np.linalg.norm(equation(solved.x), ord=np.inf) > 1e-8:
                continue
            state = expansion @ solved.x
            signature = np.r_[np.sort(np.abs(state)), np.linalg.norm(state)]
            if not any(np.linalg.norm(signature - old[0]) < 1e-5 for old in solutions):
                power = float(state @ state)
                solutions.append((signature, float(np.sum(state**4) / power**2) if power > 1e-14 else 0.0))
        rows.append({"kappa": kappa, "distinct_reflection_reduced_solutions_found": len(solutions), "IPRs": [item[1] for item in solutions]})
    return {
        "program": "ST90",
        "object": "Reflection-Reduced Disconnected-Branch Multistart Search",
        "starts_per_kappa": 80,
        "rows": rows,
        "status": "moderate_numerical_branch_inventory_not_exhaustive",
        "boundary": "A finite Levenberg-Marquardt multistart can miss branches and does not prove existence, uniqueness, or absence at any kappa.",
    }


def st91_dynamic_backreaction(a: np.ndarray) -> dict:
    lambda1 = float(np.linalg.eigvalsh(a)[1])
    mu = 0.35
    gc = mu / lambda1
    rows = []
    for g in [0.0, 0.5 * gc, 0.999 * gc, 1.001 * gc, 2 * gc]:
        origin_rate = mu - g * lambda1
        rows.append({"coupling": g, "origin_linear_growth_rate": origin_rate, "origin_stable": origin_rate < 0, "twelve_branch_side": g < gc})
    return {
        "program": "ST91",
        "object": "Dynamic Exchange of Stability in the Constructed Twelve-Branch Backreaction",
        "critical_coupling": gc,
        "rows": rows,
        "theorem": (
            "For gradient dynamics zdot=-dV_g/dzbar, the origin has radial growth rate mu-g lambda1. "
            "It is unstable below g_c=mu/lambda1 and stable above it. The twelve nonzero minima exist only below g_c, "
            "so the constructed model exhibits an exchange of stability at the same threshold as ST78."
        ),
        "status": "proven_conditional_local_dynamic_bifurcation_for_constructed_model",
        "boundary": "The potential, coupling, metric and gradient-flow clock are added; this is not strict FIN time evolution.",
    }


def st92_reflection_odd_search(a: np.ndarray) -> dict:
    reflection = np.zeros((N, N))
    for x in range(N):
        reflection[(-x) % N, x] = 1
    candidates = [np.eye(N)]
    current = np.eye(N)
    for _ in range(1, 13):
        current = current @ a
        candidates.append(current.copy())
    diagonal_candidates = [np.diag(np.diag(item)) for item in candidates]
    candidates.extend(diagonal_candidates)
    commutators = [left @ right - right @ left for left in candidates[:13] for right in candidates[:13]]
    all_candidates = candidates + commutators
    relative_odd_norms = [
        float(np.linalg.norm(0.5 * (item - reflection @ item @ reflection)) / max(np.linalg.norm(item), 1.0))
        for item in candidates
    ]
    commutator_relative_residuals = [
        float(np.linalg.norm(left @ right - right @ left) / max(np.linalg.norm(left) * np.linalg.norm(right), 1.0))
        for left in candidates[:13] for right in candidates[:13]
    ]
    exact_rank = len({round(float(value), 12) for value in np.linalg.eigvalsh(a)})
    return {
        "program": "ST92",
        "object": "Finite Strict-Generated Reflection-Odd Object Search",
        "declared_class": "powers A^0..A^12, their diagonal contractions, and all pairwise commutators of powers",
        "candidate_count": len(all_candidates),
        "generated_span_dimension_exact_from_distinct_spectral_values": exact_rank,
        "maximum_relative_reflection_odd_component_norm": max(relative_odd_norms),
        "maximum_relative_floating_commutator_residual": max(commutator_relative_residuals),
        "theorem": (
            "Because RAR=A, every polynomial in A is reflection even; its diagonal is constant for circulant A; "
            "and all polynomial commutators vanish. The entire declared finite class therefore contains no nonzero "
            "reflection-odd source and cannot discharge QW-2191."
        ),
        "status": "proven_no_go_in_declared_strict_generated_class",
        "boundary": "State-dependent nonlinear observables, boundaries, records or new typed tensors lie outside the exhausted class.",
    }


def st93_projection_schema() -> dict:
    nodes = {
        "W0": {"kind": "strict mathematical core", "dependencies": []},
        "PF": {"kind": "projection fiber", "dependencies": ["W0"], "evidence": "ST86"},
        "IR": {"kind": "irreversible quotient", "dependencies": ["W0"], "evidence": "ST75/ST88"},
        "SP": {"kind": "spin/projective lift", "dependencies": ["W0", "BC"], "evidence": "ST85"},
        "BC": {"kind": "boundary/state choice", "dependencies": []},
        "SEL": {"kind": "selector/symmetry breaking", "dependencies": ["SRC_SEL"]},
        "SRC_SEL": {"kind": "missing strict source", "dependencies": []},
        "TIME": {"kind": "dimensionful clock/order", "dependencies": ["CLK", "W0"]},
        "CLK": {"kind": "missing clock calibration", "dependencies": []},
        "SCALE": {"kind": "dimensionful energy/length scale", "dependencies": ["UNIT", "W0"]},
        "UNIT": {"kind": "missing physical unit source", "dependencies": []},
        "OBS": {"kind": "operational observables", "dependencies": ["PREP", "INST", "REC"]},
        "PREP": {"kind": "missing physical preparations", "dependencies": []},
        "INST": {"kind": "missing calibrated instrument/apparatus", "dependencies": []},
        "REC": {"kind": "missing independent physical record", "dependencies": []},
        "H_PROJ": {"kind": "counterfactual interpretation", "dependencies": ["PF", "IR", "SP", "SEL", "TIME", "SCALE", "OBS"]},
    }
    missing = [key for key, value in nodes.items() if value["kind"].startswith("missing")]
    # deterministic acyclicity check
    visiting, visited = set(), set()
    def visit(key: str) -> None:
        if key in visiting:
            raise RuntimeError("cycle")
        if key in visited:
            return
        visiting.add(key)
        for dep in nodes[key]["dependencies"]:
            visit(dep)
        visiting.remove(key); visited.add(key)
    for key in nodes:
        visit(key)
    schema = {
        "schema_version": "1.0",
        "hypothesis": "H_PROJ: observable universe is an operational projection/quotient of a deeper FIN object",
        "hypothesis_status": "counterfactual_not_established_physics",
        "nodes": nodes,
        "missing_source_nodes": missing,
        "missing_source_count": len(missing),
        "acyclic": len(visited) == len(nodes),
        "central_verdict": (
            "ST86 and ST88 prove exact finite projection/quotient mechanisms. They establish mathematical possibility "
            "and nonidentifiability, not a derivation of the physical projection map. H_PROJ still depends on independent "
            "selector, boundary, clock, unit, preparation, instrument and record resources."
        ),
    }
    SCHEMA.write_text(json.dumps(schema, indent=2, sort_keys=True), encoding="utf-8")
    return {
        "program": "ST93",
        "object": "Machine-Readable FIN Projection-Hypothesis Dependency Graph",
        "schema_file": SCHEMA.name,
        "schema_sha256": digest(SCHEMA),
        **schema,
        "status": "proven_dependency_DAG_and_open_projection_source_inventory",
        "boundary": "A consistent dependency graph is an audit artifact, not evidence for the counterfactual physical hypothesis.",
    }


def make_figures(results: dict) -> None:
    FIG_DIR.mkdir(exist_ok=True)
    fig, ax = plt.subplots(figsize=(7, 4)); st84 = results["ST84"]
    for n in range(5): ax.plot([r["time"] for r in st84["sampled_complete_monotonicity_lower_bounds"]], [r["lower_bounds_for_signed_derivatives_0_to_4"][n] for r in st84["sampled_complete_monotonicity_lower_bounds"]], "o-", label=f"n={n}")
    ax.set(yscale="log", xlabel="dimensionless time", ylabel="certified positive lower expression", title="ST84 passive complete monotonicity"); ax.legend(); fig.tight_layout(); fig.savefig(FIG_DIR / "st84_passive_memory.png", dpi=190); plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(10, 4)); st86 = results["ST86"]
    axes[0].plot([r["q"] for r in st86["rows"]], [r["coarse_generator_residual"] for r in st86["rows"]], "o-"); axes[0].set(yscale="log", xlabel="fine parameter q", ylabel="coarse residual", title="coarse generator is invariant")
    axes[1].plot([r["q"] for r in st86["rows"]], [r["fine_trace"] for r in st86["rows"]], "o-"); axes[1].set(xlabel="fine parameter q", ylabel="fine-block trace", title="hidden fiber varies")
    fig.tight_layout(); fig.savefig(FIG_DIR / "st86_projection_fiber.png", dpi=190); plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 4)); st83 = results["ST83"]
    ax.loglog(st83["finite_difference_scales"], st83["augmented_J_minimum_singular_values"], "o-"); ax.set(xlabel="finite-difference step", ylabel="minimum singular value", title="ST83 numerical transversality (not interval proof)"); fig.tight_layout(); fig.savefig(FIG_DIR / "st83_fold_attempt.png", dpi=190); plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 4)); st90 = results["ST90"]
    ax.bar([str(r["kappa"]) for r in st90["rows"]], [r["distinct_reflection_reduced_solutions_found"] for r in st90["rows"]]); ax.set(xlabel="kappa", ylabel="distinct signatures found", title="ST90 finite multistart inventory"); fig.tight_layout(); fig.savefig(FIG_DIR / "st90_branch_inventory.png", dpi=190); plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 4)); st91 = results["ST91"]
    ax.axhline(0, color="black", lw=0.8); ax.plot([r["coupling"] for r in st91["rows"]], [r["origin_linear_growth_rate"] for r in st91["rows"]], "o-"); ax.axvline(st91["critical_coupling"], color="red", ls="--"); ax.set(xlabel="coupling g", ylabel="origin growth rate", title="ST91 conditional exchange of stability"); fig.tight_layout(); fig.savefig(FIG_DIR / "st91_dynamic_threshold.png", dpi=190); plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 4)); labels = ["full M12", "pinched quotient"]; values = [0, 1]; ax.bar(labels, values); ax.set_ylim(0, 1.2); ax.set(ylabel="reversible channels exist", title="ST88 reversibility only after quotienting"); fig.tight_layout(); fig.savefig(FIG_DIR / "st88_reversible_faces.png", dpi=190); plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 4)); st87 = results["ST87"]; ax.bar(["nominal", "worst nuisance"], [st87["nominal_chernoff_information_per_event"], st87["worst_case_chernoff_information_per_event"]]); ax.set(ylabel="Chernoff information/event", title="ST87 nuisance weakens synthetic discrimination"); fig.tight_layout(); fig.savefig(FIG_DIR / "st87_robust_likelihood.png", dpi=190); plt.close(fig)

    fig, ax = plt.subplots(figsize=(8, 4)); missing = results["ST93"]["missing_source_nodes"]; ax.bar(missing, np.ones(len(missing))); ax.tick_params(axis="x", rotation=30); ax.set(ylabel="open (1=yes)", title="H_PROJ still requires independent source nodes"); fig.tight_layout(); fig.savefig(FIG_DIR / "st93_missing_sources.png", dpi=190); plt.close(fig)


def write_summary(results: dict) -> None:
    with SUMMARY.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle); writer.writerow(["program", "object", "status"])
        for number in range(82, 94):
            row = results[f"ST{number}"]; writer.writerow([row["program"], row["object"], row["status"]])


def main() -> None:
    _, a, _ = strict_operator()
    results: dict[str, Any] = {"metadata": {"programs": "ST82-ST93", "date": "2026-08-11", "seed": SEED, "python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__, "hypothesis": "H_PROJ counterfactual"}}
    results["ST82"] = st82_independent_certificate()
    results["ST83"] = st83_fold_certification_attempt(a)
    results["ST84"] = st84_passive_memory_obstruction()
    results["ST85"] = st85_spin_lift_classification(a)
    results["ST86"] = st86_projection_fiber_theorem(a)
    results["ST87"] = st87_robust_likelihood(a)
    results["ST88"] = st88_reversible_faces(a)
    results["ST89"] = st89_gibbs_preserving_separation(a)
    results["ST90"] = st90_disconnected_branch_search(a)
    results["ST91"] = st91_dynamic_backreaction(a)
    results["ST92"] = st92_reflection_odd_search(a)
    results["ST93"] = st93_projection_schema()
    results["recommended_next_programs"] = [
        {"id": "ST94", "priority": 1, "study": "prove or refute uniqueness of a projection/conditional expectation from a new operational sufficiency axiom"},
        {"id": "ST95", "priority": 2, "study": "classify all coarse-preserving A24 lifts, not only the dyadic q-family"},
        {"id": "ST96", "priority": 3, "study": "derive a Petz-recovery criterion and quantify exactly which projected coherences are recoverable"},
        {"id": "ST97", "priority": 4, "study": "construct a state-dependent selector candidate and test equivariant bifurcation plus QW-2191"},
        {"id": "ST98", "priority": 5, "study": "seek a nonequilibrium signed-residue memory sourced internally, or strengthen the ST84 no-go"},
        {"id": "ST99", "priority": 6, "study": "build a certified interval Jacobian/preconditioner for the ST83 fold"},
        {"id": "ST100", "priority": 7, "study": "classify projective spin lifts compatible with refinement 12 to 24 to 48"},
        {"id": "ST101", "priority": 8, "study": "derive identifiability bounds for hidden fine parameters from deliberately non-coarse observables"},
        {"id": "ST102", "priority": 9, "study": "formalize projection composition as a finite C*-algebra tower and test Markov commuting squares"},
        {"id": "ST103", "priority": 10, "study": "replace the ST87 differential evolution by interval global optimization of the nuisance box"},
        {"id": "ST104", "priority": 11, "study": "separate thermal operations, enhanced thermal operations, and Gibbs-preserving maps for strict spectral degeneracies"},
        {"id": "ST105", "priority": 12, "study": "build a projection-specific falsification protocol with synthetic alternatives sharing every coarse channel"},
    ]
    results["central_verdict"] = (
        "The present FIN mathematics contains exact projection and irreversible quotient mechanisms, most sharply ST86. "
        "It therefore supports H_PROJ as a coherent mathematical interpretation. The same theorem proves fatal nonuniqueness: "
        "infinitely many deeper lifts generate identical coarse dynamics. No current strict theorem selects the projection, hidden lift, "
        "spin sector, time arrow, physical scale, apparatus or record. H_PROJ is not derived physics and has no empirical confirmation."
    )
    results["epistemic_boundary"] = "No QW-2191 closure, legacy-to-strict role transfer, laboratory evidence, dimensional source, Standard Model, gravity, L_total or ToE closure is claimed."
    make_figures(results); write_summary(results)
    RESULTS.write_text(json.dumps(native(results), indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({"results": RESULTS.name, "programs": 12, "figures": 8}, indent=2))


if __name__ == "__main__":
    main()
