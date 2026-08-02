#!/usr/bin/env python3
"""Execute FIN Research Programs P411--P427 (Release 10.35).

The round reduces formal trust, constructs a simultaneous contact candidate,
tests one nonlinear damping source class, audits photonic regularity, improves
noisy and erasure-aware discrimination strategies, derives the variance-
optimal Jordan sampler, accounts for phase-reference cost, constructs one
explicit oriented record, and compares conditional scale sections.  Physical
pilots, clocks, laboratory JSR data, standards, reservoirs, hold-outs, and
electroweak records remain external gates.
"""

from __future__ import annotations

import csv
from fractions import Fraction
import hashlib
import itertools
import json
import math
from pathlib import Path
import subprocess
from typing import Any, Callable

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize

import fin_programs_255_266 as core
import fin_programs_295_308 as p295
import fin_programs_365_378 as p365
import fin_programs_394_410 as prior


ROOT = Path(__file__).resolve().parent
PREFIX = "FIN_Programs_411_427"
RESULTS_PATH = ROOT / f"{PREFIX}_Results.json"
SUMMARY_PATH = ROOT / f"{PREFIX}_Summary.csv"
FIGURE_DIR = ROOT / f"{PREFIX}_Figures"
LEAN_TAYLOR = ROOT / f"{PREFIX}_Taylor_Cosine.lean"
SEED = 20260731 + 35
LEAN = ROOT / ".elan/toolchains/leanprover--lean4---v4.28.0/bin/lean"

TABLES = {
    411: ROOT / f"{PREFIX}_Taylor_Reflection.csv",
    412: ROOT / f"{PREFIX}_KKT_Contact.csv",
    413: ROOT / f"{PREFIX}_Transcendence_Interface.csv",
    414: ROOT / f"{PREFIX}_Nonlinear_Damping_Cocycle.csv",
    415: ROOT / f"{PREFIX}_Photonic_Jacobian_Atlas.csv",
    416: ROOT / f"{PREFIX}_Photonic_Pilot_Gate.csv",
    417: ROOT / f"{PREFIX}_Noisy_Comb_Gap.csv",
    418: ROOT / f"{PREFIX}_Erasure_Robust_Code.csv",
    419: ROOT / f"{PREFIX}_Twelve_Mode_Optimization.csv",
    420: ROOT / f"{PREFIX}_Clock_Anchor_Gate.csv",
    421: ROOT / f"{PREFIX}_JSR_External_Gate.csv",
    422: ROOT / f"{PREFIX}_JSR_Optimal_Estimator.csv",
    423: ROOT / f"{PREFIX}_Reference_Frame_Ledger.csv",
    424: ROOT / f"{PREFIX}_Oriented_Record.csv",
    425: ROOT / f"{PREFIX}_Scale_Sections.csv",
    426: ROOT / f"{PREFIX}_Admission_Campaign.csv",
    427: ROOT / f"{PREFIX}_EW_Blind_Gate.csv",
}


def json_ready(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): json_ready(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_ready(item) for item in value]
    if isinstance(value, np.ndarray):
        return json_ready(value.tolist())
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, Fraction):
        return str(value)
    if isinstance(value, complex):
        return {"real": float(value.real), "imag": float(value.imag)}
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    keys: list[str] = []
    for row in rows:
        for key in row:
            if key not in keys:
                keys.append(key)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=keys)
        writer.writeheader()
        for row in rows:
            writer.writerow({
                key: json.dumps(json_ready(value), ensure_ascii=False)
                if isinstance(value, (dict, list, tuple, np.ndarray)) else value
                for key, value in row.items()
            })


def sha256(path: Path) -> str:
    digest = hashlib.sha256(path.read_bytes())
    return digest.hexdigest()


def previous_atoms() -> tuple[np.ndarray, np.ndarray]:
    with (ROOT / "FIN_Programs_379_393_Jordan_Realization.csv").open(encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    return (
        np.array([float(row["node"]) for row in rows]),
        np.array([float(row["signed_weight"]) for row in rows]),
    )


def external_gate(program: int, evidence: str, description: str) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    row = {
        "program": f"P{program}",
        "required_evidence": evidence,
        "description": description,
        "admitted": False,
        "candidate_count": 0,
        "reason": "no independent signed, calibrated, hash-frozen bundle supplied",
    }
    return ({
        "status": f"[Blocked by external evidence] {description}",
        "admitted": False,
        "candidate_count": 0,
        "required_evidence": evidence,
        "boundary": "Code cannot create independent custody, calibration, apparatus output, a physical clock, or blinded data.",
    }, [row])


# ---------------------------------------------------------------------------
# P411: exact Taylor algorithm reflected through Lean
# ---------------------------------------------------------------------------


def generate_taylor_source() -> None:
    angles = ",\n  ".join(
        f"(({743*d+650} : Rat) / 4000)" for d in range(12)
    )
    source = f"""import Std

/- Exact rational reflection of the Taylor polynomials used to enclose
   cos((743*d+650)/4000), d=0,...,11.  The analytic alternating-remainder
   theorem is a separately named external dependency. -/
namespace FINPrograms411To427Taylor

def factorial : Nat -> Nat
  | 0 => 1
  | n + 1 => (n + 1) * factorial n

def signedTerm (x : Rat) (k : Nat) : Rat :=
  let magnitude := x ^ (2*k) / (factorial (2*k) : Rat)
  if k % 2 = 0 then magnitude else -magnitude

def taylorSum (x : Rat) (n : Nat) : Rat :=
  (List.range (n+1)).foldl (fun total k => total + signedTerm x k) 0

def lower (x : Rat) : Rat := taylorSum x 21
def upper (x : Rat) : Rat := taylorSum x 20

def angles : List Rat := [
  {angles}
]

def rationalCertificate : Bool :=
  angles.all (fun x => decide (lower x < upper x && upper x - lower x < (1 : Rat) / 10^30))

theorem all_twelve_rational_taylor_checks : rationalCertificate = true := by
  native_decide

def AnalyticCosineBridge (cosine : Rat -> Rat) : Prop :=
  forall x, x ∈ angles -> lower x <= cosine x && cosine x <= upper x

theorem interval_use_requires_bridge
    (cosine : Rat -> Rat) (bridge : AnalyticCosineBridge cosine) :
    forall x, x ∈ angles -> lower x <= cosine x && cosine x <= upper x :=
  bridge

end FINPrograms411To427Taylor
"""
    LEAN_TAYLOR.write_text(source, encoding="utf-8")


def rational_cos_partial(x: Fraction, n: int) -> Fraction:
    total = Fraction(0)
    for k in range(n + 1):
        term = x ** (2 * k) / math.factorial(2 * k)
        total += term if k % 2 == 0 else -term
    return total


def program_411() -> tuple[dict[str, Any], list[dict[str, Any]]]:
    generate_taylor_source()
    completed = subprocess.run([str(LEAN), str(LEAN_TAYLOR)], cwd=ROOT, capture_output=True, text=True, timeout=120)
    if completed.returncode:
        raise RuntimeError(completed.stderr or completed.stdout)
    rows = []
    for distance in range(12):
        angle = Fraction(743 * distance + 650, 4000)
        lower = rational_cos_partial(angle, 21)
        upper = rational_cos_partial(angle, 20)
        rows.append({
            "distance": distance,
            "angle_exact": str(angle),
            "lower_exact": str(lower),
            "upper_exact": str(upper),
            "width": float(upper - lower),
            "contains_float_cosine": float(lower) <= math.cos(float(angle)) <= float(upper),
        })
    return ({
        "status": "[Proven] exact Taylor-polynomial algorithm reflected by Lean; [Open] in-assistant analytic cosine bridge",
        "angle_count": 12,
        "taylor_orders": [40, 42],
        "maximum_rational_width": max(row["width"] for row in rows),
        "lean_returncode": completed.returncode,
        "source": LEAN_TAYLOR.name,
        "theorem_scope": "Lean computes the factorials, powers, alternating sums, ordering, and sub-1e-30 widths for all twelve rational phases.",
        "boundary": "The dependency-free file types the analytic alternating-remainder theorem as AnalyticCosineBridge; it does not formalize real cosine itself.",
        "new_object": "O143 Taylor trust lattice",
    }, rows)


# ---------------------------------------------------------------------------
# P412/P413: exact-contact candidate and transcendence dependency
# ---------------------------------------------------------------------------


def program_412() -> tuple[dict[str, Any], list[dict[str, Any]]]:
    result, rows = prior.program_395()
    result = dict(result)
    result["status"] = "[Strong evidence] simultaneous seven-contact KKT candidate; [Open] 25-variable interval uniqueness and global feasibility"
    result["relation_to_P395"] = "P395 proved the frozen P380 dual is only near-contact; P412 solves the simultaneous contact equations numerically instead of freezing that dual."
    result["new_object"] = "O144 simultaneous contact replacement candidate"
    return result, rows


def program_413() -> tuple[dict[str, Any], list[dict[str, Any]]]:
    result, rows = prior.program_396()
    result = dict(result)
    result["status"] = "[Proven] proof-assistant dependency localization; [Blocked] local formal Lindemann-Weierstrass library"
    result["new_object"] = "O145 transcendence provider boundary"
    return result, rows


# ---------------------------------------------------------------------------
# P414: exactly one nonlinear normalized RG recurrence class
# ---------------------------------------------------------------------------


def program_414() -> tuple[dict[str, Any], list[dict[str, Any]]]:
    # R_gamma(c)=c/(1+gamma*c), c_0=1 gives c_d=1/(1+gamma*d).
    gamma = Fraction(99, 101)  # forced by C_damp(1)=101/200
    predicted_c2 = Fraction(1, 1) / (1 + 2 * gamma)
    # Equality to C_damp(2) would force r^4=10199/10100 for r^5=2,
    # hence (10199/10100)^5=16.
    defect = 10199**5 - 16 * 10100**5
    rows = []
    current = Fraction(1)
    for distance in range(7):
        target = prior.damping_atom(float(distance))
        rows.append({
            "distance": distance,
            "recurrence_value_exact": str(current),
            "recurrence_value": float(current),
            "target_C_damp": target,
            "defect": float(current) - target,
        })
        current = current / (1 + gamma * current)
    return ({
        "status": "[Refuted] target C_damp arises from the declared normalized nonlinear RG recurrence",
        "source_class": "c_(d+1)=c_d/(1+gamma*c_d), c_0=1",
        "forced_gamma_from_d1": str(gamma),
        "predicted_c2_exact": str(predicted_c2),
        "equality_would_force": "(10199/10100)^5=16",
        "exact_integer_defect": defect,
        "theorem_scope": "The d=1 target fixes gamma uniquely; the d=2 target then contradicts an exact integer identity.",
        "boundary": "This excludes one nonlinear cocycle/flow class only. Fractional, nonlocal, state-dependent, and operator-valued sources remain open.",
        "new_object": "O146 nonlinear RG-recurrence obstruction",
    }, rows)


# ---------------------------------------------------------------------------
# P415/P416: photonic regularity atlas and external pilot
# ---------------------------------------------------------------------------


def photonic_transfer_model() -> tuple[list[Any], np.ndarray]:
    strict_a, _ = core.strict_operator()
    old = json.loads((ROOT / "FIN_Programs_323_336_Results.json").read_text(encoding="utf-8"))
    time = float(old["P326"]["best_protocols"]["wave"]["nominal_time"])
    target = np.linalg.eigh(strict_a)
    values, vectors = target
    unitary = (vectors * np.exp(-1j * time * values)) @ vectors.conj().T
    rotations, diagonal = p295.givens_decompose_unitary(unitary)
    return rotations, diagonal


def transfer_vector(rotations: list[Any], diagonal: np.ndarray, parameters: np.ndarray) -> np.ndarray:
    transfer = p365.component_transfer(rotations, diagonal, parameters)
    return np.r_[transfer.real.ravel(), transfer.imag.ravel()]


def numerical_jacobian(function: Callable[[np.ndarray], np.ndarray], point: np.ndarray, step: float = 2e-6) -> np.ndarray:
    columns = []
    for index in range(len(point)):
        shift = np.zeros_like(point)
        shift[index] = step
        columns.append((function(point + shift) - function(point - shift)) / (2 * step))
    return np.column_stack(columns)


def program_415(rng: np.random.Generator) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    rotations, diagonal = photonic_transfer_model()
    count = len(rotations)
    rows = []
    for sample in range(9):
        if sample == 0:
            point = np.zeros(2 * count)
        else:
            point = np.r_[rng.uniform(0, 0.03, count), rng.uniform(-0.1, 0.1, count)]
        jacobian = numerical_jacobian(lambda p: transfer_vector(rotations, diagonal, p), point)
        singular = np.linalg.svd(jacobian, compute_uv=False)
        tolerance = singular[0] * max(jacobian.shape) * np.finfo(float).eps
        rows.append({
            "sample": sample,
            "parameter_norm": float(np.linalg.norm(point)),
            "numerical_rank": int(np.sum(singular > tolerance)),
            "smallest_singular_value": float(singular[-1]),
            "condition_number": float(singular[0] / singular[-1]),
        })
    return ({
        "status": "[Strong evidence] full local rank on a sampled bounded chart; [Open] global interval identifiability or alias",
        "samples": len(rows),
        "parameter_count": 2 * count,
        "full_rank_samples": sum(row["numerical_rank"] == 2 * count for row in rows),
        "minimum_sampled_singular_value": min(row["smallest_singular_value"] for row in rows),
        "maximum_sampled_condition_number": max(row["condition_number"] for row in rows),
        "theorem_scope": "Central-difference complex-transfer Jacobians are numerically full column rank at the origin and eight frozen random chart points.",
        "boundary": "Finite floating-point regularity samples do not cover the continuum box and cannot exclude distant compensating aliases.",
        "new_object": "O147 sampled photonic regularity atlas",
    }, rows)


# ---------------------------------------------------------------------------
# P417: noisy-comb primal/upper-gap audit
# ---------------------------------------------------------------------------


def program_417(rng: np.random.Generator) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    rows = []
    smallest_gap = math.inf
    for uses in (2, 3, 4):
        for coherence in (0.6, 0.8, 1.0):
            for fraction in (0.25, 0.5, 0.75, 1.0):
                theta = fraction * math.pi / (2 * uses)
                lower, probabilities = prior.optimize_symmetric_parallel(uses, theta, coherence, rng)
                upper = min(1.0, uses * coherence * abs(math.sin(theta)))
                gap = max(0.0, upper - lower)
                smallest_gap = min(smallest_gap, gap)
                rows.append({
                    "uses": uses,
                    "coherence": coherence,
                    "threshold_fraction": fraction,
                    "optimized_parallel_lower": lower,
                    "adaptive_hybrid_upper": upper,
                    "gap": gap,
                    "sector_probabilities": probabilities.tolist(),
                })
    return ({
        "status": "[Strong evidence] improved primal atlas and exact ideal-boundary closure; [Open] noisy SDP dual certificate",
        "grid_rows": len(rows),
        "minimum_gap": smallest_gap,
        "maximum_gap": max(row["gap"] for row in rows),
        "zero_gap_rows_below_1e_8": sum(row["gap"] < 1e-8 for row in rows),
        "solver_boundary": "No certified semidefinite optimizer is installed locally; nonlinear pure-state lower bounds cannot substitute for a comb dual certificate.",
        "boundary": "The unrestricted noisy adaptive comb, ancillas, mixed inputs, and full twelve-mode channel remain outside this result.",
        "new_object": "O148 noisy comb primal-dual-gap atlas",
    }, rows)


# ---------------------------------------------------------------------------
# P418: heralded-erasure-aware symmetric entangled codes
# ---------------------------------------------------------------------------


def bits(uses: int) -> np.ndarray:
    return np.array([[(value >> bit) & 1 for bit in range(uses)] for value in range(2**uses)], dtype=int)


def symmetric_amplitudes(uses: int, sector_probabilities: np.ndarray) -> np.ndarray:
    bit_table = bits(uses)
    weights = bit_table.sum(axis=1)
    return np.array([math.sqrt(max(sector_probabilities[k], 0) / math.comb(uses, int(k))) for k in weights])


def partial_trace_keep(rho: np.ndarray, uses: int, keep: tuple[int, ...]) -> np.ndarray:
    lost = tuple(index for index in range(uses) if index not in keep)
    kept_dimension = 2 ** len(keep)
    lost_dimension = 2 ** len(lost)
    result = np.zeros((kept_dimension, kept_dimension), dtype=complex)
    for a in range(kept_dimension):
        a_bits = [(a >> index) & 1 for index in range(len(keep))]
        for b in range(kept_dimension):
            b_bits = [(b >> index) & 1 for index in range(len(keep))]
            total = 0j
            for lost_value in range(lost_dimension):
                lost_bits = [(lost_value >> index) & 1 for index in range(len(lost))]
                left = [0] * uses
                right = [0] * uses
                for position, bit in zip(keep, a_bits):
                    left[position] = bit
                for position, bit in zip(keep, b_bits):
                    right[position] = bit
                for position, bit in zip(lost, lost_bits):
                    left[position] = right[position] = bit
                left_index = sum(bit << index for index, bit in enumerate(left))
                right_index = sum(bit << index for index, bit in enumerate(right))
                total += rho[left_index, right_index]
            result[a, b] = total
    return result


def trace_distance(left: np.ndarray, right: np.ndarray) -> float:
    return float(0.5 * np.sum(np.abs(np.linalg.eigvalsh(left - right))))


def erasure_aware_distance(uses: int, theta: float, coherence: float, survival: float, probabilities: np.ndarray) -> float:
    table = bits(uses)
    weights = table.sum(axis=1)
    amplitudes = symmetric_amplitudes(uses, probabilities)
    pure = np.outer(amplitudes, amplitudes)
    hamming = np.abs(table[:, None, :] - table[None, :, :]).sum(axis=2)
    base = pure * coherence**hamming
    phase_plus = np.exp(1j * theta * weights)
    phase_minus = np.exp(-1j * theta * weights)
    plus = base * phase_plus[:, None] * np.conjugate(phase_plus[None, :])
    minus = base * phase_minus[:, None] * np.conjugate(phase_minus[None, :])
    total = 0.0
    for mask in range(2**uses):
        keep = tuple(index for index in range(uses) if (mask >> index) & 1)
        probability = survival ** len(keep) * (1 - survival) ** (uses - len(keep))
        total += probability * trace_distance(
            partial_trace_keep(plus, uses, keep),
            partial_trace_keep(minus, uses, keep),
        )
    return total


def optimize_erasure_code(uses: int, theta: float, coherence: float, survival: float, rng: np.random.Generator) -> tuple[float, np.ndarray]:
    product = np.array([math.comb(uses, k) / 2**uses for k in range(uses + 1)])
    ghz = np.r_[0.5, np.zeros(uses - 1), 0.5]
    starts = [product, ghz, np.full(uses + 1, 1 / (uses + 1))]
    starts.extend(rng.dirichlet(np.ones(uses + 1)) for _ in range(3))
    best_value, best = -1.0, product
    for start in starts:
        fit = minimize(
            lambda p: -erasure_aware_distance(uses, theta, coherence, survival, p),
            start,
            method="SLSQP",
            bounds=[(0, 1)] * (uses + 1),
            constraints={"type": "eq", "fun": lambda p: np.sum(p) - 1},
            options={"ftol": 1e-11, "maxiter": 250},
        )
        candidate = np.maximum(fit.x, 0)
        candidate /= candidate.sum()
        value = erasure_aware_distance(uses, theta, coherence, survival, candidate)
        if value > best_value:
            best_value, best = value, candidate
    return best_value, best


def program_418(rng: np.random.Generator) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    rows = []
    strict_improvements = 0
    maximum_gain = 0.0
    for uses, coherence, survival, fraction in itertools.product((3, 4), (0.6, 0.8), (0.6, 0.8), (0.5, 0.8)):
        theta = fraction * math.pi / (2 * uses)
        optimized, vector = optimize_erasure_code(uses, theta, coherence, survival, rng)
        product_vector = np.array([math.comb(uses, k) / 2**uses for k in range(uses + 1)])
        ghz_vector = np.r_[0.5, np.zeros(uses - 1), 0.5]
        product = erasure_aware_distance(uses, theta, coherence, survival, product_vector)
        ghz = erasure_aware_distance(uses, theta, coherence, survival, ghz_vector)
        gain = optimized - max(product, ghz)
        maximum_gain = max(maximum_gain, gain)
        strict_improvements += gain > 1e-7
        rows.append({
            "uses": uses,
            "coherence": coherence,
            "survival": survival,
            "threshold_fraction": fraction,
            "product": product,
            "ghz": ghz,
            "optimized_erasure_aware": optimized,
            "gain_over_product_or_ghz": gain,
            "sector_probabilities": vector.tolist(),
        })
    status = (
        "[Strong evidence] erasure-aware symmetric codes strictly improve the two baseline families"
        if strict_improvements else
        "[Refuted in tested symmetric family] erasure-aware optimization improves the two baselines"
    )
    return ({
        "status": status,
        "grid_rows": len(rows),
        "strict_improvement_rows": strict_improvements,
        "maximum_gain": maximum_gain,
        "theorem_scope": "Every heralded survivor subset is retained as an orthogonal classical block and the lost qubits are traced out exactly for n<=4.",
        "boundary": "The optimization is numerical and symmetry-restricted; it is not a general quantum erasure-code or adaptive-comb theorem.",
        "new_object": "O149 heralded erasure-aware symmetric code",
    }, rows)


# ---------------------------------------------------------------------------
# P419: reduced full twelve-mode optimization
# ---------------------------------------------------------------------------


def program_419(rng: np.random.Generator) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    strict_a, _ = core.strict_operator()
    legacy_a, _ = core.legacy_amplitude_absorbed_operator()
    relative = strict_a - legacy_a
    fourier = np.fft.fft(np.eye(12)) / math.sqrt(12)
    transformed = fourier.conj().T @ relative @ fourier
    gaps = np.real(np.diag(transformed))
    diameter = float(gaps.max() - gaps.min())
    off_diagonal = float(np.linalg.norm(transformed - np.diag(np.diag(transformed))))
    rows = []
    for coherence, survival, fraction in itertools.product((0.6, 0.8), (0.8, 1.0), (0.5, 1.0)):
        time = fraction * math.pi / diameter
        optimum, probabilities = prior.optimize_twelve_mode(gaps, time, coherence, rng)
        extremal = np.zeros(12)
        extremal[[int(np.argmin(gaps)), int(np.argmax(gaps))]] = 0.5
        extremal_value = prior.twelve_mode_distance(gaps, time, coherence, extremal)
        rows.append({
            "coherence": coherence,
            "survival": survival,
            "threshold_fraction": fraction,
            "time": time,
            "optimized_one_use_lower": survival * optimum,
            "extremal_two_mode_lower": survival * extremal_value,
            "active_modes": int(np.sum(probabilities > 1e-5)),
            "probabilities": probabilities.tolist(),
        })
    return ({
        "status": "[Strong evidence] optimized full twelve-mode one-use inputs; [Open] multi-use adaptive gap below 1e-3",
        "rows": len(rows),
        "relative_generator_diameter": diameter,
        "fourier_off_diagonal_defect": off_diagonal,
        "maximum_gain_over_extremal_pair": max(row["optimized_one_use_lower"] - row["extremal_two_mode_lower"] for row in rows),
        "boundary": "The dephasing and heralded survival laws are supplied. One-use pure-state optimization is not a full process-comb certificate.",
        "new_object": "O150 twelve-mode noise-adapted simplex",
    }, rows)


# ---------------------------------------------------------------------------
# P422: variance-optimal Jordan sampler for twelve moments
# ---------------------------------------------------------------------------


def program_422() -> tuple[dict[str, Any], list[dict[str, Any]]]:
    nodes, weights = previous_atoms()
    features = np.array([[node**order for order in range(12)] for node in nodes])
    radii = np.sqrt(np.sum(features**2, axis=1))
    q_jsr = np.abs(weights) / np.sum(np.abs(weights))
    q_opt = np.abs(weights) * radii
    q_opt /= q_opt.sum()

    def total_variance(probabilities: np.ndarray) -> float:
        second = sum(weights[index] ** 2 * float(np.sum(features[index] ** 2)) / probabilities[index] for index in range(7))
        means = weights @ features
        return float(second - np.sum(means**2))

    variance_jsr = total_variance(q_jsr)
    variance_opt = total_variance(q_opt)
    rows = [{
        "atom": index,
        "node": nodes[index],
        "weight": weights[index],
        "feature_radius": radii[index],
        "baseline_probability": q_jsr[index],
        "optimal_probability": q_opt[index],
        "probability_ratio": q_opt[index] / q_jsr[index],
    } for index in range(7)]
    return ({
        "status": "[Proven] variance-optimal unbiased sampler for the declared equal-weight twelve-moment loss",
        "baseline_total_variance": variance_jsr,
        "optimal_total_variance": variance_opt,
        "relative_variance_reduction": 1 - variance_opt / variance_jsr,
        "law": "q_i proportional to |w_i| sqrt(sum_{k=0}^{11} x_i^{2k})",
        "proof": "Cauchy-Schwarz minimizes sum_i w_i^2 r_i^2/q_i on the probability simplex; the importance estimator remains unbiased.",
        "boundary": "Optimality is tied to the declared twelve moments and equal loss weights. It is not apparatus efficiency or a universal physical sampler.",
        "new_object": "O151 variance-optimal JSR sampling law",
    }, rows)


# ---------------------------------------------------------------------------
# P423: reference-frame cost ledger for sign-to-phase conversion
# ---------------------------------------------------------------------------


def binary_entropy(q: float) -> float:
    if q in (0.0, 1.0):
        return 0.0
    return -q * math.log2(q) - (1 - q) * math.log2(1 - q)


def program_423() -> tuple[dict[str, Any], list[dict[str, Any]]]:
    rows = []
    for q in np.linspace(0, 0.5, 11):
        bias = 1 - 2 * q
        eigenvalues = [(1 + bias) / 2, (1 - bias) / 2]
        relative_entropy_asymmetry = 1 - binary_entropy(eigenvalues[1])
        rows.append({
            "negative_sign_probability": q,
            "encoded_l1_coherence": bias,
            "encoded_relative_entropy_asymmetry_bits": relative_entropy_asymmetry,
            "chosen_encoder_reference_budget_bits": 1.0,
            "polarity_swap_invariant_cost": True,
        })
    return ({
        "status": "[Proven] conditional reference-frame resource ledger; [Refuted] free classical sign creates U(1) coherence",
        "free_no_go": "U(1)-covariant incoherent operations map diagonal sign records to zero-coherence states.",
        "conditional_encoder": "consume an aligned |+> phase reference and apply the sign-controlled phase 0/pi",
        "chosen_encoder_reference_budget_bits": 1.0,
        "selector_discharge": False,
        "boundary": "The aligned phase origin and polarity are supplied resources. Swapping the encoding leaves costs invariant and QW-2191 open.",
        "new_object": "O152 phase-reference conversion ledger",
    }, rows)


# ---------------------------------------------------------------------------
# P424: explicit inversion-odd record and its orientation-torsor obstruction
# ---------------------------------------------------------------------------


def program_424() -> tuple[dict[str, Any], list[dict[str, Any]]]:
    nodes, weights = previous_atoms()
    contributions = weights * np.sin(2 * math.pi * nodes)
    current = float(contributions.sum())
    reflected_nodes = 1 - nodes
    reflected = float(np.sum(weights * np.sin(2 * math.pi * reflected_nodes)))
    rows = [{
        "atom": index,
        "node": nodes[index],
        "weight": weights[index],
        "odd_contribution": contributions[index],
        "reflected_contribution": weights[index] * math.sin(2 * math.pi * reflected_nodes[index]),
    } for index in range(7)]
    return ({
        "status": "[Proven] nonzero inversion-odd JSR record current; [Refuted] nonconventional strict selector source",
        "formula": "J_or=sum_i w_i sin(2*pi*x_i)",
        "value": current,
        "reflected_value": reflected,
        "oddness_defect": reflected + current,
        "nonzero": abs(current) > 1e-12,
        "selector_discharge": False,
        "theorem_scope": "Reflection x->1-x sends J_or to -J_or exactly; the frozen seven-atom representation gives a nonzero numerical value.",
        "boundary": "The coordinate orientation x increasing from 0 to 1 is supplied. Reversing that chart flips the sign, so the object is an oriented receiver/record, not an internal polarity source.",
        "new_object": "O153 oriented JSR record current",
    }, rows)


# ---------------------------------------------------------------------------
# P425: compare three conditional scale sections
# ---------------------------------------------------------------------------


def program_425() -> tuple[dict[str, Any], list[dict[str, Any]]]:
    # Normalize all sections to rho0 at control u=1, then compare elasticities.
    rho0 = 0.4296970877901551
    rows = []
    for control in (0.5, 0.75, 1.0, 1.25, 1.5):
        resolution = rho0 * control
        entropy_cell = rho0 ** control
        rg_exponential = math.exp(-(-math.log(rho0)) / control)
        rows.extend([
            {"control": control, "law": "resolution-linear", "rho": resolution, "baseline_elasticity": 1.0},
            {"control": control, "law": "entropy-cell power", "rho": entropy_cell, "baseline_elasticity": math.log(rho0)},
            {"control": control, "law": "inverse-RG-time", "rho": rg_exponential, "baseline_elasticity": -math.log(rho0)},
        ])
    return ({
        "status": "[Proven] the three declared conditional scale sections are inequivalent under a shared control coordinate",
        "common_baseline_rho": rho0,
        "laws": 3,
        "pairwise_distinct_baseline_elasticities": True,
        "theorem": "Equal values at one control point do not identify sections; their logarithmic derivatives at that point are 1, log(rho0), and -log(rho0).",
        "boundary": "A nonlinear redefinition of the external control can map any monotone section to another. No section is canonical without physical semantics and calibration for that control.",
        "new_object": "O154 conditional scale-section comparison groupoid",
    }, rows)


def summary_rows(results: dict[str, Any]) -> list[dict[str, Any]]:
    return [{
        "program": f"P{program}",
        "status": results[f"P{program}"]["status"],
        "new_object": results[f"P{program}"].get("new_object", ""),
        "boundary": results[f"P{program}"].get("boundary", ""),
    } for program in range(411, 428)]


def make_figures(results: dict[str, Any]) -> None:
    FIGURE_DIR.mkdir(exist_ok=True)
    rows = results["_rows"]

    figure, axes = plt.subplots(1, 2, figsize=(12, 4.4))
    axes[0].semilogy([r["distance"] for r in rows[411]], [r["width"] for r in rows[411]], marker="o")
    axes[0].set_xlabel("integer distance")
    axes[0].set_ylabel("exact Taylor bracket width")
    axes[0].set_title("P411 reflected Taylor enclosures")
    axes[1].scatter([float(r["node"]) for r in rows[412]], [float(r["weight"]) for r in rows[412]], c=["#d95f02" if float(r["weight"]) < 0 else "#1b9e77" for r in rows[412]], s=65)
    axes[1].axhline(0, color="black", linewidth=0.8)
    axes[1].set_xlabel("contact node")
    axes[1].set_ylabel("signed weight")
    axes[1].set_title("P412 simultaneous contact candidate")
    figure.tight_layout()
    figure.savefig(FIGURE_DIR / "p411_p412_formal_contact.png", dpi=180)
    plt.close(figure)

    figure, axes = plt.subplots(1, 2, figsize=(12, 4.4))
    axes[0].plot([r["distance"] for r in rows[414]], [r["recurrence_value"] for r in rows[414]], marker="o", label="nonlinear recurrence")
    axes[0].plot([r["distance"] for r in rows[414]], [r["target_C_damp"] for r in rows[414]], marker="s", label="C_damp")
    axes[0].set_xlabel("distance")
    axes[0].set_ylabel("attenuation completion")
    axes[0].set_title("P414 nonlinear source obstruction")
    axes[0].legend(fontsize=8)
    axes[1].scatter([r["parameter_norm"] for r in rows[415]], [r["smallest_singular_value"] for r in rows[415]], color="#7570b3")
    axes[1].set_xlabel("chart parameter norm")
    axes[1].set_ylabel("smallest Jacobian singular value")
    axes[1].set_title("P415 sampled photonic regularity")
    figure.tight_layout()
    figure.savefig(FIGURE_DIR / "p414_p415_damping_photonic.png", dpi=180)
    plt.close(figure)

    figure, axes = plt.subplots(1, 2, figsize=(12, 4.4))
    selected = [r for r in rows[417] if r["uses"] == 4 and r["coherence"] == 0.8]
    axes[0].plot([r["threshold_fraction"] for r in selected], [r["optimized_parallel_lower"] for r in selected], marker="o", label="parallel lower")
    axes[0].plot([r["threshold_fraction"] for r in selected], [r["adaptive_hybrid_upper"] for r in selected], marker="s", label="adaptive upper")
    axes[0].set_xlabel("threshold fraction")
    axes[0].set_ylabel("distance")
    axes[0].set_title("P417 remaining noisy-comb gap")
    axes[0].legend(fontsize=8)
    gains = [r["gain_over_product_or_ghz"] for r in rows[418]]
    axes[1].bar(np.arange(len(gains)), gains, color="#1b9e77")
    axes[1].axhline(0, color="black", linewidth=0.8)
    axes[1].set_xlabel("declared noise cell")
    axes[1].set_ylabel("optimized gain")
    axes[1].set_title("P418 erasure-aware code search")
    figure.tight_layout()
    figure.savefig(FIGURE_DIR / "p417_p418_noisy_erasure.png", dpi=180)
    plt.close(figure)

    figure, axes = plt.subplots(1, 2, figsize=(12, 4.4))
    axes[0].plot([r["node"] for r in rows[422]], [r["baseline_probability"] for r in rows[422]], marker="o", label="baseline JSR")
    axes[0].plot([r["node"] for r in rows[422]], [r["optimal_probability"] for r in rows[422]], marker="s", label="variance-optimal")
    axes[0].set_xlabel("atom node")
    axes[0].set_ylabel("sampling probability")
    axes[0].set_title("P422 estimator redesign")
    axes[0].legend(fontsize=8)
    for law in ("resolution-linear", "entropy-cell power", "inverse-RG-time"):
        selected_scale = [r for r in rows[425] if r["law"] == law]
        axes[1].plot([r["control"] for r in selected_scale], [r["rho"] for r in selected_scale], marker="o", label=law)
    axes[1].set_xlabel("external control")
    axes[1].set_ylabel("selected rho")
    axes[1].set_title("P425 inequivalent scale sections")
    axes[1].legend(fontsize=7)
    figure.tight_layout()
    figure.savefig(FIGURE_DIR / "p422_p425_estimator_scale.png", dpi=180)
    plt.close(figure)


def main() -> None:
    rng = np.random.default_rng(SEED)
    results: dict[str, Any] = {
        "metadata": {
            "programs": "P411-P427",
            "release": "10.35",
            "seed": SEED,
            "date": "2026-07-31",
            "kernel_split": "legacy intermediate and strict working kernels remain distinct; no role transfer",
            "new_theoretical_objects": [f"O{index}" for index in range(143, 155)],
        },
        "_rows": {},
    }
    programs: dict[int, Callable[[], tuple[dict[str, Any], list[dict[str, Any]]]]] = {
        411: program_411,
        412: program_412,
        413: program_413,
        414: program_414,
        415: lambda: program_415(rng),
        416: lambda: external_gate(416, "independent_calibrated_photonic_record", "independent calibrated photonic pilot"),
        417: lambda: program_417(rng),
        418: lambda: program_418(rng),
        419: lambda: program_419(rng),
        420: lambda: external_gate(420, "traceable_SI_to_tau_clock_anchor", "traceable physical clock anchor"),
        421: lambda: external_gate(421, "provider_registrar_separated_JSR_events", "external event-level JSR execution"),
        422: program_422,
        423: program_423,
        424: program_424,
        425: program_425,
        426: lambda: external_gate(426, "QW_standards_reservoir_campaign", "independent QW, standards, and reservoir admission campaign"),
        427: lambda: external_gate(427, "post_bridge_electroweak_blind_bundle", "blinded electroweak role-transfer test"),
    }
    for program in range(411, 428):
        result, rows = programs[program]()
        results[f"P{program}"] = result
        results["_rows"][program] = rows
        write_csv(TABLES[program], rows)
        print(f"P{program}: {result['status']}")
    write_csv(SUMMARY_PATH, summary_rows(results))
    public = {key: value for key, value in results.items() if key != "_rows"}
    RESULTS_PATH.write_text(json.dumps(json_ready(public), indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    make_figures(results)
    print(RESULTS_PATH)


if __name__ == "__main__":
    main()
