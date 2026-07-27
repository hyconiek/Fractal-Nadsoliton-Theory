#!/usr/bin/env python3
"""FIN Programs 217--229: operational sections and spectral tomography.

The suite executes the thirteen studies recommended by Release 10.19.
Analytic theorems, exact finite checks, machine-compiled Lean fragments,
synthetic method evidence, unavailable Arb infrastructure, and external-data
gates are reported as distinct claim classes.
"""

from __future__ import annotations

from fractions import Fraction
from pathlib import Path
import hashlib
import importlib.util
import json
import math
import os
import platform
import shutil
import subprocess

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.linalg import expm
from scipy.stats import beta as beta_distribution
from scipy.stats import chi2, norm

import fin_programs_204_216_categorical_catalytic_external as previous


ROOT = Path(__file__).resolve().parent
FIG = ROOT / "FIN_Programs_217_229_Operational_Section_Spectral_Tomography_Figures"
RESULTS = ROOT / "FIN_Programs_217_229_Operational_Section_Spectral_Tomography_Results.json"
ARB_ATTESTATION = ROOT / "FIN_Programs_217_229_Arb_Execution_Attestation.json"
LEAN_BUILD = ROOT / "FIN_Programs_217_229_Lean_Build_Record.json"
PHASE_FORMAL = ROOT / "FIN_Programs_217_229_Phase_Formalization_Record.json"
EXTERNAL_REQUEST = ROOT / "FIN_Programs_217_229_Independent_Custody_Bundle_Request.json"
PREDICTION_LOCK = ROOT / "FIN_Programs_217_229_Semigroup_Holdout_Preregistration.json"
PREVIOUS_RESULTS = ROOT / "FIN_Programs_204_216_Categorical_Catalytic_External_Results.json"
PREVIOUS_ARB = ROOT / "FIN_Programs_204_216_Arb_Environment_Contract.json"
GENERAL_LEAN = ROOT / "FIN_Programs_204_216_Dirichlet_Library.lean"
DIRICHLET_CORE = ROOT / "FIN_Programs_217_229_Dirichlet_Core.lean"
PHASE_CORE = ROOT / "FIN_Programs_217_229_Phase_Torsion_Core.lean"

SEED = 20260728


def rng_for(offset: int) -> np.random.Generator:
    return np.random.default_rng(SEED + offset)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def canonical_digest(obj: object) -> str:
    payload = json.dumps(
        obj, sort_keys=True, separators=(",", ":"), ensure_ascii=False
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def write_canonical(path: Path, core: dict) -> str:
    record = {"core": core, "canonical_core_sha256": canonical_digest(core)}
    path.write_text(json.dumps(record, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return sha256(path)


def program217_operational_central_state() -> dict:
    rng = rng_for(217)
    m = 5
    uniform = np.full(m, 1.0 / m)
    max_invariance_residual = 0.0
    trials = 400
    for _ in range(trials):
        coefficients = rng.dirichlet(np.ones(12))
        matrix = np.zeros((m, m))
        for coefficient in coefficients:
            permutation = rng.permutation(m)
            matrix += coefficient * np.eye(m)[permutation]
        max_invariance_residual = max(
            max_invariance_residual,
            float(np.max(np.abs(matrix @ uniform - uniform))),
        )

    # Preparation frequencies supply a section operationally but not
    # internally.  The simultaneous Hoeffding band is distribution free.
    true_weights = np.asarray([0.10, 0.15, 0.20, 0.25, 0.30])
    sample_size = 5000
    delta = 0.05
    epsilon = math.sqrt(math.log(2 * m / delta) / (2 * sample_size))
    replicates = 600
    covered = 0
    max_errors = []
    for _ in range(replicates):
        frequency = rng.multinomial(sample_size, true_weights) / sample_size
        error = float(np.max(np.abs(frequency - true_weights)))
        max_errors.append(error)
        covered += error <= epsilon

    dimensions = np.asarray([1, 2, 2, 2, 2], dtype=float)
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.7), constrained_layout=True)
    axes[0].bar(
        np.arange(m) - 0.17,
        uniform,
        width=0.34,
        label="bistochastic invariant",
        color="#1F5A99",
    )
    axes[0].bar(
        np.arange(m) + 0.17,
        true_weights,
        width=0.34,
        label="prepared frequency law",
        color="#D55E00",
    )
    axes[0].set_xticks(range(m), [f"sector {i}" for i in range(m)])
    axes[0].tick_params(axis="x", rotation=25)
    axes[0].set_ylabel("central weight")
    axes[0].legend(fontsize=8)
    axes[1].hist(max_errors, bins=35, color="#6A3D9A", alpha=0.8)
    axes[1].axvline(epsilon, color="black", ls="--", label="Hoeffding radius")
    axes[1].set_xlabel("maximum preparation-frequency error")
    axes[1].legend(fontsize=8)
    fig.suptitle("Program 217: categorical uniform state versus operational frequency section")
    fig.savefig(FIG / "program217_operational_central_state.png", dpi=190)
    plt.close(fig)

    return {
        "constructed_object": {
            "name": "BistochasticCentralStateAndFrequencySection",
            "categorical_part": "uniform state fixed by every bistochastic sector channel",
            "operational_part": "empirical multinomial preparation frequencies with a simultaneous Hoeffding band",
        },
        "number_of_sectors": m,
        "bistochastic_trials": trials,
        "maximum_uniform_invariance_residual": max_invariance_residual,
        "uniqueness_theorem": (
            "Invariance under all permutation channels, a subset of the "
            "bistochastic channels, forces all central coordinates equal."
        ),
        "uniform_dimension_expectation_exact": "9/5",
        "uniform_dimension_expectation": float(uniform @ dimensions),
        "prepared_true_weights": true_weights.tolist(),
        "prepared_dimension_expectation": float(true_weights @ dimensions),
        "sample_size": sample_size,
        "hoeffding_radius": epsilon,
        "replicates": replicates,
        "empirical_simultaneous_coverage": covered / replicates,
        "strict_central_state_source_exported": False,
        "status": "Operational morphism category classified and frequency section constructed",
        "confidence": "Proven theorem and distribution-free band; strong simulation evidence",
        "claim_boundary": (
            "Bistochastic naturality conditionally selects the uniform state. "
            "Nonuniform sector frequencies are supplied by a preparation "
            "record and are not generated by the strict operator."
        ),
    }


def program218_affine_eta_no_go() -> dict:
    degree = 8
    # Polynomial scaling naturality F(2x)=2F(x) has diagonal coefficients
    # (2^j-2)c_j.  Exactly the linear coefficient is unconstrained.
    diagonal = [Fraction(2**j - 2) for j in range(degree + 1)]
    rank = sum(value != 0 for value in diagonal)
    nullity = degree + 1 - rank
    surviving = [j for j, value in enumerate(diagonal) if value == 0]
    x_target = Fraction(4, 5)
    slopes = np.linspace(-4.0, 4.0, 801)
    vector_at_target = slopes * float(x_target)
    fixed = np.isclose(vector_at_target, 0.0, atol=1e-14)

    fig, ax = plt.subplots(figsize=(8.8, 4.8), constrained_layout=True)
    ax.plot(slopes, vector_at_target, color="#1F5A99")
    ax.scatter(slopes[fixed], vector_at_target[fixed], color="#A61B1B", zorder=4)
    ax.axhline(0, color="black", ls="--")
    ax.set_xlabel(r"allowed natural slope $\kappa$")
    ax.set_ylabel(r"$F(4/5)=4\kappa/5$")
    ax.set_title("Program 218: affine-natural eta flow cannot select the nonzero increment")
    ax.grid(True, alpha=0.25)
    fig.savefig(FIG / "program218_affine_eta_no_go.png", dpi=190)
    plt.close(fig)

    return {
        "constructed_object": {
            "name": "AffineNaturalEtaVectorFieldNoGo",
            "coordinate": "x=eta-1",
            "naturality": "F(a*x)=a*F(x) for every a>0, with reflection oddness",
        },
        "theorem": (
            "Every continuous reflection-compatible positively homogeneous "
            "vector field on the eta increment line is F(x)=kappa*x. If "
            "kappa is nonzero its only fixed point is x=0; if kappa=0 every "
            "point is fixed. No unique nonzero fixed point can be selected."
        ),
        "polynomial_degree_audited": degree,
        "polynomial_constraint_rank": rank,
        "polynomial_constraint_nullity": nullity,
        "surviving_monomial_degrees": surviving,
        "strict_increment_exact": "4/5",
        "sampled_slopes": len(slopes),
        "sampled_slopes_fixing_target": int(np.sum(fixed)),
        "target_independent_eta_source_exported": False,
        "status": "Affine-natural eta-flow no-go proved",
        "confidence": "Proven",
        "claim_boundary": (
            "A flow can stabilize 4/5 only after a nonzero offset or equivalent "
            "scale-bearing datum is supplied. The no-go does not exclude laws "
            "with genuinely new strict-derived non-affine data."
        ),
    }


def _docker_probe() -> tuple[bool, str]:
    docker = shutil.which("docker")
    if docker is None:
        return False, "docker CLI absent"
    try:
        proc = subprocess.run(
            [docker, "version", "--format", "{{.Server.Version}}"],
            text=True,
            capture_output=True,
            timeout=10,
            check=False,
        )
        return proc.returncode == 0, (proc.stdout + proc.stderr)[-2000:]
    except subprocess.TimeoutExpired:
        return False, "docker server probe timed out"


def program219_arb_execution() -> dict:
    flint_spec = importlib.util.find_spec("flint")
    arb_binary = shutil.which("arb")
    docker_ok, docker_tail = _docker_probe()
    previous_contract = json.loads(PREVIOUS_ARB.read_text(encoding="utf-8"))
    core = {
        "attestation_id": "FIN-P219-ARB-EXECUTION-001",
        "contract_sha256": sha256(PREVIOUS_ARB),
        "contract_canonical_core_sha256": previous_contract["canonical_core_sha256"],
        "python": platform.python_version(),
        "platform": platform.platform(),
        "python_flint_origin": None if flint_spec is None else flint_spec.origin,
        "arb_binary": arb_binary,
        "docker_server_accessible": docker_ok,
        "docker_probe_tail": docker_tail,
        "all_five_nodes_same_engine": False,
        "formal_run_executed": False,
        "acceptance_passed": False,
        "reason": "no authorized locally callable Arb/FLINT engine",
    }
    attestation_sha = write_canonical(ARB_ATTESTATION, core)

    fig, ax = plt.subplots(figsize=(8.9, 4.8), constrained_layout=True)
    labels = ["python-flint", "Arb binary", "Docker server", "five-node run", "<3% pass"]
    values = [
        flint_spec is not None,
        arb_binary is not None,
        docker_ok,
        False,
        False,
    ]
    ax.bar(labels, [int(v) for v in values], color=["#19733A" if v else "#A61B1B" for v in values])
    ax.set_ylim(0, 1.1)
    ax.set_yticks([0, 1], ["no", "yes"])
    ax.tick_params(axis="x", rotation=18)
    ax.set_title("Program 219: formal execution gate remains closed")
    fig.savefig(FIG / "program219_arb_execution.png", dpi=190)
    plt.close(fig)

    return {
        "constructed_object": {
            "name": "ArbExecutionAttestation",
            "path": ARB_ATTESTATION.name,
            "sha256": attestation_sha,
        },
        "python_flint_available": flint_spec is not None,
        "arb_binary_available": arb_binary is not None,
        "docker_server_accessible": docker_ok,
        "all_five_nodes_same_engine": False,
        "formal_run_executed": False,
        "formal_sub_three_percent_certificate": False,
        "status": "Execution attempted at the environment gate; formal engine unavailable",
        "confidence": "Proven environment attestation",
        "claim_boundary": (
            "No dependency was downloaded and no ordinary-float substitute "
            "was promoted. Program 219 remains an infrastructure blocker, "
            "not a mathematical impossibility theorem."
        ),
    }


def _direct_lean() -> str | None:
    candidates = sorted((ROOT / ".elan" / "toolchains").glob("*/bin/lean"))
    return str(candidates[0]) if candidates else None


def _compile_lean(source: Path) -> tuple[bool, str]:
    lean = _direct_lean()
    if lean is None:
        return False, "no local Lean compiler"
    environment = dict(os.environ)
    environment["LEAN_PATH"] = ""
    try:
        proc = subprocess.run(
            [lean, str(source)],
            cwd=ROOT,
            env=environment,
            text=True,
            capture_output=True,
            timeout=30,
            check=False,
        )
        return proc.returncode == 0, (proc.stdout + proc.stderr)[-4000:]
    except subprocess.TimeoutExpired:
        return False, "compilation timed out after 30 seconds"


def program220_lean_build() -> dict:
    core_ok, core_output = _compile_lean(DIRICHLET_CORE)
    general_ok, general_output = _compile_lean(GENERAL_LEAN)
    core = {
        "build_id": "FIN-P220-LEAN-DIRICHLET-001",
        "lean_binary": _direct_lean(),
        "dependency_free_source": DIRICHLET_CORE.name,
        "dependency_free_sha256": sha256(DIRICHLET_CORE),
        "dependency_free_compiled": core_ok,
        "dependency_free_output": core_output,
        "general_source": GENERAL_LEAN.name,
        "general_sha256": sha256(GENERAL_LEAN),
        "general_compiled": general_ok,
        "general_output": general_output,
        "scope": "C4 exact certificate versus general Mathlib theorem",
    }
    build_sha = write_canonical(LEAN_BUILD, core)

    fig, ax = plt.subplots(figsize=(8.5, 4.6), constrained_layout=True)
    labels = ["local compiler", "C4 exact core", "general Mathlib library"]
    values = [_direct_lean() is not None, core_ok, general_ok]
    ax.bar(labels, [int(v) for v in values], color=["#19733A" if v else "#A61B1B" for v in values])
    ax.set_ylim(0, 1.1)
    ax.set_yticks([0, 1], ["no", "yes"])
    ax.set_title("Program 220: machine checking advances, general dependency remains")
    fig.savefig(FIG / "program220_lean_build.png", dpi=190)
    plt.close(fig)

    return {
        "constructed_object": {
            "name": "ScopedLeanDirichletBuild",
            "path": LEAN_BUILD.name,
            "sha256": build_sha,
        },
        "dependency_free_source": DIRICHLET_CORE.name,
        "dependency_free_source_sha256": sha256(DIRICHLET_CORE),
        "dependency_free_machine_compiled": core_ok,
        "dependency_free_output_tail": core_output,
        "machine_checked_statements": [
            "constant edge-energy null mode",
            "constant Laplacian-pairing null mode",
            "archived equality of independent C4 formulas",
            "archived energy value 30",
            "archived nonnegativity",
        ],
        "general_mathlib_library_compiled": general_ok,
        "general_output_tail": general_output,
        "status": "Concrete C4 formal certificate compiled; general library still dependency-blocked",
        "confidence": "Machine checked in stated finite scope",
        "claim_boundary": (
            "The successful core compilation is not a proof of the general "
            "finite weighted-graph theorem, semigroup positivity, or physical adequacy."
        ),
    }


def _binary_entropy(probability: float) -> float:
    if probability <= 0.0 or probability >= 1.0:
        return 0.0
    return -probability * math.log2(probability) - (1.0 - probability) * math.log2(1.0 - probability)


def _correlated_reference_mutual_information(r: float, target: float) -> float:
    q = target / r
    px = np.asarray([(1.0 + r) / 2.0, (1.0 - r) / 2.0])
    conditional = np.asarray(
        [[(1.0 + q) / 2.0, (1.0 - q) / 2.0],
         [(1.0 - q) / 2.0, (1.0 + q) / 2.0]]
    )
    joint = px[:, None] * conditional
    py = joint.sum(axis=0)
    mutual = 0.0
    for i in range(2):
        for j in range(2):
            if joint[i, j] > 0:
                mutual += joint[i, j] * math.log2(joint[i, j] / (px[i] * py[j]))
    return mutual


def program221_reference_cost() -> dict:
    target = 0.6
    radii = np.linspace(target, 1.0, 401)
    mutual = np.asarray(
        [_correlated_reference_mutual_information(float(r), target) for r in radii]
    )
    input_fidelity = 1.0 - radii**2
    ideal_product_fidelity = input_fidelity * (1.0 - target**2)
    violation = input_fidelity - ideal_product_fidelity

    example_r = 0.8
    example_q = target / example_r
    example_mi = _correlated_reference_mutual_information(example_r, target)
    target_marginal = example_r * example_q
    reference_marginal = example_r

    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.7), constrained_layout=True)
    axes[0].plot(radii, violation, color="#A61B1B")
    axes[0].axhline(0, color="black", ls="--")
    axes[0].set_xlabel("reference asymmetry r")
    axes[0].set_ylabel("fidelity monotonicity violation")
    axes[0].set_title("uncorrelated exact return is impossible for r<1")
    axes[1].plot(radii, mutual, color="#1F5A99")
    axes[1].set_xlabel("reference asymmetry r")
    axes[1].set_ylabel("required output correlation (bits)")
    axes[1].set_title("marginal return is possible by storing correlations")
    fig.suptitle("Program 221: imperfect symmetry references have a correlation cost")
    fig.savefig(FIG / "program221_reference_cost.png", dpi=190)
    plt.close(fig)

    return {
        "constructed_object": {
            "name": "CorrelatedReferenceLedger",
            "uncorrelated_rule": "return catalyst exactly and in product with the target",
            "correlated_rule": "return the catalyst marginal exactly while recording target-reference correlation",
        },
        "no_broadcast_theorem": (
            "When no other asymmetric input is consumed, kappa_r and its "
            "reflected orbit have squared fidelity 1-r^2. An ideal product "
            "output with target transverse asymmetry t has fidelity "
            "(1-r^2)(1-t^2). CPTP fidelity monotonicity forbids this for "
            "0<=r<1 and t>0."
        ),
        "target_transverse_asymmetry": target,
        "audited_reference_radius_interval": [target, 1.0],
        "positive_uncorrelated_violation_points": int(np.sum(violation[:-1] > 0)),
        "perfect_reference_violation": float(violation[-1]),
        "correlated_example": {
            "reference_radius": example_r,
            "conditional_preparation_radius": example_q,
            "returned_reference_marginal": reference_marginal,
            "target_marginal": target_marginal,
            "mutual_information_bits": example_mi,
        },
        "strict_selector_source_exported": False,
        "status": "Exact imperfect-reference no-broadcast theorem and correlated workaround constructed",
        "confidence": "Proven",
        "claim_boundary": (
            "The no-broadcast theorem assumes no second asymmetric source is "
            "consumed; it is not a full classification of source-plus-catalyst "
            "conversion. Correlation bookkeeping quantifies a supplied "
            "reference and neither creates its sign nor discharges QW-2191."
        ),
    }


def _optimal_thinning(n: int, p_lower: float, delta_sequence: float) -> dict:
    best = None
    for lag in range(1, n + 1):
        count = (n - 1) // lag + 1
        coupling = (count - 1) * (1.0 - p_lower) ** lag
        dkw_budget = delta_sequence - coupling
        if dkw_budget <= 0:
            continue
        epsilon = math.sqrt(math.log(2.0 / dkw_budget) / (2.0 * count))
        row = {
            "lag": lag,
            "count": count,
            "coupling_budget": coupling,
            "dkw_budget": dkw_budget,
            "epsilon": epsilon,
        }
        if best is None or epsilon < best["epsilon"]:
            best = row
    if best is None:
        raise RuntimeError("no admissible thinning design")
    return best


def program222_robust_mixing_design() -> dict:
    rng = rng_for(222)
    p_true = 0.06
    calibration_size = 10000
    calibration_delta = 0.01
    observed_refreshes = int(rng.binomial(calibration_size, p_true))
    p_lower = float(
        beta_distribution.ppf(
            calibration_delta,
            observed_refreshes,
            calibration_size - observed_refreshes + 1,
        )
    )
    total_delta = 0.05
    delta_sequence = (total_delta - calibration_delta) / 2.0
    robust = _optimal_thinning(10000, p_lower, delta_sequence)
    oracle = _optimal_thinning(10000, p_true, total_delta / 2.0)
    inherited_equal_split_epsilon = json.loads(
        PREVIOUS_RESULTS.read_text(encoding="utf-8")
    )["programs"]["209"]["valid_epsilon"]

    calibration_reps = 3000
    lower_covers = 0
    lower_values = []
    for _ in range(calibration_reps):
        k = int(rng.binomial(calibration_size, p_true))
        lower = (
            0.0
            if k == 0
            else float(
                beta_distribution.ppf(
                    calibration_delta,
                    k,
                    calibration_size - k + 1,
                )
            )
        )
        lower_values.append(lower)
        lower_covers += lower <= p_true

    p_grid = np.linspace(0.035, 0.065, 61)
    eps_grid = [
        _optimal_thinning(10000, float(p), delta_sequence)["epsilon"]
        for p in p_grid
    ]
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.7), constrained_layout=True)
    axes[0].hist(lower_values, bins=40, color="#1F5A99", alpha=0.8)
    axes[0].axvline(p_true, color="black", ls="--", label="true refresh rate")
    axes[0].set_xlabel("one-sided calibrated lower bound")
    axes[0].legend(fontsize=8)
    axes[1].plot(p_grid, eps_grid, color="#6A3D9A")
    axes[1].scatter([p_lower], [robust["epsilon"]], color="#D55E00", label="realized robust design")
    axes[1].set_xlabel("certified p lower bound")
    axes[1].set_ylabel("optimized DKW radius")
    axes[1].legend(fontsize=8)
    fig.suptitle("Program 222: robust acquisition design under refresh-rate uncertainty")
    fig.savefig(FIG / "program222_robust_mixing_design.png", dpi=190)
    plt.close(fig)

    return {
        "constructed_object": {
            "name": "RobustMixingDesignOptimizer",
            "objective": "minimize the DKW radius over integer thinning lag after paying the exact Berbee coupling budget",
        },
        "record_length": 10000,
        "true_refresh_rate_for_challenge": p_true,
        "calibration_size": calibration_size,
        "observed_refreshes": observed_refreshes,
        "one_sided_calibrated_p_lower": p_lower,
        "calibration_failure_budget": calibration_delta,
        "per_sequence_failure_budget": delta_sequence,
        "robust_design": robust,
        "oracle_design": oracle,
        "inherited_equal_split_epsilon": inherited_equal_split_epsilon,
        "calibration_replicates": calibration_reps,
        "empirical_lower_bound_coverage": lower_covers / calibration_reps,
        "status": "Thinning and failure-budget allocation optimized with calibration uncertainty",
        "confidence": "Proven bound and finite optimization; strong calibration simulation",
        "claim_boundary": (
            "The design is valid only when the separate calibration experiment "
            "certifies the geometric mixing envelope. It does not infer mixing "
            "from the target record alone."
        ),
    }


def _harmonic(number: int) -> float:
    return sum(1.0 / k for k in range(1, number + 1))


def _eprocess_run(
    rng: np.random.Generator,
    streams: int,
    calibration: int,
    horizon: int,
    shift: float,
    reuse_fixed: bool,
) -> tuple[float, float | None, list[float]]:
    threshold = 20.0
    crossed = 0
    crossing_times = []
    example = None
    for stream in range(streams):
        pool = list(rng.normal(size=calibration))
        fixed = np.asarray(pool)
        log_e = 0.0
        trajectory = []
        first_cross = None
        for time in range(1, horizon + 1):
            score = float(rng.normal(loc=shift))
            if reuse_fixed:
                rank = 1 + int(np.sum(fixed >= score))
                size = calibration + 1
            else:
                rank = 1 + sum(value >= score for value in pool)
                size = len(pool) + 1
            evalue = size / (_harmonic(size) * rank)
            log_e += math.log(evalue)
            trajectory.append(math.exp(min(log_e, 50.0)))
            if not reuse_fixed:
                pool.append(score)
            if first_cross is None and log_e >= math.log(threshold):
                first_cross = time
        if first_cross is not None:
            crossed += 1
            crossing_times.append(first_cross)
        if stream == 0:
            example = trajectory
    return (
        crossed / streams,
        None if not crossing_times else float(np.mean(crossing_times)),
        example or [],
    )


def program223_reusable_eprocess() -> dict:
    rng = rng_for(223)
    calibration = 30
    horizon = 50
    streams = 1200
    null_rate, null_cross, null_example = _eprocess_run(
        rng, streams, calibration, horizon, 0.0, False
    )
    alternative_rate, alternative_cross, alternative_example = _eprocess_run(
        rng, streams, calibration, horizon, 2.5, False
    )
    naive_rate, naive_cross, naive_example = _eprocess_run(
        rng, streams, calibration, horizon, 0.0, True
    )
    # Reusing one fixed calibration block fails the conditional e-value
    # requirement on a positive-probability set of calibration samples.
    conditional_trials = 5000
    calibration_uniforms = np.sort(
        norm.cdf(rng.normal(size=(conditional_trials, calibration))), axis=1
    )
    gaps = np.diff(
        np.column_stack(
            [
                np.zeros(conditional_trials),
                calibration_uniforms,
                np.ones(conditional_trials),
            ]
        ),
        axis=1,
    )
    fixed_evalues = (calibration + 1) / (
        _harmonic(calibration + 1) * np.arange(1, calibration + 2)
    )
    conditional_mean_e = gaps[:, ::-1] @ fixed_evalues

    fig, ax = plt.subplots(figsize=(9.2, 5.0), constrained_layout=True)
    ax.semilogy(range(1, horizon + 1), null_example, label="online-rank iid null")
    ax.semilogy(range(1, horizon + 1), alternative_example, label="online-rank shifted alternative")
    ax.semilogy(range(1, horizon + 1), naive_example, label="naive fixed-calibration reuse")
    ax.axhline(20.0, color="black", ls="--", label="Ville threshold")
    ax.set_xlabel("sequential observation")
    ax.set_ylabel("example product process")
    ax.legend(fontsize=8)
    ax.set_title("Program 223: insertion ranks permit exact calibration reuse")
    ax.grid(True, which="both", alpha=0.25)
    fig.savefig(FIG / "program223_reusable_eprocess.png", dpi=190)
    plt.close(fig)

    return {
        "constructed_object": {
            "name": "SequentialInsertionRankEProcess",
            "rank": "descending insertion rank of the new score among every score observed so far",
            "e_value": "N/(H_N*r_N)",
        },
        "theorem": (
            "For iid continuous scores, sequential insertion ranks are "
            "independent and uniform on {1,...,N}. Therefore the product of "
            "N/(H_N r_N) is a nonnegative mean-one martingale in the "
            "filtration generated by the insertion ranks."
        ),
        "initial_calibration_records": calibration,
        "horizon": horizon,
        "streams": streams,
        "ville_alpha": 0.05,
        "online_rank_null_rejection_rate": null_rate,
        "online_rank_null_mean_crossing_time": null_cross,
        "online_rank_shifted_rejection_rate": alternative_rate,
        "online_rank_shifted_mean_crossing_time": alternative_cross,
        "naive_fixed_calibration_null_rejection_rate": naive_rate,
        "naive_fixed_calibration_mean_crossing_time": naive_cross,
        "naive_conditional_calibration_trials": conditional_trials,
        "naive_fraction_with_conditional_mean_e_above_one": float(
            np.mean(conditional_mean_e > 1.0)
        ),
        "naive_maximum_conditional_mean_e": float(np.max(conditional_mean_e)),
        "status": "Anytime-valid reusable-calibration e-process constructed and challenged",
        "confidence": "Proven rank-filtration martingale; strong synthetic power and conditional-failure evidence",
        "claim_boundary": (
            "The theorem assumes sequential iid exchangeability, continuous "
            "scores, and stopping rules measurable in the insertion-rank "
            "filtration. Arbitrary score-dependent stopping, adaptive drift, "
            "or dependent records require a different law."
        ),
    }


VISIBILITY_TO_PARAMETERS = np.asarray(
    [[2.0, -0.5, -0.5], [2.0, -1.0, -1.0], [0.0, -0.5, 0.5]]
)


def _visibility_counts(
    rng: np.random.Generator, visibilities: np.ndarray, shots: int
) -> tuple[np.ndarray, np.ndarray]:
    zero = np.asarray(
        [rng.binomial(shots, (1.0 + value) / 2.0) for value in visibilities]
    )
    pi = np.asarray(
        [rng.binomial(shots, (1.0 - value) / 2.0) for value in visibilities]
    )
    return zero, pi


def _wald_estimate(
    zero: np.ndarray, pi: np.ndarray, shots: int
) -> tuple[np.ndarray, np.ndarray]:
    visibility = np.clip(zero / shots - pi / shots, 1e-8, 1 - 1e-8)
    p0 = np.clip(zero / shots, 1e-6, 1 - 1e-6)
    pp = np.clip(pi / shots, 1e-6, 1 - 1e-6)
    variance_visibility = p0 * (1 - p0) / shots + pp * (1 - pp) / shots
    log_variance = variance_visibility / visibility**2
    theta = VISIBILITY_TO_PARAMETERS @ np.log(visibility)
    covariance = (
        VISIBILITY_TO_PARAMETERS
        @ np.diag(log_variance)
        @ VISIBILITY_TO_PARAMETERS.T
    )
    return theta, covariance


def program224_likelihood_region() -> dict:
    rng = rng_for(224)
    blur, variance, covariance = 0.84, 0.45, 0.20
    truth = np.asarray([math.log(blur), variance, covariance])
    visibilities = np.asarray(
        [
            blur * math.exp(-variance / 2.0),
            blur * math.exp(-variance - covariance),
            blur * math.exp(-variance + covariance),
        ]
    )
    shots = 1000
    replicates = 1200
    radius2 = float(chi2.ppf(0.95, 3))
    radius = math.sqrt(radius2)
    covered = 0
    memory = 0
    widths = []
    null_control_reject = 0
    alternative_control_reject = 0
    z_cut = float(norm.ppf(0.975))
    h_control = np.asarray([1.0, -2.0, 0.0])
    null_v4 = blur * math.exp(-2.0 * variance)
    alternative_v4 = null_v4 * math.exp(-0.20)
    for _ in range(replicates):
        zero, pi = _visibility_counts(rng, visibilities, shots)
        estimate, covariance_hat = _wald_estimate(zero, pi, shots)
        difference = estimate - truth
        try:
            statistic = float(difference @ np.linalg.solve(covariance_hat, difference))
        except np.linalg.LinAlgError:
            statistic = math.inf
        covered += statistic <= radius2
        standard = np.sqrt(np.maximum(np.diag(covariance_hat), 0.0))
        widths.append(2.0 * radius * standard)
        memory += estimate[2] - radius * standard[2] > 0

        for value, label in [(null_v4, "null"), (alternative_v4, "alternative")]:
            k0, kp = _visibility_counts(rng, np.asarray([value]), shots)
            estimate_v = float(np.clip(k0[0] / shots - kp[0] / shots, 1e-8, 1 - 1e-8))
            p0 = np.clip(k0[0] / shots, 1e-6, 1 - 1e-6)
            pp = np.clip(kp[0] / shots, 1e-6, 1 - 1e-6)
            var_v = p0 * (1 - p0) / shots + pp * (1 - pp) / shots
            residual = math.log(estimate_v) - float(h_control @ estimate)
            residual_variance = var_v / estimate_v**2 + float(
                h_control @ covariance_hat @ h_control
            )
            rejected = abs(residual) / math.sqrt(residual_variance) > z_cut
            if label == "null":
                null_control_reject += rejected
            else:
                alternative_control_reject += rejected

    mean_widths = np.mean(widths, axis=0)
    inherited = np.asarray(
        list(
            json.loads(PREVIOUS_RESULTS.read_text(encoding="utf-8"))["programs"]["211"][
                "mean_parameter_interval_widths"
            ].values()
        )
    )
    fig, ax = plt.subplots(figsize=(9.0, 5.0), constrained_layout=True)
    x = np.arange(3)
    ax.bar(x - 0.18, inherited, 0.36, label="exact conservative P211", color="#9E9E9E")
    ax.bar(x + 0.18, mean_widths, 0.36, label="likelihood-quadratic P224", color="#1F5A99")
    ax.set_xticks(x, [r"$\log b$", r"$v$", r"$c$"])
    ax.set_ylabel("mean simultaneous marginal width")
    ax.legend(fontsize=8)
    ax.set_title("Program 224: efficiency gain with an explicit asymptotic coverage boundary")
    fig.savefig(FIG / "program224_likelihood_region.png", dpi=190)
    plt.close(fig)

    return {
        "constructed_object": {
            "name": "LikelihoodQuadraticProcessRegion",
            "region": "delta-method likelihood ellipsoid with chi-square(3) radius",
            "misspecification_control": "independent fourth contrast log(V4)=log(b)-2v",
        },
        "shots_per_phase_context": shots,
        "replicates": replicates,
        "simultaneous_empirical_coverage": covered / replicates,
        "memory_detection_power": memory / replicates,
        "mean_parameter_interval_widths": {
            "log_blur": float(mean_widths[0]),
            "variance": float(mean_widths[1]),
            "covariance": float(mean_widths[2]),
        },
        "inherited_exact_region_mean_widths": {
            "log_blur": float(inherited[0]),
            "variance": float(inherited[1]),
            "covariance": float(inherited[2]),
        },
        "independent_control_null_rejection_rate": null_control_reject / replicates,
        "independent_control_misspecified_rejection_rate": alternative_control_reject / replicates,
        "status": "Efficient likelihood-quadratic region and independent misspecification control executed",
        "confidence": "Strong simulation evidence; asymptotic region, not exact finite-sample coverage",
        "claim_boundary": (
            "Unlike Program 211, the ellipsoid is not an exact finite-sample "
            "confidence set. It is retained only with its empirical coverage "
            "and independent model-control diagnostic."
        ),
    }


def program225_higher_moment() -> dict:
    c1 = Fraction(4, 5)
    c2 = Fraction(1, 2)
    lower = Fraction(-11, 180)
    upper = Fraction(11, 20)
    # Endpoint measures in x=cos(theta).
    lower_mass_minus_one = Fraction(11, 335)
    lower_mass_r = Fraction(324, 335)
    lower_r = Fraction(31, 36)
    upper_mass_one = Fraction(11, 15)
    upper_mass_quarter = Fraction(4, 15)
    lower_moments = [
        lower_mass_minus_one * (-1) ** k + lower_mass_r * lower_r**k
        for k in range(1, 4)
    ]
    upper_moments = [
        upper_mass_one + upper_mass_quarter * Fraction(1, 4) ** k
        for k in range(1, 4)
    ]
    lower_c3 = 4 * lower_moments[2] - 3 * lower_moments[0]
    upper_c3 = 4 * upper_moments[2] - 3 * upper_moments[0]
    grid = np.linspace(float(lower) - 0.1, float(upper) + 0.1, 801)
    determinants = -(20 * grid - 11) * (180 * grid + 11) / 10000.0

    fig, ax = plt.subplots(figsize=(8.9, 5.0), constrained_layout=True)
    ax.plot(grid, determinants, color="#1F5A99")
    ax.axhline(0, color="black", ls="--")
    ax.axvspan(float(lower), float(upper), color="#B7D7F0", alpha=0.65)
    ax.set_xlabel(r"unmeasured third characteristic $c_3$")
    ax.set_ylabel(r"$\det T_3$")
    ax.set_title("Program 225: sharp third-moment interval from c1=4/5 and c2=1/2")
    ax.grid(True, alpha=0.25)
    fig.savefig(FIG / "program225_higher_moment.png", dpi=190)
    plt.close(fig)

    return {
        "constructed_object": {
            "name": "ThirdMomentToeplitzInterval",
            "matrix": "T3=Toeplitz(1,c1,c2,c3)",
            "determinant_exact": "-(20*c3-11)*(180*c3+11)/10000",
        },
        "given_c1_exact": "4/5",
        "given_c2_exact": "1/2",
        "sharp_c3_lower_exact": "-11/180",
        "sharp_c3_upper_exact": "11/20",
        "sharp_c3_lower": float(lower),
        "sharp_c3_upper": float(upper),
        "lower_endpoint_measure": {
            "mass_at_pi": "11/335",
            "mass_split_at_plus_minus_acos_31_over_36": "162/335 each",
        },
        "upper_endpoint_measure": {
            "mass_at_zero": "11/15",
            "mass_split_at_plus_minus_acos_1_over_4": "2/15 each",
        },
        "lower_witness_moments_x": [str(value) for value in lower_moments],
        "upper_witness_moments_x": [str(value) for value in upper_moments],
        "lower_witness_c3_exact": str(lower_c3),
        "upper_witness_c3_exact": str(upper_c3),
        "status": "Sharp order-three trigonometric moment problem solved",
        "confidence": "Proven",
        "claim_boundary": (
            "The interval constrains one further operational environment "
            "moment. It does not identify a unique microscopic phase law."
        ),
    }


def program226_phase_formalization() -> dict:
    finite_ok, finite_output = _compile_lean(PHASE_CORE)
    core = {
        "formalization_id": "FIN-P226-PHASE-TORSION-001",
        "source": PHASE_CORE.name,
        "source_sha256": sha256(PHASE_CORE),
        "machine_compiled": finite_ok,
        "output": finite_output,
        "machine_checked": [
            "complete order-eight residue orbit",
            "every image remains in an order-eight slot",
            "gcd(743,4000)=1",
        ],
        "not_machine_checked": [
            "automatic continuity of Borel homomorphisms",
            "classification of continuous U(1) endomorphisms",
            "transcendence-based non-torsion theorem",
        ],
    }
    formal_sha = write_canonical(PHASE_FORMAL, core)

    fig, ax = plt.subplots(figsize=(8.8, 4.7), constrained_layout=True)
    labels = ["finite orbit", "reduced fraction", "automatic continuity", "non-torsion theorem"]
    values = [finite_ok, finite_ok, False, False]
    ax.bar(labels, [int(value) for value in values], color=["#19733A" if value else "#A61B1B" for value in values])
    ax.set_ylim(0, 1.1)
    ax.set_yticks([0, 1], ["not checked", "machine checked"])
    ax.tick_params(axis="x", rotation=17)
    ax.set_title("Program 226: finite phase core formalized; analytic hierarchy remains")
    fig.savefig(FIG / "program226_phase_formalization.png", dpi=190)
    plt.close(fig)

    return {
        "constructed_object": {
            "name": "ScopedPhaseFormalizationRecord",
            "path": PHASE_FORMAL.name,
            "sha256": formal_sha,
        },
        "lean_source": PHASE_CORE.name,
        "lean_source_sha256": sha256(PHASE_CORE),
        "finite_torsion_core_machine_compiled": finite_ok,
        "finite_output_tail": finite_output,
        "analytic_automatic_continuity_machine_formalized": False,
        "transcendence_non_torsion_machine_formalized": False,
        "strict_phase_source_exported": False,
        "status": "Finite torsion layer machine checked; analytic no-go remains paper mathematics",
        "confidence": "Machine checked in finite scope; analytic theorem inherited",
        "claim_boundary": (
            "Compilation of the finite orbit is not formalization of automatic "
            "continuity or transcendence of pi and does not source the strict phase."
        ),
    }


def _heat_reconstruct_signature(
    transition: np.ndarray,
    duration: float,
    rng: np.random.Generator,
    shots: int,
) -> tuple[np.ndarray | None, float]:
    raw = np.column_stack(
        [
            rng.multinomial(shots, np.clip(transition[:, column], 0, 1)) / shots
            for column in range(transition.shape[1])
        ]
    )
    estimate = (raw + raw.T) / 2.0
    estimate += np.diag(1.0 - estimate.sum(axis=1))
    eigenvalues = np.linalg.eigvalsh(estimate)[::-1]
    if eigenvalues[-1] <= 0 or eigenvalues[0] <= 0:
        return None, float(np.linalg.norm(estimate - transition, ord=2))
    rates = -np.log(np.clip(eigenvalues / eigenvalues[0], 1e-15, 1.0)) / duration
    positive = np.sort(rates[1:])
    return positive / positive[-1], float(np.linalg.norm(estimate - transition, ord=2))


def program227_spectral_tomography() -> dict:
    rng = rng_for(227)
    generator, _, _ = previous._kernel_generator(
        12,
        previous.OMEGA_STRICT,
        previous.PHI_STRICT,
        previous.BETA_STRICT,
        previous.ETA_STRICT,
    )
    eigenvalues = np.linalg.eigvalsh(generator)
    duration = 1.0 / eigenvalues[-1]
    transition = expm(-duration * generator)
    target = eigenvalues[1:] / eigenvalues[-1]
    exact_from_transition = -np.log(
        np.sort(np.linalg.eigvalsh(transition))[::-1]
    ) / duration
    exact_signature = np.sort(exact_from_transition[1:]) / max(exact_from_transition[1:])
    exact_residual = float(np.max(np.abs(exact_signature - target)))

    shot_levels = [2000, 10000, 50000]
    replicates = 240
    rows = []
    distributions = {}
    for shots in shot_levels:
        distances = []
        operator_errors = []
        spd_failures = 0
        for _ in range(replicates):
            signature, error = _heat_reconstruct_signature(
                transition, duration, rng, shots
            )
            operator_errors.append(error)
            if signature is None:
                spd_failures += 1
                continue
            distances.append(float(np.max(np.abs(signature - target))))
        distributions[shots] = distances
        rows.append(
            {
                "shots_per_preparation": shots,
                "replicates": replicates,
                "spd_failures": spd_failures,
                "mean_fingerprint_distance": float(np.mean(distances)),
                "p95_fingerprint_distance": float(np.quantile(distances, 0.95)),
                "program214_acceptance_rate": float(np.mean(np.asarray(distances) <= 0.02)),
                "mean_transition_operator_error": float(np.mean(operator_errors)),
            }
        )

    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.8), constrained_layout=True)
    for shots in shot_levels:
        axes[0].hist(
            distributions[shots],
            bins=35,
            alpha=0.45,
            label=f"{shots} shots/input",
        )
    axes[0].axvline(0.02, color="black", ls="--", label="P214 threshold")
    axes[0].set_xlabel("projective fingerprint distance")
    axes[0].legend(fontsize=7)
    axes[1].semilogx(
        shot_levels,
        [row["program214_acceptance_rate"] for row in rows],
        marker="o",
        color="#19733A",
    )
    axes[1].set_ylim(-0.02, 1.02)
    axes[1].set_xlabel("shots per basis preparation")
    axes[1].set_ylabel("strict fingerprint acceptance rate")
    fig.suptitle("Program 227: an explicit heat-kernel instrument reconstructs spectral ratios")
    fig.savefig(FIG / "program227_spectral_tomography.png", dpi=190)
    plt.close(fig)

    return {
        "constructed_object": {
            "name": "HeatKernelTomographyInstrument",
            "preparations": "all twelve vertex basis states",
            "dynamics": "one calibrated positive duration of the heat semigroup",
            "measurement": "twelve-outcome vertex POVM",
            "record": "one multinomial transition column per preparation",
        },
        "identifiability_theorem": (
            "If P_tau=exp(-tau*A) is known, symmetric positive definite, and "
            "tau>0, then A=-(1/tau)log(P_tau) is unique. Without a dimensional "
            "clock, eigenvalue ratios remain uniquely recoverable from "
            "log(mu_k)/log(mu_min), modulo node permutation."
        ),
        "finite_error_theorem": (
            "If ||P_hat-P||_2<=epsilon<lambda_min(P), Weyl's theorem and the "
            "Lipschitz bound for log on [lambda_min(P)-epsilon,1+epsilon] "
            "give explicit eigenrate and projective-signature error bounds."
        ),
        "dimensionless_duration": duration,
        "transition_minimum_eigenvalue": float(np.min(np.linalg.eigvalsh(transition))),
        "exact_reconstruction_signature_residual": exact_residual,
        "acceptance_threshold": 0.02,
        "shot_rows": rows,
        "external_data_used": False,
        "status": "Instrument-to-spectrum theorem proved and finite-shot challenge executed",
        "confidence": "Proven identifiability; strong finite-shot simulation evidence",
        "claim_boundary": (
            "The instrument is a conditional operational construction. It "
            "does not provide a dimensional clock, laboratory implementation, "
            "or independent physical record."
        ),
    }


def program228_external_custody(previous_results: dict) -> dict:
    inherited = previous_results["programs"]["215"]["candidate_rows"]
    core = {
        "request_id": "FIN-P228-INDEPENDENT-CUSTODY-11-OF-11-V2",
        "accepted_resource_types": ["ordered process record", "double-slit event record"],
        "required_fields": [
            "independent provider identity and license",
            "immutable raw-event hash",
            "ordered event or detection coordinates",
            "preparation or slit-configuration labels",
            "intervention or shutter-control labels",
            "clock and dimensional calibration",
            "detector geometry, efficiency, dark-count and blur calibration",
            "environment and background record",
            "negative-control record",
            "held-out target committed before analysis",
            "provider-registrar-analyst custody separation",
        ],
        "double_slit_specialization": {
            "minimum_configurations": ["both open", "left only", "right only", "both closed"],
            "required_event_fields": ["detector coordinate", "timestamp", "configuration", "run id"],
            "forbidden_substitute": "rendered interference image without event-level data",
        },
        "admission": "all eleven fields pass before unblinding the held-out target",
        "firecrawl_or_web_used": False,
    }
    request_sha = write_canonical(EXTERNAL_REQUEST, core)
    rows = [dict(row) for row in inherited]
    admitted = [row["bundle"] for row in rows if row["admitted"]]

    fig, ax = plt.subplots(figsize=(9.0, 4.8), constrained_layout=True)
    ax.bar(
        [row["bundle"].replace("_", "\n") for row in rows],
        [row["passed_fields"] for row in rows],
        color="#A61B1B",
    )
    ax.axhline(11, color="black", ls="--", label="independent admission")
    ax.set_ylim(0, 11.7)
    ax.set_ylabel("passed custody/intake fields")
    ax.tick_params(axis="x", labelsize=7)
    ax.legend(fontsize=8)
    ax.set_title("Program 228: no independent process or double-slit bundle is admitted")
    fig.savefig(FIG / "program228_external_custody.png", dpi=190)
    plt.close(fig)

    return {
        "constructed_object": {
            "name": "IndependentCustodyBundleV2",
            "path": EXTERNAL_REQUEST.name,
            "sha256": request_sha,
        },
        "candidate_rows": rows,
        "admitted_bundles": admitted,
        "external_11_of_11_bundle_found": bool(admitted),
        "double_slit_event_bundle_found": False,
        "web_or_firecrawl_used": False,
        "status": "Acquisition schema strengthened; no independent bundle supplied",
        "confidence": "Proven for current local intake scope",
        "claim_boundary": (
            "A schema and custody split are not experimental data. No rendered "
            "double-slit image, synthetic bundle, or same-analysis record is "
            "accepted as independent physical evidence."
        ),
    }


def program229_semigroup_prediction(program228: dict) -> dict:
    core = {
        "protocol_id": "FIN-P229-EXTERNAL-SEMIGROUP-HOLDOUT-001",
        "prerequisite": "Program 228 admits one independent 11-of-11 bundle",
        "calibration": "reconstruct P_tau from the committed calibration subset",
        "held_out_prediction": "P_2tau_predicted = P_tau @ P_tau",
        "primary_statistic": "maximum column total-variation distance to held-out P_2tau",
        "spectral_secondary": "Program-214 projective fingerprint distance",
        "negative_controls": [
            "time-order permutation",
            "preparation-label permutation",
            "nearest-neighbour C12",
            "both-closed/dark-count double-slit control when applicable",
        ],
        "execution_rule": "one execution after immutable bundle hash and analyst unblinding",
        "frozen_before_external_data": True,
    }
    prediction_sha = write_canonical(PREDICTION_LOCK, core)
    gate = program228["external_11_of_11_bundle_found"]

    fig, ax = plt.subplots(figsize=(8.4, 4.5), constrained_layout=True)
    ax.bar(
        ["external bundle", "P_tau reconstruction", "P_2tau holdout", "physical validation"],
        [int(gate), 0, 0, 0],
        color=["#19733A" if gate else "#A61B1B", "#A61B1B", "#A61B1B", "#A61B1B"],
    )
    ax.set_ylim(0, 1.1)
    ax.set_yticks([0, 1], ["not completed", "completed"])
    ax.tick_params(axis="x", rotation=15)
    ax.set_title("Program 229: external semigroup prediction remains locked")
    fig.savefig(FIG / "program229_prediction_lock.png", dpi=190)
    plt.close(fig)

    return {
        "constructed_object": {
            "name": "SemigroupHoldoutPredictionLock",
            "path": PREDICTION_LOCK.name,
            "sha256": prediction_sha,
        },
        "program228_gate_passed": gate,
        "external_prediction_executed": False,
        "external_physical_validation_claimed": False,
        "status": "Semigroup holdout protocol frozen but not externally executed",
        "confidence": "Proven gate status",
        "claim_boundary": (
            "The algebraic identity P_2tau=P_tau^2 holds inside a time-homogeneous "
            "semigroup. Testing it as physics requires the independent record "
            "that Program 228 did not supply."
        ),
    }


def main() -> None:
    FIG.mkdir(exist_ok=True)
    previous_results = json.loads(PREVIOUS_RESULTS.read_text(encoding="utf-8"))
    programs: dict[str, dict] = {}
    programs["217"] = program217_operational_central_state()
    programs["218"] = program218_affine_eta_no_go()
    programs["219"] = program219_arb_execution()
    programs["220"] = program220_lean_build()
    programs["221"] = program221_reference_cost()
    programs["222"] = program222_robust_mixing_design()
    programs["223"] = program223_reusable_eprocess()
    programs["224"] = program224_likelihood_region()
    programs["225"] = program225_higher_moment()
    programs["226"] = program226_phase_formalization()
    programs["227"] = program227_spectral_tomography()
    programs["228"] = program228_external_custody(previous_results)
    programs["229"] = program229_semigroup_prediction(programs["228"])

    result = {
        "metadata": {
            "title": "FIN Programs 217-229: Operational Sections, Reference Costs, and Spectral Tomography",
            "release": "10.20",
            "version": "1.0.0",
            "date": "2026-07-27",
            "creator": "Żuchowski, Krzysztof",
            "affiliation": "Independent Researcher — Fractal Information Theory Project",
            "orcid": "0009-0002-0909-3613",
            "language": "English",
            "license": "CC BY 4.0",
            "seed": SEED,
            "firecrawl_used": False,
            "external_web_used": False,
            "external_dataset_used": False,
        },
        "programs": programs,
        "global_verdict": {
            "programs_executed": list(range(217, 230)),
            "new_theoretical_objects": [
                programs[str(number)]["constructed_object"]["name"]
                for number in range(217, 230)
            ],
            "bistochastic_central_state_classified": True,
            "strict_central_state_source_found": False,
            "affine_natural_eta_source_found": False,
            "full_arb_certificate": False,
            "finite_lean_dirichlet_core_compiled": programs["220"]["dependency_free_machine_compiled"],
            "general_lean_library_compiled": programs["220"]["general_mathlib_library_compiled"],
            "imperfect_reference_no_broadcast_proved": True,
            "robust_mixing_design_completed": True,
            "reusable_calibration_eprocess_completed": True,
            "efficient_process_region_completed": True,
            "sharp_third_moment_interval_completed": True,
            "phase_finite_core_compiled": programs["226"]["finite_torsion_core_machine_compiled"],
            "instrument_to_spectrum_theorem_completed": True,
            "external_bundle_admitted": programs["228"]["external_11_of_11_bundle_found"],
            "external_prediction_executed": False,
            "QW_2191_discharged": False,
            "strict_selector_exported": False,
            "canonical_physical_unit_exported": False,
            "legacy_strict_bridge_completed": False,
            "legacy_role_transfer_started": False,
            "L_total_or_ToE_claimed": False,
            "external_physical_validation_claimed": False,
        },
        "recommended_next_programs": {
            "230": "experimental preparation test of bistochastic versus frequency-selected central states",
            "231": "eta-flow classification with one genuinely new strict-derived non-affine datum or a broader no-go",
            "232": "authorized pinned Arb/FLINT execution with immutable container digest",
            "233": "pin Mathlib and complete the general finite Dirichlet and heat-semigroup formalization",
            "234": "multi-use finite-reference degradation, recovery, and correlation accumulation theorem",
            "235": "jointly optimal mixing calibration and target-record allocation under a fixed experimental budget",
            "236": "exchangeability-robust e-process under covariate shift and weak dependence",
            "237": "exact or anytime-valid finite-shot likelihood region for process parameters",
            "238": "order-four and order-five Toeplitz moment intervals with sparse extremal environment recovery",
            "239": "Mathlib formalization of circle endomorphisms, automatic continuity, and the non-torsion phase theorem",
            "240": "optimal heat-kernel tomography times, preparations, and nonasymptotic spectral error bounds",
            "241": "independent blind-custody acquisition of an event-level double-slit or process bundle",
            "242": "single registered execution of the semigroup holdout and projective fingerprint",
        },
    }
    RESULTS.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(
        json.dumps(
            {
                "programs_executed": result["global_verdict"]["programs_executed"],
                "new_theoretical_objects": result["global_verdict"]["new_theoretical_objects"],
                "strict_central_state_source_found": False,
                "affine_natural_eta_source_found": False,
                "full_arb_certificate": False,
                "finite_lean_dirichlet_core_compiled": result["global_verdict"]["finite_lean_dirichlet_core_compiled"],
                "general_lean_library_compiled": result["global_verdict"]["general_lean_library_compiled"],
                "imperfect_reference_no_broadcast_proved": True,
                "instrument_to_spectrum_theorem_completed": True,
                "external_bundle_admitted": result["global_verdict"]["external_bundle_admitted"],
                "external_prediction_executed": False,
                "QW_2191_discharged": False,
                "legacy_strict_bridge_completed": False,
                "L_total_or_ToE_claimed": False,
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
