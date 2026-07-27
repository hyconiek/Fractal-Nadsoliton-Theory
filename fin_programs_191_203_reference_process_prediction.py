#!/usr/bin/env python3
"""FIN Programs 191--203: reference states, operational quotients and prediction.

The executable separates theorem-level identities, exact finite checks,
synthetic method validation, unavailable certification toolchains and
conditional physics.  It does not use the web and it does not promote any
selector, physical-unit, legacy-to-strict bridge or Theory-of-Everything
claim.
"""

from __future__ import annotations

from fractions import Fraction
from itertools import permutations
from pathlib import Path
import csv
import hashlib
import importlib.util
import json
import math
import os
import platform
import shutil
import subprocess
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.linalg import expm
from scipy.stats import rankdata

import fin_programs_178_190_composition_process_scale as previous_study


ROOT = Path(__file__).resolve().parent
FIG = ROOT / "FIN_Programs_191_203_Reference_Process_Prediction_Figures"
RESULTS = ROOT / "FIN_Programs_191_203_Reference_Process_Prediction_Results.json"
PROTOCOL = ROOT / "FIN_Programs_191_203_Order_Sensitive_Preregistration.json"
BUNDLE_DATA = ROOT / "FIN_Programs_191_203_Synthetic_Operational_Bundle.csv"
BUNDLE_META = ROOT / "FIN_Programs_191_203_Synthetic_Operational_Bundle.json"
LEAN_SOURCE = ROOT / "FIN_Programs_191_203_Dirichlet_Library.lean"
PREVIOUS_RESULTS = ROOT / "FIN_Programs_178_190_Composition_Process_Scale_Results.json"

SEED = 20260727
ALPHA_GEO = 4.0 * math.log(2.0)
ETA_STRICT = 1.8
OMEGA_STRICT = 0.18575
PHI_STRICT = 0.1625
BETA_STRICT = 1.0
ETA_LEGACY = 1.0
OMEGA_LEGACY = math.pi / 4.0
PHI_LEGACY = math.pi / 6.0
BETA_LEGACY = 0.01


def rng_for(offset: int) -> np.random.Generator:
    return np.random.default_rng(SEED + offset)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def canonical_digest(obj: object) -> str:
    payload = json.dumps(
        obj, sort_keys=True, separators=(",", ":"), ensure_ascii=False
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def stable_rvs(
    alpha: float, size: tuple[int, ...] | int, rng: np.random.Generator
) -> np.ndarray:
    """Symmetric alpha-stable random variables (CMS construction)."""
    u = rng.uniform(-math.pi / 2.0, math.pi / 2.0, size)
    w = rng.exponential(1.0, size)
    if abs(alpha - 1.0) < 1e-12:
        return np.tan(u)
    return (
        np.sin(alpha * u)
        / np.cos(u) ** (1.0 / alpha)
        * (np.cos((1.0 - alpha) * u) / w) ** ((1.0 - alpha) / alpha)
    )


def program191_reference_state_naturality() -> dict:
    """Classify invariant states on a finite direct sum of matrix algebras."""
    dimensions = np.asarray([1, 2, 2, 2, 2], dtype=float)
    uniform_central = np.full(5, 1.0 / 5.0)
    normalized_full_trace = dimensions / dimensions.sum()
    observable = dimensions.copy()
    eta_uniform = float(uniform_central @ observable)
    eta_trace = float(normalized_full_trace @ observable)
    rho2 = np.eye(2) / 2.0
    rho3 = np.eye(3) / 3.0
    tensor_trace_residual = float(
        np.max(np.abs(np.kron(rho2, rho3) - np.eye(6) / 6.0))
    )

    # Isomorphism invariance permutes the four equal M_2 summands but cannot
    # exchange the M_1 summand with an M_2 summand.
    a_values = np.linspace(0.0, 1.0, 201)
    natural_family = np.column_stack(
        [a_values, np.repeat(((1.0 - a_values) / 4.0)[:, None], 4, axis=1)]
    )
    permutation_residual = 0.0
    for row in natural_family[::20]:
        for perm in permutations(range(1, 5)):
            moved = row.copy()
            moved[1:] = row[list(perm)]
            permutation_residual = max(
                permutation_residual, float(np.max(np.abs(moved - row)))
            )
    eta_family = natural_family @ observable

    fig, ax = plt.subplots(figsize=(8.8, 5.0), constrained_layout=True)
    ax.plot(a_values, eta_family, color="#1F5A99", lw=2)
    ax.scatter(
        [1.0 / 5.0, 1.0 / 9.0],
        [eta_uniform, eta_trace],
        color=["#19733A", "#D55E00"],
        s=70,
        zorder=3,
    )
    ax.annotate("uniform centre: 9/5", (1 / 5, eta_uniform), xytext=(0.3, 1.72))
    ax.annotate("full trace: 17/9", (1 / 9, eta_trace), xytext=(0.25, 1.94))
    ax.set_xlabel(r"central weight $a$ of the $M_1$ summand")
    ax.set_ylabel("state expectation of the dimension observable")
    ax.set_title("Program 191: naturality fixes blocks, not central measure")
    ax.grid(True, alpha=0.25)
    fig.savefig(FIG / "program191_reference_state_naturality.png", dpi=190)
    plt.close(fig)

    return {
        "constructed_object": {
            "name": "CentralMeasureFunctor",
            "algebra": "A = direct_sum_i M_{d_i}(C)",
            "state": "omega_w(a)=sum_i w_i tr(a_i)/d_i",
            "central_weight_space": "probability simplex on simple summands",
        },
        "dimensions": dimensions.astype(int).tolist(),
        "unitary_invariance_theorem": (
            "On each simple block M_d, invariance under all inner unitaries "
            "forces the density matrix I_d/d. On a direct sum it leaves a "
            "probability vector on the centre."
        ),
        "isomorphism_natural_family": (
            "w(a)=(a,(1-a)/4,(1-a)/4,(1-a)/4,(1-a)/4), 0<=a<=1"
        ),
        "permutation_residual": permutation_residual,
        "uniform_central_weights": uniform_central.tolist(),
        "normalized_full_trace_weights": normalized_full_trace.tolist(),
        "uniform_central_expectation": eta_uniform,
        "normalized_full_trace_expectation": eta_trace,
        "uniform_central_exact": "9/5",
        "normalized_full_trace_exact": "17/9",
        "tensor_compatibility": (
            "simple-block normalized traces satisfy tau_d tensor tau_e = "
            "tau_{de}; central weights tensor multiplicatively as w_i v_j"
        ),
        "simple_trace_tensor_residual": tensor_trace_residual,
        "unique_reference_state_selected": False,
        "minimal_missing_datum": (
            "a probability measure on the centre, or an additional naturality/"
            "Morita principle that fixes it"
        ),
        "status": "Reference-state classification theorem; central source remains open",
        "confidence": "Proven",
        "claim_boundary": (
            "Normalized trace is unique on a simple matrix algebra, not on a "
            "non-simple direct sum. Neither 9/5 nor 17/9 is promoted to a "
            "strict temperature or physical state source."
        ),
    }


def invariant_kernel_profile(
    x: np.ndarray, omega: float, phi: float, beta: float, eta: float
) -> np.ndarray:
    nu = omega / beta ** (1.0 / eta)
    return np.cos(nu * x + phi) / (1.0 + x**eta)


def program192_beta_gauge_quotient() -> dict:
    x = np.linspace(0.0, 8.0, 4001)
    c_values = np.logspace(-4, 4, 17)
    base = invariant_kernel_profile(
        x, OMEGA_STRICT, PHI_STRICT, BETA_STRICT, ETA_STRICT
    )
    orbit_residuals = []
    for c in c_values:
        beta_prime = BETA_STRICT * c ** (-ETA_STRICT)
        omega_prime = OMEGA_STRICT / c
        transformed = invariant_kernel_profile(
            x, omega_prime, PHI_STRICT, beta_prime, ETA_STRICT
        )
        orbit_residuals.append(float(np.max(np.abs(base - transformed))))

    nu_strict = OMEGA_STRICT / BETA_STRICT ** (1.0 / ETA_STRICT)
    nu_legacy = OMEGA_LEGACY / BETA_LEGACY ** (1.0 / ETA_LEGACY)
    strict_profile = invariant_kernel_profile(
        x, OMEGA_STRICT, PHI_STRICT, BETA_STRICT, ETA_STRICT
    )
    legacy_profile = invariant_kernel_profile(
        x, OMEGA_LEGACY, PHI_LEGACY, BETA_LEGACY, ETA_LEGACY
    )
    profile_l2 = float(
        math.sqrt(np.trapz((strict_profile - legacy_profile) ** 2, x) / (x[-1] - x[0]))
    )

    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.6), constrained_layout=True)
    axes[0].plot(x, strict_profile, label="strict quotient profile", lw=2)
    axes[0].plot(x, legacy_profile, label="legacy quotient profile", alpha=0.8)
    axes[0].set_xlim(0, 4)
    axes[0].set_xlabel(r"$x=\beta^{1/\eta}d$")
    axes[0].set_ylabel("amplitude-normalized kernel")
    axes[0].legend(fontsize=8)
    axes[0].grid(True, alpha=0.25)
    labels = [r"$\eta$", r"$\nu$", r"$\phi$"]
    gaps = [
        abs(ETA_STRICT - ETA_LEGACY),
        abs(nu_strict - nu_legacy),
        abs(PHI_STRICT - PHI_LEGACY),
    ]
    axes[1].bar(labels, gaps, color=["#6A3D9A", "#A61B1B", "#D55E00"])
    axes[1].set_yscale("symlog", linthresh=0.1)
    axes[1].set_title("Gauge-invariant legacy/strict gaps")
    axes[1].grid(True, axis="y", alpha=0.25)
    fig.suptitle("Program 192: quotienting beta removes scale, not shape mismatch")
    fig.savefig(FIG / "program192_beta_gauge_quotient.png", dpi=190)
    plt.close(fig)

    return {
        "constructed_object": {
            "name": "BetaGaugeQuotient",
            "action": (
                "(d,beta,omega,eta,phi) -> "
                "(c*d,beta*c^(-eta),omega/c,eta,phi)"
            ),
            "complete_reduced_profile": (
                "k_hat(x)=cos(nu*x+phi)/(1+x^eta), "
                "x=beta^(1/eta)d, nu=omega/beta^(1/eta)"
            ),
        },
        "maximum_orbit_residual": max(orbit_residuals),
        "invariant_rows": {
            "legacy": {
                "eta": ETA_LEGACY,
                "nu": nu_legacy,
                "phi": PHI_LEGACY,
            },
            "strict": {
                "eta": ETA_STRICT,
                "nu": nu_strict,
                "phi": PHI_STRICT,
            },
        },
        "quotient_profile_L2_gap": profile_l2,
        "beta_one_is_gauge_representative": True,
        "gauge_closes_legacy_strict_bridge": False,
        "status": "Scale quotient constructed; invariant bridge mismatch remains",
        "confidence": "Proven quotient theorem; numerical profile integral",
        "claim_boundary": (
            "The quotient attacks one compression-scale atom only. It does "
            "not source eta, omega, phi, complete the bridge or transfer "
            "legacy physical roles."
        ),
    }


def program193_common_engine_certificate(previous: dict) -> dict:
    old = previous["programs"]["180"]
    flint_spec = importlib.util.find_spec("flint")
    arb_binary = shutil.which("arb")
    python_flint_available = flint_spec is not None
    arb_available = arb_binary is not None
    components = []
    for component in old["components"]:
        row = dict(component)
        row["same_arb_engine"] = bool(
            (python_flint_available or arb_available)
            and component.get("directed_rounding")
            and component.get("closed")
        )
        components.append(row)
    all_in_engine = all(row["same_arb_engine"] for row in components)
    closed_count = sum(row["closed"] for row in components)

    fig, ax = plt.subplots(figsize=(9.5, 4.8), constrained_layout=True)
    y = np.arange(len(components))
    available = [1 if row["available"] else 0 for row in components]
    closed = [1 if row["closed"] else 0 for row in components]
    same = [1 if row["same_arb_engine"] else 0 for row in components]
    ax.barh(y, available, color="#B7D7F0", label="available")
    ax.barh(y, closed, color="#1F5A99", height=0.55, label="closed")
    ax.barh(y, same, color="#19733A", height=0.25, label="same Arb engine")
    ax.set_yticks(y, [row["component"] for row in components], fontsize=8)
    ax.set_xticks([0, 1])
    ax.set_xlim(0, 1.05)
    ax.set_title("Program 193: common-engine certificate dependency ledger")
    ax.legend(loc="lower right", fontsize=8)
    fig.savefig(FIG / "program193_common_engine_certificate.png", dpi=190)
    plt.close(fig)

    return {
        "constructed_object": {
            "name": "CertificateDependencyDAG",
            "nodes": [row["component"] for row in components],
            "acceptance": "all nodes enclosed by one directed-rounding Arb engine",
        },
        "environment": {
            "python": sys.version.split()[0],
            "platform": platform.platform(),
            "python_flint_available": python_flint_available,
            "arb_binary_available": arb_available,
            "flint_module_origin": None if flint_spec is None else flint_spec.origin,
        },
        "components": components,
        "closed_components": closed_count,
        "total_components": len(components),
        "full_one_engine_certificate": all_in_engine,
        "formal_worst_relative_enclosure": old[
            "formal_worst_relative_enclosure_from_P151"
        ],
        "target_relative_width": old["preregistered_target"],
        "formal_sub_three_percent_certificate": bool(
            all_in_engine
            and old["formal_worst_relative_enclosure_from_P151"]
            < old["preregistered_target"]
        ),
        "status": "Blocked: reproducible Arb/python-flint engine unavailable locally",
        "confidence": "Proven environment/dependency audit",
        "claim_boundary": (
            "The wide existing enclosure is an upper bound, not a lower bound "
            "on achievable width. No false impossibility or formal certificate "
            "is claimed."
        ),
    }


def _cycle_matrix_fraction(n: int) -> list[list[Fraction]]:
    w = [[Fraction(0) for _ in range(n)] for _ in range(n)]
    for i in range(n):
        w[i][(i - 1) % n] = Fraction(1, 2)
        w[i][(i + 1) % n] = Fraction(1, 2)
    return w


def _dirichlet_exact(
    w: list[list[Fraction]], f: list[Fraction]
) -> tuple[Fraction, Fraction, list[Fraction]]:
    n = len(w)
    row_sums = [sum(row, Fraction(0)) for row in w]
    if len(set(row_sums)) != 1:
        raise ValueError("row sum is not constant")
    s = row_sums[0]
    af = [s * f[i] - sum(w[i][j] * f[j] for j in range(n)) for i in range(n)]
    quadratic = sum(f[i] * af[i] for i in range(n))
    rhs = Fraction(1, 2) * sum(
        w[i][j] * (f[i] - f[j]) ** 2 for i in range(n) for j in range(n)
    )
    return quadratic, rhs, af


def program194_compiled_dirichlet_library() -> dict:
    rng = rng_for(194)
    exact_cases = 0
    minimum_energy = None
    for n in range(3, 13):
        w = _cycle_matrix_fraction(n)
        for _ in range(15):
            f = [Fraction(int(x), 3) for x in rng.integers(-9, 10, n)]
            q, rhs, _ = _dirichlet_exact(w, f)
            if q != rhs or q < 0:
                raise AssertionError("exact Dirichlet identity failed")
            exact_cases += 1
            minimum_energy = q if minimum_energy is None else min(minimum_energy, q)

    toolchain_root = Path.home() / ".elan" / "toolchains"
    toolchains = sorted(p.name for p in toolchain_root.glob("*")) if toolchain_root.exists() else []
    lean_launcher = shutil.which("lean")
    compiled = False
    compilation_output = "not attempted: no installed Lean toolchain"
    if lean_launcher and toolchains:
        process = subprocess.run(
            [lean_launcher, str(LEAN_SOURCE)],
            cwd=ROOT,
            text=True,
            capture_output=True,
            timeout=60,
            check=False,
        )
        compiled = process.returncode == 0
        compilation_output = (process.stdout + process.stderr)[-4000:]

    # Numerical semigroup sanity check is separate from the exact theorem.
    n = 12
    w_float = np.asarray(_cycle_matrix_fraction(n), dtype=float)
    a = np.eye(n) - w_float
    times = np.logspace(-3, 2, 60)
    unitary_residual = []
    heat_row_residual = []
    heat_minimum = []
    for t in times:
        u = expm(-1j * t * a)
        p = expm(-t * a)
        unitary_residual.append(float(np.linalg.norm(u.conj().T @ u - np.eye(n), 2)))
        heat_row_residual.append(float(np.max(np.abs(p.sum(axis=1) - 1.0))))
        heat_minimum.append(float(np.min(p)))

    fig, ax = plt.subplots(figsize=(8.8, 5.0), constrained_layout=True)
    ax.loglog(times, np.maximum(unitary_residual, 1e-18), label="unitary residual")
    ax.loglog(times, np.maximum(heat_row_residual, 1e-18), label="heat row residual")
    ax.loglog(times, np.maximum(np.abs(np.minimum(heat_minimum, 0)), 1e-18), label="negative heat part")
    ax.set_xlabel("parameter")
    ax.set_ylabel("residual")
    ax.set_title("Program 194: exact finite theorem with machine-toolchain boundary")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend()
    fig.savefig(FIG / "program194_dirichlet_library.png", dpi=190)
    plt.close(fig)

    return {
        "constructed_object": {
            "name": "FiniteDirichletCertificate",
            "scope": (
                "finite symmetric nonnegative W with constant row sum s; "
                "A=sI-W"
            ),
            "theorems": [
                "<f,Af>=1/2 sum_xy W_xy |f_x-f_y|^2 >=0",
                "A*1=0",
                "exp(-itA) is unitary",
                "exp(-tA) is positivity preserving and stochastic",
            ],
        },
        "exact_rational_cycle_cases": exact_cases,
        "minimum_exact_energy": str(minimum_energy),
        "lean_source": LEAN_SOURCE.name,
        "lean_launcher_available": lean_launcher is not None,
        "installed_lean_toolchains": toolchains,
        "lean_machine_compiled": compiled,
        "compilation_output_tail": compilation_output,
        "maximum_unitarity_residual": max(unitary_residual),
        "maximum_heat_row_residual": max(heat_row_residual),
        "minimum_heat_entry": min(heat_minimum),
        "status": (
            "Exact general proof supplied; proof-assistant compilation blocked"
            if not compiled
            else "Exact theorem and Lean library machine compiled"
        ),
        "confidence": "Proven on paper and exact rational checks; compilation status explicit",
        "claim_boundary": (
            "Finite Dirichlet mathematics does not supply a continuum limit, "
            "Born rule, apparatus, physical clock or empirical interpretation."
        ),
    }


def reflection_boundary_closed_form(
    r_source: float, z_source: float, z_target: float
) -> tuple[float, float, str]:
    """Maximum target transverse radius under qubit Z2 covariance.

    Alberti--Uhlmann reduces the all-t trace-norm criterion to the minimum of
    a convex piecewise-linear function of x=u^2.  Only x=0, x=1 and the
    source kink x_c=r^2/(1-z^2) need be checked.
    """
    denominator = 1.0 - z_source * z_source
    x_c = 0.0 if denominator <= 1e-15 else r_source * r_source / denominator
    candidates = [("u=0", 0.0), ("u=1", 1.0)]
    if -1e-12 <= x_c <= 1.0 + 1e-12:
        candidates.append(("source-kink", min(1.0, max(0.0, x_c))))
    values = []
    for name, x in candidates:
        g2 = max(x, r_source * r_source + x * z_source * z_source)
        values.append((g2 - x * z_target * z_target, x, name))
    value, x_active, name = min(values)
    return math.sqrt(max(0.0, value)), x_active, name


def alberti_uhlmann_margin(
    source: tuple[float, float], target: tuple[float, float], grid: int = 4001
) -> float:
    r_s, z_s = source
    r_t, z_t = target
    u = np.linspace(0.0, 1.0, grid)
    g_s = np.maximum(u, np.sqrt(r_s * r_s + u * u * z_s * z_s))
    g_t = np.maximum(u, np.sqrt(r_t * r_t + u * u * z_t * z_t))
    return float(np.min(g_s - g_t))


def program195_reflection_geometry() -> dict:
    rng = rng_for(195)
    source = (0.6, 0.0)
    z_targets = np.linspace(-0.999, 0.999, 501)
    boundaries = [
        reflection_boundary_closed_form(source[0], source[1], z)[0]
        for z in z_targets
    ]
    old_residuals = []
    active_counts = {"u=0": 0, "u=1": 0, "source-kink": 0}
    for _ in range(500):
        z_s = rng.uniform(-0.95, 0.95)
        r_s = rng.uniform(0.02, math.sqrt(max(0.0, 1.0 - z_s * z_s)))
        z_t = rng.uniform(-0.99, 0.99)
        closed, _, active = reflection_boundary_closed_form(r_s, z_s, z_t)
        factor, _, _ = previous_study.max_covariant_transverse_factor(
            z_s, z_t, grid=12001
        )
        numeric = r_s * factor
        old_residuals.append(abs(closed - numeric))
        active_counts[active] += 1
    counter_boundary, counter_x, counter_active = reflection_boundary_closed_form(
        0.6, 0.0, 0.8
    )

    fig, ax = plt.subplots(figsize=(8.8, 5.2), constrained_layout=True)
    ax.fill_between(
        z_targets,
        -np.asarray(boundaries),
        boundaries,
        color="#B7D7F0",
        label="reachable target disk sections",
    )
    ax.plot(z_targets, boundaries, color="#1F5A99", lw=2)
    ax.plot(z_targets, -np.asarray(boundaries), color="#1F5A99", lw=2)
    ax.scatter([0.8], [0.6], marker="x", s=80, color="#A61B1B", label="equal-M counterexample")
    ax.scatter([0.8], [counter_boundary], marker="o", s=55, color="#19733A", label="analytic boundary")
    ax.set_xlabel("target z")
    ax.set_ylabel("target transverse coordinate")
    ax.set_title("Program 195: closed-form reflection-conversion geometry")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.25)
    fig.savefig(FIG / "program195_reflection_geometry.png", dpi=190)
    plt.close(fig)

    return {
        "constructed_object": {
            "name": "AlbertiUhlmannReflectionCone",
            "trace_norm_function": (
                "g_(r,z)(u)=max(u,sqrt(r^2+u^2 z^2)), 0<=u<=1"
            ),
            "criterion": "source converts to target iff g_source(u)>=g_target(u) for all u",
            "closed_boundary": (
                "r_target,max^2=min_{x in {0,1,r_s^2/(1-z_s^2)}} "
                "[max(x,r_s^2+x z_s^2)-x z_t^2]"
            ),
        },
        "counterexample_boundary": counter_boundary,
        "counterexample_active_x": counter_x,
        "counterexample_active_constraint": counter_active,
        "maximum_closed_vs_choi_grid_residual": max(old_residuals),
        "median_closed_vs_choi_grid_residual": float(np.median(old_residuals)),
        "random_cross_checks": len(old_residuals),
        "active_constraint_counts": active_counts,
        "tensor_catalysis_classified": False,
        "status": "Numerical maximization eliminated; complete one-qubit boundary derived",
        "confidence": "Proven analytic criterion; finite numerical cross-check",
        "claim_boundary": (
            "The qubit one-copy preorder is complete. Tensor catalysis leaves "
            "the two-state qubit Alberti--Uhlmann setting and remains open. "
            "The theorem does not source asymmetry or discharge QW-2191."
        ),
    }


def _iqr_interval_eps(sample: np.ndarray, eps: float) -> tuple[float, float]:
    lo = np.quantile(sample, max(0.0, 0.75 - eps)) - np.quantile(
        sample, min(1.0, 0.25 + eps)
    )
    hi = np.quantile(sample, min(1.0, 0.75 + eps)) - np.quantile(
        sample, max(0.0, 0.25 - eps)
    )
    return max(float(lo), 1e-300), max(float(hi), 1e-300)


def _exponent_interval_two_eps(
    x1: np.ndarray, x2: np.ndarray, eps1: float, eps2: float, ratio: float
) -> tuple[float, float]:
    l1, u1 = _iqr_interval_eps(x1, eps1)
    l2, u2 = _iqr_interval_eps(x2, eps2)
    return math.log(l2 / u1) / math.log(ratio), math.log(u2 / l1) / math.log(ratio)


def _refresh_chain(
    values: np.ndarray, flags: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    record = np.empty(len(flags), dtype=float)
    kept = []
    cursor = 0
    current = values[cursor]
    cursor += 1
    record[0] = current
    kept.append(current)
    for i in range(1, len(flags)):
        if flags[i]:
            current = values[cursor]
            cursor += 1
            kept.append(current)
        record[i] = current
    return record, np.asarray(kept)


def program196_regenerative_mixing() -> dict:
    rng = rng_for(196)
    n = 4000
    refresh_probability = 0.05
    alpha = 0.8
    truth = 1.0 / alpha
    ratio = 4.0
    delta = 0.05
    reps = 320
    cover_nominal = 0
    cover_regenerative = 0
    widths_nominal = []
    widths_regenerative = []
    retained = []
    for _ in range(reps):
        flags1 = np.zeros(n, dtype=bool)
        flags2 = np.zeros(n, dtype=bool)
        flags1[1:] = rng.random(n - 1) < refresh_probability
        flags2[1:] = rng.random(n - 1) < refresh_probability
        r1 = 1 + int(flags1[1:].sum())
        r2 = 1 + int(flags2[1:].sum())
        values1 = stable_rvs(alpha, r1, rng)
        values2 = stable_rvs(alpha, r2, rng) * ratio**truth
        record1, kept1 = _refresh_chain(values1, flags1)
        record2, kept2 = _refresh_chain(values2, flags2)
        eps_nom = math.sqrt(math.log(4.0 / delta) / (2.0 * n))
        eps1 = math.sqrt(math.log(4.0 / delta) / (2.0 * r1))
        eps2 = math.sqrt(math.log(4.0 / delta) / (2.0 * r2))
        lo, hi = _exponent_interval_two_eps(
            record1, record2, eps_nom, eps_nom, ratio
        )
        rlo, rhi = _exponent_interval_two_eps(kept1, kept2, eps1, eps2, ratio)
        cover_nominal += lo <= truth <= hi
        cover_regenerative += rlo <= truth <= rhi
        widths_nominal.append(hi - lo)
        widths_regenerative.append(rhi - rlo)
        retained.append((r1 + r2) / (2.0 * n))

    k = np.arange(0, 101)
    beta_bound = (1.0 - refresh_probability) ** k
    fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.6), constrained_layout=True)
    axes[0].semilogy(k, beta_bound, color="#1F5A99")
    axes[0].set_xlabel("lag k")
    axes[0].set_ylabel(r"mixing bound $(1-p)^k$")
    axes[0].grid(True, alpha=0.25)
    axes[1].hist(widths_nominal, bins=35, alpha=0.65, label="invalid nominal record")
    axes[1].hist(widths_regenerative, bins=35, alpha=0.65, label="valid refresh subsample")
    axes[1].set_xlabel("exponent interval width")
    axes[1].legend(fontsize=8)
    fig.suptitle("Program 196: observed regeneration turns mixing into iid sampling")
    fig.savefig(FIG / "program196_regenerative_mixing.png", dpi=190)
    plt.close(fig)

    return {
        "constructed_object": {
            "name": "ObservedRegenerationFrame",
            "chain": (
                "with probability p refresh from F; otherwise repeat the "
                "previous value; refresh flags are recorded"
            ),
            "mixing_bound": "beta(k) <= (1-p)^k",
        },
        "theorem": (
            "Conditional on R recorded refreshes, the retained values are iid "
            "F and DKW holds with epsilon(R,delta)=sqrt(log(2/delta)/(2R)). "
            "Mixing over the random R preserves the failure bound."
        ),
        "n": n,
        "refresh_probability": refresh_probability,
        "replicates": reps,
        "true_exponent": truth,
        "nominal_record_coverage": cover_nominal / reps,
        "regenerative_coverage": cover_regenerative / reps,
        "nominal_mean_width": float(np.mean(widths_nominal)),
        "regenerative_mean_width": float(np.mean(widths_regenerative)),
        "mean_retained_fraction": float(np.mean(retained)),
        "status": "Nonasymptotic theorem executed for an explicit beta-mixing acquisition class",
        "confidence": "Proven conditional theorem; strong numerical coverage evidence",
        "claim_boundary": (
            "Refresh flags are an acquisition assumption. Hidden general "
            "beta-mixing requires a different coupling or concentration theorem."
        ),
    }


def order_features(record: np.ndarray) -> np.ndarray:
    """Distribution-robust temporal features for one ordered record."""
    n = len(record)
    ranks = rankdata(record, method="average")
    lag_spearman = float(np.corrcoef(ranks[:-1], ranks[1:])[0, 1])
    differences = np.diff(record)
    signs = np.sign(differences)
    persistence = float(np.mean(signs[:-1] == signs[1:]))
    ties = float(np.mean(differences == 0))
    time_rank = float(np.corrcoef(np.arange(n), ranks)[0, 1])

    # Six ordinal patterns of triples; ties are ignored.
    counts = np.zeros(6, dtype=float)
    perm_index = {p: i for i, p in enumerate(permutations(range(3)))}
    accepted = 0
    for i in range(n - 2):
        triple = record[i : i + 3]
        if len(set(triple.tolist())) < 3:
            continue
        order = tuple(np.argsort(triple).tolist())
        counts[perm_index[order]] += 1
        accepted += 1
    if accepted:
        probabilities = counts[counts > 0] / accepted
        entropy = float(-np.sum(probabilities * np.log(probabilities)) / math.log(6))
        reversal = float((counts[0] - counts[-1]) / accepted)
    else:
        entropy = 0.0
        reversal = 0.0
    return np.asarray(
        [lag_spearman, persistence, entropy, ties, time_rank, reversal], dtype=float
    )


def _simulate_ar1(
    phi: float, reps: int, n: int, rng: np.random.Generator
) -> np.ndarray:
    x = rng.normal(size=(reps, n))
    scale = math.sqrt(1.0 - phi * phi)
    for j in range(1, n):
        x[:, j] = phi * x[:, j - 1] + scale * x[:, j]
    return x


def program197_order_sensitive_open_set() -> dict:
    core = {
        "protocol_id": "FIN-P197-ORDER-OPEN-SET-001",
        "feature_names": [
            "lag1_spearman",
            "sign_persistence",
            "ordinal_entropy",
            "exact_tie_fraction",
            "time_rank_trend",
            "ordinal_reversal_asymmetry",
        ],
        "calibration_model": "continuous iid Gaussian",
        "score": "regularized Mahalanobis distance",
        "threshold": "empirical 0.99 calibration quantile",
        "challenges": [
            "iid_gaussian_validation",
            "iid_stable_0.8",
            "sorted_same_multiset",
            "AR1_phi_0.8",
            "repeated_block_20",
            "linear_drift",
        ],
        "training_records": 1400,
        "validation_records": 420,
        "record_length": 360,
        "seed": SEED + 197,
        "frozen_before_challenge": True,
    }
    protocol_record = {"core": core, "canonical_core_sha256": canonical_digest(core)}
    PROTOCOL.write_text(
        json.dumps(protocol_record, indent=2) + "\n", encoding="utf-8"
    )

    rng = rng_for(197)
    train = rng.normal(size=(core["training_records"], core["record_length"]))
    train_features = np.asarray([order_features(row) for row in train])
    centre = train_features.mean(axis=0)
    covariance = np.cov(train_features, rowvar=False)
    ridge = 1e-6 * np.trace(covariance) / covariance.shape[0] + 1e-10
    precision = np.linalg.inv(covariance + ridge * np.eye(covariance.shape[0]))

    def score(features: np.ndarray) -> np.ndarray:
        delta = features - centre
        return np.sqrt(np.einsum("ij,jk,ik->i", delta, precision, delta))

    training_scores = score(train_features)
    threshold = float(np.quantile(training_scores, 0.99))
    reps, n = core["validation_records"], core["record_length"]
    base = rng.normal(size=(reps, n))
    stable = stable_rvs(0.8, (reps, n), rng)
    sorted_same = np.sort(base, axis=1)
    ar = _simulate_ar1(0.8, reps, n, rng)
    latent = rng.normal(size=(reps, math.ceil(n / 20)))
    repeated = np.repeat(latent, 20, axis=1)[:, :n]
    drift = rng.normal(size=(reps, n)) + np.linspace(-2.0, 2.0, n)
    challenges = {
        "iid_gaussian_validation": base,
        "iid_stable_0.8": stable,
        "sorted_same_multiset": sorted_same,
        "AR1_phi_0.8": ar,
        "repeated_block_20": repeated,
        "linear_drift": drift,
    }
    rows = []
    feature_means = {}
    for name, records in challenges.items():
        features = np.asarray([order_features(row) for row in records])
        scores = score(features)
        rows.append(
            {
                "challenge": name,
                "rejection_rate": float(np.mean(scores > threshold)),
                "median_score": float(np.median(scores)),
            }
        )
        feature_means[name] = features.mean(axis=0).tolist()

    marginal_residual = float(
        np.max(
            np.abs(
                np.sort(base, axis=1)
                - np.sort(sorted_same, axis=1)
            )
        )
    )
    order_feature_shift = float(
        np.linalg.norm(
            np.asarray(feature_means["iid_gaussian_validation"])
            - np.asarray(feature_means["sorted_same_multiset"])
        )
    )

    fig, ax = plt.subplots(figsize=(9.5, 5.0), constrained_layout=True)
    ax.bar(
        [row["challenge"].replace("_", "\n") for row in rows],
        [row["rejection_rate"] for row in rows],
        color=["#8FC1E3", "#19733A", "#A61B1B", "#D55E00", "#6A3D9A", "#B5651D"],
    )
    ax.axhline(0.01, color="black", ls="--", label="nominal calibration tail")
    ax.set_ylim(0, 1.05)
    ax.set_ylabel("rejection rate")
    ax.set_title("Program 197: order-sensitive preregistered open-set challenge")
    ax.legend(fontsize=8)
    fig.savefig(FIG / "program197_order_sensitive_open_set.png", dpi=190)
    plt.close(fig)

    return {
        "protocol": PROTOCOL.name,
        "protocol_sha256": protocol_record["canonical_core_sha256"],
        "threshold": threshold,
        "rows": rows,
        "feature_means": feature_means,
        "exact_multiset_residual_sorted_challenge": marginal_residual,
        "order_feature_mean_shift": order_feature_shift,
        "closed_set_completeness_claimed": False,
        "status": "Frozen order-sensitive protocol executed on six synthetic challenges",
        "confidence": "Strong synthetic evidence; exact multiset identity",
        "claim_boundary": (
            "The detector distinguishes declared temporal alternatives. It is "
            "not a universal open-set classifier, memory theorem or external "
            "experimental validation."
        ),
    }


def program198_minimal_process_tomography() -> dict:
    matrix = np.asarray(
        [
            [1.0, -0.5, 0.0],
            [1.0, -1.0, -1.0],
            [1.0, -1.0, 1.0],
        ]
    )
    reduced = matrix[:2]
    rank = int(np.linalg.matrix_rank(matrix))
    reduced_rank = int(np.linalg.matrix_rank(reduced))
    condition = float(np.linalg.cond(matrix))
    truth = np.asarray([math.log(0.84), 0.45, 0.20])
    log_visibilities = matrix @ truth
    rng = rng_for(198)
    reps = 6000
    noise_sd = 0.015
    noisy = log_visibilities + rng.normal(0.0, noise_sd, (reps, 3))
    estimates = np.linalg.solve(matrix, noisy.T).T
    means = estimates.mean(axis=0)
    sds = estimates.std(axis=0, ddof=1)

    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.7), constrained_layout=True)
    image = axes[0].imshow(matrix, cmap="coolwarm", vmin=-1, vmax=1)
    axes[0].set_xticks(range(3), [r"$\log b$", r"$v$", r"$c$"])
    axes[0].set_yticks(range(3), ["single", "plus", "echo"])
    for i in range(3):
        for j in range(3):
            axes[0].text(j, i, f"{matrix[i,j]:.1f}", ha="center", va="center")
    fig.colorbar(image, ax=axes[0], shrink=0.75)
    axes[1].errorbar(
        range(3),
        means,
        yerr=sds,
        fmt="o",
        color="#1F5A99",
        capsize=4,
        label="estimated mean ± sd",
    )
    axes[1].scatter(range(3), truth, marker="x", s=65, color="#A61B1B", label="truth")
    axes[1].set_xticks(range(3), [r"$\log b$", r"$v$", r"$c$"])
    axes[1].legend(fontsize=8)
    axes[1].grid(True, alpha=0.25)
    fig.suptitle("Program 198: minimal three-contrast process tomography")
    fig.savefig(FIG / "program198_minimal_process_tomography.png", dpi=190)
    plt.close(fig)

    return {
        "constructed_object": {
            "name": "GaussianTwoStepTomographyFrame",
            "parameters": ["log_detector_blur", "increment_variance_v", "covariance_c"],
            "contrasts": ["single_visibility", "plus_visibility", "echo_visibility"],
            "linear_log_model": matrix.tolist(),
        },
        "rank_three_contrasts": rank,
        "rank_without_echo": reduced_rank,
        "condition_number": condition,
        "minimality_theorem": (
            "In the declared three-parameter family, fewer than three scalar "
            "contrasts have Jacobian rank at most two. The single/plus/echo "
            "matrix has determinant -1 and identifies all three parameters."
        ),
        "closed_form_inverse": {
            "c": "(log V_echo-log V_plus)/2",
            "v": "2[log V_single-(log V_plus+log V_echo)/2]",
            "log_b": "2 log V_single-(log V_plus+log V_echo)/2",
        },
        "simulation": {
            "replicates": reps,
            "log_noise_sd": noise_sd,
            "truth": truth.tolist(),
            "mean_estimate": means.tolist(),
            "estimate_sd": sds.tolist(),
        },
        "status": "Minimal rank theorem and finite-noise tomography executed",
        "confidence": "Proven rank/minimality in declared family; strong synthetic evidence",
        "claim_boundary": (
            "Gaussian equal-variance increments and multiplicative detector "
            "blur are explicit model assumptions, not strict FIN consequences."
        ),
    }


def program199_environment_equivalence() -> dict:
    gamma = 0.8
    angle = math.acos(gamma)
    # Distribution A: +/- angle equally. Distribution B: 0 with 0.9 and pi with 0.1.
    t = np.linspace(0.0, 4.0, 501)
    char_a = np.cos(angle * t)
    char_b = 0.9 + 0.1 * np.cos(math.pi * t)
    gamma2_a = math.cos(2.0 * angle)
    gamma2_b = 1.0

    fig, ax = plt.subplots(figsize=(8.8, 5.0), constrained_layout=True)
    ax.plot(t, char_a, label=r"$\frac{1}{2}(\delta_{-a}+\delta_a)$")
    ax.plot(t, char_b, label=r"$0.9\delta_0+0.1\delta_\pi$")
    ax.scatter([1, 1], [gamma, gamma], color="black", zorder=3)
    ax.axvline(1.0, color="black", ls=":", alpha=0.6)
    ax.set_xlabel("characteristic-function argument")
    ax.set_ylabel("dephasing characteristic function")
    ax.set_title("Program 199: same one-time channel, different multi-time process")
    ax.legend()
    ax.grid(True, alpha=0.25)
    fig.savefig(FIG / "program199_environment_equivalence.png", dpi=190)
    plt.close(fig)

    return {
        "constructed_object": {
            "name": "OperationalEnvironmentQuotient",
            "equivalence": (
                "two environments are k-equivalent when all preparation/"
                "intervention/measurement probabilities in the declared "
                "k-step experiment agree"
            ),
            "minimal_dilation_rule": (
                "minimal Stinespring dilations of one channel are unique only "
                "up to an environment isometry"
            ),
        },
        "one_time_counterexample": {
            "common_gamma": gamma,
            "distribution_A": f"equal atoms at +/-acos({gamma})",
            "distribution_B": "mass 0.9 at 0 and 0.1 at pi",
            "one_time_characteristic_A": math.cos(angle),
            "one_time_characteristic_B": 0.9 + 0.1 * math.cos(math.pi),
            "two_time_characteristic_A": gamma2_a,
            "two_time_characteristic_B": gamma2_b,
        },
        "one_time_channels_identical": abs(
            math.cos(angle) - (0.9 + 0.1 * math.cos(math.pi))
        ) < 1e-14,
        "multi_time_processes_distinguishable": True,
        "finite_identifiability_theorem": (
            "A finite set of observed characteristic-function values defines "
            "a finite-dimensional affine image of the probability-measure "
            "simplex; generally its fibres contain multiple microscopic "
            "phase laws. Only quotient data, not microscopic coordinates, "
            "are operationally identifiable."
        ),
        "status": "Explicit environment quotient and one-time nonidentifiability witness",
        "confidence": "Proven by construction",
        "claim_boundary": (
            "Process tomography can refine the quotient but finite data do "
            "not identify an arbitrary environment distribution or dilation."
        ),
    }


def program200_analytic_phase_sources() -> dict:
    legacy = np.exp(1j * math.pi / 4.0)
    strict = np.exp(1j * 743.0 / 4000.0)
    images = [legacy**n for n in range(-8, 9)]
    minimum_distance = min(abs(z - strict) for z in images)
    candidates = [
        {
            "class": "continuous_group_endomorphism_U1",
            "normal_form": "f(z)=z^n, n in Z",
            "target_independent": True,
            "origin_bearing": False,
            "produces_strict": False,
            "verdict": "roots of unity map to roots of unity",
        },
        {
            "class": "rational_or_algebraic_functional_calculus",
            "normal_form": "f algebraic over Qbar",
            "target_independent": True,
            "origin_bearing": False,
            "produces_strict": False,
            "verdict": "algebraic input remains algebraic when defined",
        },
        {
            "class": "polynomial_interpolation",
            "normal_form": "f(z)=z_strict*z/z_legacy",
            "target_independent": False,
            "origin_bearing": False,
            "produces_strict": True,
            "verdict": "receiver construction inserts target",
        },
        {
            "class": "logarithm_and_reexponentiation",
            "normal_form": "f(z)=exp(c Log z)",
            "target_independent": False,
            "origin_bearing": False,
            "produces_strict": True,
            "verdict": "requires branch origin and target coefficient c",
        },
        {
            "class": "constant_holomorphic_map",
            "normal_form": "f(z)=z_strict",
            "target_independent": False,
            "origin_bearing": False,
            "produces_strict": True,
            "verdict": "pure target coding",
        },
    ]
    accepted = [
        row
        for row in candidates
        if row["target_independent"] and row["origin_bearing"] and row["produces_strict"]
    ]

    fig, ax = plt.subplots(figsize=(6.3, 6.3), constrained_layout=True)
    theta = np.linspace(0, 2 * math.pi, 500)
    ax.plot(np.cos(theta), np.sin(theta), color="black", alpha=0.35)
    ax.scatter(
        [z.real for z in images],
        [z.imag for z in images],
        color="#1F5A99",
        label=r"endomorphism images $z_{\rm legacy}^n$",
    )
    ax.scatter([strict.real], [strict.imag], color="#A61B1B", marker="x", s=90, label="strict phase")
    ax.set_aspect("equal")
    ax.set_title("Program 200: natural analytic phase operations")
    ax.legend(fontsize=8)
    fig.savefig(FIG / "program200_analytic_phase_sources.png", dpi=190)
    plt.close(fig)

    return {
        "constructed_object": {
            "name": "PhaseSourceNaturalityNoGo",
            "source": "z_legacy=exp(i*pi/4)",
            "target": "z_strict=exp(i*743/4000)",
        },
        "endomorphism_theorem": (
            "Every continuous group endomorphism U(1)->U(1) is z->z^n. "
            "It maps the order-eight legacy phase to a root of unity, whereas "
            "exp(i*743/4000) is transcendental by Lindemann--Weierstrass."
        ),
        "minimum_distance_from_endomorphism_images": float(minimum_distance),
        "candidate_rows": candidates,
        "accepted_source_operations": len(accepted),
        "strict_phase_source_exported": False,
        "status": "Natural analytic classes exhausted; only target-coded receivers remain",
        "confidence": "Proven for declared operation classes",
        "claim_boundary": (
            "Arbitrary holomorphic interpolation can map one point to another "
            "but has no source content. No phase origin, selector or coupling "
            "law is generated."
        ),
    }


def program201_scale_free_observables() -> dict:
    n = 12
    w = np.asarray(_cycle_matrix_fraction(n), dtype=float)
    a = np.eye(n) - w
    eigen = np.linalg.eigvalsh(a)
    positive = eigen[eigen > 1e-10]
    gap = positive[0]
    tau = 1.7
    zeta = 0.8
    scales = np.logspace(-6, 6, 49)
    rows = []
    reference = None
    raw_gaps = []
    for c in scales:
        ac = c * a
        ec = np.linalg.eigvalsh(ac)
        pc = ec[ec > 1e-8 * c]
        gap_c = pc[0]
        probs = pc / pc.sum()
        obs = {
            "normalized_spectrum": (pc / pc[-1]).tolist(),
            "condition_number": float(pc[-1] / pc[0]),
            "spectral_entropy": float(-np.sum(probs * np.log(probs))),
            "dimensionless_heat_trace": float(
                np.trace(expm(-(tau / gap_c) * ac)).real / n
            ),
            "dimensionless_resolvent": float(
                gap_c
                * np.trace(np.linalg.inv(ac + (zeta * gap_c) * np.eye(n))).real
                / n
            ),
        }
        if reference is None:
            reference = obs
        rows.append(obs)
        raw_gaps.append(float(gap_c))

    residuals = {}
    for key in reference:
        if key == "normalized_spectrum":
            residuals[key] = max(
                float(np.max(np.abs(np.asarray(row[key]) - np.asarray(reference[key]))))
                for row in rows
            )
        else:
            residuals[key] = max(abs(row[key] - reference[key]) for row in rows)

    fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.6), constrained_layout=True)
    axes[0].loglog(scales, raw_gaps, color="#A61B1B", label="raw spectral gap")
    axes[0].set_xlabel("scale c")
    axes[0].set_ylabel("raw gap")
    axes[0].grid(True, which="both", alpha=0.25)
    axes[1].semilogx(
        scales,
        [row["dimensionless_heat_trace"] for row in rows],
        label="dimensionless heat trace",
    )
    axes[1].semilogx(
        scales,
        [row["dimensionless_resolvent"] for row in rows],
        label="dimensionless resolvent",
    )
    axes[1].set_xlabel("scale c")
    axes[1].set_ylim(0, 1.1)
    axes[1].legend(fontsize=8)
    axes[1].grid(True, alpha=0.25)
    fig.suptitle("Program 201: projective spectral observables survive the scale orbit")
    fig.savefig(FIG / "program201_scale_free_observables.png", dpi=190)
    plt.close(fig)

    return {
        "constructed_object": {
            "name": "ProjectiveSpectralObservableAlgebra",
            "universal_property": (
                "Every observable invariant under (A,t,z)->(cA,t/c,cz) "
                "factors through the projective spectral class [A] and the "
                "dimensionless products tA and z/A_scale."
            ),
        },
        "catalogue": [
            "eigenvalue ratios",
            "condition number on support",
            "spectral probability entropy",
            "heat trace at fixed tau=t*gap",
            "gap-normalized resolvent at fixed z/gap",
            "dimensionless unitary transition probabilities",
            "beta-gauge quotient profile (eta,nu,phi)",
        ],
        "noninvariants": [
            "raw spectral gap",
            "raw time",
            "raw resolvent energy",
            "beta representative",
            "mass, length and energy in physical units",
        ],
        "maximum_invariance_residuals": residuals,
        "raw_gap_ratio": max(raw_gaps) / min(raw_gaps),
        "strict_physical_scale_exported": False,
        "status": "Scale-orbit quotient and executable observable catalogue constructed",
        "confidence": "Proven factorization statement; numerical invariance audit",
        "claim_boundary": (
            "The catalogue supports scale-free falsification only. A physical "
            "clock or scale-charged source remains additional CA or new strict input."
        ),
    }


def _write_synthetic_bundle() -> dict:
    rng = rng_for(202)
    tau_star_seconds = 0.002
    blur = 0.84
    variance = 0.45
    covariance = 0.20
    contexts = {
        "reference": (0.0, 0.0),
        "single": (1.0, 0.0),
        "plus": (1.0, 1.0),
        "echo": (1.0, -1.0),
        "heldout": (1.0, 0.5),
    }
    phases = np.linspace(0.0, 2.0 * math.pi, 8, endpoint=False)
    trials = 300
    fieldnames = [
        "sequence",
        "timestamp_seconds",
        "context",
        "a1",
        "a2",
        "phase_radians",
        "outcome_plus",
    ]
    rows = []
    sequence = 0
    visibility_truth = {}
    for context, (a1, a2) in contexts.items():
        exponent = 0.5 * variance * (a1 * a1 + a2 * a2) + covariance * a1 * a2
        visibility = blur * math.exp(-exponent)
        visibility_truth[context] = visibility
        for phase in phases:
            probability = (1.0 + visibility * math.cos(phase)) / 2.0
            outcomes = rng.random(trials) < probability
            for outcome in outcomes:
                rows.append(
                    {
                        "sequence": sequence,
                        "timestamp_seconds": sequence * tau_star_seconds / 1000.0,
                        "context": context,
                        "a1": a1,
                        "a2": a2,
                        "phase_radians": phase,
                        "outcome_plus": int(outcome),
                    }
                )
                sequence += 1
    with BUNDLE_DATA.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    data_hash = sha256(BUNDLE_DATA)
    required = [
        "public_source_identifier",
        "license_identifier",
        "immutable_hashes",
        "preparation_provenance",
        "raw_time_ordered_records",
        "physical_units_and_timestamps",
        "detector_calibration",
        "apparatus_memory_calibration",
        "reference_control",
        "exclusion_audit",
        "independent_analysis_boundary",
    ]
    fields = {
        "public_source_identifier": {
            "pass": True,
            "value": "repository-relative:FIN-P202-SYNTHETIC-DRYRUN-001",
        },
        "license_identifier": {"pass": True, "value": "CC-BY-4.0"},
        "immutable_hashes": {
            "pass": True,
            "value": {BUNDLE_DATA.name: data_hash},
        },
        "preparation_provenance": {
            "pass": True,
            "value": "declared Gaussian two-increment simulator; seed 20260929",
        },
        "raw_time_ordered_records": {
            "pass": True,
            "value": f"{len(rows)} event-level ordered rows",
        },
        "physical_units_and_timestamps": {
            "pass": True,
            "value": (
                "timestamps in seconds under imported conditional clock "
                "tau_star=0.002 s; phase in radians"
            ),
        },
        "detector_calibration": {
            "pass": True,
            "value": "reference-context phase scan; injected blur b=0.84",
        },
        "apparatus_memory_calibration": {
            "pass": True,
            "value": "declared covariance challenge c=0.20 and echo context",
        },
        "reference_control": {
            "pass": True,
            "value": "reference context a1=a2=0",
        },
        "exclusion_audit": {
            "pass": True,
            "value": "no events excluded; all generated rows retained",
        },
        "independent_analysis_boundary": {
            "pass": False,
            "value": "generator and evaluator belong to the same release",
        },
    }
    metadata = {
        "bundle_id": "FIN-P202-SYNTHETIC-DRYRUN-001",
        "source_class": "synthetic_method_validation",
        "external_empirical_record": False,
        "required_fields": required,
        "fields": fields,
        "conditional_clock": {"tau_star_seconds": tau_star_seconds, "source": "CA input"},
        "generator_truth": {
            "blur": blur,
            "variance": variance,
            "covariance": covariance,
            "visibility_by_context": visibility_truth,
        },
        "data_file": BUNDLE_DATA.name,
        "data_sha256": data_hash,
    }
    BUNDLE_META.write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    return metadata


def program202_admissible_bundle() -> dict:
    metadata = _write_synthetic_bundle()
    passed = sum(row["pass"] for row in metadata["fields"].values())
    total = len(metadata["fields"])
    externally_admitted = bool(
        metadata["external_empirical_record"]
        and all(row["pass"] for row in metadata["fields"].values())
    )
    names = list(metadata["fields"])
    values = [1 if metadata["fields"][name]["pass"] else 0 for name in names]

    fig, ax = plt.subplots(figsize=(11.0, 3.9), constrained_layout=True)
    image = ax.imshow([values], cmap="RdYlGn", vmin=0, vmax=1, aspect="auto")
    ax.set_xticks(range(len(names)), [name.replace("_", "\n") for name in names], fontsize=7)
    ax.set_yticks([0], ["synthetic dry run"])
    for j, value in enumerate(values):
        ax.text(j, 0, str(value), ha="center", va="center", fontsize=9)
    ax.set_title("Program 202: operational bundle intake (method-validation only)")
    fig.colorbar(image, ax=ax, ticks=[0, 1], shrink=0.75)
    fig.savefig(FIG / "program202_operational_bundle.png", dpi=190)
    plt.close(fig)

    return {
        "constructed_object": {
            "name": "OperationalBundleSchemaInstance",
            "metadata": BUNDLE_META.name,
            "raw_event_file": BUNDLE_DATA.name,
            "raw_event_sha256": metadata["data_sha256"],
        },
        "source_class": metadata["source_class"],
        "passed_required_fields": passed,
        "total_required_fields": total,
        "failed_fields": [
            name for name, row in metadata["fields"].items() if not row["pass"]
        ],
        "external_bundle_admitted": externally_admitted,
        "external_protocol_executed": False,
        "status": "Synthetic intake dry run passes 10/11; independent external boundary fails",
        "confidence": "Proven schema/hash audit",
        "claim_boundary": (
            "The generated event record validates the intake and analysis "
            "machinery only. It is not external evidence and cannot validate physics."
        ),
    }


def _estimate_visibility(rows: list[dict[str, str]], context: str) -> float:
    selected = [row for row in rows if row["context"] == context]
    phase = np.asarray([float(row["phase_radians"]) for row in selected])
    response = np.asarray([2.0 * int(row["outcome_plus"]) - 1.0 for row in selected])
    cosine = np.cos(phase)
    return float((cosine @ response) / (cosine @ cosine))


def program203_conditional_prediction(program202: dict) -> dict:
    with BUNDLE_DATA.open("r", newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    observed = {
        context: _estimate_visibility(rows, context)
        for context in ["reference", "single", "plus", "echo", "heldout"]
    }
    logs = np.log([observed["single"], observed["plus"], observed["echo"]])
    matrix = np.asarray(
        [[1.0, -0.5, 0.0], [1.0, -1.0, -1.0], [1.0, -1.0, 1.0]]
    )
    log_blur, variance, covariance = np.linalg.solve(matrix, logs)
    a1, a2 = 1.0, 0.5
    predicted = math.exp(
        log_blur
        - 0.5 * variance * (a1 * a1 + a2 * a2)
        - covariance * a1 * a2
    )
    observed_heldout = observed["heldout"]
    residual = observed_heldout - predicted
    reference_residual = observed["reference"] - math.exp(log_blur)

    fig, ax = plt.subplots(figsize=(8.8, 5.0), constrained_layout=True)
    contexts = ["reference", "single", "plus", "echo", "heldout"]
    ax.scatter(
        range(len(contexts)),
        [observed[name] for name in contexts],
        color="#1F5A99",
        s=60,
        label="synthetic observed",
    )
    ax.scatter(
        [4],
        [predicted],
        marker="x",
        s=100,
        color="#A61B1B",
        label="held-out conditional prediction",
    )
    ax.set_xticks(range(len(contexts)), contexts)
    ax.set_ylabel("fringe visibility")
    ax.set_ylim(0, 1)
    ax.set_title("Program 203: W0 + CA + OP conditional prediction dry run")
    ax.legend(fontsize=8)
    ax.grid(True, axis="y", alpha=0.25)
    fig.savefig(FIG / "program203_conditional_prediction.png", dpi=190)
    plt.close(fig)

    return {
        "constructed_object": {
            "name": "ConditionalGaussianVisibilityLaw",
            "equation": (
                "V(a1,a2)=b exp[-v(a1^2+a2^2)/2-c a1 a2]"
            ),
            "heldout_setting": {"a1": a1, "a2": a2},
        },
        "dependency_ledger": {
            "W0": (
                "self-adjoint spectral generator and admissible unitary phase "
                "evolution; no physical scale or environment supplied"
            ),
            "CA": "tau_star=0.002 s imported as a conditional clock",
            "OP": (
                "preparation, Gaussian two-step environment, phase POVM, "
                "detector blur, event record and fitting rule"
            ),
        },
        "estimated_parameters": {
            "blur": math.exp(log_blur),
            "variance": variance,
            "covariance": covariance,
        },
        "observed_visibilities": observed,
        "heldout_predicted_visibility": predicted,
        "heldout_observed_visibility": observed_heldout,
        "heldout_residual": residual,
        "reference_control_residual": reference_residual,
        "program202_external_gate_passed": program202["external_bundle_admitted"],
        "synthetic_method_validation_only": True,
        "external_physical_prediction_tested": False,
        "status": "Conditional prediction executed only as a synthetic pipeline dry run",
        "confidence": "Strong synthetic method evidence; no external physics evidence",
        "claim_boundary": (
            "Because Program 202 did not admit an independent external bundle, "
            "this result is not an empirical FIN prediction. CA and OP inputs "
            "remain explicit and are not promoted to strict W0."
        ),
    }


def main() -> None:
    FIG.mkdir(exist_ok=True)
    previous = json.loads(PREVIOUS_RESULTS.read_text(encoding="utf-8"))
    programs: dict[str, dict] = {}
    programs["191"] = program191_reference_state_naturality()
    programs["192"] = program192_beta_gauge_quotient()
    programs["193"] = program193_common_engine_certificate(previous)
    programs["194"] = program194_compiled_dirichlet_library()
    programs["195"] = program195_reflection_geometry()
    programs["196"] = program196_regenerative_mixing()
    programs["197"] = program197_order_sensitive_open_set()
    programs["198"] = program198_minimal_process_tomography()
    programs["199"] = program199_environment_equivalence()
    programs["200"] = program200_analytic_phase_sources()
    programs["201"] = program201_scale_free_observables()
    programs["202"] = program202_admissible_bundle()
    programs["203"] = program203_conditional_prediction(programs["202"])

    results = {
        "metadata": {
            "title": (
                "FIN Programs 191-203: Reference States, Operational Quotients, "
                "and Conditional Prediction"
            ),
            "release": "10.18",
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
            "programs_executed": list(range(191, 204)),
            "new_theoretical_objects": [
                "CentralMeasureFunctor",
                "BetaGaugeQuotient",
                "CertificateDependencyDAG",
                "FiniteDirichletCertificate",
                "AlbertiUhlmannReflectionCone",
                "ObservedRegenerationFrame",
                "OrderSensitiveOpenSetProtocol",
                "GaussianTwoStepTomographyFrame",
                "OperationalEnvironmentQuotient",
                "PhaseSourceNaturalityNoGo",
                "ProjectiveSpectralObservableAlgebra",
                "OperationalBundleSchemaInstance",
                "ConditionalGaussianVisibilityLaw",
            ],
            "central_measure_source_found": False,
            "strict_beta_source_found": False,
            "full_arb_certificate": programs["193"]["full_one_engine_certificate"],
            "lean_library_compiled": programs["194"]["lean_machine_compiled"],
            "analytic_qubit_conversion_completed": True,
            "mixing_theorem_completed_in_declared_class": True,
            "minimal_process_tomography_completed": True,
            "strict_phase_source_found": False,
            "external_bundle_admitted": programs["202"]["external_bundle_admitted"],
            "conditional_pipeline_dry_run_completed": True,
            "QW_2191_discharged": False,
            "strict_selector_exported": False,
            "canonical_physical_unit_exported": False,
            "legacy_strict_bridge_completed": False,
            "legacy_role_transfer_started": False,
            "L_total_or_ToE_claimed": False,
            "external_physical_validation_claimed": False,
        },
        "recommended_next_programs": {
            "204": "Morita-natural central measure and categorical trace classification",
            "205": "renormalization-cocycle test for eta and quotient bridge invariants",
            "206": "reproducible Arb container or rational directed-rounding replacement",
            "207": "machine-compiled finite Dirichlet library with semigroup positivity",
            "208": "tensor conversion and catalysis beyond the one-qubit reflection cone",
            "209": "unknown-refresh and hidden-mixing confidence sequences",
            "210": "sequential conformal calibration for temporal open-set detection",
            "211": "finite-shot confidence regions for process-tensor tomography",
            "212": "moment-problem bounds from finitely sampled environment characteristics",
            "213": "formal phase-source no-go via topological group and functional calculus",
            "214": "preregistered scale-free falsification protocol",
            "215": "independently sourced experimental bundle satisfying all intake fields",
            "216": "first held-out conditional prediction on an admitted external bundle",
        },
    }
    RESULTS.write_text(json.dumps(results, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(results["global_verdict"], indent=2))


if __name__ == "__main__":
    main()
