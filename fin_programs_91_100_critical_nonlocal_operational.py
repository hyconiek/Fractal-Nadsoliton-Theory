#!/usr/bin/env python3
"""Execute FIN research Programs 91--100.

The suite closes several questions left open by Programs 81--90:

* a tail-controlled certificate that the strict lattice kernel is not a
  scalar-normalized fixed point of alternating-site Schur reduction;
* the exact critical mass and field rescaling needed for the nearest-neighbour
  local continuum route;
* identification of the coordinate-rescaled dense limit as a bounded
  convolution/graphon operator rather than a local Laplacian;
* polynomial propagation bounds for the strict long-range tail;
* a distinction between a projected operator-space gradient flow and a true
  quotient flow, together with a recheck of the P2772 learning candidate;
* a provenance intake for the post-P2721 selector/source artifacts;
* an optimized noisy two-slot process-tensor experiment;
* an apparatus-inclusive feedback ledger;
* a preregistered external-data acquisition calculation; and
* a necessary damping-completion atom for the legacy-to-strict bridge.

Every physical clock, detector, bath, preparation, measurement channel, unit,
and dataset remains explicit conditioning data.  The calculations do not
close QW-2191, source a dimensional standard, prove Lorentz invariance,
complete the legacy-to-strict bridge, transfer legacy physical roles, promote
L_total, or validate a physical theory.
"""

from __future__ import annotations

import hashlib
import itertools
import json
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.linalg import expm


ROOT = Path(__file__).resolve().parent
OUT = ROOT / "FIN_Programs_91_100_Critical_Nonlocal_Operational_Results.json"
INTAKE = ROOT / "FIN_Programs_91_100_External_Data_Intake_Template.json"
FIG = ROOT / "FIN_Programs_91_100_Critical_Nonlocal_Operational_Figures"
FIG.mkdir(exist_ok=True)

ALPHA_GEO = 4.0 * math.log(2.0)
BETA_TORS = 0.01
OMEGA_L = math.pi / 4.0
PHI_L = math.pi / 6.0
OMEGA_S = 0.18575
PHI_S = 0.16250
BETA_S = 1.0
ETA_S = 1.8
MASS2 = 0.25
SEED = 20260727


def k_legacy(d):
    d = np.asarray(d, dtype=float)
    return ALPHA_GEO * np.cos(OMEGA_L * d + PHI_L) / (
        1.0 + BETA_TORS * d
    )


def k_strict(d):
    d = np.asarray(d, dtype=float)
    return np.cos(OMEGA_S * d + PHI_S) / (
        1.0 + BETA_S * d**ETA_S
    )


def radial_matrix(kernel, n: int, *, absolute: bool = False) -> np.ndarray:
    delta = (
        np.arange(n)[:, None] - np.arange(n, dtype=int)[None, :]
    ) % n
    d = np.minimum(delta, n - delta)
    w = np.asarray(kernel(d), dtype=float)
    np.fill_diagonal(w, 0.0)
    return np.abs(w) if absolute else w


def laplacian(w: np.ndarray) -> np.ndarray:
    return np.diag(w.sum(axis=1)) - w


def sha256(path: Path) -> str | None:
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else None


def read_json(path: Path) -> dict:
    if not path.exists():
        return {"missing": True, "path": str(path.relative_to(ROOT))}
    return json.loads(path.read_text(encoding="utf-8"))


def interval_ratio_certificate(
    truncation: int = 1_000_000, mass2: float = MASS2
) -> dict:
    """Tail-controlled two-mode obstruction for the strict lattice profile.

    The normalized infinite-lattice jump weights are

        p_{+/-d}=a_d/Z,  a_d=|K_strict(d)|, d>=1.

    Alternating-site Schur reduction aliases q and q+pi by the harmonic mean.
    A scalar normalization cannot alter the ratio of the q=0 and q=pi
    eigenvalues.  The native and reduced ratios are enclosed with a rigorous
    absolute tail bound using a_d <= d^{-eta}.
    """

    d_int = np.arange(1, truncation + 1, dtype=np.int64)
    d = d_int.astype(float)
    a = np.abs(np.cos(OMEGA_S * d + PHI_S)) / (1.0 + d**ETA_S)
    partial_z = float(2.0 * a.sum())
    tail_z_bound = float(
        2.0 * truncation ** (1.0 - ETA_S) / (ETA_S - 1.0)
    )
    parity = np.where(d_int % 2 == 0, 1.0, -1.0)
    quarter = np.cos(0.5 * math.pi * d_int)
    p_pi_partial = float(2.0 * np.dot(a, parity) / partial_z)
    p_half_partial = float(2.0 * np.dot(a, quarter) / partial_z)

    # If |tail numerator| <= t and 0 <= tail denominator <= t, then
    # |N_inf/Z_inf - N_D/Z_D| <= 2t/Z_D.
    symbol_error = float(2.0 * tail_z_bound / partial_z)
    ppi = (p_pi_partial - symbol_error, p_pi_partial + symbol_error)
    phalf = (p_half_partial - symbol_error, p_half_partial + symbol_error)

    native_values = []
    schur_values = []
    differences = []
    for p_pi, p_h in itertools.product(ppi, phalf):
        lam_pi = mass2 + 1.0 - p_pi
        lam_half = mass2 + 1.0 - p_h
        native = mass2 / lam_pi
        harmonic_zero = 2.0 * mass2 * lam_pi / (mass2 + lam_pi)
        reduced = harmonic_zero / lam_half
        native_values.append(native)
        schur_values.append(reduced)
        differences.append(reduced - native)
    native_interval = [min(native_values), max(native_values)]
    schur_interval = [min(schur_values), max(schur_values)]
    gap_interval = [min(differences), max(differences)]
    margin = schur_interval[0] - native_interval[1]

    return {
        "truncation": truncation,
        "eta": ETA_S,
        "mass2": mass2,
        "partial_normalization": partial_z,
        "normalization_tail_upper_bound": tail_z_bound,
        "normalized_symbol_absolute_error_bound": symbol_error,
        "p_hat_pi_partial": p_pi_partial,
        "p_hat_pi_over_2_partial": p_half_partial,
        "native_zero_to_pi_ratio_interval": native_interval,
        "schur_zero_to_pi_ratio_interval": schur_interval,
        "difference_interval": gap_interval,
        "disjoint_interval_margin": margin,
        "certificate_passes": margin > 0.0,
    }


def program91_strict_nonclosure_certificate() -> dict:
    cert = interval_ratio_certificate()
    fig, ax = plt.subplots(figsize=(8.2, 4.5), constrained_layout=True)
    centers = [
        np.mean(cert["native_zero_to_pi_ratio_interval"]),
        np.mean(cert["schur_zero_to_pi_ratio_interval"]),
    ]
    errors = [
        0.5
        * (
            cert["native_zero_to_pi_ratio_interval"][1]
            - cert["native_zero_to_pi_ratio_interval"][0]
        ),
        0.5
        * (
            cert["schur_zero_to_pi_ratio_interval"][1]
            - cert["schur_zero_to_pi_ratio_interval"][0]
        ),
    ]
    ax.errorbar(
        ["native strict", "Schur retained"],
        centers,
        yerr=errors,
        fmt="o",
        capsize=7,
        color="#1f5a99",
    )
    ax.set_ylabel("scale-invariant zero/high-mode ratio")
    ax.set_title("Tail-certified strict-kernel nonclosure")
    ax.grid(True, axis="y", alpha=0.25)
    fig.savefig(FIG / "program91_strict_nonclosure_certificate.png", dpi=190)
    plt.close(fig)
    return {
        "certificate": cert,
        "theorem": (
            "For alternating-site elimination, the retained precision symbol "
            "is the harmonic mean of the aliased fine symbols. The ratio of "
            "two symbol values is invariant under every scalar field/operator "
            "normalization. The certified native and retained ratios are "
            "disjoint; hence the normalized infinite strict lattice family is "
            "not closed under one Schur step at mass2=1/4."
        ),
        "status": "Computer-assisted theorem with explicit analytic tail bound",
        "scope": (
            "This excludes scalar-normalized closure of this fixed-mass "
            "infinite-lattice realization. It is not a no-go for parameter "
            "running, additional counterterms, or a different coarse map."
        ),
    }


def nn_schur_parameters(mass: float, coupling: float) -> tuple[float, float]:
    denominator = mass + 2.0 * coupling
    coarse_coupling = coupling * coupling / denominator
    coarse_mass = mass * (mass + 4.0 * coupling) / denominator
    return coarse_mass, coarse_coupling


def program92_critical_tuning_and_field_rescaling() -> dict:
    rows = []
    for mu in [0.5, 1.0, 2.0]:
        for n in [64, 128, 256, 512, 1024, 2048, 4096]:
            coupling = 1.0
            ratio = mu * mu / (n * n)
            mass = coupling * ratio
            m2, c2 = nn_schur_parameters(mass, coupling)
            ratio2 = m2 / c2
            target = 4.0 * ratio
            rows.append(
                {
                    "mu": mu,
                    "fine_N": n,
                    "fine_mass_to_coupling": ratio,
                    "coarse_mass_to_coupling": ratio2,
                    "critical_target_for_N_over_2": target,
                    "absolute_tuning_defect": ratio2 - target,
                    "relative_tuning_defect": (ratio2 - target) / target,
                }
            )

    fig, ax = plt.subplots(figsize=(8.4, 4.7), constrained_layout=True)
    for mu in [0.5, 1.0, 2.0]:
        selected = [r for r in rows if r["mu"] == mu]
        ax.loglog(
            [r["fine_N"] for r in selected],
            [r["relative_tuning_defect"] for r in selected],
            "o-",
            label=f"mu={mu:g}",
        )
    ax.set_xlabel("fine lattice size N")
    ax.set_ylabel("relative one-step critical-tuning defect")
    ax.set_title("Nearest-neighbour critical tuning converges as N^-2")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend()
    fig.savefig(FIG / "program92_critical_tuning.png", dpi=190)
    plt.close(fig)

    m0, c0 = nn_schur_parameters(0.0, 1.0)
    return {
        "exact_parameter_map": {
            "mass_prime": "m(m+4c)/(m+2c)",
            "coupling_prime": "c^2/(m+2c)",
            "ratio_map": "r'=r(r+4), r=m/c",
        },
        "massless_check": {
            "coarse_mass": m0,
            "coarse_coupling": c0,
            "operator_rescaling_needed_to_restore_c": 2.0,
        },
        "critical_scaling": (
            "r_N=mu^2/N^2 gives r' = 4 r_N + r_N^2, while "
            "r_(N/2)=4 r_N. The relative defect is exactly mu^2/(4N^2)."
        ),
        "rows": rows,
        "status": "Exact nearest-neighbour theorem and finite verification",
        "scope": (
            "Critical mass tuning and the factor-two field/operator "
            "renormalization are additional scale-dependent prescriptions; "
            "the fixed FIN strict tuple does not source them."
        ),
    }


def coordinate_profile(x: np.ndarray) -> np.ndarray:
    d = np.minimum(x, 1.0 - x)
    return np.abs(k_strict(d))


def continuum_fourier_coefficients(max_mode: int, grid: int = 1 << 20):
    x = (np.arange(grid, dtype=float) + 0.5) / grid
    g = coordinate_profile(x)
    fft = np.fft.rfft(g)
    modes = np.arange(max_mode + 1, dtype=float)
    midpoint_phase = np.exp(-1j * math.pi * modes / grid)
    coeff = (fft[: max_mode + 1] * midpoint_phase / fft[0]).real
    return coeff


def discrete_coordinate_coefficients(n: int, max_mode: int) -> np.ndarray:
    x = np.arange(n, dtype=float) / n
    g = coordinate_profile(x)
    g[0] = 0.0
    return (np.fft.rfft(g)[: max_mode + 1] / g.sum()).real


def program93_graphon_continuum_limit() -> dict:
    max_mode = 4096
    continuum = continuum_fourier_coefficients(max_mode)
    rows = []
    for n in [128, 256, 512, 1024, 2048, 4096, 8192]:
        discrete = discrete_coordinate_coefficients(n, 32)
        error = float(np.max(np.abs(discrete - continuum[:33])))
        rows.append({"N": n, "first_33_mode_sup_error": error})
    slope = float(
        np.polyfit(
            np.log([r["N"] for r in rows]),
            np.log([r["first_33_mode_sup_error"] for r in rows]),
            1,
        )[0]
    )
    precision = MASS2 + 1.0 - continuum
    kinetic_ratio = [
        float((precision[k] - MASS2) / ((2.0 * math.pi * k) ** 2))
        for k in [1, 4, 16, 64, 256, 1024, 4096]
    ]

    fig, axes = plt.subplots(1, 2, figsize=(10.6, 4.4), constrained_layout=True)
    axes[0].loglog(
        [r["N"] for r in rows],
        [r["first_33_mode_sup_error"] for r in rows],
        "o-",
    )
    axes[0].set_xlabel("Nystrom size N")
    axes[0].set_ylabel("sup error, modes 0..32")
    axes[0].set_title("Convergence to bounded integral operator")
    axes[1].semilogx(
        np.arange(1, 257), precision[1:257], label="continuum precision"
    )
    axes[1].axhline(MASS2 + 1.0, color="black", ls="--", label="high-mode limit")
    axes[1].set_xlabel("Fourier mode")
    axes[1].set_ylabel("eigenvalue")
    axes[1].set_title("Bounded spectrum, unlike a local Laplacian")
    for ax in axes:
        ax.grid(True, which="both", alpha=0.25)
    axes[1].legend(fontsize=8)
    fig.savefig(FIG / "program93_graphon_continuum.png", dpi=190)
    plt.close(fig)

    return {
        "continuum_operator": (
            "(Tf)(x)=Z^{-1} integral_0^1 "
            "|K_strict(dist_circle(x,y))| f(y) dy; A=m^2 I+I-T"
        ),
        "nystrom_rows": rows,
        "observed_error_loglog_slope": slope,
        "selected_precision_eigenvalues": {
            str(k): float(precision[k])
            for k in [0, 1, 2, 4, 16, 64, 256, 1024, 4096]
        },
        "selected_kinetic_over_k_squared": dict(
            zip(
                [str(k) for k in [1, 4, 16, 64, 256, 1024, 4096]],
                kinetic_ratio,
            )
        ),
        "theorem": (
            "The coordinate-rescaled dense matrices are Nystrom "
            "discretizations of a bounded circle-convolution operator. Its "
            "Fourier coefficients tend to zero, so A_k tends to m^2+1; a "
            "local second-order Laplacian instead has unbounded k^2 growth."
        ),
        "status": "Strong numerical convergence plus standard Fourier theorem",
        "scope": (
            "This identifies a genuine nonlocal continuum limit, not a local "
            "spacetime or Lorentzian field theory."
        ),
    }


def lattice_tail_records(
    truncation: int = 1_000_000,
) -> tuple[list[dict], float, float]:
    d = np.arange(1, truncation + 1, dtype=float)
    a = np.abs(k_strict(d))
    partial_z = float(2.0 * a.sum())
    trunc_tail = float(
        2.0 * truncation ** (1.0 - ETA_S) / (ETA_S - 1.0)
    )
    reverse = 2.0 * np.cumsum(a[::-1])[::-1]
    records = []
    radii = [16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768]
    for radius in radii:
        numerator_partial = float(reverse[radius])
        lower = numerator_partial / (partial_z + trunc_tail)
        upper = (numerator_partial + trunc_tail) / partial_z
        analytic_upper = float(
            2.0
            * radius ** (1.0 - ETA_S)
            / ((ETA_S - 1.0) * partial_z)
        )
        records.append(
            {
                "R": radius,
                "tail_probability_lower": lower,
                "tail_probability_upper": upper,
                "analytic_power_upper": analytic_upper,
                "duhamel_wave_heat_norm_bound_t_0_5": min(
                    2.0, 2.0 * 0.5 * upper
                ),
            }
        )
    return records, partial_z, trunc_tail


def program94_long_range_propagation() -> dict:
    rows, z_partial, z_tail = lattice_tail_records()
    slope = float(
        np.polyfit(
            np.log([r["R"] for r in rows[-8:]]),
            np.log(
                [
                    0.5
                    * (r["tail_probability_lower"] + r["tail_probability_upper"])
                    for r in rows[-8:]
                ]
            ),
            1,
        )[0]
    )

    # A direct finite-cycle propagation witness computed spectrally.
    n = 8192
    d = np.minimum(np.arange(n), n - np.arange(n)).astype(float)
    weights = np.abs(k_strict(d))
    weights[0] = 0.0
    weights /= weights.sum()
    lam = 1.0 - np.fft.fft(weights).real
    t = 0.2
    wave = np.fft.ifft(np.exp(-1j * t * lam))
    heat = np.fft.ifft(np.exp(-t * lam)).real
    propagation = []
    for radius in [8, 16, 32, 64, 128, 256, 512]:
        far = d > radius
        propagation.append(
            {
                "R": radius,
                "wave_probability_outside_R": float(np.sum(np.abs(wave[far]) ** 2)),
                "heat_mass_outside_R": float(np.sum(heat[far])),
            }
        )

    fig, ax = plt.subplots(figsize=(8.4, 4.7), constrained_layout=True)
    ax.loglog(
        [r["R"] for r in rows],
        [r["tail_probability_upper"] for r in rows],
        "o-",
        label="certified tail upper bound",
    )
    ax.loglog(
        [r["R"] for r in rows],
        [r["analytic_power_upper"] for r in rows],
        "--",
        label="R^(1-eta) envelope",
    )
    ax.set_xlabel("truncation radius R")
    ax.set_ylabel("coupling mass beyond R")
    ax.set_title("Strict lattice propagation has a polynomial long-range tail")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend()
    fig.savefig(FIG / "program94_long_range_propagation.png", dpi=190)
    plt.close(fig)

    return {
        "partial_infinite_normalization": z_partial,
        "normalization_tail_bound": z_tail,
        "tail_rows": rows,
        "observed_tail_loglog_slope": slope,
        "predicted_tail_exponent": 1.0 - ETA_S,
        "finite_cycle_propagation_t_0_2": propagation,
        "theorem": (
            "For eta>1, the normalized coupling mass beyond R is at most "
            "2/[Z(eta-1)] R^(1-eta). If L_R truncates all longer edges, "
            "||exp(-zL)-exp(-zL_R)|| <= |z| ||L-L_R|| <= "
            "2|z| tail(R) for z=t or z=it. For eta=1.8 this is polynomial "
            "R^-0.8, not a strict finite propagation cone."
        ),
        "status": "Analytic propagation bound and finite spectral check",
    }


def project_symmetric_circulant_zero_diagonal(matrix: np.ndarray) -> np.ndarray:
    n = matrix.shape[0]
    row = np.array(
        [
            np.mean([matrix[i, (i + k) % n] for i in range(n)])
            for k in range(n)
        ],
        dtype=float,
    )
    for k in range(1, n):
        value = 0.5 * (row[k] + row[-k])
        row[k] = value
        row[-k] = value
    row[0] = 0.0
    return np.array([np.roll(row, i) for i in range(n)])


def moment_functional(k: np.ndarray, c: np.ndarray, a: float, b: float) -> float:
    return float(
        0.5 * a * np.trace(k @ k)
        + 0.25 * b * np.trace(k @ k @ k @ k)
        - np.trace(k @ c)
    )


def finite_cgeo_loss(name: str, parameters: dict[str, float]) -> float:
    n = 13
    shells = range(7)

    def kv(d: int) -> float:
        if name == "legacy":
            return (
                parameters["alpha_geo"]
                * math.cos(parameters["omega"] * d + parameters["phi"])
                / (1.0 + parameters["beta_tors"] * d)
            )
        return math.cos(parameters["omega"] * d + parameters["phi"]) / (
            1.0 + parameters["beta"] * d ** parameters["eta"]
        )

    source = np.array([kv(d) for d in shells])
    cyclic = np.array([source[min(x, n - x)] for x in range(n)])
    coupled_full = np.fft.ifft(np.fft.fft(cyclic) ** 2).real
    coupled = np.array(
        [
            np.mean(
                [coupled_full[x] for x in range(n) if min(x, n - x) == d]
            )
            for d in shells
        ]
    )
    gamma = float(np.dot(source, coupled) / np.dot(source, source))
    residual = coupled - gamma * source
    return float(0.5 * np.dot(residual, residual) / np.dot(coupled, coupled))


def finite_parameter_gradient(name: str, parameters: dict[str, float]) -> dict:
    gradient = {}
    for key, base in parameters.items():
        step = 1e-5 * max(1.0, abs(base))
        plus = dict(parameters)
        minus = dict(parameters)
        plus[key] += step
        minus[key] -= step
        gradient[key] = (
            finite_cgeo_loss(name, plus) - finite_cgeo_loss(name, minus)
        ) / (2.0 * step)
    return gradient


def program95_adaptive_gradient_and_quotient_audit() -> dict:
    n = 12
    rng = np.random.default_rng(SEED)
    k = project_symmetric_circulant_zero_diagonal(rng.normal(size=(n, n)))
    c = radial_matrix(k_strict, n)
    c = project_symmetric_circulant_zero_diagonal(c)
    a, b = 0.5, 0.1
    gradient = a * k + b * (k @ k @ k) - c
    projected = project_symmetric_circulant_zero_diagonal(gradient)
    velocity = -projected
    analytic = -float(np.sum(projected * projected))
    eps = 1e-7
    finite = (
        moment_functional(k + eps * velocity, c, a, b)
        - moment_functional(k - eps * velocity, c, a, b)
    ) / (2.0 * eps)
    gauge_shifts = {
        str(shift): moment_functional(k + shift * np.eye(n), c, a, b)
        - moment_functional(k, c, a, b)
        for shift in [-1.0, -0.25, 0.25, 1.0]
    }

    parameter_sets = {
        "legacy": {
            "alpha_geo": ALPHA_GEO,
            "omega": OMEGA_L,
            "phi": PHI_L,
            "beta_tors": BETA_TORS,
        },
        "strict": {
            "omega": OMEGA_S,
            "phi": PHI_S,
            "beta": BETA_S,
            "eta": ETA_S,
        },
    }
    p2772_rows = []
    for name, params in parameter_sets.items():
        grad = finite_parameter_gradient(name, params)
        norm = math.sqrt(sum(value * value for value in grad.values()))
        lr = 1e-3
        updated = {key: params[key] - lr * grad[key] for key in params}
        p2772_rows.append(
            {
                "kernel": name,
                "initial_loss": finite_cgeo_loss(name, params),
                "gradient": grad,
                "gradient_norm": norm,
                "one_step_loss": finite_cgeo_loss(name, updated),
                "stationary_at_1e_minus_9": norm < 1e-9,
            }
        )

    p2772_path = (
        ROOT
        / "fundamental_action_reconstruction/generated/"
        "p2772_s1722_self_learning_kernel_update_law_stationarity_witness.json"
    )
    return {
        "operator_space_identity": {
            "functional": "F(K)=a/2 tr(K^2)+b/4 tr(K^4)-tr(KC)",
            "flow": "dot K=-Pi_S(aK+bK^3-C)",
            "analytic_dF_dt": analytic,
            "finite_difference_dF_dt": finite,
            "identity_residual": abs(analytic - finite),
            "all_tested_gauge_shift_defects": gauge_shifts,
        },
        "p2772_reproduction": p2772_rows,
        "p2772_source_sha256": sha256(p2772_path),
        "theorem": (
            "The stated adaptive law is exactly a gradient flow after choosing "
            "an admissible linear subspace S and its orthogonal projection. It "
            "does not descend to the quotient K~K+cI unless F is scalar-shift "
            "invariant. The explicit nonzero gauge-shift defects disprove that "
            "invariance for the tested nonlinear moment functional."
        ),
        "status": "Exact projected-gradient identity; quotient claim refuted",
        "scope": (
            "The P2772 C_geo loss is a computational candidate, not an "
            "ontologically sourced FIN learning law; neither frozen kernel "
            "tuple is stationary for it."
        ),
    }


def program96_strict_chiral_source_intake() -> dict:
    base = ROOT / "fundamental_action_reconstruction/generated"
    names = [
        "p2749_s1699_minimal_inversion_odd_source_coupling_polarity_audit.json",
        "p2750_s1700_concrete_odd_source_sign_value_inventory_no_go.json",
        "p2759_s1709_post_p2758_no_new_live_frontier_reconciliation.json",
        "p2777_s1727_symmetry_source_selector_geometry_audit.json",
    ]
    rows = []
    for name in names:
        path = base / name
        payload = read_json(path)
        rows.append(
            {
                "artifact": name,
                "sha256": sha256(path),
                "status": payload.get("status", "missing"),
                "exists": path.exists(),
            }
        )
    all_present = all(row["exists"] for row in rows)
    return {
        "audited_post_p2721_artifacts": rows,
        "all_declared_inputs_present": all_present,
        "accepted_new_strict_signed_source_count": 0,
        "required_acceptance_tuple": [
            "explicit strict formula",
            "computable nonzero signed value",
            "non-premise provenance",
            "one selected P2721 coupling polarity",
            "robustness under the declared automorphism/gauge action",
        ],
        "decision": (
            "The new artifacts sharpen the representation and geometry "
            "boundaries but export no object satisfying all five criteria. "
            "The paired chiral state construction therefore remains "
            "conditional and QW-2191 remains open."
        ),
        "status": "Provenance revalidation; no admitted source",
    }


def js_divergence(p: np.ndarray, q: np.ndarray) -> float:
    p = np.asarray(p, dtype=float).ravel()
    q = np.asarray(q, dtype=float).ravel()
    p /= p.sum()
    q /= q.sum()
    m = 0.5 * (p + q)
    total = 0.0
    for distribution in [p, q]:
        keep = distribution > 0.0
        total += 0.5 * float(
            np.sum(distribution[keep] * np.log(distribution[keep] / m[keep]))
        )
    return total


def chernoff_information(p: np.ndarray, q: np.ndarray) -> tuple[float, float]:
    p = np.maximum(np.asarray(p, dtype=float).ravel(), 1e-300)
    q = np.maximum(np.asarray(q, dtype=float).ravel(), 1e-300)
    p /= p.sum()
    q /= q.sum()
    grid = np.linspace(0.0, 1.0, 1001)
    overlaps = np.array(
        [np.sum(np.exp(s * np.log(p) + (1.0 - s) * np.log(q))) for s in grid]
    )
    index = int(np.argmin(overlaps))
    return float(-math.log(overlaps[index])), float(grid[index])


def coarse_map(kind: str, n: int) -> np.ndarray:
    if kind == "full_sites":
        return np.eye(n)
    if kind == "adjacent_pairs":
        result = np.zeros((n // 2, n))
        for i in range(n):
            result[i // 2, i] = 1.0
        return result
    if kind == "parity":
        result = np.zeros((2, n))
        for i in range(n):
            result[i % 2, i] = 1.0
        return result
    raise ValueError(kind)


def confusion_matrix(size: int, error: float) -> np.ndarray:
    if size == 1:
        return np.ones((1, 1))
    return (1.0 - error) * np.eye(size) + error / (size - 1) * (
        np.ones((size, size)) - np.eye(size)
    )


def two_time_joint(l: np.ndarray, half_time: float) -> tuple[np.ndarray, np.ndarray]:
    n = len(l)
    u = expm(-1j * half_time * l)
    p = expm(-half_time * l)
    wave = np.zeros((n, n))
    diffusion = np.zeros((n, n))
    for y in range(n):
        wave[y, :] = abs(u[y, 0]) ** 2 * abs(u[:, y]) ** 2
        diffusion[y, :] = p[y, 0] * p[:, y]
    return wave, diffusion


def observed_joint(
    joint: np.ndarray, instrument: str, detector_error: float
) -> np.ndarray:
    g = coarse_map(instrument, joint.shape[0])
    coarse = g @ joint @ g.T
    c = confusion_matrix(coarse.shape[0], detector_error)
    observed = c @ coarse @ c.T
    observed /= observed.sum()
    return observed


def monte_carlo_likelihood_error(
    p: np.ndarray, q: np.ndarray, shots: int, replicates: int, rng
) -> float:
    p = np.maximum(p.ravel(), 1e-300)
    q = np.maximum(q.ravel(), 1e-300)
    p /= p.sum()
    q /= q.sum()
    log_ratio = np.log(p / q)
    error_wave = 0
    error_diff = 0
    for _ in range(replicates):
        count_p = rng.multinomial(shots, p)
        count_q = rng.multinomial(shots, q)
        error_wave += float(np.dot(count_p, log_ratio)) < 0.0
        error_diff += float(np.dot(count_q, log_ratio)) >= 0.0
    return float((error_wave + error_diff) / (2.0 * replicates))


def program97_noisy_process_tensor_design() -> dict:
    n = 12
    w = radial_matrix(k_strict, n, absolute=True)
    w /= w.sum(axis=1)[0]
    l = laplacian(w)
    times = np.geomspace(0.025, 2.0, 80)
    errors = [0.0, 0.05, 0.10]
    instruments = ["full_sites", "adjacent_pairs", "parity"]
    rows = []
    for instrument in instruments:
        for half_time in times:
            wave, diffusion = two_time_joint(l, float(half_time))
            js_values = []
            for error in errors:
                pw = observed_joint(wave, instrument, error)
                pd = observed_joint(diffusion, instrument, error)
                js_values.append(js_divergence(pw, pd))
            rows.append(
                {
                    "instrument": instrument,
                    "half_time": float(half_time),
                    "js_by_detector_error": dict(
                        zip([str(e) for e in errors], js_values)
                    ),
                    "worst_case_js": min(js_values),
                }
            )
    best = max(rows, key=lambda row: row["worst_case_js"])
    wave, diffusion = two_time_joint(l, best["half_time"])
    p_best = observed_joint(wave, best["instrument"], 0.10)
    q_best = observed_joint(diffusion, best["instrument"], 0.10)
    chernoff, s_star = chernoff_information(p_best, q_best)
    rng = np.random.default_rng(SEED)
    monte_carlo = [
        {
            "shots": shots,
            "equal_prior_likelihood_error": monte_carlo_likelihood_error(
                p_best, q_best, shots, 2000, rng
            ),
        }
        for shots in [20, 50, 100, 250]
    ]

    fig, ax = plt.subplots(figsize=(8.8, 4.8), constrained_layout=True)
    for instrument in instruments:
        selected = [row for row in rows if row["instrument"] == instrument]
        ax.semilogx(
            [row["half_time"] for row in selected],
            [row["worst_case_js"] for row in selected],
            label=instrument.replace("_", " "),
        )
    ax.axvline(best["half_time"], color="black", ls="--", alpha=0.6)
    ax.set_xlabel("half-time (dimensionless)")
    ax.set_ylabel("worst-case Jensen-Shannon divergence")
    ax.set_title("Robust two-time instrument design under detector confusion")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend()
    fig.savefig(FIG / "program97_process_tensor_design.png", dpi=190)
    plt.close(fig)

    return {
        "search_size": len(rows),
        "detector_error_grid": errors,
        "best_maximin_design": best,
        "chernoff_information_at_error_0_10": chernoff,
        "chernoff_s_star": s_star,
        "monte_carlo_at_error_0_10": monte_carlo,
        "best_wave_probabilities": p_best.tolist(),
        "best_diffusion_probabilities": q_best.tolist(),
        "top_five_designs": sorted(
            rows, key=lambda row: row["worst_case_js"], reverse=True
        )[:5],
        "status": "Conditional optimized finite operational experiment",
        "scope": (
            "The preparation, projective intervention, detector confusion, "
            "dimensionless clock, and decision rule are imported. The result "
            "distinguishes two operator functions, not a physical FIN signal."
        ),
    }


def binary_entropy(error: float) -> float:
    if error <= 0.0 or error >= 1.0:
        return 0.0
    return -error * math.log(error) - (1.0 - error) * math.log(1.0 - error)


def complete_feedback_row(error: float) -> dict:
    delta = math.log((1.0 - error) / error)
    conditional = binary_entropy(error)
    mutual = math.log(2.0) - conditional
    delta_f = math.log(2.0) - math.log1p(math.exp(-delta))
    system_work = error * delta
    reset = math.log(2.0)
    total = system_work + reset
    return {
        "measurement_error": error,
        "mutual_information": mutual,
        "conditional_entropy": conditional,
        "system_feedback_work": system_work,
        "free_energy_difference": delta_f,
        "system_equality_residual": abs(system_work - (delta_f - mutual)),
        "unconditional_memory_reset_cost": reset,
        "apparatus_inclusive_work": total,
        "excess_over_deltaF": total - delta_f,
        "excess_equals_conditional_entropy_residual": abs(
            total - delta_f - conditional
        ),
    }


def program98_apparatus_inclusive_feedback() -> dict:
    rows = [
        complete_feedback_row(error)
        for error in [0.001, 0.01, 0.05, 0.10, 0.20, 0.40, 0.49]
    ]
    fig, ax = plt.subplots(figsize=(8.5, 4.7), constrained_layout=True)
    ax.plot(
        [row["measurement_error"] for row in rows],
        [row["system_feedback_work"] for row in rows],
        "o-",
        label="system feedback work",
    )
    ax.plot(
        [row["measurement_error"] for row in rows],
        [row["apparatus_inclusive_work"] for row in rows],
        "s-",
        label="plus memory reset",
    )
    ax.plot(
        [row["measurement_error"] for row in rows],
        [row["free_energy_difference"] for row in rows],
        "--",
        label="Delta F baseline",
    )
    ax.set_xlabel("binary measurement error")
    ax.set_ylabel("dimensionless work / information")
    ax.set_title("The complete memory ledger removes apparent free work")
    ax.grid(True, alpha=0.25)
    ax.legend()
    fig.savefig(FIG / "program98_feedback_ledger.png", dpi=190)
    plt.close(fig)
    return {
        "rows": rows,
        "exact_identity": (
            "W_system=DeltaF-I and W_reset=H(Y)=ln2 imply "
            "W_complete-DeltaF=ln2-I=H(Y|X)>=0."
        ),
        "status": "Exact apparatus-inclusive information/work identity",
        "scope": (
            "beta=1, the Hamiltonian protocol, bath, measurement, controller, "
            "and erasure convention are external operational inputs. FIN does "
            "not derive their dimensional energy scale."
        ),
    }


def write_external_intake_template() -> dict:
    record = {
        "schema": "FIN Programs 91-100 independent dataset intake v1",
        "required_fields": {
            "dataset_id": None,
            "immutable_source_uri": None,
            "retrieval_timestamp_utc": None,
            "license": None,
            "raw_sha256": None,
            "calibration_record": None,
            "instrument_definition": None,
            "preparation_definition": None,
            "clock_and_units": None,
            "pre_registered_primary_statistic": (
                "two-slot wave-versus-diffusion log-likelihood ratio"
            ),
            "pre_registered_exclusion_rule": None,
            "blind_holdout_identifier": None,
        },
        "admission_rule": (
            "No evidential claim before immutable raw bytes, provenance, "
            "calibration, declared instrument, and a blind holdout are present."
        ),
        "current_admitted_datasets": [],
    }
    canonical = json.dumps(
        record, sort_keys=True, separators=(",", ":")
    ).encode("utf-8")
    record["canonical_core_sha256"] = hashlib.sha256(canonical).hexdigest()
    INTAKE.write_text(
        json.dumps(record, indent=2, ensure_ascii=False, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return record


def program99_external_acquisition_design(program97: dict) -> dict:
    chernoff = program97["chernoff_information_at_error_0_10"]
    requirements = []
    for alpha in [0.05, 0.01, 0.001, 0.0001]:
        shots = math.ceil(math.log(0.5 / alpha) / chernoff)
        requirements.append(
            {
                "target_equal_prior_error_upper_bound": alpha,
                "chernoff_sufficient_shots": shots,
            }
        )
    intake = write_external_intake_template()
    return {
        "best_design_from_program97": program97["best_maximin_design"],
        "chernoff_information": chernoff,
        "sample_size_bounds": requirements,
        "intake_template": INTAKE.name,
        "intake_digest": intake["canonical_core_sha256"],
        "admitted_external_datasets": [],
        "status": "Acquisition and power design completed; no external data admitted",
        "scope": (
            "Chernoff counts are sufficient model-based bounds under i.i.d. "
            "records and the declared detector model. They are not empirical "
            "power until an independent calibrated dataset is admitted."
        ),
    }


def program100_damping_completion_atom() -> dict:
    d = np.unique(
        np.rint(np.geomspace(1.0, 1_000_000.0, 1000)).astype(int)
    ).astype(float)
    legacy_envelope = 1.0 / (1.0 + BETA_TORS * d)
    strict_envelope = 1.0 / (1.0 + BETA_S * d**ETA_S)
    completion = strict_envelope / legacy_envelope
    reconstructed = legacy_envelope * completion
    residual = float(
        np.max(np.abs(reconstructed - strict_envelope))
        / np.max(np.abs(strict_envelope))
    )
    tail = d >= 10_000
    completion_slope = float(
        np.polyfit(np.log(d[tail]), np.log(completion[tail]), 1)[0]
    )
    legacy_slope = float(
        np.polyfit(np.log(d[tail]), np.log(legacy_envelope[tail]), 1)[0]
    )
    strict_slope = float(
        np.polyfit(np.log(d[tail]), np.log(strict_envelope[tail]), 1)[0]
    )

    gap_files = [
        "p2760_s1710_foundation_kernel_lagrangian_gap_matrix.json",
        "p2761_s1711_kernel_moment_coupling_provenance_obstruction.json",
        "p2766_s1716_post_moment_provenance_state_reconciliation.json",
    ]
    gap_rows = []
    base = ROOT / "fundamental_action_reconstruction/generated"
    for name in gap_files:
        path = base / name
        payload = read_json(path)
        gap_rows.append(
            {
                "artifact": name,
                "status": payload.get("status", "missing"),
                "sha256": sha256(path),
            }
        )

    fig, ax = plt.subplots(figsize=(8.5, 4.8), constrained_layout=True)
    ax.loglog(d, legacy_envelope, label="legacy damping envelope")
    ax.loglog(d, strict_envelope, label="strict damping envelope")
    ax.loglog(d, completion, "--", label="necessary completion multiplier")
    ax.set_xlabel("distance d")
    ax.set_ylabel("amplitude")
    ax.set_title("The legacy-to-strict damping bridge requires a d^-0.8 atom")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend()
    fig.savefig(FIG / "program100_damping_completion_atom.png", dpi=190)
    plt.close(fig)

    return {
        "exact_typed_atom": (
            "C_damp(d)=(1+beta_tors d)/(1+beta d^(9/5)); "
            "D_strict(d)=C_damp(d) D_legacy(d)"
        ),
        "relative_reconstruction_residual": residual,
        "tail_slopes": {
            "legacy_envelope": legacy_slope,
            "strict_envelope": strict_slope,
            "completion_multiplier": completion_slope,
            "predicted_completion": 1.0 - ETA_S,
        },
        "provenance_gap_inputs": gap_rows,
        "theorem": (
            "A constant amplitude change and a linear distance rescaling "
            "preserve the power-law exponent, so they cannot map a d^-1 "
            "legacy envelope to a d^-9/5 strict envelope. Any multiplicative "
            "completion must have asymptotic degree -4/5. The displayed atom "
            "is the unique pointwise positive multiplier for the two frozen "
            "denominators."
        ),
        "status": "Necessary bridge atom constructed; source theorem absent",
        "guardrail": (
            "C_damp is target-coded from both kernels. It does not source "
            "beta=1 or eta=9/5, does not bridge the phase/frequency or "
            "amplitude layers, and licenses no legacy physical-role transfer."
        ),
    }


def main() -> dict:
    results = {
        "metadata": {
            "release": "10.10",
            "executed_programs": list(range(91, 101)),
            "date": "2026-07-27",
            "seed": SEED,
            "creator": "Żuchowski, Krzysztof",
            "orcid": "0009-0002-0909-3613",
        },
        "program91": program91_strict_nonclosure_certificate(),
        "program92": program92_critical_tuning_and_field_rescaling(),
        "program93": program93_graphon_continuum_limit(),
        "program94": program94_long_range_propagation(),
        "program95": program95_adaptive_gradient_and_quotient_audit(),
        "program96": program96_strict_chiral_source_intake(),
    }
    results["program97"] = program97_noisy_process_tensor_design()
    results["program98"] = program98_apparatus_inclusive_feedback()
    results["program99"] = program99_external_acquisition_design(
        results["program97"]
    )
    results["program100"] = program100_damping_completion_atom()
    results["global_guardrail"] = (
        "No result discharges QW-2191, sources a dimensional standard, proves "
        "Lorentz invariance or local spacetime, identifies information with "
        "physical entropy without a conversion package, completes the "
        "legacy-to-strict bridge, transfers legacy physical roles, promotes "
        "role-bearing L_total, or constitutes external physical evidence. "
        "The nadsoliton remains the primordial information in a solitonic "
        "state, with no separate informational layer beneath it."
    )
    OUT.write_text(
        json.dumps(results, indent=2, ensure_ascii=False, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(f"Wrote {OUT.name}")
    for key in [f"program{i}" for i in range(91, 101)]:
        print(key, ":", results[key]["status"])
    return results


if __name__ == "__main__":
    main()
