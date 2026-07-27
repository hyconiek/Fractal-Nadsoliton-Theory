#!/usr/bin/env python3
"""Execute FIN Programs 61--70.

The suite tests continuum/projective compatibility, regularizer independence,
Schur-compression composition, signed/Krein alternatives, operational
tomography, causal-order compatibility, conversion-scale identifiability,
one explicit chiral source law, a Landauer protocol, and a blinded
held-out model challenge.

All FIN-core inputs are dimensionless.  Imported operational, bath, unit, or
sector data are marked as conditioning data and are not strict exports.
"""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.linalg import expm
from scipy.optimize import minimize_scalar


ROOT = Path(__file__).resolve().parent
OUT = ROOT / "FIN_Programs_61_70_Continuum_Operational_Physics_Results.json"
FIG = ROOT / "FIN_Programs_61_70_Continuum_Operational_Physics_Figures"
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
    return ALPHA_GEO * np.cos(OMEGA_L * d + PHI_L) / (1.0 + BETA_TORS * d)


def k_strict(d):
    d = np.asarray(d, dtype=float)
    return np.cos(OMEGA_S * d + PHI_S) / (1.0 + BETA_S * d**ETA_S)


def cycle_distance(i: int, j: int, n: int) -> int:
    return min((i - j) % n, (j - i) % n)


def radial_matrix(kernel, n: int) -> np.ndarray:
    w = np.zeros((n, n), dtype=float)
    for i in range(n):
        for j in range(n):
            if i != j:
                w[i, j] = float(kernel(cycle_distance(i, j, n)))
    return w


def laplacian(w: np.ndarray) -> np.ndarray:
    return np.diag(w.sum(axis=1)) - w


def reflection(n: int) -> np.ndarray:
    r = np.zeros((n, n))
    for i in range(n):
        r[i, (-i) % n] = 1.0
    return r


def shift_to_spd(a: np.ndarray, floor: float = MASS2) -> tuple[np.ndarray, float]:
    mineig = float(np.linalg.eigvalsh(a).min())
    shift = max(0.0, floor - mineig)
    return a + shift * np.eye(a.shape[0]), shift


def schur(a: np.ndarray, keep: np.ndarray) -> np.ndarray:
    keep = np.asarray(keep, dtype=int)
    keep_set = set(map(int, keep))
    drop = np.array([i for i in range(a.shape[0]) if i not in keep_set], dtype=int)
    akk = a[np.ix_(keep, keep)]
    ako = a[np.ix_(keep, drop)]
    aoo = a[np.ix_(drop, drop)]
    return akk - ako @ np.linalg.solve(aoo, ako.conj().T)


def normalize_precision(a: np.ndarray) -> np.ndarray:
    scale = float(np.trace(a).real / a.shape[0])
    return a / scale


def positive_normalized_precision(kernel, n: int, repair: str) -> np.ndarray:
    w = radial_matrix(kernel, n)
    if repair == "absolute":
        w = np.abs(w)
    elif repair == "positive_part":
        w = np.maximum(w, 0.0)
    else:
        raise ValueError(repair)
    row_sum = float(w.sum(axis=1)[0])
    w = w / row_sum
    return MASS2 * np.eye(n) + laplacian(w)


def projective_defect(kernel, repair: str, n: int) -> dict:
    parent = positive_normalized_precision(kernel, 2 * n, repair)
    compressed = schur(parent, np.arange(0, 2 * n, 2))
    native = positive_normalized_precision(kernel, n, repair)
    cn = normalize_precision(compressed)
    nn = normalize_precision(native)
    operator_defect = float(np.linalg.norm(cn - nn, "fro") / np.linalg.norm(nn, "fro"))
    green_defect = float(
        np.linalg.norm(np.linalg.inv(cn) - np.linalg.inv(nn), "fro")
        / np.linalg.norm(np.linalg.inv(nn), "fro")
    )
    return {
        "n": n,
        "operator_projective_defect": operator_defect,
        "green_projective_defect": green_defect,
    }


def program61_continuum_functor() -> dict:
    sizes = [6, 12, 24, 48]
    families = {
        "strict_absolute": (k_strict, "absolute"),
        "legacy_absolute": (k_legacy, "absolute"),
    }
    rows = {}
    fig, axes = plt.subplots(1, 2, figsize=(10.2, 4.2), constrained_layout=True)
    for name, (kernel, repair) in families.items():
        data = [projective_defect(kernel, repair, n) for n in sizes]
        op_values = np.array([x["operator_projective_defect"] for x in data])
        green_values = np.array([x["green_projective_defect"] for x in data])
        rows[name] = {
            "rows": data,
            "operator_defect_monotone_decreasing": bool(
                np.all(np.diff(op_values) < 0.0)
            ),
            "green_defect_monotone_decreasing": bool(
                np.all(np.diff(green_values) < 0.0)
            ),
            "operator_loglog_slope": float(
                np.polyfit(np.log(sizes), np.log(op_values), 1)[0]
            ),
            "green_loglog_slope": float(
                np.polyfit(np.log(sizes), np.log(green_values), 1)[0]
            ),
        }
        axes[0].loglog(
            sizes,
            [x["operator_projective_defect"] for x in data],
            "o-",
            label=name.replace("_", " "),
        )
        axes[1].loglog(
            sizes,
            [x["green_projective_defect"] for x in data],
            "s-",
            label=name.replace("_", " "),
        )
    for ax, title in zip(
        axes,
        ["operator naturality defect", "Green naturality defect"],
    ):
        ax.set_xlabel("retained cycle size N")
        ax.set_ylabel("relative defect")
        ax.set_title(title)
        ax.legend(fontsize=8)
        ax.grid(True, which="both", alpha=0.25)
    fig.savefig(FIG / "program61_projective_continuum.png", dpi=190)
    plt.close(fig)
    return {
        "fixed_normalization": "mean diagonal equals one",
        "positive_realization": "absolute-value repair, row-normalized weights",
        "sizes": sizes,
        "families": rows,
        "verdict": (
            "Exact Schur compression exists at every finite level, but a "
            "continuum/projective functor additionally requires the compressed "
            "operator to agree with the native lower-resolution operator. "
            "The measured naturality defects test, rather than assume, that claim."
        ),
    }


def resistance(l: np.ndarray) -> np.ndarray:
    g = np.linalg.pinv(l, hermitian=True)
    d = np.diag(g)
    return np.maximum(d[:, None] + d[None, :] - 2.0 * g, 0.0)


def shifted_resistance(l: np.ndarray, delta: float) -> tuple[np.ndarray, float]:
    g = np.linalg.inv(l + delta * np.eye(l.shape[0]))
    d = np.diag(g)
    r = np.maximum(d[:, None] + d[None, :] - 2.0 * g, 0.0)
    return r, float(np.linalg.norm(g, "fro"))


def regularizer_rows(l: np.ndarray) -> list[dict]:
    target = resistance(l)
    rows = []
    for delta in np.logspace(-1, -7, 7):
        r, gnorm = shifted_resistance(l, float(delta))
        rows.append(
            {
                "delta": float(delta),
                "resistance_relative_error": float(
                    np.linalg.norm(r - target, "fro") / np.linalg.norm(target, "fro")
                ),
                "green_fro_norm": gnorm,
                "delta_times_green_fro_norm": float(delta * gnorm),
            }
        )
    return rows


def program62_regularizer_independence() -> dict:
    ls = laplacian(radial_matrix(k_strict, 12))
    ll = laplacian(np.abs(radial_matrix(k_legacy, 12)))
    strict_rows = regularizer_rows(ls)
    legacy_rows = regularizer_rows(ll)

    fig, ax = plt.subplots(figsize=(8.2, 4.7), constrained_layout=True)
    for name, rows in [("strict", strict_rows), ("legacy absolute", legacy_rows)]:
        ax.loglog(
            [x["delta"] for x in rows],
            [x["resistance_relative_error"] for x in rows],
            "o-",
            label=name,
        )
    ax.invert_xaxis()
    ax.set_xlabel(r"regularizer $\delta\to0^+$")
    ax.set_ylabel("relative resistance error")
    ax.set_title("Zero-mode divergence cancels in resistance differences")
    ax.legend()
    ax.grid(True, which="both", alpha=0.25)
    fig.savefig(FIG / "program62_regularizer_limit.png", dpi=190)
    plt.close(fig)
    return {
        "strict": strict_rows,
        "legacy_absolute": legacy_rows,
        "theorem": (
            "For a connected graph Laplacian, the divergent delta^-1 zero-mode "
            "term in (L+delta I)^-1 cancels from pairwise resistance. The shifted "
            "resistance converges to the Moore-Penrose resistance, while the full "
            "Green norm diverges."
        ),
        "scope": (
            "This limit is available for declared positive graph realizations. "
            "It does not make the raw signed legacy Laplacian a Markov generator."
        ),
    }


def program63_compression_semigroup() -> dict:
    a48 = positive_normalized_precision(k_legacy, 48, "absolute")
    direct = schur(a48, np.arange(0, 48, 4))
    first = schur(a48, np.arange(0, 48, 2))
    sequential = schur(first, np.arange(0, 24, 2))
    composition_residual = float(
        np.linalg.norm(direct - sequential, "fro") / np.linalg.norm(direct, "fro")
    )

    defects = [projective_defect(k_legacy, "absolute", n) for n in [6, 12, 24, 48]]
    return {
        "direct_48_to_12_vs_sequential_48_to_24_to_12_relative_residual": composition_residual,
        "semigroup_theorem": (
            "Nested Schur elimination is transitive: if E is contained in F, "
            "Schur_E(Schur_F(A))=Schur_E(A), provided eliminated blocks are invertible."
        ),
        "legacy_projective_defects": defects,
        "fixed_point_status": (
            "Exact composition does not imply a fixed point. A fixed point also "
            "requires equality with the independently defined native lower-scale "
            "operator under one normalization."
        ),
    }


def program64_krein_signed_alternative() -> dict:
    n = 12
    l_signed = laplacian(radial_matrix(k_legacy, n))
    a = l_signed + 0.05 * np.eye(n)
    eig, vec = np.linalg.eigh(a)
    signs = np.sign(eig)
    j = vec @ np.diag(signs) @ vec.T
    abs_a = j @ a
    t = 0.1
    u = expm(-1j * t * a)
    diffusion = expm(-t * a)
    stable = expm(-t * abs_a)
    g = np.linalg.inv(1j * np.eye(n) - a)
    g0 = np.linalg.inv(a)
    green_d = np.diag(g0)[:, None] + np.diag(g0)[None, :] - 2 * g0
    return {
        "inertia": {
            "negative": int(np.sum(eig < -1e-10)),
            "positive": int(np.sum(eig > 1e-10)),
            "zero": int(np.sum(abs(eig) <= 1e-10)),
        },
        "minimum_eigenvalue": float(eig.min()),
        "maximum_eigenvalue": float(eig.max()),
        "J_involution_residual": float(np.linalg.norm(j @ j - np.eye(n), "fro")),
        "Krein_selfadjoint_residual": float(np.linalg.norm(a.T @ j - j @ a, "fro")),
        "JA_minimum_eigenvalue": float(np.linalg.eigvalsh(abs_a).min()),
        "unitarity_residual": float(np.linalg.norm(u.conj().T @ u - np.eye(n), "fro")),
        "raw_diffusion_operator_norm": float(np.linalg.norm(diffusion, 2)),
        "absolute_generator_diffusion_norm": float(np.linalg.norm(stable, 2)),
        "raw_diffusion_minimum_entry": float(diffusion.min()),
        "nonreal_resolvent_inverse_residual": float(
            np.linalg.norm((1j * np.eye(n) - a) @ g - np.eye(n), "fro")
        ),
        "indefinite_green_distance_negative_entries": int(np.sum(green_d < -1e-10)),
        "verdict": (
            "The signed Hermitian generator already gives unitary dynamics and "
            "off-spectrum resolvents. A Krein fundamental symmetry converts it "
            "to |A| for a positive diffusion, but that changes the generator. "
            "Krein structure does not export a Markov process or physical metric "
            "from the raw signed legacy operator."
        ),
    }


def channel_records(l: np.ndarray, t: float) -> tuple[np.ndarray, np.ndarray]:
    u = expm(-1j * t * l)
    q = np.abs(u) ** 2
    p = expm(-t * l)
    return q, p


def circulant_channel_residual(c: np.ndarray) -> float:
    base = c[:, 0]
    return float(max(np.linalg.norm(c[:, j] - np.roll(base, j)) for j in range(c.shape[1])))


def program65_operational_tomography() -> dict:
    n = 12
    l = laplacian(radial_matrix(k_strict, n))
    t = 1.0
    q, p = channel_records(l, t)
    d = q - p
    column = d[:, 0]
    tv = float(0.5 * np.sum(np.abs(column)))
    optimal_subset = np.flatnonzero(column > 0.0)
    binary_difference = float(abs(column[optimal_subset].sum()))
    coarse = np.zeros((n // 2, n))
    for i in range(n):
        coarse[i // 2, i] = 1.0
    tv_pair = float(0.5 * np.sum(np.abs(coarse @ column)))

    times = np.linspace(0.0, 3.0, 241)
    tvs = []
    pair_tvs = []
    for ti in times:
        qi, pi = channel_records(l, float(ti))
        col = qi[:, 0] - pi[:, 0]
        tvs.append(0.5 * np.sum(abs(col)))
        pair_tvs.append(0.5 * np.sum(abs(coarse @ col)))

    fig, ax = plt.subplots(figsize=(8.2, 4.7), constrained_layout=True)
    ax.plot(times, tvs, label="full site instrument")
    ax.plot(times, pair_tvs, label="pair-coarse instrument")
    ax.set_xlabel("dimensionless time")
    ax.set_ylabel("single-preparation total variation")
    ax.set_title("Operational distinguishability of unitary and diffusive calculi")
    ax.legend()
    fig.savefig(FIG / "program65_operational_tomography.png", dpi=190)
    plt.close(fig)
    return {
        "dimensionless_time": t,
        "channel_difference_rank": int(np.linalg.matrix_rank(d, tol=1e-10)),
        "unitary_circulant_residual": circulant_channel_residual(q),
        "diffusive_circulant_residual": circulant_channel_residual(p),
        "one_localized_preparation_full_detector_TV": tv,
        "optimal_binary_detector_sites": optimal_subset.tolist(),
        "one_preparation_optimal_binary_difference": binary_difference,
        "one_preparation_pair_coarse_TV": tv_pair,
        "minimum_exact_discrimination_under_known_model_class": (
            "one localized preparation plus one nontrivial binary instrument at t>0"
        ),
        "tomography_scope": (
            "One preparation reconstructs a translation-covariant circulant "
            "transition kernel. Without the circulant prior, full channel "
            "tomography requires a spanning preparation set."
        ),
    }


def nearest_neighbor_laplacian(n: int, weight: float) -> np.ndarray:
    w = np.zeros((n, n))
    for i in range(n):
        w[i, (i + 1) % n] = weight
        w[i, (i - 1) % n] = weight
    return laplacian(w)


def first_nonzero_power(a: np.ndarray, source: int, target: int) -> tuple[int, float]:
    power = np.eye(a.shape[0])
    for m in range(1, a.shape[0] + 1):
        power = power @ a
        value = float(power[target, source])
        if abs(value) > 1e-12:
            return m, value
    raise RuntimeError("no connecting power")


def program66_causal_order() -> dict:
    n = 12
    source, target = 0, 6
    full = laplacian(radial_matrix(k_strict, n))
    local = nearest_neighbor_laplacian(n, float(k_strict(1)))
    m_full, coeff_full = first_nonzero_power(full, source, target)
    m_local, coeff_local = first_nonzero_power(local, source, target)
    times = np.logspace(-2, 0, 160)
    curves = {}
    for name, a in [("strict_all_distance", full), ("nearest_neighbor_control", local)]:
        wave = []
        diff = []
        for t in times:
            wave.append(abs(expm(-1j * t * a)[target, source]) ** 2)
            diff.append(max(float(expm(-t * a)[target, source]), 1e-300))
        curves[name] = {"wave": wave, "diffusion": diff}

    fig, axes = plt.subplots(1, 2, figsize=(10.2, 4.2), constrained_layout=True)
    for name, vals in curves.items():
        axes[0].loglog(times, np.maximum(vals["wave"], 1e-300), label=name.replace("_", " "))
        axes[1].loglog(times, vals["diffusion"], label=name.replace("_", " "))
    axes[0].set_title("far-site wave probability")
    axes[1].set_title("far-site diffusion probability")
    for ax in axes:
        ax.set_xlabel("dimensionless time")
        ax.set_ylabel("response at opposite site")
        ax.legend(fontsize=8)
        ax.grid(True, which="both", alpha=0.25)
    fig.savefig(FIG / "program66_causal_tails.png", dpi=190)
    plt.close(fig)
    return {
        "cycle_distance": 6,
        "strict_all_distance_first_nonzero_power": m_full,
        "strict_wave_leading_probability_exponent": 2 * m_full,
        "strict_diffusion_leading_exponent": m_full,
        "nearest_neighbor_first_nonzero_power": m_local,
        "nearest_neighbor_wave_leading_probability_exponent": 2 * m_local,
        "nearest_neighbor_diffusion_leading_exponent": m_local,
        "strict_leading_matrix_coefficient": coeff_full,
        "nearest_neighbor_leading_matrix_coefficient": coeff_local,
        "verdict": (
            "The full strict C12 kernel couples the opposite site at first order, "
            "so it has no exact finite-speed cone relative to cycle distance. "
            "The local control is suppressed to graph-distance order but matrix "
            "exponentials still have nonzero tails for every t>0. A Lorentzian "
            "causal order is an additional structure."
        ),
    }


def program67_calibration_identifiability() -> dict:
    # Log-parameter order: log ell, log tau, log hbar.
    j_length = np.array([[1.0, 0.0, 0.0]])
    j_time = np.array([[0.0, 1.0, 0.0]])
    j_energy = np.array([[0.0, -1.0, 1.0]])
    full_j = np.vstack([j_length, j_time, j_energy])
    true = np.log(np.array([2.5e-6, 4.0e-9, 1.7e-34]))
    d_hat, theta_hat, lambda_hat = 3.0, 1.25, 0.7541211542070796
    obs = np.array(
        [
            true[0] + math.log(d_hat),
            true[1] + math.log(theta_hat),
            true[2] - true[1] + math.log(lambda_hat),
        ]
    )
    offsets = np.log([d_hat, theta_hat, lambda_hat])
    offsets[2] = math.log(lambda_hat)
    rhs = np.array([obs[0] - offsets[0], obs[1] - offsets[1], obs[2] - offsets[2]])
    recovered = np.linalg.solve(full_j, rhs)
    return {
        "rank_length_only": int(np.linalg.matrix_rank(j_length)),
        "rank_length_plus_time": int(np.linalg.matrix_rank(np.vstack([j_length, j_time]))),
        "rank_length_plus_energy": int(np.linalg.matrix_rank(np.vstack([j_length, j_energy]))),
        "rank_time_plus_energy": int(np.linalg.matrix_rank(np.vstack([j_time, j_energy]))),
        "rank_length_time_energy": int(np.linalg.matrix_rank(full_j)),
        "true_conversion_triple": np.exp(true).tolist(),
        "recovered_conversion_triple": np.exp(recovered).tolist(),
        "maximum_relative_recovery_error": float(
            np.max(np.abs(np.exp(recovered) / np.exp(true) - 1.0))
        ),
        "minimal_experiment_classes": [
            "one calibrated length record",
            "one calibrated clock record",
            "one calibrated energy/action record",
        ],
        "verdict": (
            "For independent length, time, and action units the log-Jacobian has "
            "rank three only when all three experiment classes are present. "
            "More dimensionless FIN observations cannot replace a missing "
            "calibration direction."
        ),
    }


def chiral_loop(d: np.ndarray) -> float:
    n = d.shape[0]
    value = 0.0
    for i in range(n):
        value += (
            d[i, (i + 1) % n]
            * d[(i + 1) % n, (i + 2) % n]
            * d[(i + 2) % n, i]
        )
        value -= (
            d[i, (i - 1) % n]
            * d[(i - 1) % n, (i - 2) % n]
            * d[(i - 2) % n, i]
        )
    return float(value)


def program68_chiral_source_law() -> dict:
    n = 12
    r = reflection(n)
    translation = np.roll(np.eye(n), 1, axis=0)
    b = np.zeros((n, n))
    for i in range(n):
        b[i, (i + 1) % n] = 1.0
        b[(i + 1) % n, i] = -1.0
    rows = {}
    for name, w in [
        ("strict", radial_matrix(k_strict, n)),
        ("legacy", radial_matrix(k_legacy, n)),
    ]:
        epsilon = 0.1
        dp = w + epsilon * b
        dm = w - epsilon * b
        rows[name] = {
            "radial_candidate": chiral_loop(w),
            "directed_plus": chiral_loop(dp),
            "directed_minus": chiral_loop(dm),
            "sign_pair_residual": abs(chiral_loop(dp) + chiral_loop(dm)),
            "reflection_odd_residual": abs(chiral_loop(r @ dp @ r.T) + chiral_loop(dp)),
            "translation_invariance_residual": abs(
                chiral_loop(translation @ dp @ translation.T) - chiral_loop(dp)
            ),
        }
    return {
        "candidate_formula": (
            "sum_i[D_i,i+1 D_i+1,i+2 D_i+2,i - "
            "D_i,i-1 D_i-1,i-2 D_i-2,i]"
        ),
        "results": rows,
        "general_no_go": (
            "If F(RWR)=-F(W) and the radial kernel obeys RWR=W, then "
            "F(W)=-F(W), hence F(W)=0."
        ),
        "verdict": (
            "The explicit loop is translation invariant and reflection odd, but "
            "vanishes on both radial kernels. It becomes nonzero only after an "
            "oriented antisymmetric datum is inserted; its sign then follows the "
            "inserted orientation. This is a receiver, not a strict source law."
        ),
    }


def binary_entropy(p: float) -> float:
    if p <= 0.0 or p >= 1.0:
        return 0.0
    return -p * math.log(p) - (1.0 - p) * math.log(1.0 - p)


def program69_landauer_protocol() -> dict:
    rows = []
    xs = [0.5, 1.0, 2.0, 5.0, 10.0, 20.0]
    for x in xs:
        error = 1.0 / (1.0 + math.exp(x))
        sf = binary_entropy(error)
        beta_work = math.log(2.0) - math.log1p(math.exp(-x))
        beta_delta_u = x * error
        beta_heat_bath = math.log(2.0) - sf
        rows.append(
            {
                "beta_gap": x,
                "reset_error": error,
                "final_Shannon_entropy_nats": sf,
                "beta_reversible_work": beta_work,
                "beta_heat_to_bath": beta_heat_bath,
                "beta_internal_energy_change": beta_delta_u,
                "first_law_residual": beta_work - beta_heat_bath - beta_delta_u,
                "Clausius_equality_residual": beta_heat_bath - (math.log(2.0) - sf),
            }
        )

    fig, ax1 = plt.subplots(figsize=(8.2, 4.7), constrained_layout=True)
    ax1.plot(xs, [x["beta_reversible_work"] for x in rows], "o-", label=r"$\beta W_{\rm rev}$")
    ax1.plot(xs, [x["beta_heat_to_bath"] for x in rows], "s-", label=r"$\beta Q_{\rm bath}$")
    ax1.axhline(math.log(2.0), color="black", linestyle="--", label=r"$\ln2$")
    ax1.set_xlabel(r"dimensionless final gap $\beta\Delta$")
    ax1.set_ylabel("dimensionless work/heat")
    ax2 = ax1.twinx()
    ax2.semilogy(xs, [x["reset_error"] for x in rows], "^-", color="tab:red", label="reset error")
    ax2.set_ylabel("reset error")
    lines = ax1.get_lines() + ax2.get_lines()
    ax1.legend(lines, [line.get_label() for line in lines], fontsize=8, loc="center right")
    ax1.set_title("Explicit reversible erasure protocol")
    fig.savefig(FIG / "program69_landauer_protocol.png", dpi=190)
    plt.close(fig)
    return {
        "conditioning_data": [
            "two-state memory",
            "thermal bath",
            "temperature T",
            "quasistatic level-raising protocol",
            "work and heat sign convention",
        ],
        "rows": rows,
        "limit_beta_gap_to_infinity": {
            "reset_error": 0.0,
            "beta_work": math.log(2.0),
            "beta_heat_to_bath": math.log(2.0),
        },
        "verdict": (
            "Shannon information becomes thermodynamic only after k_B, T, a bath, "
            "a Hamiltonian protocol, and work/heat instruments are supplied. "
            "The Landauer limit is derived inside that conditioned protocol; "
            "it is not generated by the dimensionless FIN kernel alone."
        ),
    }


def model_probability(a: np.ndarray, prep: int, t: float) -> np.ndarray:
    u = expm(-1j * t * a)
    p = np.abs(u[:, prep]) ** 2
    return np.maximum(p / p.sum(), 1e-15)


def empirical_kl(counts: np.ndarray, prediction: np.ndarray) -> float:
    empirical = counts / counts.sum()
    mask = empirical > 0
    return float(np.sum(empirical[mask] * np.log(empirical[mask] / prediction[mask])))


def compressed_legacy_candidate() -> np.ndarray:
    l24 = laplacian(radial_matrix(k_legacy, 24))
    a24, _ = shift_to_spd(l24)
    s12 = schur(a24, np.arange(0, 24, 2))
    # A scalar identity shift is invisible to unitary site probabilities.
    return s12 - np.trace(s12).real / 12.0 * np.eye(12)


def program70_blinded_challenge() -> dict:
    n = 12
    rng = np.random.default_rng(SEED)
    w = radial_matrix(k_strict, n)
    modulation = 1.0 + 0.08 * np.cos(2.0 * np.pi * np.arange(n) / n + 0.37)
    hidden_w = w * np.outer(modulation, modulation)
    hidden_a = laplacian(hidden_w)
    hidden_digest = hashlib.sha256(hidden_a.tobytes()).hexdigest()

    candidates = {
        "strict": laplacian(radial_matrix(k_strict, n)),
        "legacy_signed": laplacian(radial_matrix(k_legacy, n)),
        "legacy_absolute": laplacian(np.abs(radial_matrix(k_legacy, n))),
        "legacy_green_schur": compressed_legacy_candidate(),
        "nearest_neighbor_null": nearest_neighbor_laplacian(n, float(k_strict(1))),
    }
    preparations = [0, 3, 5]
    train_times = [0.35, 0.70, 1.05]
    test_times = [1.40, 1.75]
    shots = 20000
    records = {}
    for prep in preparations:
        for t in train_times + test_times:
            p = model_probability(hidden_a, prep, t)
            records[(prep, t)] = rng.multinomial(shots, p)

    def score(a: np.ndarray, scale: float, times: list[float]) -> float:
        vals = []
        for prep in preparations:
            for t in times:
                pred = model_probability(a, prep, scale * t)
                vals.append(empirical_kl(records[(prep, t)], pred))
        return float(np.mean(vals))

    ranking = []
    for name, a in candidates.items():
        fit = minimize_scalar(
            lambda x: score(a, float(x), train_times),
            bounds=(0.35, 1.65),
            method="bounded",
            options={"xatol": 1e-10},
        )
        scale = float(fit.x)
        ranking.append(
            {
                "model": name,
                "fitted_time_scale": scale,
                "training_mean_empirical_KL": score(a, scale, train_times),
                "heldout_mean_empirical_KL": score(a, scale, test_times),
            }
        )
    ranking.sort(key=lambda x: x["heldout_mean_empirical_KL"])

    fig, ax = plt.subplots(figsize=(8.5, 4.8), constrained_layout=True)
    names = [x["model"].replace("_", " ") for x in ranking]
    values = [x["heldout_mean_empirical_KL"] for x in ranking]
    ax.barh(names, values)
    ax.invert_yaxis()
    ax.set_xlabel("held-out mean empirical KL (lower is better)")
    ax.set_title("Blinded challenge: exact hidden generator excluded from candidates")
    fig.savefig(FIG / "program70_blinded_challenge.png", dpi=190)
    plt.close(fig)
    return {
        "preregistered_seed": SEED,
        "hidden_generator_sha256": hidden_digest,
        "hidden_generator_in_candidate_set": False,
        "preparations": preparations,
        "training_times": train_times,
        "heldout_times": test_times,
        "shots_per_record": shots,
        "candidate_ranking": ranking,
        "winner": ranking[0]["model"],
        "post_score_reveal": (
            "The hidden generator was a non-circulant 8% site-modulated strict "
            "weight matrix. No exact hidden model was offered to the scorer."
        ),
        "scope": (
            "This is a methodology and identifiability test on synthetic data, "
            "not evidence that nature selected the winning FIN kernel."
        ),
    }


def main() -> None:
    programs = {
        "program_61_continuum_functor": program61_continuum_functor(),
        "program_62_regularizer_independence": program62_regularizer_independence(),
        "program_63_compression_semigroup": program63_compression_semigroup(),
        "program_64_krein_signed_alternative": program64_krein_signed_alternative(),
        "program_65_operational_tomography": program65_operational_tomography(),
        "program_66_causal_order": program66_causal_order(),
        "program_67_calibration_identifiability": program67_calibration_identifiability(),
        "program_68_chiral_source_law": program68_chiral_source_law(),
        "program_69_landauer_protocol": program69_landauer_protocol(),
        "program_70_blinded_challenge": program70_blinded_challenge(),
    }
    payload = {
        "release": "10.7",
        "programs": "61-70",
        "seed": SEED,
        "scope": "finite dimensionless FIN core plus explicitly conditioned physics interfaces",
        "programs_data": programs,
        "global_verdict": (
            "Schur elimination is an exact compression semigroup and shifted "
            "resistance has a regularizer-independent limit, but the audited "
            "native kernel families do not automatically form a continuum "
            "projective functor. Signed/Krein structure does not create Markov "
            "physics; causal order, calibration, bath thermodynamics, and chiral "
            "sign provenance remain additional typed interfaces. Operational "
            "tomography and blinded held-out scoring are feasible once those "
            "interfaces are declared."
        ),
    }
    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps({"output": str(OUT), "figures": sorted(p.name for p in FIG.iterdir())}, indent=2))


if __name__ == "__main__":
    main()
