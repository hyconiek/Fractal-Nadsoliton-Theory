#!/usr/bin/env python3
"""Reproducible checks and figures for the FIN dual-dynamics note.

The program constructs the 12-site cyclic-distance matrix sampled from the
strict FIN profile, forms its weighted graph Laplacian A = sI - W, audits the
unitary and heat functional calculi, and simulates a finite two-path
(double-slit) operational experiment.

Outputs
-------
FIN_Dual_Dynamics_reference_results.json
FIN_One_Spectrum_Two_Dynamics.png
FIN_Double_Slit_Operational_Simulation.png
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


N = 12
OMEGA = 0.18575
PHI = 0.16250
BETA = 1.0
ETA = 1.8
OUT = Path(__file__).resolve().parent


def strict_profile(distance: int | np.ndarray) -> np.ndarray:
    d = np.asarray(distance, dtype=float)
    return np.cos(OMEGA * d + PHI) / (1.0 + BETA * d**ETA)


def cyclic_distance(x: int, y: int) -> int:
    raw = abs(x - y)
    return min(raw, N - raw)


def construct_operator() -> tuple[np.ndarray, float, np.ndarray]:
    weights = strict_profile(np.arange(1, N // 2 + 1))
    W = np.zeros((N, N), dtype=float)
    for x in range(N):
        for y in range(N):
            if x != y:
                W[x, y] = weights[cyclic_distance(x, y) - 1]
    s = float(W[0].sum())
    A = s * np.eye(N) - W
    return W, s, A


W, S, A = construct_operator()
EVALS, EVECS = np.linalg.eigh(A)
EVALS[np.abs(EVALS) < 1e-13] = 0.0


def unitary(time: float) -> np.ndarray:
    return (EVECS * np.exp(-1j * time * EVALS)) @ EVECS.T


def heat(time: float) -> np.ndarray:
    return (EVECS * np.exp(-time * EVALS)) @ EVECS.T


def shannon(probability: np.ndarray) -> float:
    p = np.clip(np.real(probability), 0.0, None)
    p = p / p.sum()
    nz = p > 0.0
    return float(-np.sum(p[nz] * np.log(p[nz])))


def localized_profiles(time: float) -> tuple[np.ndarray, np.ndarray]:
    U = unitary(time)
    P = heat(time)
    return np.abs(U[:, 0]) ** 2, np.real(P[:, 0])


def double_slit_profiles(
    time: float, phase: float = 0.0, coherence: float = 1.0
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return coherent/dephased quantum, incoherent quantum, and heat profiles.

    Slits are the two vertices -2 and +2.  ``coherence`` is the real overlap
    of the corresponding environment records: 1 erases which-path data and
    0 makes the two records orthogonal.
    """
    slit_a, slit_b = (-2) % N, 2
    U = unitary(time)
    ua, ub = U[:, slit_a], U[:, slit_b]
    incoherent = 0.5 * (np.abs(ua) ** 2 + np.abs(ub) ** 2)
    cross = np.real(np.exp(-1j * phase) * ua * np.conj(ub))
    coherent = incoherent + coherence * cross
    entrance = np.zeros(N)
    entrance[[slit_a, slit_b]] = 0.5
    diffusive = np.real(heat(time) @ entrance)
    return coherent, incoherent, diffusive, cross


def first_interference_peak(tmax: float = 5.0) -> tuple[float, float]:
    times = np.linspace(0.0, tmax, 10001)
    strengths = np.empty_like(times)
    for j, time in enumerate(times):
        _, _, _, cross = double_slit_profiles(float(time))
        strengths[j] = 0.5 * np.sum(np.abs(cross))
    candidates = np.where(
        (strengths[1:-1] > strengths[:-2])
        & (strengths[1:-1] >= strengths[2:])
    )[0] + 1
    candidates = candidates[times[candidates] > 1e-3]
    if candidates.size == 0:
        index = int(np.argmax(strengths))
    else:
        index = int(candidates[0])
    return float(times[index]), float(strengths[index])


def bisection_mixing_time(epsilon: float) -> float:
    uniform = np.ones(N) / N

    def tv(time: float) -> float:
        return float(0.5 * np.sum(np.abs(heat(time)[:, 0] - uniform)))

    lo, hi = 0.0, 1.0
    while tv(hi) > epsilon:
        hi *= 2.0
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if tv(mid) > epsilon:
            lo = mid
        else:
            hi = mid
    return hi


def quantum_cesaro_limit() -> np.ndarray:
    # Sum squared projections onto each distinct eigenspace.
    source = np.zeros(N)
    source[0] = 1.0
    limit = np.zeros(N)
    used = np.zeros(N, dtype=bool)
    for j, value in enumerate(EVALS):
        if used[j]:
            continue
        indices = np.where(np.isclose(EVALS, value, atol=1e-11, rtol=0.0))[0]
        used[indices] = True
        projector = EVECS[:, indices] @ EVECS[:, indices].T
        amplitudes = projector @ source
        limit += np.abs(amplitudes) ** 2
    return limit


def make_dynamics_figure() -> None:
    times = np.linspace(0.0, 20.0, 1201)
    return_u, return_p, entropy_u, entropy_p = [], [], [], []
    for time in times:
        pu, pp = localized_profiles(float(time))
        return_u.append(pu[0])
        return_p.append(pp[0])
        entropy_u.append(shannon(pu))
        entropy_p.append(shannon(pp))

    plt.style.use("seaborn-v0_8-whitegrid")
    fig, axes = plt.subplots(2, 1, figsize=(10.2, 7.3), sharex=True)
    blue, orange = "#1f5a99", "#d55e00"
    axes[0].plot(times, return_u, color=blue, lw=2.0, label=r"unitary $|U_t(0,0)|^2$")
    axes[0].plot(times, return_p, color=orange, lw=2.0, label=r"diffusive $P_t(0,0)$")
    axes[0].axhline(1 / N, color="0.35", ls="--", lw=1.1, label=r"$1/12$")
    axes[0].set_ylabel("return probability")
    axes[0].legend(ncol=3, fontsize=9, loc="upper right")
    axes[0].set_title("One spectrum, two temporal calculi")

    axes[1].plot(times, entropy_u, color=blue, lw=2.0, label="unitary populations")
    axes[1].plot(times, entropy_p, color=orange, lw=2.0, label="diffusive law")
    axes[1].axhline(np.log(N), color="0.35", ls="--", lw=1.1, label=r"$\log 12$")
    axes[1].set_xlabel("dimensionless time $t$")
    axes[1].set_ylabel("Shannon entropy (nats)")
    axes[1].legend(ncol=3, fontsize=9, loc="lower right")
    fig.tight_layout()
    fig.savefig(OUT / "FIN_One_Spectrum_Two_Dynamics.png", dpi=240, bbox_inches="tight")
    plt.close(fig)


def make_double_slit_figure(peak_time: float) -> dict[str, float | int]:
    coherent, incoherent, diffusive, cross = double_slit_profiles(peak_time)
    sites = np.arange(N)

    # Interference through time.
    times = np.linspace(0.0, 5.0, 401)
    interference = np.empty((N, times.size))
    strength = np.empty(times.size)
    for j, time in enumerate(times):
        _, _, _, term = double_slit_profiles(float(time))
        interference[:, j] = term
        strength[j] = 0.5 * np.sum(np.abs(term))

    # Phase scan at the symmetry detector x=0.
    phases = np.linspace(0.0, 2.0 * np.pi, 361)
    p_eta = {}
    for eta in (1.0, 0.5, 0.0):
        p_eta[eta] = np.array(
            [double_slit_profiles(peak_time, phase=float(phi), coherence=eta)[0][0] for phi in phases]
        )

    plt.style.use("seaborn-v0_8-whitegrid")
    fig = plt.figure(figsize=(12.0, 9.0))
    grid = fig.add_gridspec(2, 2, width_ratios=(1.08, 1.0), hspace=0.34, wspace=0.25)

    ax0 = fig.add_subplot(grid[0, 0])
    width = 0.25
    ax0.bar(sites - width, coherent, width, color="#1f5a99", label=r"coherent, $\eta=1$")
    ax0.bar(sites, incoherent, width, color="#8a8a8a", label=r"which-path, $\eta=0$")
    ax0.bar(sites + width, diffusive, width, color="#d55e00", label=r"diffusion $P_t$")
    ax0.set_xticks(sites)
    ax0.set_xlabel("detector vertex $x$")
    ax0.set_ylabel("detection probability")
    ax0.set_title(rf"Twelve-detector screen at $t_*={peak_time:.3f}$")
    ax0.legend(fontsize=8.5)

    ax1 = fig.add_subplot(grid[0, 1])
    ax1.bar(sites, cross, color=np.where(cross >= 0.0, "#0072b2", "#cc3311"))
    ax1.axhline(0.0, color="black", lw=0.8)
    ax1.set_xticks(sites)
    ax1.set_xlabel("detector vertex $x$")
    ax1.set_ylabel(r"$p_{\rm coh}-p_{\rm mix}$")
    ax1.set_title("Signed interference term")

    ax2 = fig.add_subplot(grid[1, 0])
    vmax = float(np.max(np.abs(interference)))
    image = ax2.imshow(
        interference,
        origin="lower",
        aspect="auto",
        extent=(times[0], times[-1], -0.5, N - 0.5),
        cmap="RdBu_r",
        vmin=-vmax,
        vmax=vmax,
        interpolation="nearest",
    )
    ax2.axvline(peak_time, color="black", ls="--", lw=1.0)
    ax2.set_yticks(sites)
    ax2.set_xlabel("dimensionless propagation time $t$")
    ax2.set_ylabel("detector vertex $x$")
    ax2.set_title(r"Interference map $p_{\rm coh}(x,t)-p_{\rm mix}(x,t)$")
    fig.colorbar(image, ax=ax2, fraction=0.046, pad=0.04)

    ax3 = fig.add_subplot(grid[1, 1])
    colors = {1.0: "#1f5a99", 0.5: "#7b4ab5", 0.0: "#8a8a8a"}
    for eta in (1.0, 0.5, 0.0):
        ax3.plot(phases, p_eta[eta], lw=2.0, color=colors[eta], label=rf"record overlap $\eta={eta:g}$")
    ax3.set_xticks([0.0, np.pi / 2, np.pi, 3 * np.pi / 2, 2 * np.pi])
    ax3.set_xticklabels(["0", r"$\pi/2$", r"$\pi$", r"$3\pi/2$", r"$2\pi$"])
    ax3.set_xlabel(r"relative slit phase $\phi$")
    ax3.set_ylabel(r"probability at symmetry detector $x=0$")
    ax3.set_title("Phase fringe and which-path record")
    ax3.legend(fontsize=8.5)

    fig.suptitle(
        "Finite FIN double-slit analogue: preparation, propagation, instrument, and record",
        fontsize=14,
        y=0.985,
    )
    fig.savefig(OUT / "FIN_Double_Slit_Operational_Simulation.png", dpi=240, bbox_inches="tight")
    plt.close(fig)

    detector_mix = float(incoherent[0])
    pmax = float(np.max(p_eta[1.0]))
    pmin = float(np.min(p_eta[1.0]))
    visibility = (pmax - pmin) / (pmax + pmin) if pmax + pmin else 0.0
    return {
        "slit_a": int((-2) % N),
        "slit_b": 2,
        "peak_time": peak_time,
        "total_variation_interference": float(0.5 * np.sum(np.abs(cross))),
        "symmetry_detector": 0,
        "symmetry_detector_mixture_probability": detector_mix,
        "symmetry_detector_phase_max": pmax,
        "symmetry_detector_phase_min": pmin,
        "symmetry_detector_visibility": float(visibility),
        "coherent_probability_sum": float(np.sum(coherent)),
        "incoherent_probability_sum": float(np.sum(incoherent)),
        "diffusive_probability_sum": float(np.sum(diffusive)),
    }


def main() -> None:
    peak_time, peak_strength = first_interference_peak()
    make_dynamics_figure()
    slit = make_double_slit_figure(peak_time)

    rng = np.random.default_rng(20260719)
    f = rng.normal(size=N) + 1j * rng.normal(size=N)
    lhs = float(np.real(np.vdot(f, A @ f)))
    rhs = float(
        0.5
        * sum(W[x, y] * abs(f[x] - f[y]) ** 2 for x in range(N) for y in range(N))
    )

    snapshots = {}
    for time in (0.1, 1.0, 5.0):
        pu, pp = localized_profiles(time)
        U, P = unitary(time), heat(time)
        snapshots[str(time)] = {
            "unitary_return": float(pu[0]),
            "unitary_population_entropy_nats": shannon(pu),
            "diffusive_return": float(pp[0]),
            "diffusive_entropy_nats": shannon(pp),
            "unitarity_residual_frobenius": float(np.linalg.norm(U.conj().T @ U - np.eye(N))),
            "stochastic_row_sum_residual_max": float(np.max(np.abs(P.sum(axis=1) - 1.0))),
            "stochastic_minimum_entry": float(P.min()),
            "chapman_kolmogorov_heat_frobenius": float(np.linalg.norm(heat(2 * time) - P @ P)),
            "chapman_kolmogorov_born_frobenius": float(
                np.linalg.norm(np.abs(unitary(2 * time)) ** 2 - (np.abs(U) ** 2) @ (np.abs(U) ** 2))
            ),
        }

    qbar = quantum_cesaro_limit()
    raw_w_evals, raw_w_evecs = np.linalg.eigh(W)
    raw_heat = (raw_w_evecs * np.exp(-raw_w_evals)) @ raw_w_evecs.T
    profile_extension = {str(d): float(strict_profile(d)) for d in range(1, 13)}

    results = {
        "model": {
            "vertices": N,
            "omega": OMEGA,
            "phi": PHI,
            "beta": BETA,
            "eta": ETA,
            "weights_d1_to_d6": [float(v) for v in strict_profile(np.arange(1, 7))],
            "row_sum_s": S,
            "s_display_15_decimals": f"{S:.15f}",
            "A_eigenvalues": [float(v) for v in EVALS],
            "spectral_gap": float(EVALS[1]),
            "largest_A_eigenvalue": float(EVALS[-1]),
            "raw_W_minimum_eigenvalue": float(raw_w_evals[0]),
        },
        "audits": {
            "A_symmetry_residual_frobenius": float(np.linalg.norm(A - A.T)),
            "A_constant_mode_residual": float(np.linalg.norm(A @ np.ones(N))),
            "A_minimum_eigenvalue": float(EVALS[0]),
            "dirichlet_lhs": lhs,
            "dirichlet_rhs": rhs,
            "dirichlet_identity_absolute_residual": abs(lhs - rhs),
            "short_time_diffusive_escape_rate": S,
            "short_time_unitary_variance": float(np.sum(W[:, 0] ** 2)),
        },
        "snapshots": snapshots,
        "mixing_times_total_variation": {
            str(eps): bisection_mixing_time(eps) for eps in (0.1, 0.05, 0.01, 0.001)
        },
        "unitary_cesaro_position_law": [float(v) for v in qbar],
        "unitary_cesaro_tv_from_uniform": float(0.5 * np.sum(np.abs(qbar - 1.0 / N))),
        "raw_W_heat_at_t1": {
            "row_sum": float(raw_heat[0].sum()),
            "minimum_entry": float(raw_heat.min()),
            "operator_norm": float(np.max(np.exp(-raw_w_evals))),
        },
        "profile_extension_d1_to_d12": profile_extension,
        "first_negative_profile_distance": next(d for d in range(1, 13) if strict_profile(d) < 0),
        "double_slit": slit,
        "double_slit_peak_strength_crosscheck": peak_strength,
    }
    (OUT / "FIN_Dual_Dynamics_reference_results.json").write_text(
        json.dumps(results, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(results, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
