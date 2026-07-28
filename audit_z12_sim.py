#!/usr/bin/env python3
"""Independent numerical audit for ``Z12 sim.html``.

The audit evaluates the mathematical model directly with NumPy.  It does not
execute the page and therefore cannot be fooled by its canvas presentation.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np


N = 12
K = np.array(
    [
        0.46998567264502017,
        0.19204355169010282,
        0.09142861427792497,
        0.04702916874565040,
        0.02413122336363006,
        0.011070817321442113,
    ],
    dtype=float,
)
S_REFERENCE = 1.660307278766099
HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "Z12_sim_methodological_audit.json"


def matrices() -> tuple[np.ndarray, np.ndarray]:
    w = np.zeros((N, N), dtype=float)
    for x in range(N):
        for y in range(N):
            delta = abs(x - y) % N
            d = min(delta, N - delta)
            if d:
                w[x, y] = K[d - 1]
    s = float(w[0].sum())
    return w, s * np.eye(N) - w


W, A = matrices()
EIGENVALUES, EIGENVECTORS = np.linalg.eigh(A)


def state(phi: float, eta: float, time: float, a: int, b: int) -> np.ndarray:
    """Return the evolved partially coherent two-source density matrix."""
    basis = np.eye(N)
    ket_a, ket_b = basis[:, a], basis[:, b]
    rho0 = 0.5 * (
        np.outer(ket_a, ket_a)
        + np.outer(ket_b, ket_b)
        + eta * np.exp(-1j * phi) * np.outer(ket_a, ket_b)
        + eta * np.exp(+1j * phi) * np.outer(ket_b, ket_a)
    )
    unitary = (EIGENVECTORS * np.exp(-1j * time * EIGENVALUES)) @ EIGENVECTORS.T
    return unitary @ rho0 @ unitary.conj().T


def observables(
    phi: float, eta: float, time: float, a: int = 10, b: int = 2
) -> dict[str, float | np.ndarray]:
    rho = state(phi, eta, time, a, b)
    # J[x,y] is the current directed from x to y.
    current = 2.0 * W * np.imag(rho.T)
    c1 = float(sum(current[x, (x + 1) % N] for x in range(N)))
    c_chi = float(
        sum(
            d * current[x, (x + d) % N]
            for x in range(N)
            for d in range(1, N // 2)
        )
    )
    rho_dot = -1j * (A @ rho - rho @ A)
    probability = np.real(np.diag(rho))
    continuity_residual = float(
        np.max(np.abs(np.real(np.diag(rho_dot)) + current.sum(axis=1)))
    )
    return {
        "rho": rho,
        "current": current,
        "probability": probability,
        "C1": c1,
        "C_chi": c_chi,
        "normalization_error": abs(float(probability.sum()) - 1.0),
        "continuity_residual": continuity_residual,
    }


def main() -> None:
    t = 1.5525
    phi = 0.7
    plus = observables(+phi, 1.0, t)
    minus = observables(-phi, 1.0, t)
    half = observables(+phi, 0.5, t)
    incoherent = observables(+phi, 0.0, t)
    symmetric = observables(0.0, 1.0, t)
    pi_state = observables(np.pi, 1.0, t)
    separation_control = observables(0.0, 1.0, t, 10, 3)

    tests = {
        "row_sum_matches_declared_s": abs(float(W[0].sum()) - S_REFERENCE) < 1e-12,
        "A_is_positive_semidefinite": float(EIGENVALUES.min()) > -1e-12,
        "A_has_constant_zero_mode": abs(float(EIGENVALUES.min())) < 1e-12,
        "probability_is_normalized": plus["normalization_error"] < 1e-12,
        "continuity_equation_holds": plus["continuity_residual"] < 1e-12,
        "current_is_antisymmetric": float(
            np.max(np.abs(plus["current"] + plus["current"].T))
        )
        < 1e-12,
        "C_chi_is_odd_in_phi": abs(float(plus["C_chi"] + minus["C_chi"])) < 1e-12,
        "C_chi_zero_at_phi_zero": abs(float(symmetric["C_chi"])) < 1e-12,
        "C_chi_zero_at_phi_pi": abs(float(pi_state["C_chi"])) < 1e-12,
        "C_chi_is_linear_in_eta": abs(
            float(half["C_chi"]) - 0.5 * float(plus["C_chi"])
        )
        < 1e-12,
        "eta_zero_kills_C_chi": abs(float(incoherent["C_chi"])) < 1e-12,
        "eta_zero_does_not_kill_all_local_currents": float(
            np.max(np.abs(incoherent["current"]))
        )
        > 1e-6,
        "changed_separation_is_not_a_selector": abs(
            float(separation_control["C_chi"])
        )
        < 1e-12,
        "nearest_neighbour_control_is_blind_for_10_2": max(
            abs(float(plus["C1"])), abs(float(minus["C1"]))
        )
        < 1e-12,
        "long_range_C_chi_detects_phase_for_10_2": abs(float(plus["C_chi"]))
        > 1e-6,
    }

    report = {
        "status": "PASS" if all(tests.values()) else "FAIL",
        "tests": tests,
        "numbers": {
            "declared_s": S_REFERENCE,
            "computed_s": float(W[0].sum()),
            "lambda_min": float(EIGENVALUES.min()),
            "lambda_max": float(EIGENVALUES.max()),
            "C_chi_phi_plus_0_7": float(plus["C_chi"]),
            "C_chi_phi_minus_0_7": float(minus["C_chi"]),
            "C1_phi_plus_0_7": float(plus["C1"]),
            "max_local_current_eta_zero": float(
                np.max(np.abs(incoherent["current"]))
            ),
            "continuity_residual": float(plus["continuity_residual"]),
            "normalization_error": float(plus["normalization_error"]),
        },
        "methodological_boundary": {
            "C1": "nearest-neighbour first-harmonic control; blind to the 10--2 preparation phase",
            "C_chi": "orientation- and lift-dependent long-range group-current observable",
            "claim": "a detector of an inserted chiral phase, not an internal FIN selector theorem",
        },
    }
    OUTPUT.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps(report, indent=2, ensure_ascii=False))
    if report["status"] != "PASS":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
