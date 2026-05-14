from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path


def coeffs_from_params(omega: float, phi: float, beta: float, eta: float) -> tuple[float, float, float]:
    """Strict-only surrogate map: kernel params -> effective coefficient bundle."""
    lam = abs(math.cos(phi)) / (1.0 + 0.3 * beta)
    kap = (omega * omega + 0.5 * eta) / (0.2 + abs(math.cos(phi)))
    eps = (beta * eta + 0.1 * omega) / (1.0 + abs(math.sin(phi)))
    return (lam, kap, eps)


def jacobian(params: tuple[float, float, float, float], h: float = 1e-6) -> list[list[float]]:
    base = coeffs_from_params(*params)
    cols: list[list[float]] = []
    for j, p in enumerate(params):
        perturbed = list(params)
        perturbed[j] = p + h
        f1 = coeffs_from_params(*perturbed)
        cols.append([(f1[i] - base[i]) / h for i in range(3)])
    return [[cols[c][r] for c in range(4)] for r in range(3)]


def det3(m: list[list[float]]) -> float:
    return (
        m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
    )


def main() -> None:
    root = Path(__file__).resolve().parent
    generated = root / "generated"
    generated.mkdir(parents=True, exist_ok=True)

    params = (0.18575, 0.16250, 1.0, 1.8)
    lam, kap, eps = coeffs_from_params(*params)
    j = jacobian(params)

    minors = {
        "freeze_omega_use_phi_beta_eta": [j[r][1] for r in range(3)],
        "freeze_phi_use_omega_beta_eta": [j[r][0] for r in range(3)],
        "freeze_beta_use_omega_phi_eta": [j[r][0] for r in range(3)],
        "freeze_eta_use_omega_phi_beta": [j[r][0] for r in range(3)],
    }
    dets = {
        "omega_frozen": det3([[j[r][1], j[r][2], j[r][3]] for r in range(3)]),
        "phi_frozen": det3([[j[r][0], j[r][2], j[r][3]] for r in range(3)]),
        "beta_frozen": det3([[j[r][0], j[r][1], j[r][3]] for r in range(3)]),
        "eta_frozen": det3([[j[r][0], j[r][1], j[r][2]] for r in range(3)]),
    }

    stable_slice_threshold = 1e-4
    identifiable_slices = [k for k, v in dets.items() if abs(v) > stable_slice_threshold]

    summary = {
        "checkpoint": "P1567_S517",
        "date_utc": "2026-05-14",
        "status": "PARTIAL_PASS_LOCAL_SLICE_ONLY",
        "strict_only": True,
        "legacy_bridge_used": False,
        "forward_coefficients": {
            "lambda_sm_eff": lam,
            "kappa_gr_eff": kap,
            "epsilon_mix_eff": eps,
        },
        "jacobian_3x4": j,
        "slice_3x3_determinants": dets,
        "stable_slice_threshold": stable_slice_threshold,
        "locally_identifiable_slices": identifiable_slices,
        "full_4param_identifiable_from_3observables": False,
        "qw2191_closed": False,
        "toe_closed": False,
        "missing_exports_witnesses_theorems": [
            "W1567A_underdetermined_orbit_witness",
            "W1567B_minimal_4th_observable_theorem",
            "E1567C_analytic_kernel_to_eom_gradient_hessian_export",
            "T1567D_inverse_stability_theorem",
        ],
        "next_honest_step": "P1568_construct_underdetermined_orbit_witness_and_4th_observable_candidate",
        "lay_summary": "3 obserwowalne nie domykają 4 parametrów; lokalne przekroje są odwracalne, pełna rekonstrukcja jeszcze nie.",
    }

    summary["audit_digest"] = hashlib.sha256(json.dumps(summary, sort_keys=True).encode("utf-8")).hexdigest()

    out = generated / "p1567_s517_strict_kernel_parameter_inversion_and_identifiability_summary.json"
    out.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"[P1567] wrote {out} status={summary['status']}")


if __name__ == "__main__":
    main()
