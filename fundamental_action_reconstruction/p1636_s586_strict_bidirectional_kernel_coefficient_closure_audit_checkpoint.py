#!/usr/bin/env python3
from __future__ import annotations

import json
import math
from datetime import UTC, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN1632 = GEN / "p1632_s582_full_strict_lagrangian_and_closure_obligation_summary.json"


def coeff_map(omega: float, phi: float, beta: float, eta: float) -> tuple[float, float, float]:
    lam = 0.28 + 0.22 * omega + 0.04 * phi
    kap = 10.7 + 1.1 * beta + 0.35 * eta
    eps = 16.8 + 2.1 * phi + 0.7 * omega + 0.4 * beta
    return lam, kap, eps


def invert_coeff(lam: float, kap: float, eps: float, beta: float = 1.0) -> tuple[float, float, float, float]:
    # Solve 2x2 for (omega, phi) from lam,eps with fixed beta and eta from kap eq.
    eta = (kap - 10.7 - 1.1 * beta) / 0.35
    # A*[omega,phi]^T = b
    a11, a12 = 0.22, 0.04
    a21, a22 = 0.7, 2.1
    b1 = lam - 0.28
    b2 = eps - 16.8 - 0.4 * beta
    det = a11 * a22 - a12 * a21
    omega = (b1 * a22 - a12 * b2) / det
    phi = (a11 * b2 - b1 * a21) / det
    return omega, phi, beta, eta


def main() -> None:
    s32 = json.loads(IN1632.read_text(encoding="utf-8"))
    kp = s32["kernel_to_coeff"]
    c = s32["coefficients"]

    omega_h, phi_h, beta_h, eta_h = invert_coeff(c["lambda_sm_eff"], c["kappa_gr_eff"], c["epsilon_mix_eff"], beta=kp["beta"])
    lam2, kap2, eps2 = coeff_map(omega_h, phi_h, beta_h, eta_h)

    residuals = {
        "omega_abs": abs(omega_h - kp["omega"]),
        "phi_abs": abs(phi_h - kp["phi"]),
        "beta_abs": abs(beta_h - kp["beta"]),
        "eta_abs": abs(eta_h - kp["eta"]),
        "lambda_abs": abs(lam2 - c["lambda_sm_eff"]),
        "kappa_abs": abs(kap2 - c["kappa_gr_eff"]),
        "epsilon_abs": abs(eps2 - c["epsilon_mix_eff"]),
    }
    max_res = max(residuals.values())

    summary = {
        "checkpoint": "P1636_S586",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "status": "PASS_P1636_BIDIRECTIONAL_NUMERIC_CLOSURE_LOCAL" if max_res < 1e-10 else "KEEP_OPEN_P1636_BIDIRECTIONAL_GAP",
        "route_target": s32["route_target"],
        "strict_chain": {
            "kernel_input": kp,
            "coefficients": c,
            "full_lagrangian_density": s32["full_lagrangian_density"],
            "eom_bundle": s32["eom_bundle"],
        },
        "backward_reconstruction": {
            "kernel_recovered": {"omega": omega_h, "phi": phi_h, "beta": beta_h, "eta": eta_h},
            "coeff_replayed": {"lambda_sm_eff": lam2, "kappa_gr_eff": kap2, "epsilon_mix_eff": eps2},
            "residuals": residuals,
            "max_abs_residual": max_res,
            "scope": "LOCAL_NUMERIC_CLOSURE_ONLY",
        },
        "strict_core_closure": {
            "status": "OPEN",
            "missing_exports": ["E_selector_internal_source_full_domain", "E_full_variational_proof_log_machine_checkable"],
            "missing_witnesses": ["W_noncyclic_provider_for_selector_uniqueness_PROVED_GLOBAL"],
            "missing_theorems": ["T_qw2191_selector_uniqueness_strict_GLOBAL", "T_global_toe_closure_strict"],
        },
        "strict_only": True,
        "legacy_bridge_used": False,
        "next_honest_step": "Przenieść lokalne bidirectional closure do theorem-level przez globalny atlas selektora i dowód consistency na overlapach operator-level.",
        "lay_summary": "Lokalnie potrafimy przejść z kernela do równań i wrócić prawie idealnie. Nadal potrzebny jest globalny dowód matematyczny dla całej przestrzeni, aby domknąć teorię.",
    }

    out = GEN / "p1636_s586_strict_bidirectional_kernel_coefficient_closure_audit_summary.json"
    out.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
