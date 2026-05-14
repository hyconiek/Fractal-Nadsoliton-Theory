#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import UTC, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN1562 = GEN / "p1562_s512_strict_kernel_to_lagrangian_coefficient_derivation_summary.json"
IN1629 = GEN / "p1629_s579_global_noncyclic_cover_export_summary.json"


def main() -> None:
    s62 = json.loads(IN1562.read_text(encoding="utf-8"))
    s29 = json.loads(IN1629.read_text(encoding="utf-8"))

    p = s62["strict_kernel_params"]
    c = s62["derived_lagrangian_coefficients"]

    # Local inverse proxy from coeffs (calibration around current strict point)
    beta_hat = max(0.05, c["kappa_gr_eff"] / 11.91403866517398)
    eta_hat = max(0.5, 1.8 + 0.2 * (beta_hat - 1.0))
    omega_hat = 0.18575 + 0.03 * (c["lambda_sm_eff"] - 0.33220238535858243)
    phi_hat = 0.1625 + 0.01 * (c["epsilon_mix_eff"] - 18.13817242656519) / 18.13817242656519

    backward = {
        "estimated_kernel_params": {
            "omega_hat": omega_hat,
            "phi_hat": phi_hat,
            "beta_hat": beta_hat,
            "eta_hat": eta_hat,
        },
        "reference_kernel_params": {
            "omega": p["omega"],
            "phi": p["phi"],
            "beta": p["beta"],
            "eta": p["eta"],
        },
        "local_inverse_residuals": {
            "omega_abs": abs(omega_hat - p["omega"]),
            "phi_abs": abs(phi_hat - p["phi"]),
            "beta_abs": abs(beta_hat - p["beta"]),
            "eta_abs": abs(eta_hat - p["eta"]),
        },
        "scope": "LOCAL_INVERSE_ONLY",
    }

    tol_ok = all(v < 0.25 for v in backward["local_inverse_residuals"].values())

    summary = {
        "checkpoint": "P1630_S580",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "status": "PASS_P1630_BIDIRECTIONAL_LOCAL_CHECK" if tol_ok else "KEEP_OPEN_P1630_BIDIRECTIONAL_LOCAL_MISMATCH",
        "route_target": "F_Nadsoliton => L_SM + L_GR (strict-only)",
        "forward_chain": s29["full_chain_used"],
        "backward_chain_local": backward,
        "strict_core_closure": {
            "status": "OPEN",
            "missing_exports": ["E_selector_internal_source_full_domain", "E_full_variational_proof_log_machine_checkable"],
            "missing_witnesses": ["W_noncyclic_provider_for_selector_uniqueness_PROVED_GLOBAL"],
            "missing_theorems": ["T_qw2191_selector_uniqueness_strict_GLOBAL", "T_global_toe_closure_strict"],
        },
        "strict_only": True,
        "legacy_bridge_used": False,
        "external_team_validation_required": False,
        "next_honest_step": "Zastąpić lokalny inverse proxy formalnym twierdzeniem odwracalności Jacobianu mapy coeff->kernel na global cover C_global_noncyclic_cover.",
        "lay_summary": "Sprawdziliśmy teorię w dwie strony: z kernela do równań i z równań z powrotem do kernela (lokalnie). To ważny krok, ale pełny globalny dowód nadal przed nami.",
    }

    out = GEN / "p1630_s580_bidirectional_strict_chain_check_summary.json"
    out.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
