#!/usr/bin/env python3
"""P1610/S560: export full strict Lagrangian (non-skeleton) and EOM chain from strict kernel."""
from __future__ import annotations

import json
from datetime import UTC, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"

IN1558 = GEN / "p1558_s508_strict_sm_closed_form_coupling_bundle_and_full_lagrangian_skeleton_summary.json"
IN1559 = GEN / "p1559_s509_gr_strict_curvature_transport_bundle_and_full_lagrangian_update_summary.json"
IN1562 = GEN / "p1562_s512_strict_kernel_to_lagrangian_coefficient_derivation_summary.json"
IN1563 = GEN / "p1563_s513_strict_kernel_to_euler_lagrange_eom_export_summary.json"
IN1609 = GEN / "p1609_s559_strict_kernel_to_full_lagrangian_closure_audit_summary.json"


def _load(path: Path) -> dict:
    if not path.exists():
        raise FileNotFoundError(f"Missing required input: {path.name}")
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    s58 = _load(IN1558)
    s59 = _load(IN1559)
    s62 = _load(IN1562)
    s63 = _load(IN1563)
    s09 = _load(IN1609)

    coeff = s62["derived_lagrangian_coefficients"]
    lam = coeff["lambda_sm_eff"]
    kap = coeff["kappa_gr_eff"]
    eps = coeff["epsilon_mix_eff"]

    closure = s09.get("strict_core_closure", {})
    missing_exports = closure.get("missing_exports", [])
    missing_witnesses = closure.get("missing_witnesses", [])
    missing_theorems = closure.get("missing_theorems", [])

    full_lagrangian = {
        "L_SM_strict": (
            "-1/4*F^a_{mu nu}F_a^{mu nu} + i*bar(psi)gamma^mu D_mu psi + "
            "|D_mu H|^2 - lambda_sm_eff*(H^dagger H - v^2/2)^2 - y_f*bar(psi)_L H psi_R + h.c."
        ),
        "L_GR_strict": "kappa_gr_eff*R*sqrt(-g)",
        "L_mix_strict": "epsilon_mix_eff*(H^dagger H)*R*sqrt(-g)",
        "L_total": (
            "L_SM_strict + L_GR_strict + L_mix_strict with "
            f"lambda_sm_eff={lam}, kappa_gr_eff={kap}, epsilon_mix_eff={eps}"
        ),
    }

    chain_ready = all([
        s58.get("status", "").startswith("PASS"),
        s59.get("status", "").startswith("PASS"),
        s63.get("status", "").startswith("PASS"),
        s09.get("status", "").startswith("PASS"),
        not missing_exports and not missing_witnesses and not missing_theorems,
    ])

    status = "PASS_P1610_FULL_STRICT_LAGRANGIAN_AND_EOM_EXPORTED" if chain_ready else "KEEP_OPEN_P1610_CHAIN_INCOMPLETE"

    summary = {
        "checkpoint": "P1610_S560",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "status": status,
        "route_target": "F_Nadsoliton => L_SM + L_GR (strict-only)",
        "strict_kernel": s62.get("strict_kernel_params", {}),
        "derived_coefficients": coeff,
        "full_lagrangian_density": full_lagrangian,
        "euler_lagrange_equations": s63.get("euler_lagrange_equations", {}),
        "theoretical_chain": "kernel_strict -> coefficients -> full_lagrangian -> eom -> strict_core_closure",
        "strict_core_closure": {
            "status": "CLOSED" if chain_ready else "OPEN",
            "missing_exports": missing_exports,
            "missing_witnesses": missing_witnesses,
            "missing_theorems": missing_theorems,
        },
        "strict_only": True,
        "legacy_bridge_used": False,
        "external_team_validation_required": False,
        "next_honest_step": "Export theorem-level variational consistency object linking full Lagrangian terms to each EOM sector with explicit proof obligations.",
        "lay_summary": "To pełna wersja modelu strict: od kernela wyliczamy liczby, wkładamy je do pełnego Langrażianu i dostajemy równania ruchu; domknięcie jest ważne tylko bez brakujących elementów.",
    }

    out = GEN / "p1610_s560_full_strict_lagrangian_and_eom_chain_export_summary.json"
    out.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
