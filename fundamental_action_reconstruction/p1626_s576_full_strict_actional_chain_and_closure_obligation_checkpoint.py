#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import UTC, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"

IN1562 = GEN / "p1562_s512_strict_kernel_to_lagrangian_coefficient_derivation_summary.json"
IN1622 = GEN / "p1622_s572_full_strict_lagrangian_density_and_eom_summary.json"
IN1625 = GEN / "p1625_s575_previous_candidate_failure_audit_and_next_strict_move_summary.json"


def load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    s62 = load(IN1562)
    s22 = load(IN1622)
    s25 = load(IN1625)

    params = s62["strict_kernel_params"]
    coeff = s62["derived_lagrangian_coefficients"]

    full_action = {
        "kernel": "K_strict(d)=cos(omega*d+phi)/(1+beta*d^eta)",
        "kernel_params": params,
        "coeff_mapping": {
            "lambda_sm_eff": coeff["lambda_sm_eff"],
            "kappa_gr_eff": coeff["kappa_gr_eff"],
            "epsilon_mix_eff": coeff["epsilon_mix_eff"],
        },
        "lagrangian_density": s22["full_lagrangian_density"],
        "action": "S = ∫ d^4x sqrt(-g) [L_strict_scalar + L_SM + L_GR + L_mix]",
        "eom": s22["euler_lagrange_eom"],
    }

    closure_obligations = {
        "missing_exports": s25["strict_core_closure"]["missing_exports"],
        "missing_witnesses": s25["strict_core_closure"]["missing_witnesses"],
        "missing_theorems": s25["strict_core_closure"]["missing_theorems"],
    }

    summary = {
        "checkpoint": "P1626_S576",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "status": "KEEP_OPEN_P1626_ACTION_CHAIN_COMPLETE_BUT_CLOSURE_OPEN",
        "route_target": "F_Nadsoliton => L_SM + L_GR (strict-only)",
        "full_strict_chain": full_action,
        "strict_core_closure": {
            "status": "OPEN",
            **closure_obligations,
        },
        "strict_only": True,
        "legacy_bridge_used": False,
        "external_team_validation_required": False,
        "next_honest_step": "Dowód theorem-level: T_qw2191_selector_uniqueness_strict z witness PROVED + export E_selector_internal_source_full_domain.",
        "lay_summary": "Mamy pełny przepis od kernela do równań ruchu; brakują już tylko formalne dowody domykające teorię, nie kolejne warianty modelu.",
    }

    out = GEN / "p1626_s576_full_strict_actional_chain_and_closure_obligation_summary.json"
    out.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
