#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import UTC, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN1641 = GEN / "p1641_s591_theorem_level_closure_requirement_matrix_summary.json"
IN1563 = GEN / "p1563_s513_strict_kernel_to_euler_lagrange_eom_export_summary.json"


def main() -> None:
    s41 = json.loads(IN1641.read_text(encoding="utf-8"))
    s63 = json.loads(IN1563.read_text(encoding="utf-8"))

    full_lagrangian = {
        "L_strict_scalar": "1/2 (∂_μ φ)(∂^μ φ) - V_strict(φ)",
        "L_SM_gauge": "-1/4 G^A_{μν}G_A^{μν} - 1/4 W^I_{μν}W_I^{μν} - 1/4 B_{μν}B^{μν}",
        "L_SM_fermions": "Σ_f i ψ̄_f γ^μ D_μ ψ_f - Σ_f y_f(ψ̄_{fL} H ψ_{fR} + h.c.)",
        "L_SM_higgs": "(D_μ H)^†(D^μ H) - μ_H^2 H^†H - λ_H(H^†H)^2",
        "L_GR": "M_Pl^2/2·R - Λ + c1 R^2 + c2 R_{μν}R^{μν} + c3 R_{μναβ}R^{μναβ}",
        "L_mix": "ξ_HR H^†H R + ξ_φR φ^2 R + λ_{φH} φ^2 H^†H",
        "L_total": "L_strict_scalar + L_SM_gauge + L_SM_fermions + L_SM_higgs + L_GR + L_mix",
    }

    eom_system = {
        "E_phi": "δS/δφ = 0",
        "E_H": "δS/δH† = 0",
        "E_psi_f": "δS/δψ̄_f = 0",
        "E_SU3": "D_μ G^{A μν} = J^{A ν}",
        "E_SU2": "D_μ W^{I μν} = J^{I ν}",
        "E_U1": "∂_μ B^{μν} = J_Y^ν",
        "E_metric": "δS/δg^{μν} = 0",
    }

    bidirectional_obligations = [
        {
            "id": "B1_forward_kernel_to_coeff",
            "direction": "K_strict -> coefficients",
            "status": "OPEN",
            "gap": "global identifiability theorem eliminating null directions",
        },
        {
            "id": "B2_forward_coeff_to_lagrangian",
            "direction": "coefficients -> L_total",
            "status": "PARTIAL",
            "gap": "theorem-level uniqueness of coefficient placement in full operator basis",
        },
        {
            "id": "B3_forward_lagrangian_to_eom",
            "direction": "L_total -> EOM",
            "status": "PARTIAL",
            "gap": "machine-checkable global variational log beyond local CAS proxy",
        },
        {
            "id": "B4_reverse_eom_to_lagrangian",
            "direction": "EOM -> L_total",
            "status": "OPEN",
            "gap": "global Helmholtz-type integrability witness for strict operator class",
        },
        {
            "id": "B5_reverse_lagrangian_to_coeff",
            "direction": "L_total -> coefficients",
            "status": "OPEN",
            "gap": "injective recovery map with overlap compatibility theorem",
        },
        {
            "id": "B6_reverse_coeff_to_kernel",
            "direction": "coefficients -> K_strict",
            "status": "OPEN",
            "gap": "strict-core uniqueness witness removing selector-source ambiguity (QW-2191)",
        },
    ]

    summary = {
        "checkpoint": "P1643_S593",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "status": "KEEP_OPEN_P1643_FULL_LAGRANGIAN_AND_BIDIRECTIONAL_OBLIGATIONS_EXPORTED",
        "route_target": s41["route_target"],
        "strict_only": True,
        "legacy_bridge_used": False,
        "full_lagrangian_density": full_lagrangian,
        "eom_system": eom_system,
        "forward_chain_reference": s63.get("chain", "K_strict -> coefficients -> L_total -> EOM"),
        "bidirectional_obligations": bidirectional_obligations,
        "strict_core_closure": {
            "status": "OPEN",
            "reason": "Bidirectional theorem-level witnesses B1..B6 not fully discharged",
        },
        "next_honest_step": "Wyeksportować theorem-level witness B4 (EOM -> L_total) dla pełnej klasy operatorów strict, następnie spiąć go z B5 i B6.",
        "lay_summary": "Mamy już pełny zapis lagranżianu i równania ruchu, ale formalny dowód działania w obie strony nadal wymaga kilku brakujących twierdzeń. To uczciwie trzyma status OPEN.",
    }

    out = GEN / "p1643_s593_strict_full_lagrangian_and_bidirectional_eom_obligation_summary.json"
    out.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
