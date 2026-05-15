#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import UTC, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN1707 = GEN / "p1707_s657_strict_full_nonproxy_eom_bundle_export_contract_checkpoint.json"
IN1662 = GEN / "p1662_s612_strict_full_lagrangian_explicit_density_summary.json"
OUT = GEN / "p1708_s658_strict_nonproxy_covariant_eom_first_explicit_formula_export_checkpoint.json"


def main() -> None:
    p1707 = json.loads(IN1707.read_text(encoding="utf-8"))
    p1662 = json.loads(IN1662.read_text(encoding="utf-8"))

    fullL = p1662.get("full_lagrangian_density_explicit", {})

    explicit_eom_formula_pack = {
        "metric_eom_template": (
            "E_{μν} := (M_Pl^2/2)G_{μν} + Λ g_{μν} + H_{μν}[c1 R^2 + c2 R_{αβ}R^{αβ} + c3 R_{αβγδ}R^{αβγδ}] "
            "- T^{(gauge+higgs+fermion+scalar+mix)}_{μν} = 0"
        ),
        "gauge_eom_template": "D_μ(Z_a F^{a μν}) + χ_RG ∇_μ(R F^{a μν}) - J^{a ν}_{matter} = 0",
        "higgs_eom_template": "D_μD^μH + μ_H^2 H + 2 λ_H(H†H)H + ξ_HR R H + δ_Yukawa/δH† + δ_mix/δH† = 0",
        "fermion_eom_template": "i γ^a e_a^{ μ}D_μ ψ_f - y_f H ψ_{fR/L} - m_f ψ_f = 0, and adjoint equation for ψ̄_f",
    }

    payload = {
        "checkpoint": "P1708_S658",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "strict_only": True,
        "legacy_bridge_used": False,
        "route_target": "F_Nadsoliton => L_SM + L_GR (strict-only)",
        "chain": "K_strict -> coefficients -> full explicit L_total -> first explicit nonproxy covariant EOM formulas",
        "contract_anchor": p1707.get("full_nonproxy_eom_bundle_contract", {}),
        "full_lagrangian_explicit_anchor": fullL,
        "first_explicit_nonproxy_formula_pack": explicit_eom_formula_pack,
        "formula_scope_note": "First explicit covariant formula export packet (template-level, convention-fixed). Full symbolic tensor expansion remains OPEN.",
        "bidirectional_status": {
            "forward_formula_export": "FIRST_EXPLICIT_PACK_EXPORTED",
            "reverse_local_identity": p1707.get("bidirectional_contract", {}).get("reverse_local_identity_status", "UNKNOWN"),
            "reverse_global_nonproxy": "OPEN_THEOREM_REQUIRED",
        },
        "status": "KEEP_OPEN_QG_THEOREM_LEVEL_REQUIRED",
        "open_obligations": [
            "full_symbolic_tensor_expansion_metric_sector",
            "full_symbolic_curved_gamma_spin_connection_expansion",
            "global_helmholtz_integrability_nonproxy",
            "brst_nilpotency_and_cohomology_nonproxy",
            "cutkosky_unitarity_full_sector",
            "counterterm_flow_renormalization_closure",
            "background_independence_family_theorem",
            "qw2191_selector_source_or_nonclosure_theorem",
        ],
        "next_honest_step": "Rozwinąć first_explicit_nonproxy_formula_pack do jawnej symbolicznej postaci tensorowej/spinorowej z jedną zamrożoną konwencją indeksową i podpisać local EL residual tests dla wszystkich sektorów.",
        "lay_summary": "Po raz pierwszy zapisaliśmy docelowe równania w pełnej kowariantnej formie sektorowej (jeszcze jako pakiet jawnych szablonów). To zbliża teorię do finalnej postaci obliczeniowej potrzebnej do ścisłych dowodów.",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
