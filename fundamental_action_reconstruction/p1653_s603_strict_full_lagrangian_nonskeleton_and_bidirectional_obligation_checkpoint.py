#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import UTC, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN1652 = GEN / "p1652_s602_strict_full_chain_bidirectional_consistency_dossier_summary.json"


def main() -> None:
    s52 = json.loads(IN1652.read_text(encoding="utf-8"))

    full_lagrangian_nonskeleton = {
        "L_strict_scalar": {
            "kinetic": "(1/2) g^{μν}(∇_μ φ)(∇_ν φ)",
            "potential": "V_strict(φ)=m_φ^2 φ^2/2 + λ_3 φ^3/3! + λ_4 φ^4/4!",
        },
        "L_SM_gauge": {
            "SU3": "-(1/4) G^A_{μν}G_A^{μν}",
            "SU2": "-(1/4) W^I_{μν}W_I^{μν}",
            "U1": "-(1/4) B_{μν}B^{μν}",
        },
        "L_SM_fermions": "Σ_f i ψ̄_f γ^μ D_μ ψ_f - Σ_f y_f(ψ̄_{fL} H ψ_{fR} + h.c.)",
        "L_SM_higgs": "(D_μ H)^†(D^μ H) - μ_H^2(H^†H) - λ_H(H^†H)^2",
        "L_GR": "(M_Pl^2/2)R - Λ + c1 R^2 + c2 R_{μν}R^{μν} + c3 R_{μναβ}R^{μναβ}",
        "L_mix": "ξ_HR(H^†H)R + ξ_φR φ^2R + λ_{φH}φ^2(H^†H)",
        "L_total": "L_strict_scalar + L_SM_gauge + L_SM_fermions + L_SM_higgs + L_GR + L_mix",
    }

    reverse_obligation_backlog = [
        "Helmholtz integrability H1..H4 for coupled gauge-metric-scalar-fermion system",
        "Global injective recovery theorem from EL equations to coefficient chart",
        "Selector-nullspace removal theorem closing coefficients -> K_strict (QW-2191 dependent)",
    ]

    summary = {
        "checkpoint": "P1653_S603",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "status": "KEEP_OPEN_P1653_FULL_LAGRANGIAN_NONSKELETON_EXPORTED_WITH_REVERSE_BACKLOG",
        "strict_only": True,
        "legacy_bridge_used": False,
        "route_target": s52["route_target"],
        "kernel_definition": s52["kernel_definition"],
        "forward_chain": "K_strict -> coefficients -> full L_total -> EOM",
        "reverse_chain": "EOM -> full L_total -> coefficients -> K_strict",
        "full_lagrangian_density_nonskeleton": full_lagrangian_nonskeleton,
        "reverse_theorem_obligation_backlog": reverse_obligation_backlog,
        "final_strict_core_closure": {
            "status": "OPEN",
            "reason": "Reverse-chain theorem exports still missing; no false-pass admitted",
        },
        "next_honest_step": "Wyeksportować theorem-level Helmholtz witness (H1..H4) dla pełnego układu sprzężonego jako pierwszy twardy krok toru odwrotnego.",
        "lay_summary": "Mamy już pełny zapis lagranżianu strict i wiemy, jakie dokładnie twierdzenia są jeszcze potrzebne, żeby odtworzyć teorię także wstecz z równań ruchu.",
    }

    out = GEN / "p1653_s603_strict_full_lagrangian_nonskeleton_and_bidirectional_obligation_summary.json"
    out.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
