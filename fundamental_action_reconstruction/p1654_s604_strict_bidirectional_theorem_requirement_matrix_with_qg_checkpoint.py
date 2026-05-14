#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import UTC, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN1653 = GEN / "p1653_s603_strict_full_lagrangian_nonskeleton_and_bidirectional_obligation_summary.json"


def main() -> None:
    s53 = json.loads(IN1653.read_text(encoding="utf-8"))

    theorem_matrix = {
        "forward_chain": [
            {"gate": "K_strict -> coefficients global identifiability theorem", "status": "PARTIAL"},
            {"gate": "coefficients -> full L_total completeness theorem", "status": "PARTIAL"},
            {"gate": "full L_total -> full EOM variational theorem", "status": "PARTIAL"},
        ],
        "reverse_chain": [
            {"gate": "EOM -> L_total Helmholtz H1..H4 discharge", "status": "OPEN"},
            {"gate": "L_total -> coefficients global injective recovery", "status": "OPEN"},
            {"gate": "coefficients -> K_strict selector-nullspace closure (QW-2191)", "status": "OPEN"},
        ],
    }

    qg_closure_gates = {
        "renormalization": {
            "status": "OPEN",
            "required_export": "strict UV counterterm-closure theorem for gravity-coupled sectors",
        },
        "unitarity": {
            "status": "OPEN",
            "required_export": "ghost-free + optical-theorem compatible amplitude theorem in strict sector",
        },
        "background_independence": {
            "status": "OPEN",
            "required_export": "background-independent observable algebra and evolution consistency theorem",
        },
    }

    summary = {
        "checkpoint": "P1654_S604",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "status": "KEEP_OPEN_P1654_BIDIRECTIONAL_AND_QG_THEOREM_MATRIX_EXPORTED",
        "strict_only": True,
        "legacy_bridge_used": False,
        "route_target": "F_Nadsoliton => L_SM + L_GR (strict-only)",
        "kernel_definition": s53["kernel_definition"],
        "full_chain_forward": s53["forward_chain"],
        "full_chain_reverse": s53["reverse_chain"],
        "theorem_requirement_matrix": theorem_matrix,
        "quantum_gravity_closure_gates": qg_closure_gates,
        "final_strict_core_closure": {
            "status": "OPEN",
            "reason": "Reverse-chain and QG theorem-level gates remain OPEN",
        },
        "next_honest_step": "Zbudować pierwszy formalny witness H1 (lokalna warunkowalność Helmholtza) dla pełnego układu sprzężonego i podpiąć go do macierzy bramek P1654.",
        "lay_summary": "Mamy już pełną mapę, czego brakuje do domknięcia teorii: nie tylko matematyczne cofanie z równań ruchu, ale też trzy twarde warunki kwantowej grawitacji: renormalizacja, unitarność i niezależność od tła.",
    }

    out = GEN / "p1654_s604_strict_bidirectional_theorem_requirement_matrix_with_qg_summary.json"
    out.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
