#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import UTC, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN1651 = GEN / "p1651_s601_strict_detj_lower_bound_certificate_summary.json"
IN1643 = GEN / "p1643_s593_strict_full_lagrangian_and_bidirectional_eom_obligation_summary.json"


def main() -> None:
    s51 = json.loads(IN1651.read_text(encoding="utf-8"))
    s43 = json.loads(IN1643.read_text(encoding="utf-8"))

    forward_chain = [
        {"step": "K_strict -> coefficients", "status": "PARTIAL", "evidence": "P1647 + P1648 + P1649 + P1651"},
        {"step": "coefficients -> full L_total", "status": "PARTIAL", "evidence": "operator-basis map + full Lagrangian sectors exported"},
        {"step": "full L_total -> EOM", "status": "PARTIAL", "evidence": "EOM bundle exported; theorem-level global variational log still open"},
    ]

    reverse_chain = [
        {"step": "EOM -> full L_total", "status": "OPEN", "blocker": "Helmholtz H1..H4 global discharge missing"},
        {"step": "full L_total -> coefficients", "status": "OPEN", "blocker": "global injective recovery theorem missing"},
        {"step": "coefficients -> K_strict", "status": "OPEN", "blocker": "global injectivity + selector nullspace theorem (QW-2191) missing"},
    ]

    summary = {
        "checkpoint": "P1652_S602",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "status": "KEEP_OPEN_P1652_FULL_CHAIN_BIDIRECTIONAL_DOSSIER_EXPORTED",
        "strict_only": True,
        "legacy_bridge_used": False,
        "route_target": s43["route_target"],
        "kernel_definition": s51["kernel_definition"],
        "full_lagrangian_density": s43["full_lagrangian_density"],
        "eom_system": s43["eom_system"],
        "forward_chain_matrix": forward_chain,
        "reverse_chain_matrix": reverse_chain,
        "o1_progress": s51["o1_certificate"],
        "strict_core_closure": {
            "status": "OPEN",
            "reason": "Reverse chain theorem-level obligations remain unresolved",
        },
        "next_honest_step": "Zacząć O2: wyeksportować theorem-level overlap cocycle consistency dla patchy parametrów i połączyć z O1/O3.",
        "lay_summary": "Mamy już pełny zapis toru do przodu (kernel→współczynniki→lagranżian→równania), a także listę braków dla drogi wstecz. To porządkuje dokładnie, co trzeba jeszcze udowodnić do domknięcia teorii.",
    }

    out = GEN / "p1652_s602_strict_full_chain_bidirectional_consistency_dossier_summary.json"
    out.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
