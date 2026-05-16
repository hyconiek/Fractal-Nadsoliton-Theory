#!/usr/bin/env python3
"""P1835 S785 strict forward->reverse unlock matrix checkpoint."""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"


def load(name: str) -> dict:
    p = GEN / name
    if not p.exists():
        return {"_missing": name, "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(p.read_text(encoding="utf-8"))


def main() -> None:
    p1834 = load("p1834_s784_strict_full_lagrangian_sector_export_gate_checkpoint.json")
    p1833 = load("p1833_s783_strict_full_lagrangian_and_reverse_closure_worklist_checkpoint.json")

    sector_gate = p1834.get("full_lagrangian_sector_gate", [])
    all_forward_ready = bool(sector_gate) and all(item.get("sector_gate_pass", False) for item in sector_gate)

    reverse = p1833.get("reverse_qg_worklist", [])
    unlock_matrix = []
    for item in reverse:
        gate = item.get("theorem_gate", "UNKNOWN")
        unlock_matrix.append(
            {
                "theorem_gate": gate,
                "forward_prerequisite": "ALL_FULL_LAGRANGIAN_SECTORS_PASS",
                "forward_prerequisite_met": all_forward_ready,
                "current_status": item.get("status", "OPEN_OBSTRUCTION_WITH_TRACE"),
                "unlock_state": "UNLOCKED_FOR_EXECUTION" if all_forward_ready else "LOCKED_BY_FORWARD_CHAIN",
            }
        )

    out = {
        "packet_id": "P1835",
        "stage_id": "S785",
        "status": "OPEN_OBSTRUCTION_WITH_TRACE",
        "route": "strict_only",
        "forward_chain_ready": all_forward_ready,
        "reverse_unlock_matrix": unlock_matrix,
        "technical_progress": "Forward-completeness gate is now wired into an explicit reverse-theorem unlock matrix.",
        "proven": "Reverse theorem gates remain locked until full non-skeleton Lagrangian sector exports are complete.",
        "open": "Current forward sector gate is not passed, so reverse QG closure execution stays blocked.",
        "false_pass_risk": "Attempting reverse theorem claims before forward sector completeness would decouple EOM/QG from declared L_total content.",
        "next_honest_step": "Complete sector export gates in P1834, then execute reverse theorem worklist in this unlock order.",
        "lay_explanation": "Najpierw trzeba domknąć pełny Lagrangian; dopiero wtedy odblokowują się końcowe testy QG.",
    }

    path = GEN / "p1835_s785_strict_forward_reverse_unlock_matrix_checkpoint.json"
    path.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(path)


if __name__ == "__main__":
    main()
