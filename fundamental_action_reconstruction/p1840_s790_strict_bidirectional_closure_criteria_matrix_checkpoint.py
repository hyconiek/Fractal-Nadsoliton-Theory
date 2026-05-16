#!/usr/bin/env python3
"""P1840 S790 strict bidirectional closure criteria matrix checkpoint."""

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
    p1839 = load("p1839_s789_strict_full_lagrangian_term_delivery_contract_checkpoint.json")
    p1838 = load("p1838_s788_strict_kernel_to_eom_to_qg_theorem_gate_map_checkpoint.json")
    p1835 = load("p1835_s785_strict_forward_reverse_unlock_matrix_checkpoint.json")

    sectors = p1839.get("term_delivery_contract", [])
    tg_map = p1838.get("theorem_gate_map", {})
    unlock = p1835.get("reverse_unlock_matrix", [])

    sector_ready_map = {s.get("sector"): bool(s.get("sector_ready", False)) for s in sectors}
    reverse_unlock_map = {u.get("theorem_gate"): u.get("unlock_state", "LOCKED_BY_FORWARD_CHAIN") for u in unlock}

    matrix = {}
    for gate, req in tg_map.items():
        forward_req = req.get("requires_forward", [])
        reverse_req = req.get("requires_reverse", [])

        forward_ok = all(sector_ready_map.get(s, False) for s in forward_req)
        reverse_ok = all(reverse_unlock_map.get(r, "LOCKED_BY_FORWARD_CHAIN") == "UNLOCKED_FOR_EXECUTION" for r in reverse_req)

        matrix[gate] = {
            "forward_requirements": forward_req,
            "forward_ready": forward_ok,
            "reverse_requirements": reverse_req,
            "reverse_unlocked": reverse_ok,
            "gate_ready": forward_ok and reverse_ok,
        }

    strict_core_ready = all(v.get("gate_ready", False) for v in matrix.values()) and bool(matrix)

    out = {
        "packet_id": "P1840",
        "stage_id": "S790",
        "status": "PASS_ZERO" if strict_core_ready else "OPEN_OBSTRUCTION_WITH_TRACE",
        "route": "strict_only",
        "bidirectional_closure_matrix": matrix,
        "strict_core_ready": strict_core_ready,
        "technical_progress": "Forward term-delivery and reverse unlock signals are fused into a single strict-core closure criteria matrix.",
        "proven": "A theorem gate is ready only when both forward sector content and reverse unlock requirements are jointly satisfied.",
        "open": "At least one TG criterion remains unsatisfied, so strict-core closure is not achieved.",
        "false_pass_risk": "Evaluating TG readiness without bidirectional criteria matrix can yield one-sided false-pass claims.",
        "next_honest_step": "Satisfy missing forward sector deliveries and reverse unlock states per gate row, then rerun P1840.",
        "lay_explanation": "To tabela warunków domknięcia: każda bramka potrzebuje jednocześnie kompletu wzoru i odblokowanych dowodów QG.",
    }

    path = GEN / "p1840_s790_strict_bidirectional_closure_criteria_matrix_checkpoint.json"
    path.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(path)


if __name__ == "__main__":
    main()
