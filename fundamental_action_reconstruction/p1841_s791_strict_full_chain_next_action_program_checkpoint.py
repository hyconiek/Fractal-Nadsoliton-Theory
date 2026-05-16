#!/usr/bin/env python3
"""P1841 S791 strict full-chain next-action program checkpoint."""

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
    p1840 = load("p1840_s790_strict_bidirectional_closure_criteria_matrix_checkpoint.json")

    terms = p1839.get("term_delivery_contract", [])
    matrix = p1840.get("bidirectional_closure_matrix", {})

    high_value_order = ["gravity_sector", "gauge_sector", "fermion_sector", "higgs_sector", "interaction_terms", "covariant_structures", "scalar_mix_sector"]
    term_map = {t.get("sector"): t for t in terms}

    program = []
    for sector in high_value_order:
        t = term_map.get(sector)
        if not t:
            continue
        if t.get("sector_ready", False):
            continue
        program.append(
            {
                "sector": sector,
                "density_symbol": t.get("density_symbol", "MISSING_DENSITY_SYMBOL"),
                "missing_exports": t.get("missing_exports", []),
                "completion_target": "sector_ready == true",
            }
        )

    out = {
        "packet_id": "P1841",
        "stage_id": "S791",
        "status": "OPEN_OBSTRUCTION_WITH_TRACE",
        "route": "strict_only",
        "next_action_program": program,
        "gate_readiness_snapshot": matrix,
        "strict_core_ready": p1840.get("strict_core_ready", False),
        "technical_progress": "A concrete sector-by-sector next-action program is extracted from full-Lagrangian term contract and TG readiness matrix.",
        "proven": "Current strict-core closure remains blocked by incomplete sector exports and unsatisfied gate rows.",
        "open": "Program remains non-empty; theorem gates cannot be promoted yet.",
        "false_pass_risk": "Skipping sector completion program and jumping to gate evaluation would break strict bidirectional closure discipline.",
        "next_honest_step": "Execute program in listed order and rerun P1839->P1840 after each sector delivery batch.",
        "lay_explanation": "To konkretna lista kolejnych ruchów: domykamy sektory jeden po drugim, a po każdym kroku sprawdzamy czy bramki teorii się odblokowały.",
    }

    path = GEN / "p1841_s791_strict_full_chain_next_action_program_checkpoint.json"
    path.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(path)


if __name__ == "__main__":
    main()
