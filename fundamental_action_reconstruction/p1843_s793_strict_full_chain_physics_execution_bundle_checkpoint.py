#!/usr/bin/env python3
"""P1843 S793 strict full-chain physics execution bundle checkpoint."""

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
    p1841 = load("p1841_s791_strict_full_chain_next_action_program_checkpoint.json")
    p1842 = load("p1842_s792_strict_kernel_to_full_lagrangian_formula_pack_checkpoint.json")
    p1840 = load("p1840_s790_strict_bidirectional_closure_criteria_matrix_checkpoint.json")

    action_program = p1841.get("next_action_program", [])
    formula_pack = p1842.get("strict_formula_pack", {})
    closure_matrix = p1840.get("bidirectional_closure_matrix", {})

    # Build a compact execution bundle: sector term deliveries first, then gate-level theorem actions
    sector_bundle = [
        {
            "sector": item.get("sector", "UNKNOWN"),
            "density_symbol": item.get("density_symbol", "MISSING_DENSITY_SYMBOL"),
            "missing_exports": item.get("missing_exports", []),
            "required_output": "explicit nonproxy sector terms + covariant EOM residual trace",
        }
        for item in action_program
    ]

    theorem_bundle = []
    for gate, row in closure_matrix.items():
        theorem_bundle.append(
            {
                "gate": gate,
                "forward_ready": row.get("forward_ready", False),
                "reverse_unlocked": row.get("reverse_unlocked", False),
                "gate_ready": row.get("gate_ready", False),
                "required_output": "theorem witness object with PASS_ZERO or obstruction trace",
            }
        )

    out = {
        "packet_id": "P1843",
        "stage_id": "S793",
        "status": "OPEN_OBSTRUCTION_WITH_TRACE",
        "route": "strict_only",
        "kernel_anchor": formula_pack.get("kernel_anchor", "K_strict_anchor_missing"),
        "coefficient_map_schema": formula_pack.get("coefficient_map_schema", {}),
        "execution_bundle": {
            "phase_1_sector_term_delivery": sector_bundle,
            "phase_2_theorem_gate_witnesses": theorem_bundle,
        },
        "technical_progress": "Strict chain is now compiled into a physics execution bundle from kernel/coefficient scaffolding to theorem-gate witness outputs.",
        "proven": "The path kernel->coefficients->L_total->EOM->TG/QG is executable only as a two-phase bundle with sector terms before theorem witnesses.",
        "open": "Sector exports and theorem witnesses remain incomplete; strict-core closure is still open.",
        "false_pass_risk": "Skipping phase-1 sector term delivery and jumping to theorem gates would produce nonphysical closure claims.",
        "next_honest_step": "Deliver phase-1 sectors in listed order, then rerun closure matrix and execute phase-2 theorem witness exports.",
        "lay_explanation": "To plan zamknięcia teorii: najpierw kompletujemy składniki pełnego wzoru, potem dopiero robimy finalne dowody bramek QG.",
    }

    path = GEN / "p1843_s793_strict_full_chain_physics_execution_bundle_checkpoint.json"
    path.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(path)


if __name__ == "__main__":
    main()
