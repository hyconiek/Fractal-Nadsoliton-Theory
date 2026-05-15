#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import UTC, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN1760 = GEN / "p1760_s710_strict_shared_background_family_contract_delivery_checkpoint.json"
OUT = GEN / "p1761_s711_strict_boundary_term_control_clause_delivery_checkpoint.json"


def main() -> None:
    p1760 = json.loads(IN1760.read_text(encoding="utf-8"))

    clause = {
        "clause_name": "boundary_term_control_clause_v1",
        "formal_statement": "Boundary contributions from integration-by-parts in H1 cross-variation must be explicitly exported and constrained before any weak-form PASS claim.",
        "required_exports": [
            "surface_term_symbolic_form",
            "boundary_family_assumption",
            "vanishing_or_compensated_boundary_condition",
        ],
        "promotion_rule": "weak-form pass alone cannot promote strict-core closure",
        "status": "EXPORTED_CONTRACT_TEMPLATE",
    }

    payload = {
        "checkpoint": "P1761_S711",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "strict_only": True,
        "legacy_bridge_used": False,
        "chain": "K_strict -> coefficients -> full non-skeleton L_total -> EOM -> D4 delivery (boundary term control clause)",
        "input_anchor": "p1760_s710",
        "d4_delivery": clause,
        "execution_sequence_status": {
            "D1_export_E_A_mu_nonproxy": p1760.get("execution_sequence_status", {}).get("D1_export_E_A_mu_nonproxy", "PARTIAL_TEMPLATE_DELIVERED"),
            "D2_export_E_H_nonproxy": p1760.get("execution_sequence_status", {}).get("D2_export_E_H_nonproxy", "PARTIAL_TEMPLATE_DELIVERED"),
            "D3_export_shared_background_family_contract": p1760.get("execution_sequence_status", {}).get("D3_export_shared_background_family_contract", "TEMPLATE_DELIVERED"),
            "D4_export_boundary_term_control_clause": "TEMPLATE_DELIVERED",
            "D5_finalize_boundary_control_contract": "PENDING",
            "D6_run_nonproxy_H1_4D": "BLOCKED",
        },
        "status": "KEEP_OPEN_QG_THEOREM_LEVEL_REQUIRED",
        "no_false_pass_claim": True,
        "next_honest_step": "Dowieźć D5: final boundary-control contract spinający D3+D4, a następnie uruchomić D6 (nonproxy H1 4D).",
        "lay_summary": "Dopisaliśmy twarde reguły dla wyrazów brzegowych. To zabezpiecza przed sztucznym PASS w sytuacji, gdy wynik zależy od granic obszaru.",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
