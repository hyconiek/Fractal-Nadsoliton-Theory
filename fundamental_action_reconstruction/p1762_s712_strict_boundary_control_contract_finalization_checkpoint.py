#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import UTC, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN1761 = GEN / "p1761_s711_strict_boundary_term_control_clause_delivery_checkpoint.json"
OUT = GEN / "p1762_s712_strict_boundary_control_contract_finalization_checkpoint.json"


def main() -> None:
    p1761 = json.loads(IN1761.read_text(encoding="utf-8"))

    final_contract = {
        "contract_name": "boundary_control_contract_v1",
        "binds": [
            "shared_background_family_contract_v1",
            "boundary_term_control_clause_v1",
            "E_A_mu_nonproxy_template_v1",
            "E_H_nonproxy_template_v1",
        ],
        "admissible_verdicts": {
            "strict_local": ["PASS_ZERO", "OBSTRUCTION"],
            "weak_form": ["PASS_WEAK_FORM_WITH_BOUNDARY_CLAUSE", "OBSTRUCTION"],
        },
        "promotion_rule": "STRICT_CORE_CLOSURE cannot be claimed from weak_form pass alone",
        "status": "FINALIZED_CONTRACT_TEMPLATE",
    }

    payload = {
        "checkpoint": "P1762_S712",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "strict_only": True,
        "legacy_bridge_used": False,
        "chain": "K_strict -> coefficients -> full non-skeleton L_total -> EOM -> D5 final boundary-control contract",
        "input_anchor": "p1761_s711",
        "d5_delivery": final_contract,
        "execution_sequence_status": {
            "D1_export_E_A_mu_nonproxy": p1761.get("execution_sequence_status", {}).get("D1_export_E_A_mu_nonproxy", "PARTIAL_TEMPLATE_DELIVERED"),
            "D2_export_E_H_nonproxy": p1761.get("execution_sequence_status", {}).get("D2_export_E_H_nonproxy", "PARTIAL_TEMPLATE_DELIVERED"),
            "D3_export_shared_background_family_contract": p1761.get("execution_sequence_status", {}).get("D3_export_shared_background_family_contract", "TEMPLATE_DELIVERED"),
            "D4_export_boundary_term_control_clause": p1761.get("execution_sequence_status", {}).get("D4_export_boundary_term_control_clause", "TEMPLATE_DELIVERED"),
            "D5_finalize_boundary_control_contract": "TEMPLATE_DELIVERED",
            "D6_run_nonproxy_H1_4D": "READY_FOR_FIRST_EXECUTION_ATTEMPT",
        },
        "status": "KEEP_OPEN_QG_THEOREM_LEVEL_REQUIRED",
        "no_false_pass_claim": True,
        "next_honest_step": "Uruchomić D6: pierwszy nonproxy H1(A_mu,H) 4D z dualnym raportem strict_local + weak_form i werdyktem PASS_ZERO albo OBSTRUCTION.",
        "lay_summary": "Domknęliśmy zasady bezpieczeństwa dla testu 4D. Teraz można wykonać pierwszy pełny bieg testu bez ryzyka metodologicznych skrótów.",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
