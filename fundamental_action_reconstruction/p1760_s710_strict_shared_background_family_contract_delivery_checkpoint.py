#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import UTC, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN1759 = GEN / "p1759_s709_strict_nonproxy_e_h_template_delivery_checkpoint.json"
OUT = GEN / "p1760_s710_strict_shared_background_family_contract_delivery_checkpoint.json"


def main() -> None:
    p1759 = json.loads(IN1759.read_text(encoding="utf-8"))

    contract = {
        "contract_name": "shared_background_family_contract_v1",
        "scope": "common 4D background family for E_A^mu and E_H nonproxy runs",
        "required_fields": [
            "metric_family_label",
            "chart_and_index_convention",
            "boundary_hypothesis",
            "field_regularization_domain",
        ],
        "status": "EXPORTED_CONTRACT_TEMPLATE",
    }

    payload = {
        "checkpoint": "P1760_S710",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "strict_only": True,
        "legacy_bridge_used": False,
        "chain": "K_strict -> coefficients -> full non-skeleton L_total -> EOM -> D3 delivery (shared background contract)",
        "input_anchor": "p1759_s709",
        "d3_delivery": contract,
        "execution_sequence_status": {
            "D1_export_E_A_mu_nonproxy": p1759.get("execution_sequence_status", {}).get("D1_export_E_A_mu_nonproxy", "PARTIAL_TEMPLATE_DELIVERED"),
            "D2_export_E_H_nonproxy": p1759.get("execution_sequence_status", {}).get("D2_export_E_H_nonproxy", "PARTIAL_TEMPLATE_DELIVERED"),
            "D3_export_shared_background_family_contract": "TEMPLATE_DELIVERED",
            "D4_export_boundary_term_control_clause": "PENDING",
            "D5_finalize_boundary_control_contract": "PENDING",
            "D6_run_nonproxy_H1_4D": "BLOCKED",
        },
        "status": "KEEP_OPEN_QG_THEOREM_LEVEL_REQUIRED",
        "no_false_pass_claim": True,
        "next_honest_step": "Dowieźć D4 (boundary term control clause) i D5 (boundary control contract), potem uruchomić D6.",
        "lay_summary": "Ustaliliśmy wspólne ramy tła obliczeń dla obu równań. To warunek, by test 4D był porównywalny i uczciwy.",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
