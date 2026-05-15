#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import UTC, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN1769 = GEN / "p1769_s719_strict_gbw_componentwise_first_execution_attempt_checkpoint.json"
OUT = GEN / "p1770_s720_strict_curvature_tensor_notation_sanitization_and_gbw_retry_contract_checkpoint.json"


def main() -> None:
    p1769 = json.loads(IN1769.read_text(encoding="utf-8"))

    sanitized_index_lock = {
        "H_R2_munu": "H^(R2)_{mu nu} = 2 R R_{mu nu} - (1/2) g_{mu nu} R^2 + 2(g_{mu nu} Box - nabla_mu nabla_nu)R",
        "H_Ric2_munu": "H^(Ric2)_{mu nu} = 2 R_{mu alpha}R_{nu}^alpha - (1/2)g_{mu nu}R_{alpha beta}R^{alpha beta} + Box R_{mu nu} + g_{mu nu} nabla_alpha nabla_beta R^{alpha beta} - 2 nabla_alpha nabla_(mu R_{nu)}^alpha",
        "H_Riem2_munu": "H^(Riem2)_{mu nu} = 2 R_{mu alpha beta gamma}R_{nu}^{alpha beta gamma} - (1/2)g_{mu nu}R_{alpha beta gamma delta}R^{alpha beta gamma delta} - 4 nabla_alpha nabla_beta R_{mu nu}^{alpha beta}",
    }

    payload = {
        "checkpoint": "P1770_S720",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "strict_only": True,
        "legacy_bridge_used": False,
        "input_anchor": "p1769_s719",
        "issue_detected": "Curvature-tensor notation encoding ambiguity in componentwise lock string (risk for parser/backends).",
        "sanitized_curvature_tensor_index_lock": sanitized_index_lock,
        "gbw_retry_contract": {
            "precondition": [
                "Use sanitized curvature tensor notation (ASCII-safe, parser-safe).",
                "Keep same background family and same residual basis B1/B2/B3/C1/C2.",
                "Keep no-false-pass verdict policy unchanged.",
            ],
            "must_output": [
                "componentwise_EL_g_minus_E_munu_B1_B2_B3_C1_C2",
                "componentwise_nabla_mu_E_total_munu_trace",
                "verdict in {PASS_ZERO, OBSTRUCTION_WITH_DIVERGENCE_TRACE}",
            ],
        },
        "gate_update": {
            "G_BW": "RETRY_REQUIRED_AFTER_NOTATION_SANITIZATION",
            "G_BRST": "BLOCKED_UNTIL_G_BW_FINAL_VERDICT",
            "G_CUT": "BLOCKED_UNTIL_G_BW_AND_G_BRST",
        },
        "status": "KEEP_OPEN_QG_THEOREM_LEVEL_REQUIRED",
        "no_false_pass_claim": True,
        "next_honest_step": "Wykonać ponowną próbę G_BW z oczyszczoną notacją tensorową i opublikować pełny trace residual/divergence bez zmiany kryteriów werdyktu.",
        "lay_summary": "Najpierw naprawiamy zapis matematyczny tak, by narzędzia liczące rozumiały go jednoznacznie. Potem dopiero powtarzamy test i sprawdzamy, czy blokada znika.",
        "previous_obstruction_context": p1769.get("componentwise_execution_result", {}),
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
