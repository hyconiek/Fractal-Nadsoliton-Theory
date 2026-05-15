#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import UTC, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN1716 = GEN / "p1716_s666_strict_metric_index_convention_normalization_audit_checkpoint.json"
IN1727 = GEN / "p1727_s677_strict_componentwise_metric_residual_coefficient_vector_stub_checkpoint.json"
IN1767 = GEN / "p1767_s717_strict_bianchi_ward_to_brst_cutkosky_gate_sequencing_checkpoint.json"
OUT = GEN / "p1768_s718_strict_bw_componentwise_execution_input_pack_checkpoint.json"


def main() -> None:
    p1716 = json.loads(IN1716.read_text(encoding="utf-8"))
    p1727 = json.loads(IN1727.read_text(encoding="utf-8"))
    p1767 = json.loads(IN1767.read_text(encoding="utf-8"))

    payload = {
        "checkpoint": "P1768_S718",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "strict_only": True,
        "legacy_bridge_used": False,
        "goal": "Prepare concrete input pack for G_BW componentwise execution without symbolic placeholders",
        "input_anchors": ["p1716_s666", "p1727_s677", "p1767_s717"],
        "gbw_execution_input_pack": {
            "background_family": "shared_background_family_contract_v1",
            "index_sign_lock": p1716.get("normalized_h_terms", {}),
            "residual_basis": ["B1", "B2", "B3", "C1", "C2"],
            "coefficient_vector_slots": p1727.get("residual_vector_format", {}),
            "required_nonproxy_operators": [
                "E_A_mu_nonproxy_explicit_v2",
                "E_H_nonproxy_explicit_v2",
                "EL_g_nonproxy_explicit_v1",
            ],
            "required_outputs": {
                "metric_residual": "componentwise_EL_g_minus_E_munu_B1_B2_B3_C1_C2",
                "divergence_trace": "componentwise_nabla_mu_E_total_munu_trace",
                "verdict": ["PASS_ZERO", "OBSTRUCTION_WITH_DIVERGENCE_TRACE"],
            },
        },
        "gate_dependency": p1767.get("gate_sequencing_contract", {}),
        "status": "READY_FOR_COMPONENTWISE_EXECUTION_ATTEMPT",
        "no_false_pass_claim": True,
        "next_honest_step": "Uruchomić wykonanie G_BW na tym input-packu i opublikować wyłącznie wynik PASS_ZERO albo OBSTRUCTION_WITH_DIVERGENCE_TRACE z pełnym wektorem współczynników.",
        "lay_summary": "Spakowano wszystkie potrzebne dane do pierwszej pełnej próby testu zgodności Bianchi/Ward. Teraz trzeba wykonać obliczenie i uczciwie pokazać wynik.",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
