#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import UTC, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN1764 = GEN / "p1764_s714_strict_covariant_nonproxy_ea_eh_explicit_export_checkpoint.json"
IN1765 = GEN / "p1765_s715_strict_nonproxy_metric_elg_explicit_export_checkpoint.json"
IN1781 = GEN / "p1781_s731_strict_reverse_chain_execution_readiness_delta_checkpoint.json"
OUT = GEN / "p1782_s732_strict_priority_closure_gap_matrix_checkpoint.json"


def main() -> None:
    p1764 = json.loads(IN1764.read_text(encoding="utf-8"))
    p1765 = json.loads(IN1765.read_text(encoding="utf-8"))
    p1781 = json.loads(IN1781.read_text(encoding="utf-8"))

    payload = {
        "checkpoint": "P1782_S732",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "strict_only": True,
        "legacy_bridge_used": False,
        "input_anchors": ["p1764_s714", "p1765_s715", "p1781_s731"],
        "priority_closure_gap_matrix": {
            "E_A_mu_nonproxy": {
                "state": p1764.get("d1_d2_upgrade", {}).get("E_A_mu_nonproxy_explicit_v2", {}).get("componentwise_state", "OPEN"),
                "class": "NONPROXY_OPERATOR_EXPORT_DONE_COMPONENTWISE_OPEN",
            },
            "E_H_nonproxy": {
                "state": p1764.get("d1_d2_upgrade", {}).get("E_H_nonproxy_explicit_v2", {}).get("componentwise_state", "OPEN"),
                "class": "NONPROXY_OPERATOR_EXPORT_DONE_COMPONENTWISE_OPEN",
            },
            "EL_g_nonproxy": {
                "state": p1765.get("explicit_nonproxy_EL_g_munu_v1", {}).get("classification", {}).get("strict_local", "OPEN"),
                "class": "NONPROXY_OPERATOR_EXPORT_DONE_COMPONENTWISE_OPEN",
            },
            "H1_weak_form": {"state": "UPGRADED_OPERATOR_LEVEL", "class": "WEAK_FORM_READY_STRICT_LOCAL_OPEN"},
            "Bianchi_Ward": {"state": "BLOCKED", "class": "PENDING_JOINT_COMPONENTWISE_WITNESS"},
            "BRST_Cutkosky": {"state": "BLOCKED", "class": "PENDING_BW_PASS_AND_THEOREM_WITNESSES"},
        },
        "readiness_link": p1781.get("reverse_chain_readiness_delta", {}),
        "proofs_and_nonproofs": {
            "proven_now": [
                "Priority closure gaps are explicitly classified with NONPROXY/LOCAL/GLOBAL semantics and no false pass.",
            ],
            "still_open": [
                "Componentwise H1 and EL_g-E_munu residual witnesses.",
                "Theorem-level BW/BRST/Cutkosky closures.",
            ],
        },
        "status": "KEEP_OPEN_QG_THEOREM_LEVEL_REQUIRED",
        "no_false_pass_claim": True,
        "next_honest_step": "Dowieźć W1 FULL_EXPORT, wykonać joint componentwise witness run, a następnie zaktualizować tylko te klasy, które mają jawny PASS_ZERO albo OBSTRUCTION_WITH_DIVERGENCE_TRACE.",
        "lay_summary": "Mamy tabelę braków zamknięcia: wiemy dokładnie, które elementy są już jawnie zapisane, a które nadal czekają na twarde testy.",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
