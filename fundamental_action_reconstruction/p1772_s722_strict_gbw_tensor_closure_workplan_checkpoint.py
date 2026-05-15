#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import UTC, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN1771 = GEN / "p1771_s721_strict_gbw_componentwise_retry_after_sanitization_checkpoint.json"
OUT = GEN / "p1772_s722_strict_gbw_tensor_closure_workplan_checkpoint.json"


def main() -> None:
    p1771 = json.loads(IN1771.read_text(encoding="utf-8"))

    payload = {
        "checkpoint": "P1772_S722",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "strict_only": True,
        "legacy_bridge_used": False,
        "input_anchor": "p1771_s721",
        "objective": "Close remaining tensor-expansion blockers for G_BW without changing background family or verdict policy.",
        "tensor_closure_workplan": {
            "W1_H_R2_componentwise": {
                "target": "Expand full componentwise form of H_R2 and its divergence contribution.",
                "deliverable": "explicit_H_R2_component_basis_export",
                "status": "OPEN",
            },
            "W2_H_Ric2_componentwise": {
                "target": "Resolve mixed-covariant derivative commutator terms for H_Ric2.",
                "deliverable": "explicit_H_Ric2_component_basis_export",
                "status": "OPEN",
            },
            "W3_H_Riem2_componentwise": {
                "target": "Resolve double-divergence transport and index contractions in H_Riem2.",
                "deliverable": "explicit_H_Riem2_component_basis_export",
                "status": "OPEN",
            },
            "W4_CT_tensor_basis": {
                "target": "Export full componentwise CT tensor basis compatible with B1/B2/B3/C1/C2.",
                "deliverable": "explicit_T_CT_component_basis_export",
                "status": "OPEN",
            },
        },
        "gbw_reexecution_contract": {
            "fixed_background_family": True,
            "fixed_residual_basis": ["B1", "B2", "B3", "C1", "C2"],
            "fixed_verdict_policy": ["PASS_ZERO", "OBSTRUCTION_WITH_DIVERGENCE_TRACE"],
            "forbidden_actions": [
                "No BRST/Cutkosky gate promotion before G_BW verdict.",
                "No reduced-sector substitution for nonproxy closure.",
            ],
        },
        "status": "KEEP_OPEN_QG_THEOREM_LEVEL_REQUIRED",
        "no_false_pass_claim": True,
        "next_honest_step": "Dowieźć W1-W4 i dopiero wtedy uruchomić kolejne wykonanie G_BW z pełnym residual/divergence trace.",
        "lay_summary": "Mamy już listę dokładnych braków matematycznych. Następny krok to je po kolei uzupełnić i dopiero potem powtórzyć test końcowy.",
        "previous_retry_context": p1771.get("componentwise_retry_result", {}),
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
