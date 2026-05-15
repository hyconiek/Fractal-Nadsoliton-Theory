#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import UTC, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN1738 = GEN / "p1738_s688_strict_kernel_to_full_lagrangian_to_eom_and_qg_closure_sequence_checkpoint.json"
OUT = GEN / "p1739_s689_strict_full_chain_execution_frontier_and_minimal_nonproxy_delivery_checkpoint.json"


def main() -> None:
    p1738 = json.loads(IN1738.read_text(encoding="utf-8"))
    seq = p1738.get("strict_closure_sequence", {})
    contract = p1738.get("next_action_contract_anchor", {})

    step1 = contract.get("step_1_h1_nonproxy_export", {})
    step2 = contract.get("step_2_metric_nonproxy_export", {})

    frontier = {
        "frontier_id": "F_EXEC_STRICT_REV_01",
        "priority_order": ["S2", "S3", "S4"],
        "minimal_nonproxy_delivery": {
            "for_S2_H1": step1.get("must_export", []),
            "for_S3_metric": step2.get("must_export", []),
        },
        "decision_policy": {
            "S2": step1.get("decision", "PASS_ZERO_or_OBSTRUCTION_ONLY"),
            "S3": step2.get("decision", "PASS_ZERO_or_OBSTRUCTION_ONLY"),
            "S4": "UPDATE_ONLY_AFTER_S2_S3_RESULTS",
        },
    }

    payload = {
        "checkpoint": "P1739_S689",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "strict_only": True,
        "legacy_bridge_used": False,
        "route_target": "F_Nadsoliton => L_SM + L_GR (strict-only)",
        "chain": "K_strict -> coefficients -> full nonskeleton L_total -> EOM -> execution frontier with minimal nonproxy delivery",
        "full_lagrangian_density_nonskeleton_instantiated": p1738.get(
            "full_lagrangian_density_nonskeleton_instantiated", {}
        ),
        "strict_closure_sequence_anchor": seq,
        "execution_frontier": frontier,
        "status": "KEEP_OPEN_QG_THEOREM_LEVEL_REQUIRED",
        "no_false_pass_claim": True,
        "next_honest_step": "Dowieźć minimal_nonproxy_delivery.for_S2_H1 i wykonać S2; po wyniku PASS_ZERO/OBSTRUCTION dowieźć minimal_nonproxy_delivery.for_S3_metric i wykonać S3, a następnie zaktualizować S4.",
        "lay_summary": "To precyzyjna lista minimalnych rzeczy do zrobienia teraz. Najpierw dostarczamy brakujące elementy do testu H1, potem do testu grawitacyjnego, i dopiero wtedy ruszamy duże bramki QG.",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
