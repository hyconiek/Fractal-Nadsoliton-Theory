#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import UTC, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN1739 = GEN / "p1739_s689_strict_full_chain_execution_frontier_and_minimal_nonproxy_delivery_checkpoint.json"
OUT = GEN / "p1740_s690_strict_kernel_to_eom_reverse_execution_trigger_contract_checkpoint.json"


def main() -> None:
    p1739 = json.loads(IN1739.read_text(encoding="utf-8"))
    frontier = p1739.get("execution_frontier", {})
    delivery = frontier.get("minimal_nonproxy_delivery", {})

    trigger_contract = {
        "trigger_T1_for_S2_H1": {
            "required_delivery": delivery.get("for_S2_H1", []),
            "on_delivery_action": "execute_H1_cross_variation_immediately",
            "result_policy": "PASS_ZERO_or_OBSTRUCTION_ONLY",
        },
        "trigger_T2_for_S3_metric": {
            "required_delivery": delivery.get("for_S3_metric", []),
            "on_delivery_action": "execute_metric_ELg_minus_Emunu_componentwise_immediately",
            "result_policy": "PASS_ZERO_or_OBSTRUCTION_ONLY",
        },
        "trigger_T3_for_S4_qg_update": {
            "required_preconditions": [
                "T1_result_published",
                "T2_result_published",
                "both_results_classified_PASS_ZERO_or_OBSTRUCTION",
            ],
            "on_preconditions_action": "allow_qg_gate_update",
            "result_policy": "NO_CLOSURE_UPGRADE_WITHOUT_THEOREM_WITNESS",
        },
    }

    payload = {
        "checkpoint": "P1740_S690",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "strict_only": True,
        "legacy_bridge_used": False,
        "route_target": "F_Nadsoliton => L_SM + L_GR (strict-only)",
        "chain": "K_strict -> coefficients -> full nonskeleton L_total -> EOM -> execution frontier -> trigger contract",
        "full_lagrangian_density_nonskeleton_instantiated": p1739.get(
            "full_lagrangian_density_nonskeleton_instantiated", {}
        ),
        "execution_frontier_anchor": frontier,
        "reverse_execution_trigger_contract": trigger_contract,
        "status": "KEEP_OPEN_QG_THEOREM_LEVEL_REQUIRED",
        "no_false_pass_claim": True,
        "next_honest_step": "Dowieźć required_delivery dla T1 i uruchomić H1 natychmiast; następnie analogicznie T2 dla metryki; T3 dopuścić dopiero po publikacji obu wyników i bez sztucznego upgrade closure.",
        "lay_summary": "To kontrakt uruchomień: gdy brakujące elementy są gotowe, obliczenia mają ruszyć od razu. Dzięki temu teoria nie stoi w miejscu i nie ma pozornych deklaracji sukcesu.",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
