#!/usr/bin/env python3
"""P1827 S777 strict S1 target execution order checkpoint."""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"


def load(name: str) -> dict:
    return json.loads((GEN / name).read_text(encoding="utf-8"))


def main() -> None:
    p1820 = load("p1820_s770_strict_current_priority_bottleneck_to_execution_contract_checkpoint.json")
    p1824 = load("p1824_s774_strict_s1_missing_witness_worklist_checkpoint.json")
    p1825 = load("p1825_s775_strict_s1_worklist_delivery_tracker_checkpoint.json")

    ranking = p1820.get("bottleneck_ranking", [])
    missing_targets = {w.get("target") for w in p1824.get("s1_missing_witness_worklist", [])}

    delivery = {d.get("target"): d for d in p1825.get("delivery_tracker", [])}

    ordered_targets = []
    for row in ranking:
        item = row.get("item")
        if item in missing_targets:
            d = delivery.get(item, {})
            artifacts = d.get("artifact_delivery", [])
            delivered = sum(1 for a in artifacts if a.get("status", "").startswith("PASS"))
            total = len(artifacts)
            ordered_targets.append(
                {
                    "target": item,
                    "severity": row.get("severity", 0),
                    "delivery_progress": {
                        "delivered": delivered,
                        "total": total,
                    },
                    "next_required_action": "deliver_next_OPEN_UNDELIVERED_artifact",
                }
            )

    out = {
        "packet_id": "P1827",
        "stage_id": "S777",
        "status": "OPEN_OBSTRUCTION_WITH_TRACE",
        "route": "strict_only",
        "ordered_s1_execution_targets": ordered_targets,
        "technical_progress": "S1 missing targets are now ordered by strict bottleneck severity and delivery progress.",
        "proven": "Execution order is deterministic and tied to current strict priority ranking.",
        "open": "All ordered targets remain unfinished until artifact delivery fields move beyond OPEN_UNDELIVERED.",
        "false_pass_risk": "Unordered parallel claims can hide unresolved high-severity targets and create misleading readiness narratives.",
        "next_honest_step": "Execute targets in listed order, update P1825 delivery statuses, then rerun P1826/P1823.",
        "lay_explanation": "To kolejka prac: najpierw robimy najważniejsze brakujące cele i śledzimy postęp artefakt po artefakcie.",
    }

    path = GEN / "p1827_s777_strict_s1_target_execution_order_checkpoint.json"
    path.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(path)


if __name__ == "__main__":
    main()
