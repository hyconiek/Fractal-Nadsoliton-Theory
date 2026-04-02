#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

AS_OF = "2026-03-27"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P1063 = GENERATED / "p1063_current_strict_t173_t176_post_f961_negative_3_cycle_reference_to_existing_t216_t218_pair12_provider_frontier_route_decision_audit_probe_summary.json"
OUT_JSON = GENERATED / "f962_current_strict_t173_t176_post_f961_negative_3_cycle_reference_to_existing_t216_t218_pair12_provider_frontier_route_decision_packet.json"
OUT_SUMMARY = GENERATED / "f962_current_strict_t173_t176_post_f961_negative_3_cycle_reference_to_existing_t216_t218_pair12_provider_frontier_route_decision_packet_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def main() -> None:
    if not IN_P1063.exists():
        artifact = {
            "stage": "F962",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": [str(IN_P1063.relative_to(REPO))],
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    p1063 = load_json(IN_P1063)
    admitted = bool(p1063.get("rejoin_to_existing_t216_t218_frontier") is True)
    status = (
        "F962_EXECUTED_CURRENT_STRICT_T173_T176_POST_F961_NEGATIVE_3_CYCLE_REFERENCE_TO_EXISTING_T216_T218_PAIR12_PROVIDER_FRONTIER_ROUTE_DECISION_PACKET_NO_FALSE_PASS"
        if admitted
        else "F962_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_ROUTE_DECISION_STATE"
    )

    artifact = {
        "stage": "F962",
        "status": status,
        "as_of": AS_OF,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "depends_on": {
            "p1063_route_decision_probe_summary": str(IN_P1063.relative_to(REPO)),
        },
        "route_decision": {
            "object_id": "PostF961Negative3CycleReferenceRejoinToExistingT216T218Pair12ProviderFrontier_v1",
            "negative_3_cycle_reference_stays_reference_only": True,
            "primary_continuation_target": "W_strict_t173_pair12_source_side_branch_selection_provider_actual_realization_attempt_v1",
            "exact_interface_frontier_target": "W_strict_t173_pair12_source_side_branch_selection_provider_actual_realization_attempt_chart_sensitive_pair12_typed_descent_interface_target_v1",
            "rejoin_to_existing_t216_t218_frontier": admitted,
        },
        "hard_limits": [
            "Does not discharge T216.",
            "Does not discharge T218.",
            "Does not discharge T183.",
            "Does not discharge T176.",
            "Does not discharge QW-2191.",
            "Does not claim ToE closure.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": artifact["stage"],
        "status": artifact["status"],
        "as_of": artifact["as_of"],
        "route_object_id": artifact["route_decision"]["object_id"],
        "negative_3_cycle_reference_stays_reference_only": True,
        "primary_continuation_target": artifact["route_decision"]["primary_continuation_target"],
        "exact_interface_frontier_target": artifact["route_decision"]["exact_interface_frontier_target"],
        "rejoin_to_existing_t216_t218_frontier": admitted,
        "no_false_pass": True,
    }

    write_json(OUT_JSON, artifact)
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
