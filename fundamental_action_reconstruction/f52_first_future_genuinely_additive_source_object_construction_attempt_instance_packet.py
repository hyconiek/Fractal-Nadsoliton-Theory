#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT_JSON = (
    GENERATED
    / "f52_first_future_genuinely_additive_source_object_construction_attempt_instance_packet.json"
)
OUT_SUMMARY = (
    GENERATED
    / "f52_first_future_genuinely_additive_source_object_construction_attempt_instance_packet_summary.json"
)


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    f51 = load_json(
        "fundamental_action_reconstruction/generated/f51_first_future_genuinely_additive_source_object_construction_attempt_packet_summary.json"
    )
    n151 = load_json(
        "fundamental_action_reconstruction/generated/n151_only_remaining_positive_move_class_reduced_to_one_first_additive_construction_attempt_theorem_summary.json"
    )

    checks_spec = [
        {
            "id": "n151_remaining_positive_move_class",
            "actual": n151["theorem_result"]["only_remaining_positive_move_class"],
            "expected": "future_genuinely_additive_new_strict_core_source_object_construction",
            "meaning": "N151 keeps the same sole remaining positive move class",
        },
        {
            "id": "f51_target_id",
            "actual": f51["future_additive_attempt_target"]["target_id"],
            "expected": "S_sel_int_additive_attempt_target_v1",
            "meaning": "F51 fixes the target identity for the future attempt",
        },
        {
            "id": "f51_future_admissibility_target",
            "actual": f51["future_additive_attempt_target"]["future_admissibility_target"],
            "expected": True,
            "meaning": "F51 already keeps the target compatible with later admissibility work",
        },
    ]

    checks: list[dict[str, Any]] = []
    for item in checks_spec:
        checks.append(
            {
                "id": item["id"],
                "actual": item["actual"],
                "expected": item["expected"],
                "pass": item["actual"] == item["expected"],
                "meaning": item["meaning"],
            }
        )

    artifact = {
        "stage": "F52",
        "lane": "future_additive_construction_attempt_instance_only",
        "goal": "freeze_one_explicit_future_additive_construction_attempt_instance",
        "status": "F52_EXECUTED_FIRST_FUTURE_GENUINELY_ADDITIVE_SOURCE_OBJECT_CONSTRUCTION_ATTEMPT_INSTANCE_PACKET_NO_FALSE_PASS",
        "future_additive_attempt_instance": {
            "attempt_id": "construct_attempt_v1(S_sel_int_additive_attempt_target_v1)",
            "target_id": "S_sel_int_additive_attempt_target_v1",
            "constructed_source_object": False,
            "success_verdict_present": False,
            "failure_verdict_present": False,
            "admissible_S_sel_int_present": False,
        },
        "checks": checks,
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": "F52",
        "status": artifact["status"],
        "lane": artifact["lane"],
        "future_additive_attempt_instance": artifact["future_additive_attempt_instance"],
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    OUT_JSON.write_text(
        json.dumps(artifact, indent=2, ensure_ascii=True) + "\n",
        encoding="ascii",
    )
    OUT_SUMMARY.write_text(
        json.dumps(summary, indent=2, ensure_ascii=True) + "\n",
        encoding="ascii",
    )
    print(OUT_JSON)


if __name__ == "__main__":
    main()
