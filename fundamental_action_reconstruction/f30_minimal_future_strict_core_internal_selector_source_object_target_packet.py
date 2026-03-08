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
    / "f30_minimal_future_strict_core_internal_selector_source_object_target_packet.json"
)
OUT_SUMMARY = (
    GENERATED
    / "f30_minimal_future_strict_core_internal_selector_source_object_target_packet_summary.json"
)


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    n126 = load_json(
        "fundamental_action_reconstruction/generated/n126_current_repo_exports_no_admissible_strict_core_internal_selector_source_object_theorem_summary.json"
    )
    p2 = load_json(
        "fundamental_action_reconstruction/generated/p2_strict_core_sigma_int_to_a1_pair1_probe_summary.json"
    )

    checks_spec = [
        {
            "id": "n126_no_admissible_object_present",
            "actual": n126["theorem_result"][
                "admissible_strict_core_internal_selector_source_object_present"
            ],
            "expected": False,
            "meaning": "N126 already excludes all current objects from admissible source status",
        },
        {
            "id": "n126_future_move_must_add_new_source_object",
            "actual": n126["theorem_result"][
                "future_positive_move_must_add_new_source_object"
            ],
            "expected": True,
            "meaning": "N126 already says the next positive move must add a new source object",
        },
        {
            "id": "p2_no_downstream_operator_reachability",
            "actual": p2["status"],
            "expected": "NOT_COMPUTABLE_FROM_CURRENT_STRICT_CORE_ROUTE",
            "meaning": "P2 still keeps downstream strict-core operator reachability absent",
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
        "stage": "F30",
        "lane": "future_strict_core_internal_selector_source_target_only",
        "goal": "reduce_the_last_positive_branch_to_one_minimal_construction_target",
        "status": "F30_EXECUTED_MINIMAL_FUTURE_STRICT_CORE_INTERNAL_SELECTOR_SOURCE_OBJECT_TARGET_PACKET_NO_FALSE_PASS",
        "minimal_target_chain": {
            "strict_core_source_object": "S_sel_int",
            "internal_orientation_export": "E_orient",
            "strict_core_bridge": "B_sel",
            "selector_reduction": "R_sel",
            "downstream_operator_reachability": "O_sel",
        },
        "checks": checks,
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": "F30",
        "status": artifact["status"],
        "lane": artifact["lane"],
        "minimal_target_chain": artifact["minimal_target_chain"],
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
