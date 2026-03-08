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
    / "f31_initial_future_strict_core_selector_source_seed_target_packet.json"
)
OUT_SUMMARY = (
    GENERATED
    / "f31_initial_future_strict_core_selector_source_seed_target_packet_summary.json"
)


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    n127 = load_json(
        "fundamental_action_reconstruction/generated/n127_last_positive_branch_reduced_to_one_minimal_future_source_object_target_theorem_summary.json"
    )

    target_chain = n127["theorem_result"]["future_source_object_target_chain"]

    checks_spec = [
        {
            "id": "n127_last_branch_reduced_to_one_target",
            "actual": n127["theorem_result"]["last_positive_branch_reduced_to_one_target"],
            "expected": True,
            "meaning": "N127 already reduces the last positive branch to one explicit target chain",
        },
        {
            "id": "target_begins_with_source_object",
            "actual": target_chain["strict_core_source_object"],
            "expected": "S_sel_int",
            "meaning": "the chain begins with the new strict-core source object",
        },
        {
            "id": "target_second_node_is_orientation_export",
            "actual": target_chain["internal_orientation_export"],
            "expected": "E_orient",
            "meaning": "the second node is the internal orientation export",
        },
        {
            "id": "downstream_nodes_remain_later",
            "actual": [
                target_chain["strict_core_bridge"],
                target_chain["selector_reduction"],
                target_chain["downstream_operator_reachability"],
            ],
            "expected": ["B_sel", "R_sel", "O_sel"],
            "meaning": "the bridge, reduction, and downstream operator nodes remain downstream of the source seed",
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
        "stage": "F31",
        "lane": "future_strict_core_selector_source_seed_target_only",
        "goal": "reduce_the_last_positive_branch_to_the_narrowest_forced_first_construction_subtarget",
        "status": "F31_EXECUTED_INITIAL_FUTURE_STRICT_CORE_SELECTOR_SOURCE_SEED_TARGET_PACKET_NO_FALSE_PASS",
        "initial_seed_target": {
            "strict_core_source_object": "S_sel_int",
            "internal_orientation_export": "E_orient",
        },
        "downstream_chain_left_open": {
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
        "stage": "F31",
        "status": artifact["status"],
        "lane": artifact["lane"],
        "initial_seed_target": artifact["initial_seed_target"],
        "downstream_chain_left_open": artifact["downstream_chain_left_open"],
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
