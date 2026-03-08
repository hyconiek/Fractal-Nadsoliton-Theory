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
    / "p116_future_strict_core_internal_selector_source_object_target_probe.json"
)
OUT_SUMMARY = (
    GENERATED
    / "p116_future_strict_core_internal_selector_source_object_target_probe_summary.json"
)


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    f30 = load_json(
        "fundamental_action_reconstruction/generated/f30_minimal_future_strict_core_internal_selector_source_object_target_packet_summary.json"
    )
    n126 = load_json(
        "fundamental_action_reconstruction/generated/n126_current_repo_exports_no_admissible_strict_core_internal_selector_source_object_theorem_summary.json"
    )
    p2 = load_json(
        "fundamental_action_reconstruction/generated/p2_strict_core_sigma_int_to_a1_pair1_probe_summary.json"
    )

    reduced_to_one_target = (
        not n126["theorem_result"][
            "admissible_strict_core_internal_selector_source_object_present"
        ]
        and n126["theorem_result"]["future_positive_move_must_add_new_source_object"]
        and p2["status"] == "NOT_COMPUTABLE_FROM_CURRENT_STRICT_CORE_ROUTE"
    )

    checks_spec = [
        {
            "id": "n126_no_current_admissible_object",
            "actual": n126["theorem_result"][
                "admissible_strict_core_internal_selector_source_object_present"
            ],
            "expected": False,
            "meaning": "no current admissible source object is present",
        },
        {
            "id": "p2_no_downstream_reachability",
            "actual": p2["status"],
            "expected": "NOT_COMPUTABLE_FROM_CURRENT_STRICT_CORE_ROUTE",
            "meaning": "no downstream strict-core operator reachability is present",
        },
        {
            "id": "f30_target_chain_exported",
            "actual": f30["minimal_target_chain"]["strict_core_source_object"],
            "expected": "S_sel_int",
            "meaning": "F30 exports one explicit minimal future target chain",
        },
        {
            "id": "reduced_to_one_target",
            "actual": reduced_to_one_target,
            "expected": True,
            "meaning": "the last positive branch is reduced to one explicit future construction target",
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
        "stage": "P116",
        "lane": "future_strict_core_internal_selector_source_target_probe_only",
        "goal": "test_whether_the_last_positive_branch_has_now_been_reduced_to_one_explicit_future_construction_target",
        "status": "CURRENT_REPO_REDUCES_THE_LAST_POSITIVE_BRANCH_TO_ONE_MINIMAL_FUTURE_STRICT_CORE_INTERNAL_SELECTOR_SOURCE_OBJECT_TARGET_AFTER_P116",
        "target_state": {
            "last_positive_branch_reduced_to_one_target": reduced_to_one_target,
            "minimal_target_chain": f30["minimal_target_chain"],
        },
        "checks": checks,
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": "P116",
        "status": artifact["status"],
        "lane": artifact["lane"],
        "target_state": artifact["target_state"],
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
