#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT_JSON = (
    GENERATED / "p119_first_source_seed_construction_target_probe.json"
)
OUT_SUMMARY = (
    GENERATED / "p119_first_source_seed_construction_target_probe_summary.json"
)


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    n129 = load_json(
        "fundamental_action_reconstruction/generated/n129_last_positive_branch_reduced_to_one_initial_source_plus_orientation_package_theorem_summary.json"
    )
    f33 = load_json(
        "fundamental_action_reconstruction/generated/f33_future_strict_core_source_seed_construction_target_packet_summary.json"
    )

    reduced_to_one_first_source_seed_target = (
        n129["theorem_result"]["initial_package"]["strict_core_source_object"] == "S_sel_int"
        and n129["remaining_open_branches"][0] == "future_construction_of_admissible_S_sel_int"
        and f33["source_seed_construction_target"]["strict_core_source_object"] == "S_sel_int"
    )

    checks_spec = [
        {
            "id": "n129_source_construction_is_first_open_branch",
            "actual": n129["remaining_open_branches"][0],
            "expected": "future_construction_of_admissible_S_sel_int",
            "meaning": "source construction is the first remaining open branch",
        },
        {
            "id": "f33_source_seed_target_name",
            "actual": f33["source_seed_construction_target"]["strict_core_source_object"],
            "expected": "S_sel_int",
            "meaning": "the source-seed construction target is anchored on S_sel_int",
        },
        {
            "id": "reduced_to_one_first_source_seed_target",
            "actual": reduced_to_one_first_source_seed_target,
            "expected": True,
            "meaning": "the last positive branch is reduced to one explicit first source-seed construction target",
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
        "stage": "P119",
        "lane": "first_source_seed_construction_target_only",
        "goal": "test_whether_the_last_positive_branch_is_now_reduced_to_one_explicit_first_source_seed_construction_target",
        "status": "CURRENT_REPO_REDUCES_THE_LAST_POSITIVE_BRANCH_TO_ONE_FIRST_SOURCE_SEED_CONSTRUCTION_TARGET_AFTER_P119",
        "target_state": {
            "last_positive_branch_reduced_to_one_first_source_seed_target": reduced_to_one_first_source_seed_target,
            "first_source_seed_construction_target": f33["source_seed_construction_target"],
            "later_open_branches": [
                "future_derivation_of_admissible_E_orient_from_S_sel_int",
                "future_completion_of_B_sel_R_sel_O_sel_after_seed_package",
            ],
        },
        "checks": checks,
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": "P119",
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
