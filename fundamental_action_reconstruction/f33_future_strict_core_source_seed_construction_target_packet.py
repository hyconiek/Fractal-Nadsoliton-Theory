#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT_JSON = (
    GENERATED / "f33_future_strict_core_source_seed_construction_target_packet.json"
)
OUT_SUMMARY = (
    GENERATED
    / "f33_future_strict_core_source_seed_construction_target_packet_summary.json"
)


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    n129 = load_json(
        "fundamental_action_reconstruction/generated/n129_last_positive_branch_reduced_to_one_initial_source_plus_orientation_package_theorem_summary.json"
    )

    checks_spec = [
        {
            "id": "n129_initial_package_reduced",
            "actual": n129["theorem_result"][
                "last_positive_branch_reduced_to_one_initial_package"
            ],
            "expected": True,
            "meaning": "N129 already reduces the last branch to one initial package",
        },
        {
            "id": "n129_first_remaining_branch_is_source_construction",
            "actual": n129["remaining_open_branches"][0],
            "expected": "future_construction_of_admissible_S_sel_int",
            "meaning": "the first remaining open branch is future source construction",
        },
        {
            "id": "n129_source_name",
            "actual": n129["theorem_result"]["initial_package"]["strict_core_source_object"],
            "expected": "S_sel_int",
            "meaning": "the initial package is anchored on the future source object",
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

    source_seed_construction_target = {
        "strict_core_source_object": "S_sel_int",
        "genuinely_new_strict_core_object_required": True,
        "source_carrying_enough_for_later_E_orient_export_required": True,
        "orientation_not_yet_counted_as_exported": True,
        "strict_core_only_required": True,
        "silent_legacy_to_strict_substitution_forbidden": True,
        "selector_acceptance_outside_strict_core_may_not_count_as_source_seed": True,
    }

    artifact = {
        "stage": "F33",
        "lane": "future_strict_core_source_seed_construction_target_only",
        "goal": "freeze_the_narrowest_first_source_construction_target_inside_the_last_positive_branch",
        "status": "F33_EXECUTED_FUTURE_STRICT_CORE_SOURCE_SEED_CONSTRUCTION_TARGET_PACKET_NO_FALSE_PASS",
        "source_seed_construction_target": source_seed_construction_target,
        "checks": checks,
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": "F33",
        "status": artifact["status"],
        "lane": artifact["lane"],
        "source_seed_construction_target": source_seed_construction_target,
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
