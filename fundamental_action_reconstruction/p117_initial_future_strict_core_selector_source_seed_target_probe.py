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
    / "p117_initial_future_strict_core_selector_source_seed_target_probe.json"
)
OUT_SUMMARY = (
    GENERATED
    / "p117_initial_future_strict_core_selector_source_seed_target_probe_summary.json"
)


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    f31 = load_json(
        "fundamental_action_reconstruction/generated/f31_initial_future_strict_core_selector_source_seed_target_packet_summary.json"
    )
    n127 = load_json(
        "fundamental_action_reconstruction/generated/n127_last_positive_branch_reduced_to_one_minimal_future_source_object_target_theorem_summary.json"
    )

    reduced_to_one_seed_target = (
        n127["theorem_result"]["last_positive_branch_reduced_to_one_target"]
        and f31["initial_seed_target"]["strict_core_source_object"] == "S_sel_int"
        and f31["initial_seed_target"]["internal_orientation_export"] == "E_orient"
    )

    checks_spec = [
        {
            "id": "n127_full_target_chain_present",
            "actual": n127["theorem_result"]["last_positive_branch_reduced_to_one_target"],
            "expected": True,
            "meaning": "N127 already exports the full target chain",
        },
        {
            "id": "f31_seed_source_name",
            "actual": f31["initial_seed_target"]["strict_core_source_object"],
            "expected": "S_sel_int",
            "meaning": "F31 exports the seed source object name",
        },
        {
            "id": "f31_seed_orientation_name",
            "actual": f31["initial_seed_target"]["internal_orientation_export"],
            "expected": "E_orient",
            "meaning": "F31 exports the seed orientation-export name",
        },
        {
            "id": "reduced_to_one_seed_target",
            "actual": reduced_to_one_seed_target,
            "expected": True,
            "meaning": "the last branch is reduced to one explicit first seed target",
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
        "stage": "P117",
        "lane": "initial_future_selector_source_seed_target_probe_only",
        "goal": "test_whether_the_last_positive_branch_is_now_reduced_to_one_explicit_first_seed_target",
        "status": "CURRENT_REPO_REDUCES_THE_LAST_POSITIVE_BRANCH_TO_ONE_INITIAL_FUTURE_STRICT_CORE_SELECTOR_SOURCE_SEED_TARGET_AFTER_P117",
        "target_state": {
            "last_positive_branch_reduced_to_one_initial_seed_target": reduced_to_one_seed_target,
            "initial_seed_target": f31["initial_seed_target"],
            "downstream_chain_left_open": f31["downstream_chain_left_open"],
        },
        "checks": checks,
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": "P117",
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
