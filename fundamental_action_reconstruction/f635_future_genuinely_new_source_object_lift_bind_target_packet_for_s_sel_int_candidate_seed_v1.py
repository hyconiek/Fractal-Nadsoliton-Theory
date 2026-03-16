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
    / "f635_future_genuinely_new_source_object_lift_bind_target_packet_for_s_sel_int_candidate_seed_v1.json"
)
OUT_SUMMARY = (
    GENERATED
    / "f635_future_genuinely_new_source_object_lift_bind_target_packet_for_s_sel_int_candidate_seed_v1_summary.json"
)


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    n525 = load_json(
        "fundamental_action_reconstruction/generated/n525_current_first_clause_obstruction_theorem_for_s_sel_int_candidate_seed_v1_summary.json"
    )
    f318 = load_json(
        "fundamental_action_reconstruction/generated/f318_first_current_strict_side_source_seed_strict_sigma_int_upgrade_candidate_instance_packet_summary.json"
    )

    checks_spec = [
        {
            "id": "n525_first_clause_unsatisfied",
            "actual": n525["theorem_result"]["first_clause_currently_satisfied"],
            "expected": False,
            "meaning": "N525 blocks the upgraded seed on the first clause",
        },
        {
            "id": "f318_precursor_topological_layer",
            "actual": f318["candidate_seed_instance"]["construction_route"][
                "local_topological_protection_layer"
            ],
            "expected": "QW-2206_local_topological_protection_layer",
            "meaning": "the precursor topological layer remains fixed",
        },
        {
            "id": "f318_precursor_sigma_input",
            "actual": f318["candidate_seed_instance"]["construction_route"][
                "internal_binary_candidate"
            ],
            "expected": "sigma_int_strict_derived_v1",
            "meaning": "the precursor internal datum is the strict sigma-int source-upgrade object",
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

    future_lift_bind_target = {
        "future_target_name": "S_sel_int_new_object_target_v1",
        "construction_shape": "strict_core_single_object_lift_bind",
        "input_materials": [
            "QW-2206_local_topological_protection_layer",
            "sigma_int_strict_derived_v1",
        ],
        "target_output": "future_S_sel_int",
        "counts_only_as": "future_genuinely_new_source_object_lift_bind_target",
        "does_not_count_as": [
            "constructed_source_object",
            "admissible_S_sel_int",
            "admissible_E_orient",
            "strict_core_selector_closure",
        ],
    }

    artifact = {
        "stage": "F635",
        "lane": "future_genuinely_new_source_object_lift_bind_target_only",
        "goal": "freeze_the_narrowest_future_target_for_recovering_the_first_admissibility_clause_beyond_the_failed_seed_v1_attempt",
        "status": "F635_EXECUTED_FUTURE_GENUINELY_NEW_SOURCE_OBJECT_LIFT_BIND_TARGET_PACKET_FOR_S_SEL_INT_CANDIDATE_SEED_V1_NO_FALSE_PASS",
        "future_lift_bind_target": future_lift_bind_target,
        "checks": checks,
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": "F635",
        "status": artifact["status"],
        "lane": artifact["lane"],
        "future_lift_bind_target": artifact["future_lift_bind_target"],
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_JSON)


if __name__ == "__main__":
    main()

