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
    / "f51_first_future_genuinely_additive_source_object_construction_attempt_packet.json"
)
OUT_SUMMARY = (
    GENERATED
    / "f51_first_future_genuinely_additive_source_object_construction_attempt_packet_summary.json"
)


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    f50 = load_json(
        "fundamental_action_reconstruction/generated/f50_genuinely_additive_source_object_construction_contract_packet_summary.json"
    )
    n150 = load_json(
        "fundamental_action_reconstruction/generated/n150_only_remaining_positive_move_class_reduced_to_one_additive_construction_contract_theorem_summary.json"
    )

    checks_spec = [
        {
            "id": "n150_only_remaining_positive_move_class",
            "actual": n150["theorem_result"]["only_remaining_positive_move_class"],
            "expected": "future_genuinely_additive_new_strict_core_source_object_construction",
            "meaning": "N150 already fixes the only remaining positive move class",
        },
        {
            "id": "f50_new_exported_object_identity",
            "actual": f50["minimal_additive_construction_contract"]["new_exported_object_identity"],
            "expected": True,
            "meaning": "F50 requires a future new exported object identity",
        },
        {
            "id": "f50_not_just_packaged_reuse",
            "actual": f50["minimal_additive_construction_contract"]["not_just_packaged_reuse"],
            "expected": True,
            "meaning": "F50 excludes packaged reuse as a positive move",
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
        "stage": "F51",
        "lane": "future_additive_attempt_target_only",
        "goal": "freeze_one_explicit_future_attempt_target_for_the_only_remaining_positive_move_class",
        "status": "F51_EXECUTED_FIRST_FUTURE_GENUINELY_ADDITIVE_SOURCE_OBJECT_CONSTRUCTION_ATTEMPT_PACKET_NO_FALSE_PASS",
        "future_additive_attempt_target": {
            "target_id": "S_sel_int_additive_attempt_target_v1",
            "new_exported_object_identity": True,
            "not_identified_with_any_current_export": True,
            "not_just_packaged_reuse": True,
            "strict_core_only": True,
            "no_external_selector_import": True,
            "kernel_split_safe": True,
            "future_admissibility_target": True,
        },
        "checks": checks,
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": "F51",
        "status": artifact["status"],
        "lane": artifact["lane"],
        "future_additive_attempt_target": artifact["future_additive_attempt_target"],
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
