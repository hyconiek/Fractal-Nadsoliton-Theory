#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = (
    GENERATED
    / "p656_current_exported_s_sel_int_strict_core_source_object_selector_reduction_operator_probe_summary.json"
)


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    n546 = load_json(
        "fundamental_action_reconstruction/generated/n546_current_exported_s_sel_int_strict_core_source_object_admissible_orientation_export_theorem_summary.json"
    )
    n547 = load_json(
        "fundamental_action_reconstruction/generated/n547_current_exported_s_sel_int_strict_core_source_object_selector_bridge_operator_theorem_summary.json"
    )
    f656 = load_json(
        "fundamental_action_reconstruction/generated/f656_first_exported_s_sel_int_strict_core_source_object_selector_reduction_operator_packet_summary.json"
    )

    props = f656["selector_reduction_properties"]
    response = f656["source_selector_response"]
    checks_spec = [
        {
            "id": "orientation_input_is_admissible",
            "actual": bool(n546["theorem_result"]["admissible_E_orient"]),
            "expected": True,
        },
        {
            "id": "selector_bridge_is_admissible",
            "actual": bool(n547["theorem_result"]["admissible_B_sel"]),
            "expected": True,
        },
        {
            "id": "derived_only_from_orientation_and_bridge",
            "actual": props["derived_only_from_orientation_and_bridge"],
            "expected": True,
        },
        {
            "id": "strict_core_only",
            "actual": props["strict_core_only"],
            "expected": True,
        },
        {
            "id": "internal_selector_sector_exported",
            "actual": props["internal_selector_sector_exported"],
            "expected": True,
        },
        {
            "id": "positive_plus_channel",
            "actual": response["positive_plus_channel"],
            "expected": True,
        },
        {
            "id": "vanishing_minus_channel",
            "actual": response["vanishing_minus_channel"],
            "expected": True,
        },
        {
            "id": "bridge_ready_for_O_sel",
            "actual": props["bridge_ready_for_O_sel"],
            "expected": True,
        },
        {
            "id": "kernel_split_safe",
            "actual": props["kernel_split_safe"],
            "expected": True,
        },
    ]

    checks = []
    mismatches = []
    for item in checks_spec:
        ok = item["actual"] == item["expected"]
        checks.append(
            {
                "id": item["id"],
                "actual": item["actual"],
                "expected": item["expected"],
                "pass": ok,
            }
        )
        if not ok:
            mismatches.append(item["id"])

    if mismatches:
        summary = {
            "stage": "P656",
            "lane": "current_exported_s_sel_int_strict_core_source_object_selector_reduction_operator_only",
            "status": "P656_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_SELECTOR_REDUCTION_STATE",
            "checks": checks,
            "blocking_mismatches": mismatches,
            "no_false_pass": True,
        }
    else:
        summary = {
            "stage": "P656",
            "lane": "current_exported_s_sel_int_strict_core_source_object_selector_reduction_operator_only",
            "status": "CURRENT_REPO_EXPORTS_ONE_ADMISSIBLE_SELECTOR_REDUCTION_OPERATOR_FROM_B_SEL_S_SEL_INT_SOURCE_OBJECT_V1_AFTER_P656",
            "source_object": "S_sel_int_strict_core_source_object_v1",
            "orientation_input": "E_orient_s_sel_int_source_object_v1",
            "selector_bridge_operator": "B_sel_s_sel_int_source_object_v1",
            "selector_reduction_operator": "R_sel_s_sel_int_source_object_v1",
            "checks": checks,
            "actual_R_sel_constructed": True,
            "actual_O_sel_constructed": False,
            "strict_core_selector_closure": False,
            "no_false_pass": True,
        }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()

