#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = (
    GENERATED
    / "p657_current_exported_s_sel_int_strict_core_source_object_selector_output_operator_probe_summary.json"
)


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    n548 = load_json(
        "fundamental_action_reconstruction/generated/n548_current_exported_s_sel_int_strict_core_source_object_selector_reduction_operator_theorem_summary.json"
    )
    f657 = load_json(
        "fundamental_action_reconstruction/generated/f657_first_exported_s_sel_int_strict_core_source_object_selector_output_operator_packet_summary.json"
    )

    props = f657["selector_output_properties"]
    response = f657["source_selector_output_response"]
    checks_spec = [
        {
            "id": "selector_reduction_is_admissible",
            "actual": bool(n548["theorem_result"]["admissible_R_sel"]),
            "expected": True,
        },
        {
            "id": "derived_only_from_selector_reduction",
            "actual": props["derived_only_from_selector_reduction"],
            "expected": True,
        },
        {
            "id": "strict_core_only",
            "actual": props["strict_core_only"],
            "expected": True,
        },
        {
            "id": "selector_output_sector_exported",
            "actual": props["selector_output_sector_exported"],
            "expected": True,
        },
        {
            "id": "positive_plus_output",
            "actual": response["positive_plus_output"],
            "expected": True,
        },
        {
            "id": "vanishing_minus_output",
            "actual": response["vanishing_minus_output"],
            "expected": True,
        },
        {
            "id": "bridge_ready_for_downstream_completion",
            "actual": props["bridge_ready_for_downstream_completion"],
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
            "stage": "P657",
            "lane": "current_exported_s_sel_int_strict_core_source_object_selector_output_operator_only",
            "status": "P657_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_SELECTOR_OUTPUT_STATE",
            "checks": checks,
            "blocking_mismatches": mismatches,
            "no_false_pass": True,
        }
    else:
        summary = {
            "stage": "P657",
            "lane": "current_exported_s_sel_int_strict_core_source_object_selector_output_operator_only",
            "status": "CURRENT_REPO_EXPORTS_ONE_ADMISSIBLE_SELECTOR_OUTPUT_OPERATOR_FROM_R_SEL_S_SEL_INT_SOURCE_OBJECT_V1_AFTER_P657",
            "source_object": "S_sel_int_strict_core_source_object_v1",
            "selector_reduction_operator": "R_sel_s_sel_int_source_object_v1",
            "selector_output_operator": "O_sel_s_sel_int_source_object_v1",
            "checks": checks,
            "actual_O_sel_constructed": True,
            "emergent_observer_constructed": False,
            "strict_core_selector_closure": False,
            "no_false_pass": True,
        }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()

