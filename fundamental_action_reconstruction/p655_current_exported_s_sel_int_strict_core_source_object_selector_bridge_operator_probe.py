#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = (
    GENERATED
    / "p655_current_exported_s_sel_int_strict_core_source_object_selector_bridge_operator_probe_summary.json"
)


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    n546 = load_json(
        "fundamental_action_reconstruction/generated/n546_current_exported_s_sel_int_strict_core_source_object_admissible_orientation_export_theorem_summary.json"
    )
    f655 = load_json(
        "fundamental_action_reconstruction/generated/f655_first_exported_s_sel_int_strict_core_source_object_selector_bridge_operator_packet_summary.json"
    )

    props = f655["operator_properties"]
    checks_spec = [
        {
            "id": "orientation_input_is_admissible",
            "actual": bool(n546["theorem_result"]["admissible_E_orient"]),
            "expected": True,
        },
        {
            "id": "derived_only_from_orientation_export",
            "actual": props["derived_only_from_orientation_export"],
            "expected": True,
        },
        {
            "id": "strict_core_only",
            "actual": props["strict_core_only"],
            "expected": True,
        },
        {
            "id": "symmetric",
            "actual": props["symmetric"],
            "expected": True,
        },
        {
            "id": "traceless_on_pair1_plane",
            "actual": props["traceless_on_pair1_plane"],
            "expected": True,
        },
        {
            "id": "exports_selector_decomposition",
            "actual": props["exports_internal_signed_selector_decomposition"],
            "expected": True,
        },
        {
            "id": "selector_bearing_without_external_anchor",
            "actual": props["selector_bearing_without_external_anchor"]
            and not props["uses_imported_psi0"]
            and not props["uses_external_selector_control"]
            and not props["uses_observer_information"],
            "expected": True,
        },
        {
            "id": "positive_source_alignment_witness",
            "actual": f655["source_alignment_witness"]["positive_signed_selector_response"],
            "expected": True,
        },
        {
            "id": "bridge_ready_for_R_sel",
            "actual": props["bridge_ready_for_R_sel"],
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
            "stage": "P655",
            "lane": "current_exported_s_sel_int_strict_core_source_object_selector_bridge_operator_only",
            "status": "P655_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_SELECTOR_BRIDGE_STATE",
            "checks": checks,
            "blocking_mismatches": mismatches,
            "no_false_pass": True,
        }
    else:
        summary = {
            "stage": "P655",
            "lane": "current_exported_s_sel_int_strict_core_source_object_selector_bridge_operator_only",
            "status": "CURRENT_REPO_EXPORTS_ONE_ADMISSIBLE_SELECTOR_BRIDGE_OPERATOR_FROM_E_ORIENT_S_SEL_INT_SOURCE_OBJECT_V1_AFTER_P655",
            "source_object": "S_sel_int_strict_core_source_object_v1",
            "orientation_input": "E_orient_s_sel_int_source_object_v1",
            "selector_bridge_operator": "B_sel_s_sel_int_source_object_v1",
            "checks": checks,
            "actual_B_sel_constructed": True,
            "actual_R_sel_constructed": False,
            "actual_O_sel_constructed": False,
            "strict_core_selector_closure": False,
            "no_false_pass": True,
        }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()

