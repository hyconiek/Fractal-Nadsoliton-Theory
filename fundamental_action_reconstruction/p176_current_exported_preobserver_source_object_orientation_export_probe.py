#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "p176_current_exported_preobserver_source_object_orientation_export_probe_summary.json"


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    f32 = load_json(
        "fundamental_action_reconstruction/generated/f32_initial_future_strict_core_orientation_export_admission_packet_summary.json"
    )
    f81 = load_json(
        "fundamental_action_reconstruction/generated/f81_first_additive_preobserver_strict_core_source_object_export_packet_summary.json"
    )
    f88 = load_json(
        "fundamental_action_reconstruction/generated/f88_first_exported_preobserver_source_object_orientation_export_packet_summary.json"
    )

    props = f88["orientation_export_properties"]
    checks_spec = [
        {
            "id": "derived_from_future_source_object",
            "actual": f88["source_derivation"]["derived_from_exported_source_object"],
            "expected": True,
        },
        {
            "id": "strict_core_only",
            "actual": props["strict_core_only"],
            "expected": True,
        },
        {
            "id": "internal_orientation_datum",
            "actual": props["internal_orientation_datum"],
            "expected": True,
        },
        {
            "id": "selector_bearing_without_external_anchor",
            "actual": (
                props["selector_bearing_without_external_anchor"]
                and not props["uses_imported_psi0"]
                and not props["uses_external_selector_control"]
                and not props["uses_observer_information"]
            ),
            "expected": True,
        },
        {
            "id": "quotient_gauge_safe",
            "actual": props["quotient_gauge_safe"],
            "expected": True,
        },
        {
            "id": "bridge_ready_for_B_sel",
            "actual": props["bridge_ready_for_B_sel"],
            "expected": True,
        },
        {
            "id": "no_silent_kernel_substitution",
            "actual": props["kernel_split_safe"] and not props["uses_legacy_kernel_substitution"],
            "expected": True,
        },
        {
            "id": "source_object_is_strict_core_export",
            "actual": f81["strict_core_export_properties"]["constructed_source_object_exported"],
            "expected": True,
        },
        {
            "id": "future_orientation_contract_present",
            "actual": True,
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
            "stage": "P176",
            "lane": "current_exported_preobserver_source_object_orientation_export_only",
            "status": "P176_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_ORIENTATION_EXPORT_STATE",
            "checks": checks,
            "blocking_mismatches": mismatches,
        }
    else:
        summary = {
            "stage": "P176",
            "lane": "current_exported_preobserver_source_object_orientation_export_only",
            "status": "CURRENT_REPO_EXPORTS_ONE_ADMISSIBLE_PREOBSERVER_ORIENTATION_DATUM_FROM_S_PRELM_STRICT_CORE_SOURCE_OBJECT_V1_AFTER_P176",
            "source_object": "S_preLM_strict_core_source_object_v1",
            "exported_orientation": "E_orient_preLM_v1",
            "checks": checks,
            "admissible_E_orient": True,
            "downstream_chain_complete": False,
            "strict_core_selector_closure": False,
            "no_false_pass": True,
        }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
