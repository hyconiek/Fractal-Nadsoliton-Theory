#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = (
    GENERATED
    / "p654_current_exported_s_sel_int_strict_core_source_object_orientation_export_probe_summary.json"
)


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    f32 = load_json(
        "fundamental_action_reconstruction/generated/f32_initial_future_strict_core_orientation_export_admission_packet_summary.json"
    )
    n545 = load_json(
        "fundamental_action_reconstruction/generated/n545_current_sixth_admissibility_clause_discharge_theorem_for_s_sel_int_strict_core_source_object_v1_summary.json"
    )
    f654 = load_json(
        "fundamental_action_reconstruction/generated/f654_first_exported_s_sel_int_strict_core_source_object_orientation_export_packet_summary.json"
    )

    props = f654["orientation_export_properties"]
    checks_spec = [
        {
            "id": "derived_from_future_source_object",
            "actual": f654["source_derivation"]["derived_from_exported_source_object"],
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
            "id": "source_object_full_admissibility_pass",
            "actual": bool(n545["theorem_result"]["full_admissibility_pass"]),
            "expected": True,
        },
        {
            "id": "future_orientation_contract_present",
            "actual": bool(f32.get("orientation_export_contract")),
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
            "stage": "P654",
            "lane": "current_exported_s_sel_int_strict_core_source_object_orientation_export_only",
            "status": "P654_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_ORIENTATION_EXPORT_STATE",
            "checks": checks,
            "blocking_mismatches": mismatches,
            "no_false_pass": True,
        }
    else:
        summary = {
            "stage": "P654",
            "lane": "current_exported_s_sel_int_strict_core_source_object_orientation_export_only",
            "status": "CURRENT_REPO_EXPORTS_ONE_ADMISSIBLE_ORIENTATION_DATUM_FROM_S_SEL_INT_STRICT_CORE_SOURCE_OBJECT_V1_AFTER_P654",
            "source_object": "S_sel_int_strict_core_source_object_v1",
            "exported_orientation": f654["exported_orientation"]["object"],
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
