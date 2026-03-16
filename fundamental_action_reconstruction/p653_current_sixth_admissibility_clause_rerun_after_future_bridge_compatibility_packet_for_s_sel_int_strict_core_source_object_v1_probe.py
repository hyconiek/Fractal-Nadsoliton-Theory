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
    / "p653_current_sixth_admissibility_clause_rerun_after_future_bridge_compatibility_packet_for_s_sel_int_strict_core_source_object_v1_probe.json"
)
OUT_SUMMARY = (
    GENERATED
    / "p653_current_sixth_admissibility_clause_rerun_after_future_bridge_compatibility_packet_for_s_sel_int_strict_core_source_object_v1_probe_summary.json"
)


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    f34 = load_json(
        "fundamental_action_reconstruction/generated/f34_minimal_admissible_strict_core_source_seed_construction_contract_packet_summary.json"
    )
    n540 = load_json(
        "fundamental_action_reconstruction/generated/n540_current_first_admissibility_clause_discharge_theorem_for_s_sel_int_strict_core_source_object_v1_summary.json"
    )
    n541 = load_json(
        "fundamental_action_reconstruction/generated/n541_current_second_admissibility_clause_discharge_theorem_for_s_sel_int_strict_core_source_object_v1_summary.json"
    )
    n542 = load_json(
        "fundamental_action_reconstruction/generated/n542_current_third_admissibility_clause_discharge_theorem_for_s_sel_int_strict_core_source_object_v1_summary.json"
    )
    n543 = load_json(
        "fundamental_action_reconstruction/generated/n543_current_fourth_admissibility_clause_discharge_theorem_for_s_sel_int_strict_core_source_object_v1_summary.json"
    )
    n544 = load_json(
        "fundamental_action_reconstruction/generated/n544_current_fifth_admissibility_clause_discharge_theorem_for_s_sel_int_strict_core_source_object_v1_summary.json"
    )
    f649 = load_json(
        "fundamental_action_reconstruction/generated/f649_first_exported_s_sel_int_strict_core_source_object_second_clause_typing_packet_summary.json"
    )
    f650 = load_json(
        "fundamental_action_reconstruction/generated/f650_first_exported_s_sel_int_strict_core_source_object_source_seed_only_packet_summary.json"
    )
    f652 = load_json(
        "fundamental_action_reconstruction/generated/f652_first_exported_s_sel_int_strict_core_source_object_selector_acceptance_independence_packet_summary.json"
    )
    f653 = load_json(
        "fundamental_action_reconstruction/generated/f653_first_exported_s_sel_int_strict_core_source_object_future_bridge_compatibility_packet_summary.json"
    )

    exported_object = "S_sel_int_strict_core_source_object_v1"
    sixth_clause = "future_bridge_compatible_required"

    checks_spec = [
        {
            "id": "f34_sixth_clause_required",
            "actual": f34["minimal_source_seed_construction_contract"][sixth_clause],
            "expected": True,
            "meaning": "F34 keeps the future-bridge-compatibility clause active for admissible S_sel_int",
        },
        {"id": "first_clause_discharged", "actual": bool(n540["theorem_result"]["discharged"]), "expected": True},
        {"id": "second_clause_discharged", "actual": bool(n541["theorem_result"]["discharged"]), "expected": True},
        {"id": "third_clause_discharged", "actual": bool(n542["theorem_result"]["discharged"]), "expected": True},
        {"id": "fourth_clause_discharged", "actual": bool(n543["theorem_result"]["discharged"]), "expected": True},
        {"id": "fifth_clause_discharged", "actual": bool(n544["theorem_result"]["discharged"]), "expected": True},
        {
            "id": "later_E_orient_export_meaningful_target_frame",
            "actual": f649["future_orientation_export_target_frame"],
            "expected": ["c1", "s1"],
        },
        {
            "id": "E_orient_not_exported",
            "actual": bool(f649["E_orient_exported"]),
            "expected": False,
        },
        {
            "id": "source_seed_only",
            "actual": bool(f650["source_seed_only"]),
            "expected": True,
        },
        {
            "id": "no_E_orient_smuggling",
            "actual": bool(f650["E_orient_exported"]),
            "expected": False,
        },
        {
            "id": "kernel_split_safe",
            "actual": bool(f653["kernel_split_safe"]),
            "expected": True,
        },
        {
            "id": "no_external_selector_import",
            "actual": bool(f653["no_external_selector_import"]),
            "expected": True,
        },
        {
            "id": "selector_acceptance_independence_axiom_lane_free",
            "actual": bool(f652["uses_axiom_lane_artifact"]),
            "expected": False,
        },
        {
            "id": "selector_acceptance_independence_scope",
            "actual": f652["accepted_scope"],
            "expected": "axiom_augmented_only",
        },
        {
            "id": "strict_core_unchanged_by_selector_acceptance",
            "actual": bool(f652["strict_core_changed"]),
            "expected": False,
        },
        {
            "id": "upstream_of_observer",
            "actual": bool(f653["upstream_of_observer"]),
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
                "meaning": item.get("meaning") or "",
            }
        )
        if not ok:
            mismatches.append(item["id"])

    if mismatches:
        status = "P653_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_SIXTH_CLAUSE_STATE_FOR_S_SEL_INT_AFTER_FUTURE_BRIDGE_COMPATIBILITY_PACKET"
        artifact: dict[str, Any] = {
            "stage": "P653",
            "lane": "sixth_admissibility_clause_rerun_after_future_bridge_compatibility_only",
            "status": status,
            "exported_object": exported_object,
            "sixth_clause": sixth_clause,
            "checks": checks,
            "blocking_mismatches": mismatches,
            "no_false_pass": True,
        }
        summary: dict[str, Any] = {
            "stage": "P653",
            "lane": artifact["lane"],
            "status": status,
            "checks": checks,
            "blocking_mismatches": mismatches,
            "no_false_pass": True,
        }
    else:
        status = "CURRENT_REPO_EXPORTS_ONE_STRICT_CORE_SOURCE_OBJECT_SATISFYING_THE_SIXTH_ADMISSIBILITY_CLAUSE_FOR_S_SEL_INT_AFTER_P653"
        artifact = {
            "stage": "P653",
            "lane": "sixth_admissibility_clause_rerun_after_future_bridge_compatibility_only",
            "status": status,
            "exported_object": exported_object,
            "sixth_clause": sixth_clause,
            "checks": checks,
            "remaining_admissibility_clauses_unresolved": [],
            "no_false_pass": True,
        }
        summary = {
            "stage": "P653",
            "lane": artifact["lane"],
            "status": status,
            "exported_object": exported_object,
            "sixth_clause": sixth_clause,
            "checks": checks,
            "remaining_admissibility_clauses_unresolved": [],
            "no_false_pass": True,
        }

    OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()

