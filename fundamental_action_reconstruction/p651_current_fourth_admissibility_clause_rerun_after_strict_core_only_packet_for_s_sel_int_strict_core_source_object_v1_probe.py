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
    / "p651_current_fourth_admissibility_clause_rerun_after_strict_core_only_packet_for_s_sel_int_strict_core_source_object_v1_probe.json"
)
OUT_SUMMARY = (
    GENERATED
    / "p651_current_fourth_admissibility_clause_rerun_after_strict_core_only_packet_for_s_sel_int_strict_core_source_object_v1_probe_summary.json"
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
    f651 = load_json(
        "fundamental_action_reconstruction/generated/f651_first_exported_s_sel_int_strict_core_source_object_strict_core_only_packet_summary.json"
    )

    exported_object = "S_sel_int_strict_core_source_object_v1"
    fourth_clause = "strict_core_only_required"

    checks_spec = [
        {
            "id": "f34_fourth_clause_required",
            "actual": f34["minimal_source_seed_construction_contract"][fourth_clause],
            "expected": True,
            "meaning": "F34 keeps the strict-core-only clause active for admissible S_sel_int",
        },
        {
            "id": "first_clause_discharged",
            "actual": bool(n540["theorem_result"]["discharged"]),
            "expected": True,
            "meaning": "N540 discharges the first admissibility clause for S_sel_int on current repo state",
        },
        {
            "id": "second_clause_discharged",
            "actual": bool(n541["theorem_result"]["discharged"]),
            "expected": True,
            "meaning": "N541 discharges the second admissibility clause for S_sel_int on current repo state",
        },
        {
            "id": "third_clause_discharged",
            "actual": bool(n542["theorem_result"]["discharged"]),
            "expected": True,
            "meaning": "N542 discharges the third admissibility clause for S_sel_int on current repo state",
        },
        {
            "id": "strict_core_only",
            "actual": bool(f651["strict_core_only"]),
            "expected": True,
            "meaning": "F651 explicitly freezes strict_core_only = true for the exported source object",
        },
        {
            "id": "uses_axiom_lane_artifact",
            "actual": bool(f651["uses_axiom_lane_artifact"]),
            "expected": False,
            "meaning": "F651 explicitly freezes uses_axiom_lane_artifact = false",
        },
        {
            "id": "upstream_of_observer",
            "actual": bool(f651["upstream_of_observer"]),
            "expected": True,
            "meaning": "F651 explicitly freezes upstream_of_observer = true",
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
                "meaning": item["meaning"],
            }
        )
        if not ok:
            mismatches.append(item["id"])

    if mismatches:
        status = "P651_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_FOURTH_CLAUSE_STATE_FOR_S_SEL_INT_AFTER_STRICT_CORE_ONLY_PACKET"
        artifact: dict[str, Any] = {
            "stage": "P651",
            "lane": "fourth_admissibility_clause_rerun_after_strict_core_only_only",
            "status": status,
            "exported_object": exported_object,
            "fourth_clause": fourth_clause,
            "checks": checks,
            "blocking_mismatches": mismatches,
            "no_false_pass": True,
        }
        summary: dict[str, Any] = {
            "stage": "P651",
            "lane": artifact["lane"],
            "status": status,
            "checks": checks,
            "blocking_mismatches": mismatches,
            "no_false_pass": True,
        }
    else:
        status = "CURRENT_REPO_EXPORTS_ONE_STRICT_CORE_SOURCE_OBJECT_SATISFYING_THE_FOURTH_ADMISSIBILITY_CLAUSE_FOR_S_SEL_INT_AFTER_P651"
        artifact = {
            "stage": "P651",
            "lane": "fourth_admissibility_clause_rerun_after_strict_core_only_only",
            "status": status,
            "exported_object": exported_object,
            "fourth_clause": fourth_clause,
            "checks": checks,
            "remaining_admissibility_clauses_unresolved": [
                "selector_acceptance_outside_strict_core_may_not_count_as_source_construction",
                "future_bridge_compatible_required",
            ],
            "no_false_pass": True,
        }
        summary = {
            "stage": "P651",
            "lane": artifact["lane"],
            "status": status,
            "exported_object": exported_object,
            "fourth_clause": fourth_clause,
            "checks": checks,
            "remaining_admissibility_clauses_unresolved": artifact["remaining_admissibility_clauses_unresolved"],
            "no_false_pass": True,
        }

    OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()

