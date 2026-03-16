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
    / "p652_current_fifth_admissibility_clause_rerun_after_selector_acceptance_independence_packet_for_s_sel_int_strict_core_source_object_v1_probe.json"
)
OUT_SUMMARY = (
    GENERATED
    / "p652_current_fifth_admissibility_clause_rerun_after_selector_acceptance_independence_packet_for_s_sel_int_strict_core_source_object_v1_probe_summary.json"
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
    f652 = load_json(
        "fundamental_action_reconstruction/generated/f652_first_exported_s_sel_int_strict_core_source_object_selector_acceptance_independence_packet_summary.json"
    )

    exported_object = "S_sel_int_strict_core_source_object_v1"
    fifth_clause = "selector_acceptance_outside_strict_core_may_not_count_as_source_construction"

    checks_spec = [
        {
            "id": "f34_fifth_clause_required",
            "actual": f34["minimal_source_seed_construction_contract"][fifth_clause],
            "expected": True,
            "meaning": "F34 keeps the selector-acceptance-independence clause active for admissible S_sel_int",
        },
        {
            "id": "first_clause_discharged",
            "actual": bool(n540["theorem_result"]["discharged"]),
            "expected": True,
        },
        {
            "id": "second_clause_discharged",
            "actual": bool(n541["theorem_result"]["discharged"]),
            "expected": True,
        },
        {
            "id": "third_clause_discharged",
            "actual": bool(n542["theorem_result"]["discharged"]),
            "expected": True,
        },
        {
            "id": "fourth_clause_discharged",
            "actual": bool(n543["theorem_result"]["discharged"]),
            "expected": True,
        },
        {
            "id": "uses_axiom_lane_artifact",
            "actual": bool(f652["uses_axiom_lane_artifact"]),
            "expected": False,
        },
        {
            "id": "strict_core_only",
            "actual": bool(f652["strict_core_only"]),
            "expected": True,
        },
        {
            "id": "selector_requirement_accepted_at_theory_level",
            "actual": bool(f652["selector_requirement_accepted_at_theory_level"]),
            "expected": True,
        },
        {
            "id": "accepted_scope",
            "actual": f652["accepted_scope"],
            "expected": "axiom_augmented_only",
        },
        {
            "id": "strict_core_changed",
            "actual": bool(f652["strict_core_changed"]),
            "expected": False,
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
        status = "P652_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_FIFTH_CLAUSE_STATE_FOR_S_SEL_INT_AFTER_SELECTOR_ACCEPTANCE_INDEPENDENCE_PACKET"
        artifact: dict[str, Any] = {
            "stage": "P652",
            "lane": "fifth_admissibility_clause_rerun_after_selector_acceptance_independence_only",
            "status": status,
            "exported_object": exported_object,
            "fifth_clause": fifth_clause,
            "checks": checks,
            "blocking_mismatches": mismatches,
            "no_false_pass": True,
        }
        summary: dict[str, Any] = {
            "stage": "P652",
            "lane": artifact["lane"],
            "status": status,
            "checks": checks,
            "blocking_mismatches": mismatches,
            "no_false_pass": True,
        }
    else:
        status = "CURRENT_REPO_EXPORTS_ONE_STRICT_CORE_SOURCE_OBJECT_SATISFYING_THE_FIFTH_ADMISSIBILITY_CLAUSE_FOR_S_SEL_INT_AFTER_P652"
        artifact = {
            "stage": "P652",
            "lane": "fifth_admissibility_clause_rerun_after_selector_acceptance_independence_only",
            "status": status,
            "exported_object": exported_object,
            "fifth_clause": fifth_clause,
            "checks": checks,
            "remaining_admissibility_clauses_unresolved": [
                "future_bridge_compatible_required",
            ],
            "no_false_pass": True,
        }
        summary = {
            "stage": "P652",
            "lane": artifact["lane"],
            "status": status,
            "exported_object": exported_object,
            "fifth_clause": fifth_clause,
            "checks": checks,
            "remaining_admissibility_clauses_unresolved": artifact["remaining_admissibility_clauses_unresolved"],
            "no_false_pass": True,
        }

    OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()

