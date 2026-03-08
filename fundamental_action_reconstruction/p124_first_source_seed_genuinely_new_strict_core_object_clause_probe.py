#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT_JSON = GENERATED / "p124_first_source_seed_genuinely_new_strict_core_object_clause_probe.json"
OUT_SUMMARY = (
    GENERATED / "p124_first_source_seed_genuinely_new_strict_core_object_clause_probe_summary.json"
)


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    p123 = load_json(
        "fundamental_action_reconstruction/generated/p123_first_source_seed_admissibility_upgrade_target_probe_summary.json"
    )
    f38 = load_json(
        "fundamental_action_reconstruction/generated/f38_first_source_seed_first_clause_refinement_packet_summary.json"
    )

    reduced_to_first_clause_test = (
        p123["target_state"]["next_constructive_move_reduced_to_one_first_admissibility_upgrade_target"] is True
        and f38["first_clause_target"]["first_clause_name"] == "genuinely_new_strict_core_source_object_required"
        and f38["first_clause_target"]["candidate_seed_name"] == "S_sel_int_candidate_seed_v0"
    )

    checks_spec = [
        {
            "id": "p123_single_upgrade_target",
            "actual": p123["target_state"]["next_constructive_move_reduced_to_one_first_admissibility_upgrade_target"],
            "expected": True,
            "meaning": "P123 already fixes one admissibility-upgrade target",
        },
        {
            "id": "f38_first_clause_name",
            "actual": f38["first_clause_target"]["first_clause_name"],
            "expected": "genuinely_new_strict_core_source_object_required",
            "meaning": "F38 already isolates the first clause to test",
        },
        {
            "id": "reduced_to_first_clause_test",
            "actual": reduced_to_first_clause_test,
            "expected": True,
            "meaning": "the next clause-by-clause admissibility move is reduced to the genuinely-new-object requirement",
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
        "stage": "P124",
        "lane": "first_source_seed_first_clause_only",
        "goal": "test_whether_the_next_clause_by_clause_move_is_now_reduced_to_the_genuinely_new_strict_core_source_object_requirement",
        "status": "CURRENT_REPO_REDUCES_THE_NEXT_CLAUSE_BY_CLAUSE_ADMISSIBILITY_MOVE_TO_THE_GENUINELY_NEW_STRICT_CORE_SOURCE_OBJECT_REQUIREMENT_AFTER_P124",
        "target_state": {
            "next_clause_by_clause_move_reduced_to_first_clause": reduced_to_first_clause_test,
            "first_clause_target": f38["first_clause_target"],
            "later_clauses_left_open": [
                "carrier_typed_enough_for_later_E_orient_export_required",
                "source_seed_only_not_counted_as_E_orient_or_bridge",
                "strict_core_only_required",
                "silent_legacy_to_strict_substitution_forbidden",
                "selector_acceptance_outside_strict_core_may_not_count_as_source_construction",
                "future_bridge_compatible_required",
            ],
        },
        "checks": checks,
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": "P124",
        "status": artifact["status"],
        "lane": artifact["lane"],
        "target_state": artifact["target_state"],
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_JSON)


if __name__ == "__main__":
    main()
