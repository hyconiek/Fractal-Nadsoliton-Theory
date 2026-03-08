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
    / "p125_genuinely_new_strict_core_source_object_clause_probe_for_s_sel_int_candidate_seed_v0.json"
)
OUT_SUMMARY = (
    GENERATED
    / "p125_genuinely_new_strict_core_source_object_clause_probe_for_s_sel_int_candidate_seed_v0_summary.json"
)


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    f34 = load_json(
        "fundamental_action_reconstruction/generated/f34_minimal_admissible_strict_core_source_seed_construction_contract_packet_summary.json"
    )
    f36 = load_json(
        "fundamental_action_reconstruction/generated/f36_first_source_seed_candidate_construction_instance_packet_summary.json"
    )
    p124 = load_json(
        "fundamental_action_reconstruction/generated/p124_first_source_seed_genuinely_new_strict_core_object_clause_probe_summary.json"
    )

    clause_satisfied = False

    checks_spec = [
        {
            "id": "p124_first_clause_fixed",
            "actual": p124["target_state"]["first_clause_target"]["first_clause_name"],
            "expected": "genuinely_new_strict_core_source_object_required",
            "meaning": "P124 already fixes the first clause to be tested",
        },
        {
            "id": "f34_first_clause_required",
            "actual": f34["minimal_source_seed_construction_contract"]["genuinely_new_strict_core_source_object_required"],
            "expected": True,
            "meaning": "F34 requires a genuinely new strict-core source object",
        },
        {
            "id": "f36_candidate_built_from_current_sigma_candidate",
            "actual": f36["candidate_seed_instance"]["construction_route"]["internal_binary_candidate"],
            "expected": "sigma_int_candidate",
            "meaning": "F36 candidate seed still reuses the current sigma_int candidate",
        },
        {
            "id": "f36_candidate_built_from_current_topological_layer",
            "actual": f36["candidate_seed_instance"]["construction_route"]["local_topological_protection_layer"],
            "expected": "QW-2206_local_topological_protection_layer",
            "meaning": "F36 candidate seed still reuses the current QW-2206 local protection layer",
        },
        {
            "id": "first_clause_currently_satisfied",
            "actual": clause_satisfied,
            "expected": False,
            "meaning": "current repo does not yet show that the candidate seed is a genuinely new strict-core source object",
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
        "stage": "P125",
        "lane": "first_source_seed_first_clause_test_only",
        "goal": "test_whether_S_sel_int_candidate_seed_v0_already_satisfies_the_genuinely_new_strict_core_source_object_clause",
        "status": "CURRENT_REPO_DOES_NOT_YET_SHOW_THAT_S_SEL_INT_CANDIDATE_SEED_V0_SATISFIES_THE_GENUINELY_NEW_STRICT_CORE_SOURCE_OBJECT_CLAUSE_AFTER_P125",
        "clause_test_result": {
            "candidate_seed_name": "S_sel_int_candidate_seed_v0",
            "first_clause_name": "genuinely_new_strict_core_source_object_required",
            "currently_satisfied": clause_satisfied,
            "current_reason": "candidate_seed_is_still_packaged_from_current_route_ingredients_rather_than_exported_as_a_genuinely_new_strict_core_source_object",
        },
        "checks": checks,
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": "P125",
        "status": artifact["status"],
        "lane": artifact["lane"],
        "clause_test_result": artifact["clause_test_result"],
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_JSON)


if __name__ == "__main__":
    main()
