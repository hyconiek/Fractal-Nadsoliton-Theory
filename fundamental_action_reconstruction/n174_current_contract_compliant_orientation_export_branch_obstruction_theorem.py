#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "n174_current_contract_compliant_orientation_export_branch_obstruction_theorem_summary.json"


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    p158 = load_json(
        "fundamental_action_reconstruction/generated/p158_current_contract_compliant_orientation_export_branch_discharge_probe_summary.json"
    )

    checks_spec = [
        {
            "id": "p158_fixed_attempt_instance",
            "actual": p158["fixed_attempt_instance"],
            "expected": "construct_attempt_v2(S_sel_int_future_additive_upstream_target_v2)",
            "meaning": "P158 fixes the contract-compliant additive construction attempt under discussion",
        },
        {
            "id": "p158_tested_lower_branch",
            "actual": p158["tested_lower_branch"],
            "expected": "future_derivation_of_admissible_E_orient_from_a_future_new_source_object_for_the_fixed_first_contract_compliant_additive_attempt",
            "meaning": "P158 fixes the first remaining orientation-export branch under test",
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
        summary = {
            "step": "N174",
            "status": "N174_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_CONTRACT_COMPLIANT_ORIENTATION_EXPORT_BRANCH_STATE",
            "scope": "current_contract_compliant_orientation_export_branch_discharge_only",
            "checks": checks,
            "blocking_mismatches": mismatches,
            "theorem_result": {"discharged": False},
        }
    else:
        summary = {
            "step": "N174",
            "status": "N174_DISCHARGED_CURRENT_CONTRACT_COMPLIANT_ORIENTATION_EXPORT_BRANCH_OBSTRUCTION_THEOREM_NO_FALSE_PASS",
            "scope": "current_contract_compliant_orientation_export_branch_discharge_only",
            "checks": checks,
            "theorem_result": {
                "discharged": True,
                "current_contract_compliant_orientation_export_branch_obstructed": True,
                "fixed_attempt_instance": "construct_attempt_v2(S_sel_int_future_additive_upstream_target_v2)",
                "tested_lower_branch": "future_derivation_of_admissible_E_orient_from_a_future_new_source_object_for_the_fixed_first_contract_compliant_additive_attempt",
                "full_closure_pass": False,
            },
            "hard_limits": [
                "admissible_E_orient_not_yet_constructed",
                "downstream_chain_not_yet_constructed",
                "constructed_source_object_not_yet_exported",
                "admissible_S_sel_int_not_yet_constructed",
                "no_strict_core_selector_closure",
                "no_QW2191_discharge",
                "no_ToE_closure",
            ],
        }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
