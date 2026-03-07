#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
OUT = (
    ROOT
    / "generated"
    / "n107_current_legacy_gravity_hierarchy_replaced_successor_subbranch_obstruction_theorem_summary.json"
)


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    f24 = load_json(
        "fundamental_action_reconstruction/generated/f24_legacy_gravity_hierarchy_replaced_branch_refinement_packet_summary.json"
    )
    p100 = load_json(
        "fundamental_action_reconstruction/generated/p100_legacy_gravity_hierarchy_replaced_successor_subbranch_probe_summary.json"
    )

    checks_spec = [
        {
            "id": "f24_object_candidate_present",
            "actual": f24["candidate_state"]["strict_object_candidate_present"],
            "expected": True,
            "meaning": "F24 confirms that the strict side exports gravity_hierarchy_beta20 as a candidate object",
        },
        {
            "id": "f24_method_candidate_present",
            "actual": f24["candidate_state"]["strict_method_candidate_present"],
            "expected": True,
            "meaning": "F24 confirms that the strict side exports qw2115_micro_supported_beta_hierarchy_bridge as a candidate method",
        },
        {
            "id": "object_replaced_verdict_absent",
            "actual": p100["object_replaced_verdict_present"],
            "expected": False,
            "meaning": "the explicit object-successor verdict is absent on the current repo state",
        },
        {
            "id": "method_replaced_verdict_absent",
            "actual": p100["method_replaced_verdict_present"],
            "expected": False,
            "meaning": "the explicit method-successor-semantics verdict is absent on the current repo state",
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
            "step": "N107",
            "status": "N107_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_LEGACY_GRAVITY_HIERARCHY_REPLACED_BRANCH_STATE",
            "scope": "current_legacy_gravity_hierarchy_replaced_successor_subbranch_question_only",
            "checks": checks,
            "blocking_mismatches": mismatches,
            "theorem_result": {"discharged": False},
        }
    else:
        summary = {
            "step": "N107",
            "status": "N107_DISCHARGED_CURRENT_LEGACY_GRAVITY_HIERARCHY_REPLACED_SUCCESSOR_SUBBRANCH_OBSTRUCTION_THEOREM_NO_FALSE_PASS",
            "scope": "current_legacy_gravity_hierarchy_replaced_successor_subbranch_question_only",
            "checks": checks,
            "theorem_result": {
                "discharged": True,
                "object_replaced_verdict_present": False,
                "method_replaced_verdict_present": False,
                "full_closure_pass": False,
            },
            "remaining_open_branch": [
                "explicit_object_successor_verdict_identifying_gravity_hierarchy_beta20_as_the_strict_side_successor_object_replacing_the_legacy_gravity_hierarchy_role",
                "explicit_method_successor_semantics_verdict_identifying_qw2115_micro_supported_beta_hierarchy_bridge_as_the_strict_side_successor_semantics_replacing_the_legacy_gravity_hierarchy_role",
            ],
            "hard_limits": [
                "no_proof_that_either_replaced_subbranch_is_impossible_forever",
                "no_selector_closure",
                "no_QW2191_discharge",
                "no_ToE_closure",
            ],
        }

    OUT.write_text(
        json.dumps(summary, indent=2, ensure_ascii=True) + "\n",
        encoding="ascii",
    )
    print(OUT)


if __name__ == "__main__":
    main()
