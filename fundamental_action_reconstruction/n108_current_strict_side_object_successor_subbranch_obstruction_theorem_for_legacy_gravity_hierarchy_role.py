#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
OUT = (
    ROOT
    / "generated"
    / "n108_current_strict_side_object_successor_subbranch_obstruction_theorem_for_legacy_gravity_hierarchy_role_summary.json"
)


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    f25 = load_json(
        "fundamental_action_reconstruction/generated/f25_legacy_gravity_hierarchy_object_successor_refinement_packet_summary.json"
    )
    p101 = load_json(
        "fundamental_action_reconstruction/generated/p101_strict_side_object_successor_subbranch_probe_for_legacy_gravity_hierarchy_role_summary.json"
    )

    checks_spec = [
        {
            "id": "f25_object_chain_present",
            "actual": f25["candidate_state"]["object_chain_present"],
            "expected": True,
            "meaning": "F25 confirms that the repo exports the gravity_hierarchy_beta20 object chain",
        },
        {
            "id": "textual_object_successor_verdict_absent",
            "actual": p101["textual_object_successor_verdict_present"],
            "expected": False,
            "meaning": "the explicit textual object-successor verdict is absent on the current repo state",
        },
        {
            "id": "object_lineage_upgrade_verdict_absent",
            "actual": p101["object_lineage_upgrade_verdict_present"],
            "expected": False,
            "meaning": "the explicit object-lineage-upgrade verdict is absent on the current repo state",
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
            "step": "N108",
            "status": "N108_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_LEGACY_GRAVITY_HIERARCHY_OBJECT_SUCCESSOR_STATE",
            "scope": "current_legacy_gravity_hierarchy_object_successor_subbranch_question_only",
            "checks": checks,
            "blocking_mismatches": mismatches,
            "theorem_result": {"discharged": False},
        }
    else:
        summary = {
            "step": "N108",
            "status": "N108_DISCHARGED_CURRENT_STRICT_SIDE_OBJECT_SUCCESSOR_SUBBRANCH_OBSTRUCTION_THEOREM_FOR_LEGACY_GRAVITY_HIERARCHY_ROLE_NO_FALSE_PASS",
            "scope": "current_legacy_gravity_hierarchy_object_successor_subbranch_question_only",
            "checks": checks,
            "theorem_result": {
                "discharged": True,
                "textual_object_successor_verdict_present": False,
                "object_lineage_upgrade_verdict_present": False,
                "full_closure_pass": False,
            },
            "remaining_open_branch": [
                "explicit_textual_object_successor_verdict_identifying_gravity_hierarchy_beta20_as_the_strict_side_successor_object_replacing_the_legacy_gravity_hierarchy_role",
                "explicit_object_lineage_upgrade_verdict_elevating_the_existing_gravity_hierarchy_beta20_candidate_chain_into_replacement_semantics_for_the_legacy_gravity_hierarchy_role",
            ],
            "hard_limits": [
                "no_proof_that_either_object_successor_subbranch_is_impossible_forever",
                "the_method_successor_semantics_branch",
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
