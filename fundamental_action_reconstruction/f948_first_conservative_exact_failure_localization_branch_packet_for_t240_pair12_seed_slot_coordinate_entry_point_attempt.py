#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-21"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT_JSON = (
    GENERATED
    / "f948_first_conservative_exact_failure_localization_branch_packet_for_t240_pair12_seed_slot_coordinate_entry_point_attempt.json"
)
OUT_SUMMARY = (
    GENERATED
    / "f948_first_conservative_exact_failure_localization_branch_packet_for_t240_pair12_seed_slot_coordinate_entry_point_attempt_summary.json"
)

P949_JSON = (
    GENERATED
    / "p949_current_strict_t240_pair12_seed_slot_coordinate_entry_point_subsubsubsubsubinterface_actual_realization_attempt_probe.json"
)
P950_SUMMARY = (
    GENERATED
    / "p950_current_strict_t241_pair12_seed_slot_coordinate_entry_point_subsubsubsubsubinterface_actual_realization_attempt_verdict_or_exact_failure_localization_nonexport_audit_probe_summary.json"
)

FAILURE_BRANCH = (
    "future_exact_failure_localization_below_the_chart_label_retaining_pair12_typed_"
    "seed_slot_coordinate_entry_point_subsubsubsubsubinterface"
)
SUCCESS_BRANCH_PREFIX = "future_success_or_failure_verdict_for_"
ATTEMPT_NAME = (
    "W_strict_t173_pair12_seed_slot_coordinate_entry_point_subsubsubsubsubinterface_"
    "actual_realization_attempt_v1"
)


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prerequisites = [P949_JSON, P950_SUMMARY]
    missing = [str(path.relative_to(REPO)) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "F948",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    p949 = load_json(P949_JSON)
    p950 = load_json(P950_SUMMARY)
    attempt = p949.get("first_actual_t238_subsubsubsubsubinterface_realization_attempt") or {}
    later_open_branches = attempt.get("later_open_branches") or []

    checks_spec = [
        {
            "id": "p950_nonexport_boundary_passed",
            "actual": p950.get("status"),
            "expected": "PASS_STRICT_T241_PAIR12_SEED_SLOT_COORDINATE_ENTRY_POINT_SUBSUBSUBSUBSUBINTERFACE_ACTUAL_REALIZATION_ATTEMPT_VERDICT_OR_EXACT_FAILURE_LOCALIZATION_NONEXPORT_AUDITED",
            "meaning": "P950 already freezes that the exact T240 attempt still lacks both real success verdict and exact lower failure-localization export.",
        },
        {
            "id": "t240_attempt_exported",
            "actual": attempt.get("attempt_name"),
            "expected": ATTEMPT_NAME,
            "meaning": "The fixed exact T240 attempt name is exported in P949.",
        },
        {
            "id": "declared_failure_branch_present",
            "actual": FAILURE_BRANCH in later_open_branches,
            "expected": True,
            "meaning": "P949 already declares one exact future failure-localization branch on the same fixed T240 route.",
        },
        {
            "id": "success_verdict_still_absent",
            "actual": p950.get("current_repo_still_lacks_success_verdict_for_t240_exact_attempt"),
            "expected": True,
            "meaning": "The repo still lacks a real success verdict for the exact T240 attempt.",
        },
        {
            "id": "failure_localization_still_absent",
            "actual": p950.get("current_repo_still_lacks_exact_failure_localization_below_t240_exact_attempt"),
            "expected": True,
            "meaning": "The repo still lacks an actual exact failure-localization export below the exact T240 attempt.",
        },
        {
            "id": "next_move_already_reduced_to_failure_localization_freeze",
            "actual": p950.get("next_honest_move_is_freeze_exact_failure_localization_below_t240_exact_attempt"),
            "expected": True,
            "meaning": "P950 already reduces the next honest move to freezing exact failure-localization below the same T240 attempt.",
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

    blocking = [check["id"] for check in checks if not check["pass"]]
    status = (
        "F948_EXECUTED_FIRST_CONSERVATIVE_EXACT_FAILURE_LOCALIZATION_BRANCH_PACKET_FOR_T240_PAIR12_SEED_SLOT_COORDINATE_ENTRY_POINT_ATTEMPT_NO_FALSE_PASS"
        if not blocking
        else "F948_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_T240_FAILURE_BRANCH_STATE"
    )

    artifact = {
        "stage": "F948",
        "status": status,
        "as_of": AS_OF,
        "lane": "conservative_exact_failure_localization_branch_first_only",
        "goal": "freeze_the_exact_failure_localization_branch_as_the_first_branch_to_attack_under_no_false_pass_ordering_for_the_fixed_t240_attempt",
        "fixed_attempt_name": ATTEMPT_NAME,
        "first_branch_to_attack": FAILURE_BRANCH,
        "sibling_declared_branches": later_open_branches,
        "branch_ordering_basis": [
            "fixed_exact_t240_attempt_already_frozen_by_T240_P949_N782",
            "P950_N783_already_freeze_absence_of_real_success_verdict",
            "P950_N783_already_freeze_absence_of_exact_failure_localization_export",
            "failure_localization_first_avoids_inventing_a_lower_object_by_fiat",
            "failure_localization_first_is_methodologically_more_conservative_than_success_first",
        ],
        "checks": checks,
        "blocking_checks": blocking,
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": "F948",
        "status": status,
        "lane": artifact["lane"],
        "fixed_attempt_name": artifact["fixed_attempt_name"],
        "first_branch_to_attack": artifact["first_branch_to_attack"],
        "branch_ordering_basis": artifact["branch_ordering_basis"],
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
