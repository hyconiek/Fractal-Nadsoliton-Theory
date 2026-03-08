#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = (
    GENERATED
    / "n126_current_repo_exports_no_admissible_strict_core_internal_selector_source_object_theorem_summary.json"
)


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    p115 = load_json(
        "fundamental_action_reconstruction/generated/p115_current_admissible_strict_core_internal_selector_source_object_probe_summary.json"
    )

    checks_spec = [
        {
            "id": "p115_no_admissible_object_present",
            "actual": p115["admission_state"][
                "admissible_strict_core_internal_selector_source_object_present"
            ],
            "expected": False,
            "meaning": "P115 already proves that no admissible strict-core selector source object is currently exported",
        },
        {
            "id": "future_move_must_add_new_source_object",
            "actual": p115["admission_state"][
                "future_positive_move_must_add_new_source_object"
            ],
            "expected": True,
            "meaning": "the next positive move must add a genuinely new source object",
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
            "step": "N126",
            "status": "N126_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_ADMISSIBLE_SOURCE_OBJECT_STATE",
            "scope": "current_repo_admissible_strict_core_internal_selector_source_object_question_only",
            "checks": checks,
            "blocking_mismatches": mismatches,
            "theorem_result": {"discharged": False},
        }
    else:
        summary = {
            "step": "N126",
            "status": "N126_DISCHARGED_CURRENT_REPO_EXPORTS_NO_ADMISSIBLE_STRICT_CORE_INTERNAL_SELECTOR_SOURCE_OBJECT_THEOREM_NO_FALSE_PASS",
            "scope": "current_repo_admissible_strict_core_internal_selector_source_object_question_only",
            "checks": checks,
            "theorem_result": {
                "discharged": True,
                "admissible_strict_core_internal_selector_source_object_present": False,
                "future_positive_move_must_add_new_source_object": True,
                "full_closure_pass": False,
            },
            "remaining_open_branches": [
                "future_new_strict_core_internal_selector_source_object",
            ],
            "hard_limits": [
                "no_proof_that_no_future_admissible_source_object_can_exist",
                "no_strict_core_selector_closure",
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
