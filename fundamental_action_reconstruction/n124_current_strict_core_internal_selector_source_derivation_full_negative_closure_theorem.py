#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = (
    GENERATED
    / "n124_current_strict_core_internal_selector_source_derivation_full_negative_closure_theorem_summary.json"
)


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    p113 = load_json(
        "fundamental_action_reconstruction/generated/p113_current_strict_core_internal_selector_source_derivation_discharge_probe_summary.json"
    )

    checks_spec = [
        {
            "id": "p113_package_non_discharge_supported",
            "actual": p113["support_state"][
                "package_level_non_discharge_supported"
            ],
            "expected": True,
            "meaning": "P113 already supports package-level non-discharge for the strict-core internal selector source frontier",
        },
        {
            "id": "generic_hidden_source_branch_closed",
            "actual": p113["support_state"][
                "generic_hidden_source_branch_closed_negatively"
            ],
            "expected": True,
            "meaning": "the generic hidden-source branch is closed negatively",
        },
        {
            "id": "psi0_branch_closed",
            "actual": p113["support_state"]["psi0_branch_closed_negatively"],
            "expected": True,
            "meaning": "the current psi0 branch is closed negatively",
        },
        {
            "id": "fr_branch_closed",
            "actual": p113["support_state"]["fr_branch_closed_negatively"],
            "expected": True,
            "meaning": "the current FR branch is closed negatively",
        },
        {
            "id": "sigma_int_branch_closed",
            "actual": p113["support_state"][
                "sigma_int_bridge_branch_closed_negatively"
            ],
            "expected": True,
            "meaning": "the current sigma-int bridge branch is closed negatively",
        },
        {
            "id": "no_downstream_a1_pair1_reachability",
            "actual": p113["support_state"][
                "strict_core_downstream_A1_pair1_reachability_present"
            ],
            "expected": False,
            "meaning": "no strict-core downstream A1(pair1) reachability is present",
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
            "step": "N124",
            "status": "N124_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_STRICT_CORE_INTERNAL_SELECTOR_SOURCE_STATE",
            "scope": "current_strict_core_internal_selector_source_frontier_only",
            "checks": checks,
            "blocking_mismatches": mismatches,
            "theorem_result": {"discharged": False},
        }
    else:
        summary = {
            "step": "N124",
            "status": "N124_DISCHARGED_CURRENT_STRICT_CORE_INTERNAL_SELECTOR_SOURCE_DERIVATION_FULL_NEGATIVE_CLOSURE_THEOREM_NO_FALSE_PASS",
            "scope": "current_strict_core_internal_selector_source_frontier_only",
            "checks": checks,
            "theorem_result": {
                "discharged": True,
                "strict_core_internal_selector_source_derivation_discharge_present": False,
                "strict_core_internal_selector_source_frontier_closed_negatively_on_current_repo_state": True,
                "full_closure_pass": False,
            },
            "remaining_open_branches": [],
            "hard_limits": [
                "no_proof_that_no_future_strict_core_selector_source_can_exist",
                "no_strict_core_selector_closure",
                "no_QW2191_discharge",
                "no_theory_level_acceptance_of_selector_requirement",
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
