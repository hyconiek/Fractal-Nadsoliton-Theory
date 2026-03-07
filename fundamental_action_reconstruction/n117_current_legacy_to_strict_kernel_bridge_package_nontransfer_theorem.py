#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = (
    GENERATED
    / "n117_current_legacy_to_strict_kernel_bridge_package_nontransfer_theorem_summary.json"
)


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    n50 = load_json(
        "fundamental_action_reconstruction/generated/n50_current_legacy_ontological_kernel_to_strict_gate_kernel_nonidentification_theorem_summary.json"
    )
    n116 = load_json(
        "fundamental_action_reconstruction/generated/n116_current_legacy_physical_role_package_full_negative_closure_theorem_summary.json"
    )
    p107 = load_json(
        "fundamental_action_reconstruction/generated/p107_legacy_to_strict_kernel_bridge_and_role_package_probe_summary.json"
    )

    checks_spec = [
        {
            "id": "n50_kernel_bridge_absent",
            "actual": n50["theorem_result"]["rigorous_bridge_present"],
            "expected": False,
            "meaning": "N50 already keeps rigorous kernel identification absent",
        },
        {
            "id": "n116_legacy_package_closed",
            "actual": n116["theorem_result"]["legacy_physical_role_package_closed_negatively_on_current_repo_state"],
            "expected": True,
            "meaning": "N116 already closes the full legacy physical-role package negatively on the current repo state",
        },
        {
            "id": "p107_package_bridge_negative",
            "actual": p107["kernel_bridge_present"],
            "expected": False,
            "meaning": "P107 confirms that no package-level bridge is currently exported",
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
            "step": "N117",
            "status": "N117_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_LEGACY_TO_STRICT_PACKAGE_STATE",
            "scope": "current_legacy_to_strict_package_question_only",
            "checks": checks,
            "blocking_mismatches": mismatches,
            "theorem_result": {"discharged": False},
        }
    else:
        summary = {
            "step": "N117",
            "status": "N117_DISCHARGED_CURRENT_LEGACY_TO_STRICT_KERNEL_BRIDGE_PACKAGE_NONTRANSFER_THEOREM_NO_FALSE_PASS",
            "scope": "current_legacy_to_strict_package_question_only",
            "checks": checks,
            "theorem_result": {
                "discharged": True,
                "rigorous_legacy_to_strict_kernel_identification_present": False,
                "legacy_physical_role_package_closed_negatively_on_current_repo_state": True,
                "package_level_bridge_present": False,
                "legacy_to_strict_package_nontransfer_on_current_repo_state": True,
                "full_closure_pass": False,
            },
            "remaining_open_branches": [
                "explicit_legacy_to_strict_kernel_bridge_or_nonbridge_theorem_with_package_level_scope",
                "explicit_selector_or_symmetry_breaking_requirement_resolving_the_QW2191_frontier"
            ],
            "hard_limits": [
                "no_proof_that_no_bridge_can_ever_exist",
                "no_proof_that_strict_side_successor_semantics_can_never_exist",
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
