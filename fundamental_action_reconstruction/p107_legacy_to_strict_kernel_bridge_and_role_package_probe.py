#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT_JSON = (
    GENERATED
    / "p107_legacy_to_strict_kernel_bridge_and_role_package_probe.json"
)
OUT_SUMMARY = (
    GENERATED
    / "p107_legacy_to_strict_kernel_bridge_and_role_package_probe_summary.json"
)


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    n50 = load_json(
        "fundamental_action_reconstruction/generated/n50_current_legacy_ontological_kernel_to_strict_gate_kernel_nonidentification_theorem_summary.json"
    )
    n116 = load_json(
        "fundamental_action_reconstruction/generated/n116_current_legacy_physical_role_package_full_negative_closure_theorem_summary.json"
    )

    checks_spec = [
        {
            "id": "n50_rigorous_kernel_bridge_absent",
            "actual": n50["theorem_result"]["rigorous_bridge_present"],
            "expected": False,
            "meaning": "N50 already keeps rigorous legacy-to-strict kernel identification absent",
        },
        {
            "id": "n116_legacy_package_closed_negatively",
            "actual": n116["theorem_result"]["legacy_physical_role_package_closed_negatively_on_current_repo_state"],
            "expected": True,
            "meaning": "N116 already closes the full legacy physical-role package negatively on the current repo state",
        },
    ]

    checks = []
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
        "stage": "P107",
        "lane": "legacy_to_strict_kernel_bridge_and_role_package_probe_current_repo_state_only",
        "goal": "test_whether_the_current_repo_exports_a_package_level_bridge_identifying_the_strict_side_as_the_rigorous_carrier_of_the_full_legacy_physical_role_package",
        "status": "CURRENT_REPO_EXPORTS_NEITHER_A_RIGOROUS_LEGACY_TO_STRICT_KERNEL_IDENTIFICATION_NOR_A_FULL_LEGACY_PHYSICAL_ROLE_PACKAGE_TRANSFER_TO_THE_STRICT_SIDE_AFTER_P107",
        "reason": "N50 already closes rigorous kernel identification negatively and N116 already closes the full legacy physical-role package negatively on the current repo state; therefore the current repo exports no package-level bridge from the legacy kernel/package to the strict side",
        "kernel_bridge_present": n50["theorem_result"]["rigorous_bridge_present"],
        "legacy_role_package_closed_negatively_on_current_repo_state": n116["theorem_result"]["legacy_physical_role_package_closed_negatively_on_current_repo_state"],
        "remaining_missing_objects": [
            "explicit_legacy_to_strict_kernel_bridge_or_nonbridge_theorem_with_package_level_scope",
            "explicit_selector_or_symmetry_breaking_requirement_resolving_the_QW2191_frontier"
        ],
        "checks": checks,
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": "P107",
        "status": artifact["status"],
        "lane": artifact["lane"],
        "kernel_bridge_present": artifact["kernel_bridge_present"],
        "legacy_role_package_closed_negatively_on_current_repo_state": artifact["legacy_role_package_closed_negatively_on_current_repo_state"],
        "remaining_missing_objects": artifact["remaining_missing_objects"],
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    OUT_JSON.write_text(
        json.dumps(artifact, indent=2, ensure_ascii=True) + "\n",
        encoding="ascii",
    )
    OUT_SUMMARY.write_text(
        json.dumps(summary, indent=2, ensure_ascii=True) + "\n",
        encoding="ascii",
    )
    print(OUT_JSON)


if __name__ == "__main__":
    main()
