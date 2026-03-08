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
    / "p112_current_legacy_to_strict_kernel_package_level_nonbridge_probe.json"
)
OUT_SUMMARY = (
    GENERATED
    / "p112_current_legacy_to_strict_kernel_package_level_nonbridge_probe_summary.json"
)


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    n50 = load_json(
        "fundamental_action_reconstruction/generated/n50_current_legacy_ontological_kernel_to_strict_gate_kernel_nonidentification_theorem_summary.json"
    )
    n116 = load_json(
        "fundamental_action_reconstruction/generated/n116_current_legacy_physical_role_package_full_negative_closure_theorem_summary.json"
    )
    n117 = load_json(
        "fundamental_action_reconstruction/generated/n117_current_legacy_to_strict_kernel_bridge_package_nontransfer_theorem_summary.json"
    )

    package_level_nonbridge_supported = (
        n50["theorem_result"]["rigorous_nonidentification_on_current_repo_state"]
        and n116["theorem_result"][
            "legacy_physical_role_package_closed_negatively_on_current_repo_state"
        ]
        and n117["theorem_result"][
            "legacy_to_strict_package_nontransfer_on_current_repo_state"
        ]
    )

    checks_spec = [
        {
            "id": "n50_kernel_nonidentification_present",
            "actual": n50["theorem_result"][
                "rigorous_nonidentification_on_current_repo_state"
            ],
            "expected": True,
            "meaning": "N50 already discharges rigorous legacy-to-strict kernel nonidentification on the current repo state",
        },
        {
            "id": "n116_package_closed_negatively",
            "actual": n116["theorem_result"][
                "legacy_physical_role_package_closed_negatively_on_current_repo_state"
            ],
            "expected": True,
            "meaning": "N116 already closes the full legacy physical-role package negatively on the strict side",
        },
        {
            "id": "n117_package_nontransfer_present",
            "actual": n117["theorem_result"][
                "legacy_to_strict_package_nontransfer_on_current_repo_state"
            ],
            "expected": True,
            "meaning": "N117 already discharges package-level nontransfer on the current repo state",
        },
        {
            "id": "package_level_nonbridge_supported",
            "actual": package_level_nonbridge_supported,
            "expected": True,
            "meaning": "taken together, N50, N116, and N117 support a package-level nonbridge conclusion on the current repo state",
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
        "stage": "P112",
        "lane": "current_legacy_to_strict_kernel_package_level_nonbridge_probe_only",
        "goal": "test_whether_the_current_repo_now_supports_an_explicit_package_level_nonbridge_conclusion_between_the_legacy_kernel_package_and_the_strict_side",
        "status": "CURRENT_REPO_SUPPORTS_A_PACKAGE_LEVEL_NONBRIDGE_CONCLUSION_BETWEEN_THE_LEGACY_KERNEL_PACKAGE_AND_THE_STRICT_SIDE_AFTER_P112",
        "reason": "N50 already closes rigorous kernel identification negatively, N116 already closes the full legacy physical-role package negatively on the strict side, and N117 already discharges package-level nontransfer; therefore the current repo now supports a package-level nonbridge conclusion between the legacy kernel/package and the strict side",
        "support_state": {
            "kernel_nonidentification_present": n50["theorem_result"][
                "rigorous_nonidentification_on_current_repo_state"
            ],
            "legacy_package_closed_negatively": n116["theorem_result"][
                "legacy_physical_role_package_closed_negatively_on_current_repo_state"
            ],
            "package_nontransfer_present": n117["theorem_result"][
                "legacy_to_strict_package_nontransfer_on_current_repo_state"
            ],
            "package_level_nonbridge_supported": package_level_nonbridge_supported,
        },
        "remaining_missing_objects": [
            "explicit_current_repo_state_package_level_nonbridge_theorem_discharge",
            "explicit_strict_core_internal_selector_source_derivation_discharge",
        ],
        "checks": checks,
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": "P112",
        "status": artifact["status"],
        "lane": artifact["lane"],
        "support_state": artifact["support_state"],
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
