#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = (
    GENERATED
    / "n123_current_legacy_to_strict_kernel_package_level_nonbridge_theorem_summary.json"
)


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    p112 = load_json(
        "fundamental_action_reconstruction/generated/p112_current_legacy_to_strict_kernel_package_level_nonbridge_probe_summary.json"
    )

    checks_spec = [
        {
            "id": "p112_package_nonbridge_supported",
            "actual": p112["support_state"]["package_level_nonbridge_supported"],
            "expected": True,
            "meaning": "P112 already supports the package-level nonbridge conclusion",
        },
        {
            "id": "kernel_nonidentification_present",
            "actual": p112["support_state"]["kernel_nonidentification_present"],
            "expected": True,
            "meaning": "the kernel nonidentification ingredient is present",
        },
        {
            "id": "legacy_package_closed_negatively",
            "actual": p112["support_state"]["legacy_package_closed_negatively"],
            "expected": True,
            "meaning": "the full legacy package closure ingredient is present",
        },
        {
            "id": "package_nontransfer_present",
            "actual": p112["support_state"]["package_nontransfer_present"],
            "expected": True,
            "meaning": "the package-level nontransfer ingredient is present",
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
            "step": "N123",
            "status": "N123_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_PACKAGE_LEVEL_NONBRIDGE_STATE",
            "scope": "current_legacy_to_strict_package_level_nonbridge_question_only",
            "checks": checks,
            "blocking_mismatches": mismatches,
            "theorem_result": {"discharged": False},
        }
    else:
        summary = {
            "step": "N123",
            "status": "N123_DISCHARGED_CURRENT_LEGACY_TO_STRICT_KERNEL_PACKAGE_LEVEL_NONBRIDGE_THEOREM_NO_FALSE_PASS",
            "scope": "current_legacy_to_strict_package_level_nonbridge_question_only",
            "checks": checks,
            "theorem_result": {
                "discharged": True,
                "package_level_nonbridge_on_current_repo_state": True,
                "legacy_to_strict_bridge_present": False,
                "full_closure_pass": False,
            },
            "remaining_open_branches": [
                "explicit_strict_core_internal_selector_source_derivation_discharge",
            ],
            "hard_limits": [
                "no_proof_that_no_bridge_can_ever_exist",
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
