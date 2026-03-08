#!/usr/bin/env python3

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "p242_current_legacy_to_strict_kernel_nonbridge_strengthening_target_probe_summary.json"


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    f152 = load_json(
        GENERATED / "f152_first_legacy_to_strict_kernel_nonbridge_strengthening_target_packet_summary.json"
    )

    checks_spec = [
        {
            "id": "nonbridge_target_exported",
            "actual": f152["nonbridge_strengthening_target"] == "NB_legacy_strict_strengthening_target_v1",
            "expected": True,
        },
        {
            "id": "future_only_target_exported",
            "actual": f152["future_only_target_exported"],
            "expected": True,
        },
        {
            "id": "package_level_nonbridge_base_present",
            "actual": f152["package_level_nonbridge_base_present"],
            "expected": True,
        },
        {
            "id": "positive_bridge_branch_still_open",
            "actual": f152["positive_bridge_branch_still_open"],
            "expected": True,
        },
        {
            "id": "permanent_no_bridge_NOT_claimed",
            "actual": f152["permanent_no_bridge_claimed"],
            "expected": False,
        },
        {
            "id": "kernel_split_safe",
            "actual": f152["kernel_split_safe"],
            "expected": True,
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
            }
        )
        if not ok:
            mismatches.append(item["id"])

    if mismatches:
        status = "P242_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_NONBRIDGE_TARGET_STATE"
    else:
        status = "CURRENT_REPO_EXPORTS_ONE_FUTURE_ONLY_LEGACY_TO_STRICT_KERNEL_NONBRIDGE_STRENGTHENING_TARGET_BELOW_ACTUAL_NONBRIDGE_DISCHARGE_AFTER_P242"

    summary = {
        "stage": "P242",
        "lane": "kernel_nonbridge_strengthening",
        "status": status,
        "checks": checks,
        "blocking_mismatches": mismatches,
        "nonbridge_target_exported": f152["nonbridge_strengthening_target"] == "NB_legacy_strict_strengthening_target_v1",
        "future_only_target_exported": f152["future_only_target_exported"],
        "package_level_nonbridge_base_present": f152["package_level_nonbridge_base_present"],
        "positive_bridge_branch_still_open": f152["positive_bridge_branch_still_open"],
        "actual_nonbridge_strengthening_discharged": f152["actual_nonbridge_strengthening_discharged"],
        "permanent_no_bridge_claimed": f152["permanent_no_bridge_claimed"],
        "kernel_split_safe": f152["kernel_split_safe"],
        "no_false_pass": True,
    }
    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
