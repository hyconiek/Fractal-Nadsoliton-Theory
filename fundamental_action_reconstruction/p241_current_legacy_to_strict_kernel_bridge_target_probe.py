#!/usr/bin/env python3

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "p241_current_legacy_to_strict_kernel_bridge_target_probe_summary.json"

def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))

def main() -> None:
    f151 = load_json(
        GENERATED / "f151_first_legacy_to_strict_kernel_bridge_target_packet_summary.json"
    )

    checks_spec = [
        {
            "id": "bridge_target_exported",
            "actual": f151["bridge_target"] == "B_legacy_strict_bridge_target_v1",
            "expected": True,
        },
        {
            "id": "split_safety_maintained",
            "actual": f151["kernel_split_safe"],
            "expected": True,
        },
        {
            "id": "future_only_target_exported",
            "actual": f151["future_only_target_exported"],
            "expected": True,
        },
        {
            "id": "nonbridge_branch_still_open",
            "actual": f151["nonbridge_branch_still_open"],
            "expected": True,
        },
        {
            "id": "legacy_to_strict_bridge_NOT_claimed",
            "actual": f151["legacy_to_strict_bridge_claimed"],
            "expected": False,
        },
        {
            "id": "legacy_physical_role_transfer_NOT_claimed",
            "actual": f151["legacy_physical_role_transfer_claimed"],
            "expected": False,
        },
        {
            "id": "bridge_sufficiency_for_global_qw2191_NOT_claimed",
            "actual": f151["bridge_sufficiency_for_global_qw2191_claimed"],
            "expected": False,
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
        status = "P241_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_BRIDGE_TARGET_STATE"
    else:
        status = "CURRENT_REPO_EXPORTS_ONE_FUTURE_ONLY_LEGACY_TO_STRICT_KERNEL_BRIDGE_TARGET_BELOW_ACTUAL_BRIDGE_DISCHARGE_AFTER_P241"

    summary = {
        "stage": "P241",
        "lane": "kernel_bridge",
        "status": status,
        "checks": checks,
        "blocking_mismatches": mismatches,
        "bridge_target_exported": f151["bridge_target"] == "B_legacy_strict_bridge_target_v1",
        "future_only_target_exported": f151["future_only_target_exported"],
        "kernel_split_safe": f151["kernel_split_safe"],
        "nonbridge_branch_still_open": f151["nonbridge_branch_still_open"],
        "legacy_to_strict_bridge_claimed": f151["legacy_to_strict_bridge_claimed"],
        "actual_bridge_discharged": f151["actual_bridge_discharged"],
        "legacy_physical_role_transfer_claimed": f151["legacy_physical_role_transfer_claimed"],
        "current_strict_core_selector_closure": f151["current_strict_core_selector_closure"],
        "current_global_selector_closure": f151["current_global_selector_closure"],
        "current_global_qw2191_discharge": f151["current_global_qw2191_discharge"],
        "bridge_sufficiency_for_global_qw2191_claimed": f151["bridge_sufficiency_for_global_qw2191_claimed"],
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)

if __name__ == "__main__":
    main()
