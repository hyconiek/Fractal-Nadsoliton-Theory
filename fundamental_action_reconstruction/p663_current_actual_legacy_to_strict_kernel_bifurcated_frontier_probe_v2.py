#!/usr/bin/env python3

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
GENERATED = ROOT / "generated"

OUT = (
    GENERATED
    / "p663_current_actual_legacy_to_strict_kernel_bifurcated_frontier_probe_v2_summary.json"
)


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    f663 = load_json(
        GENERATED
        / "f663_current_actual_legacy_to_strict_kernel_bifurcated_frontier_packet_v2_summary.json"
    )

    checks_spec = [
        {
            "id": "bridge_target_exported",
            "actual": f663["bridge_target_exported"],
            "expected": True,
        },
        {
            "id": "bridge_branch_future_only",
            "actual": f663["bridge_branch_future_only"],
            "expected": True,
        },
        {
            "id": "nonbridge_target_exported",
            "actual": f663["nonbridge_target_exported"],
            "expected": True,
        },
        {
            "id": "nonbridge_branch_NOT_future_only",
            "actual": f663["nonbridge_branch_future_only"],
            "expected": False,
        },
        {
            "id": "actual_bridge_NOT_discharged",
            "actual": f663["actual_bridge_discharged"],
            "expected": False,
        },
        {
            "id": "actual_nonbridge_strengthening_discharged",
            "actual": f663["actual_nonbridge_strengthening_discharged"],
            "expected": True,
        },
        {
            "id": "branch_selection_NOT_justified",
            "actual": f663["branch_selection_justified_on_current_repo_state"],
            "expected": False,
        },
        {
            "id": "kernel_split_safe",
            "actual": f663["kernel_split_safe"],
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
        status = "P663_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_KERNEL_FRONTIER_V2_STATE"
    else:
        status = "CURRENT_REPO_EXPORTS_ONE_ACTUAL_LEGACY_TO_STRICT_KERNEL_BIFURCATED_FRONTIER_PACKET_V2_POST_N554_AFTER_P663"

    summary = {
        "stage": "P663",
        "lane": "legacy_strict_kernel_frontier_status_v2",
        "status": status,
        "checks": checks,
        "blocking_mismatches": mismatches,
        "actual_nonbridge_strengthening_discharged": f663["actual_nonbridge_strengthening_discharged"],
        "branch_selection_justified_on_current_repo_state": f663["branch_selection_justified_on_current_repo_state"],
        "kernel_split_safe": f663["kernel_split_safe"],
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()

