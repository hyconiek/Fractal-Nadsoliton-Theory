#!/usr/bin/env python3

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "p243_current_actual_legacy_to_strict_kernel_bifurcated_frontier_probe_summary.json"


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    f153 = load_json(
        GENERATED / "f153_first_actual_legacy_to_strict_kernel_bifurcated_frontier_packet_summary.json"
    )

    checks_spec = [
        {
            "id": "bridge_target_exported",
            "actual": f153["bridge_target_exported"],
            "expected": True,
        },
        {
            "id": "bridge_branch_future_only",
            "actual": f153["bridge_branch_future_only"],
            "expected": True,
        },
        {
            "id": "nonbridge_target_exported",
            "actual": f153["nonbridge_target_exported"],
            "expected": True,
        },
        {
            "id": "nonbridge_branch_future_only",
            "actual": f153["nonbridge_branch_future_only"],
            "expected": True,
        },
        {
            "id": "actual_bridge_NOT_discharged",
            "actual": f153["actual_bridge_discharged"],
            "expected": False,
        },
        {
            "id": "actual_nonbridge_strengthening_NOT_discharged",
            "actual": f153["actual_nonbridge_strengthening_discharged"],
            "expected": False,
        },
        {
            "id": "branch_selection_NOT_justified",
            "actual": f153["branch_selection_justified_on_current_repo_state"],
            "expected": False,
        },
        {
            "id": "kernel_split_safe",
            "actual": f153["kernel_split_safe"],
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
        status = "P243_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_LEGACY_STRICT_KERNEL_FRONTIER_STATE"
    else:
        status = "CURRENT_REPO_EXPORTS_ONE_ACTUAL_LEGACY_TO_STRICT_KERNEL_BIFURCATED_FRONTIER_PACKET_WITH_NO_JUSTIFIED_BRANCH_SELECTION_AFTER_P243"

    summary = {
        "stage": "P243",
        "lane": "legacy_strict_kernel_frontier_state_only",
        "status": status,
        "checks": checks,
        "blocking_mismatches": mismatches,
        "actual_bifurcated_frontier_packet_exported": f153["actual_bifurcated_frontier_packet_exported"],
        "bridge_target_exported": f153["bridge_target_exported"],
        "nonbridge_target_exported": f153["nonbridge_target_exported"],
        "bridge_branch_future_only": f153["bridge_branch_future_only"],
        "nonbridge_branch_future_only": f153["nonbridge_branch_future_only"],
        "branch_selection_justified_on_current_repo_state": f153["branch_selection_justified_on_current_repo_state"],
        "actual_bridge_discharged": f153["actual_bridge_discharged"],
        "actual_nonbridge_strengthening_discharged": f153["actual_nonbridge_strengthening_discharged"],
        "kernel_split_safe": f153["kernel_split_safe"],
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
