#!/usr/bin/env python3

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
GENERATED = ROOT / "generated"
OUT = (
    GENERATED
    / "p662_current_actual_legacy_to_strict_kernel_nonbridge_strengthening_discharge_witness_probe_summary.json"
)


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    f662 = load_json(
        GENERATED
        / "f662_current_actual_legacy_to_strict_kernel_nonbridge_strengthening_discharge_witness_packet_summary.json"
    )

    comps = f662.get("components", {})
    checks_spec = [
        {
            "id": "package_level_nonbridge_present",
            "actual": bool(comps.get("package_level_nonbridge_on_current_repo_state")),
            "expected": True,
        },
        {
            "id": "package_nontransfer_present",
            "actual": bool(comps.get("legacy_to_strict_package_nontransfer_on_current_repo_state")),
            "expected": True,
        },
        {
            "id": "a_abs_obstruction_discharged",
            "actual": bool(comps.get("a_abs_nonbridge_obstruction_discharged")),
            "expected": True,
        },
        {
            "id": "r_damp_obstruction_discharged",
            "actual": bool(comps.get("r_damp_nonbridge_obstruction_discharged")),
            "expected": True,
        },
        {
            "id": "phase_frequency_obstruction_discharged",
            "actual": bool(comps.get("actual_phase_frequency_nonconformal_obstruction_discharged")),
            "expected": True,
        },
        {
            "id": "actual_nonbridge_strengthening_discharged",
            "actual": bool(f662.get("actual_nonbridge_strengthening_discharged")),
            "expected": True,
        },
        {
            "id": "permanent_no_bridge_NOT_claimed",
            "actual": bool(f662.get("permanent_no_bridge_claimed")),
            "expected": False,
        },
        {
            "id": "kernel_split_safe",
            "actual": bool(f662.get("kernel_split_safe")),
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
        status = "P662_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_NONBRIDGE_STRENGTHENING_DISCHARGE_STATE"
    else:
        status = "CURRENT_REPO_EXPORTS_ONE_ACTUAL_LEGACY_TO_STRICT_KERNEL_NONBRIDGE_STRENGTHENING_DISCHARGE_WITNESS_AFTER_P662"

    summary = {
        "stage": "P662",
        "lane": "legacy_strict_kernel_nonbridge_strengthening",
        "status": status,
        "checks": checks,
        "blocking_mismatches": mismatches,
        "actual_nonbridge_strengthening_discharged": bool(
            f662.get("actual_nonbridge_strengthening_discharged")
        ),
        "kernel_split_safe": bool(f662.get("kernel_split_safe")),
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()

