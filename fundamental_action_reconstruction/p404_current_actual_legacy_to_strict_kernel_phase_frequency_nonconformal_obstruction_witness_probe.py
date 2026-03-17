#!/usr/bin/env python3

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
GENERATED = ROOT / "generated"
OUT = (
    GENERATED
    / "p404_current_actual_legacy_to_strict_kernel_phase_frequency_nonconformal_obstruction_witness_probe_summary.json"
)


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    f326 = load_json(
        GENERATED
        / "f326_first_actual_legacy_to_strict_kernel_phase_frequency_nonconformal_obstruction_witness_packet_summary.json"
    )

    checks_spec = [
        {
            "id": "explicit_phase_frequency_bridge_present",
            "actual": f326["explicit_phase_frequency_bridge_present"],
            "expected": False,
        },
        {
            "id": "package_nontransfer_present",
            "actual": f326["legacy_to_strict_package_nontransfer_on_current_repo_state"],
            "expected": True,
        },
        {
            "id": "a_abs_obstruction_discharged",
            "actual": f326["a_abs_nonbridge_obstruction_discharged"],
            "expected": True,
        },
        {
            "id": "r_damp_obstruction_discharged",
            "actual": f326["r_damp_nonbridge_obstruction_discharged"],
            "expected": True,
        },
        {
            "id": "p_shift_obstruction_discharged",
            "actual": f326["p_shift_nonbridge_obstruction_discharged"],
            "expected": True,
        },
        {
            "id": "kernel_split_safe",
            "actual": f326["kernel_split_safe"],
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
        status = "P404_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_PHASE_FREQUENCY_OBSTRUCTION_STATE"
    else:
        status = "CURRENT_REPO_EXPORTS_ONE_ACTUAL_LEGACY_TO_STRICT_KERNEL_PHASE_FREQUENCY_NONCONFORMAL_OBSTRUCTION_WITNESS_AFTER_P404"

    summary = {
        "stage": "P404",
        "lane": "legacy_strict_kernel_nonbridge_phase_frequency_only",
        "status": status,
        "checks": checks,
        "blocking_mismatches": mismatches,
        "explicit_phase_frequency_bridge_present": f326["explicit_phase_frequency_bridge_present"],
        "p_shift_nonbridge_obstruction_discharged": f326["p_shift_nonbridge_obstruction_discharged"],
        "kernel_split_safe": f326["kernel_split_safe"],
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()

