#!/usr/bin/env python3

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "p247_current_actual_legacy_to_strict_kernel_full_amplitude_nonabsorption_obstruction_witness_probe_summary.json"


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    f157 = load_json(
        GENERATED / "f157_first_actual_legacy_to_strict_kernel_full_amplitude_nonabsorption_obstruction_witness_packet_summary.json"
    )

    checks_spec = [
        {
            "id": "explicit_amplitude_absorption_map_for_alpha_geo_absent",
            "actual": f157["explicit_amplitude_absorption_map_for_alpha_geo_present"],
            "expected": False,
        },
        {
            "id": "package_nontransfer_present",
            "actual": f157["legacy_to_strict_package_nontransfer_on_current_repo_state"],
            "expected": True,
        },
        {
            "id": "amplitude_coverage_packet_exported",
            "actual": f157["actual_amplitude_nonabsorption_coverage_packet_exported"],
            "expected": True,
        },
        {
            "id": "a_abs_obstruction_discharged",
            "actual": f157["a_abs_nonbridge_obstruction_discharged"],
            "expected": True,
        },
        {
            "id": "damping_not_discharged",
            "actual": f157["actual_damping_nonrenormalization_obstruction_discharged"],
            "expected": False,
        },
        {
            "id": "phase_not_discharged",
            "actual": f157["actual_phase_frequency_nonconformal_obstruction_discharged"],
            "expected": False,
        },
        {
            "id": "kernel_split_safe",
            "actual": f157["kernel_split_safe"],
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
        status = "P247_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_FULL_A_ABS_OBSTRUCTION_STATE"
    else:
        status = "CURRENT_REPO_EXPORTS_ONE_ACTUAL_LEGACY_TO_STRICT_KERNEL_FULL_AMPLITUDE_NONABSORPTION_OBSTRUCTION_WITNESS_AFTER_P247"

    summary = {
        "stage": "P247",
        "lane": "legacy_strict_kernel_nonbridge_full_amplitude_only",
        "status": status,
        "checks": checks,
        "blocking_mismatches": mismatches,
        "a_abs_nonbridge_obstruction_discharged": f157["a_abs_nonbridge_obstruction_discharged"],
        "legacy_to_strict_package_nontransfer_on_current_repo_state": f157["legacy_to_strict_package_nontransfer_on_current_repo_state"],
        "actual_amplitude_nonabsorption_coverage_packet_exported": f157["actual_amplitude_nonabsorption_coverage_packet_exported"],
        "actual_damping_nonrenormalization_obstruction_discharged": f157["actual_damping_nonrenormalization_obstruction_discharged"],
        "actual_phase_frequency_nonconformal_obstruction_discharged": f157["actual_phase_frequency_nonconformal_obstruction_discharged"],
        "actual_nonbridge_strengthening_discharged": f157["actual_nonbridge_strengthening_discharged"],
        "kernel_split_safe": f157["kernel_split_safe"],
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
