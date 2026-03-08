#!/usr/bin/env python3

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "p248_current_actual_legacy_to_strict_kernel_damping_nonrenormalization_obstruction_witness_probe_summary.json"


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    f158 = load_json(
        GENERATED / "f158_first_actual_legacy_to_strict_kernel_damping_nonrenormalization_obstruction_witness_packet_summary.json"
    )

    checks_spec = [
        {
            "id": "explicit_beta_tors_to_beta_eta_translation_absent",
            "actual": f158["explicit_beta_tors_to_beta_eta_translation_present"],
            "expected": False,
        },
        {
            "id": "package_nontransfer_present",
            "actual": f158["legacy_to_strict_package_nontransfer_on_current_repo_state"],
            "expected": True,
        },
        {
            "id": "a_abs_obstruction_discharged",
            "actual": f158["a_abs_nonbridge_obstruction_discharged"],
            "expected": True,
        },
        {
            "id": "r_damp_obstruction_discharged",
            "actual": f158["r_damp_nonbridge_obstruction_discharged"],
            "expected": True,
        },
        {
            "id": "phase_not_discharged",
            "actual": f158["actual_phase_frequency_nonconformal_obstruction_discharged"],
            "expected": False,
        },
        {
            "id": "kernel_split_safe",
            "actual": f158["kernel_split_safe"],
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
        status = "P248_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_FULL_R_DAMP_OBSTRUCTION_STATE"
    else:
        status = "CURRENT_REPO_EXPORTS_ONE_ACTUAL_LEGACY_TO_STRICT_KERNEL_DAMPING_NONRENORMALIZATION_OBSTRUCTION_WITNESS_AFTER_P248"

    summary = {
        "stage": "P248",
        "lane": "legacy_strict_kernel_nonbridge_full_damping_only",
        "status": status,
        "checks": checks,
        "blocking_mismatches": mismatches,
        "r_damp_nonbridge_obstruction_discharged": f158["r_damp_nonbridge_obstruction_discharged"],
        "legacy_to_strict_package_nontransfer_on_current_repo_state": f158["legacy_to_strict_package_nontransfer_on_current_repo_state"],
        "a_abs_nonbridge_obstruction_discharged": f158["a_abs_nonbridge_obstruction_discharged"],
        "actual_phase_frequency_nonconformal_obstruction_discharged": f158["actual_phase_frequency_nonconformal_obstruction_discharged"],
        "actual_nonbridge_strengthening_discharged": f158["actual_nonbridge_strengthening_discharged"],
        "kernel_split_safe": f158["kernel_split_safe"],
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
