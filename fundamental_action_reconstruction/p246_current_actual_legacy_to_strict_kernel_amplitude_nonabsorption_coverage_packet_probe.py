#!/usr/bin/env python3

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "p246_current_actual_legacy_to_strict_kernel_amplitude_nonabsorption_coverage_packet_probe_summary.json"


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    f156 = load_json(
        GENERATED / "f156_first_actual_legacy_to_strict_kernel_amplitude_nonabsorption_coverage_packet_summary.json"
    )

    checks_spec = [
        {
            "id": "legacy_weinberg_frontier_closed",
            "actual": f156["legacy_weinberg_claim_specific_frontier_closed"],
            "expected": True,
        },
        {
            "id": "legacy_fine_structure_frontier_closed",
            "actual": f156["legacy_fine_structure_claim_specific_frontier_closed"],
            "expected": True,
        },
        {
            "id": "legacy_gravity_hierarchy_frontier_closed",
            "actual": f156["legacy_gravity_hierarchy_claim_specific_frontier_closed"],
            "expected": True,
        },
        {
            "id": "legacy_physical_role_package_closed",
            "actual": f156["legacy_physical_role_package_closed_negatively"],
            "expected": True,
        },
        {
            "id": "actual_claim_specific_amplitude_witness_exported",
            "actual": f156["actual_claim_specific_amplitude_witness_exported"],
            "expected": True,
        },
        {
            "id": "coverage_packet_exported",
            "actual": f156["actual_amplitude_nonabsorption_coverage_packet_exported"],
            "expected": True,
        },
        {
            "id": "full_a_abs_obstruction_NOT_discharged",
            "actual": f156["full_a_abs_obstruction_discharged"],
            "expected": False,
        },
        {
            "id": "kernel_split_safe",
            "actual": f156["kernel_split_safe"],
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
        status = "P246_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_AMPLITUDE_COVERAGE_STATE"
    else:
        status = "CURRENT_REPO_EXPORTS_ONE_ACTUAL_LEGACY_TO_STRICT_KERNEL_AMPLITUDE_NONABSORPTION_COVERAGE_PACKET_BELOW_FULL_A_ABS_OBSTRUCTION_AFTER_P246"

    summary = {
        "stage": "P246",
        "lane": "legacy_strict_kernel_nonbridge_amplitude_coverage_only",
        "status": status,
        "checks": checks,
        "blocking_mismatches": mismatches,
        "actual_amplitude_nonabsorption_coverage_packet_exported": f156["actual_amplitude_nonabsorption_coverage_packet_exported"],
        "legacy_physical_role_package_closed_negatively": f156["legacy_physical_role_package_closed_negatively"],
        "actual_claim_specific_amplitude_witness_exported": f156["actual_claim_specific_amplitude_witness_exported"],
        "full_a_abs_obstruction_discharged": f156["full_a_abs_obstruction_discharged"],
        "actual_nonbridge_strengthening_discharged": f156["actual_nonbridge_strengthening_discharged"],
        "kernel_split_safe": f156["kernel_split_safe"],
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
