#!/usr/bin/env python3

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "p245_current_actual_legacy_to_strict_kernel_claim_specific_amplitude_nonabsorption_witness_probe_summary.json"


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    f155 = load_json(
        GENERATED / "f155_first_actual_legacy_to_strict_kernel_claim_specific_amplitude_nonabsorption_witness_packet_summary.json"
    )

    checks_spec = [
        {
            "id": "actual_component_witness_exported",
            "actual": f155["actual_component_witness_exported"],
            "expected": True,
        },
        {
            "id": "claim_specific_witness_exported",
            "actual": f155["actual_claim_specific_amplitude_nonabsorption_witness_exported"],
            "expected": True,
        },
        {
            "id": "claim_specific_scope_preserved",
            "actual": f155["claim_specific_scope"],
            "expected": "legacy_weinberg_amplitude_role_only",
        },
        {
            "id": "frontier_undecided",
            "actual": f155["bridge_nonbridge_frontier_undecided"],
            "expected": True,
        },
        {
            "id": "full_amplitude_obstruction_NOT_discharged",
            "actual": f155["full_amplitude_nonabsorption_obstruction_discharged"],
            "expected": False,
        },
        {
            "id": "kernel_split_safe",
            "actual": f155["kernel_split_safe"],
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
        status = "P245_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_CLAIM_SPECIFIC_AMPLITUDE_WITNESS_STATE"
    else:
        status = "CURRENT_REPO_EXPORTS_ONE_ACTUAL_CLAIM_SPECIFIC_LEGACY_TO_STRICT_KERNEL_AMPLITUDE_NONABSORPTION_WITNESS_BELOW_FULL_AMPLITUDE_OBSTRUCTION_AFTER_P245"

    summary = {
        "stage": "P245",
        "lane": "legacy_strict_kernel_nonbridge_claim_specific_amplitude_only",
        "status": status,
        "checks": checks,
        "blocking_mismatches": mismatches,
        "actual_claim_specific_amplitude_nonabsorption_witness_exported": f155["actual_claim_specific_amplitude_nonabsorption_witness_exported"],
        "claim_specific_scope": f155["claim_specific_scope"],
        "bridge_nonbridge_frontier_undecided": f155["bridge_nonbridge_frontier_undecided"],
        "full_amplitude_nonabsorption_obstruction_discharged": f155["full_amplitude_nonabsorption_obstruction_discharged"],
        "actual_nonbridge_strengthening_discharged": f155["actual_nonbridge_strengthening_discharged"],
        "kernel_split_safe": f155["kernel_split_safe"],
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
