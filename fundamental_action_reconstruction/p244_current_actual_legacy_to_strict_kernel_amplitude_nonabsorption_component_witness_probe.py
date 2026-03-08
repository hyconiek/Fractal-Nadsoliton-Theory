#!/usr/bin/env python3

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "p244_current_actual_legacy_to_strict_kernel_amplitude_nonabsorption_component_witness_probe_summary.json"


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    f154 = load_json(
        GENERATED / "f154_first_actual_legacy_to_strict_kernel_amplitude_nonabsorption_component_witness_packet_summary.json"
    )

    checks_spec = [
        {
            "id": "actual_component_witness_exported",
            "actual": f154["actual_amplitude_nonabsorption_component_witness_exported"],
            "expected": True,
        },
        {
            "id": "legacy_weinberg_role_transfer_absent",
            "actual": f154["legacy_weinberg_role_transfer_present"],
            "expected": False,
        },
        {
            "id": "strict_side_candidate_object_present",
            "actual": f154["strict_side_candidate_object_present"],
            "expected": True,
        },
        {
            "id": "explicit_role_equivalence_verdict_absent",
            "actual": f154["explicit_role_equivalence_verdict_present"],
            "expected": False,
        },
        {
            "id": "legacy_weinberg_frontier_closed_negatively",
            "actual": f154["legacy_weinberg_claim_specific_frontier_closed_negatively"],
            "expected": True,
        },
        {
            "id": "frontier_undecided",
            "actual": f154["bridge_nonbridge_frontier_undecided"],
            "expected": True,
        },
        {
            "id": "full_amplitude_obstruction_NOT_discharged",
            "actual": f154["actual_amplitude_nonabsorption_obstruction_discharged"],
            "expected": False,
        },
        {
            "id": "kernel_split_safe",
            "actual": f154["kernel_split_safe"],
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
        status = "P244_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_AMPLITUDE_NONABSORPTION_COMPONENT_STATE"
    else:
        status = "CURRENT_REPO_EXPORTS_ONE_ACTUAL_LEGACY_TO_STRICT_KERNEL_AMPLITUDE_NONABSORPTION_COMPONENT_WITNESS_BELOW_FULL_AMPLITUDE_OBSTRUCTION_AFTER_P244"

    summary = {
        "stage": "P244",
        "lane": "legacy_strict_kernel_nonbridge_amplitude_component_only",
        "status": status,
        "checks": checks,
        "blocking_mismatches": mismatches,
        "actual_amplitude_nonabsorption_component_witness_exported": f154["actual_amplitude_nonabsorption_component_witness_exported"],
        "legacy_weinberg_role_transfer_present": f154["legacy_weinberg_role_transfer_present"],
        "strict_side_candidate_object_present": f154["strict_side_candidate_object_present"],
        "explicit_role_equivalence_verdict_present": f154["explicit_role_equivalence_verdict_present"],
        "legacy_weinberg_claim_specific_frontier_closed_negatively": f154["legacy_weinberg_claim_specific_frontier_closed_negatively"],
        "bridge_nonbridge_frontier_undecided": f154["bridge_nonbridge_frontier_undecided"],
        "actual_amplitude_nonabsorption_obstruction_discharged": f154["actual_amplitude_nonabsorption_obstruction_discharged"],
        "actual_nonbridge_strengthening_discharged": f154["actual_nonbridge_strengthening_discharged"],
        "kernel_split_safe": f154["kernel_split_safe"],
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
