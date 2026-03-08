#!/usr/bin/env python3

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "p249_current_nadsoliton_macroscopic_identification_role_separation_probe_summary.json"


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    f159 = load_json(
        GENERATED / "f159_first_actual_nadsoliton_macroscopic_identification_role_separation_packet_summary.json"
    )

    checks_spec = [
        {
            "id": "single_nadsoliton_ontology_respected",
            "actual": f159["single_nadsoliton_ontology_respected"],
            "expected": True,
        },
        {
            "id": "legacy_kernel_role_macro_identification_only",
            "actual": f159["legacy_kernel_role"],
            "expected": "macroscopic_nadsoliton_identification_tool_only",
        },
        {
            "id": "strict_kernel_role_source_topology_only",
            "actual": f159["strict_kernel_role"],
            "expected": "strict_source_topology_working_kernel_only",
        },
        {
            "id": "cross_kernel_absorption_not_required_for_toe_closure",
            "actual": f159["cross_kernel_absorption_required_for_toe_closure"],
            "expected": False,
        },
        {
            "id": "lack_of_cross_kernel_absorption_not_toe_failure",
            "actual": f159["lack_of_cross_kernel_absorption_counts_as_toe_failure"],
            "expected": False,
        },
        {
            "id": "t15_bridge_not_required_for_t14_closure",
            "actual": f159["t15_bridge_required_for_t14_selector_closure"],
            "expected": False,
        },
        {
            "id": "t16_nonbridge_not_required_for_t14_closure",
            "actual": f159["t16_nonbridge_required_for_t14_selector_closure"],
            "expected": False,
        },
        {
            "id": "kernel_identity_not_claimed",
            "actual": f159["kernel_identity_claimed"],
            "expected": False,
        },
        {
            "id": "legacy_role_transfer_not_claimed",
            "actual": f159["legacy_physical_role_transfer_claimed"],
            "expected": False,
        },
        {
            "id": "bridge_nonbridge_frontier_still_undecided",
            "actual": f159["bridge_nonbridge_frontier_undecided"],
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
        status = "P249_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_NADSOLITON_ROLE_SEPARATION_STATE"
    else:
        status = (
            "CURRENT_REPO_EXPORTS_ONE_ACTUAL_NADSOLITON_MACROSCOPIC_IDENTIFICATION_ROLE_SEPARATION_PACKET_WITHDRAWING_THE_T15_T16_DEADLOCK_AS_A_MANDATORY_T14_CLOSURE_GATE_AFTER_P249"
        )

    summary = {
        "stage": "P249",
        "lane": "nadsoliton_role_separation_boundary_only",
        "status": status,
        "checks": checks,
        "blocking_mismatches": mismatches,
        "new_closure_level_ingredient_added_after_n260": f159["new_closure_level_ingredient_added_after_n260"],
        "bridge_nonbridge_frontier_explicit": f159["bridge_nonbridge_frontier_explicit"],
        "bridge_nonbridge_frontier_undecided": f159["bridge_nonbridge_frontier_undecided"],
        "t15_bridge_required_for_t14_selector_closure": f159["t15_bridge_required_for_t14_selector_closure"],
        "t16_nonbridge_required_for_t14_selector_closure": f159["t16_nonbridge_required_for_t14_selector_closure"],
        "lack_of_cross_kernel_absorption_counts_as_toe_failure": f159["lack_of_cross_kernel_absorption_counts_as_toe_failure"],
        "kernel_identity_claimed": f159["kernel_identity_claimed"],
        "legacy_physical_role_transfer_claimed": f159["legacy_physical_role_transfer_claimed"],
        "current_strict_core_selector_closure": f159["current_strict_core_selector_closure"],
        "current_global_qw2191_discharge": f159["current_global_qw2191_discharge"],
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
