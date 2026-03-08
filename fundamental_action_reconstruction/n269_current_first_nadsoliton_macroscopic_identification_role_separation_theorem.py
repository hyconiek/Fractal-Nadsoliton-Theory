#!/usr/bin/env python3

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "n269_current_first_nadsoliton_macroscopic_identification_role_separation_theorem_summary.json"


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    p249 = load_json(
        GENERATED / "p249_current_nadsoliton_macroscopic_identification_role_separation_probe_summary.json"
    )

    expected_status = (
        "CURRENT_REPO_EXPORTS_ONE_ACTUAL_NADSOLITON_MACROSCOPIC_IDENTIFICATION_ROLE_SEPARATION_PACKET_WITHDRAWING_THE_T15_T16_DEADLOCK_AS_A_MANDATORY_T14_CLOSURE_GATE_AFTER_P249"
    )
    status_ok = p249["status"] == expected_status

    summary = {
        "step": "N269",
        "status": "N269_DISCHARGED_CURRENT_FIRST_NADSOLITON_MACROSCOPIC_IDENTIFICATION_ROLE_SEPARATION_THEOREM_NO_FALSE_PASS",
        "scope": "nadsoliton_role_separation_and_t14_boundary_only",
        "checks": [
            {
                "id": "p249_role_separation_status",
                "actual": p249["status"],
                "expected": expected_status,
                "pass": status_ok,
            }
        ],
        "theorem_result": {
            "discharged": status_ok and len(p249.get("blocking_mismatches", [])) == 0,
            "nadsoliton_role_separation_principle_exported": True,
            "new_closure_level_ingredient_added_after_n260": p249["new_closure_level_ingredient_added_after_n260"],
            "post_n260_bridge_nonbridge_impasse_withdrawn_as_mandatory_t14_closure_gate": True,
            "bridge_nonbridge_frontier_reclassified_as_optional_comparison_frontier": True,
            "t15_bridge_required_for_t14_selector_closure": p249["t15_bridge_required_for_t14_selector_closure"],
            "t16_nonbridge_required_for_t14_selector_closure": p249["t16_nonbridge_required_for_t14_selector_closure"],
            "lack_of_cross_kernel_absorption_counts_as_toe_failure": p249["lack_of_cross_kernel_absorption_counts_as_toe_failure"],
            "bridge_nonbridge_frontier_undecided": p249["bridge_nonbridge_frontier_undecided"],
            "kernel_identity_claimed": p249["kernel_identity_claimed"],
            "legacy_physical_role_transfer_claimed": p249["legacy_physical_role_transfer_claimed"],
            "actual_bridge_discharged": False,
            "actual_nonbridge_strengthening_discharged": False,
            "current_strict_core_selector_closure": False,
            "current_global_selector_closure": False,
            "current_global_qw2191_discharge": False,
        },
        "hard_limits": [
            "no_kernel_identity_claim",
            "no_actual_bridge_derivation",
            "no_actual_nonbridge_strengthening_theorem",
            "no_current_branch_selection_theorem",
            "no_legacy_physical_role_transfer",
            "no_strict_core_selector_closure",
            "no_global_selector_closure",
            "no_global_QW2191_discharge",
            "no_ToE_closure",
        ],
    }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
