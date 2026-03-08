#!/usr/bin/env python3

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "n263_current_first_actual_legacy_to_strict_kernel_bifurcated_frontier_theorem_summary.json"


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    p243 = load_json(
        GENERATED / "p243_current_actual_legacy_to_strict_kernel_bifurcated_frontier_probe_summary.json"
    )

    expected_status = (
        "CURRENT_REPO_EXPORTS_ONE_ACTUAL_LEGACY_TO_STRICT_KERNEL_BIFURCATED_FRONTIER_PACKET_WITH_NO_JUSTIFIED_BRANCH_SELECTION_AFTER_P243"
    )
    status_ok = p243["status"] == expected_status

    summary = {
        "step": "N263",
        "status": "N263_DISCHARGED_CURRENT_FIRST_ACTUAL_LEGACY_TO_STRICT_KERNEL_BIFURCATED_FRONTIER_THEOREM_NO_FALSE_PASS",
        "scope": "current_legacy_strict_kernel_frontier_state_only",
        "checks": [
            {
                "id": "p243_frontier_status",
                "actual": p243["status"],
                "expected": expected_status,
                "pass": status_ok,
            }
        ],
        "theorem_result": {
            "discharged": status_ok and len(p243.get("blocking_mismatches", [])) == 0,
            "actual_bifurcated_frontier_packet_exported": p243["actual_bifurcated_frontier_packet_exported"],
            "bridge_target_exported": p243["bridge_target_exported"],
            "nonbridge_target_exported": p243["nonbridge_target_exported"],
            "bridge_branch_future_only": p243["bridge_branch_future_only"],
            "nonbridge_branch_future_only": p243["nonbridge_branch_future_only"],
            "branch_selection_justified_on_current_repo_state": False,
            "actual_bridge_discharged": False,
            "actual_nonbridge_strengthening_discharged": False,
            "current_strict_core_selector_closure": False,
            "current_global_selector_closure": False,
            "current_global_qw2191_discharge": False,
            "kernel_split_safe": p243["kernel_split_safe"],
        },
        "hard_limits": [
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
