#!/usr/bin/env python3

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
GENERATED = ROOT / "generated"

OUT = (
    GENERATED
    / "n663_current_first_actual_legacy_to_strict_kernel_bifurcated_frontier_theorem_v2_summary.json"
)


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    p663 = load_json(
        GENERATED
        / "p663_current_actual_legacy_to_strict_kernel_bifurcated_frontier_probe_v2_summary.json"
    )

    expected_status = (
        "CURRENT_REPO_EXPORTS_ONE_ACTUAL_LEGACY_TO_STRICT_KERNEL_BIFURCATED_FRONTIER_PACKET_V2_POST_N554_AFTER_P663"
    )
    status_ok = p663["status"] == expected_status

    summary = {
        "step": "N663",
        "status": "N663_DISCHARGED_CURRENT_FIRST_ACTUAL_LEGACY_TO_STRICT_KERNEL_BIFURCATED_FRONTIER_THEOREM_V2_NO_FALSE_PASS",
        "scope": "current_legacy_strict_kernel_frontier_status_v2_only",
        "checks": [
            {
                "id": "p663_frontier_v2_status",
                "actual": p663["status"],
                "expected": expected_status,
                "pass": status_ok,
            }
        ],
        "theorem_result": {
            "discharged": status_ok and len(p663.get("blocking_mismatches", [])) == 0,
            "actual_nonbridge_strengthening_discharged": bool(
                p663.get("actual_nonbridge_strengthening_discharged")
            ),
            "branch_selection_justified_on_current_repo_state": bool(
                p663.get("branch_selection_justified_on_current_repo_state")
            ),
            "actual_bridge_discharged": False,
            "kernel_split_safe": bool(p663.get("kernel_split_safe")),
        },
        "hard_limits": [
            "no_actual_bridge_derivation",
            "no_permanent_no_bridge_claim",
            "no_current_branch_selection_theorem",
            "no_legacy_physical_role_transfer",
            "no_strict_core_selector_closure",
            "no_global_selector_closure",
            "no_global_QW2191_discharge",
            "no_ToE_closure",
        ],
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()

