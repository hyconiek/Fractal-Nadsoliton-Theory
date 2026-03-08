#!/usr/bin/env python3

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
IN_N123 = ROOT / "generated" / "n123_current_legacy_to_strict_kernel_package_level_nonbridge_theorem_summary.json"
IN_N261 = ROOT / "generated" / "n261_current_first_legacy_to_strict_kernel_bridge_target_theorem_summary.json"
OUT = ROOT / "generated" / "f152_first_legacy_to_strict_kernel_nonbridge_strengthening_target_packet_summary.json"


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    n123 = load_json(IN_N123)
    n261 = load_json(IN_N261)

    summary = {
        "packet_id": "F152",
        "status": "F152_EXECUTED_FIRST_LEGACY_TO_STRICT_KERNEL_NONBRIDGE_STRENGTHENING_TARGET_PACKET_FUTURE_ROUTE_ONLY_NO_FALSE_PASS",
        "as_of": "2026-03-08",
        "package_level_nonbridge_base_present": n123["theorem_result"]["package_level_nonbridge_on_current_repo_state"],
        "input_legacy_kernel": "K_legacy_ont",
        "input_strict_kernel": "K_strict_gate",
        "nonbridge_strengthening_target": "NB_legacy_strict_strengthening_target_v1",
        "future_only_target_exported": True,
        "nonbridge_scope": "structural_kernel_comparison_negative_branch_only",
        "components_target": {
            "name": "Delta_nonbridge_components_target_v1",
            "components": [
                "A_abs_nonbridge_obstruction_target_v1",
                "R_damp_nonbridge_obstruction_target_v1",
                "P_shift_nonbridge_obstruction_target_v1",
            ],
        },
        "positive_bridge_branch_still_open": n261["theorem_result"]["bridge_target_exported"],
        "actual_nonbridge_strengthening_discharged": False,
        "permanent_no_bridge_claimed": False,
        "actual_bridge_discharged": False,
        "current_strict_core_selector_closure": False,
        "current_global_selector_closure": False,
        "current_global_qw2191_discharge": False,
        "kernel_split_safe": True,
        "no_false_pass": True,
    }
    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
