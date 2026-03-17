#!/usr/bin/env python3

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
GENERATED = ROOT / "generated"

IN_N261 = GENERATED / "n261_current_first_legacy_to_strict_kernel_bridge_target_theorem_summary.json"
IN_N262 = GENERATED / "n262_current_first_legacy_to_strict_kernel_nonbridge_strengthening_target_theorem_summary.json"
IN_N554 = (
    GENERATED
    / "n554_current_first_actual_legacy_to_strict_kernel_nonbridge_strengthening_discharge_witness_theorem_summary.json"
)

OUT = (
    GENERATED
    / "f663_current_actual_legacy_to_strict_kernel_bifurcated_frontier_packet_v2_summary.json"
)


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    n261 = load_json(IN_N261)
    n262 = load_json(IN_N262)
    n554 = load_json(IN_N554)

    bridge_target_exported = bool(n261["theorem_result"]["bridge_target_exported"])
    nonbridge_target_exported = bool(n262["theorem_result"]["nonbridge_target_exported"])

    actual_nonbridge_strengthening_discharged = bool(
        n554["theorem_result"]["actual_nonbridge_strengthening_discharged"]
    )

    summary = {
        "packet_id": "F663",
        "status": "F663_EXECUTED_CURRENT_ACTUAL_LEGACY_TO_STRICT_KERNEL_BIFURCATED_FRONTIER_PACKET_V2_NO_FALSE_PASS",
        "as_of": "2026-03-17",
        "frontier_packet_name": "Xi_legacy_strict_frontier_status_packet_v2",
        "bridge_target_exported": bridge_target_exported,
        "bridge_branch_future_only": True,
        "nonbridge_target_exported": nonbridge_target_exported,
        "nonbridge_branch_future_only": bool(not actual_nonbridge_strengthening_discharged),
        "actual_bifurcated_frontier_packet_exported": bool(
            bridge_target_exported and nonbridge_target_exported
        ),
        "branch_selection_justified_on_current_repo_state": False,
        "actual_bridge_discharged": False,
        "actual_nonbridge_strengthening_discharged": actual_nonbridge_strengthening_discharged,
        "kernel_split_safe": True,
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()

