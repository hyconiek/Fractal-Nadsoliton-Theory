#!/usr/bin/env python3

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "f153_first_actual_legacy_to_strict_kernel_bifurcated_frontier_packet_summary.json"


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    n261 = load_json(
        GENERATED / "n261_current_first_legacy_to_strict_kernel_bridge_target_theorem_summary.json"
    )
    n262 = load_json(
        GENERATED / "n262_current_first_legacy_to_strict_kernel_nonbridge_strengthening_target_theorem_summary.json"
    )

    summary = {
        "packet_id": "F153",
        "status": "F153_EXECUTED_FIRST_ACTUAL_LEGACY_TO_STRICT_KERNEL_BIFURCATED_FRONTIER_PACKET_NO_FALSE_PASS",
        "as_of": "2026-03-08",
        "frontier_packet_name": "Xi_legacy_strict_frontier_bifurcation_packet_v1",
        "bridge_target_exported": n261["theorem_result"]["bridge_target_exported"],
        "bridge_branch_future_only": n261["theorem_result"]["future_only_target_exported"],
        "nonbridge_target_exported": n262["theorem_result"]["nonbridge_target_exported"],
        "nonbridge_branch_future_only": n262["theorem_result"]["future_only_target_exported"],
        "frontier_packet": {
            "positive_branch": "B_legacy_strict_bridge_target_v1",
            "negative_branch": "NB_legacy_strict_strengthening_target_v1",
        },
        "actual_bifurcated_frontier_packet_exported": True,
        "branch_selection_justified_on_current_repo_state": False,
        "actual_bridge_discharged": False,
        "actual_nonbridge_strengthening_discharged": False,
        "legacy_physical_role_transfer_claimed": False,
        "permanent_no_bridge_claimed": False,
        "kernel_split_safe": (
            n261["theorem_result"]["kernel_split_safe"]
            and n262["theorem_result"]["kernel_split_safe"]
        ),
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
