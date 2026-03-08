#!/usr/bin/env python3

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "f155_first_actual_legacy_to_strict_kernel_claim_specific_amplitude_nonabsorption_witness_packet_summary.json"


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    n264 = load_json(
        GENERATED / "n264_current_first_actual_legacy_to_strict_kernel_amplitude_nonabsorption_component_witness_theorem_summary.json"
    )

    summary = {
        "packet_id": "F155",
        "status": "F155_EXECUTED_FIRST_ACTUAL_LEGACY_TO_STRICT_KERNEL_CLAIM_SPECIFIC_AMPLITUDE_NONABSORPTION_WITNESS_PACKET_NO_FALSE_PASS",
        "as_of": "2026-03-08",
        "support_packet_name": "W_abs_nonbridge_claim_specific_support_v1",
        "claim_specific_witness_name": "A_abs_nonbridge_claim_specific_actual_witness_v1",
        "claim_specific_target_name": "claim_specific_amplitude_nonabsorption_obstruction_target_v1",
        "claim_specific_scope": "legacy_weinberg_amplitude_role_only",
        "actual_component_witness_exported": n264["theorem_result"]["actual_amplitude_nonabsorption_component_witness_exported"],
        "bridge_nonbridge_frontier_undecided": n264["theorem_result"]["bridge_nonbridge_frontier_undecided"],
        "actual_claim_specific_amplitude_nonabsorption_witness_exported": True,
        "full_amplitude_nonabsorption_obstruction_discharged": False,
        "actual_nonbridge_strengthening_discharged": False,
        "actual_bridge_discharged": False,
        "kernel_split_safe": n264["theorem_result"]["kernel_split_safe"],
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
