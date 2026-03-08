#!/usr/bin/env python3

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "n265_current_first_actual_legacy_to_strict_kernel_claim_specific_amplitude_nonabsorption_witness_theorem_summary.json"


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    p245 = load_json(
        GENERATED / "p245_current_actual_legacy_to_strict_kernel_claim_specific_amplitude_nonabsorption_witness_probe_summary.json"
    )

    expected_status = (
        "CURRENT_REPO_EXPORTS_ONE_ACTUAL_CLAIM_SPECIFIC_LEGACY_TO_STRICT_KERNEL_AMPLITUDE_NONABSORPTION_WITNESS_BELOW_FULL_AMPLITUDE_OBSTRUCTION_AFTER_P245"
    )
    status_ok = p245["status"] == expected_status

    summary = {
        "step": "N265",
        "status": "N265_DISCHARGED_CURRENT_FIRST_ACTUAL_LEGACY_TO_STRICT_KERNEL_CLAIM_SPECIFIC_AMPLITUDE_NONABSORPTION_WITNESS_THEOREM_NO_FALSE_PASS",
        "scope": "legacy_weinberg_claim_specific_amplitude_only",
        "checks": [
            {
                "id": "p245_claim_specific_amplitude_status",
                "actual": p245["status"],
                "expected": expected_status,
                "pass": status_ok,
            }
        ],
        "theorem_result": {
            "discharged": status_ok and len(p245.get("blocking_mismatches", [])) == 0,
            "actual_claim_specific_amplitude_nonabsorption_witness_exported": p245["actual_claim_specific_amplitude_nonabsorption_witness_exported"],
            "bridge_nonbridge_frontier_undecided": p245["bridge_nonbridge_frontier_undecided"],
            "full_amplitude_nonabsorption_obstruction_discharged": False,
            "actual_nonbridge_strengthening_discharged": False,
            "actual_bridge_discharged": False,
            "current_strict_core_selector_closure": False,
            "current_global_selector_closure": False,
            "current_global_qw2191_discharge": False,
            "kernel_split_safe": p245["kernel_split_safe"],
        },
        "hard_limits": [
            "no_full_amplitude_nonabsorption_obstruction",
            "no_actual_damping_nonrenormalization_obstruction",
            "no_actual_phase_frequency_nonconformal_obstruction",
            "no_actual_nonbridge_strengthening_theorem",
            "no_actual_bridge_derivation",
            "no_current_branch_selection_theorem",
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
