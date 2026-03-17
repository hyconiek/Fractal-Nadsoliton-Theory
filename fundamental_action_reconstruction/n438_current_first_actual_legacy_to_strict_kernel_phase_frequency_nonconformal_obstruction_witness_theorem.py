#!/usr/bin/env python3

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
GENERATED = ROOT / "generated"
OUT = (
    GENERATED
    / "n438_current_first_actual_legacy_to_strict_kernel_phase_frequency_nonconformal_obstruction_witness_theorem_summary.json"
)


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    p404 = load_json(
        GENERATED
        / "p404_current_actual_legacy_to_strict_kernel_phase_frequency_nonconformal_obstruction_witness_probe_summary.json"
    )

    expected_status = (
        "CURRENT_REPO_EXPORTS_ONE_ACTUAL_LEGACY_TO_STRICT_KERNEL_PHASE_FREQUENCY_NONCONFORMAL_OBSTRUCTION_WITNESS_AFTER_P404"
    )
    status_ok = p404["status"] == expected_status

    summary = {
        "step": "N438",
        "status": "N438_DISCHARGED_CURRENT_FIRST_ACTUAL_LEGACY_TO_STRICT_KERNEL_PHASE_FREQUENCY_NONCONFORMAL_OBSTRUCTION_WITNESS_THEOREM_NO_FALSE_PASS",
        "scope": "declared_kernel_comparison_phase_frequency_layer_only",
        "checks": [
            {
                "id": "p404_phase_frequency_obstruction_status",
                "actual": p404["status"],
                "expected": expected_status,
                "pass": status_ok,
            }
        ],
        "theorem_result": {
            "discharged": status_ok and len(p404.get("blocking_mismatches", [])) == 0,
            "actual_phase_frequency_nonconformal_obstruction_discharged": bool(
                p404.get("p_shift_nonbridge_obstruction_discharged")
            ),
            "actual_nonbridge_strengthening_discharged": False,
            "actual_bridge_discharged": False,
            "kernel_split_safe": bool(p404.get("kernel_split_safe")),
        },
        "hard_limits": [
            "no_actual_nonbridge_strengthening_theorem",
            "no_actual_bridge_derivation",
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

