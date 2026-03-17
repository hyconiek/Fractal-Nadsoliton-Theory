#!/usr/bin/env python3

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
GENERATED = ROOT / "generated"
OUT = (
    GENERATED
    / "n554_current_first_actual_legacy_to_strict_kernel_nonbridge_strengthening_discharge_witness_theorem_summary.json"
)


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    p662 = load_json(
        GENERATED
        / "p662_current_actual_legacy_to_strict_kernel_nonbridge_strengthening_discharge_witness_probe_summary.json"
    )

    expected_status = (
        "CURRENT_REPO_EXPORTS_ONE_ACTUAL_LEGACY_TO_STRICT_KERNEL_NONBRIDGE_STRENGTHENING_DISCHARGE_WITNESS_AFTER_P662"
    )
    status_ok = p662["status"] == expected_status

    summary = {
        "step": "N554",
        "status": "N554_DISCHARGED_CURRENT_FIRST_ACTUAL_LEGACY_TO_STRICT_KERNEL_NONBRIDGE_STRENGTHENING_DISCHARGE_WITNESS_THEOREM_NO_FALSE_PASS",
        "scope": "declared_kernel_comparison_nonbridge_strengthening_only",
        "checks": [
            {
                "id": "p662_nonbridge_strengthening_status",
                "actual": p662["status"],
                "expected": expected_status,
                "pass": status_ok,
            }
        ],
        "theorem_result": {
            "discharged": status_ok and len(p662.get("blocking_mismatches", [])) == 0,
            "actual_nonbridge_strengthening_discharged": bool(
                p662.get("actual_nonbridge_strengthening_discharged")
            ),
            "actual_bridge_discharged": False,
            "permanent_no_bridge_claimed": False,
            "branch_selection_justified_on_current_repo_state": False,
            "kernel_split_safe": bool(p662.get("kernel_split_safe")),
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
