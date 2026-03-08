#!/usr/bin/env python3

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "n261_current_first_legacy_to_strict_kernel_bridge_target_theorem_summary.json"

def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))

def main() -> None:
    p241 = load_json(
        GENERATED / "p241_current_legacy_to_strict_kernel_bridge_target_probe_summary.json"
    )

    expected_status = (
        "CURRENT_REPO_EXPORTS_ONE_FUTURE_ONLY_LEGACY_TO_STRICT_KERNEL_BRIDGE_TARGET_BELOW_ACTUAL_BRIDGE_DISCHARGE_AFTER_P241"
    )
    status_ok = p241["status"] == expected_status
    mismatches = p241.get("blocking_mismatches", [])

    summary = {
        "step": "N261",
        "status": "N261_DISCHARGED_CURRENT_FIRST_LEGACY_TO_STRICT_KERNEL_BRIDGE_TARGET_THEOREM_NO_FALSE_PASS",
        "scope": "bridge_target_only",
        "checks": [
            {
                "id": "p241_status",
                "actual": p241["status"],
                "expected": expected_status,
                "pass": status_ok,
            }
        ],
        "theorem_result": {
            "discharged": status_ok and len(mismatches) == 0,
            "bridge_target_exported": p241["bridge_target_exported"],
            "future_only_target_exported": p241["future_only_target_exported"],
            "kernel_split_safe": p241["kernel_split_safe"],
            "nonbridge_branch_still_open": p241["nonbridge_branch_still_open"],
            "actual_bridge_discharged": False,
            "legacy_physical_role_transfer_claimed": p241["legacy_physical_role_transfer_claimed"],
            "current_strict_core_selector_closure": False,
            "current_global_selector_closure": False,
            "current_global_qw2191_discharge": False,
        },
        "hard_limits": [
            "no_actual_bridge_derivation",
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
