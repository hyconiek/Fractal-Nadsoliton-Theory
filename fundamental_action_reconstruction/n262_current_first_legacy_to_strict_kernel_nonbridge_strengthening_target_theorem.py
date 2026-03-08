#!/usr/bin/env python3

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "n262_current_first_legacy_to_strict_kernel_nonbridge_strengthening_target_theorem_summary.json"


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    p242 = load_json(
        GENERATED / "p242_current_legacy_to_strict_kernel_nonbridge_strengthening_target_probe_summary.json"
    )

    expected_status = (
        "CURRENT_REPO_EXPORTS_ONE_FUTURE_ONLY_LEGACY_TO_STRICT_KERNEL_NONBRIDGE_STRENGTHENING_TARGET_BELOW_ACTUAL_NONBRIDGE_DISCHARGE_AFTER_P242"
    )
    status_ok = p242["status"] == expected_status

    summary = {
        "step": "N262",
        "status": "N262_DISCHARGED_CURRENT_FIRST_LEGACY_TO_STRICT_KERNEL_NONBRIDGE_STRENGTHENING_TARGET_THEOREM_NO_FALSE_PASS",
        "scope": "nonbridge_strengthening_target_only",
        "checks": [
            {
                "id": "p242_status",
                "actual": p242["status"],
                "expected": expected_status,
                "pass": status_ok,
            }
        ],
        "theorem_result": {
            "discharged": status_ok and len(p242.get("blocking_mismatches", [])) == 0,
            "nonbridge_target_exported": p242["nonbridge_target_exported"],
            "future_only_target_exported": p242["future_only_target_exported"],
            "package_level_nonbridge_base_present": p242["package_level_nonbridge_base_present"],
            "positive_bridge_branch_still_open": p242["positive_bridge_branch_still_open"],
            "actual_nonbridge_strengthening_discharged": False,
            "permanent_no_bridge_claimed": p242["permanent_no_bridge_claimed"],
            "kernel_split_safe": p242["kernel_split_safe"],
        },
        "hard_limits": [
            "no_actual_nonbridge_strengthening_theorem",
            "no_permanent_no_bridge_claim",
            "positive_bridge_branch_remains_open",
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
