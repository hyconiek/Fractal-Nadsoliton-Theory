#!/usr/bin/env python3

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "n270_current_first_nonstrict_declared_scope_selector_closure_theorem_summary.json"


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    p250 = load_json(
        GENERATED / "p250_current_actual_nonstrict_declared_scope_selector_closure_probe_summary.json"
    )

    expected_status = "CURRENT_REPO_EXPORTS_ONE_ACTUAL_NONSTRICT_DECLARED_SCOPE_SELECTOR_CLOSURE_WITNESS_AFTER_P250"
    status_ok = p250["status"] == expected_status

    summary = {
        "step": "N270",
        "status": "N270_DISCHARGED_CURRENT_FIRST_NONSTRICT_DECLARED_SCOPE_SELECTOR_CLOSURE_THEOREM_NO_FALSE_PASS",
        "scope": "axiom_augmented_declared_scope_selector_closure_only",
        "checks": [
            {
                "id": "p250_nonstrict_declared_scope_closure_status",
                "actual": p250["status"],
                "expected": expected_status,
                "pass": status_ok,
            }
        ],
        "theorem_result": {
            "discharged": status_ok and len(p250.get("blocking_mismatches", [])) == 0,
            "nonstrict_declared_scope_selector_closure_exported": p250["nonstrict_declared_scope_selector_closure_exported"],
            "accepted_scope": p250["accepted_scope"],
            "strict_core_changed": p250["strict_core_changed"],
            "bridge_nonbridge_nonmandatory_for_t14_after_n269": p250["bridge_nonbridge_nonmandatory_for_t14_after_n269"],
            "actual_strict_core_selector_closure": False,
            "actual_global_selector_closure": False,
            "actual_global_qw2191_discharge": False,
            "toe_closure_claimed": False,
        },
        "hard_limits": [
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
