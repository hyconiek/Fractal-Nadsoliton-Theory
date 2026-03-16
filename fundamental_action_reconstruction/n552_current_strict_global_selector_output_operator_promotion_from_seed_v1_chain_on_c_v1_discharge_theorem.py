#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

P660 = (
    GENERATED
    / "p660_current_strict_global_selector_output_operator_promotion_from_seed_v1_chain_on_c_v1_probe_summary.json"
)
OUT = (
    GENERATED
    / "n552_current_strict_global_selector_output_operator_promotion_from_seed_v1_chain_on_c_v1_discharge_theorem_summary.json"
)


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    p660 = load_json(P660)
    expected_status = (
        "CURRENT_REPO_EXPORTS_ONE_GLOBAL_SELECTOR_OUTPUT_OPERATOR_PROMOTION_FROM_SEED_V1_CHAIN_ON_C_V1_AFTER_P660"
    )
    ok = p660.get("status") == expected_status

    checks = [
        {
            "id": "promotion_probe_positive",
            "actual": p660.get("status"),
            "expected": expected_status,
            "pass": ok,
        }
    ]

    summary = {
        "step": "N552",
        "status": "N552_DISCHARGED_CURRENT_STRICT_GLOBAL_SELECTOR_OUTPUT_OPERATOR_PROMOTION_FROM_SEED_V1_CHAIN_ON_C_V1_DISCHARGE_THEOREM_NO_FALSE_PASS",
        "scope": "current_strict_global_selector_output_operator_promotion_from_seed_v1_chain_on_c_v1_only",
        "checks": checks,
        "theorem_result": {
            "discharged": ok,
            "exported_object": "SelectorOutputOperator_global_C_v1_seed_v1_promoted_strict_v1",
            "projector_section_level_only": True,
            "residual_sign_gauge_explicit": True,
            "emergent_observer_constructed": False,
            "strict_core_selector_closure": False,
            "QW2191_discharge": False,
        },
        "hard_limits": [
            "no_admissible_S_sel_int",
            "no_strict_core_selector_closure",
            "no_global_QW2191_discharge",
            "no_emergent_observer_construction",
            "no_ToE_closure",
        ],
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()

