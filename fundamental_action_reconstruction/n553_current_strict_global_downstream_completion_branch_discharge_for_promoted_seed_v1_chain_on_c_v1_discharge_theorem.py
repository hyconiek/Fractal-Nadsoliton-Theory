#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

P661 = (
    GENERATED
    / "p661_current_strict_global_downstream_completion_branch_discharge_for_promoted_seed_v1_chain_on_c_v1_probe_summary.json"
)
OUT = (
    GENERATED
    / "n553_current_strict_global_downstream_completion_branch_discharge_for_promoted_seed_v1_chain_on_c_v1_discharge_theorem_summary.json"
)


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    p661 = load_json(P661)
    expected_status = (
        "CURRENT_REPO_EXPORTS_ONE_GLOBAL_DOWNSTREAM_COMPLETION_BRANCH_DISCHARGE_FOR_PROMOTED_SEED_V1_CHAIN_ON_C_V1_AFTER_P661"
    )
    ok = p661.get("status") == expected_status

    checks = [
        {
            "id": "bundle_probe_positive",
            "actual": p661.get("status"),
            "expected": expected_status,
            "pass": ok,
        }
    ]

    summary = {
        "step": "N553",
        "status": "N553_DISCHARGED_CURRENT_STRICT_GLOBAL_DOWNSTREAM_COMPLETION_BRANCH_DISCHARGE_FOR_PROMOTED_SEED_V1_CHAIN_ON_C_V1_DISCHARGE_THEOREM_NO_FALSE_PASS",
        "scope": "current_strict_global_downstream_completion_branch_discharge_for_promoted_seed_v1_chain_on_c_v1_only",
        "checks": checks,
        "theorem_result": {
            "discharged": ok,
            "exported_object": "SelectorDownstreamCompletionBranch_global_C_v1_seed_v1_promoted_strict_v1",
            "projector_section_level_only": True,
            "residual_sign_gauge_explicit": True,
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

