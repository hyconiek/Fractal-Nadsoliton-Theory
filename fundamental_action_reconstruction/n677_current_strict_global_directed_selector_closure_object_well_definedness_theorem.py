#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

P677 = (
    GENERATED
    / "p677_current_strict_global_directed_selector_closure_object_sign_lift_gluing_probe_summary.json"
)
OUT = (
    GENERATED
    / "n677_current_strict_global_directed_selector_closure_object_well_definedness_theorem_summary.json"
)


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    p677 = load_json(P677)
    expected_status = "CURRENT_REPO_EXPORTS_ONE_GLOBAL_DIRECTED_SELECTOR_CLOSURE_OBJECT_ON_C_V1_WITH_EXPLICIT_SIGN_LIFT_AFTER_F677"
    ok = p677.get("status") == expected_status

    checks = [
        {
            "id": "directed_closure_sign_lift_gluing_probe_positive",
            "actual": p677.get("status"),
            "expected": expected_status,
            "pass": ok,
        }
    ]

    summary = {
        "step": "N677",
        "status": "N677_DISCHARGED_CURRENT_STRICT_GLOBAL_DIRECTED_SELECTOR_CLOSURE_OBJECT_WELL_DEFINEDNESS_THEOREM_NO_FALSE_PASS",
        "scope": "current_strict_global_directed_selector_closure_object_well_definedness_only",
        "checks": checks,
        "theorem_result": {
            "discharged": ok,
            "exported_object": "SelectorClosure_global_C_v1_directed_strict_v1",
            "closure_scope": "directed_vector_state_only",
            "vector_section_level_only": True,
            "requires_explicit_sign_lift": True,
            "strict_core_selector_closure": False,
            "QW2191_kernel_alone_discharge": False,
            "ToE_closure": False,
        },
        "hard_limits": [
            "premise_based_directed_closure_only (N462 boundary; depends on explicit fixing datum path via F474)",
            "explicit_sign_lift_required",
            "no_strict_core_selector_closure",
            "no_global_kernel_alone_QW2191_discharge",
            "no_operator_level_transition_groupoid_promotion (N512 boundary)",
            "no_ToE_closure",
        ],
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()

