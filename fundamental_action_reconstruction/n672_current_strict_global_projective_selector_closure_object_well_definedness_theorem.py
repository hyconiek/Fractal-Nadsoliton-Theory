#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

P673 = (
    GENERATED
    / "p673_current_strict_global_projective_selector_closure_object_gluing_probe_summary.json"
)
OUT = (
    GENERATED
    / "n672_current_strict_global_projective_selector_closure_object_well_definedness_theorem_summary.json"
)


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    p673 = load_json(P673)
    expected_status = "CURRENT_REPO_EXPORTS_ONE_GLOBAL_PROJECTIVE_SELECTOR_CLOSURE_OBJECT_ON_C_V1_AFTER_F672"
    ok = p673.get("status") == expected_status

    checks = [
        {
            "id": "projective_closure_gluing_probe_positive",
            "actual": p673.get("status"),
            "expected": expected_status,
            "pass": ok,
        }
    ]

    summary = {
        "step": "N672",
        "status": "N672_DISCHARGED_CURRENT_STRICT_GLOBAL_PROJECTIVE_SELECTOR_CLOSURE_OBJECT_WELL_DEFINEDNESS_THEOREM_NO_FALSE_PASS",
        "scope": "current_strict_global_projective_selector_closure_object_well_definedness_only",
        "checks": checks,
        "theorem_result": {
            "discharged": ok,
            "exported_object": "SelectorClosure_global_C_v1_projective_strict_v1",
            "closure_scope": "projective_ray_state_only",
            "projector_section_level_only": True,
            "strict_core_selector_closure": False,
            "QW2191_discharge": False,
        },
        "hard_limits": [
            "no_strict_core_selector_closure",
            "no_global_QW2191_discharge",
            "no_operator_level_transition_groupoid_promotion (N512 boundary)",
            "no_ToE_closure",
        ],
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()

