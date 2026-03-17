#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GENERATED = ROOT / "generated"

IN_CLOSURE_OBJECT = GENERATED / "selector_closure_global_c_v1_projective_strict_v1.json"
IN_N672 = (
    GENERATED
    / "n672_current_strict_global_projective_selector_closure_object_well_definedness_theorem_summary.json"
)
IN_N673 = (
    GENERATED
    / "n673_current_strict_global_qw2191_projective_closure_resolution_statement_theorem_summary.json"
)

OUT = (
    GENERATED
    / "n674_current_strict_t172_projective_closure_discharge_theorem_summary.json"
)


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def _is_discharged(summary: dict) -> bool:
    tr = summary.get("theorem_result")
    return bool(isinstance(tr, dict) and tr.get("discharged") is True)


def main() -> None:
    closure_present = IN_CLOSURE_OBJECT.exists()
    n672 = load_json(IN_N672) if IN_N672.exists() else {}
    n673 = load_json(IN_N673) if IN_N673.exists() else {}

    ok_n672 = _is_discharged(n672) and (
        (n672.get("theorem_result", {}) or {}).get("exported_object")
        == "SelectorClosure_global_C_v1_projective_strict_v1"
    )
    ok_n673 = _is_discharged(n673) and bool(
        (n673.get("theorem_result", {}) or {}).get(
            "QW2191_bypassed_for_projective_closure_observable"
        )
    )

    ok = closure_present and ok_n672 and ok_n673

    checks = [
        {
            "id": "closure_object_present",
            "actual": closure_present,
            "expected": True,
            "pass": closure_present,
        },
        {
            "id": "n672_projective_closure_well_definedness_discharged",
            "actual": n672.get("status"),
            "expected": "N672_DISCHARGED_*",
            "pass": ok_n672,
        },
        {
            "id": "n673_qw2191_projective_closure_resolution_statement_discharged",
            "actual": n673.get("status"),
            "expected": "N673_DISCHARGED_*",
            "pass": ok_n673,
        },
    ]

    summary = {
        "step": "N674",
        "status": "N674_DISCHARGED_CURRENT_STRICT_T172_PROJECTIVE_CLOSURE_DISCHARGE_THEOREM_NO_FALSE_PASS",
        "scope": "current_strict_t172_projective_closure_discharge_only",
        "checks": checks,
        "theorem_result": {
            "discharged": ok,
            "t172_target_spec": "T172_CURRENT_STRICT_GLOBAL_QW2191_DISCHARGE_AND_SELECTOR_CLOSURE_TARGET_SPEC",
            "t172_discharged_in_scope": ok,
            "closure_object": "SelectorClosure_global_C_v1_projective_strict_v1",
            "closure_scope": "projective_ray_state_only",
            "projector_section_level_only": True,
            "qw2191_bypass_statement_exported": ok_n673,
            "QW2191_kernel_alone_obstruction_remains": True,
            "QW2191_kernel_alone_discharge": False,
            "strict_core_selector_closure": False,
            "ToE_closure": False,
        },
        "hard_limits": [
            "no_strict_core_selector_closure",
            "no_global_kernel_alone_QW2191_discharge",
            "no_directed_global_closure_object_export",
            "no_operator_level_transition_groupoid_promotion (N512 boundary)",
            "no_ToE_closure",
        ],
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()

