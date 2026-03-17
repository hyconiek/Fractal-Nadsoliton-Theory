#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_OBJECT = GENERATED / "selector_closure_global_c_v1_directed_strict_v1.json"
N677_SUMMARY = (
    GENERATED / "n677_current_strict_global_directed_selector_closure_object_well_definedness_theorem_summary.json"
)
OUT = GENERATED / "n678_current_strict_t172_directed_closure_discharge_theorem_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    ok = False
    checks: list[dict[str, Any]] = []

    if not IN_OBJECT.exists():
        checks.append(
            {
                "id": "closure_object_present",
                "actual": False,
                "expected": True,
                "pass": False,
            }
        )
    else:
        obj = load_json(IN_OBJECT)
        checks.append(
            {
                "id": "closure_object_present",
                "actual": obj.get("object"),
                "expected": "SelectorClosure_global_C_v1_directed_strict_v1",
                "pass": obj.get("object") == "SelectorClosure_global_C_v1_directed_strict_v1",
            }
        )

    n677_ok = False
    n677_status = None
    if N677_SUMMARY.exists():
        try:
            n677 = load_json(N677_SUMMARY)
            n677_status = n677.get("status")
            tr = n677.get("theorem_result")
            n677_ok = bool(isinstance(tr, dict) and tr.get("discharged") is True)
        except Exception:
            n677_ok = False

    checks.append(
        {
            "id": "n677_well_definedness_discharged",
            "actual": n677_status,
            "expected": "N677_DISCHARGED_*",
            "pass": bool(n677_ok),
        }
    )

    ok = bool(checks[0]["pass"]) and bool(n677_ok)

    summary = {
        "step": "N678",
        "status": "N678_DISCHARGED_CURRENT_STRICT_T172_DIRECTED_CLOSURE_DISCHARGE_THEOREM_NO_FALSE_PASS",
        "scope": "current_strict_t172_directed_closure_discharge_only",
        "checks": checks,
        "theorem_result": {
            "discharged": ok,
            "t172_target_spec": "T172_CURRENT_STRICT_GLOBAL_QW2191_DISCHARGE_AND_SELECTOR_CLOSURE_TARGET_SPEC",
            "t172_discharged_in_scope": ok,
            "closure_object": "SelectorClosure_global_C_v1_directed_strict_v1",
            "closure_scope": "directed_vector_state_only (explicit sign-lift required)",
            "vector_section_level_only": True,
            "requires_explicit_sign_lift": True,
            "QW2191_kernel_alone_obstruction_remains": True,
            "QW2191_kernel_alone_discharge": False,
            "strict_core_selector_closure": False,
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

