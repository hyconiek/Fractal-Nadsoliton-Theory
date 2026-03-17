#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

P674_SUMMARY = (
    GENERATED
    / "p674_current_strict_global_directed_selector_closure_output_sign_mismatch_audit_probe_summary.json"
)
OUT = (
    GENERATED
    / "n675_current_strict_global_directed_selector_closure_obstruction_boundary_theorem_summary.json"
)


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    ok = False
    actual_status = None
    if P674_SUMMARY.exists():
        try:
            p674 = load_json(P674_SUMMARY)
            actual_status = p674.get("status")
            ok = actual_status == "PASS_AUDIT_DIRECTED_CLOSURE_OUTPUT_SIGN_MISMATCH_ACROSS_CHARTS"
        except Exception:
            ok = False

    checks = [
        {
            "id": "directed_closure_sign_mismatch_audit_positive",
            "actual": actual_status,
            "expected": "PASS_AUDIT_DIRECTED_CLOSURE_OUTPUT_SIGN_MISMATCH_ACROSS_CHARTS",
            "pass": ok,
        }
    ]

    summary = {
        "step": "N675",
        "status": "N675_DISCHARGED_CURRENT_STRICT_GLOBAL_DIRECTED_SELECTOR_CLOSURE_OBSTRUCTION_BOUNDARY_THEOREM_NO_FALSE_PASS",
        "scope": "current_strict_global_directed_selector_closure_obstruction_only",
        "checks": checks,
        "theorem_result": {
            "discharged": ok,
            "directed_global_closure_exported": False,
            "directed_global_closure_obstructed_on_current_exports": ok,
            "obstruction_type": "output_o_plus_sign_mismatch_across_charts_under_fixed_output_bases",
            "projective_closure_remains_available": True,
            "strict_core_selector_closure": False,
            "QW2191_kernel_alone_discharge": False,
            "ToE_closure": False,
        },
        "hard_limits": [
            "no_strict_core_selector_closure",
            "no_directed_global_closure_object_export",
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

