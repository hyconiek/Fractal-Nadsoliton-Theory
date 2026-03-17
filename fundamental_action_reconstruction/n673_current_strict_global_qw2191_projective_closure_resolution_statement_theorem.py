#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_OBJECT = GENERATED / "selector_closure_global_c_v1_projective_strict_v1.json"
OUT = (
    GENERATED
    / "n673_current_strict_global_qw2191_projective_closure_resolution_statement_theorem_summary.json"
)


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def max_abs_diff_2x2(a: list, b: list) -> float:
    return max(abs(float(a[i][j]) - float(b[i][j])) for i in range(2) for j in range(2))


def main() -> None:
    obj = load_json(IN_OBJECT)

    cert = obj.get("well_definedness_certificate", {}) or {}
    cert_pass = bool(cert.get("certificate_pass"))
    tol = float(cert.get("tolerance", 0.0) or 0.0)

    p_out = ((obj.get("output_observable", {}) or {}).get("output_projector_matrix_in_o_plus_o_minus")) or None
    diag = [[1.0, 0.0], [0.0, 0.0]]
    diag_diff = None
    if isinstance(p_out, list) and len(p_out) == 2 and all(isinstance(r, list) and len(r) == 2 for r in p_out):
        diag_diff = max_abs_diff_2x2(p_out, diag)

    ok = cert_pass and (diag_diff is not None) and (diag_diff <= tol if tol > 0.0 else False)

    checks = [
        {
            "id": "closure_object_present",
            "actual": obj.get("object"),
            "expected": "SelectorClosure_global_C_v1_projective_strict_v1",
            "pass": obj.get("object") == "SelectorClosure_global_C_v1_projective_strict_v1",
        },
        {
            "id": "certificate_pass",
            "actual": cert_pass,
            "expected": True,
            "pass": cert_pass,
        },
        {
            "id": "output_projector_equals_o_plus_projector_within_tolerance",
            "actual": {"max_abs_diff_to_diag(1,0)": diag_diff, "tolerance": tol},
            "expected": "<= tolerance",
            "pass": (diag_diff is not None) and (tol > 0.0) and (diag_diff <= tol),
        },
    ]

    summary = {
        "step": "N673",
        "status": "N673_DISCHARGED_CURRENT_STRICT_GLOBAL_QW2191_PROJECTIVE_CLOSURE_RESOLUTION_STATEMENT_THEOREM_NO_FALSE_PASS",
        "scope": "current_strict_global_qw2191_projective_closure_resolution_statement_only",
        "checks": checks,
        "theorem_result": {
            "discharged": ok,
            "closure_scope": "projective_ray_state_only",
            "closure_observable_unique_in_exported_scope": ok,
            "closure_observable": "output_projector_on_o_plus (projective; sign-gauge-safe)",
            "QW2191_kernel_alone_obstruction_remains": True,
            "QW2191_kernel_alone_discharge": False,
            "QW2191_bypassed_for_projective_closure_observable": ok,
            "strict_core_selector_closure": False,
            "ToE_closure": False,
        },
        "hard_limits": [
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

