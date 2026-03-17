#!/usr/bin/env python3

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

AS_OF = "2026-03-17"

IN_P699 = (
    GENERATED
    / "p699_current_strict_projective_emergent_observer_chain_from_global_projective_selector_closure_output_probe.json"
)
OUT = GENERATED / "n699_current_first_strict_projective_emergent_observer_chain_computability_theorem_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def is_finite_number(x: Any) -> bool:
    return isinstance(x, (int, float)) and math.isfinite(float(x))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    if not IN_P699.exists():
        summary = {
            "step": "N699",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": [str(IN_P699.relative_to(REPO))],
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT)
        return

    p699 = load_json(IN_P699)
    results = p699.get("results") or {}

    required = ["z_commit", "z_residual", "w_commit", "w_residual", "x_commit", "x_residual", "u_commit", "u_residual", "c_closure"]
    missing_fields = [k for k in required if k not in results]

    status_ok = (
        p699.get("no_false_pass") is True
        and p699.get("stage") == "P699"
        and p699.get("status")
        == "PASS_PROJECTIVE_EMERGENT_OBSERVER_CHAIN_COMPUTABLE_FROM_PROJECTIVE_SELECTOR_CLOSURE_OUTPUT"
        and not missing_fields
    )

    numeric_ok = all(is_finite_number(results.get(k)) for k in required)

    tol = results.get("tolerance_abs")
    if not is_finite_number(tol):
        tol = 1e-12
    tol = float(tol)

    residuals_ok = bool(
        (float(results.get("z_commit")) > 0.0)
        and (abs(float(results.get("z_residual"))) <= tol)
        and (float(results.get("w_commit")) > 0.0)
        and (abs(float(results.get("w_residual"))) <= tol)
        and (float(results.get("x_commit")) > 0.0)
        and (abs(float(results.get("x_residual"))) <= tol)
        and (float(results.get("u_commit")) > 0.0)
        and (abs(float(results.get("u_residual"))) <= tol)
        and (float(results.get("c_closure")) > 0.0)
    )

    discharged = bool(status_ok and numeric_ok and residuals_ok)

    checks = [
        {
            "id": "p699_status_is_pass",
            "actual": p699.get("status"),
            "expected": "PASS_PROJECTIVE_EMERGENT_OBSERVER_CHAIN_COMPUTABLE_FROM_PROJECTIVE_SELECTOR_CLOSURE_OUTPUT",
            "pass": bool(
                p699.get("status")
                == "PASS_PROJECTIVE_EMERGENT_OBSERVER_CHAIN_COMPUTABLE_FROM_PROJECTIVE_SELECTOR_CLOSURE_OUTPUT"
            ),
        },
        {
            "id": "p699_no_false_pass_true",
            "actual": p699.get("no_false_pass"),
            "expected": True,
            "pass": bool(p699.get("no_false_pass") is True),
        },
        {
            "id": "p699_results_required_fields_present",
            "actual": missing_fields,
            "expected": [],
            "pass": bool(len(missing_fields) == 0),
        },
        {
            "id": "p699_results_numeric",
            "actual": numeric_ok,
            "expected": True,
            "pass": bool(numeric_ok),
        },
        {
            "id": "p699_residuals_within_tolerance",
            "actual": {
                "tolerance_abs": tol,
                "z_residual": results.get("z_residual"),
                "w_residual": results.get("w_residual"),
                "x_residual": results.get("x_residual"),
                "u_residual": results.get("u_residual"),
            },
            "expected": "all_residuals_abs<=tolerance_abs",
            "pass": bool(residuals_ok),
        },
    ]

    status = "N699_DERIVABLE_CURRENT_FIRST_STRICT_PROJECTIVE_EMERGENT_OBSERVER_CHAIN_COMPUTABILITY_THEOREM_NO_FALSE_PASS"
    if not discharged:
        status = "N699_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_PROJECTIVE_EMERGENT_OBSERVER_CHAIN_STATE"

    theorem_result = {
        "discharged": bool(discharged),
        "projective_emergent_observer_chain_computable_from_projective_selector_closure_output_projector": bool(discharged),
        "z_commit": results.get("z_commit"),
        "z_residual": results.get("z_residual"),
        "w_commit": results.get("w_commit"),
        "w_residual": results.get("w_residual"),
        "x_commit": results.get("x_commit"),
        "x_residual": results.get("x_residual"),
        "u_commit": results.get("u_commit"),
        "u_residual": results.get("u_residual"),
        "c_closure_candidate": results.get("c_closure"),
        "kernel_alone_qw2191_discharge": False,
        "ToE_closure": False,
        "actual_emergent_observer_closure": False,
        "evidence": {"P699": str(IN_P699.relative_to(REPO))},
    }

    summary = {
        "step": "N699",
        "status": status,
        "as_of": AS_OF,
        "scope": "strict_projective_emergent_observer_chain_computability_only",
        "checks": checks,
        "theorem_result": theorem_result,
        "hard_limits": [
            "no_kernel_alone_global_QW2191_discharge",
            "no_directed_sign_sensitive_physical_orientation_claim",
            "no_ToE_closure",
            "no_actual_emergent_observer_closure",
        ],
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()

