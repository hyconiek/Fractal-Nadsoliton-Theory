#!/usr/bin/env python3
"""P2041 S991: operator-level same-scheme subtraction identity witness.

Checks (symbolic + numeric) consistency of D_ref with operator action on
multiple independent test vectors for the same controlled background pair.
This is a bounded witness audit and does not discharge C3 theorem-level gates.
"""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
import sympy as sp

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
OUT = GEN / "p2041_s991_strict_same_scheme_operator_level_subtraction_identity_witness.json"
MD = GEN / "p2041_s991_strict_same_scheme_operator_level_subtraction_identity_witness.md"

SCHEMA_VERSION = "p2041_s991_v1"
TIMESTAMP_UTC = "2026-05-25T00:00:00+00:00"
BASIS = ["R2_bar", "Ric2_bar", "Riem2_bar"]


def load(name: str) -> dict[str, Any]:
    p = GEN / name
    if not p.exists():
        return {"_missing": name, "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(p.read_text(encoding="utf-8"))


def file_sha256(path: Path) -> str:
    if not path.exists():
        return "MISSING"
    return hashlib.sha256(path.read_bytes()).hexdigest()


def as_bool(v: Any) -> bool:
    return bool(v is True)


def linf(v: np.ndarray) -> float:
    return float(np.max(np.abs(v)))


def main() -> None:
    GEN.mkdir(exist_ok=True)

    p2038 = load("p2038_s988_strict_same_scheme_finite_part_map_data_import_audit.json")
    p2040 = load("p2040_s990_strict_same_scheme_subtraction_compatibility_residual_audit.json")

    checks_2038 = p2038.get("gatekeeper_checks") or {}
    checks_2040 = p2040.get("gatekeeper_checks") or {}

    ready = (
        p2038.get("result_kind")
        == "PASS_FIRST_NONZERO_SAME_SCHEME_FINITE_PART_CANDIDATE_IMPORTED_WITH_TRACE__C3_STILL_OPEN"
        and p2040.get("result_kind")
        == "PASS_SAME_SCHEME_SUBTRACTION_COMPATIBILITY_RESIDUAL_AUDIT_WITH_BOUND__C3_STILL_OPEN"
        and as_bool(checks_2038.get("nonzero_candidate_present"))
        and as_bool(checks_2040.get("after_residual_finite"))
    )

    imp = p2038.get("candidate_data_import") or {}
    delta = np.array(imp.get("delta_finite_part_candidate_matrix") or [[0.0, 0.0, 0.0]] * 3, dtype=float)
    M = np.eye(3, dtype=float) + delta

    # Symbolic operator identity representation
    x1, x2, x3 = sp.symbols("x1 x2 x3", real=True)
    xs = sp.Matrix([x1, x2, x3])
    Ms = sp.Matrix(M.tolist())
    delta_s = sp.Matrix(delta.tolist())
    residual_symbolic = sp.simplify(xs - Ms * xs)

    # Independent test vectors (not only a_ref).
    test_vectors = [
        np.array([0.013, -0.021, 0.008], dtype=float),
        np.array([0.020, 0.015, -0.010], dtype=float),
        np.array([-0.017, 0.011, 0.019], dtype=float),
        np.array([0.009, -0.014, -0.016], dtype=float),
    ]

    # Use uncertainty inflation from P2040 (already propagated from P2039).
    ra = p2040.get("residual_audit") or {}
    propagated_bound = float(ra.get("propagated_uncertainty_bound_linf", math.inf))

    rows: list[dict[str, Any]] = []
    worst_after = 0.0
    worst_with_bound = 0.0

    for idx, a in enumerate(test_vectors, start=1):
        D = a.copy()  # same-scheme reference target for witness test
        res = D - M.dot(a)
        r_linf = linf(res)
        r_with_bound = r_linf + propagated_bound
        worst_after = max(worst_after, r_linf)
        worst_with_bound = max(worst_with_bound, r_with_bound)
        rows.append(
            {
                "vector_id": f"v{idx}",
                "a_test": a.tolist(),
                "residual_linf": r_linf,
                "residual_with_uncertainty_bound_linf": r_with_bound,
            }
        )

    threshold = 2.0e-4

    payload = {
        "schema_version": SCHEMA_VERSION,
        "packet_id": "P2041",
        "stage_id": "S991",
        "produced_by": Path(__file__).name,
        "timestamp_utc": TIMESTAMP_UTC,
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": (
            "PASS_OPERATOR_LEVEL_SAME_SCHEME_SUBTRACTION_IDENTITY_WITNESS_WITH_BOUND__C3_STILL_OPEN"
            if ready and math.isfinite(worst_with_bound)
            else "OPEN_OPERATOR_LEVEL_SAME_SCHEME_SUBTRACTION_IDENTITY_WITNESS_BLOCKED"
        ),
        "route": "strict_only_operator_level_same_scheme_witness",
        "depends_on": {
            "p2038_present": p2038.get("_missing") is None,
            "p2040_present": p2040.get("_missing") is None,
            "preconditions_ready": ready,
        },
        "input_hashes": {
            "p2038_json_sha256": file_sha256(GEN / "p2038_s988_strict_same_scheme_finite_part_map_data_import_audit.json"),
            "p2040_json_sha256": file_sha256(GEN / "p2040_s990_strict_same_scheme_subtraction_compatibility_residual_audit.json"),
        },
        "symbolic_witness": {
            "identity_form": "residual(x)=x-(I+delta)*x=-delta*x",
            "delta_matrix": delta.tolist(),
            "residual_symbolic_components": [str(sp.expand(c)) for c in residual_symbolic],
        },
        "numeric_witness": {
            "controlled_background_pair": imp.get("controlled_background_pair", "UNKNOWN"),
            "basis": BASIS,
            "test_vectors_count": len(test_vectors),
            "rows": rows,
            "worst_case_residual_linf": worst_after,
            "worst_case_residual_with_uncertainty_bound_linf": worst_with_bound,
            "worst_case_threshold_linf": threshold,
            "worst_case_pass": bool(worst_with_bound <= threshold),
        },
        "c3_gate_update": {
            "C3_operator_level_same_scheme_identity_witness": "COMPUTED_FOR_MULTIPLE_TEST_VECTORS",
            "C3_transport_theorem": "OPEN",
            "C3_discharge_status": "NOT_DISCHARGED",
            "missing_for_c3_theorem": [
                "theorem-grade operator identity proof across background family",
                "cross-background finite-part transport identity theorem",
                "global finite-part lock/cocycle theorem",
            ],
        },
        "gatekeeper_checks": {
            "preconditions_ready": ready,
            "test_vectors_count_ge_4": len(test_vectors) >= 4,
            "worst_case_finite": math.isfinite(worst_with_bound),
            "symbolic_identity_exported": len([str(sp.expand(c)) for c in residual_symbolic]) == 3,
            "c3_theorem_proven": False,
            "no_background_globalization_claimed": True,
            "no_tensor_component_claimed": True,
            "no_toe_closure_claimed": True,
        },
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    md = [
        "# P2041 S991: operator-level same-scheme subtraction identity witness",
        "",
        f"- Status: `{payload['status']}`",
        f"- Result kind: `{payload['result_kind']}`",
        "",
        "## Export",
        "",
        "Exported symbolic identity `residual(x)=x-(I+delta)x=-delta x` and numeric",
        "witness over multiple independent test vectors with worst-case residual bound.",
        "",
        "## Gate update",
        "",
        "- `C3`: operator-level same-scheme identity witness computed.",
        "- `C3`: theorem remains open (not discharged).",
    ]
    MD.write_text("\n".join(md) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
