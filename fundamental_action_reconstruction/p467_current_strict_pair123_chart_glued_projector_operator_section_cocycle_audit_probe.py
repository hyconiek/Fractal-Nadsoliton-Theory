#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_A1 = GENERATED / "a_1_pair1_orientation_projector_operator_strict_core_v1.json"
IN_O12 = GENERATED / "o12_pair1_pair2_selector_chart_transport_operator_strict_derived_from_sigma_int_alpha12_v1.json"
IN_A2 = GENERATED / "a_2_pair2_orientation_projector_operator_strict_core_v1.json"

IN_A3 = GENERATED / "a_3_pair3_orientation_projector_operator_strict_core_v1.json"
IN_O23 = GENERATED / "o23_pair2_pair3_selector_chart_transport_operator_axis_only_alpha23_mod_pi_strict_core_v1.json"
IN_O13 = GENERATED / "o13_pair1_pair3_selector_chart_transport_operator_axis_only_alpha13_mod_pi_strict_core_v1.json"

OUT_JSON = GENERATED / "p467_current_strict_pair123_chart_glued_projector_operator_section_cocycle_audit_probe.json"
OUT_SUMMARY = GENERATED / "p467_current_strict_pair123_chart_glued_projector_operator_section_cocycle_audit_probe_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def max_abs(a: np.ndarray) -> float:
    return float(np.max(np.abs(a)))


def as_vector(x: Any, n: int, label: str) -> np.ndarray:
    if not (isinstance(x, list) and len(x) == n and all(isinstance(v, (int, float)) for v in x)):
        raise ValueError(f"{label} must be a length-{n} numeric list")
    v = np.array([float(v) for v in x], dtype=float)
    if not np.isfinite(v).all():
        raise ValueError(f"{label} must contain only finite numbers")
    return v


def as_square_matrix(x: Any, n: int, label: str) -> np.ndarray:
    if not (isinstance(x, list) and len(x) == n and all(isinstance(row, list) and len(row) == n for row in x)):
        raise ValueError(f"{label} must be a {n}x{n} nested list")
    a = np.array(x, dtype=float)
    if not np.isfinite(a).all():
        raise ValueError(f"{label} must contain only finite numbers")
    return a


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    required = (IN_A1, IN_O12, IN_A2, IN_A3, IN_O23, IN_O13)
    missing = [str(p.relative_to(REPO)) for p in required if not p.exists()]
    if missing:
        summary = {
            "stage": "P467",
            "status": "NOT_COMPUTABLE_MISSING_INPUTS",
            "missing": missing,
            "theorem_level_pass": False,
            "no_false_pass": True,
        }
        OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    a1_obj = load_json(IN_A1)
    o12_obj = load_json(IN_O12)
    a2_obj = load_json(IN_A2)
    a3_obj = load_json(IN_A3)
    o23_obj = load_json(IN_O23)
    o13_obj = load_json(IN_O13)

    try:
        u1 = as_vector((a1_obj.get("data") or {}).get("u_1"), n=12, label="F456.A_1.u_1")
        u2 = as_vector((a2_obj.get("data") or {}).get("u_2"), n=12, label="F462.A_2.u_2")
        u3 = as_vector((a3_obj.get("data") or {}).get("u_3"), n=12, label="F464.A_3.u_3")
    except Exception as exc:
        raise SystemExit(json.dumps({"status": "INVALID_INPUT_VECTORS", "error": str(exc)}, ensure_ascii=True))

    try:
        O12 = as_square_matrix((o12_obj.get("outputs") or {}).get("O12"), n=12, label="F461.O12")
    except Exception as exc:
        raise SystemExit(json.dumps({"status": "INVALID_F461_O12_SHAPE", "error": str(exc)}, ensure_ascii=True))

    try:
        O23 = as_square_matrix((o23_obj.get("outputs") or {}).get("O23"), n=12, label="F464.O23")
    except Exception as exc:
        raise SystemExit(json.dumps({"status": "INVALID_F464_O23_SHAPE", "error": str(exc)}, ensure_ascii=True))

    try:
        O13 = as_square_matrix((o13_obj.get("outputs") or {}).get("O13"), n=12, label="F464.O13")
    except Exception as exc:
        raise SystemExit(json.dumps({"status": "INVALID_F464_O13_SHAPE", "error": str(exc)}, ensure_ascii=True))

    n = 12
    I = np.eye(n, dtype=float)

    A1 = np.outer(u1, u1)
    A2 = np.outer(u2, u2)
    A3 = np.outer(u3, u3)

    # Orthogonality / involution checks.
    o12_orth_res = max_abs(O12.T @ O12 - I)
    o23_orth_res = max_abs(O23.T @ O23 - I)
    o13_orth_res = max_abs(O13.T @ O13 - I)

    o12_invol_res = max_abs(O12 @ O12 - I)
    o23_invol_res = max_abs(O23 @ O23 - I)
    o13_invol_res = max_abs(O13 @ O13 - I)

    # Gluing checks.
    gluing12_res = max_abs(A2 - (O12 @ A1 @ O12.T))
    gluing23_res = max_abs(A3 - (O23 @ A2 @ O23.T))
    gluing13_res = max_abs(A3 - (O13 @ A1 @ O13.T))

    # Cocycle/path-independence on the exported projector section:
    # A3 via pair2 equals A3 via direct pair1->pair3.
    O23O12 = O23 @ O12
    cocycle_1_to_3_res = max_abs((O23O12 @ A1 @ O23O12.T) - (O13 @ A1 @ O13.T))

    # Sign-gauge invariance (projector-level).
    a1_sign_res = max_abs(A1 - np.outer(-u1, -u1))
    a2_sign_res = max_abs(A2 - np.outer(-u2, -u2))
    a3_sign_res = max_abs(A3 - np.outer(-u3, -u3))

    artifact = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "stage": "P467",
        "goal": "audit_pair123_projector_section_gluing_and_cocycle_path_independence",
        "inputs": {
            "a1_pair1": str(IN_A1.relative_to(REPO)),
            "o12": str(IN_O12.relative_to(REPO)),
            "a2_pair2": str(IN_A2.relative_to(REPO)),
            "o23": str(IN_O23.relative_to(REPO)),
            "o13": str(IN_O13.relative_to(REPO)),
            "a3_pair3": str(IN_A3.relative_to(REPO)),
        },
        "audits": {
            "orthogonality": {
                "o12_max_abs_residual": o12_orth_res,
                "o23_max_abs_residual": o23_orth_res,
                "o13_max_abs_residual": o13_orth_res,
            },
            "involution": {
                "o12_max_abs_residual": o12_invol_res,
                "o23_max_abs_residual": o23_invol_res,
                "o13_max_abs_residual": o13_invol_res,
            },
            "gluing": {
                "gluing12_max_abs_residual_A2_minus_O12_A1_O12T": gluing12_res,
                "gluing23_max_abs_residual_A3_minus_O23_A2_O23T": gluing23_res,
                "gluing13_max_abs_residual_A3_minus_O13_A1_O13T": gluing13_res,
            },
            "cocycle": {
                "max_abs_residual_A3_via_O23O12_minus_A3_via_O13": cocycle_1_to_3_res,
                "meaning": "path-independence on the exported projector section (sign-free)",
            },
            "sign_gauge": {
                "A1_projector_sign_invariance_max_abs": a1_sign_res,
                "A2_projector_sign_invariance_max_abs": a2_sign_res,
                "A3_projector_sign_invariance_max_abs": a3_sign_res,
            },
        },
        "hard_limits": [
            "Probe-level numeric audit only; does not export a global selector atlas nor global cocycle data on the full strict domain.",
            "Does not discharge QW-2191.",
            "Projector-level only: does not derive a sign-sensitive physical orientation datum.",
            "Does not claim strict-core selector closure / admissible S_sel_int.",
            "Does not claim ToE closure.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "P467",
        "status": "PASS_PROBE_READY",
        "generated_utc": artifact["generated_utc"],
        "audits": {
            "o23_orthogonality_max_abs_residual": o23_orth_res,
            "o13_orthogonality_max_abs_residual": o13_orth_res,
            "gluing23_max_abs_residual": gluing23_res,
            "gluing13_max_abs_residual": gluing13_res,
            "cocycle_1_to_3_max_abs_residual": cocycle_1_to_3_res,
        },
        "no_false_pass": True,
    }

    OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()

