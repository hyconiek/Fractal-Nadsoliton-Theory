#!/usr/bin/env python3
from __future__ import annotations

import json
import math
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

OUT_JSON = GENERATED / "p466_current_strict_pair12_chart_glued_projector_operator_section_audit_probe.json"
OUT_SUMMARY = GENERATED / "p466_current_strict_pair12_chart_glued_projector_operator_section_audit_probe_summary.json"


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


def rotation_so2(alpha: float) -> np.ndarray:
    return np.array([[math.cos(alpha), -math.sin(alpha)], [math.sin(alpha), math.cos(alpha)]], dtype=float)


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    required = (IN_A1, IN_O12, IN_A2)
    missing = [str(p.relative_to(REPO)) for p in required if not p.exists()]
    if missing:
        summary = {
            "stage": "P466",
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

    try:
        u1 = as_vector((a1_obj.get("data") or {}).get("u_1"), n=12, label="F456.A_1.u_1")
        a1_2x2 = as_square_matrix((a1_obj.get("data") or {}).get("A_1_pair1_matrix_in_c1_s1"), n=2, label="F456.A_1 2x2")
    except Exception as exc:
        raise SystemExit(json.dumps({"status": "INVALID_F456_A1_SHAPE", "error": str(exc)}, ensure_ascii=True))

    try:
        u2 = as_vector((a2_obj.get("data") or {}).get("u_2"), n=12, label="F462.A_2.u_2")
        a2_2x2 = as_square_matrix((a2_obj.get("data") or {}).get("A_2_pair2_matrix_in_c2_s2"), n=2, label="F462.A_2 2x2")
    except Exception as exc:
        raise SystemExit(json.dumps({"status": "INVALID_F462_A2_SHAPE", "error": str(exc)}, ensure_ascii=True))

    try:
        outputs = o12_obj.get("outputs") or {}
        alpha12 = float(outputs.get("alpha12_mod_2pi"))
        O12 = as_square_matrix(outputs.get("O12"), n=12, label="F461.O12")
        G12 = as_square_matrix(outputs.get("G12_so2"), n=2, label="F461.G12_so2")
    except Exception as exc:
        raise SystemExit(json.dumps({"status": "INVALID_F461_O12_SHAPE", "error": str(exc)}, ensure_ascii=True))

    n = 12
    I = np.eye(n, dtype=float)

    A1_full = np.outer(u1, u1)
    A2_full = np.outer(u2, u2)

    # Full-carrier gluing check.
    A2_from_transport = O12 @ A1_full @ O12.T
    gluing_full_res = max_abs(A2_full - A2_from_transport)

    # Pair-plane 2x2 gluing check.
    G12_from_alpha = rotation_so2(alpha12)
    gluing_2x2_res = max_abs(a2_2x2 - (G12 @ a1_2x2 @ G12.T))
    g12_diff = max_abs(G12 - G12_from_alpha)

    # Orthogonality and involution.
    orth_res = max_abs(O12.T @ O12 - I)
    invol_res = max_abs(O12 @ O12 - I)

    # Sign-gauge checks.
    a1_sign_res = max_abs(A1_full - np.outer(-u1, -u1))
    a2_sign_res = max_abs(A2_full - np.outer(-u2, -u2))

    artifact = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "stage": "P466",
        "goal": "audit_consistency_of_pair12_chart_glued_projector_operator_section_A2_equals_O12_A1_O12T",
        "inputs": {
            "a1_pair1": str(IN_A1.relative_to(REPO)),
            "o12": str(IN_O12.relative_to(REPO)),
            "a2_pair2": str(IN_A2.relative_to(REPO)),
        },
        "audits": {
            "o12": {
                "alpha12_mod_2pi": alpha12,
                "orthogonality_max_abs_residual": orth_res,
                "involution_max_abs_residual": invol_res,
                "g12_max_abs_diff_vs_rotation(alpha12)": g12_diff,
            },
            "gluing": {
                "full_carrier_max_abs_residual_A2_minus_O12_A1_O12T": gluing_full_res,
                "pair_plane_2x2_max_abs_residual_A2_minus_G12_A1_G12T": gluing_2x2_res,
            },
            "sign_gauge": {
                "A1_projector_sign_invariance_max_abs": a1_sign_res,
                "A2_projector_sign_invariance_max_abs": a2_sign_res,
            },
        },
        "hard_limits": [
            "Lane-scoped audit only; does not export a global selector atlas nor global cocycle data.",
            "Does not discharge QW-2191.",
            "Projector-level only: does not derive a sign-sensitive physical orientation datum.",
            "Does not claim strict-core selector closure / admissible S_sel_int.",
            "Does not claim ToE closure.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "P466",
        "status": "PASS_PROBE_READY",
        "generated_utc": artifact["generated_utc"],
        "audits": {
            "orthogonality_max_abs_residual": orth_res,
            "gluing_full_max_abs_residual": gluing_full_res,
            "gluing_pair_plane_2x2_max_abs_residual": gluing_2x2_res,
            "a1_sign_invariance_max_abs": a1_sign_res,
            "a2_sign_invariance_max_abs": a2_sign_res,
        },
        "no_false_pass": True,
    }

    OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()

