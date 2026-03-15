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
IN_A2 = GENERATED / "a_2_pair2_orientation_projector_operator_strict_core_v1.json"
IN_A3 = GENERATED / "a_3_pair3_orientation_projector_operator_strict_core_v1.json"
IN_A4 = GENERATED / "a_4_pair4_orientation_projector_operator_strict_core_v1.json"
IN_A5 = GENERATED / "a_5_pair5_orientation_projector_operator_strict_core_v1.json"

IN_O12 = GENERATED / "o12_pair1_pair2_selector_chart_transport_operator_strict_derived_from_sigma_int_alpha12_v1.json"
IN_O23 = GENERATED / "o23_pair2_pair3_selector_chart_transport_operator_axis_only_alpha23_mod_pi_strict_core_v1.json"
IN_O13 = GENERATED / "o13_pair1_pair3_selector_chart_transport_operator_axis_only_alpha13_mod_pi_strict_core_v1.json"
IN_O34 = GENERATED / "o34_pair3_pair4_selector_chart_transport_operator_axis_only_alpha34_mod_pi_strict_core_v1.json"
IN_O45 = GENERATED / "o45_pair4_pair5_selector_chart_transport_operator_axis_only_alpha45_mod_pi_strict_core_v1.json"
IN_O24 = GENERATED / "o24_pair2_pair4_selector_chart_transport_operator_axis_only_alpha24_mod_pi_strict_core_v1.json"
IN_O35 = GENERATED / "o35_pair3_pair5_selector_chart_transport_operator_axis_only_alpha35_mod_pi_strict_core_v1.json"

OUT_JSON = GENERATED / "p468_current_strict_pair12345_chart_glued_projector_operator_section_cocycle_audit_probe.json"
OUT_SUMMARY = GENERATED / "p468_current_strict_pair12345_chart_glued_projector_operator_section_cocycle_audit_probe_summary.json"


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


def load_O(obj: dict[str, Any], *, key_candidates: list[str], label: str) -> np.ndarray:
    outputs = obj.get("outputs") or {}
    for k in key_candidates:
        if k in outputs:
            return as_square_matrix(outputs.get(k), n=12, label=f"{label}.{k}")
    raise ValueError(f"{label} missing any of keys: {key_candidates}")


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    required = (IN_A1, IN_A2, IN_A3, IN_A4, IN_A5, IN_O12, IN_O23, IN_O13, IN_O34, IN_O45, IN_O24, IN_O35)
    missing = [str(p.relative_to(REPO)) for p in required if not p.exists()]
    if missing:
        summary = {
            "stage": "P468",
            "status": "NOT_COMPUTABLE_MISSING_INPUTS",
            "missing": missing,
            "theorem_level_pass": False,
            "no_false_pass": True,
        }
        OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    a1_obj = load_json(IN_A1)
    a2_obj = load_json(IN_A2)
    a3_obj = load_json(IN_A3)
    a4_obj = load_json(IN_A4)
    a5_obj = load_json(IN_A5)

    o12_obj = load_json(IN_O12)
    o23_obj = load_json(IN_O23)
    o13_obj = load_json(IN_O13)
    o34_obj = load_json(IN_O34)
    o45_obj = load_json(IN_O45)
    o24_obj = load_json(IN_O24)
    o35_obj = load_json(IN_O35)

    try:
        u1 = as_vector((a1_obj.get("data") or {}).get("u_1"), n=12, label="A1.u_1")
        u2 = as_vector((a2_obj.get("data") or {}).get("u_2"), n=12, label="A2.u_2")
        u3 = as_vector((a3_obj.get("data") or {}).get("u_3"), n=12, label="A3.u_3")
        u4 = as_vector((a4_obj.get("data") or {}).get("u_4"), n=12, label="A4.u_4")
        u5 = as_vector((a5_obj.get("data") or {}).get("u_5"), n=12, label="A5.u_5")
    except Exception as exc:
        raise SystemExit(json.dumps({"status": "INVALID_INPUT_VECTORS", "error": str(exc)}, ensure_ascii=True))

    try:
        O12 = load_O(o12_obj, key_candidates=["O12"], label="O12")
        O23 = load_O(o23_obj, key_candidates=["O23"], label="O23")
        O13 = load_O(o13_obj, key_candidates=["O13"], label="O13")
        O34 = load_O(o34_obj, key_candidates=["O34", "O"], label="O34")
        O45 = load_O(o45_obj, key_candidates=["O45", "O"], label="O45")
        O24 = load_O(o24_obj, key_candidates=["O24", "O"], label="O24")
        O35 = load_O(o35_obj, key_candidates=["O35", "O"], label="O35")
    except Exception as exc:
        raise SystemExit(json.dumps({"status": "INVALID_TRANSPORT_OPERATOR_SHAPE", "error": str(exc)}, ensure_ascii=True))

    n = 12
    I = np.eye(n, dtype=float)

    A1 = np.outer(u1, u1)
    A2 = np.outer(u2, u2)
    A3 = np.outer(u3, u3)
    A4 = np.outer(u4, u4)
    A5 = np.outer(u5, u5)

    # Orthogonality / involution checks.
    orth = {
        "o12": max_abs(O12.T @ O12 - I),
        "o23": max_abs(O23.T @ O23 - I),
        "o13": max_abs(O13.T @ O13 - I),
        "o34": max_abs(O34.T @ O34 - I),
        "o45": max_abs(O45.T @ O45 - I),
        "o24": max_abs(O24.T @ O24 - I),
        "o35": max_abs(O35.T @ O35 - I),
    }
    invol = {
        "o12": max_abs(O12 @ O12 - I),
        "o23": max_abs(O23 @ O23 - I),
        "o13": max_abs(O13 @ O13 - I),
        "o34": max_abs(O34 @ O34 - I),
        "o45": max_abs(O45 @ O45 - I),
        "o24": max_abs(O24 @ O24 - I),
        "o35": max_abs(O35 @ O35 - I),
    }

    # Gluing checks.
    gluing = {
        "gluing12_A2_minus_O12_A1_O12T": max_abs(A2 - (O12 @ A1 @ O12.T)),
        "gluing23_A3_minus_O23_A2_O23T": max_abs(A3 - (O23 @ A2 @ O23.T)),
        "gluing13_A3_minus_O13_A1_O13T": max_abs(A3 - (O13 @ A1 @ O13.T)),
        "gluing34_A4_minus_O34_A3_O34T": max_abs(A4 - (O34 @ A3 @ O34.T)),
        "gluing45_A5_minus_O45_A4_O45T": max_abs(A5 - (O45 @ A4 @ O45.T)),
        "gluing24_A4_minus_O24_A2_O24T": max_abs(A4 - (O24 @ A2 @ O24.T)),
        "gluing35_A5_minus_O35_A3_O35T": max_abs(A5 - (O35 @ A3 @ O35.T)),
    }

    # Cocycle/path-independence checks on the exported projector section.
    O23O12 = O23 @ O12
    cocycle_1_to_3 = max_abs((O23O12 @ A1 @ O23O12.T) - (O13 @ A1 @ O13.T))

    O34O23 = O34 @ O23
    cocycle_2_to_4 = max_abs((O34O23 @ A2 @ O34O23.T) - (O24 @ A2 @ O24.T))

    O45O34 = O45 @ O34
    cocycle_3_to_5 = max_abs((O45O34 @ A3 @ O45O34.T) - (O35 @ A3 @ O35.T))

    # Sign-gauge invariance (projector-level).
    sign = {
        "a1": max_abs(A1 - np.outer(-u1, -u1)),
        "a2": max_abs(A2 - np.outer(-u2, -u2)),
        "a3": max_abs(A3 - np.outer(-u3, -u3)),
        "a4": max_abs(A4 - np.outer(-u4, -u4)),
        "a5": max_abs(A5 - np.outer(-u5, -u5)),
    }

    artifact = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "stage": "P468",
        "goal": "audit_pair12345_projector_section_gluing_and_local_cocycle_path_independence",
        "inputs": {
            "a1_pair1": str(IN_A1.relative_to(REPO)),
            "a2_pair2": str(IN_A2.relative_to(REPO)),
            "a3_pair3": str(IN_A3.relative_to(REPO)),
            "a4_pair4": str(IN_A4.relative_to(REPO)),
            "a5_pair5": str(IN_A5.relative_to(REPO)),
            "o12": str(IN_O12.relative_to(REPO)),
            "o23": str(IN_O23.relative_to(REPO)),
            "o13": str(IN_O13.relative_to(REPO)),
            "o34": str(IN_O34.relative_to(REPO)),
            "o45": str(IN_O45.relative_to(REPO)),
            "o24": str(IN_O24.relative_to(REPO)),
            "o35": str(IN_O35.relative_to(REPO)),
        },
        "audits": {
            "orthogonality": orth,
            "involution": invol,
            "gluing": gluing,
            "cocycle": {
                "pair1_pair2_pair3": cocycle_1_to_3,
                "pair2_pair3_pair4": cocycle_2_to_4,
                "pair3_pair4_pair5": cocycle_3_to_5,
                "meaning": "local path-independence on the exported projector section (sign-free)",
            },
            "sign_gauge": sign,
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
        "stage": "P468",
        "status": "PASS_PROBE_READY",
        "generated_utc": artifact["generated_utc"],
        "audits": {
            "gluing34_max_abs_residual": gluing["gluing34_A4_minus_O34_A3_O34T"],
            "gluing45_max_abs_residual": gluing["gluing45_A5_minus_O45_A4_O45T"],
            "cocycle_1_to_3_max_abs_residual": cocycle_1_to_3,
            "cocycle_2_to_4_max_abs_residual": cocycle_2_to_4,
            "cocycle_3_to_5_max_abs_residual": cocycle_3_to_5,
        },
        "no_false_pass": True,
    }

    OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()

