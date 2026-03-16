#!/usr/bin/env python3
from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np

AS_OF = "2026-03-16"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

STRICT_PROJECTIVE_STATE = GENERATED / "selector_state_global_c_v1_projective_strict_v1.json"
EXT_ORIENTED_STATE = GENERATED / "strict_extension_lane_selector_state_global_c_v1_oriented_vector_v1.json"

OUT_JSON = (
    GENERATED
    / "p473_current_strict_extension_lane_global_oriented_selector_state_projector_consistency_audit_probe.json"
)
OUT_SUMMARY = (
    GENERATED
    / "p473_current_strict_extension_lane_global_oriented_selector_state_projector_consistency_audit_probe_summary.json"
)


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def fourier_pair_basis(n: int, m: int) -> tuple[np.ndarray, np.ndarray]:
    scale = math.sqrt(2.0 / float(n))
    c = np.array(
        [scale * math.cos(2.0 * math.pi * m * k / float(n)) for k in range(n)],
        dtype=float,
    )
    s = np.array(
        [scale * math.sin(2.0 * math.pi * m * k / float(n)) for k in range(n)],
        dtype=float,
    )
    return c, s


def projector_from_coords(coords: np.ndarray) -> np.ndarray:
    v = coords.reshape(2, 1)
    return v @ v.T


def max_abs(a: np.ndarray) -> float:
    return float(np.max(np.abs(a)).item())


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    missing: list[str] = []
    for p in (STRICT_PROJECTIVE_STATE, EXT_ORIENTED_STATE):
        if not p.exists():
            missing.append(str(p.relative_to(REPO)))
    if missing:
        raise SystemExit(
            json.dumps(
                {
                    "stage": "P473",
                    "status": "NOT_COMPUTABLE_MISSING_INPUTS",
                    "as_of": AS_OF,
                    "missing": missing,
                    "expected": [
                        str(STRICT_PROJECTIVE_STATE.relative_to(REPO)),
                        str(EXT_ORIENTED_STATE.relative_to(REPO)),
                    ],
                    "no_false_pass": True,
                },
                ensure_ascii=True,
            )
        )

    strict_state = load_json(STRICT_PROJECTIVE_STATE)
    section_ref = strict_state.get("projector_section_ref")
    if not isinstance(section_ref, str) or not section_ref:
        raise SystemExit(
            json.dumps(
                {
                    "stage": "P473",
                    "status": "INVALID_STRICT_PROJECTIVE_STATE_MISSING_PROJECTOR_SECTION_REF",
                    "as_of": AS_OF,
                    "strict_state_ref": str(STRICT_PROJECTIVE_STATE.relative_to(REPO)),
                    "no_false_pass": True,
                },
                ensure_ascii=True,
            )
        )

    section_path = REPO / section_ref
    if not section_path.exists():
        raise SystemExit(
            json.dumps(
                {
                    "stage": "P473",
                    "status": "NOT_COMPUTABLE_MISSING_PROJECTOR_SECTION_OBJECT",
                    "as_of": AS_OF,
                    "projector_section_ref": section_ref,
                    "missing": [section_ref],
                    "no_false_pass": True,
                },
                ensure_ascii=True,
            )
        )

    section = load_json(section_path)
    section_inputs = section.get("inputs") or {}
    if not isinstance(section_inputs, dict):
        section_inputs = {}

    ext_state = load_json(EXT_ORIENTED_STATE)
    outputs = ext_state.get("outputs") or {}
    u_vectors = (outputs.get("u_vectors_oriented_extension_lane") if isinstance(outputs, dict) else None) or {}
    if not isinstance(u_vectors, dict):
        u_vectors = {}

    n = 12
    tol_norm = 1e-9
    tol_proj = 1e-12

    pair_results: dict[str, Any] = {}
    overall_pass = True
    for m in range(1, 6):
        pair_id = f"pair{m}"
        u_key = f"u_{m}"
        u_raw = u_vectors.get(u_key)
        if not (
            isinstance(u_raw, list)
            and len(u_raw) == n
            and all(isinstance(v, (int, float)) and math.isfinite(float(v)) for v in u_raw)
        ):
            pair_results[pair_id] = {
                "status": "INVALID_EXTENSION_VECTOR_SHAPE",
                "u_key": u_key,
                "expected": f"length-{n} numeric list",
            }
            overall_pass = False
            continue

        u = np.array([float(v) for v in u_raw], dtype=float)
        u_norm = float(np.linalg.norm(u).item())

        c, s = fourier_pair_basis(n=n, m=m)
        coords = np.array([float(np.dot(u, c)), float(np.dot(u, s))], dtype=float)
        coords_norm = float(np.linalg.norm(coords).item())
        P_ext = projector_from_coords(coords)
        P_ext_minus = projector_from_coords(-coords)
        minus_invariance_max_abs = max_abs(P_ext - P_ext_minus)

        a_ref_key = f"a{m}_pair{m}_operator_ref"
        a_ref = section_inputs.get(a_ref_key)
        if not isinstance(a_ref, str) or not a_ref:
            pair_results[pair_id] = {
                "status": "MISSING_STRICT_PROJECTOR_OPERATOR_REF",
                "projector_section_ref": section_ref,
                "missing_key": a_ref_key,
            }
            overall_pass = False
            continue

        a_path = REPO / a_ref
        if not a_path.exists():
            pair_results[pair_id] = {
                "status": "MISSING_STRICT_PROJECTOR_OPERATOR_OBJECT",
                "a_ref": a_ref,
            }
            overall_pass = False
            continue

        a_obj = load_json(a_path)
        matrix_key = f"A_{m}_pair{m}_matrix_in_c{m}_s{m}"
        A_raw = ((a_obj.get("data") or {}) if isinstance(a_obj, dict) else {}).get(matrix_key)
        if not (
            isinstance(A_raw, list)
            and len(A_raw) == 2
            and all(isinstance(row, list) and len(row) == 2 for row in A_raw)
            and all(isinstance(v, (int, float)) and math.isfinite(float(v)) for row in A_raw for v in row)
        ):
            pair_results[pair_id] = {
                "status": "INVALID_STRICT_PROJECTOR_MATRIX_SHAPE",
                "a_ref": a_ref,
                "matrix_key": matrix_key,
            }
            overall_pass = False
            continue

        P_strict = np.array([[float(A_raw[0][0]), float(A_raw[0][1])], [float(A_raw[1][0]), float(A_raw[1][1])]])
        diff = P_ext - P_strict
        diff_max_abs = max_abs(diff)

        pair_pass = (
            abs(u_norm - 1.0) <= tol_norm
            and abs(coords_norm - 1.0) <= tol_norm
            and diff_max_abs <= tol_proj
            and minus_invariance_max_abs <= tol_proj
        )
        overall_pass = overall_pass and pair_pass

        pair_results[pair_id] = {
            "status": "PASS" if pair_pass else "FAIL",
            "u_norm": float(u_norm),
            "coords_in_c_s": [float(coords[0]), float(coords[1])],
            "coords_norm": float(coords_norm),
            "projector_ext": [[float(P_ext[0, 0]), float(P_ext[0, 1])], [float(P_ext[1, 0]), float(P_ext[1, 1])]],
            "projector_strict_ref": a_ref,
            "projector_strict": [
                [float(P_strict[0, 0]), float(P_strict[0, 1])],
                [float(P_strict[1, 0]), float(P_strict[1, 1])],
            ],
            "projector_diff_max_abs": float(diff_max_abs),
            "minus_invariance_projector_diff_max_abs": float(minus_invariance_max_abs),
            "tolerances": {"norm_abs_tol": tol_norm, "projector_max_abs_tol": tol_proj},
        }

    artifact = {
        "stage": "P473",
        "as_of": AS_OF,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "goal": (
            "audit_that_the_extension_lane_global_oriented_vector_selector_state_object (AX29) "
            "is_projector_consistent_with_the_strict_global_projective_selector_state_object (F470) "
            "on_each_pair_m_chart_basis (c_m,s_m), so_the_extension_lane_only_fixes_a_sign_gauge_representative "
            "without_changing_the_underlying_strict_projective_state"
        ),
        "inputs": {
            "strict_projective_state_ref": str(STRICT_PROJECTIVE_STATE.relative_to(REPO)),
            "projector_section_ref": section_ref,
            "extension_oriented_vector_state_ref": str(EXT_ORIENTED_STATE.relative_to(REPO)),
        },
        "result": {
            "overall_pass": bool(overall_pass),
            "pairs": pair_results,
        },
        "no_false_pass": True,
        "hard_limits": [
            "Probe-only: does not discharge H37 and does not export a strict sign-sensitive physical orientation datum.",
            "Audits only projector consistency in the declared (c_m,s_m) chart bases; does not promote any physical interpretation.",
            "Does not claim strict-core selector closure / admissible S_sel_int.",
            "Does not claim global discharge of QW-2191.",
            "Does not claim ToE closure.",
        ],
    }

    summary = {
        "stage": "P473",
        "status": "PASS_PROJECTOR_CONSISTENT" if overall_pass else "FAIL_PROJECTOR_INCONSISTENT",
        "overall_pass": bool(overall_pass),
        "no_false_pass": True,
    }

    OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_JSON)


if __name__ == "__main__":
    main()

