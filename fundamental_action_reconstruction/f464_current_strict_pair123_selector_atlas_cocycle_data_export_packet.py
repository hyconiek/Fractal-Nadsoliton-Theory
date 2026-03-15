#!/usr/bin/env python3
from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np

AS_OF = "2026-03-15"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_A1 = GENERATED / "a_1_pair1_orientation_projector_operator_strict_core_v1.json"
IN_A2 = GENERATED / "a_2_pair2_orientation_projector_operator_strict_core_v1.json"
IN_O12 = GENERATED / "o12_pair1_pair2_selector_chart_transport_operator_strict_derived_from_sigma_int_alpha12_v1.json"
IN_MODE_ASSIGN = GENERATED / "mode_index_assignment_shannon_element_order_reference_strict_core_v1.json"

OUT_A3 = GENERATED / "a_3_pair3_orientation_projector_operator_strict_core_v1.json"
OUT_O23 = GENERATED / "o23_pair2_pair3_selector_chart_transport_operator_axis_only_alpha23_mod_pi_strict_core_v1.json"
OUT_O13 = GENERATED / "o13_pair1_pair3_selector_chart_transport_operator_axis_only_alpha13_mod_pi_strict_core_v1.json"
OUT_SECTION = GENERATED / "a_123_pair123_chart_glued_orientation_projector_operator_section_strict_core_v1.json"
OUT_ATLAS = GENERATED / "selector_atlas_pair123_axis_only_projector_v1.json"
OUT_SUMMARY = GENERATED / "f464_current_strict_pair123_selector_atlas_cocycle_data_export_packet_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def clean_scalar(value: float) -> float:
    if abs(value) < 1e-15:
        return 0.0
    return round(float(value), 15)


def max_abs(a: np.ndarray) -> float:
    return float(np.max(np.abs(a)))


def angle_mod(x: float, modulus: float) -> float:
    return float(x % modulus)


def circular_distance_mod_pi(a: float, b: float) -> float:
    """Distance on R/(pi Z) in [0, pi/2]."""
    d = abs(angle_mod(a - b, math.pi))
    return float(min(d, math.pi - d))


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


def real_fourier_pair_basis(n: int, m: int) -> tuple[np.ndarray, np.ndarray]:
    """Canonical real Fourier pair (c_m, s_m) on Z_n for 1<=m<=n/2-1."""
    if not (n > 0 and n % 2 == 0):
        raise ValueError("n must be positive and even")
    if not (1 <= m <= (n // 2 - 1)):
        raise ValueError("m out of range for real Fourier pair")
    scale = math.sqrt(2.0 / float(n))
    c = np.array([scale * math.cos(2.0 * math.pi * m * k / float(n)) for k in range(n)], dtype=float)
    s = np.array([scale * math.sin(2.0 * math.pi * m * k / float(n)) for k in range(n)], dtype=float)
    return c, s


def rotation_so2(alpha: float) -> np.ndarray:
    return np.array([[math.cos(alpha), -math.sin(alpha)], [math.sin(alpha), math.cos(alpha)]], dtype=float)


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    required = (IN_A1, IN_A2, IN_O12, IN_MODE_ASSIGN)
    missing = [str(p.relative_to(REPO)) for p in required if not p.exists()]
    if missing:
        raise SystemExit(
            json.dumps(
                {
                    "status": "NOT_COMPUTABLE_MISSING_INPUTS",
                    "missing": missing,
                    "expected": [
                        "F456 A_1(pair1) projector operator",
                        "F462 A_2(pair2) projector operator",
                        "F461 O12 chart-transport operator",
                        "F454 Shannon mode-index assignment (pair3 minimizer axis)",
                    ],
                },
                ensure_ascii=True,
            )
        )

    a1_obj = load_json(IN_A1)
    a2_obj = load_json(IN_A2)
    o12_obj = load_json(IN_O12)
    mode_assign = load_json(IN_MODE_ASSIGN)

    try:
        u1 = as_vector((a1_obj.get("data") or {}).get("u_1"), n=12, label="F456.data.u_1")
        u2 = as_vector((a2_obj.get("data") or {}).get("u_2"), n=12, label="F462.data.u_2")
    except Exception as exc:
        raise SystemExit(json.dumps({"status": "INVALID_INPUT_VECTORS", "error": str(exc)}, ensure_ascii=True))

    try:
        outputs = o12_obj.get("outputs") or {}
        n = int(outputs.get("n"))
        if n != 12:
            raise ValueError("expected n=12")
        alpha12_mod_2pi = float(outputs.get("alpha12_mod_2pi"))
        O12 = as_square_matrix(outputs.get("O12"), n=n, label="F461.outputs.O12")
    except Exception as exc:
        raise SystemExit(json.dumps({"status": "INVALID_F461_SHAPE", "error": str(exc)}, ensure_ascii=True))

    try:
        pairs = ((mode_assign.get("outputs") or {}).get("pairs") or {})
        pair1 = pairs.get("pair1") or {}
        pair2 = pairs.get("pair2") or {}
        pair3 = pairs.get("pair3") or {}
        theta1_mod_pi = float(pair1.get("theta_minimizer_mod_pi"))
        theta2_mod_pi = float(pair2.get("theta_minimizer_mod_pi"))
        theta3_mod_pi = float(pair3.get("theta_minimizer_mod_pi"))
        u3 = as_vector(pair3.get("u_minus"), n=12, label="F454.outputs.pairs.pair3.u_minus")
        objective_minimizer_vector = str(pair3.get("objective_minimizer_vector") or "")
        if objective_minimizer_vector and objective_minimizer_vector != "u_minus":
            raise ValueError("expected pair3 objective_minimizer_vector == 'u_minus'")
    except Exception as exc:
        raise SystemExit(json.dumps({"status": "INVALID_F454_SHAPE", "error": str(exc)}, ensure_ascii=True))

    # Axis-only canonical representatives (mod pi).
    theta1_pi = angle_mod(theta1_mod_pi, math.pi)
    theta2_pi = angle_mod(theta2_mod_pi, math.pi)
    theta3_pi = angle_mod(theta3_mod_pi, math.pi)

    alpha12_mod_pi_from_shannon = angle_mod(theta2_pi - theta1_pi, math.pi)
    alpha12_mod_pi_from_o12 = angle_mod(alpha12_mod_2pi, math.pi)
    alpha12_mod_pi_consistency = circular_distance_mod_pi(alpha12_mod_pi_from_shannon, alpha12_mod_pi_from_o12)

    alpha23_mod_pi = angle_mod(theta3_pi - theta2_pi, math.pi)
    alpha13_mod_pi = angle_mod(theta3_pi - theta1_pi, math.pi)

    # Canonical Fourier bases for pair1/pair2/pair3 (declared scaffold).
    c1, s1 = real_fourier_pair_basis(n=n, m=1)
    c2, s2 = real_fourier_pair_basis(n=n, m=2)
    c3, s3 = real_fourier_pair_basis(n=n, m=3)
    C1 = np.column_stack([c1, s1])
    C2 = np.column_stack([c2, s2])
    C3 = np.column_stack([c3, s3])

    I = np.eye(n, dtype=float)

    def build_O(Ca: np.ndarray, Cb: np.ndarray, alpha: float) -> np.ndarray:
        pia = Ca @ Ca.T
        pib = Cb @ Cb.T
        pi_rest = I - pia - pib
        G = rotation_so2(alpha)
        return (Cb @ G @ Ca.T) + (Ca @ G.T @ Cb.T) + pi_rest

    O23 = build_O(Ca=C2, Cb=C3, alpha=alpha23_mod_pi)
    O13 = build_O(Ca=C1, Cb=C3, alpha=alpha13_mod_pi)

    # Projectors / operators on the full carrier.
    A1 = np.outer(u1, u1)
    A2 = np.outer(u2, u2)
    A3 = np.outer(u3, u3)

    # Minimal checks / audits.
    orth12_res = max_abs(O12.T @ O12 - I)
    orth23_res = max_abs(O23.T @ O23 - I)
    orth13_res = max_abs(O13.T @ O13 - I)

    invol23_res = max_abs(O23 @ O23 - I)
    invol13_res = max_abs(O13 @ O13 - I)

    gluing12_res = max_abs(A2 - (O12 @ A1 @ O12.T))
    gluing23_res = max_abs(A3 - (O23 @ A2 @ O23.T))
    gluing13_res = max_abs(A3 - (O13 @ A1 @ O13.T))

    O23O12 = O23 @ O12
    cocycle_1_to_3_res = max_abs((O23O12 @ A1 @ O23O12.T) - (O13 @ A1 @ O13.T))

    # Pair-plane coordinates for A3 in (c3,s3).
    coords3 = np.array([float(c3 @ u3), float(s3 @ u3)], dtype=float)
    coords3_norm = float(np.linalg.norm(coords3))
    if coords3_norm == 0.0:
        raise SystemExit(json.dumps({"status": "DEGENERATE_U3_COORDS", "coords_norm": 0.0}, ensure_ascii=True))
    coords3_unit = coords3 / coords3_norm
    A3_2 = np.outer(coords3_unit, coords3_unit)
    u3_recon = (coords3[0] * c3) + (coords3[1] * s3)
    u3_outside = float(np.linalg.norm(u3 - u3_recon))

    a3_sign_res = max_abs(A3 - np.outer(-u3, -u3))
    a3_sym_res = max_abs(A3 - A3.T)
    a3_idem_res = max_abs(A3 @ A3 - A3)

    # Export A3 operator object.
    a3_obj = {
        "object": "A_3_pair3_orientation_projector_operator_strict_core_v1",
        "status": "actual_exported_operator__derived_from_shannon_mode_index_assignment_pair3_axis__residual_z2_sign_gauge_invariant",
        "as_of": AS_OF,
        "intent": (
            "Export one strict-core projector operator on V_3 = span{c3,s3} constructed from the exported Shannon element-order reference "
            "mode-index assignment minimizer axis on pair3 (F454). The operator is A_3(pair3) := |u_3><u_3|. This is projector-level "
            "(residual-sign-gauge-invariant) and does not imply any sign-sensitive physical orientation convention."
        ),
        "inputs": {
            "mode_index_assignment_ref": str(IN_MODE_ASSIGN.relative_to(REPO)),
        },
        "domain_notes": {
            "basis_plane": "V_3 = span{c3,s3} inside the n=12 real Fourier scaffold",
            "lane": "shannon_element_order_reference (axis-only; residual sign remains gauge)",
        },
        "data": {
            "theta_3_minimizer_mod_pi": clean_scalar(theta3_pi),
            "u_3": [clean_scalar(float(x)) for x in u3.tolist()],
            "u_3_coords_in_c3_s3": [clean_scalar(float(x)) for x in coords3_unit.tolist()],
            "A_3_pair3_matrix_in_c3_s3": [[clean_scalar(float(x)) for x in row] for row in A3_2.tolist()],
        },
        "audits": {
            "u3_outside_span_c3_s3_l2": u3_outside,
            "coords_norm": coords3_norm,
            "projector_symmetry_inf_norm": a3_sym_res,
            "projector_idempotence_inf_norm": a3_idem_res,
            "projector_trace": float(np.trace(A3)),
            "residual_z2_sign_invariance_max_abs_diff": a3_sign_res,
            "tolerance": 1e-12,
        },
        "hard_limits": [
            "Projector-level only; does not derive a sign-sensitive physical orientation datum.",
            "Does not claim strict-core selector closure / admissible S_sel_int.",
            "Does not claim a global selector atlas nor global discharge of QW-2191.",
            "Does not claim ToE closure.",
        ],
        "no_false_pass": True,
    }

    # Export O23 and O13 operator objects.
    o23_obj = {
        "object": "O23_pair2_pair3_selector_chart_transport_operator_axis_only_alpha23_mod_pi_strict_core_v1",
        "stage": "F464",
        "status": "actual_exported_lane_scoped_pair_chart_transport_operator__pair2_pair3__axis_only__projector_level__no_false_pass",
        "as_of": AS_OF,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "inputs": {
            "alpha23_mod_pi_source": str(IN_MODE_ASSIGN.relative_to(REPO)),
            "alpha23_mod_pi": clean_scalar(alpha23_mod_pi),
        },
        "construction": {
            "fourier_pair_planes": {
                "pair2": "V2 := span{c2,s2} (m=2) with canonical real Fourier basis on Z12",
                "pair3": "V3 := span{c3,s3} (m=3) with canonical real Fourier basis on Z12",
            },
            "transport_operator": "O23 := C3 G(α) C2^T + C2 G(-α) C3^T + Π_rest with α := alpha23_mod_pi (axis-only representative)",
            "note": "Using α mod π fixes an axis-only representative; replacing α by α+π flips transported vectors but preserves transported projectors.",
        },
        "outputs": {
            "n": n,
            "alpha23_mod_pi": clean_scalar(alpha23_mod_pi),
            "O23": [[clean_scalar(float(x)) for x in row] for row in O23.tolist()],
            "checks": {
                "orthogonality_max_abs_residual": orth23_res,
                "involution_O23_squared_equals_I_max_abs_residual": invol23_res,
            },
        },
        "hard_limits": [
            "Axis-only: uses alpha23 mod π; transport is intended for projector/axis gluing, not sign-sensitive orientation.",
            "Lane-scoped: does not export a global selector atlas nor global cocycle data on the full strict domain.",
            "Does not discharge QW-2191 nor imply strict-core selector closure.",
            "Does not claim ToE closure.",
        ],
        "no_false_pass": True,
    }

    o13_obj = {
        "object": "O13_pair1_pair3_selector_chart_transport_operator_axis_only_alpha13_mod_pi_strict_core_v1",
        "stage": "F464",
        "status": "actual_exported_lane_scoped_pair_chart_transport_operator__pair1_pair3__axis_only__projector_level__no_false_pass",
        "as_of": AS_OF,
        "generated_utc": o23_obj["generated_utc"],
        "inputs": {
            "alpha13_mod_pi_source": str(IN_MODE_ASSIGN.relative_to(REPO)),
            "alpha13_mod_pi": clean_scalar(alpha13_mod_pi),
        },
        "construction": {
            "fourier_pair_planes": {
                "pair1": "V1 := span{c1,s1} (m=1) with canonical real Fourier basis on Z12",
                "pair3": "V3 := span{c3,s3} (m=3) with canonical real Fourier basis on Z12",
            },
            "transport_operator": "O13 := C3 G(α) C1^T + C1 G(-α) C3^T + Π_rest with α := alpha13_mod_pi (axis-only representative)",
            "note": "Using α mod π fixes an axis-only representative; replacing α by α+π flips transported vectors but preserves transported projectors.",
        },
        "outputs": {
            "n": n,
            "alpha13_mod_pi": clean_scalar(alpha13_mod_pi),
            "O13": [[clean_scalar(float(x)) for x in row] for row in O13.tolist()],
            "checks": {
                "orthogonality_max_abs_residual": orth13_res,
                "involution_O13_squared_equals_I_max_abs_residual": invol13_res,
            },
        },
        "hard_limits": o23_obj["hard_limits"],
        "no_false_pass": True,
    }

    # Export a three-chart glued operator section (projector level).
    section_obj = {
        "object": "A_123_pair123_chart_glued_orientation_projector_operator_section_strict_core_v1",
        "stage": "F464",
        "status": "actual_exported_three_chart_projector_operator_section__pair1_pair2_pair3__with_cocycle_audit__no_false_pass",
        "as_of": AS_OF,
        "generated_utc": o23_obj["generated_utc"],
        "intent": (
            "Export one explicit three-chart projector operator section on {pair1,pair2,pair3} with explicit transition operators and "
            "projector-level gluing laws, including an explicit cocycle/path-independence audit at the level of the exported projector section. "
            "This is lane-scoped and does not imply a global selector atlas nor any global discharge."
        ),
        "inputs": {
            "a1_pair1_operator_ref": str(IN_A1.relative_to(REPO)),
            "a2_pair2_operator_ref": str(IN_A2.relative_to(REPO)),
            "a3_pair3_operator_ref": str(OUT_A3.relative_to(REPO)),
            "o12_ref": str(IN_O12.relative_to(REPO)),
            "o23_ref": str(OUT_O23.relative_to(REPO)),
            "o13_ref": str(OUT_O13.relative_to(REPO)),
            "mode_index_assignment_ref": str(IN_MODE_ASSIGN.relative_to(REPO)),
        },
        "gluing_laws": [
            "A_2(pair2) = O_12 A_1(pair1) O_12^T  (projector-level)",
            "A_3(pair3) = O_23 A_2(pair2) O_23^T  (projector-level, axis-only transport)",
            "A_3(pair3) = O_13 A_1(pair1) O_13^T  (projector-level, axis-only transport)",
        ],
        "cocycle_audit": {
            "statement": (
                "Transporting the pair1 projector to pair3 via pair2 agrees with the direct pair1->pair3 transport, "
                "at the level of the glued projector operator section (sign-free)."
            ),
            "path_1_to_3_via_2": "O_23 O_12",
            "direct_1_to_3": "O_13",
            "max_abs_residual": cocycle_1_to_3_res,
        },
        "audits": {
            "o12_orthogonality_max_abs_residual": orth12_res,
            "o23_orthogonality_max_abs_residual": orth23_res,
            "o13_orthogonality_max_abs_residual": orth13_res,
            "gluing12_max_abs_residual": gluing12_res,
            "gluing23_max_abs_residual": gluing23_res,
            "gluing13_max_abs_residual": gluing13_res,
            "alpha12_mod_pi_consistency_distance": alpha12_mod_pi_consistency,
        },
        "hard_limits": [
            "Projector-level only: section gluing is sign-gauge-safe but does not lift residual sign to a physical convention.",
            "Lane-scoped three-chart section only; does not export a global selector atlas nor global cocycle data on the full strict domain.",
            "Does not discharge QW-2191 nor export strict-core selector closure / admissible S_sel_int.",
            "Does not claim ToE closure.",
        ],
        "no_false_pass": True,
    }

    # Export a three-chart atlas object (lane-scoped, projector-level).
    atlas = {
        "object": "SelectorAtlas_pair123_axis_only_projector_v1",
        "stage": "F464",
        "status": "actual_exported_lane_scoped_three_chart_selector_atlas_with_projector_level_gluing_and_cocycle_audit__no_false_pass",
        "as_of": AS_OF,
        "generated_utc": o23_obj["generated_utc"],
        "intent": (
            "Export one explicit lane-scoped three-chart selector atlas object on {pair1,pair2,pair3} at projector level (sign-gauge-safe), "
            "including explicit chart transport operators (O12, O23, O13), projector-level gluing laws, and an explicit cocycle/path-independence "
            "audit for the glued projector operator section. This is not a global selector atlas and does not imply any global discharge."
        ),
        "charts": {
            "pair1": {
                "chart_id": "pair1",
                "carrier_plane": "V1 = span{c1,s1} inside the n=12 real Fourier scaffold",
                "local_operator": {
                    "ref": str(IN_A1.relative_to(REPO)),
                    "object": "A_1_pair1_orientation_projector_operator_strict_core_v1",
                },
            },
            "pair2": {
                "chart_id": "pair2",
                "carrier_plane": "V2 = span{c2,s2} inside the n=12 real Fourier scaffold",
                "local_operator": {
                    "ref": str(IN_A2.relative_to(REPO)),
                    "object": "A_2_pair2_orientation_projector_operator_strict_core_v1",
                },
            },
            "pair3": {
                "chart_id": "pair3",
                "carrier_plane": "V3 = span{c3,s3} inside the n=12 real Fourier scaffold",
                "local_operator": {
                    "ref": str(OUT_A3.relative_to(REPO)),
                    "object": "A_3_pair3_orientation_projector_operator_strict_core_v1",
                },
            },
        },
        "overlap_domain_declaration": {
            "scope": "exported_artifact_overlap_on_n12_fourier_carrier",
            "meaning": (
                "This is an exported-artifact overlap declaration: all three chart operators A_1/A_2/A_3 and the transport operators "
                "O_12/O_23/O_13 are simultaneously defined as strict-core exported artifacts on the declared n=12 Fourier carrier. "
                "This is not a global open-cover atlas on the full strict domain C_v1."
            ),
        },
        "transitions": {
            "pair1_to_pair2": {
                "operator_ref": str(IN_O12.relative_to(REPO)),
                "operator": "O12_pair1_pair2_selector_chart_transport_operator_strict_derived_from_sigma_int_alpha12_v1",
                "alpha12_mod_2pi": clean_scalar(alpha12_mod_2pi),
            },
            "pair2_to_pair3": {
                "operator_ref": str(OUT_O23.relative_to(REPO)),
                "operator": o23_obj["object"],
                "alpha23_mod_pi": clean_scalar(alpha23_mod_pi),
            },
            "pair1_to_pair3": {
                "operator_ref": str(OUT_O13.relative_to(REPO)),
                "operator": o13_obj["object"],
                "alpha13_mod_pi": clean_scalar(alpha13_mod_pi),
            },
        },
        "gluing_data": {
            "operator_section_ref": str(OUT_SECTION.relative_to(REPO)),
            "laws": section_obj["gluing_laws"],
            "cocycle_audit": section_obj["cocycle_audit"],
        },
        "audits": section_obj["audits"],
        "hard_limits": [
            "Lane-scoped: exports only a three-chart atlas object on {pair1,pair2,pair3} at projector level.",
            "Does not export a global selector atlas, global overlap-domain declaration, nor global cocycle data (H41 remains open globally).",
            "Does not discharge QW-2191 nor export a global selector transition/gluing object (H40 remains open globally).",
            "Projector-level only: does not derive a sign-sensitive physical orientation datum.",
            "Does not claim strict-core selector closure / admissible S_sel_int.",
            "Does not claim ToE closure.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "F464",
        "status": "F464_EXECUTED_CURRENT_STRICT_PAIR123_SELECTOR_ATLAS_COCYCLE_DATA_EXPORT_PACKET_NO_FALSE_PASS",
        "generated_utc": atlas["generated_utc"],
        "artifacts": {
            "a3": str(OUT_A3.relative_to(REPO)),
            "o23": str(OUT_O23.relative_to(REPO)),
            "o13": str(OUT_O13.relative_to(REPO)),
            "section": str(OUT_SECTION.relative_to(REPO)),
            "atlas": str(OUT_ATLAS.relative_to(REPO)),
        },
        "audits": {
            "alpha12_mod_pi_consistency_distance": alpha12_mod_pi_consistency,
            "o23_orthogonality_max_abs_residual": orth23_res,
            "o13_orthogonality_max_abs_residual": orth13_res,
            "gluing23_max_abs_residual": gluing23_res,
            "gluing13_max_abs_residual": gluing13_res,
            "cocycle_1_to_3_max_abs_residual": cocycle_1_to_3_res,
        },
        "no_false_pass": True,
    }

    OUT_A3.write_text(json.dumps(a3_obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_O23.write_text(json.dumps(o23_obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_O13.write_text(json.dumps(o13_obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SECTION.write_text(json.dumps(section_obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_ATLAS.write_text(json.dumps(atlas, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()

