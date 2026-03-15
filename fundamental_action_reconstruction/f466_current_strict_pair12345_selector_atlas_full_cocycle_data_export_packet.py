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

OUT_O14 = GENERATED / "o14_pair1_pair4_selector_chart_transport_operator_axis_only_alpha14_mod_pi_strict_core_v1.json"
OUT_O15 = GENERATED / "o15_pair1_pair5_selector_chart_transport_operator_axis_only_alpha15_mod_pi_strict_core_v1.json"
OUT_O25 = GENERATED / "o25_pair2_pair5_selector_chart_transport_operator_axis_only_alpha25_mod_pi_strict_core_v1.json"

OUT_SECTION = GENERATED / "a_12345_pair12345_chart_glued_orientation_projector_operator_section_strict_core_v2.json"
OUT_ATLAS = GENERATED / "selector_atlas_pair12345_axis_only_projector_v2.json"
OUT_SUMMARY = GENERATED / "f466_current_strict_pair12345_selector_atlas_full_cocycle_data_export_packet_summary.json"


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


def build_swap_transport(I: np.ndarray, Ca: np.ndarray, Cb: np.ndarray, alpha: float) -> np.ndarray:
    pia = Ca @ Ca.T
    pib = Cb @ Cb.T
    pi_rest = I - pia - pib
    G = rotation_so2(alpha)
    return (Cb @ G @ Ca.T) + (Ca @ G.T @ Cb.T) + pi_rest


def theta_mod_pi_from_u(u: np.ndarray, c: np.ndarray, s: np.ndarray) -> float:
    coords = np.array([float(c @ u), float(s @ u)], dtype=float)
    norm = float(np.linalg.norm(coords))
    if norm == 0.0:
        raise ValueError("degenerate coordinates (zero)")
    theta = math.atan2(float(coords[1]), float(coords[0]))
    return angle_mod(theta, math.pi)


def export_axis_only_transport_operator(
    *,
    object_name: str,
    m_from: int,
    m_to: int,
    alpha_mod_pi: float,
    Ca: np.ndarray,
    Cb: np.ndarray,
    I: np.ndarray,
    out_path: Path,
    generated_utc: str,
    input_ref: str,
) -> tuple[dict[str, Any], np.ndarray, float, float]:
    O = build_swap_transport(I=I, Ca=Ca, Cb=Cb, alpha=alpha_mod_pi)
    orth_res = max_abs(O.T @ O - I)
    invol_res = max_abs(O @ O - I)

    obj = {
        "object": object_name,
        "stage": "F466",
        "status": "actual_exported_lane_scoped_pair_chart_transport_operator__axis_only__projector_level__no_false_pass",
        "as_of": AS_OF,
        "generated_utc": generated_utc,
        "inputs": {
            "alpha_mod_pi_source_ref": input_ref,
            "alpha_mod_pi": clean_scalar(alpha_mod_pi),
        },
        "construction": {
            "fourier_pair_planes": {
                f"pair{m_from}": f"V{m_from} := span{{c{m_from},s{m_from}}} (m={m_from}) with canonical real Fourier basis on Z12",
                f"pair{m_to}": f"V{m_to} := span{{c{m_to},s{m_to}}} (m={m_to}) with canonical real Fourier basis on Z12",
            },
            "transport_operator": (
                f"O := C_to G(α) C_from^T + C_from G(-α) C_to^T + Π_rest with α := alpha_mod_pi (axis-only representative); "
                "this transports projectors/spans and is invariant under α -> α+π at projector level."
            ),
        },
        "outputs": {
            "n": int(I.shape[0]),
            "alpha_mod_pi": clean_scalar(alpha_mod_pi),
            "O": [[clean_scalar(float(x)) for x in row] for row in O.tolist()],
            "checks": {
                "orthogonality_max_abs_residual": orth_res,
                "involution_O_squared_equals_I_max_abs_residual": invol_res,
            },
        },
        "hard_limits": [
            "Axis-only: uses alpha mod π; transport is intended for projector/axis gluing, not sign-sensitive orientation.",
            "Lane-scoped: does not export a global selector atlas nor global cocycle data on the full strict domain.",
            "Does not discharge QW-2191 nor imply strict-core selector closure.",
            "Does not claim ToE closure.",
        ],
        "no_false_pass": True,
    }
    out_path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    return obj, O, orth_res, invol_res


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    required = (
        IN_A1,
        IN_A2,
        IN_A3,
        IN_A4,
        IN_A5,
        IN_O12,
        IN_O23,
        IN_O13,
        IN_O34,
        IN_O45,
        IN_O24,
        IN_O35,
    )
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
                        "F464 A_3(pair3) projector operator",
                        "F465 A_4(pair4) and A_5(pair5) projector operators",
                        "F461 O12 chart-transport operator",
                        "F464 O23 and O13 chart-transport operators",
                        "F465 O34/O45/O24/O35 chart-transport operators",
                    ],
                },
                ensure_ascii=True,
            )
        )

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
        u1 = as_vector((a1_obj.get("data") or {}).get("u_1"), n=12, label="F456.data.u_1")
        u2 = as_vector((a2_obj.get("data") or {}).get("u_2"), n=12, label="F462.data.u_2")
        u3 = as_vector((a3_obj.get("data") or {}).get("u_3"), n=12, label="F464.data.u_3")
        u4 = as_vector((a4_obj.get("data") or {}).get("u_4"), n=12, label="F465.data.u_4")
        u5 = as_vector((a5_obj.get("data") or {}).get("u_5"), n=12, label="F465.data.u_5")
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
        out23 = o23_obj.get("outputs") or {}
        O23 = as_square_matrix(out23.get("O23"), n=n, label="F464.outputs.O23")
        alpha23_mod_pi = float(out23.get("alpha23_mod_pi"))
    except Exception as exc:
        raise SystemExit(json.dumps({"status": "INVALID_F464_O23_SHAPE", "error": str(exc)}, ensure_ascii=True))

    try:
        out13 = o13_obj.get("outputs") or {}
        O13 = as_square_matrix(out13.get("O13"), n=n, label="F464.outputs.O13")
        alpha13_mod_pi = float(out13.get("alpha13_mod_pi"))
    except Exception as exc:
        raise SystemExit(json.dumps({"status": "INVALID_F464_O13_SHAPE", "error": str(exc)}, ensure_ascii=True))

    def parse_axis_only_O(obj: dict[str, Any], *, label: str) -> tuple[np.ndarray, float]:
        outputs = obj.get("outputs") or {}
        matrix = outputs.get("O")
        if matrix is None:
            matrix = outputs.get(label)
        alpha = outputs.get("alpha_mod_pi")
        if matrix is None or alpha is None:
            raise ValueError(f"missing outputs.O or outputs.alpha_mod_pi in {label}")
        return as_square_matrix(matrix, n=n, label=f"{label}.outputs.O"), float(alpha)

    try:
        O34, alpha34_mod_pi = parse_axis_only_O(o34_obj, label="F465.outputs.O34")
        O45, alpha45_mod_pi = parse_axis_only_O(o45_obj, label="F465.outputs.O45")
        O24, alpha24_mod_pi = parse_axis_only_O(o24_obj, label="F465.outputs.O24")
        O35, alpha35_mod_pi = parse_axis_only_O(o35_obj, label="F465.outputs.O35")
    except Exception as exc:
        raise SystemExit(json.dumps({"status": "INVALID_F465_AXIS_ONLY_O_SHAPE", "error": str(exc)}, ensure_ascii=True))

    # Canonical Fourier bases for pair1..pair5 (declared scaffold).
    c1, s1 = real_fourier_pair_basis(n=n, m=1)
    c2, s2 = real_fourier_pair_basis(n=n, m=2)
    c3, s3 = real_fourier_pair_basis(n=n, m=3)
    c4, s4 = real_fourier_pair_basis(n=n, m=4)
    c5, s5 = real_fourier_pair_basis(n=n, m=5)

    C1 = np.column_stack([c1, s1])
    C2 = np.column_stack([c2, s2])
    C3 = np.column_stack([c3, s3])
    C4 = np.column_stack([c4, s4])
    C5 = np.column_stack([c5, s5])

    I = np.eye(n, dtype=float)

    # Axis-only angles extracted directly from exported u vectors (mod pi).
    theta1_pi = theta_mod_pi_from_u(u1, c1, s1)
    theta2_pi = theta_mod_pi_from_u(u2, c2, s2)
    theta3_pi = theta_mod_pi_from_u(u3, c3, s3)
    theta4_pi = theta_mod_pi_from_u(u4, c4, s4)
    theta5_pi = theta_mod_pi_from_u(u5, c5, s5)

    alpha12_mod_pi_from_u = angle_mod(theta2_pi - theta1_pi, math.pi)
    alpha12_mod_pi_from_o12 = angle_mod(alpha12_mod_2pi, math.pi)
    alpha12_mod_pi_consistency = circular_distance_mod_pi(alpha12_mod_pi_from_u, alpha12_mod_pi_from_o12)

    alpha23_mod_pi_from_u = angle_mod(theta3_pi - theta2_pi, math.pi)
    alpha13_mod_pi_from_u = angle_mod(theta3_pi - theta1_pi, math.pi)
    alpha23_mod_pi_consistency = circular_distance_mod_pi(alpha23_mod_pi_from_u, angle_mod(alpha23_mod_pi, math.pi))
    alpha13_mod_pi_consistency = circular_distance_mod_pi(alpha13_mod_pi_from_u, angle_mod(alpha13_mod_pi, math.pi))

    # New long-edge angles (axis-only).
    alpha14_mod_pi = angle_mod(theta4_pi - theta1_pi, math.pi)
    alpha15_mod_pi = angle_mod(theta5_pi - theta1_pi, math.pi)
    alpha25_mod_pi = angle_mod(theta5_pi - theta2_pi, math.pi)

    generated_utc = datetime.now(timezone.utc).isoformat()

    # Export additional axis-only transport operators to complete the full triple cocycle set.
    o14_obj, O14, o14_orth_res, o14_invol_res = export_axis_only_transport_operator(
        object_name="O14_pair1_pair4_selector_chart_transport_operator_axis_only_alpha14_mod_pi_strict_core_v1",
        m_from=1,
        m_to=4,
        alpha_mod_pi=alpha14_mod_pi,
        Ca=C1,
        Cb=C4,
        I=I,
        out_path=OUT_O14,
        generated_utc=generated_utc,
        input_ref=f"{IN_A1.relative_to(REPO)} + {IN_A4.relative_to(REPO)}",
    )
    o15_obj, O15, o15_orth_res, o15_invol_res = export_axis_only_transport_operator(
        object_name="O15_pair1_pair5_selector_chart_transport_operator_axis_only_alpha15_mod_pi_strict_core_v1",
        m_from=1,
        m_to=5,
        alpha_mod_pi=alpha15_mod_pi,
        Ca=C1,
        Cb=C5,
        I=I,
        out_path=OUT_O15,
        generated_utc=generated_utc,
        input_ref=f"{IN_A1.relative_to(REPO)} + {IN_A5.relative_to(REPO)}",
    )
    o25_obj, O25, o25_orth_res, o25_invol_res = export_axis_only_transport_operator(
        object_name="O25_pair2_pair5_selector_chart_transport_operator_axis_only_alpha25_mod_pi_strict_core_v1",
        m_from=2,
        m_to=5,
        alpha_mod_pi=alpha25_mod_pi,
        Ca=C2,
        Cb=C5,
        I=I,
        out_path=OUT_O25,
        generated_utc=generated_utc,
        input_ref=f"{IN_A2.relative_to(REPO)} + {IN_A5.relative_to(REPO)}",
    )

    # Projectors / operators on the full carrier.
    A1 = np.outer(u1, u1)
    A2 = np.outer(u2, u2)
    A3 = np.outer(u3, u3)
    A4 = np.outer(u4, u4)
    A5 = np.outer(u5, u5)

    def orth_invol(O: np.ndarray) -> tuple[float, float]:
        return max_abs(O.T @ O - I), max_abs(O @ O - I)

    # Orthogonality/involution audits for reused operators.
    o12_orth_res, o12_invol_res = orth_invol(O12)
    o23_orth_res, o23_invol_res = orth_invol(O23)
    o13_orth_res, o13_invol_res = orth_invol(O13)
    o34_orth_res, o34_invol_res = orth_invol(O34)
    o45_orth_res, o45_invol_res = orth_invol(O45)
    o24_orth_res, o24_invol_res = orth_invol(O24)
    o35_orth_res, o35_invol_res = orth_invol(O35)

    # Gluing residuals (projector level).
    gluing12_res = max_abs(A2 - (O12 @ A1 @ O12.T))
    gluing23_res = max_abs(A3 - (O23 @ A2 @ O23.T))
    gluing13_res = max_abs(A3 - (O13 @ A1 @ O13.T))
    gluing34_res = max_abs(A4 - (O34 @ A3 @ O34.T))
    gluing24_res = max_abs(A4 - (O24 @ A2 @ O24.T))
    gluing14_res = max_abs(A4 - (O14 @ A1 @ O14.T))
    gluing45_res = max_abs(A5 - (O45 @ A4 @ O45.T))
    gluing35_res = max_abs(A5 - (O35 @ A3 @ O35.T))
    gluing25_res = max_abs(A5 - (O25 @ A2 @ O25.T))
    gluing15_res = max_abs(A5 - (O15 @ A1 @ O15.T))

    def cocycle_residual(*, path: np.ndarray, direct: np.ndarray, A: np.ndarray) -> float:
        return max_abs((path @ A @ path.T) - (direct @ A @ direct.T))

    # Full cocycle/path-independence audits on the glued projector section (all triples on {pair1..pair5}).
    cocycle = [
        {
            "triple": "pair1-pair2-pair3",
            "statement": "O_23 O_12 transports A_1 to the same A_3 as the direct transport O_13 (projector section).",
            "path": "O_23 O_12",
            "direct": "O_13",
            "max_abs_residual": cocycle_residual(path=(O23 @ O12), direct=O13, A=A1),
        },
        {
            "triple": "pair1-pair2-pair4",
            "statement": "O_24 O_12 transports A_1 to the same A_4 as the direct transport O_14 (projector section).",
            "path": "O_24 O_12",
            "direct": "O_14",
            "max_abs_residual": cocycle_residual(path=(O24 @ O12), direct=O14, A=A1),
        },
        {
            "triple": "pair1-pair2-pair5",
            "statement": "O_25 O_12 transports A_1 to the same A_5 as the direct transport O_15 (projector section).",
            "path": "O_25 O_12",
            "direct": "O_15",
            "max_abs_residual": cocycle_residual(path=(O25 @ O12), direct=O15, A=A1),
        },
        {
            "triple": "pair1-pair3-pair4",
            "statement": "O_34 O_13 transports A_1 to the same A_4 as the direct transport O_14 (projector section).",
            "path": "O_34 O_13",
            "direct": "O_14",
            "max_abs_residual": cocycle_residual(path=(O34 @ O13), direct=O14, A=A1),
        },
        {
            "triple": "pair1-pair3-pair5",
            "statement": "O_35 O_13 transports A_1 to the same A_5 as the direct transport O_15 (projector section).",
            "path": "O_35 O_13",
            "direct": "O_15",
            "max_abs_residual": cocycle_residual(path=(O35 @ O13), direct=O15, A=A1),
        },
        {
            "triple": "pair1-pair4-pair5",
            "statement": "O_45 O_14 transports A_1 to the same A_5 as the direct transport O_15 (projector section).",
            "path": "O_45 O_14",
            "direct": "O_15",
            "max_abs_residual": cocycle_residual(path=(O45 @ O14), direct=O15, A=A1),
        },
        {
            "triple": "pair2-pair3-pair4",
            "statement": "O_34 O_23 transports A_2 to the same A_4 as the direct transport O_24 (projector section).",
            "path": "O_34 O_23",
            "direct": "O_24",
            "max_abs_residual": cocycle_residual(path=(O34 @ O23), direct=O24, A=A2),
        },
        {
            "triple": "pair2-pair3-pair5",
            "statement": "O_35 O_23 transports A_2 to the same A_5 as the direct transport O_25 (projector section).",
            "path": "O_35 O_23",
            "direct": "O_25",
            "max_abs_residual": cocycle_residual(path=(O35 @ O23), direct=O25, A=A2),
        },
        {
            "triple": "pair2-pair4-pair5",
            "statement": "O_45 O_24 transports A_2 to the same A_5 as the direct transport O_25 (projector section).",
            "path": "O_45 O_24",
            "direct": "O_25",
            "max_abs_residual": cocycle_residual(path=(O45 @ O24), direct=O25, A=A2),
        },
        {
            "triple": "pair3-pair4-pair5",
            "statement": "O_45 O_34 transports A_3 to the same A_5 as the direct transport O_35 (projector section).",
            "path": "O_45 O_34",
            "direct": "O_35",
            "max_abs_residual": cocycle_residual(path=(O45 @ O34), direct=O35, A=A3),
        },
    ]

    # Residual sign gauge invariance at projector level.
    a1_sign_res = max_abs(A1 - np.outer(-u1, -u1))
    a2_sign_res = max_abs(A2 - np.outer(-u2, -u2))
    a3_sign_res = max_abs(A3 - np.outer(-u3, -u3))
    a4_sign_res = max_abs(A4 - np.outer(-u4, -u4))
    a5_sign_res = max_abs(A5 - np.outer(-u5, -u5))

    section_obj = {
        "object": "A_12345_pair12345_chart_glued_orientation_projector_operator_section_strict_core_v2",
        "stage": "F466",
        "status": "actual_exported_five_chart_projector_operator_section__pair1_to_pair5__with_full_cocycle_audits__no_false_pass",
        "as_of": AS_OF,
        "generated_utc": generated_utc,
        "intent": (
            "Export one explicit five-chart projector operator section on {pair1..pair5} on the declared n=12 Fourier carrier, "
            "upgrading the previous local-cocycle stage by exporting additional axis-only long-edge chart transport operators and "
            "recording explicit cocycle/path-independence audits for all triple overlaps on {pair1..pair5} at the level of the glued projector section. "
            "This remains lane-scoped and does not imply a global selector atlas nor any global discharge."
        ),
        "inputs": {
            "a1_pair1_operator_ref": str(IN_A1.relative_to(REPO)),
            "a2_pair2_operator_ref": str(IN_A2.relative_to(REPO)),
            "a3_pair3_operator_ref": str(IN_A3.relative_to(REPO)),
            "a4_pair4_operator_ref": str(IN_A4.relative_to(REPO)),
            "a5_pair5_operator_ref": str(IN_A5.relative_to(REPO)),
            "o12_ref": str(IN_O12.relative_to(REPO)),
            "o23_ref": str(IN_O23.relative_to(REPO)),
            "o13_ref": str(IN_O13.relative_to(REPO)),
            "o34_ref": str(IN_O34.relative_to(REPO)),
            "o45_ref": str(IN_O45.relative_to(REPO)),
            "o24_ref": str(IN_O24.relative_to(REPO)),
            "o35_ref": str(IN_O35.relative_to(REPO)),
            "o14_ref": str(OUT_O14.relative_to(REPO)),
            "o15_ref": str(OUT_O15.relative_to(REPO)),
            "o25_ref": str(OUT_O25.relative_to(REPO)),
        },
        "gluing_laws": [
            "A_2(pair2) = O_12 A_1(pair1) O_12^T  (projector-level)",
            "A_3(pair3) = O_23 A_2(pair2) O_23^T  (projector-level, axis-only transport)",
            "A_3(pair3) = O_13 A_1(pair1) O_13^T  (projector-level, axis-only transport)",
            "A_4(pair4) = O_34 A_3(pair3) O_34^T  (projector-level, axis-only transport)",
            "A_4(pair4) = O_24 A_2(pair2) O_24^T  (projector-level, axis-only transport)",
            "A_4(pair4) = O_14 A_1(pair1) O_14^T  (projector-level, axis-only transport)",
            "A_5(pair5) = O_45 A_4(pair4) O_45^T  (projector-level, axis-only transport)",
            "A_5(pair5) = O_35 A_3(pair3) O_35^T  (projector-level, axis-only transport)",
            "A_5(pair5) = O_25 A_2(pair2) O_25^T  (projector-level, axis-only transport)",
            "A_5(pair5) = O_15 A_1(pair1) O_15^T  (projector-level, axis-only transport)",
        ],
        "cocycle_audits": cocycle,
        "audits": {
            "alpha_mod_pi_consistency": {
                "alpha12_mod_pi_distance_u_vs_o12": alpha12_mod_pi_consistency,
                "alpha23_mod_pi_distance_u_vs_o23": alpha23_mod_pi_consistency,
                "alpha13_mod_pi_distance_u_vs_o13": alpha13_mod_pi_consistency,
            },
            "orthogonality": {
                "o12_max_abs_residual": o12_orth_res,
                "o23_max_abs_residual": o23_orth_res,
                "o13_max_abs_residual": o13_orth_res,
                "o34_max_abs_residual": o34_orth_res,
                "o45_max_abs_residual": o45_orth_res,
                "o24_max_abs_residual": o24_orth_res,
                "o35_max_abs_residual": o35_orth_res,
                "o14_max_abs_residual": o14_orth_res,
                "o15_max_abs_residual": o15_orth_res,
                "o25_max_abs_residual": o25_orth_res,
            },
            "involution": {
                "o12_max_abs_residual": o12_invol_res,
                "o23_max_abs_residual": o23_invol_res,
                "o13_max_abs_residual": o13_invol_res,
                "o34_max_abs_residual": o34_invol_res,
                "o45_max_abs_residual": o45_invol_res,
                "o24_max_abs_residual": o24_invol_res,
                "o35_max_abs_residual": o35_invol_res,
                "o14_max_abs_residual": o14_invol_res,
                "o15_max_abs_residual": o15_invol_res,
                "o25_max_abs_residual": o25_invol_res,
            },
            "gluing": {
                "gluing12_max_abs_residual": gluing12_res,
                "gluing23_max_abs_residual": gluing23_res,
                "gluing13_max_abs_residual": gluing13_res,
                "gluing34_max_abs_residual": gluing34_res,
                "gluing24_max_abs_residual": gluing24_res,
                "gluing14_max_abs_residual": gluing14_res,
                "gluing45_max_abs_residual": gluing45_res,
                "gluing35_max_abs_residual": gluing35_res,
                "gluing25_max_abs_residual": gluing25_res,
                "gluing15_max_abs_residual": gluing15_res,
            },
            "sign_gauge": {
                "a1_projector_sign_invariance_max_abs": a1_sign_res,
                "a2_projector_sign_invariance_max_abs": a2_sign_res,
                "a3_projector_sign_invariance_max_abs": a3_sign_res,
                "a4_projector_sign_invariance_max_abs": a4_sign_res,
                "a5_projector_sign_invariance_max_abs": a5_sign_res,
            },
        },
        "hard_limits": [
            "Projector-level only: section gluing is sign-gauge-safe but does not lift residual sign to a physical convention.",
            "Lane-scoped five-chart section only; does not export a global selector atlas nor global cocycle data on the full strict domain.",
            "Does not discharge QW-2191 nor export strict-core selector closure / admissible S_sel_int.",
            "Does not claim ToE closure.",
        ],
        "no_false_pass": True,
    }

    atlas = {
        "object": "SelectorAtlas_pair12345_axis_only_projector_v2",
        "stage": "F466",
        "status": "actual_exported_lane_scoped_five_chart_selector_atlas_with_projector_level_gluing_and_full_cocycle_audits__no_false_pass",
        "as_of": AS_OF,
        "generated_utc": generated_utc,
        "intent": (
            "Export one explicit lane-scoped five-chart selector atlas object on {pair1..pair5} at projector level (sign-gauge-safe), "
            "upgrading the previous local-cocycle stage by exporting additional long-edge chart transport operators and explicit cocycle/path-independence "
            "audits for all triple overlaps on {pair1..pair5} (projector section). This is not a global selector atlas and does not imply any global discharge."
        ),
        "charts": {
            "pair1": {
                "chart_id": "pair1",
                "carrier_plane": "V1 = span{c1,s1} inside the n=12 real Fourier scaffold",
                "local_operator": {"ref": str(IN_A1.relative_to(REPO)), "object": a1_obj.get("object")},
            },
            "pair2": {
                "chart_id": "pair2",
                "carrier_plane": "V2 = span{c2,s2} inside the n=12 real Fourier scaffold",
                "local_operator": {"ref": str(IN_A2.relative_to(REPO)), "object": a2_obj.get("object")},
            },
            "pair3": {
                "chart_id": "pair3",
                "carrier_plane": "V3 = span{c3,s3} inside the n=12 real Fourier scaffold",
                "local_operator": {"ref": str(IN_A3.relative_to(REPO)), "object": a3_obj.get("object")},
            },
            "pair4": {
                "chart_id": "pair4",
                "carrier_plane": "V4 = span{c4,s4} inside the n=12 real Fourier scaffold",
                "local_operator": {"ref": str(IN_A4.relative_to(REPO)), "object": a4_obj.get("object")},
            },
            "pair5": {
                "chart_id": "pair5",
                "carrier_plane": "V5 = span{c5,s5} inside the n=12 real Fourier scaffold",
                "local_operator": {"ref": str(IN_A5.relative_to(REPO)), "object": a5_obj.get("object")},
            },
        },
        "overlap_domain_declaration": {
            "scope": "exported_artifact_overlap_on_n12_fourier_carrier",
            "meaning": (
                "This is an exported-artifact overlap declaration: all five chart operators A_1..A_5 and the chart transport operators "
                "are simultaneously defined as strict-core exported artifacts on the declared n=12 Fourier carrier. "
                "This is not a global open-cover atlas on the full strict domain C_v1."
            ),
        },
        "transitions": {
            "pair1_to_pair2": {
                "operator_ref": str(IN_O12.relative_to(REPO)),
                "operator": o12_obj.get("object"),
                "alpha12_mod_2pi": clean_scalar(alpha12_mod_2pi),
            },
            "pair2_to_pair3": {"operator_ref": str(IN_O23.relative_to(REPO)), "operator": o23_obj.get("object"), "alpha23_mod_pi": clean_scalar(alpha23_mod_pi)},
            "pair1_to_pair3": {"operator_ref": str(IN_O13.relative_to(REPO)), "operator": o13_obj.get("object"), "alpha13_mod_pi": clean_scalar(alpha13_mod_pi)},
            "pair3_to_pair4": {"operator_ref": str(IN_O34.relative_to(REPO)), "operator": o34_obj.get("object"), "alpha34_mod_pi": clean_scalar(alpha34_mod_pi)},
            "pair4_to_pair5": {"operator_ref": str(IN_O45.relative_to(REPO)), "operator": o45_obj.get("object"), "alpha45_mod_pi": clean_scalar(alpha45_mod_pi)},
            "pair2_to_pair4": {"operator_ref": str(IN_O24.relative_to(REPO)), "operator": o24_obj.get("object"), "alpha24_mod_pi": clean_scalar(alpha24_mod_pi)},
            "pair3_to_pair5": {"operator_ref": str(IN_O35.relative_to(REPO)), "operator": o35_obj.get("object"), "alpha35_mod_pi": clean_scalar(alpha35_mod_pi)},
            "pair1_to_pair4": {"operator_ref": str(OUT_O14.relative_to(REPO)), "operator": o14_obj.get("object"), "alpha14_mod_pi": clean_scalar(alpha14_mod_pi)},
            "pair1_to_pair5": {"operator_ref": str(OUT_O15.relative_to(REPO)), "operator": o15_obj.get("object"), "alpha15_mod_pi": clean_scalar(alpha15_mod_pi)},
            "pair2_to_pair5": {"operator_ref": str(OUT_O25.relative_to(REPO)), "operator": o25_obj.get("object"), "alpha25_mod_pi": clean_scalar(alpha25_mod_pi)},
        },
        "gluing_data": {
            "operator_section_ref": str(OUT_SECTION.relative_to(REPO)),
            "laws": section_obj["gluing_laws"],
            "cocycle_audits": section_obj["cocycle_audits"],
        },
        "audits": section_obj["audits"],
        "hard_limits": [
            "Lane-scoped: exports only a five-chart atlas object on {pair1..pair5} at projector level.",
            "Does not export a global selector atlas, global overlap-domain declaration, nor global cocycle data (H41 remains open globally).",
            "Does not discharge QW-2191 nor export a global selector transition/gluing object (H40 remains open globally).",
            "Projector-level only: does not derive a sign-sensitive physical orientation datum.",
            "Does not claim strict-core selector closure / admissible S_sel_int.",
            "Does not claim ToE closure.",
        ],
        "no_false_pass": True,
    }

    cocycle_max = max((float(item["max_abs_residual"]) for item in cocycle), default=0.0)
    gluing_max = float(max(gluing12_res, gluing23_res, gluing13_res, gluing34_res, gluing24_res, gluing14_res, gluing45_res, gluing35_res, gluing25_res, gluing15_res))

    summary = {
        "stage": "F466",
        "status": "F466_EXECUTED_CURRENT_STRICT_PAIR12345_SELECTOR_ATLAS_FULL_COCYCLE_DATA_EXPORT_PACKET_NO_FALSE_PASS",
        "generated_utc": generated_utc,
        "artifacts": {
            "o14": str(OUT_O14.relative_to(REPO)),
            "o15": str(OUT_O15.relative_to(REPO)),
            "o25": str(OUT_O25.relative_to(REPO)),
            "section_v2": str(OUT_SECTION.relative_to(REPO)),
            "atlas_v2": str(OUT_ATLAS.relative_to(REPO)),
        },
        "audits": {
            "alpha12_mod_pi_consistency_distance": alpha12_mod_pi_consistency,
            "gluing_max_abs_residual": gluing_max,
            "cocycle_max_abs_residual": cocycle_max,
        },
        "no_false_pass": True,
    }

    OUT_SECTION.write_text(json.dumps(section_obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_ATLAS.write_text(json.dumps(atlas, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
