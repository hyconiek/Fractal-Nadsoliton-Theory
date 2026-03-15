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

IN_A3 = GENERATED / "a_3_pair3_orientation_projector_operator_strict_core_v1.json"
IN_O23 = GENERATED / "o23_pair2_pair3_selector_chart_transport_operator_axis_only_alpha23_mod_pi_strict_core_v1.json"
IN_O13 = GENERATED / "o13_pair1_pair3_selector_chart_transport_operator_axis_only_alpha13_mod_pi_strict_core_v1.json"

IN_MODE_ASSIGN = GENERATED / "mode_index_assignment_shannon_element_order_reference_strict_core_v1.json"

OUT_A4 = GENERATED / "a_4_pair4_orientation_projector_operator_strict_core_v1.json"
OUT_A5 = GENERATED / "a_5_pair5_orientation_projector_operator_strict_core_v1.json"

OUT_O34 = GENERATED / "o34_pair3_pair4_selector_chart_transport_operator_axis_only_alpha34_mod_pi_strict_core_v1.json"
OUT_O45 = GENERATED / "o45_pair4_pair5_selector_chart_transport_operator_axis_only_alpha45_mod_pi_strict_core_v1.json"
OUT_O24 = GENERATED / "o24_pair2_pair4_selector_chart_transport_operator_axis_only_alpha24_mod_pi_strict_core_v1.json"
OUT_O35 = GENERATED / "o35_pair3_pair5_selector_chart_transport_operator_axis_only_alpha35_mod_pi_strict_core_v1.json"

OUT_SECTION = GENERATED / "a_12345_pair12345_chart_glued_orientation_projector_operator_section_strict_core_v1.json"
OUT_ATLAS = GENERATED / "selector_atlas_pair12345_axis_only_projector_v1.json"
OUT_SUMMARY = GENERATED / "f465_current_strict_pair12345_selector_atlas_cocycle_data_export_packet_summary.json"


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


def theta_mod_pi_from_u(u: np.ndarray, c: np.ndarray, s: np.ndarray) -> tuple[float, np.ndarray, float]:
    coords = np.array([float(c @ u), float(s @ u)], dtype=float)
    norm = float(np.linalg.norm(coords))
    if norm == 0.0:
        raise ValueError("degenerate coordinates (zero)")
    theta = math.atan2(float(coords[1]), float(coords[0]))
    return angle_mod(theta, math.pi), coords / norm, norm


def export_pair_projector_operator(
    *,
    pair_label: str,
    pair_m: int,
    u: np.ndarray,
    c: np.ndarray,
    s: np.ndarray,
    input_ref: str,
    out_path: Path,
) -> dict[str, Any]:
    theta_pi, coords_unit, coords_norm = theta_mod_pi_from_u(u, c, s)
    u_recon = (float(coords_unit[0]) * c) + (float(coords_unit[1]) * s)
    u_outside = float(np.linalg.norm(u - u_recon))

    A_full = np.outer(u, u)
    A_2 = np.outer(coords_unit, coords_unit)

    sym_res = max_abs(A_full - A_full.T)
    idem_res = max_abs(A_full @ A_full - A_full)
    sign_res = max_abs(A_full - np.outer(-u, -u))

    obj = {
        "object": f"A_{pair_m}_{pair_label}_orientation_projector_operator_strict_core_v1",
        "status": "actual_exported_operator__derived_from_axis_only_pair_plane_orientation_vector__residual_z2_sign_gauge_invariant",
        "as_of": AS_OF,
        "intent": (
            f"Export one strict-core projector operator on V_{pair_m} = span{{c{pair_m},s{pair_m}}} constructed from the exported axis-only "
            f"pair-plane orientation vector u_{pair_m}. The operator is A_{pair_m}({pair_label}) := |u_{pair_m}><u_{pair_m}|. "
            "This is projector-level (residual-sign-gauge-invariant) and does not imply any sign-sensitive physical orientation convention."
        ),
        "inputs": {
            "axis_only_source_ref": input_ref,
        },
        "domain_notes": {
            "basis_plane": f"V_{pair_m} = span{{c{pair_m},s{pair_m}}} inside the n=12 real Fourier scaffold",
            "lane": "axis_only_projector_level",
        },
        "data": {
            f"theta_{pair_m}_mod_pi": clean_scalar(theta_pi),
            f"u_{pair_m}": [clean_scalar(float(x)) for x in u.tolist()],
            f"u_{pair_m}_coords_in_c{pair_m}_s{pair_m}": [clean_scalar(float(x)) for x in coords_unit.tolist()],
            f"A_{pair_m}_{pair_label}_matrix_in_c{pair_m}_s{pair_m}": [[clean_scalar(float(x)) for x in row] for row in A_2.tolist()],
        },
        "audits": {
            f"u{pair_m}_outside_span_c{pair_m}_s{pair_m}_l2": u_outside,
            "coords_norm": coords_norm,
            "projector_symmetry_inf_norm": sym_res,
            "projector_idempotence_inf_norm": idem_res,
            "projector_trace": float(np.trace(A_full)),
            "residual_z2_sign_invariance_max_abs_diff": sign_res,
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
    out_path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    return obj


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
        "stage": "F465",
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

    required = (IN_A1, IN_A2, IN_O12, IN_A3, IN_O23, IN_O13, IN_MODE_ASSIGN)
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
                        "F464 A_3(pair3) projector operator",
                        "F464 O23 and O13 chart-transport operators",
                        "F454 Shannon mode-index assignment (pair4/pair5 minimizer axes)",
                    ],
                },
                ensure_ascii=True,
            )
        )

    a1_obj = load_json(IN_A1)
    a2_obj = load_json(IN_A2)
    o12_obj = load_json(IN_O12)
    a3_obj = load_json(IN_A3)
    o23_obj = load_json(IN_O23)
    o13_obj = load_json(IN_O13)
    mode_assign = load_json(IN_MODE_ASSIGN)

    try:
        u1 = as_vector((a1_obj.get("data") or {}).get("u_1"), n=12, label="F456.data.u_1")
        u2 = as_vector((a2_obj.get("data") or {}).get("u_2"), n=12, label="F462.data.u_2")
        u3 = as_vector((a3_obj.get("data") or {}).get("u_3"), n=12, label="F464.data.u_3")
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
        O23 = as_square_matrix((o23_obj.get("outputs") or {}).get("O23"), n=12, label="F464.outputs.O23")
        alpha23_mod_pi = float((o23_obj.get("outputs") or {}).get("alpha23_mod_pi"))
    except Exception as exc:
        raise SystemExit(json.dumps({"status": "INVALID_F464_O23_SHAPE", "error": str(exc)}, ensure_ascii=True))

    try:
        O13 = as_square_matrix((o13_obj.get("outputs") or {}).get("O13"), n=12, label="F464.outputs.O13")
        alpha13_mod_pi = float((o13_obj.get("outputs") or {}).get("alpha13_mod_pi"))
    except Exception as exc:
        raise SystemExit(json.dumps({"status": "INVALID_F464_O13_SHAPE", "error": str(exc)}, ensure_ascii=True))

    try:
        pairs = ((mode_assign.get("outputs") or {}).get("pairs") or {})
        u4 = as_vector((pairs.get("pair4") or {}).get("u_minus"), n=12, label="F454.outputs.pairs.pair4.u_minus")
        u5 = as_vector((pairs.get("pair5") or {}).get("u_minus"), n=12, label="F454.outputs.pairs.pair5.u_minus")
    except Exception as exc:
        raise SystemExit(json.dumps({"status": "INVALID_F454_SHAPE", "error": str(exc)}, ensure_ascii=True))

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

    # Axis-only angles extracted directly from the exported u vectors (mod pi).
    theta1_pi, _, _ = theta_mod_pi_from_u(u1, c1, s1)
    theta2_pi, _, _ = theta_mod_pi_from_u(u2, c2, s2)
    theta3_pi, _, _ = theta_mod_pi_from_u(u3, c3, s3)
    theta4_pi, _, _ = theta_mod_pi_from_u(u4, c4, s4)
    theta5_pi, _, _ = theta_mod_pi_from_u(u5, c5, s5)

    alpha12_mod_pi_from_u = angle_mod(theta2_pi - theta1_pi, math.pi)
    alpha12_mod_pi_from_o12 = angle_mod(alpha12_mod_2pi, math.pi)
    alpha12_mod_pi_consistency = circular_distance_mod_pi(alpha12_mod_pi_from_u, alpha12_mod_pi_from_o12)

    alpha23_mod_pi_from_u = angle_mod(theta3_pi - theta2_pi, math.pi)
    alpha13_mod_pi_from_u = angle_mod(theta3_pi - theta1_pi, math.pi)
    alpha23_mod_pi_consistency = circular_distance_mod_pi(alpha23_mod_pi_from_u, angle_mod(alpha23_mod_pi, math.pi))
    alpha13_mod_pi_consistency = circular_distance_mod_pi(alpha13_mod_pi_from_u, angle_mod(alpha13_mod_pi, math.pi))

    alpha34_mod_pi = angle_mod(theta4_pi - theta3_pi, math.pi)
    alpha45_mod_pi = angle_mod(theta5_pi - theta4_pi, math.pi)
    alpha24_mod_pi = angle_mod(theta4_pi - theta2_pi, math.pi)
    alpha35_mod_pi = angle_mod(theta5_pi - theta3_pi, math.pi)

    generated_utc = datetime.now(timezone.utc).isoformat()

    # Export A4/A5 projector operators.
    a4_obj = export_pair_projector_operator(
        pair_label="pair4",
        pair_m=4,
        u=u4,
        c=c4,
        s=s4,
        input_ref=str(IN_MODE_ASSIGN.relative_to(REPO)),
        out_path=OUT_A4,
    )
    a5_obj = export_pair_projector_operator(
        pair_label="pair5",
        pair_m=5,
        u=u5,
        c=c5,
        s=s5,
        input_ref=str(IN_MODE_ASSIGN.relative_to(REPO)),
        out_path=OUT_A5,
    )

    # Export additional axis-only transport operators.
    o34_obj, O34, o34_orth_res, o34_invol_res = export_axis_only_transport_operator(
        object_name="O34_pair3_pair4_selector_chart_transport_operator_axis_only_alpha34_mod_pi_strict_core_v1",
        m_from=3,
        m_to=4,
        alpha_mod_pi=alpha34_mod_pi,
        Ca=C3,
        Cb=C4,
        I=I,
        out_path=OUT_O34,
        generated_utc=generated_utc,
        input_ref=str(IN_MODE_ASSIGN.relative_to(REPO)),
    )
    o45_obj, O45, o45_orth_res, o45_invol_res = export_axis_only_transport_operator(
        object_name="O45_pair4_pair5_selector_chart_transport_operator_axis_only_alpha45_mod_pi_strict_core_v1",
        m_from=4,
        m_to=5,
        alpha_mod_pi=alpha45_mod_pi,
        Ca=C4,
        Cb=C5,
        I=I,
        out_path=OUT_O45,
        generated_utc=generated_utc,
        input_ref=str(IN_MODE_ASSIGN.relative_to(REPO)),
    )
    o24_obj, O24, o24_orth_res, o24_invol_res = export_axis_only_transport_operator(
        object_name="O24_pair2_pair4_selector_chart_transport_operator_axis_only_alpha24_mod_pi_strict_core_v1",
        m_from=2,
        m_to=4,
        alpha_mod_pi=alpha24_mod_pi,
        Ca=C2,
        Cb=C4,
        I=I,
        out_path=OUT_O24,
        generated_utc=generated_utc,
        input_ref=str(IN_MODE_ASSIGN.relative_to(REPO)),
    )
    o35_obj, O35, o35_orth_res, o35_invol_res = export_axis_only_transport_operator(
        object_name="O35_pair3_pair5_selector_chart_transport_operator_axis_only_alpha35_mod_pi_strict_core_v1",
        m_from=3,
        m_to=5,
        alpha_mod_pi=alpha35_mod_pi,
        Ca=C3,
        Cb=C5,
        I=I,
        out_path=OUT_O35,
        generated_utc=generated_utc,
        input_ref=str(IN_MODE_ASSIGN.relative_to(REPO)),
    )

    # Projectors / operators on the full carrier.
    A1 = np.outer(u1, u1)
    A2 = np.outer(u2, u2)
    A3 = np.outer(u3, u3)
    A4 = np.outer(u4, u4)
    A5 = np.outer(u5, u5)

    # Orthogonality/involution audits for reused operators.
    orth12_res = max_abs(O12.T @ O12 - I)
    orth23_res = max_abs(O23.T @ O23 - I)
    orth13_res = max_abs(O13.T @ O13 - I)
    invol12_res = max_abs(O12 @ O12 - I)
    invol23_res = max_abs(O23 @ O23 - I)
    invol13_res = max_abs(O13 @ O13 - I)

    # Gluing residuals (projector level).
    gluing12_res = max_abs(A2 - (O12 @ A1 @ O12.T))
    gluing23_res = max_abs(A3 - (O23 @ A2 @ O23.T))
    gluing13_res = max_abs(A3 - (O13 @ A1 @ O13.T))
    gluing34_res = max_abs(A4 - (O34 @ A3 @ O34.T))
    gluing45_res = max_abs(A5 - (O45 @ A4 @ O45.T))
    gluing24_res = max_abs(A4 - (O24 @ A2 @ O24.T))
    gluing35_res = max_abs(A5 - (O35 @ A3 @ O35.T))

    # Cocycle/path-independence audits on the glued projector section.
    O23O12 = O23 @ O12
    cocycle_1_to_3_res = max_abs((O23O12 @ A1 @ O23O12.T) - (O13 @ A1 @ O13.T))

    O34O23 = O34 @ O23
    cocycle_2_to_4_res = max_abs((O34O23 @ A2 @ O34O23.T) - (O24 @ A2 @ O24.T))

    O45O34 = O45 @ O34
    cocycle_3_to_5_res = max_abs((O45O34 @ A3 @ O45O34.T) - (O35 @ A3 @ O35.T))

    # Residual sign gauge invariance at projector level.
    a1_sign_res = max_abs(A1 - np.outer(-u1, -u1))
    a2_sign_res = max_abs(A2 - np.outer(-u2, -u2))
    a3_sign_res = max_abs(A3 - np.outer(-u3, -u3))
    a4_sign_res = max_abs(A4 - np.outer(-u4, -u4))
    a5_sign_res = max_abs(A5 - np.outer(-u5, -u5))

    section_obj = {
        "object": "A_12345_pair12345_chart_glued_orientation_projector_operator_section_strict_core_v1",
        "stage": "F465",
        "status": "actual_exported_five_chart_projector_operator_section__pair1_to_pair5__with_local_cocycle_audits__no_false_pass",
        "as_of": AS_OF,
        "generated_utc": generated_utc,
        "intent": (
            "Export one explicit five-chart projector operator section on {pair1..pair5} on the declared n=12 Fourier carrier, "
            "with explicit chart transport operators and projector-level gluing laws, including explicit local cocycle/path-independence "
            "audits for adjacent triple overlaps (1-2-3, 2-3-4, 3-4-5) at the level of the glued projector section. "
            "This is lane-scoped and does not imply a global selector atlas nor any global discharge."
        ),
        "inputs": {
            "a1_pair1_operator_ref": str(IN_A1.relative_to(REPO)),
            "a2_pair2_operator_ref": str(IN_A2.relative_to(REPO)),
            "a3_pair3_operator_ref": str(IN_A3.relative_to(REPO)),
            "a4_pair4_operator_ref": str(OUT_A4.relative_to(REPO)),
            "a5_pair5_operator_ref": str(OUT_A5.relative_to(REPO)),
            "o12_ref": str(IN_O12.relative_to(REPO)),
            "o23_ref": str(IN_O23.relative_to(REPO)),
            "o13_ref": str(IN_O13.relative_to(REPO)),
            "o34_ref": str(OUT_O34.relative_to(REPO)),
            "o45_ref": str(OUT_O45.relative_to(REPO)),
            "o24_ref": str(OUT_O24.relative_to(REPO)),
            "o35_ref": str(OUT_O35.relative_to(REPO)),
            "mode_index_assignment_ref": str(IN_MODE_ASSIGN.relative_to(REPO)),
        },
        "gluing_laws": [
            "A_2(pair2) = O_12 A_1(pair1) O_12^T  (projector-level)",
            "A_3(pair3) = O_23 A_2(pair2) O_23^T  (projector-level, axis-only transport)",
            "A_4(pair4) = O_34 A_3(pair3) O_34^T  (projector-level, axis-only transport)",
            "A_5(pair5) = O_45 A_4(pair4) O_45^T  (projector-level, axis-only transport)",
            "A_3(pair3) = O_13 A_1(pair1) O_13^T  (projector-level, axis-only transport)",
            "A_4(pair4) = O_24 A_2(pair2) O_24^T  (projector-level, axis-only transport)",
            "A_5(pair5) = O_35 A_3(pair3) O_35^T  (projector-level, axis-only transport)",
        ],
        "local_cocycle_audits": [
            {
                "triple": "pair1-pair2-pair3",
                "statement": "O_23 O_12 transports A_1 to the same A_3 as the direct transport O_13 (projector section).",
                "path": "O_23 O_12",
                "direct": "O_13",
                "max_abs_residual": cocycle_1_to_3_res,
            },
            {
                "triple": "pair2-pair3-pair4",
                "statement": "O_34 O_23 transports A_2 to the same A_4 as the direct transport O_24 (projector section).",
                "path": "O_34 O_23",
                "direct": "O_24",
                "max_abs_residual": cocycle_2_to_4_res,
            },
            {
                "triple": "pair3-pair4-pair5",
                "statement": "O_45 O_34 transports A_3 to the same A_5 as the direct transport O_35 (projector section).",
                "path": "O_45 O_34",
                "direct": "O_35",
                "max_abs_residual": cocycle_3_to_5_res,
            },
        ],
        "audits": {
            "alpha_mod_pi_consistency": {
                "alpha12_mod_pi_distance_u_vs_o12": alpha12_mod_pi_consistency,
                "alpha23_mod_pi_distance_u_vs_o23": alpha23_mod_pi_consistency,
                "alpha13_mod_pi_distance_u_vs_o13": alpha13_mod_pi_consistency,
            },
            "orthogonality": {
                "o12_max_abs_residual": orth12_res,
                "o23_max_abs_residual": orth23_res,
                "o13_max_abs_residual": orth13_res,
                "o34_max_abs_residual": o34_orth_res,
                "o45_max_abs_residual": o45_orth_res,
                "o24_max_abs_residual": o24_orth_res,
                "o35_max_abs_residual": o35_orth_res,
            },
            "involution": {
                "o12_max_abs_residual": invol12_res,
                "o23_max_abs_residual": invol23_res,
                "o13_max_abs_residual": invol13_res,
                "o34_max_abs_residual": o34_invol_res,
                "o45_max_abs_residual": o45_invol_res,
                "o24_max_abs_residual": o24_invol_res,
                "o35_max_abs_residual": o35_invol_res,
            },
            "gluing": {
                "gluing12_max_abs_residual": gluing12_res,
                "gluing23_max_abs_residual": gluing23_res,
                "gluing13_max_abs_residual": gluing13_res,
                "gluing34_max_abs_residual": gluing34_res,
                "gluing45_max_abs_residual": gluing45_res,
                "gluing24_max_abs_residual": gluing24_res,
                "gluing35_max_abs_residual": gluing35_res,
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
        "object": "SelectorAtlas_pair12345_axis_only_projector_v1",
        "stage": "F465",
        "status": "actual_exported_lane_scoped_five_chart_selector_atlas_with_projector_level_gluing_and_local_cocycle_audits__no_false_pass",
        "as_of": AS_OF,
        "generated_utc": generated_utc,
        "intent": (
            "Export one explicit lane-scoped five-chart selector atlas object on {pair1..pair5} at projector level (sign-gauge-safe), "
            "including explicit chart transport operators, projector-level gluing laws, and explicit local cocycle/path-independence audits "
            "for adjacent triple overlaps on the glued projector section. This is not a global selector atlas and does not imply any global discharge."
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
                    "ref": str(IN_A3.relative_to(REPO)),
                    "object": "A_3_pair3_orientation_projector_operator_strict_core_v1",
                },
            },
            "pair4": {
                "chart_id": "pair4",
                "carrier_plane": "V4 = span{c4,s4} inside the n=12 real Fourier scaffold",
                "local_operator": {
                    "ref": str(OUT_A4.relative_to(REPO)),
                    "object": "A_4_pair4_orientation_projector_operator_strict_core_v1",
                },
            },
            "pair5": {
                "chart_id": "pair5",
                "carrier_plane": "V5 = span{c5,s5} inside the n=12 real Fourier scaffold",
                "local_operator": {
                    "ref": str(OUT_A5.relative_to(REPO)),
                    "object": "A_5_pair5_orientation_projector_operator_strict_core_v1",
                },
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
                "operator": "O12_pair1_pair2_selector_chart_transport_operator_strict_derived_from_sigma_int_alpha12_v1",
                "alpha12_mod_2pi": clean_scalar(alpha12_mod_2pi),
            },
            "pair2_to_pair3": {"operator_ref": str(IN_O23.relative_to(REPO)), "operator": o23_obj.get("object"), "alpha23_mod_pi": clean_scalar(alpha23_mod_pi)},
            "pair1_to_pair3": {"operator_ref": str(IN_O13.relative_to(REPO)), "operator": o13_obj.get("object"), "alpha13_mod_pi": clean_scalar(alpha13_mod_pi)},
            "pair3_to_pair4": {"operator_ref": str(OUT_O34.relative_to(REPO)), "operator": o34_obj.get("object"), "alpha34_mod_pi": clean_scalar(alpha34_mod_pi)},
            "pair4_to_pair5": {"operator_ref": str(OUT_O45.relative_to(REPO)), "operator": o45_obj.get("object"), "alpha45_mod_pi": clean_scalar(alpha45_mod_pi)},
            "pair2_to_pair4": {"operator_ref": str(OUT_O24.relative_to(REPO)), "operator": o24_obj.get("object"), "alpha24_mod_pi": clean_scalar(alpha24_mod_pi)},
            "pair3_to_pair5": {"operator_ref": str(OUT_O35.relative_to(REPO)), "operator": o35_obj.get("object"), "alpha35_mod_pi": clean_scalar(alpha35_mod_pi)},
        },
        "gluing_data": {
            "operator_section_ref": str(OUT_SECTION.relative_to(REPO)),
            "laws": section_obj["gluing_laws"],
            "local_cocycle_audits": section_obj["local_cocycle_audits"],
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

    summary = {
        "stage": "F465",
        "status": "F465_EXECUTED_CURRENT_STRICT_PAIR12345_SELECTOR_ATLAS_COCYCLE_DATA_EXPORT_PACKET_NO_FALSE_PASS",
        "generated_utc": generated_utc,
        "artifacts": {
            "a4": str(OUT_A4.relative_to(REPO)),
            "a5": str(OUT_A5.relative_to(REPO)),
            "o34": str(OUT_O34.relative_to(REPO)),
            "o45": str(OUT_O45.relative_to(REPO)),
            "o24": str(OUT_O24.relative_to(REPO)),
            "o35": str(OUT_O35.relative_to(REPO)),
            "section": str(OUT_SECTION.relative_to(REPO)),
            "atlas": str(OUT_ATLAS.relative_to(REPO)),
        },
        "audits": {
            "alpha12_mod_pi_consistency_distance": alpha12_mod_pi_consistency,
            "gluing34_max_abs_residual": gluing34_res,
            "gluing45_max_abs_residual": gluing45_res,
            "cocycle_1_to_3_max_abs_residual": cocycle_1_to_3_res,
            "cocycle_2_to_4_max_abs_residual": cocycle_2_to_4_res,
            "cocycle_3_to_5_max_abs_residual": cocycle_3_to_5_res,
        },
        "no_false_pass": True,
    }

    OUT_SECTION.write_text(json.dumps(section_obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_ATLAS.write_text(json.dumps(atlas, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()

