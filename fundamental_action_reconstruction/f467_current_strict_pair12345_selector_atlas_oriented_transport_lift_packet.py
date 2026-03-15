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

OUT_THETA_FAMILY = GENERATED / "theta_family_pair12345_oriented_mod_2pi_strict_convention_v1.json"

OUT_O13 = GENERATED / "o13_pair1_pair3_selector_chart_transport_operator_oriented_alpha13_mod_2pi_strict_convention_v1.json"
OUT_O14 = GENERATED / "o14_pair1_pair4_selector_chart_transport_operator_oriented_alpha14_mod_2pi_strict_convention_v1.json"
OUT_O15 = GENERATED / "o15_pair1_pair5_selector_chart_transport_operator_oriented_alpha15_mod_2pi_strict_convention_v1.json"
OUT_O23 = GENERATED / "o23_pair2_pair3_selector_chart_transport_operator_oriented_alpha23_mod_2pi_strict_convention_v1.json"
OUT_O24 = GENERATED / "o24_pair2_pair4_selector_chart_transport_operator_oriented_alpha24_mod_2pi_strict_convention_v1.json"
OUT_O25 = GENERATED / "o25_pair2_pair5_selector_chart_transport_operator_oriented_alpha25_mod_2pi_strict_convention_v1.json"
OUT_O34 = GENERATED / "o34_pair3_pair4_selector_chart_transport_operator_oriented_alpha34_mod_2pi_strict_convention_v1.json"
OUT_O35 = GENERATED / "o35_pair3_pair5_selector_chart_transport_operator_oriented_alpha35_mod_2pi_strict_convention_v1.json"
OUT_O45 = GENERATED / "o45_pair4_pair5_selector_chart_transport_operator_oriented_alpha45_mod_2pi_strict_convention_v1.json"

OUT_SECTION = GENERATED / "u_12345_pair12345_chart_glued_orientation_vector_section_oriented_mod_2pi_strict_convention_v1.json"
OUT_ATLAS = GENERATED / "selector_atlas_pair12345_oriented_transport_mod_2pi_strict_convention_v1.json"
OUT_SUMMARY = GENERATED / "f467_current_strict_pair12345_selector_atlas_oriented_transport_lift_packet_summary.json"


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


def circular_distance_mod_2pi(a: float, b: float) -> float:
    d = abs(angle_mod(a - b, 2.0 * math.pi))
    return float(min(d, 2.0 * math.pi - d))


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


def theta_mod_2pi_from_u(u: np.ndarray, c: np.ndarray, s: np.ndarray) -> float:
    coords = np.array([float(c @ u), float(s @ u)], dtype=float)
    norm = float(np.linalg.norm(coords))
    if norm == 0.0:
        raise ValueError("degenerate coordinates (zero)")
    theta = math.atan2(float(coords[1]), float(coords[0]))
    return angle_mod(theta, 2.0 * math.pi)


def export_oriented_transport_operator(
    *,
    object_name: str,
    alpha_mod_2pi: float,
    Ca: np.ndarray,
    Cb: np.ndarray,
    I: np.ndarray,
    out_path: Path,
    generated_utc: str,
    input_ref: str,
) -> tuple[dict[str, Any], np.ndarray, float, float]:
    O = build_swap_transport(I=I, Ca=Ca, Cb=Cb, alpha=alpha_mod_2pi)
    orth_res = max_abs(O.T @ O - I)
    invol_res = max_abs(O @ O - I)

    obj = {
        "object": object_name,
        "stage": "F467",
        "status": "actual_exported_lane_scoped_pair_chart_transport_operator__oriented_alpha_mod_2pi__sign_tracked_convention__no_false_pass",
        "as_of": AS_OF,
        "generated_utc": generated_utc,
        "inputs": {
            "alpha_mod_2pi_source_ref": input_ref,
            "alpha_mod_2pi": clean_scalar(alpha_mod_2pi),
        },
        "construction": {
            "transport_operator": (
                "O := C_to G(α) C_from^T + C_from G(-α) C_to^T + Π_rest with α := alpha_mod_2pi (oriented lift). "
                "This lift is sign-tracked: changing θ_m -> θ_m+π changes α by π and flips the transported oriented vector, "
                "while the transported projector is unchanged."
            ),
        },
        "outputs": {
            "n": int(I.shape[0]),
            "alpha_mod_2pi": clean_scalar(alpha_mod_2pi),
            "O": [[clean_scalar(float(x)) for x in row] for row in O.tolist()],
            "checks": {
                "orthogonality_max_abs_residual": orth_res,
                "involution_O_squared_equals_I_max_abs_residual": invol_res,
            },
        },
        "hard_limits": [
            "Oriented alpha mod 2π is a tracked gauge/convention lift, not a sign-sensitive physical orientation datum.",
            "Lane-scoped: does not export a global selector atlas nor global cocycle data on the full strict domain C_v1.",
            "Does not discharge QW-2191 nor imply strict-core selector closure.",
            "Does not claim ToE closure.",
        ],
        "no_false_pass": True,
    }
    out_path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    return obj, O, orth_res, invol_res


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    required = (IN_A1, IN_A2, IN_A3, IN_A4, IN_A5, IN_O12)
    missing = [str(p.relative_to(REPO)) for p in required if not p.exists()]
    if missing:
        raise SystemExit(
            json.dumps(
                {
                    "status": "NOT_COMPUTABLE_MISSING_INPUTS",
                    "missing": missing,
                    "expected": [
                        "F456 A_1(pair1) projector operator with u_1",
                        "F462 A_2(pair2) projector operator with u_2",
                        "F464 A_3(pair3) projector operator with u_3",
                        "F465 A_4(pair4) projector operator with u_4",
                        "F465 A_5(pair5) projector operator with u_5",
                        "F461 O12 transport operator with alpha12_mod_2pi",
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

    try:
        u1 = as_vector((a1_obj.get("data") or {}).get("u_1"), n=12, label="A1.u_1")
        u2 = as_vector((a2_obj.get("data") or {}).get("u_2"), n=12, label="A2.u_2")
        u3 = as_vector((a3_obj.get("data") or {}).get("u_3"), n=12, label="A3.u_3")
        u4 = as_vector((a4_obj.get("data") or {}).get("u_4"), n=12, label="A4.u_4")
        u5 = as_vector((a5_obj.get("data") or {}).get("u_5"), n=12, label="A5.u_5")
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

    # Oriented angles (mod 2π) induced by the exported oriented u vectors.
    theta1_2pi = theta_mod_2pi_from_u(u1, c1, s1)
    theta2_2pi = theta_mod_2pi_from_u(u2, c2, s2)
    theta3_2pi = theta_mod_2pi_from_u(u3, c3, s3)
    theta4_2pi = theta_mod_2pi_from_u(u4, c4, s4)
    theta5_2pi = theta_mod_2pi_from_u(u5, c5, s5)

    alpha12_mod_2pi_from_u = angle_mod(theta2_2pi - theta1_2pi, 2.0 * math.pi)
    alpha12_distance = circular_distance_mod_2pi(alpha12_mod_2pi_from_u, angle_mod(alpha12_mod_2pi, 2.0 * math.pi))

    # Remaining oriented edge angles (mod 2π).
    def a(theta_i: float, theta_j: float) -> float:
        return angle_mod(theta_j - theta_i, 2.0 * math.pi)

    alpha13 = a(theta1_2pi, theta3_2pi)
    alpha14 = a(theta1_2pi, theta4_2pi)
    alpha15 = a(theta1_2pi, theta5_2pi)
    alpha23 = a(theta2_2pi, theta3_2pi)
    alpha24 = a(theta2_2pi, theta4_2pi)
    alpha25 = a(theta2_2pi, theta5_2pi)
    alpha34 = a(theta3_2pi, theta4_2pi)
    alpha35 = a(theta3_2pi, theta5_2pi)
    alpha45 = a(theta4_2pi, theta5_2pi)

    generated_utc = datetime.now(timezone.utc).isoformat()

    # Export oriented transport operators (mod 2π lift) for all non-(1,2) edges.
    o13_obj, O13, o13_orth, o13_invol = export_oriented_transport_operator(
        object_name="O13_pair1_pair3_selector_chart_transport_operator_oriented_alpha13_mod_2pi_strict_convention_v1",
        alpha_mod_2pi=alpha13,
        Ca=C1,
        Cb=C3,
        I=I,
        out_path=OUT_O13,
        generated_utc=generated_utc,
        input_ref=str(OUT_THETA_FAMILY.relative_to(REPO)),
    )
    o14_obj, O14, o14_orth, o14_invol = export_oriented_transport_operator(
        object_name="O14_pair1_pair4_selector_chart_transport_operator_oriented_alpha14_mod_2pi_strict_convention_v1",
        alpha_mod_2pi=alpha14,
        Ca=C1,
        Cb=C4,
        I=I,
        out_path=OUT_O14,
        generated_utc=generated_utc,
        input_ref=str(OUT_THETA_FAMILY.relative_to(REPO)),
    )
    o15_obj, O15, o15_orth, o15_invol = export_oriented_transport_operator(
        object_name="O15_pair1_pair5_selector_chart_transport_operator_oriented_alpha15_mod_2pi_strict_convention_v1",
        alpha_mod_2pi=alpha15,
        Ca=C1,
        Cb=C5,
        I=I,
        out_path=OUT_O15,
        generated_utc=generated_utc,
        input_ref=str(OUT_THETA_FAMILY.relative_to(REPO)),
    )
    o23_obj, O23, o23_orth, o23_invol = export_oriented_transport_operator(
        object_name="O23_pair2_pair3_selector_chart_transport_operator_oriented_alpha23_mod_2pi_strict_convention_v1",
        alpha_mod_2pi=alpha23,
        Ca=C2,
        Cb=C3,
        I=I,
        out_path=OUT_O23,
        generated_utc=generated_utc,
        input_ref=str(OUT_THETA_FAMILY.relative_to(REPO)),
    )
    o24_obj, O24, o24_orth, o24_invol = export_oriented_transport_operator(
        object_name="O24_pair2_pair4_selector_chart_transport_operator_oriented_alpha24_mod_2pi_strict_convention_v1",
        alpha_mod_2pi=alpha24,
        Ca=C2,
        Cb=C4,
        I=I,
        out_path=OUT_O24,
        generated_utc=generated_utc,
        input_ref=str(OUT_THETA_FAMILY.relative_to(REPO)),
    )
    o25_obj, O25, o25_orth, o25_invol = export_oriented_transport_operator(
        object_name="O25_pair2_pair5_selector_chart_transport_operator_oriented_alpha25_mod_2pi_strict_convention_v1",
        alpha_mod_2pi=alpha25,
        Ca=C2,
        Cb=C5,
        I=I,
        out_path=OUT_O25,
        generated_utc=generated_utc,
        input_ref=str(OUT_THETA_FAMILY.relative_to(REPO)),
    )
    o34_obj, O34, o34_orth, o34_invol = export_oriented_transport_operator(
        object_name="O34_pair3_pair4_selector_chart_transport_operator_oriented_alpha34_mod_2pi_strict_convention_v1",
        alpha_mod_2pi=alpha34,
        Ca=C3,
        Cb=C4,
        I=I,
        out_path=OUT_O34,
        generated_utc=generated_utc,
        input_ref=str(OUT_THETA_FAMILY.relative_to(REPO)),
    )
    o35_obj, O35, o35_orth, o35_invol = export_oriented_transport_operator(
        object_name="O35_pair3_pair5_selector_chart_transport_operator_oriented_alpha35_mod_2pi_strict_convention_v1",
        alpha_mod_2pi=alpha35,
        Ca=C3,
        Cb=C5,
        I=I,
        out_path=OUT_O35,
        generated_utc=generated_utc,
        input_ref=str(OUT_THETA_FAMILY.relative_to(REPO)),
    )
    o45_obj, O45, o45_orth, o45_invol = export_oriented_transport_operator(
        object_name="O45_pair4_pair5_selector_chart_transport_operator_oriented_alpha45_mod_2pi_strict_convention_v1",
        alpha_mod_2pi=alpha45,
        Ca=C4,
        Cb=C5,
        I=I,
        out_path=OUT_O45,
        generated_utc=generated_utc,
        input_ref=str(OUT_THETA_FAMILY.relative_to(REPO)),
    )

    # Vector transport audits (oriented).
    def v_res(O: np.ndarray, a: np.ndarray, b: np.ndarray) -> float:
        return float(np.linalg.norm((O @ a) - b))

    transport = {
        "12": {"u2_minus_O12_u1_l2": v_res(O12, u1, u2), "u1_minus_O12_u2_l2": v_res(O12, u2, u1)},
        "13": {"u3_minus_O13_u1_l2": v_res(O13, u1, u3), "u1_minus_O13_u3_l2": v_res(O13, u3, u1)},
        "14": {"u4_minus_O14_u1_l2": v_res(O14, u1, u4), "u1_minus_O14_u4_l2": v_res(O14, u4, u1)},
        "15": {"u5_minus_O15_u1_l2": v_res(O15, u1, u5), "u1_minus_O15_u5_l2": v_res(O15, u5, u1)},
        "23": {"u3_minus_O23_u2_l2": v_res(O23, u2, u3), "u2_minus_O23_u3_l2": v_res(O23, u3, u2)},
        "24": {"u4_minus_O24_u2_l2": v_res(O24, u2, u4), "u2_minus_O24_u4_l2": v_res(O24, u4, u2)},
        "25": {"u5_minus_O25_u2_l2": v_res(O25, u2, u5), "u2_minus_O25_u5_l2": v_res(O25, u5, u2)},
        "34": {"u4_minus_O34_u3_l2": v_res(O34, u3, u4), "u3_minus_O34_u4_l2": v_res(O34, u4, u3)},
        "35": {"u5_minus_O35_u3_l2": v_res(O35, u3, u5), "u3_minus_O35_u5_l2": v_res(O35, u5, u3)},
        "45": {"u5_minus_O45_u4_l2": v_res(O45, u4, u5), "u4_minus_O45_u5_l2": v_res(O45, u5, u4)},
    }

    # Full triple cocycle/path-independence audits on oriented vectors.
    def cocycle_residual(*, path: np.ndarray, direct: np.ndarray, u: np.ndarray) -> float:
        return float(np.linalg.norm((path @ u) - (direct @ u)))

    cocycle = [
        {
            "triple": "pair1-pair2-pair3",
            "path": "O_23 O_12",
            "direct": "O_13",
            "u_base": "u_1",
            "l2_residual": cocycle_residual(path=(O23 @ O12), direct=O13, u=u1),
        },
        {
            "triple": "pair1-pair2-pair4",
            "path": "O_24 O_12",
            "direct": "O_14",
            "u_base": "u_1",
            "l2_residual": cocycle_residual(path=(O24 @ O12), direct=O14, u=u1),
        },
        {
            "triple": "pair1-pair2-pair5",
            "path": "O_25 O_12",
            "direct": "O_15",
            "u_base": "u_1",
            "l2_residual": cocycle_residual(path=(O25 @ O12), direct=O15, u=u1),
        },
        {
            "triple": "pair1-pair3-pair4",
            "path": "O_34 O_13",
            "direct": "O_14",
            "u_base": "u_1",
            "l2_residual": cocycle_residual(path=(O34 @ O13), direct=O14, u=u1),
        },
        {
            "triple": "pair1-pair3-pair5",
            "path": "O_35 O_13",
            "direct": "O_15",
            "u_base": "u_1",
            "l2_residual": cocycle_residual(path=(O35 @ O13), direct=O15, u=u1),
        },
        {
            "triple": "pair1-pair4-pair5",
            "path": "O_45 O_14",
            "direct": "O_15",
            "u_base": "u_1",
            "l2_residual": cocycle_residual(path=(O45 @ O14), direct=O15, u=u1),
        },
        {
            "triple": "pair2-pair3-pair4",
            "path": "O_34 O_23",
            "direct": "O_24",
            "u_base": "u_2",
            "l2_residual": cocycle_residual(path=(O34 @ O23), direct=O24, u=u2),
        },
        {
            "triple": "pair2-pair3-pair5",
            "path": "O_35 O_23",
            "direct": "O_25",
            "u_base": "u_2",
            "l2_residual": cocycle_residual(path=(O35 @ O23), direct=O25, u=u2),
        },
        {
            "triple": "pair2-pair4-pair5",
            "path": "O_45 O_24",
            "direct": "O_25",
            "u_base": "u_2",
            "l2_residual": cocycle_residual(path=(O45 @ O24), direct=O25, u=u2),
        },
        {
            "triple": "pair3-pair4-pair5",
            "path": "O_45 O_34",
            "direct": "O_35",
            "u_base": "u_3",
            "l2_residual": cocycle_residual(path=(O45 @ O34), direct=O35, u=u3),
        },
    ]

    theta_family = {
        "object": "ThetaFamily_pair12345_oriented_mod_2pi_strict_convention_v1",
        "stage": "F467",
        "status": "actual_exported_lane_scoped_theta_family__oriented_lift_induced_by_exported_u_vectors__sign_tracked_convention__no_false_pass",
        "as_of": AS_OF,
        "generated_utc": generated_utc,
        "intent": (
            "Export one lane-scoped oriented theta-family lift (angles mod 2π) on {pair1..pair5} induced by the currently exported "
            "oriented representative vectors u_1..u_5. This is a tracked gauge/convention layer (sign-tracked); it does not claim a "
            "sign-sensitive physical orientation datum and does not discharge QW-2191."
        ),
        "inputs": {
            "a1_pair1": str(IN_A1.relative_to(REPO)),
            "a2_pair2": str(IN_A2.relative_to(REPO)),
            "a3_pair3": str(IN_A3.relative_to(REPO)),
            "a4_pair4": str(IN_A4.relative_to(REPO)),
            "a5_pair5": str(IN_A5.relative_to(REPO)),
            "o12_pair1_pair2": str(IN_O12.relative_to(REPO)),
        },
        "outputs": {
            "n": n,
            "pair1": {"m": 1, "theta_1_mod_2pi": clean_scalar(theta1_2pi)},
            "pair2": {"m": 2, "theta_2_mod_2pi": clean_scalar(theta2_2pi)},
            "pair3": {"m": 3, "theta_3_mod_2pi": clean_scalar(theta3_2pi)},
            "pair4": {"m": 4, "theta_4_mod_2pi": clean_scalar(theta4_2pi)},
            "pair5": {"m": 5, "theta_5_mod_2pi": clean_scalar(theta5_2pi)},
            "derived_alphas_mod_2pi": {
                "alpha12_mod_2pi": clean_scalar(alpha12_mod_2pi_from_u),
                "alpha13_mod_2pi": clean_scalar(alpha13),
                "alpha14_mod_2pi": clean_scalar(alpha14),
                "alpha15_mod_2pi": clean_scalar(alpha15),
                "alpha23_mod_2pi": clean_scalar(alpha23),
                "alpha24_mod_2pi": clean_scalar(alpha24),
                "alpha25_mod_2pi": clean_scalar(alpha25),
                "alpha34_mod_2pi": clean_scalar(alpha34),
                "alpha35_mod_2pi": clean_scalar(alpha35),
                "alpha45_mod_2pi": clean_scalar(alpha45),
            },
        },
        "audits": {
            "alpha12_mod_2pi_distance_u_vs_exported_o12": alpha12_distance,
        },
        "hard_limits": [
            "Oriented angles mod 2π are a tracked gauge/convention lift, not a sign-sensitive physical orientation datum.",
            "Lane-scoped: does not export a global selector atlas on C_v1.",
            "Does not discharge QW-2191 nor imply strict-core selector closure.",
            "Does not claim ToE closure.",
        ],
        "no_false_pass": True,
    }

    # Write theta family first so transport operators can reference it.
    OUT_THETA_FAMILY.write_text(json.dumps(theta_family, indent=2, ensure_ascii=True) + "\n", encoding="ascii")

    section_obj = {
        "object": "U_12345_pair12345_chart_glued_orientation_vector_section_oriented_mod_2pi_strict_convention_v1",
        "stage": "F467",
        "status": "actual_exported_lane_scoped_oriented_vector_section_with_full_triple_cocycle_audits__sign_tracked_convention__no_false_pass",
        "as_of": AS_OF,
        "generated_utc": generated_utc,
        "intent": (
            "Export one explicit oriented (mod 2π) vector-level lift of the already exported projector-level five-chart atlas ingredient on {pair1..pair5}. "
            "This packages the currently exported oriented representative vectors u_1..u_5 together with oriented chart-transport operators (α mod 2π) "
            "so that transport and full triple cocycle/path-independence holds at the level of the exported vectors. "
            "This is a tracked gauge/convention layer and does not claim a sign-sensitive physical orientation datum."
        ),
        "inputs": {
            "theta_family_ref": str(OUT_THETA_FAMILY.relative_to(REPO)),
            "o12_ref": str(IN_O12.relative_to(REPO)),
            "o13_ref": str(OUT_O13.relative_to(REPO)),
            "o14_ref": str(OUT_O14.relative_to(REPO)),
            "o15_ref": str(OUT_O15.relative_to(REPO)),
            "o23_ref": str(OUT_O23.relative_to(REPO)),
            "o24_ref": str(OUT_O24.relative_to(REPO)),
            "o25_ref": str(OUT_O25.relative_to(REPO)),
            "o34_ref": str(OUT_O34.relative_to(REPO)),
            "o35_ref": str(OUT_O35.relative_to(REPO)),
            "o45_ref": str(OUT_O45.relative_to(REPO)),
        },
        "outputs": {
            "u_vectors": {
                "u_1": [clean_scalar(float(x)) for x in u1.tolist()],
                "u_2": [clean_scalar(float(x)) for x in u2.tolist()],
                "u_3": [clean_scalar(float(x)) for x in u3.tolist()],
                "u_4": [clean_scalar(float(x)) for x in u4.tolist()],
                "u_5": [clean_scalar(float(x)) for x in u5.tolist()],
            },
        },
        "audits": {
            "transport_vector_level": transport,
            "cocycle_vector_level": cocycle,
            "orthogonality_max_abs_residual": {
                "o12": max_abs(O12.T @ O12 - I),
                "o13": o13_orth,
                "o14": o14_orth,
                "o15": o15_orth,
                "o23": o23_orth,
                "o24": o24_orth,
                "o25": o25_orth,
                "o34": o34_orth,
                "o35": o35_orth,
                "o45": o45_orth,
            },
            "involution_max_abs_residual": {
                "o12": max_abs(O12 @ O12 - I),
                "o13": o13_invol,
                "o14": o14_invol,
                "o15": o15_invol,
                "o23": o23_invol,
                "o24": o24_invol,
                "o25": o25_invol,
                "o34": o34_invol,
                "o35": o35_invol,
                "o45": o45_invol,
            },
        },
        "hard_limits": [
            "Vector-level lift is a tracked gauge/convention layer; it does not derive a sign-sensitive physical orientation datum.",
            "Lane-scoped: does not export a global selector atlas nor global overlap-domain declaration on C_v1.",
            "Does not discharge QW-2191 nor imply strict-core selector closure / admissible S_sel_int.",
            "Does not claim ToE closure.",
        ],
        "no_false_pass": True,
    }

    atlas = {
        "object": "SelectorAtlas_pair12345_oriented_transport_mod_2pi_strict_convention_v1",
        "stage": "F467",
        "status": "actual_exported_lane_scoped_five_chart_selector_atlas_with_oriented_transport_mod_2pi__full_triple_cocycle_vector_section__sign_tracked_convention__no_false_pass",
        "as_of": AS_OF,
        "generated_utc": generated_utc,
        "intent": (
            "Export one explicit lane-scoped five-chart selector-atlas ingredient on {pair1..pair5} whose transport operators are lifted to "
            "oriented angles α mod 2π (sign-tracked convention) so that the exported oriented representative vectors u_1..u_5 satisfy full triple "
            "cocycle/path-independence at vector level. This does not export a global selector atlas on C_v1 and does not discharge QW-2191."
        ),
        "overlap_domain_declaration": {
            "scope": "exported_artifact_overlap_on_n12_fourier_carrier",
            "meaning": (
                "This is an exported-artifact overlap declaration: all five chart operators and oriented transport operators are simultaneously defined "
                "as strict exported artifacts on the declared n=12 Fourier carrier. This is not a global open-cover atlas on the full strict domain C_v1."
            ),
        },
        "inputs": {
            "projector_level_atlas_ref": "fundamental_action_reconstruction/generated/selector_atlas_pair12345_axis_only_projector_v2.json",
            "theta_family_ref": str(OUT_THETA_FAMILY.relative_to(REPO)),
            "vector_section_ref": str(OUT_SECTION.relative_to(REPO)),
        },
        "transitions": {
            "pair1_to_pair2": {
                "operator_ref": str(IN_O12.relative_to(REPO)),
                "operator": o12_obj.get("object"),
                "alpha12_mod_2pi": clean_scalar(alpha12_mod_2pi),
            },
            "pair1_to_pair3": {"operator_ref": str(OUT_O13.relative_to(REPO)), "operator": o13_obj.get("object"), "alpha13_mod_2pi": clean_scalar(alpha13)},
            "pair1_to_pair4": {"operator_ref": str(OUT_O14.relative_to(REPO)), "operator": o14_obj.get("object"), "alpha14_mod_2pi": clean_scalar(alpha14)},
            "pair1_to_pair5": {"operator_ref": str(OUT_O15.relative_to(REPO)), "operator": o15_obj.get("object"), "alpha15_mod_2pi": clean_scalar(alpha15)},
            "pair2_to_pair3": {"operator_ref": str(OUT_O23.relative_to(REPO)), "operator": o23_obj.get("object"), "alpha23_mod_2pi": clean_scalar(alpha23)},
            "pair2_to_pair4": {"operator_ref": str(OUT_O24.relative_to(REPO)), "operator": o24_obj.get("object"), "alpha24_mod_2pi": clean_scalar(alpha24)},
            "pair2_to_pair5": {"operator_ref": str(OUT_O25.relative_to(REPO)), "operator": o25_obj.get("object"), "alpha25_mod_2pi": clean_scalar(alpha25)},
            "pair3_to_pair4": {"operator_ref": str(OUT_O34.relative_to(REPO)), "operator": o34_obj.get("object"), "alpha34_mod_2pi": clean_scalar(alpha34)},
            "pair3_to_pair5": {"operator_ref": str(OUT_O35.relative_to(REPO)), "operator": o35_obj.get("object"), "alpha35_mod_2pi": clean_scalar(alpha35)},
            "pair4_to_pair5": {"operator_ref": str(OUT_O45.relative_to(REPO)), "operator": o45_obj.get("object"), "alpha45_mod_2pi": clean_scalar(alpha45)},
        },
        "audits": {
            "alpha12_mod_2pi_consistency_distance_u_vs_exported_o12": alpha12_distance,
            "transport_vector_level": transport,
            "cocycle_vector_level": cocycle,
        },
        "hard_limits": [
            "Oriented transport is a tracked gauge/convention lift; it does not derive a sign-sensitive physical orientation datum.",
            "Lane-scoped: does not export a global selector atlas nor global overlap-domain declaration on C_v1 (H41 remains open globally).",
            "Does not discharge QW-2191 nor export strict-core selector closure / admissible S_sel_int (H40 remains open globally).",
            "Does not claim ToE closure.",
        ],
        "no_false_pass": True,
    }

    # Summary.
    cocycle_max = max((float(item["l2_residual"]) for item in cocycle), default=0.0)
    transport_max = max(float(v) for edge in transport.values() for v in edge.values())

    summary = {
        "stage": "F467",
        "status": "F467_EXECUTED_CURRENT_STRICT_PAIR12345_SELECTOR_ATLAS_ORIENTED_TRANSPORT_LIFT_PACKET_NO_FALSE_PASS",
        "generated_utc": generated_utc,
        "artifacts": {
            "theta_family": str(OUT_THETA_FAMILY.relative_to(REPO)),
            "o13": str(OUT_O13.relative_to(REPO)),
            "o14": str(OUT_O14.relative_to(REPO)),
            "o15": str(OUT_O15.relative_to(REPO)),
            "o23": str(OUT_O23.relative_to(REPO)),
            "o24": str(OUT_O24.relative_to(REPO)),
            "o25": str(OUT_O25.relative_to(REPO)),
            "o34": str(OUT_O34.relative_to(REPO)),
            "o35": str(OUT_O35.relative_to(REPO)),
            "o45": str(OUT_O45.relative_to(REPO)),
            "vector_section": str(OUT_SECTION.relative_to(REPO)),
            "atlas": str(OUT_ATLAS.relative_to(REPO)),
        },
        "audits": {
            "alpha12_mod_2pi_consistency_distance_u_vs_exported_o12": alpha12_distance,
            "transport_vector_level_max_l2_residual": transport_max,
            "cocycle_vector_level_max_l2_residual": cocycle_max,
        },
        "no_false_pass": True,
    }

    OUT_SECTION.write_text(json.dumps(section_obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_ATLAS.write_text(json.dumps(atlas, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()

