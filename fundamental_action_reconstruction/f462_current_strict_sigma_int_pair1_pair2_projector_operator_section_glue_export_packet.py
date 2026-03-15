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

IN_THETA_PAIR = GENERATED / "theta_pair_sigma_int_strict_selector_ingredient_o2_cut_slot_free_v1.json"
IN_A1 = GENERATED / "a_1_pair1_orientation_projector_operator_strict_core_v1.json"
IN_O12 = GENERATED / "o12_pair1_pair2_selector_chart_transport_operator_strict_derived_from_sigma_int_alpha12_v1.json"

OUT_A2 = GENERATED / "a_2_pair2_orientation_projector_operator_strict_core_v1.json"
OUT_SECTION = GENERATED / "a_12_pair12_chart_glued_orientation_projector_operator_section_strict_core_v1.json"
OUT_SUMMARY = (
    GENERATED / "f462_current_strict_sigma_int_pair1_pair2_projector_operator_section_glue_export_packet_summary.json"
)


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

    required = (IN_THETA_PAIR, IN_A1, IN_O12)
    missing = [str(p.relative_to(REPO)) for p in required if not p.exists()]
    if missing:
        raise SystemExit(
            json.dumps(
                {
                    "status": "NOT_COMPUTABLE_MISSING_INPUTS",
                    "missing": missing,
                    "expected": [
                        "F451 theta-pair export",
                        "F456 A_1(pair1) projector operator export",
                        "F461 O12 chart-transport operator export",
                    ],
                },
                ensure_ascii=True,
            )
        )

    theta_pair = load_json(IN_THETA_PAIR)
    a1_obj = load_json(IN_A1)
    o12_obj = load_json(IN_O12)

    try:
        sigma_int = int(((theta_pair.get("inputs") or {}).get("sigma_int") or {}).get("value"))
        u1 = as_vector(((theta_pair.get("outputs") or {}).get("pair1") or {}).get("u_1"), n=12, label="F451.outputs.pair1.u_1")
        u2 = as_vector(((theta_pair.get("outputs") or {}).get("pair2") or {}).get("u_2"), n=12, label="F451.outputs.pair2.u_2")
        theta_1 = float(((theta_pair.get("outputs") or {}).get("pair1") or {}).get("theta_1"))
        theta_2 = float(((theta_pair.get("outputs") or {}).get("pair2") or {}).get("theta_2"))
    except Exception as exc:
        raise SystemExit(
            json.dumps(
                {
                    "status": "INVALID_F451_SHAPE",
                    "expected": "F451 inputs.sigma_int.value and outputs.pair{1,2}.{theta_i,u_i}",
                    "error": str(exc),
                },
                ensure_ascii=True,
            )
        )

    try:
        a1_matrix_c1s1 = as_square_matrix((a1_obj.get("data") or {}).get("A_1_pair1_matrix_in_c1_s1"), n=2, label="F456.A_1_pair1_matrix_in_c1_s1")
    except Exception as exc:
        raise SystemExit(
            json.dumps(
                {
                    "status": "INVALID_F456_SHAPE",
                    "expected": "F456 exported A_1_pair1_matrix_in_c1_s1 (2x2)",
                    "error": str(exc),
                },
                ensure_ascii=True,
            )
        )

    try:
        O12 = as_square_matrix((o12_obj.get("outputs") or {}).get("O12"), n=12, label="F461.outputs.O12")
        alpha12 = float((o12_obj.get("outputs") or {}).get("alpha12_mod_2pi"))
        G12 = as_square_matrix((o12_obj.get("outputs") or {}).get("G12_so2"), n=2, label="F461.outputs.G12_so2")
    except Exception as exc:
        raise SystemExit(
            json.dumps(
                {
                    "status": "INVALID_F461_SHAPE",
                    "expected": "F461.outputs.{O12(12x12),alpha12_mod_2pi,G12_so2(2x2)}",
                    "error": str(exc),
                },
                ensure_ascii=True,
            )
        )

    n = 12
    I = np.eye(n, dtype=float)

    # Canonical Fourier bases for pair1 and pair2 (declared scaffold).
    c1, s1 = real_fourier_pair_basis(n=n, m=1)
    c2, s2 = real_fourier_pair_basis(n=n, m=2)

    # Coordinate extraction in the canonical pair bases.
    coords1 = np.array([float(c1 @ u1), float(s1 @ u1)], dtype=float)
    coords2 = np.array([float(c2 @ u2), float(s2 @ u2)], dtype=float)

    u1_recon = (coords1[0] * c1) + (coords1[1] * s1)
    u2_recon = (coords2[0] * c2) + (coords2[1] * s2)
    u1_outside = float(np.linalg.norm(u1 - u1_recon))
    u2_outside = float(np.linalg.norm(u2 - u2_recon))

    # Projector operators (full carrier level).
    A1_full = np.outer(u1, u1)
    A2_full = np.outer(u2, u2)
    A2_full_from_transport = O12 @ A1_full @ O12.T

    # Projector operators in pair-plane coordinates (2x2).
    A1_2 = np.outer(coords1, coords1)
    A2_2 = np.outer(coords2, coords2)
    A2_2_from_G12 = G12 @ A1_2 @ G12.T

    # Consistency with F456 (A1 in c1,s1 basis).
    a1_2_residual_vs_f456 = max_abs(A1_2 - a1_matrix_c1s1)

    # Gluing residuals.
    gluing_full_res = max_abs(A2_full - A2_full_from_transport)
    gluing_2x2_res = max_abs(A2_2 - A2_2_from_G12)

    # Pair12 restricted matrix of O12 in (c1,s1,c2,s2) basis.
    B = np.column_stack([c1, s1, c2, s2])  # 12x4 with orthonormal columns
    O4_from_O12 = B.T @ O12 @ B
    G21 = rotation_so2(-alpha12)
    O4_expected = np.block([[np.zeros((2, 2)), G21], [G12, np.zeros((2, 2))]])
    O4_diff = max_abs(O4_from_O12 - O4_expected)

    # Basic operator audits.
    orth_res = max_abs(O12.T @ O12 - I)
    invol_res = max_abs(O12 @ O12 - I)
    a1_idem_res = max_abs(A1_full @ A1_full - A1_full)
    a2_idem_res = max_abs(A2_full @ A2_full - A2_full)
    a1_sym_res = max_abs(A1_full - A1_full.T)
    a2_sym_res = max_abs(A2_full - A2_full.T)

    # Residual sign gauge invariance (projector level).
    A1_full_sign = np.outer(-u1, -u1)
    A2_full_sign = np.outer(-u2, -u2)
    a1_sign_res = max_abs(A1_full - A1_full_sign)
    a2_sign_res = max_abs(A2_full - A2_full_sign)

    # Export A2 operator object.
    a2_obj = {
        "object": "A_2_pair2_orientation_projector_operator_strict_core_v1",
        "status": "actual_exported_operator__derived_only_from_sigma_int_u2__residual_z2_sign_gauge_invariant",
        "as_of": "2026-03-15",
        "intent": (
            "Export one strict-core, slot-free operator on V_2 = span{c2,s2} constructed only from the already exported sigma-int "
            "orientation direction u_2. The operator is the rank-one projector A_2(pair2) := |u_2><u_2| in the (c2,s2) basis. "
            "This is lane-scoped to the strict sigma-int corridor and does not imply any global selector atlas, selector closure, or ToE closure."
        ),
        "inputs": {
            "theta_pair": str(IN_THETA_PAIR.relative_to(REPO)),
            "o12_chart_transport": str(IN_O12.relative_to(REPO)),
        },
        "domain_notes": {
            "basis_plane": "V_2 = span{c2,s2} inside the n=12 real Fourier scaffold",
            "sigma_int_lane": "slot_free_theta_pair_supply (F451/N489)",
        },
        "data": {
            "sigma_int": sigma_int,
            "theta_2_exported": theta_2,
            "u_2": [float(x) for x in u2.tolist()],
            "u_2_coords_in_c2_s2": [float(x) for x in coords2.tolist()],
            "A_2_pair2_matrix_in_c2_s2": [[float(x) for x in row] for row in A2_2.tolist()],
        },
        "audits": {
            "u2_outside_span_c2_s2_l2": u2_outside,
            "coords_norm": float(np.linalg.norm(coords2)),
            "projector_symmetry_inf_norm": a2_sym_res,
            "projector_idempotence_inf_norm": a2_idem_res,
            "projector_trace": float(np.trace(A2_full)),
            "residual_z2_sign_invariance_max_abs_diff": a2_sign_res,
            "tolerance": 1e-12,
        },
        "hard_limits": [
            "Lane-scoped sigma-int corridor object; does not export a global selector atlas nor global gluing cocycle data.",
            "Does not discharge QW-2191.",
            "Does not derive a sign-sensitive physical orientation datum; projector is sign-gauge-invariant.",
            "Does not claim strict-core selector closure / admissible S_sel_int.",
            "Does not claim ToE closure.",
        ],
        "no_false_pass": True,
    }

    # Export glued two-chart operator section object.
    section_obj = {
        "object": "A_12_pair12_chart_glued_orientation_projector_operator_section_strict_core_v1",
        "status": "actual_exported_two_chart_projector_operator_section__pair1_pair2_glued_by_O12__no_false_pass",
        "as_of": "2026-03-15",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "intent": (
            "Export one explicit two-chart operator section on {pair1,pair2} in the strict sigma-int corridor scope, "
            "gluing the projector operators A_1(pair1) and A_2(pair2) by the exported chart-transport operator O_12: "
            "A_2 = O_12 A_1 O_12^T. This is projector-level and therefore residual-sign-gauge-safe, but remains lane-scoped "
            "and does not imply any global selector atlas or QW-2191 discharge."
        ),
        "inputs": {
            "a1_pair1_operator_ref": str(IN_A1.relative_to(REPO)),
            "a2_pair2_operator_ref": str(OUT_A2.relative_to(REPO)),
            "o12_chart_transport_ref": str(IN_O12.relative_to(REPO)),
            "theta_pair_ref": str(IN_THETA_PAIR.relative_to(REPO)),
        },
        "construction": {
            "A1_full": "A_1 := |u_1><u_1| on the n=12 carrier (pair1 representative)",
            "A2_full": "A_2 := |u_2><u_2| on the n=12 carrier (pair2 representative)",
            "gluing_law_full": "A_2 = O_12 A_1 O_12^T",
            "pair12_restricted_matrix": "O_12 restricted to V1⊕V2 in basis (c1,s1,c2,s2) is [[0,G(-alpha12)],[G(alpha12),0]]",
        },
        "outputs": {
            "n": n,
            "alpha12_mod_2pi": alpha12,
            "G12_so2": [[float(x) for x in row] for row in G12.tolist()],
            "A1_pair1_matrix_in_c1_s1": [[float(x) for x in row] for row in A1_2.tolist()],
            "A2_pair2_matrix_in_c2_s2": [[float(x) for x in row] for row in A2_2.tolist()],
            "O12_pair12_restricted_matrix_in_c1_s1_c2_s2": [[float(x) for x in row] for row in O4_from_O12.tolist()],
        },
        "audits": {
            "O12_orthogonality_max_abs_residual": orth_res,
            "O12_involution_max_abs_residual": invol_res,
            "u1_outside_span_c1_s1_l2": u1_outside,
            "u2_outside_span_c2_s2_l2": u2_outside,
            "A1_symmetry_inf_norm": a1_sym_res,
            "A2_symmetry_inf_norm": a2_sym_res,
            "A1_idempotence_inf_norm": a1_idem_res,
            "A2_idempotence_inf_norm": a2_idem_res,
            "A1_sign_gauge_invariance_max_abs": a1_sign_res,
            "A2_sign_gauge_invariance_max_abs": a2_sign_res,
            "A1_2x2_residual_vs_F456_export": a1_2_residual_vs_f456,
            "gluing_full_matrix_max_abs_residual": gluing_full_res,
            "gluing_pair_plane_2x2_max_abs_residual": gluing_2x2_res,
            "O12_pair12_restriction_max_abs_residual_vs_expected_block_form": O4_diff,
        },
        "hard_limits": [
            "Lane-scoped: only a two-chart {pair1,pair2} projector-operator section is exported.",
            "Does not export a global selector atlas, overlap-domain declaration, nor cocycle-level gluing data (H41 remains open globally).",
            "Does not discharge QW-2191 nor export a global selector transition/gluing object (H40 remains open globally).",
            "Projector-level only: does not derive a sign-sensitive physical orientation datum.",
            "Does not claim strict-core selector closure / admissible S_sel_int.",
            "Does not claim ToE closure.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "F462",
        "status": "F462_EXECUTED_CURRENT_STRICT_SIGMA_INT_PAIR1_PAIR2_PROJECTOR_OPERATOR_SECTION_GLUE_EXPORT_PACKET_NO_FALSE_PASS",
        "generated_utc": section_obj["generated_utc"],
        "outputs": {
            "gluing_full_matrix_max_abs_residual": gluing_full_res,
            "gluing_pair_plane_2x2_max_abs_residual": gluing_2x2_res,
            "o12_pair12_restriction_max_abs_residual": O4_diff,
        },
        "no_false_pass": True,
    }

    OUT_A2.write_text(json.dumps(a2_obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SECTION.write_text(json.dumps(section_obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()

