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

IN_ALPHA12 = GENERATED / "alpha12_pair1_pair2_transition_angle_strict_derived_from_sigma_int_slot_free_theta_pair_v1.json"

OUT_OBJECT = GENERATED / "o12_pair1_pair2_selector_chart_transport_operator_strict_derived_from_sigma_int_alpha12_v1.json"
OUT_SUMMARY = (
    GENERATED / "o12_pair1_pair2_selector_chart_transport_operator_strict_derived_from_sigma_int_alpha12_v1_summary.json"
)


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def max_abs(a: np.ndarray) -> float:
    return float(np.max(np.abs(a)))


def as_vector(x: Any, n: int, label: str) -> np.ndarray:
    if not (isinstance(x, list) and len(x) == n and all(isinstance(v, (int, float)) for v in x)):
        raise ValueError(f"{label} must be a length-{n} list of finite numbers")
    v = np.array([float(v) for v in x], dtype=float)
    if not np.isfinite(v).all():
        raise ValueError(f"{label} must contain only finite numbers")
    return v


def real_fourier_pair_basis(n: int, m: int) -> tuple[np.ndarray, np.ndarray]:
    """Return the canonical real Fourier pair (c_m, s_m) on Z_n for 1 <= m <= n/2-1."""
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

    if not IN_ALPHA12.exists():
        raise SystemExit(
            json.dumps(
                {
                    "status": "NOT_COMPUTABLE_MISSING_INPUT",
                    "missing": str(IN_ALPHA12.relative_to(REPO)),
                    "expected": "F457 alpha12 export object",
                },
                ensure_ascii=True,
            )
        )

    alpha12_obj = load_json(IN_ALPHA12)
    try:
        alpha_mod_2pi = float((alpha12_obj.get("outputs") or {}).get("alpha_12_mod_2pi"))
        theta_pair_ref = str(((alpha12_obj.get("inputs") or {}).get("theta_pair_source") or {}).get("ref"))
        if not theta_pair_ref:
            raise ValueError("missing inputs.theta_pair_source.ref")
    except Exception as exc:
        raise SystemExit(
            json.dumps(
                {
                    "status": "INVALID_INPUT_SHAPE",
                    "expected": "F457 object with inputs.theta_pair_source.ref and outputs.alpha_12_mod_2pi",
                    "error": str(exc),
                },
                ensure_ascii=True,
            )
        )

    theta_pair_path = Path(theta_pair_ref)
    if not theta_pair_path.is_absolute():
        theta_pair_path = REPO / theta_pair_path
    if not theta_pair_path.exists():
        raise SystemExit(
            json.dumps(
                {
                    "status": "NOT_COMPUTABLE_MISSING_INPUT",
                    "missing": str(theta_pair_path.relative_to(REPO)),
                    "expected": "F451 theta-pair export object referenced by F457",
                },
                ensure_ascii=True,
            )
        )

    theta_pair_obj = load_json(theta_pair_path)
    try:
        u1 = as_vector(((theta_pair_obj.get("outputs") or {}).get("pair1") or {}).get("u_1"), n=12, label="F451.outputs.pair1.u_1")
        u2 = as_vector(((theta_pair_obj.get("outputs") or {}).get("pair2") or {}).get("u_2"), n=12, label="F451.outputs.pair2.u_2")
        sigma_int = int(((theta_pair_obj.get("inputs") or {}).get("sigma_int") or {}).get("value"))
        theta_1 = float(((theta_pair_obj.get("outputs") or {}).get("pair1") or {}).get("theta_1"))
        theta_2 = float(((theta_pair_obj.get("outputs") or {}).get("pair2") or {}).get("theta_2"))
    except Exception as exc:
        raise SystemExit(
            json.dumps(
                {
                    "status": "INVALID_INPUT_SHAPE",
                    "expected": "F451 object with inputs.sigma_int.value and outputs.pair{1,2}.{theta_i,u_i}",
                    "error": str(exc),
                },
                ensure_ascii=True,
            )
        )

    n = 12
    I = np.eye(n, dtype=float)

    c1, s1 = real_fourier_pair_basis(n=n, m=1)
    c2, s2 = real_fourier_pair_basis(n=n, m=2)
    C1 = np.column_stack([c1, s1])
    C2 = np.column_stack([c2, s2])

    # Pair-plane projectors and complement projector.
    pi1 = C1 @ C1.T
    pi2 = C2 @ C2.T
    pi_rest = I - pi1 - pi2

    G12 = rotation_so2(alpha_mod_2pi)
    G21 = rotation_so2(-alpha_mod_2pi)

    # Explicit orthogonal swap-with-rotation transport on V1⊕V2, identity on orthogonal complement.
    O12 = (C2 @ G12 @ C1.T) + (C1 @ G21 @ C2.T) + pi_rest

    # Checks / audits.
    orth_res = max_abs(O12.T @ O12 - I)
    invol_res = max_abs(O12 @ O12 - I)
    u1_to_u2_res = max_abs(O12 @ u1 - u2)
    u2_to_u1_res = max_abs(O12.T @ u2 - u1)

    Pu1 = np.outer(u1, u1)
    Pu2 = np.outer(u2, u2)
    proj_transport_res = max_abs(Pu2 - (O12 @ Pu1 @ O12.T))

    artifact = {
        "object": "O12_pair1_pair2_selector_chart_transport_operator_strict_derived_from_sigma_int_alpha12_v1",
        "stage": "F461",
        "status": "actual_exported_lane_scoped_pair_chart_transport_operator__sigma_int_corridor_pair1_pair2__no_false_pass",
        "as_of": "2026-03-15",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "goal": "export_an_explicit_lane_scoped_pair1_pair2_selector_chart_transport_operator_on_n12_fourier_carrier_derived_only_from_sigma_int_theta_pair_transition_angle_alpha12",
        "scope": {
            "n": n,
            "charts": ["pair1", "pair2"],
            "lane": "sigma_int_corridor",
            "global_selector_atlas": False,
        },
        "inputs": {
            "alpha12_object_ref": str(IN_ALPHA12.relative_to(REPO)),
            "theta_pair_object_ref": str(theta_pair_path.relative_to(REPO)),
            "sigma_int": sigma_int,
            "theta_1": theta_1,
            "theta_2": theta_2,
            "note": "Uses the declared real Fourier scaffold for n=12 and the exported alpha12 from the strict sigma-int theta supply; does not introduce any new selector slots.",
        },
        "construction": {
            "fourier_pair_planes": {
                "pair1": "V1 := span{c1,s1} (m=1) with canonical real Fourier basis on Z12",
                "pair2": "V2 := span{c2,s2} (m=2) with canonical real Fourier basis on Z12",
            },
            "pair_plane_projectors": "Π1 := C1 C1^T, Π2 := C2 C2^T, Π_rest := I - Π1 - Π2",
            "so2_rotation": "G(α) := [[cos α, -sin α],[sin α, cos α]] with α := alpha12_mod_2pi",
            "transport_operator": "O12 := C2 G(α) C1^T + C1 G(-α) C2^T + Π_rest",
            "projector_transport": "P(u2) ≈ O12 P(u1) O12^T with P(u):=u u^T (sign-gauge-invariant)",
        },
        "outputs": {
            "n": n,
            "alpha12_mod_2pi": alpha_mod_2pi,
            "G12_so2": [[float(x) for x in row] for row in G12.tolist()],
            "O12": [[float(x) for x in row] for row in O12.tolist()],
            "checks": {
                "orthogonality_max_abs_residual": orth_res,
                "involution_O12_squared_equals_I_max_abs_residual": invol_res,
                "u1_to_u2_max_abs_residual": u1_to_u2_res,
                "u2_to_u1_max_abs_residual": u2_to_u1_res,
                "projector_transport_max_abs_residual": proj_transport_res,
            },
        },
        "hard_limits": [
            "Lane-scoped: exports only a pair1<->pair2 chart-transport operator on the declared n=12 Fourier carrier.",
            "Does not export a global selector atlas, overlap-domain declaration, nor cocycle/gluing data (H41 remains open).",
            "Does not discharge QW-2191 nor export a global selector transition/gluing object (H40 remains globally open).",
            "Does not derive a sign-sensitive physical orientation datum; only projector-level transport is sign-gauge-invariant.",
            "Does not assign a theorem-level physical interpretation to alpha12.",
            "Does not claim strict-core selector closure / admissible S_sel_int.",
            "Does not claim ToE closure.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": artifact["stage"],
        "status": "F461_EXECUTED_CURRENT_STRICT_SIGMA_INT_PAIR1_PAIR2_SELECTOR_CHART_TRANSITION_OBJECT_EXPORT_PACKET_NO_FALSE_PASS",
        "generated_utc": artifact["generated_utc"],
        "outputs": {
            "alpha12_mod_2pi": artifact["outputs"]["alpha12_mod_2pi"],
            "orthogonality_max_abs_residual": artifact["outputs"]["checks"]["orthogonality_max_abs_residual"],
            "involution_max_abs_residual": artifact["outputs"]["checks"]["involution_O12_squared_equals_I_max_abs_residual"],
            "u1_to_u2_max_abs_residual": artifact["outputs"]["checks"]["u1_to_u2_max_abs_residual"],
            "projector_transport_max_abs_residual": artifact["outputs"]["checks"]["projector_transport_max_abs_residual"],
        },
        "no_false_pass": True,
    }

    OUT_OBJECT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()

