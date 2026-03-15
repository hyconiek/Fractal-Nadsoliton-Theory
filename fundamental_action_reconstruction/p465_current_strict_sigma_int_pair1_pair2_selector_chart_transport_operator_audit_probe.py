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

IN_F461 = GENERATED / "o12_pair1_pair2_selector_chart_transport_operator_strict_derived_from_sigma_int_alpha12_v1.json"

OUT_JSON = GENERATED / "p465_current_strict_sigma_int_pair1_pair2_selector_chart_transport_operator_audit_probe.json"
OUT_SUMMARY = (
    GENERATED / "p465_current_strict_sigma_int_pair1_pair2_selector_chart_transport_operator_audit_probe_summary.json"
)


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def max_abs(a: np.ndarray) -> float:
    return float(np.max(np.abs(a)))


def as_square_matrix(x: Any, n: int, label: str) -> np.ndarray:
    if not (isinstance(x, list) and len(x) == n and all(isinstance(row, list) and len(row) == n for row in x)):
        raise ValueError(f"{label} must be a {n}x{n} nested list")
    a = np.array(x, dtype=float)
    if not np.isfinite(a).all():
        raise ValueError(f"{label} must contain only finite numbers")
    return a


def real_fourier_pair_basis(n: int, m: int) -> tuple[np.ndarray, np.ndarray]:
    scale = math.sqrt(2.0 / float(n))
    c = np.array([scale * math.cos(2.0 * math.pi * m * k / float(n)) for k in range(n)], dtype=float)
    s = np.array([scale * math.sin(2.0 * math.pi * m * k / float(n)) for k in range(n)], dtype=float)
    return c, s


def rotation_so2(alpha: float) -> np.ndarray:
    return np.array([[math.cos(alpha), -math.sin(alpha)], [math.sin(alpha), math.cos(alpha)]], dtype=float)


def u_theta(theta: float, c: np.ndarray, s: np.ndarray) -> np.ndarray:
    return (math.cos(theta) * c) + (math.sin(theta) * s)


def tau_theta(theta: float, c: np.ndarray, s: np.ndarray) -> np.ndarray:
    return (-math.sin(theta) * c) + (math.cos(theta) * s)


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    if not IN_F461.exists():
        raise SystemExit(
            json.dumps(
                {
                    "status": "NOT_COMPUTABLE_MISSING_INPUT",
                    "missing": str(IN_F461.relative_to(REPO)),
                    "expected": "F461 chart transport operator export object",
                },
                ensure_ascii=True,
            )
        )

    f461 = load_json(IN_F461)
    try:
        outputs = f461.get("outputs") or {}
        n = int(outputs.get("n"))
        if n != 12:
            raise ValueError("expected n=12")
        alpha12 = float(outputs.get("alpha12_mod_2pi"))
        O12 = as_square_matrix(outputs.get("O12"), n=n, label="F461.outputs.O12")
        G12 = as_square_matrix(outputs.get("G12_so2"), n=2, label="F461.outputs.G12_so2")
    except Exception as exc:
        raise SystemExit(
            json.dumps(
                {
                    "status": "INVALID_F461_SHAPE",
                    "expected": "F461.outputs.{n,alpha12_mod_2pi,O12,G12_so2}",
                    "error": str(exc),
                },
                ensure_ascii=True,
            )
        )

    I = np.eye(n, dtype=float)

    # Canonical Fourier pair-plane bases (declared scaffold).
    c1, s1 = real_fourier_pair_basis(n=n, m=1)
    c2, s2 = real_fourier_pair_basis(n=n, m=2)
    C1 = np.column_stack([c1, s1])
    C2 = np.column_stack([c2, s2])

    pi1 = C1 @ C1.T
    pi2 = C2 @ C2.T
    pi_rest = I - pi1 - pi2

    # Structural checks.
    orth_res = max_abs(O12.T @ O12 - I)
    invol_res = max_abs(O12 @ O12 - I)
    plane1_to_plane2_res = max_abs((O12 @ pi1 @ O12.T) - pi2)
    plane2_to_plane1_res = max_abs((O12 @ pi2 @ O12.T) - pi1)
    rest_fixed_res = max_abs((O12 @ pi_rest @ O12.T) - pi_rest)

    # Basis-level mapping checks implied by alpha12.
    alphaG = rotation_so2(alpha12)
    # Expected: O12 c1 = cos α c2 + sin α s2, O12 s1 = -sin α c2 + cos α s2.
    expected_c1 = (alphaG[0, 0] * c2) + (alphaG[1, 0] * s2)
    expected_s1 = (alphaG[0, 1] * c2) + (alphaG[1, 1] * s2)
    c1_map_res = max_abs(O12 @ c1 - expected_c1)
    s1_map_res = max_abs(O12 @ s1 - expected_s1)

    # Sampled direction transport checks.
    # Choose a finite theta grid; avoid claiming global completeness.
    thetas = [0.0, math.pi / 6.0, math.pi / 4.0, math.pi / 3.0, math.pi / 2.0, -math.pi / 3.0]
    dir_u_res = 0.0
    dir_tau_res = 0.0
    proj_res = 0.0
    for theta in thetas:
        u1 = u_theta(theta, c1, s1)
        t1 = tau_theta(theta, c1, s1)
        u2 = u_theta(theta + alpha12, c2, s2)
        t2 = tau_theta(theta + alpha12, c2, s2)

        dir_u_res = max(dir_u_res, max_abs(O12 @ u1 - u2))
        dir_tau_res = max(dir_tau_res, max_abs(O12 @ t1 - t2))

        Pu1 = np.outer(u1, u1)
        Pu2 = np.outer(u2, u2)
        proj_res = max(proj_res, max_abs(Pu2 - (O12 @ Pu1 @ O12.T)))

    artifact = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "stage": "P465",
        "goal": "audit_lane_scoped_pair1_pair2_selector_chart_transport_operator_properties_for_f461_on_n12_fourier_scaffold",
        "inputs": {"f461_object": str(IN_F461.relative_to(REPO))},
        "construction": {
            "pair_plane_bases": "canonical real Fourier (c1,s1) and (c2,s2) on Z12",
            "theta_grid": thetas,
        },
        "audits": {
            "orthogonality_max_abs_residual": orth_res,
            "involution_max_abs_residual": invol_res,
            "plane_transport": {
                "O12_pi1_O12T_minus_pi2_max_abs": plane1_to_plane2_res,
                "O12_pi2_O12T_minus_pi1_max_abs": plane2_to_plane1_res,
                "O12_pi_rest_O12T_minus_pi_rest_max_abs": rest_fixed_res,
            },
            "basis_mapping_max_abs_residuals": {"c1_map_res": c1_map_res, "s1_map_res": s1_map_res},
            "sampled_direction_transport_max_abs_residuals": {
                "u_theta_transport_max_abs": dir_u_res,
                "tau_theta_transport_max_abs": dir_tau_res,
                "projector_transport_max_abs": proj_res,
            },
            "g12_consistency": {
                "alpha12_mod_2pi": alpha12,
                "g12_from_alpha12": [[float(x) for x in row] for row in alphaG.tolist()],
                "g12_exported_by_f461": [[float(x) for x in row] for row in G12.tolist()],
                "g12_max_abs_diff": max_abs(alphaG - G12),
            },
        },
        "hard_limits": [
            "Lane-scoped audit only: samples a finite theta grid and checks structural properties on n=12.",
            "Does not export a global selector atlas nor global gluing cocycle data.",
            "Does not discharge QW-2191.",
            "Does not derive a sign-sensitive physical orientation datum.",
            "Does not claim strict-core selector closure / admissible S_sel_int.",
            "Does not claim ToE closure.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "P465",
        "status": "PASS_PROBE_READY",
        "generated_utc": artifact["generated_utc"],
        "audits": {
            "orthogonality_max_abs_residual": orth_res,
            "plane_transport_max_abs_residual": max(plane1_to_plane2_res, plane2_to_plane1_res, rest_fixed_res),
            "sampled_u_transport_max_abs_residual": dir_u_res,
            "sampled_projector_transport_max_abs_residual": proj_res,
            "g12_max_abs_diff": artifact["audits"]["g12_consistency"]["g12_max_abs_diff"],
        },
        "no_false_pass": True,
    }

    OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()

