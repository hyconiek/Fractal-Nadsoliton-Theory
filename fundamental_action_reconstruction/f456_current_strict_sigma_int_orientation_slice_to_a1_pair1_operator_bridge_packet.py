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

IN_THETA_PAIR = GENERATED / "theta_pair_sigma_int_strict_selector_ingredient_o2_cut_slot_free_v1.json"
IN_QW2190 = REPO / "report_qw2190_kernel_mode_representation_emergence_gate.json"

OUT_OPERATOR = GENERATED / "a_1_pair1_orientation_projector_operator_strict_core_v1.json"
OUT_SUMMARY = GENERATED / "f456_current_strict_sigma_int_orientation_slice_to_a1_pair1_operator_bridge_packet_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def clean_scalar(value: float) -> float:
    if abs(value) < 1e-15:
        return 0.0
    return round(float(value), 15)


def real_fourier_basis_pair1(n: int) -> tuple[np.ndarray, np.ndarray]:
    x = np.arange(n, dtype=float)
    c1 = math.sqrt(2.0 / n) * np.cos(2.0 * math.pi * 1.0 * x / n)
    s1 = math.sqrt(2.0 / n) * np.sin(2.0 * math.pi * 1.0 * x / n)
    return c1, s1


def max_abs_diff(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.max(np.abs(a - b)))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    missing: list[str] = []
    for p in (IN_THETA_PAIR, IN_QW2190):
        if not p.exists():
            missing.append(str(p.relative_to(REPO)))
    if missing:
        raise SystemExit(json.dumps({"status": "MISSING_INPUT_ARTIFACTS", "missing": missing}, ensure_ascii=True))

    theta_pair = load_json(IN_THETA_PAIR)
    q2190 = load_json(IN_QW2190)

    expected_theta_obj = "ThetaPair_sigma_int_strict_selector_ingredient_o2_cut_slot_free_v1"
    if theta_pair.get("object") != expected_theta_obj:
        raise SystemExit(
            json.dumps(
                {
                    "status": "UNEXPECTED_THETA_PAIR_OBJECT",
                    "expected": expected_theta_obj,
                    "actual": theta_pair.get("object"),
                },
                ensure_ascii=True,
            )
        )

    n = int((q2190.get("mode_mapping") or {}).get("n_octaves") or 0)
    if n != 12:
        raise SystemExit(json.dumps({"status": "UNEXPECTED_N", "expected": 12, "actual": n}, ensure_ascii=True))

    theta_1 = float((((theta_pair.get("outputs") or {}).get("pair1") or {}).get("theta_1")))
    u1_list = (((theta_pair.get("outputs") or {}).get("pair1") or {}).get("u_1"))
    if not (isinstance(u1_list, list) and len(u1_list) == n and all(isinstance(x, (int, float)) for x in u1_list)):
        raise SystemExit(
            json.dumps({"status": "INVALID_U1_VECTOR", "expected": f"numeric list length {n}", "actual": u1_list}, ensure_ascii=True)
        )

    u1 = np.array([float(x) for x in u1_list], dtype=float)

    c1, s1 = real_fourier_basis_pair1(n)
    x1 = float(np.dot(c1, u1))
    y1 = float(np.dot(s1, u1))
    u1_in_plane = x1 * c1 + y1 * s1
    u1_outside = float(np.linalg.norm(u1 - u1_in_plane))

    coords = np.array([x1, y1], dtype=float)
    coords_norm = float(np.linalg.norm(coords))
    if coords_norm == 0.0:
        raise SystemExit(json.dumps({"status": "DEGENERATE_U1_COORDS", "coords_norm": 0.0}, ensure_ascii=True))
    coords_unit = coords / coords_norm

    P = np.outer(coords_unit, coords_unit)  # 2x2 projector in (c1,s1)
    I2 = np.eye(2, dtype=float)
    P_sym_res = float(np.linalg.norm(P - P.T, ord=np.inf))
    P_idem_res = float(np.linalg.norm(P @ P - P, ord=np.inf))
    P_trace = float(np.trace(P))
    P_det = float(np.linalg.det(P))
    P_eigs = np.linalg.eigvalsh((P + P.T) * 0.5)

    P_from_neg_u1 = np.outer(-coords_unit, -coords_unit)
    sign_invariance_res = max_abs_diff(P, P_from_neg_u1)

    tol = 1e-12
    checks = [
        {
            "id": "u1_in_span_c1_s1",
            "pass": bool(u1_outside <= tol),
            "expected": f"<= {tol:g}",
            "actual": u1_outside,
            "meaning": "exported u_1 lies in span{c1,s1} (numerical tolerance)",
        },
        {
            "id": "coords_unit_norm_1",
            "pass": bool(abs(coords_norm - 1.0) <= 1e-9),
            "expected": "~1",
            "actual": coords_norm,
            "meaning": "u_1 coordinates in (c1,s1) have unit norm (numerical tolerance)",
        },
        {
            "id": "projector_symmetric",
            "pass": bool(P_sym_res <= tol),
            "expected": f"<= {tol:g}",
            "actual": P_sym_res,
            "meaning": "A_1(pair1) is symmetric in the (c1,s1) basis",
        },
        {
            "id": "projector_idempotent",
            "pass": bool(P_idem_res <= tol),
            "expected": f"<= {tol:g}",
            "actual": P_idem_res,
            "meaning": "A_1(pair1)^2 = A_1(pair1) (projector)",
        },
        {
            "id": "residual_z2_sign_invariant",
            "pass": bool(sign_invariance_res <= tol),
            "expected": f"<= {tol:g}",
            "actual": sign_invariance_res,
            "meaning": "u_1 -> -u_1 leaves the operator unchanged (residual Z2 sign gauge)",
        },
    ]

    audits_pass = bool(all(bool(c["pass"]) for c in checks))

    operator = {
        "object": "A_1_pair1_orientation_projector_operator_strict_core_v1",
        "status": "actual_exported_operator__derived_only_from_sigma_int_u1__residual_z2_sign_gauge_invariant",
        "as_of": AS_OF,
        "intent": (
            "Export one strict-core, slot-free operator-level bridge target on V_1 = span{c1,s1} constructed only from the "
            "already exported sigma-int orientation direction u_1. The operator is the rank-one projector A_1(pair1) := |u_1><u_1| "
            "in the (c1,s1) basis. This is a minimal downstream operator reachability witness and does not identify with any extension-only "
            "H/O-lane A_1_ext nor with any host matching claim."
        ),
        "inputs": {
            "theta_pair": str(IN_THETA_PAIR.relative_to(REPO)),
            "qw2190_scaffold": str(IN_QW2190.relative_to(REPO)),
        },
        "domain_notes": {
            "basis_plane": "V_1 = span{c1,s1} inside the n=12 real Fourier scaffold",
            "sigma_int_lane": "slot_free_theta_pair_supply (F451/N489)",
        },
        "data": {
            "theta_1_exported": theta_1,
            "u_1": [clean_scalar(float(x)) for x in u1.tolist()],
            "u_1_coords_in_c1_s1": [clean_scalar(float(x)) for x in coords_unit.tolist()],
            "A_1_pair1_matrix_in_c1_s1": [[clean_scalar(float(x)) for x in row] for row in P.tolist()],
        },
        "audits": {
            "u1_outside_span_c1_s1": u1_outside,
            "coords_norm": coords_norm,
            "projector_symmetry_inf_norm": P_sym_res,
            "projector_idempotence_inf_norm": P_idem_res,
            "projector_trace": P_trace,
            "projector_det": P_det,
            "projector_eigenvalues": [clean_scalar(float(x)) for x in P_eigs.tolist()],
            "residual_z2_sign_invariance_max_abs_diff": sign_invariance_res,
            "tolerance": tol,
        },
        "checks": checks,
        "hard_limits": [
            "Does not claim identification with the extension-only H/O-lane A_1_ext operator.",
            "Does not claim any host matching/cancellation, coefficient extraction, or K_obs factorization.",
            "Does not claim strict-core selector closure / admissible S_sel_int.",
            "Does not claim global discharge of QW-2191.",
            "Does not claim ToE closure.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "F456",
        "status": "F456_EXECUTED_CURRENT_STRICT_SIGMA_INT_ORIENTATION_SLICE_TO_A1_PAIR1_OPERATOR_BRIDGE_PACKET_NO_FALSE_PASS",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "constructed_object": str(OUT_OPERATOR.relative_to(REPO)),
        "audits_pass": bool(audits_pass),
        "no_false_pass": True,
    }

    OUT_OPERATOR.write_text(json.dumps(operator, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()

