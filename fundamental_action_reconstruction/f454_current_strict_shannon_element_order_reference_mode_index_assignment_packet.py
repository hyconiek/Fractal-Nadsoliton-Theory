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

IN_R_ORD = GENERATED / "r_ord_z12_v1_reference_distribution.json"
IN_ALPHA_GEO = GENERATED / "alpha_geo_strict_derived_v1.json"
IN_QW2190 = REPO / "report_qw2190_kernel_mode_representation_emergence_gate.json"

OUT_ASSIGNMENT = GENERATED / "mode_index_assignment_shannon_element_order_reference_strict_core_v1.json"
OUT_SUMMARY = GENERATED / "mode_index_assignment_shannon_element_order_reference_strict_core_v1_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def real_fourier_basis(n: int) -> dict[str, np.ndarray]:
    x = np.arange(n, dtype=float)
    basis: dict[str, np.ndarray] = {"e0": np.ones(n, dtype=float) / math.sqrt(n)}
    for m in range(1, n // 2):
        basis[f"c{m}"] = math.sqrt(2.0 / n) * np.cos(2.0 * math.pi * m * x / n)
        basis[f"s{m}"] = math.sqrt(2.0 / n) * np.sin(2.0 * math.pi * m * x / n)
    if n % 2 == 0:
        basis[f"e{n//2}"] = ((-1.0) ** x) / math.sqrt(n)
    return basis


def orthonormal_residual(b: np.ndarray) -> float:
    gram = b.T @ b
    return float(np.linalg.norm(gram - np.eye(gram.shape[0])))


def defect_F_2m(values: np.ndarray, m: int, n: int) -> complex:
    k = np.arange(n, dtype=float)
    return complex(np.sum(values * np.exp(1j * (4.0 * math.pi * m * k / n))))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    if not IN_R_ORD.exists():
        raise SystemExit(json.dumps({"status": "MISSING_R_ORD_REFERENCE", "missing": str(IN_R_ORD)}, ensure_ascii=True))
    if not IN_ALPHA_GEO.exists():
        raise SystemExit(json.dumps({"status": "MISSING_ALPHA_GEO", "missing": str(IN_ALPHA_GEO)}, ensure_ascii=True))

    r_ord = load_json(IN_R_ORD)
    alpha_geo = load_json(IN_ALPHA_GEO)
    q2190 = load_json(IN_QW2190)

    n = int(q2190["mode_mapping"]["n_octaves"])
    if n != 12:
        raise SystemExit(json.dumps({"status": "UNSUPPORTED_N_EXPECTED_12", "n": n}, ensure_ascii=True))

    ord_by_x = [float(x) for x in r_ord["ord_z12_by_x"]]
    ord_vec = np.array(ord_by_x, dtype=float)
    if ord_vec.shape != (n,):
        raise SystemExit(
            json.dumps({"status": "BAD_R_ORD_SHAPE", "expected_len": n, "actual_len": int(ord_vec.shape[0])}, ensure_ascii=True)
        )

    basis = real_fourier_basis(n)
    e0 = basis["e0"]
    e_half = basis.get("e6")

    pairs: dict[str, Any] = {}
    full_cols: list[np.ndarray] = [e0]

    tol = 1e-12
    all_pairs_cut = True

    for m in range(1, n // 2):
        c = basis[f"c{m}"]
        s = basis[f"s{m}"]

        F = defect_F_2m(ord_vec, m, n)
        theta_star = 0.5 * math.atan2(float(np.imag(F)), float(np.real(F)))

        u_plus = math.cos(theta_star) * c + math.sin(theta_star) * s
        u_minus = -math.sin(theta_star) * c + math.cos(theta_star) * s

        lam_plus = float(np.dot(u_plus, ord_vec * u_plus))
        lam_minus = float(np.dot(u_minus, ord_vec * u_minus))

        if abs(F) <= tol:
            all_pairs_cut = False

        pairs[f"pair{m}"] = {
            "m": m,
            "F_2m_ord": {"Re": float(np.real(F)), "Im": float(np.imag(F)), "abs": float(abs(F))},
            "theta_star": float(theta_star),
            "theta_minimizer_mod_pi": float((theta_star + math.pi / 2.0) % math.pi),
            "u_plus": [float(x) for x in u_plus.tolist()],
            "u_minus": [float(x) for x in u_minus.tolist()],
            "eigenvalues_on_diag_ord_restriction": {"lambda_plus": lam_plus, "lambda_minus": lam_minus},
            "objective_minimizer_vector": "u_minus",
        }

        full_cols.append(u_plus)
        full_cols.append(u_minus)

    if e_half is not None:
        full_cols.append(e_half)

    full_basis = np.column_stack(full_cols)
    orth_res = orthonormal_residual(full_basis)

    assignment = {
        "object": "ModeIndexAssignment_shannon_element_order_reference_strict_core_v1",
        "status": "actual_exported_strict_core_mode_index_assignment__shannon_element_order_reference_lane__no_false_pass",
        "as_of": "2026-03-15",
        "intent": (
            "Export one strict-core mode-index assignment basis object on the n=12 QW-2190 Fourier scaffold by using the "
            "strict Shannon element-order reference datum (r_ord) as an internal symmetry-breaking selector ingredient. "
            "This cuts each pair_m O(2) family down to residual Z2 (sign) and exports an explicit orthonormal basis "
            "(u_plus,u_minus) on every degenerate pair plane."
        ),
        "scope": "strict_core_shannon_element_order_reference_lane_only",
        "derived_from": {
            "r_ord_reference": str(IN_R_ORD.relative_to(REPO)),
            "alpha_geo_strict_derived_v1": str(IN_ALPHA_GEO.relative_to(REPO)),
            "qw2190_mode_scaffold": str(IN_QW2190.relative_to(REPO)),
        },
        "construction": {
            "diagonal_reference_profile": "ord_Z12(x) (stored in r_ord reference packet; direction-free under Aut(Z_12))",
            "pair_m_defect": "F_{2m}(ord) := sum_x ord(x) * exp(i * 4π m x / n)",
            "theta_star_rule": "theta_* = (1/2) atan2(Im(F_{2m}), Re(F_{2m}))",
            "u_plus_rule": "u_+ := cos(theta_*) c_m + sin(theta_*) s_m",
            "u_minus_rule": "u_- := -sin(theta_*) c_m + cos(theta_*) s_m",
            "objective_note": "For cross-entropy J(theta)=E_p[alpha_geo*ord(x)]+const, the minimizer is u_- (smaller eigenvalue).",
            "references": ["F446", "N479", "N480", "N488", "N496"],
        },
        "outputs": {
            "n": n,
            "basis_vectors_order": (
                ["e0", *[f"pair{m}:u_plus,u_minus" for m in range(1, n // 2)], (f"e{n//2}" if (e_half is not None) else None)]
            ),
            "e0": [float(x) for x in e0.tolist()],
            (f"e{n//2}" if e_half is not None else "e_half_absent"): (
                [float(x) for x in e_half.tolist()] if e_half is not None else None
            ),
            "pairs": pairs,
            "all_pairs_cut": bool(all_pairs_cut),
        },
        "audits": {
            "full_basis_shape": [int(full_basis.shape[0]), int(full_basis.shape[1])],
            "full_basis_orthonormal_residual_vs_identity": orth_res,
            "full_basis_det_gram": float(np.linalg.det(full_basis.T @ full_basis)),
        },
        "hard_limits": [
            "Does not claim strict-core selector closure nor admissible S_sel_int.",
            "Does not claim global discharge of QW-2191 (kernel-alone obstruction remains true as a canonical-representative obstruction).",
            "Does not claim ToE closure.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "F454",
        "status": "F454_EXECUTED_CURRENT_STRICT_SHANNON_ELEMENT_ORDER_REFERENCE_MODE_INDEX_ASSIGNMENT_PACKET_NO_FALSE_PASS",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "outputs": {
            "n": n,
            "all_pairs_cut": bool(all_pairs_cut),
            "full_basis_orthonormal_residual_vs_identity": orth_res,
        },
        "inputs": {
            "r_ord_reference": str(IN_R_ORD.relative_to(REPO)),
            "alpha_geo": str(IN_ALPHA_GEO.relative_to(REPO)),
        },
        "no_false_pass": True,
    }

    OUT_ASSIGNMENT.write_text(json.dumps(assignment, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()

