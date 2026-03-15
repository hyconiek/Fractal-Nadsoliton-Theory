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

IN_ALPHA_GEO = GENERATED / "alpha_geo_strict_derived_v1.json"
IN_I24 = GENERATED / "i_24_v1_index_set.json"
IN_Z24 = GENERATED / "z_24_v1_group.json"
IN_TAU_Z24 = GENERATED / "tau_z24_v1_regular_action_on_i_24_v1.json"

OUT_R_ORD = GENERATED / "r_ord_z24_v1_reference_distribution.json"
OUT_ASSIGNMENT = GENERATED / "mode_index_assignment_shannon_element_order_reference_z24_strict_core_v1.json"
OUT_SUMMARY = GENERATED / "mode_index_assignment_shannon_element_order_reference_z24_strict_core_v1_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def zn_element_order(n: int, k: int) -> int:
    kk = int(k) % int(n)
    if kk == 0:
        return 1
    return int(n) // math.gcd(kk, int(n))


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

    required = (IN_ALPHA_GEO, IN_I24, IN_Z24, IN_TAU_Z24)
    missing = [str(p.relative_to(REPO)) for p in required if not p.exists()]
    if missing:
        summary = {
            "stage": "F468",
            "status": "NOT_COMPUTABLE_MISSING_INPUTS",
            "missing": missing,
            "no_false_pass": True,
        }
        OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    alpha_geo = load_json(IN_ALPHA_GEO)
    _i24 = load_json(IN_I24)
    _z24 = load_json(IN_Z24)
    _tau = load_json(IN_TAU_Z24)

    n = 24
    tol = 1e-12

    ord_by_x = [int(zn_element_order(n, x)) for x in range(n)]
    ord_vec = np.array([float(v) for v in ord_by_x], dtype=float)

    r_ord = {
        "object_id": "r_ord_z24_v1_reference_distribution",
        "type": "z24_reference_distribution_v1",
        "definition": "r_ord(x) ∝ exp(-alpha_geo * ord_Z24(x))",
        "alpha_geo": {
            "source_object_id": str(alpha_geo.get("object_id") or "alpha_geo_strict_derived_v1"),
            "value_symbolic": "4 ln 2",
        },
        "carrier": {
            "group_object_id": "z_24_v1_group",
            "index_set_object_id": "i_24_v1_index_set",
        },
        "ord_z24_by_x": ord_by_x,
        "normalization": {
            "symbolic": "Z = Σ_{x∈Z_24} exp(-alpha_geo*ord_z24[x])",
            "notes": [
                "This artifact intentionally stores ord_z24 values; numeric normalization can be computed downstream.",
            ],
        },
        "invariance_notes": [
            "Aut(Z_24)-invariant reference shape (N503): no marked generator/direction.",
            "Not translation-invariant on the regular action: it distinguishes the identity orbit {0}.",
        ],
        "provenance": {"packet": "F468", "theorems": ["N503"]},
    }

    basis = real_fourier_basis(n)
    e0 = basis["e0"]
    e_half = basis.get(f"e{n//2}")

    pairs: dict[str, Any] = {}
    full_cols: list[np.ndarray] = [e0]

    all_pairs_cut = True
    for m in range(1, n // 2):
        c = basis[f"c{m}"]
        s = basis[f"s{m}"]

        F = defect_F_2m(ord_vec, m=m, n=n)
        absF = float(abs(F))
        if absF <= tol:
            all_pairs_cut = False

        theta_star = 0.5 * math.atan2(float(np.imag(F)), float(np.real(F)))

        u_plus = math.cos(theta_star) * c + math.sin(theta_star) * s
        u_minus = -math.sin(theta_star) * c + math.cos(theta_star) * s

        lam_plus = float(np.dot(u_plus, ord_vec * u_plus))
        lam_minus = float(np.dot(u_minus, ord_vec * u_minus))

        pairs[f"pair{m}"] = {
            "m": m,
            "F_2m_ord": {"Re": float(np.real(F)), "Im": float(np.imag(F)), "abs": absF},
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
        "object": "ModeIndexAssignment_shannon_element_order_reference_z24_strict_core_v1",
        "status": "actual_exported_strict_core_mode_index_assignment__shannon_element_order_reference_lane__z24_scope_extension__no_false_pass",
        "as_of": AS_OF,
        "intent": (
            "Export one strict, typed mode-index assignment basis object on Z_24 by using the Shannon element-order reference "
            "shape ord_{Z_24}(x) as an internal symmetry-breaking ingredient. This is a cautious scope-extension infrastructure "
            "export only: it does not promote n=24 into the QW-2190 physical scaffold."
        ),
        "scope": "strict_core_shannon_element_order_reference_lane__z24_scope_extension_only",
        "derived_from": {
            "alpha_geo_strict_derived_v1": str(IN_ALPHA_GEO.relative_to(REPO)),
            "z24_carrier": str(IN_Z24.relative_to(REPO)),
            "i24_index_set": str(IN_I24.relative_to(REPO)),
            "tau_z24_regular_action": str(IN_TAU_Z24.relative_to(REPO)),
            "r_ord_z24_reference": str(OUT_R_ORD.relative_to(REPO)),
        },
        "construction": {
            "diagonal_reference_profile": "ord_Z24(x) (Aut(Z_24)-invariant; no marked direction by N503)",
            "pair_m_defect": "F_{2m}(ord) := sum_x ord(x) * exp(i * 4π m x / n)",
            "theta_star_rule": "theta_* = (1/2) atan2(Im(F_{2m}), Re(F_{2m}))",
            "u_plus_rule": "u_+ := cos(theta_*) c_m + sin(theta_*) s_m",
            "u_minus_rule": "u_- := -sin(theta_*) c_m + cos(theta_*) s_m",
            "objective_note": "For cross-entropy J(theta)=E_p[alpha_geo*ord(x)]+const, the minimizer is u_- (smaller eigenvalue).",
            "references": ["F458", "N503", "P462"],
        },
        "outputs": {
            "n": n,
            "basis_vectors_order": [
                "e0",
                *[f"pair{m}:u_plus,u_minus" for m in range(1, n // 2)],
                f"e{n//2}",
            ],
            "e0": [float(x) for x in e0.tolist()],
            f"e{n//2}": ([float(x) for x in e_half.tolist()] if e_half is not None else None),
            "pairs": pairs,
            "all_pairs_cut": bool(all_pairs_cut),
        },
        "audits": {
            "full_basis_shape": [int(full_basis.shape[0]), int(full_basis.shape[1])],
            "full_basis_orthonormal_residual_vs_identity": orth_res,
            "full_basis_det_gram": float(np.linalg.det(full_basis.T @ full_basis)),
        },
        "hard_limits": [
            "Does not claim any strict-core promotion of n=24 into the QW-2190 physical scaffold.",
            "Does not claim global discharge of QW-2191 beyond declared n=12 lanes.",
            "Does not claim strict-core selector closure nor admissible S_sel_int.",
            "Does not claim any physical sign-sensitive orientation datum.",
            "Does not claim ToE closure.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "F468",
        "status": "F468_EXECUTED_CURRENT_STRICT_SHANNON_ELEMENT_ORDER_REFERENCE_Z24_MODE_INDEX_ASSIGNMENT_PACKET_NO_FALSE_PASS",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "outputs": {
            "n": n,
            "all_pairs_cut": bool(all_pairs_cut),
            "full_basis_orthonormal_residual_vs_identity": orth_res,
        },
        "inputs": {
            "alpha_geo": str(IN_ALPHA_GEO.relative_to(REPO)),
            "z24_carrier": str(IN_Z24.relative_to(REPO)),
        },
        "no_false_pass": True,
    }

    OUT_R_ORD.write_text(json.dumps(r_ord, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_ASSIGNMENT.write_text(json.dumps(assignment, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()

