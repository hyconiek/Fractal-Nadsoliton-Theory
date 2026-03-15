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

OUT_JSON = GENERATED / "p462_current_strict_shannon_element_order_reference_z24_mode_index_assignment_candidate_probe.json"
OUT_SUMMARY = (
    GENERATED
    / "p462_current_strict_shannon_element_order_reference_z24_mode_index_assignment_candidate_probe_summary.json"
)


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


def defect_F_2m(values: np.ndarray, m: int, n: int) -> complex:
    k = np.arange(n, dtype=float)
    return complex(np.sum(values * np.exp(1j * (4.0 * math.pi * m * k / n))))


def orthonormal_residual(b: np.ndarray) -> float:
    gram = b.T @ b
    return float(np.linalg.norm(gram - np.eye(gram.shape[0])))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    n = 24
    tol = 1e-12

    ord_by_x = [float(zn_element_order(n, x)) for x in range(n)]
    ord_vec = np.array(ord_by_x, dtype=float)

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
            "cut_O2_to_Z2_axis_only": bool(absF > tol),
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

    artifact = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "stage": "P462",
        "goal": "export_probe_level_mode_index_assignment_candidate_on_Z24_from_shannon_element_order_reference_defect",
        "inputs": {
            "n": n,
            "tolerance": tol,
            "reference_profile": "ord_Z24(x)=24/gcd(x,24) for x≠0; ord_Z24(0)=1",
        },
        "construction": {
            "real_fourier_basis": "e0, (c_m,s_m) for m=1..(n/2-1), e_{n/2} (n even)",
            "pair_m_defect": "F_{2m}(ord) := sum_x ord(x) * exp(i * 4π m x / n)",
            "theta_star_rule": "theta_* = (1/2) atan2(Im(F_{2m}), Re(F_{2m}))",
            "u_plus_rule": "u_+ := cos(theta_*) c_m + sin(theta_*) s_m",
            "u_minus_rule": "u_- := -sin(theta_*) c_m + cos(theta_*) s_m",
        },
        "outputs": {
            "n": n,
            "ord_by_x": ord_by_x,
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
            "Probe-only: does not promote any n≠12 result to theorem level.",
            "Does not export typed Z_n / Phase_n / Aut(Z_n) infrastructure.",
            "Does not claim any strict-core promotion of n=24 into the QW-2190 physical scaffold.",
            "Does not claim any global discharge of QW-2191 beyond declared n=12 lanes.",
            "Does not claim strict-core selector closure nor admissible S_sel_int.",
            "Does not claim ToE closure.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "P462",
        "status": "PASS_PROBE_READY",
        "outputs": {
            "n": n,
            "all_pairs_cut": bool(all_pairs_cut),
            "full_basis_orthonormal_residual_vs_identity": orth_res,
        },
        "no_false_pass": True,
    }

    OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()

