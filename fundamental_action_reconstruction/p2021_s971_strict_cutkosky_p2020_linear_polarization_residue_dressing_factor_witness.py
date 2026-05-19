#!/usr/bin/env python3
"""P2021 S971 strict Cutkosky P2020 proxy residue admissibility audit.

Corrected next honest step after review: P1994's Z(s) is only a proxy/mirrored
import, so it must not be accepted as a P1953 dressed-residue export.  This
script keeps the useful algebraic sandbox (a positive scalar proxy factor does
not rotate the P2020 {plus,cross} matrix), but the theorem-level verdict is an
admissibility failure for dressed Cutkosky closure until a loop-derived,
same-scheme residue or DiscM_common_basis is exported.
"""
from __future__ import annotations

import json
import platform
from pathlib import Path
from typing import Any

import numpy as np
import scipy.linalg as la
import sympy as sp

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
OUT = GEN / "p2021_s971_strict_cutkosky_p2020_linear_polarization_residue_dressing_factor_witness.json"
TS = "2026-05-18T00:00:00+00:00"
CHANNEL = "graviton->gauge_gauge"


def load(name: str) -> dict[str, Any]:
    path = GEN / name
    if not path.exists():
        return {"_missing": name, "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def sympify_matrix(rows: list[list[str]]) -> sp.Matrix:
    return sp.Matrix([[sp.sympify(item, locals={"pi": sp.pi}) for item in row] for row in rows])


def matrix_to_strings(matrix: sp.Matrix) -> list[list[str]]:
    return [
        [str(sp.factor(sp.together(sp.simplify(matrix[i, j])))) for j in range(matrix.cols)]
        for i in range(matrix.rows)
    ]


def is_zero_matrix(matrix: sp.Matrix) -> bool:
    return all(sp.simplify(entry) == 0 for entry in matrix)


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2020 = load("p2020_s970_strict_cutkosky_p2019_tree_phase_space_cut_sum_witness.json")
    p1994 = load("p1994_s944_strict_cutkosky_dressed_amplitude_first_import_witness.json")

    s = sp.symbols("s", positive=True, real=True)
    p2020_exact = p2020.get("exact_phase_space_cut_sum", {})
    tree_no_sym = sympify_matrix(
        p2020_exact.get(
            "linear_polarization_resolved_CutSum_no_identical_symmetry_over_kappa2_Zgauge2",
            [["0", "0"], ["0", "0"]],
        )
    )
    tree_identical = sympify_matrix(
        p2020_exact.get(
            "linear_polarization_resolved_CutSum_identical_final_state_over_kappa2_Zgauge2",
            [["0", "0"], ["0", "0"]],
        )
    )

    z_expr_text = ((p1994.get("dressed_import") or {}).get("Z(s)") or "0")
    z_factor = sp.factor(sp.sympify(z_expr_text, locals={"s": s, "pi": sp.pi, "log": sp.log, "ln": sp.log}))
    proxy_residue_factor = sp.factor(sp.together(z_factor**2))

    proxy_no_sym = sp.simplify(proxy_residue_factor * tree_no_sym)
    proxy_identical = sp.simplify(proxy_residue_factor * tree_identical)
    offdiag_no_sym = sp.Matrix([[0, proxy_no_sym[0, 1]], [proxy_no_sym[1, 0], 0]])
    offdiag_identical = sp.Matrix([[0, proxy_identical[0, 1]], [proxy_identical[1, 0], 0]])
    eigen_no_sym = [sp.factor(sp.together(ev)) for ev in proxy_no_sym.eigenvals().keys()]
    eigen_identical = [sp.factor(sp.together(ev)) for ev in proxy_identical.eigenvals().keys()]

    numerator, denominator = sp.fraction(proxy_residue_factor)
    proxy_positive_certificate = {
        "domain": "s > 0 real",
        "Z_proxy(s)": str(z_factor),
        "R_proxy(s)=Z_proxy(s)^2": str(proxy_residue_factor),
        "numerator": str(sp.factor(numerator)),
        "denominator": str(sp.factor(denominator)),
        "algebraic_reason": "R_proxy(s) is the square of a real rational function; for P1994 z1=1/50, Z_proxy(s)=(51*s + 50)/(50*(s + 1)) is strictly positive for s>0.",
        "admissibility_warning": "This positivity is only a scalar proxy-transport fact.  It is not a loop-derived dressed residue and is not accepted as P1953 M_dressed_common_basis or DiscM_common_basis input.",
    }

    grid = [sp.Rational(1, 2), sp.Integer(1), sp.Integer(2), sp.Integer(4), sp.Integer(8)]
    rows: list[dict[str, Any]] = []
    max_exact_numeric_error = 0.0
    for sv in grid:
        exact_matrix = sp.N(proxy_no_sym.subs(s, sv), 50)
        numeric_matrix = np.array(exact_matrix.tolist(), dtype=float)
        eig_numeric = np.linalg.eigvalsh(numeric_matrix)
        eig_exact = np.array([float(sp.N(ev.subs(s, sv), 50)) for ev in eigen_no_sym], dtype=float)
        eig_error = float(la.norm(np.sort(eig_numeric) - np.sort(eig_exact), 2))
        max_exact_numeric_error = max(max_exact_numeric_error, eig_error)
        rows.append(
            {
                "s": str(sv),
                "Z_proxy(s)": float(sp.N(z_factor.subs(s, sv), 50)),
                "R_proxy(s)=Z_proxy(s)^2": float(sp.N(proxy_residue_factor.subs(s, sv), 50)),
                "proxy_transport_no_symmetry_matrix_over_kappa2_Zgauge2": numeric_matrix.tolist(),
                "proxy_transport_no_symmetry_eigvals_numeric": eig_numeric.tolist(),
                "proxy_transport_no_symmetry_trace_numeric": float(np.trace(numeric_matrix)),
                "positive_semidefinite_numeric": bool(np.all(eig_numeric > 0.0)),
                "eig_numeric_vs_exact_l2_error": eig_error,
            }
        )

    tree_input_ok = p2020.get("result_kind") == "PASS_TREE_PHASE_SPACE_LINEAR_POLARIZATION_CUT_SUM_COMPONENT_WITNESS"
    p1994_is_proxy = p1994.get("result_kind") == "PASS_DRESSED_IMPORT_OPTICAL_ZERO_PROXY_WITNESS"
    expected_tree_no_sym = sp.eye(2) / sp.pi
    expected_tree_identical = sp.eye(2) / (2 * sp.pi)
    matrices_match_p2020 = is_zero_matrix(tree_no_sym - expected_tree_no_sym) and is_zero_matrix(tree_identical - expected_tree_identical)
    proxy_factor_positive_symbolic = sp.factor(z_factor - (51 * s + 50) / (50 * (s + 1))) == 0
    no_basis_rotation = is_zero_matrix(offdiag_no_sym) and is_zero_matrix(offdiag_identical)
    psd_symbolic = all(
        sp.factor(sp.together(ev / proxy_residue_factor)) in {sp.Rational(1, 1) / sp.pi, sp.Rational(1, 2) / sp.pi}
        for ev in eigen_no_sym + eigen_identical
    )
    numeric_grid_positive = all(row["positive_semidefinite_numeric"] for row in rows)

    local_algebra_sanity_checks = {
        "p2020_linear_polarization_matrix_available": tree_input_ok,
        "p1994_proxy_factor_available": p1994_is_proxy,
        "p2020_matrix_normalization_matches_expected": matrices_match_p2020,
        "proxy_factor_symbolically_positive_on_s_positive": bool(proxy_factor_positive_symbolic),
        "linear_polarization_basis_preserved_no_rotation": bool(no_basis_rotation),
        "proxy_transport_matrices_symbolically_psd_given_R_positive": bool(psd_symbolic),
        "numeric_grid_psd_and_exact_eigenvalue_match": bool(numeric_grid_positive and max_exact_numeric_error < 1e-14),
    }

    admissibility_checks = {
        "loop_derived_from_L_total": False,
        "same_scheme_as_P1953_MSbar_B1_seed_locked": False,
        "DiscM_common_basis_exported": False,
        "DiscM_minus_CutSum_simplified_evaluated": False,
        "BRST_physical_state_projector_exported": False,
        "proxy_not_promoted_to_dressed_residue": True,
    }
    admissible_for_p1953 = all(admissibility_checks.values())

    out = {
        "ledger_id": "P2021_S971_STRICT_CUTKOSKY_P2020_PROXY_RESIDUE_ADMISSIBILITY_AUDIT",
        "packet_id": "P2021",
        "stage_id": "S971",
        "produced_by": Path(__file__).name,
        "timestamp_utc": TS,
        "route": "strict_only",
        "legacy_bridge_used": False,
        "channel": CHANNEL,
        "depends_on": {
            "p2020": p2020.get("ledger_id"),
            "p1994": p1994.get("ledger_id"),
        },
        "corrective_review_note": "The P1994 Z(s) factor is retained only as an algebraic proxy sandbox.  P2021 now rejects it as a theorem-grade dressed-residue input because it is not loop-derived, not same-scheme locked to P1953, and not tied to BRST-projected DiscM.",
        "input_scope": {
            "tree_cut_sum_source": "P2020 exact linear-polarization-resolved tree phase-space CutSum matrix",
            "proxy_factor_source": "P1994 first dressed proxy import Z(s); P1994 itself labels the construction proxy and mirrored, not a loop-derived backend",
            "normalization": "matrices are in units of kappa^2*Z_gauge^2 and preserve the P2020 {plus,cross} real linear-polarization basis",
            "scheme_guard": "No MSbar_B1_seed loop discontinuity is locked to this proxy factor.",
        },
        "p2020_tree_linear_polarization_matrices_over_kappa2_Zgauge2": {
            "no_identical_symmetry": matrix_to_strings(tree_no_sym),
            "identical_final_state": matrix_to_strings(tree_identical),
            "basis_order": ["plus", "cross"],
        },
        "proxy_residue_factor_sandbox": proxy_positive_certificate,
        "proxy_transport_linear_polarization_matrices_over_kappa2_Zgauge2": {
            "no_identical_symmetry": matrix_to_strings(proxy_no_sym),
            "identical_final_state": matrix_to_strings(proxy_identical),
            "basis_order": ["plus", "cross"],
            "trace_no_identical_symmetry": str(sp.factor(sp.together(sp.trace(proxy_no_sym)))),
            "trace_identical_final_state": str(sp.factor(sp.together(sp.trace(proxy_identical)))),
            "eigenvalues_no_identical_symmetry": [str(ev) for ev in eigen_no_sym],
            "eigenvalues_identical_final_state": [str(ev) for ev in eigen_identical],
            "accepted_as_dressed_residue": False,
        },
        "numeric_grid_certificate": {
            "grid_s_values": [str(v) for v in grid],
            "max_exact_numeric_eigenvalue_l2_error": max_exact_numeric_error,
            "rows": rows,
            "interpretation": "Numerical PSD confirms only the proxy scalar-transport algebra, not theorem-grade unitarity.",
        },
        "p1953_contract_non_update": {
            "M_dressed_common_basis": "OPEN_NOT_EXPORTED__P2021_PROXY_FACTOR_REJECTED_AS_DRESSED_INPUT",
            "AbsM_dressed_squared_common_basis": "OPEN_NOT_EXPORTED__PROXY_TRANSPORT_SANDBOX_ONLY",
            "CutSum_common_basis": "P2020_TREE_LINEAR_POLARIZATION_CUTSUM_REMAINS_AVAILABLE_BUT_NOT_DRESSED",
            "DiscM_common_basis": "OPEN_NOT_EVALUATED",
            "DiscM_minus_CutSum_simplified": "OPEN_NOT_EVALUATED",
            "external_state_projectors including BRST physical-state projector": "OPEN_BRST_COHOMOLOGY_PROJECTOR_NOT_EXPORTED",
        },
        "local_algebra_sanity_checks": local_algebra_sanity_checks,
        "admissibility_checks_for_dressed_cutkosky": admissibility_checks,
        "admissible_for_p1953_dressed_interface": admissible_for_p1953,
        "result_kind": "OPEN_PROXY_RESIDUE_TRANSPORT_SANITY_ONLY_NOT_P1953_ADMISSIBLE",
        "status": "OPEN_OBSTRUCTION_WITH_TRACE",
        "false_pass_guard": "P2021 proves only that a scalar positive proxy factor would preserve the P2020 linear-polarization matrix algebra.  It explicitly rejects that proxy as a loop-derived dressed residue, does not update M_dressed_common_basis, does not export DiscM_common_basis, does not evaluate DiscM=CutSum, does not supply BRST projection, does not discharge QW-2191, and does not close ToE.",
        "next_honest_step": "Build P2022 as a real same-scheme source audit/derivation: either extract a loop-derived hAA self-energy/discontinuity residue from L_total with BRST projector data, or export a nonavailability theorem listing the exact missing vertices, gauges, and projector inputs.",
        "toe_progress": "Corrects the unitarity route by preventing a proxy residue factor from being mistaken for dressed Cutkosky data.  The useful local algebra survives as a sandbox, but the theorem frontier is sharpened to the missing loop-derived same-scheme DiscM/residue and BRST projector.",
        "lay_explanation": "Sprawdziliśmy, że dodatni próbny mnożnik nie psuje macierzy dwóch polaryzacji. Ale po recenzji nie wolno go uznać za prawdziwą poprawkę pętlową: nie pochodzi z pełnego rachunku pętli i nie ma projektora BRST. Dlatego wynik jest audytem blokady, nie krokiem domykającym unitarność.",
        "environment": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "scipy": __import__("scipy").__version__,
            "sympy": sp.__version__,
        },
    }

    OUT.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"[P2021] wrote proxy residue admissibility audit: {OUT}")


if __name__ == "__main__":
    main()
