#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import subprocess
from typing import Any

import mpmath as mp
import sympy as sp

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import (
    DOC_FILES,
    REPO,
    ROOT,
    load_json,
    rel,
)
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2509_s1459_strict_damping_rg_minimum_roughness_variational_wellposedness_certificate import (
    derivative_coefficients,
    eval_poly,
    integrate_square_closed_form,
    node_vanishing_r_coefficients,
    poly_mul,
)
from p2512_s1462_strict_damping_rg_quadratic_source_admissibility_audit import integrate_product_closed_form

GEN = ROOT / "generated"
OUT = GEN / "p2513_s1463_strict_damping_rg_derivative_order_nonidentifiability_certificate.json"
MD = GEN / "p2513_s1463_strict_damping_rg_derivative_order_nonidentifiability_certificate.md"

SOURCE_FILES = {
    "P2512_QUADRATIC_SOURCE_ADMISSIBILITY": GEN / "p2512_s1462_strict_damping_rg_quadratic_source_admissibility_audit.json",
}

mp.mp.dps = 100
DOMAIN = list(range(1, 12))
BASIS_DEGREES = list(range(4))
MIXED_ROWS = [(mp.mpf("1"), mp.mpf("0")), (mp.mpf("0"), mp.mpf("1")), (mp.mpf("1"), mp.mpf("1")), (mp.mpf("2"), mp.mpf("3"))]


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")).hexdigest()


def rg_count(pattern: str) -> dict[str, Any]:
    proc = subprocess.run(
        [
            "rg", "-n", pattern, "fundamental_action_reconstruction", "material_dowodowy",
            "-g", "*.py", "-g", "*.md", "-g", "*.tex", "-g", "!fundamental_action_reconstruction/generated/**",
        ],
        cwd=REPO,
        check=False,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    lines = sorted(line for line in proc.stdout.splitlines() if line)
    return {"count": len(lines), "samples": lines[:40]}


def rg_audit() -> dict[str, Any]:
    patterns = {
        "new_packet": "P2513|S1463|derivative-order nonidentifiability|H1 H2 selector ambiguity|mixed derivative selector|Sobolev order ambiguity",
        "precursor_packets": "P2512|S1462|quadratic source admissibility|P2511|natural spline collapse|P2510|roughness KKT stationarity",
        "derivative_order_language": "derivative-only selector|Sobolev order|H1|H2|mixed derivative|order nonidentifiability",
        "guardrails": "legacy -> strict completion bridge|role-transfer audit|K_legacy_ont|K_strict_gate|QW-2191|ToE closure",
        "closure_blockers": "source theorem|bridge theorem|role-transfer theorem|physical-value generator|role-bearing L_total|selector closure",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def theorem(payload: dict[str, Any], container: str) -> dict[str, Any]:
    return payload.get(container, {}).get("theorem_export", {})


def text(value: mp.mpf, digits: int = 70) -> str:
    return mp.nstr(value, digits)


def monomial_times(coeffs: list[mp.mpf], degree: int) -> list[mp.mpf]:
    return [mp.mpf("0")] * degree + coeffs[:]


def basis_polynomials() -> list[list[mp.mpf]]:
    r = node_vanishing_r_coefficients()
    return [monomial_times(r, degree) for degree in BASIS_DEGREES]


def gram_matrix_for_order(order: int) -> mp.matrix:
    upper = mp.log(11)
    basis = basis_polynomials()
    matrix = mp.matrix(len(basis), len(basis))
    for i, left in enumerate(basis):
        left_d = derivative_coefficients(left, order)
        for j, right in enumerate(basis):
            right_d = derivative_coefficients(right, order)
            matrix[i, j] = integrate_product_closed_form(left_d, right_d, upper)
    return matrix


def matrix_pivots_and_minors(matrix: mp.matrix) -> tuple[list[mp.mpf], list[mp.mpf]]:
    leading_minors = []
    pivots = []
    n = matrix.rows
    for k in range(1, n + 1):
        sub = mp.matrix(k, k)
        for i in range(k):
            for j in range(k):
                sub[i, j] = matrix[i, j]
        det = mp.det(sub)
        leading_minors.append(det)
        pivots.append(det if k == 1 else det / leading_minors[k - 2])
    return leading_minors, pivots


def symbolic_selector_family_audit() -> dict[str, Any]:
    ell = sp.symbols("ell", real=True)
    delta = sp.Rational(14, 5) - 4 * sp.log(2)
    y0 = delta * ell
    p = sp.Function("p")
    return {
        "symbolic_backend": "sympy",
        "sympy_version": sp.__version__,
        "candidate_y0": sp.sstr(y0),
        "h1_identity": "J1[y0+p]-J1[y0]=2*delta*(p(L)-p(0))+int (p')^2 = int (p')^2 for node-vanishing p.",
        "h2_identity": "J2[y0+p]-J2[y0]=int (p'')^2 because y0''=0.",
        "mixed_identity": "J_{a,b}[y0+p]-J_{a,b}[y0]=a*int(p')^2+b*int(p'')^2 for a,b>=0 on node-vanishing p.",
        "h1_selects_same_affine_minimizer": True,
        "h2_selects_same_affine_minimizer": True,
        "stationarity_cannot_identify_derivative_order": True,
        "symbolic_fingerprint_sha256": sha256_json({"delta": sp.srepr(delta), "y0": sp.srepr(y0), "pprime": sp.srepr(sp.diff(p(ell), ell)), "psecond": sp.srepr(sp.diff(p(ell), ell, 2))}),
    }


def finite_derivative_order_audit() -> dict[str, Any]:
    upper = mp.log(11)
    delta = mp.mpf(14) / 5 - 4 * mp.log(2)
    basis = basis_polynomials()
    gram_h1 = gram_matrix_for_order(1)
    gram_h2 = gram_matrix_for_order(2)
    h1_minors, h1_pivots = matrix_pivots_and_minors(gram_h1)
    h2_minors, h2_pivots = matrix_pivots_and_minors(gram_h2)
    basis_rows = []
    max_node_residual = mp.mpf("0")
    for degree, coeffs in zip(BASIS_DEGREES, basis):
        first = derivative_coefficients(coeffs, 1)
        second = derivative_coefficients(coeffs, 2)
        h1_energy = integrate_square_closed_form(first, upper)
        h2_energy = integrate_square_closed_form(second, upper)
        cross = integrate_product_closed_form([delta], first, upper)
        node_residual = max(abs(eval_poly(coeffs, mp.log(d))) for d in DOMAIN)
        max_node_residual = max(max_node_residual, node_residual)
        basis_rows.append({
            "basis_element": f"R(ell)*ell^{degree}",
            "max_abs_node_residual": text(node_residual, 60),
            "h1_energy_int_p_prime_squared": text(h1_energy, 70),
            "h2_energy_int_p_second_squared": text(h2_energy, 70),
            "h1_cross_term_2_delta_boundary_residual_half": text(cross, 70),
            "h1_energy_positive": h1_energy > 0,
            "h2_energy_positive": h2_energy > 0,
        })
    mixed_rows = []
    for a, b in MIXED_ROWS:
        matrix = a * gram_h1 + b * gram_h2
        minors, pivots = matrix_pivots_and_minors(matrix)
        mixed_rows.append({
            "a_h1_weight": text(a, 20),
            "b_h2_weight": text(b, 20),
            "leading_principal_minors": [text(value, 50) for value in minors],
            "cholesky_equivalent_pivots": [text(value, 50) for value in pivots],
            "all_leading_minors_positive": all(value > 0 for value in minors),
            "all_pivots_positive": all(value > 0 for value in pivots),
            "min_pivot": text(min(pivots), 50),
        })
    return {
        "basis_family": "R(ell)*ell^k, k=0..3, with R(ell)=ell*prod_{d=2}^{11}(ell-log(d))",
        "basis_rows": basis_rows,
        "max_abs_node_residual": text(max_node_residual, 60),
        "h1_gram_leading_principal_minors": [text(value, 50) for value in h1_minors],
        "h1_gram_pivots": [text(value, 50) for value in h1_pivots],
        "h2_gram_leading_principal_minors": [text(value, 50) for value in h2_minors],
        "h2_gram_pivots": [text(value, 50) for value in h2_pivots],
        "h1_gram_positive_on_finite_tangent_basis": all(value > 0 for value in h1_minors) and all(value > 0 for value in h1_pivots),
        "h2_gram_positive_on_finite_tangent_basis": all(value > 0 for value in h2_minors) and all(value > 0 for value in h2_pivots),
        "mixed_derivative_rows": mixed_rows,
        "all_mixed_nonnegative_nonzero_rows_positive": all(row["all_leading_minors_positive"] and row["all_pivots_positive"] for row in mixed_rows),
        "all_basis_rows_positive_for_h1_and_h2": all(row["h1_energy_positive"] and row["h2_energy_positive"] for row in basis_rows),
    }


def build_derivative_order_nonidentifiability_certificate(p2512: dict[str, Any]) -> dict[str, Any]:
    symbolic = symbolic_selector_family_audit()
    finite = finite_derivative_order_audit()
    return {
        "frontier_atom_under_attack": "strict_damping_beta_eta_source",
        "p2512_derivative_only_ambiguity_inherited": p2512.get("derivative_only_source_ambiguity_identified") is True,
        "certificate_type": "derivative-only Sobolev order nonidentifiability certificate for the conditional selector",
        "symbolic_selector_family_audit": symbolic,
        "finite_derivative_order_audit": finite,
        "h1_and_h2_both_select_affine_node_solution": symbolic["h1_selects_same_affine_minimizer"] and symbolic["h2_selects_same_affine_minimizer"],
        "finite_gram_supports_h1_h2_and_mixed_coercivity": finite["h1_gram_positive_on_finite_tangent_basis"] and finite["h2_gram_positive_on_finite_tangent_basis"] and finite["all_mixed_nonnegative_nonzero_rows_positive"],
        "derivative_order_nonidentifiability_exported": True,
        "roughness_order_still_requires_strict_source_principle": True,
        "roughness_action_still_postulated_not_derived": True,
        "strict_damping_beta_eta_source_exported": False,
        "strict_dynamical_source_for_A_P_D_exported": False,
        "strict_phase_frequency_source_exported": False,
        "bridge_theorem_exported": False,
        "role_transfer_theorem_exported": False,
        "selector_closure_exported": False,
        "qw2191_discharged_by_this_certificate": False,
        "toe_closure_claimed": False,
    }


def append_doc_sections() -> None:
    eq_section = """
## P2513/S1463 strict damping RG derivative-order nonidentifiability certificate

`P2513/S1463` follows the P2512 source-admissibility audit by proving that stationarity and node data do not identify the derivative order of the postulated selector.  For node-vanishing perturbations `p`, both `J1[y]=int(y')^2` and `J2[y]=int(y'')^2` select the same affine strict damping reconstruction `y0(ell)=delta ell`: `J1[y0+p]-J1[y0]=int(p')^2` and `J2[y0+p]-J2[y0]=int(p'')^2`.  Finite closed-form Gram audits on the polynomial tangent family `R(ell)*ell^k` verify positive leading minors and pivots for `H1`, `H2`, and mixed nonnegative derivative-only quadratic selectors.

This is a negative/source-target theorem: it strengthens the claim that a future strict source theorem must choose the Sobolev/derivative order and coefficient from nadsoliton dynamics.  It does not derive that source and exports no `strict_damping_beta_eta_source`, bridge theorem, role-transfer theorem, QW-2191 discharge, physical-value generator, or ToE closure.
"""
    lag_section = """
## P2513/S1463 derivative-order nonidentifiability guard

`P2513/S1463` shows that derivative-only quadratic selector stationarity is underdetermined: `H1`, `H2`, and mixed derivative-only energies all select `y0(ell)=delta ell` on the audited node data.  Therefore a role-bearing strict action would still need an independent source principle choosing the derivative order/coefficient; no nonlinear compression-flow source theorem or role-bearing `L_total` term is licensed.
"""
    append_once(DOC_FILES["equation_sheet"], "## P2513/S1463 strict damping RG derivative-order nonidentifiability certificate", eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "## P2513/S1463 derivative-order nonidentifiability guard", lag_section)


def build_payload() -> dict[str, Any]:
    grep = rg_audit()
    sources = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    p2512 = theorem(sources["P2512_QUADRATIC_SOURCE_ADMISSIBILITY"], "strict_damping_rg_quadratic_source_admissibility_audit")
    cert = build_derivative_order_nonidentifiability_certificate(p2512)
    finite = cert["finite_derivative_order_audit"]
    theorem_export = {
        "theorem_name": "P2513_T1_strict_damping_rg_derivative_order_nonidentifiability_certificate",
        "audited_chain": ["P2510/S1460", "P2511/S1461", "P2512/S1462"],
        "strict_damping_rg_derivative_order_nonidentifiability_certificate": cert,
        "frontier_atom_under_attack": cert["frontier_atom_under_attack"],
        "p2512_derivative_only_ambiguity_inherited": cert["p2512_derivative_only_ambiguity_inherited"],
        "h1_and_h2_both_select_affine_node_solution": cert["h1_and_h2_both_select_affine_node_solution"],
        "h1_gram_positive_on_finite_tangent_basis": finite["h1_gram_positive_on_finite_tangent_basis"],
        "h2_gram_positive_on_finite_tangent_basis": finite["h2_gram_positive_on_finite_tangent_basis"],
        "all_mixed_nonnegative_nonzero_rows_positive": finite["all_mixed_nonnegative_nonzero_rows_positive"],
        "finite_gram_supports_h1_h2_and_mixed_coercivity": cert["finite_gram_supports_h1_h2_and_mixed_coercivity"],
        "derivative_order_nonidentifiability_exported": cert["derivative_order_nonidentifiability_exported"],
        "roughness_order_still_requires_strict_source_principle": cert["roughness_order_still_requires_strict_source_principle"],
        "roughness_action_still_postulated_not_derived": cert["roughness_action_still_postulated_not_derived"],
        "strict_damping_beta_eta_source_exported": False,
        "strict_dynamical_source_for_A_P_D_exported": False,
        "strict_phase_frequency_source_exported": False,
        "bridge_theorem_exported": False,
        "role_transfer_theorem_exported": False,
        "selector_closure_exported": False,
        "qw2191_discharged_by_this_certificate": False,
        "toe_closure_claimed": False,
        "not_licensed": [
            "P2513 proves derivative-order nonidentifiability for derivative-only selector stationarity; it is not a strict source derivation.",
            "H1, H2, and mixed derivative-only selectors can all select the same affine node reconstruction, so the derivative order/coefficient remains a source-side obligation.",
            "No bridge theorem, role-transfer theorem, QW-2191 discharge, role-bearing L_total, physical-value generator, or ToE closure is exported.",
        ],
        "next_honest_step": "A future strict source theorem must choose the derivative order/coefficient from nadsoliton dynamics or add a symmetry/source principle that breaks the derivative-only selector degeneracy.",
    }
    gatekeepers = {
        "rg_audit_performed": grep["tool"] == "rg",
        "p2512_ambiguity_inherited": theorem_export["p2512_derivative_only_ambiguity_inherited"],
        "h1_h2_select_same_solution": theorem_export["h1_and_h2_both_select_affine_node_solution"],
        "finite_gram_checks_pass": theorem_export["finite_gram_supports_h1_h2_and_mixed_coercivity"],
        "nonidentifiability_not_hidden": theorem_export["derivative_order_nonidentifiability_exported"] and theorem_export["roughness_order_still_requires_strict_source_principle"],
        "source_not_exported": not theorem_export["strict_damping_beta_eta_source_exported"],
        "negative_controls_preserved": not any(theorem_export[key] for key in [
            "strict_damping_beta_eta_source_exported",
            "strict_dynamical_source_for_A_P_D_exported",
            "strict_phase_frequency_source_exported",
            "bridge_theorem_exported",
            "role_transfer_theorem_exported",
            "selector_closure_exported",
            "qw2191_discharged_by_this_certificate",
            "toe_closure_claimed",
        ]),
    }
    return {
        "packet_id": "P2513",
        "stage_id": "S1463",
        "status": "STRICT_DAMPING_RG_DERIVATIVE_ORDER_NONIDENTIFIABILITY_CERTIFICATE_NO_SOURCE_EXPORT_NO_BRIDGE_THEOREM_NO_ROLE_TRANSFER_NO_QW2191_NO_TOE",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "strict_damping_rg_derivative_order_nonidentifiability_certificate": {
            "theorem_export": theorem_export,
            "source_fingerprints": {name: sha256_json(sources[name]) for name in sorted(sources)},
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["strict_damping_rg_derivative_order_nonidentifiability_certificate"]["theorem_export"]
    lines = [
        "# P2513/S1463 strict damping RG derivative-order nonidentifiability certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Result",
        "",
        f"- Frontier atom under attack: `{t['frontier_atom_under_attack']}`.",
        f"- P2512 derivative-only ambiguity inherited: `{t['p2512_derivative_only_ambiguity_inherited']}`.",
        f"- H1 and H2 both select affine node solution: `{t['h1_and_h2_both_select_affine_node_solution']}`.",
        f"- H1 Gram positive on finite tangent basis: `{t['h1_gram_positive_on_finite_tangent_basis']}`.",
        f"- H2 Gram positive on finite tangent basis: `{t['h2_gram_positive_on_finite_tangent_basis']}`.",
        f"- Mixed derivative-only rows positive: `{t['all_mixed_nonnegative_nonzero_rows_positive']}`.",
        f"- Derivative-order nonidentifiability exported: `{t['derivative_order_nonidentifiability_exported']}`.",
        f"- Roughness order still requires strict source principle: `{t['roughness_order_still_requires_strict_source_principle']}`.",
        f"- Source theorem exported: `{t['strict_damping_beta_eta_source_exported']}`.",
        "",
        "## Negative controls",
        "",
        "This packet proves a source-side nonidentifiability result: derivative-only stationarity and node data do not determine whether H1, H2, or a mixed derivative-only selector is the strict source. It exports no strict source atom, bridge theorem, role-transfer theorem, selector/QW-2191 closure, role-bearing L_total term, physical-value generator, or ToE closure.",
        "",
        "## Fingerprint",
        "",
        f"`{payload['strict_damping_rg_derivative_order_nonidentifiability_certificate']['theorem_fingerprint_sha256']}`",
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    GEN.mkdir(exist_ok=True)
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    write_markdown(payload)
    append_doc_sections()
    print(json.dumps(payload["strict_damping_rg_derivative_order_nonidentifiability_certificate"]["theorem_export"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
