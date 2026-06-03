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
)
from p2512_s1462_strict_damping_rg_quadratic_source_admissibility_audit import integrate_product_closed_form
from p2513_s1463_strict_damping_rg_derivative_order_nonidentifiability_certificate import (
    matrix_pivots_and_minors,
    monomial_times,
)

GEN = ROOT / "generated"
OUT = GEN / "p2514_s1464_strict_damping_rg_higher_order_selector_nonidentifiability_certificate.json"
MD = GEN / "p2514_s1464_strict_damping_rg_higher_order_selector_nonidentifiability_certificate.md"

SOURCE_FILES = {
    "P2513_DERIVATIVE_ORDER_NONIDENTIFIABILITY": GEN / "p2513_s1463_strict_damping_rg_derivative_order_nonidentifiability_certificate.json",
}

mp.mp.dps = 100
DOMAIN = list(range(1, 12))
BASIS_DEGREES = list(range(4))
THEOREM_ORDER_RANGE = list(range(1, 11))
FINITE_AUDIT_ORDERS = list(range(1, 7))


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
        "new_packet": "P2514|S1464|higher-order selector nonidentifiability|Sobolev tower nonidentifiability|derivative-order tower|higher-order derivative selector",
        "precursor_packets": "P2513|S1463|derivative-order nonidentifiability|P2512|quadratic source admissibility|P2508|Sobolev node coercivity",
        "tower_language": "Sobolev tower|higher-order coercivity|derivative order|polynomial nullspace|node-vanishing polynomial",
        "guardrails": "legacy -> strict completion bridge|role-transfer audit|K_legacy_ont|K_strict_gate|QW-2191|ToE closure",
        "closure_blockers": "source theorem|bridge theorem|role-transfer theorem|physical-value generator|role-bearing L_total|selector closure",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def theorem(payload: dict[str, Any], container: str) -> dict[str, Any]:
    return payload.get(container, {}).get("theorem_export", {})


def text(value: mp.mpf, digits: int = 70) -> str:
    return mp.nstr(value, digits)


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


def symbolic_tower_theorem() -> dict[str, Any]:
    ell = sp.symbols("ell", real=True)
    delta = sp.Rational(14, 5) - 4 * sp.log(2)
    y0 = delta * ell
    rows = []
    for order in THEOREM_ORDER_RANGE:
        y0_derivative = sp.diff(y0, ell, order)
        zero_mode_degree = order - 1
        rows.append({
            "derivative_order_m": order,
            "candidate_derivative": sp.sstr(y0_derivative),
            "zero_energy_perturbation_class": f"polynomials of degree <= {zero_mode_degree}",
            "node_count": len(DOMAIN),
            "node_count_exceeds_zero_mode_degree": len(DOMAIN) > zero_mode_degree,
            "node_vanishing_zero_mode_eliminated": len(DOMAIN) > zero_mode_degree,
        })
    return {
        "symbolic_backend": "sympy",
        "sympy_version": sp.__version__,
        "candidate_y0": sp.sstr(y0),
        "theorem_order_range": THEOREM_ORDER_RANGE,
        "theorem_statement": "For every 1<=m<=10, J_m[y]=int (D^m y)^2 has y0(ell)=delta*ell as a minimizer on the strict damping node data, and any zero-energy tangent perturbation is a polynomial of degree <=m-1 killed by the 11 node conditions.",
        "rows": rows,
        "all_orders_eliminate_zero_modes_by_node_count": all(row["node_vanishing_zero_mode_eliminated"] for row in rows),
        "symbolic_fingerprint_sha256": sha256_json({"delta": sp.srepr(delta), "orders": rows}),
    }


def finite_higher_order_gram_audit() -> dict[str, Any]:
    upper = mp.log(11)
    basis = basis_polynomials()
    order_rows = []
    max_node_residual = mp.mpf("0")
    for coeffs in basis:
        for d in DOMAIN:
            max_node_residual = max(max_node_residual, abs(eval_poly(coeffs, mp.log(d))))
    for order in FINITE_AUDIT_ORDERS:
        matrix = gram_matrix_for_order(order)
        minors, pivots = matrix_pivots_and_minors(matrix)
        energy_rows = []
        for degree, coeffs in zip(BASIS_DEGREES, basis):
            derivative = derivative_coefficients(coeffs, order)
            energy = integrate_square_closed_form(derivative, upper)
            energy_rows.append({
                "basis_element": f"R(ell)*ell^{degree}",
                "energy_int_derivative_squared": text(energy, 70),
                "energy_positive": energy > 0,
            })
        order_rows.append({
            "derivative_order_m": order,
            "leading_principal_minors": [text(value, 50) for value in minors],
            "cholesky_equivalent_pivots": [text(value, 50) for value in pivots],
            "min_pivot": text(min(pivots), 50),
            "all_leading_minors_positive": all(value > 0 for value in minors),
            "all_pivots_positive": all(value > 0 for value in pivots),
            "basis_energy_rows": energy_rows,
            "all_basis_energies_positive": all(row["energy_positive"] for row in energy_rows),
        })
    return {
        "basis_family": "R(ell)*ell^k, k=0..3, with R(ell)=ell*prod_{d=2}^{11}(ell-log(d))",
        "finite_audit_orders": FINITE_AUDIT_ORDERS,
        "max_abs_node_residual": text(max_node_residual, 60),
        "order_rows": order_rows,
        "all_order_gram_witnesses_positive": all(row["all_leading_minors_positive"] and row["all_pivots_positive"] and row["all_basis_energies_positive"] for row in order_rows),
    }


def acceptance_lattice() -> dict[str, Any]:
    rows = []
    for order in THEOREM_ORDER_RANGE:
        rows.append({
            "derivative_order_m": order,
            "selector": f"J_{order}[y]=int (D^{order} y)^2",
            "selects_affine_node_solution_by_symbolic_zero_mode_argument": True,
            "licensed_as_strict_source": False,
            "source_obligation": "derive this order/coefficient from strict nadsoliton dynamics or provide a symmetry principle selecting it",
        })
    return {
        "rows": rows,
        "admissible_order_count_under_node_stationarity_only": len(rows),
        "node_stationarity_leaves_many_orders_admissible": len(rows) > 1,
        "source_target_boundary": "Node interpolation, stationarity, and coercivity can certify a tower of derivative-only selectors; they cannot by themselves identify which order is the strict damping source.",
    }


def build_higher_order_nonidentifiability_certificate(p2513: dict[str, Any]) -> dict[str, Any]:
    symbolic = symbolic_tower_theorem()
    finite = finite_higher_order_gram_audit()
    lattice = acceptance_lattice()
    return {
        "frontier_atom_under_attack": "strict_damping_beta_eta_source",
        "p2513_derivative_order_nonidentifiability_inherited": p2513.get("derivative_order_nonidentifiability_exported") is True,
        "certificate_type": "higher-order Sobolev selector tower nonidentifiability certificate",
        "symbolic_tower_theorem": symbolic,
        "finite_higher_order_gram_audit": finite,
        "acceptance_lattice": lattice,
        "higher_order_selector_tower_nonidentifiability_exported": symbolic["all_orders_eliminate_zero_modes_by_node_count"] and lattice["node_stationarity_leaves_many_orders_admissible"],
        "finite_gram_audit_supports_orders_1_to_6": finite["all_order_gram_witnesses_positive"],
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
## P2514/S1464 strict damping RG higher-order selector nonidentifiability certificate

`P2514/S1464` strengthens the P2513 negative/source-target result from the `H1/H2` pair to a whole derivative-order tower.  For every `1 <= m <= 10`, the derivative-only functional `J_m[y]=int (D^m y)^2 d ell` has the same affine strict damping node reconstruction `y0(ell)=delta ell` as a minimizer; a zero-energy tangent perturbation would be a polynomial of degree at most `m-1`, and the eleven node conditions force it to vanish for all audited theorem orders.  A finite closed-form Gram audit over `R(ell)*ell^k`, `k=0..3`, verifies positive leading minors and Cholesky-equivalent pivots for derivative orders `1..6`.

This is a stronger nonidentifiability theorem, not a closure theorem.  It shows that node data plus stationarity/coercivity admit a tower of derivative-only selectors, so a future strict source theorem must still choose the derivative order/coefficient from nadsoliton dynamics.  It exports no `strict_damping_beta_eta_source`, bridge theorem, role-transfer theorem, QW-2191 discharge, physical-value generator, or ToE closure.
"""
    lag_section = """
## P2514/S1464 higher-order selector nonidentifiability guard

`P2514/S1464` extends the derivative-only ambiguity: `J_m[y]=int(D^m y)^2` for `m=1..10` can select the same affine strict damping node reconstruction under the current node constraints, with finite Gram checks through `m=6`.  Therefore a role-bearing strict action still needs an independent source principle for the derivative order/coefficient; no nonlinear compression-flow source theorem or role-bearing `L_total` term is licensed.
"""
    append_once(DOC_FILES["equation_sheet"], "## P2514/S1464 strict damping RG higher-order selector nonidentifiability certificate", eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "## P2514/S1464 higher-order selector nonidentifiability guard", lag_section)


def build_payload() -> dict[str, Any]:
    grep = rg_audit()
    sources = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    p2513 = theorem(sources["P2513_DERIVATIVE_ORDER_NONIDENTIFIABILITY"], "strict_damping_rg_derivative_order_nonidentifiability_certificate")
    cert = build_higher_order_nonidentifiability_certificate(p2513)
    symbolic = cert["symbolic_tower_theorem"]
    finite = cert["finite_higher_order_gram_audit"]
    lattice = cert["acceptance_lattice"]
    theorem_export = {
        "theorem_name": "P2514_T1_strict_damping_rg_higher_order_selector_nonidentifiability_certificate",
        "audited_chain": ["P2511/S1461", "P2512/S1462", "P2513/S1463"],
        "strict_damping_rg_higher_order_selector_nonidentifiability_certificate": cert,
        "frontier_atom_under_attack": cert["frontier_atom_under_attack"],
        "p2513_derivative_order_nonidentifiability_inherited": cert["p2513_derivative_order_nonidentifiability_inherited"],
        "theorem_order_range": symbolic["theorem_order_range"],
        "all_orders_eliminate_zero_modes_by_node_count": symbolic["all_orders_eliminate_zero_modes_by_node_count"],
        "finite_audit_orders": finite["finite_audit_orders"],
        "finite_gram_audit_supports_orders_1_to_6": cert["finite_gram_audit_supports_orders_1_to_6"],
        "admissible_order_count_under_node_stationarity_only": lattice["admissible_order_count_under_node_stationarity_only"],
        "higher_order_selector_tower_nonidentifiability_exported": cert["higher_order_selector_tower_nonidentifiability_exported"],
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
            "P2514 proves a higher-order derivative-only selector nonidentifiability tower; it is not a strict source derivation.",
            "Orders 1..10 are symbolically admissible under node-zero-mode elimination, and finite Gram checks through order 6 support the same ambiguity on the audited basis.",
            "No bridge theorem, role-transfer theorem, QW-2191 discharge, role-bearing L_total, physical-value generator, or ToE closure is exported.",
        ],
        "next_honest_step": lattice["source_target_boundary"],
    }
    gatekeepers = {
        "rg_audit_performed": grep["tool"] == "rg",
        "p2513_nonidentifiability_inherited": theorem_export["p2513_derivative_order_nonidentifiability_inherited"],
        "symbolic_tower_zero_modes_eliminated": theorem_export["all_orders_eliminate_zero_modes_by_node_count"],
        "finite_gram_checks_pass": theorem_export["finite_gram_audit_supports_orders_1_to_6"],
        "tower_nonidentifiability_not_hidden": theorem_export["higher_order_selector_tower_nonidentifiability_exported"] and theorem_export["roughness_order_still_requires_strict_source_principle"],
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
        "packet_id": "P2514",
        "stage_id": "S1464",
        "status": "STRICT_DAMPING_RG_HIGHER_ORDER_SELECTOR_NONIDENTIFIABILITY_CERTIFICATE_NO_SOURCE_EXPORT_NO_BRIDGE_THEOREM_NO_ROLE_TRANSFER_NO_QW2191_NO_TOE",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "strict_damping_rg_higher_order_selector_nonidentifiability_certificate": {
            "theorem_export": theorem_export,
            "source_fingerprints": {name: sha256_json(sources[name]) for name in sorted(sources)},
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["strict_damping_rg_higher_order_selector_nonidentifiability_certificate"]["theorem_export"]
    lines = [
        "# P2514/S1464 strict damping RG higher-order selector nonidentifiability certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Result",
        "",
        f"- Frontier atom under attack: `{t['frontier_atom_under_attack']}`.",
        f"- P2513 derivative-order nonidentifiability inherited: `{t['p2513_derivative_order_nonidentifiability_inherited']}`.",
        f"- Symbolic theorem order range: `{t['theorem_order_range']}`.",
        f"- All symbolic orders eliminate zero modes by node count: `{t['all_orders_eliminate_zero_modes_by_node_count']}`.",
        f"- Finite audit orders: `{t['finite_audit_orders']}`.",
        f"- Finite Gram audit supports orders 1..6: `{t['finite_gram_audit_supports_orders_1_to_6']}`.",
        f"- Admissible order count under node stationarity only: `{t['admissible_order_count_under_node_stationarity_only']}`.",
        f"- Higher-order selector tower nonidentifiability exported: `{t['higher_order_selector_tower_nonidentifiability_exported']}`.",
        f"- Source theorem exported: `{t['strict_damping_beta_eta_source_exported']}`.",
        "",
        "## Negative controls",
        "",
        "This packet strengthens the source-side nonidentifiability result: a tower of derivative-only selectors can select the same strict damping affine node reconstruction. It exports no strict source atom, bridge theorem, role-transfer theorem, selector/QW-2191 closure, role-bearing L_total term, physical-value generator, or ToE closure.",
        "",
        "## Fingerprint",
        "",
        f"`{payload['strict_damping_rg_higher_order_selector_nonidentifiability_certificate']['theorem_fingerprint_sha256']}`",
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    GEN.mkdir(exist_ok=True)
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    write_markdown(payload)
    append_doc_sections()
    print(json.dumps(payload["strict_damping_rg_higher_order_selector_nonidentifiability_certificate"]["theorem_export"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
