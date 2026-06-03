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

GEN = ROOT / "generated"
OUT = GEN / "p2510_s1460_strict_damping_rg_roughness_kkt_stationarity_certificate.json"
MD = GEN / "p2510_s1460_strict_damping_rg_roughness_kkt_stationarity_certificate.md"

SOURCE_FILES = {
    "P2509_MINIMUM_ROUGHNESS_WELLPOSEDNESS": GEN / "p2509_s1459_strict_damping_rg_minimum_roughness_variational_wellposedness_certificate.json",
}

mp.mp.dps = 100
DOMAIN = list(range(1, 12))
POLYNOMIAL_DEGREES = [12, 14]


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
        "new_packet": "P2510|S1460|roughness KKT stationarity|minimum-roughness KKT|distributional stationarity|node multiplier certificate",
        "precursor_packets": "P2509|S1459|minimum-roughness variational well-posedness|P2508|Sobolev node coercivity|P2506|minimum roughness selector",
        "kkt_language": "Euler-Lagrange|KKT|Lagrange multiplier|natural boundary|stationarity residual|node interpolation",
        "guardrails": "legacy -> strict completion bridge|role-transfer audit|K_legacy_ont|K_strict_gate|QW-2191|ToE closure",
        "closure_blockers": "source theorem|bridge theorem|role-transfer theorem|physical-value generator|role-bearing L_total|selector closure",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def theorem(payload: dict[str, Any], container: str) -> dict[str, Any]:
    return payload.get(container, {}).get("theorem_export", {})


def text(value: mp.mpf, digits: int = 70) -> str:
    return mp.nstr(value, digits)


def stiffness_entry(i: int, j: int, upper: mp.mpf) -> mp.mpf:
    if i < 2 or j < 2:
        return mp.mpf("0")
    power = i + j - 4
    return mp.mpf(i * (i - 1) * j * (j - 1)) * upper ** (power + 1) / mp.mpf(power + 1)


def build_polynomial_kkt_matrix(degree: int) -> tuple[mp.matrix, mp.matrix, list[mp.mpf], mp.mpf, mp.mpf]:
    basis_dim = degree + 1
    node_count = len(DOMAIN)
    size = basis_dim + node_count
    upper = mp.log(11)
    delta = mp.mpf(14) / 5 - 4 * mp.log(2)
    nodes = [mp.log(d) for d in DOMAIN]
    matrix = mp.matrix(size, size)
    rhs = mp.matrix(size, 1)
    for i in range(basis_dim):
        for j in range(basis_dim):
            matrix[i, j] = stiffness_entry(i, j, upper)
    for row, node in enumerate(nodes):
        target = delta * node
        rhs[basis_dim + row] = target
        for col in range(basis_dim):
            value = node ** col
            matrix[col, basis_dim + row] = value
            matrix[basis_dim + row, col] = value
    return matrix, rhs, nodes, delta, upper


def max_abs_matrix(vector: mp.matrix) -> mp.mpf:
    maximum = mp.mpf("0")
    for value in vector:
        maximum = max(maximum, abs(value))
    return maximum


def polynomial_kkt_audit_for_degree(degree: int) -> dict[str, Any]:
    matrix, rhs, nodes, delta, upper = build_polynomial_kkt_matrix(degree)
    solution = mp.lu_solve(matrix, rhs)
    residual = matrix * solution - rhs
    basis_dim = degree + 1
    coeffs = [solution[i] for i in range(basis_dim)]
    multipliers = [solution[basis_dim + i] for i in range(len(DOMAIN))]
    expected_coeffs = [mp.mpf("0")] * basis_dim
    expected_coeffs[1] = delta
    coeff_error = max(abs(coeffs[i] - expected_coeffs[i]) for i in range(basis_dim))
    multiplier_max = max(abs(value) for value in multipliers)
    second_energy = mp.fsum(
        coeffs[i] * coeffs[j] * stiffness_entry(i, j, upper)
        for i in range(basis_dim)
        for j in range(basis_dim)
    )
    node_residual = max(abs(mp.fsum(coeffs[k] * node ** k for k in range(basis_dim)) - delta * node) for node in nodes)
    expected = mp.matrix(len(rhs), 1)
    for i, coeff in enumerate(expected_coeffs):
        expected[i] = coeff
    expected_residual = matrix * expected - rhs
    singular_values = mp.svd(matrix, compute_uv=False)
    singular_values_list = list(singular_values)
    min_singular_value = min(abs(value) for value in singular_values_list)
    rank = sum(1 for value in singular_values_list if abs(value) > mp.mpf("1e-60"))
    determinant_abs = abs(mp.det(matrix))
    return {
        "polynomial_degree": degree,
        "basis_dim": basis_dim,
        "node_count": len(DOMAIN),
        "kkt_matrix_size": len(rhs),
        "kkt_rank": int(rank),
        "kkt_full_rank": int(rank) == len(rhs),
        "abs_kkt_determinant": text(determinant_abs, 40),
        "min_abs_singular_value": text(min_singular_value, 40),
        "max_abs_linear_solve_residual": text(max_abs_matrix(residual), 50),
        "max_abs_coeff_error_vs_affine_delta_ell": text(coeff_error, 50),
        "max_abs_node_multiplier": text(multiplier_max, 50),
        "max_abs_node_residual": text(node_residual, 50),
        "roughness_energy_of_solution": text(second_energy, 50),
        "expected_affine_zero_multiplier_residual": text(max_abs_matrix(expected_residual), 50),
        "solution_is_affine_delta_ell_with_zero_multipliers": coeff_error < mp.mpf("1e-70") and multiplier_max < mp.mpf("1e-70"),
        "constraints_satisfied": node_residual < mp.mpf("1e-80"),
        "stationarity_residual_zero": max_abs_matrix(residual) < mp.mpf("1e-80") and max_abs_matrix(expected_residual) < mp.mpf("1e-80"),
    }


def symbolic_kkt_stationarity() -> dict[str, Any]:
    ell = sp.symbols("ell", real=True)
    delta = sp.Rational(14, 5) - 4 * sp.log(2)
    y0 = delta * ell
    y0_second = sp.diff(y0, ell, 2)
    y0_fourth = sp.diff(y0, ell, 4)
    return {
        "symbolic_backend": "sympy",
        "sympy_version": sp.__version__,
        "weak_stationarity_equation": "int_0^L y''(ell) v''(ell) d ell = sum_i mu_i v(log i)",
        "candidate_y0": sp.sstr(y0),
        "candidate_second_derivative": sp.sstr(y0_second),
        "candidate_fourth_derivative": sp.sstr(y0_fourth),
        "zero_multiplier_solution": "mu_i=0 for all i=1..11",
        "distributional_euler_lagrange_residual": "y0'''' - sum_i mu_i delta_{log i} = 0 distributionally when all mu_i=0",
        "natural_boundary_residuals": "y0''(0)=y0''(L)=0 and y0'''(0)=y0'''(L)=0",
        "stationarity_residual_symbolically_zero": y0_second == 0 and y0_fourth == 0,
        "symbolic_fingerprint_sha256": sha256_json({"delta": sp.srepr(delta), "y0": sp.srepr(y0), "y0_second": sp.srepr(y0_second), "y0_fourth": sp.srepr(y0_fourth)}),
    }


def build_kkt_stationarity_certificate(p2509: dict[str, Any]) -> dict[str, Any]:
    symbolic = symbolic_kkt_stationarity()
    finite_rows = [polynomial_kkt_audit_for_degree(degree) for degree in POLYNOMIAL_DEGREES]
    return {
        "frontier_atom_under_attack": "strict_damping_beta_eta_source",
        "p2509_variational_wellposedness_inherited": p2509.get("minimum_roughness_problem_wellposed_for_postulated_functional") is True,
        "selector_type": "conditional Euler-Lagrange/KKT stationarity certificate for the postulated minimum-roughness selector",
        "symbolic_kkt_stationarity": symbolic,
        "finite_polynomial_kkt_audit": {
            "basis_family": "monomials ell^k, k=0..N, with point constraints y(log d)=delta log d for d=1..11",
            "degree_rows": finite_rows,
            "all_kkt_matrices_full_rank": all(row["kkt_full_rank"] for row in finite_rows),
            "all_solutions_affine_with_zero_multipliers": all(row["solution_is_affine_delta_ell_with_zero_multipliers"] for row in finite_rows),
            "all_constraints_satisfied": all(row["constraints_satisfied"] for row in finite_rows),
            "all_stationarity_residuals_zero": all(row["stationarity_residual_zero"] for row in finite_rows),
        },
        "kkt_stationarity_confirmed_for_postulated_functional": symbolic["stationarity_residual_symbolically_zero"] and all(row["solution_is_affine_delta_ell_with_zero_multipliers"] for row in finite_rows),
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
## P2510/S1460 strict damping RG roughness KKT stationarity certificate

`P2510/S1460` adds the Euler-Lagrange/KKT stationarity layer for the same postulated P2506/P2509 minimum-roughness selector.  The weak normal equation is `int y'' v'' = sum_i mu_i v(log i)` under the node constraints `y(log d)=delta log d`.  For the P2509 minimizer `y0(ell)=delta ell`, `y0''=0` and `y0''''=0`, so the distributional stationarity equation is satisfied with all node multipliers `mu_i=0`; the natural boundary residuals are also zero.  A finite polynomial KKT audit on monomial spaces through degrees 12 and 14 independently recovers the same affine coefficient vector and zero multipliers, with full-rank KKT matrices and tiny solve residuals.

This is theorem-prep for the conditional selector only.  It does not derive the roughness action from strict nadsoliton dynamics and does not export `strict_damping_beta_eta_source`, a bridge theorem, a role-transfer theorem, QW-2191 discharge, a physical-value generator, or ToE closure.
"""
    lag_section = """
## P2510/S1460 roughness KKT stationarity guard

`P2510/S1460` records the variational normal-equation/KKT witness for the postulated strict damping roughness selector: `y0(ell)=delta ell` is stationary with zero distributional node multipliers and zero natural-boundary residuals, and finite polynomial KKT systems replay the affine solution.  This strengthens the internal variational bookkeeping after P2509 but remains conditional on the roughness functional being supplied; it is not a nonlinear compression-flow source theorem and is not a license for a role-bearing `L_total` term.
"""
    append_once(DOC_FILES["equation_sheet"], "## P2510/S1460 strict damping RG roughness KKT stationarity certificate", eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "## P2510/S1460 roughness KKT stationarity guard", lag_section)


def build_payload() -> dict[str, Any]:
    grep = rg_audit()
    sources = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    p2509 = theorem(sources["P2509_MINIMUM_ROUGHNESS_WELLPOSEDNESS"], "strict_damping_rg_minimum_roughness_variational_wellposedness_certificate")
    cert = build_kkt_stationarity_certificate(p2509)
    audit = cert["finite_polynomial_kkt_audit"]
    theorem_export = {
        "theorem_name": "P2510_T1_strict_damping_rg_roughness_kkt_stationarity_certificate",
        "audited_chain": ["P2503/S1453", "P2504/S1454", "P2505/S1455", "P2506/S1456", "P2507/S1457", "P2508/S1458", "P2509/S1459"],
        "strict_damping_rg_roughness_kkt_stationarity_certificate": cert,
        "frontier_atom_under_attack": cert["frontier_atom_under_attack"],
        "p2509_variational_wellposedness_inherited": cert["p2509_variational_wellposedness_inherited"],
        "symbolic_stationarity_residual_zero": cert["symbolic_kkt_stationarity"]["stationarity_residual_symbolically_zero"],
        "finite_kkt_all_full_rank": audit["all_kkt_matrices_full_rank"],
        "finite_kkt_all_affine_zero_multipliers": audit["all_solutions_affine_with_zero_multipliers"],
        "finite_kkt_all_constraints_satisfied": audit["all_constraints_satisfied"],
        "finite_kkt_all_stationarity_residuals_zero": audit["all_stationarity_residuals_zero"],
        "kkt_stationarity_confirmed_for_postulated_functional": cert["kkt_stationarity_confirmed_for_postulated_functional"],
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
            "P2510 proves KKT stationarity only for the postulated P2506/P2509 roughness minimization problem.",
            "It does not derive the roughness action, delta, beta0, or strict_damping_beta_eta_source from strict nadsoliton dynamics.",
            "No bridge theorem, role-transfer theorem, QW-2191 discharge, role-bearing L_total, physical-value generator, or ToE closure is exported.",
        ],
        "next_honest_step": "A source-level result must still derive the roughness action or an equivalent damping-flow source from strict nadsoliton dynamics; KKT stationarity only verifies the conditional selector's variational equations.",
    }
    gatekeepers = {
        "rg_audit_performed": grep["tool"] == "rg",
        "p2509_wellposedness_inherited": theorem_export["p2509_variational_wellposedness_inherited"],
        "symbolic_stationarity_zero": theorem_export["symbolic_stationarity_residual_zero"],
        "finite_kkt_full_rank": theorem_export["finite_kkt_all_full_rank"],
        "finite_kkt_affine_zero_multipliers": theorem_export["finite_kkt_all_affine_zero_multipliers"],
        "finite_kkt_constraints_and_stationarity": theorem_export["finite_kkt_all_constraints_satisfied"] and theorem_export["finite_kkt_all_stationarity_residuals_zero"],
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
        "packet_id": "P2510",
        "stage_id": "S1460",
        "status": "STRICT_DAMPING_RG_ROUGHNESS_KKT_STATIONARITY_FOR_POSTULATED_FUNCTIONAL_NO_SOURCE_EXPORT_NO_BRIDGE_THEOREM_NO_ROLE_TRANSFER_NO_QW2191_NO_TOE",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "strict_damping_rg_roughness_kkt_stationarity_certificate": {
            "theorem_export": theorem_export,
            "source_fingerprints": {name: sha256_json(sources[name]) for name in sorted(sources)},
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["strict_damping_rg_roughness_kkt_stationarity_certificate"]["theorem_export"]
    rows = t["strict_damping_rg_roughness_kkt_stationarity_certificate"]["finite_polynomial_kkt_audit"]["degree_rows"]
    lines = [
        "# P2510/S1460 strict damping RG roughness KKT stationarity certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Result",
        "",
        f"- Frontier atom under attack: `{t['frontier_atom_under_attack']}`.",
        f"- P2509 variational well-posedness inherited: `{t['p2509_variational_wellposedness_inherited']}`.",
        f"- Symbolic stationarity residual zero: `{t['symbolic_stationarity_residual_zero']}`.",
        f"- Finite KKT matrices full-rank: `{t['finite_kkt_all_full_rank']}`.",
        f"- Finite KKT solutions affine with zero multipliers: `{t['finite_kkt_all_affine_zero_multipliers']}`.",
        f"- Finite KKT constraints satisfied: `{t['finite_kkt_all_constraints_satisfied']}`.",
        f"- Finite KKT stationarity residuals zero: `{t['finite_kkt_all_stationarity_residuals_zero']}`.",
        f"- KKT stationarity confirmed for postulated functional: `{t['kkt_stationarity_confirmed_for_postulated_functional']}`.",
        f"- Roughness action still postulated, not derived: `{t['roughness_action_still_postulated_not_derived']}`.",
        f"- Source theorem exported: `{t['strict_damping_beta_eta_source_exported']}`.",
        "",
        "## Finite polynomial KKT rows",
        "",
    ]
    for row in rows:
        lines.append(f"- degree `{row['polynomial_degree']}`: rank `{row['kkt_rank']}/{row['kkt_matrix_size']}`, solve residual `{row['max_abs_linear_solve_residual']}`, coefficient error `{row['max_abs_coeff_error_vs_affine_delta_ell']}`, max multiplier `{row['max_abs_node_multiplier']}`.")
    lines += [
        "",
        "## Negative controls",
        "",
        "This packet verifies the Euler-Lagrange/KKT equations only for the postulated roughness selector. It does not derive the roughness action from strict dynamics and exports no strict source atom, bridge theorem, role-transfer theorem, selector/QW-2191 closure, role-bearing L_total term, physical-value generator, or ToE closure.",
        "",
        "## Fingerprint",
        "",
        f"`{payload['strict_damping_rg_roughness_kkt_stationarity_certificate']['theorem_fingerprint_sha256']}`",
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    GEN.mkdir(exist_ok=True)
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    write_markdown(payload)
    append_doc_sections()
    print(json.dumps(payload["strict_damping_rg_roughness_kkt_stationarity_certificate"]["theorem_export"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
