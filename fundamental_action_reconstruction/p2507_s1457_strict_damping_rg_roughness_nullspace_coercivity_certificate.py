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
OUT = GEN / "p2507_s1457_strict_damping_rg_roughness_nullspace_coercivity_certificate.json"
MD = GEN / "p2507_s1457_strict_damping_rg_roughness_nullspace_coercivity_certificate.md"

SOURCE_FILES = {
    "P2506_ROUGHNESS_SELECTOR": GEN / "p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate.json",
}

mp.mp.dps = 90
DOMAIN = list(range(1, 12))
BASIS_DEGREE = 3


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
        "new_packet": "P2507|S1457|roughness nullspace coercivity|RG roughness coercivity|polynomial nullspace coercivity|roughness Gram matrix",
        "precursor_packets": "P2506|S1456|minimum roughness selector|P2505|finite-node RG-flow nullspace|P2504|P2503",
        "coercivity_language": "coercivity|positive definite|Gram matrix|Sobolev roughness|polynomial nullspace|node-vanishing perturbation",
        "guardrails": "legacy -> strict completion bridge|role-transfer audit|K_legacy_ont|K_strict_gate|QW-2191|ToE closure",
        "closure_blockers": "source theorem|bridge theorem|role-transfer theorem|physical-value generator|role-bearing L_total|selector closure",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def theorem(payload: dict[str, Any], container: str) -> dict[str, Any]:
    return payload.get(container, {}).get("theorem_export", {})


def text(value: mp.mpf, digits: int = 60) -> str:
    return mp.nstr(value, digits)


def poly_mul(a: list[mp.mpf], b: list[mp.mpf]) -> list[mp.mpf]:
    out = [mp.mpf("0")] * (len(a) + len(b) - 1)
    for i, av in enumerate(a):
        for j, bv in enumerate(b):
            out[i + j] += av * bv
    return out


def derivative_coefficients(coeffs: list[mp.mpf], order: int) -> list[mp.mpf]:
    current = coeffs[:]
    for _ in range(order):
        current = [mp.mpf(index + 1) * current[index + 1] for index in range(len(current) - 1)]
    return current or [mp.mpf("0")]


def eval_poly(coeffs: list[mp.mpf], x: mp.mpf) -> mp.mpf:
    value = mp.mpf("0")
    for coeff in reversed(coeffs):
        value = value * x + coeff
    return value


def node_vanishing_r_coefficients() -> list[mp.mpf]:
    # Ascending coefficients for R(ell)=ell*prod_{d=2}^{11}(ell-log(d)).
    coeffs = [mp.mpf("1")]
    for root in [mp.mpf("0")] + [mp.log(d) for d in range(2, 12)]:
        coeffs = poly_mul(coeffs, [-root, mp.mpf("1")])
    return coeffs


def convolve(a: list[mp.mpf], b: list[mp.mpf]) -> list[mp.mpf]:
    return poly_mul(a, b)


def integrate_product_closed_form(coeffs_a: list[mp.mpf], coeffs_b: list[mp.mpf], upper: mp.mpf) -> mp.mpf:
    # Exact antiderivative of the polynomial product, evaluated with high-precision
    # logarithmic coefficients.  This avoids using quadrature as the primary Gram
    # source and makes the finite-dimensional audit reproducible coefficient by
    # coefficient.
    product_coeffs = convolve(coeffs_a, coeffs_b)
    return mp.fsum(coeff * upper ** (power + 1) / mp.mpf(power + 1) for power, coeff in enumerate(product_coeffs))


def integrate_product_split_quadrature(coeffs_a: list[mp.mpf], coeffs_b: list[mp.mpf], split_nodes: list[mp.mpf]) -> mp.mpf:
    def product(x: mp.mpf) -> mp.mpf:
        return eval_poly(coeffs_a, x) * eval_poly(coeffs_b, x)

    return mp.fsum(mp.quad(product, [split_nodes[i], split_nodes[i + 1]]) for i in range(len(split_nodes) - 1))


def cholesky_pivots(matrix: list[list[mp.mpf]]) -> list[mp.mpf]:
    n = len(matrix)
    lower = [[mp.mpf("0") for _ in range(n)] for _ in range(n)]
    pivots: list[mp.mpf] = []
    for i in range(n):
        for j in range(i + 1):
            accum = mp.fsum(lower[i][k] * lower[j][k] for k in range(j))
            if i == j:
                pivot = matrix[i][i] - accum
                pivots.append(pivot)
                lower[i][j] = mp.sqrt(pivot) if pivot > 0 else mp.nan
            else:
                lower[i][j] = (matrix[i][j] - accum) / lower[j][j]
    return pivots


def symbolic_coercivity_argument() -> dict[str, Any]:
    ell = sp.symbols("ell", real=True)
    delta = sp.Rational(14, 5) - 4 * sp.log(2)
    return {
        "symbolic_backend": "sympy",
        "sympy_version": sp.__version__,
        "selector_functional": "J[y]=integral_0^log(11) (y''(ell))^2 d ell",
        "affine_kernel_of_roughness": "J[p]=0 for an admitted polynomial perturbation p only when p''=0, hence p=a*ell+b.",
        "node_vanishing_affine_elimination": "Node constraints p(0)=0 and p(log(2))=0 force b=0 and a=0, so no nonzero node-vanishing affine perturbation survives.",
        "consequence": "The roughness quadratic form is positive definite on any finite-dimensional polynomial subspace of nonzero perturbations that vanish at all audited nodes.",
        "delta_reference": sp.sstr(delta),
        "symbolic_fingerprint_sha256": sha256_json({"delta": sp.srepr(delta), "basis_degree": BASIS_DEGREE}),
    }


def gram_coercivity_certificate() -> dict[str, Any]:
    r_coeffs = node_vanishing_r_coefficients()
    split_nodes = [mp.mpf("0")] + [mp.log(d) for d in range(2, 12)]
    basis_coeffs = [poly_mul(r_coeffs, [mp.mpf("0")] * k + [mp.mpf("1")]) for k in range(BASIS_DEGREE + 1)]
    basis_second = [derivative_coefficients(coeffs, 2) for coeffs in basis_coeffs]
    upper = mp.log(11)
    gram = [[integrate_product_closed_form(basis_second[i], basis_second[j], upper) for j in range(BASIS_DEGREE + 1)] for i in range(BASIS_DEGREE + 1)]
    quadrature_gram = [[integrate_product_split_quadrature(basis_second[i], basis_second[j], split_nodes) for j in range(BASIS_DEGREE + 1)] for i in range(BASIS_DEGREE + 1)]
    closed_form_quadrature_residual = max(
        abs(gram[i][j] - quadrature_gram[i][j])
        for i in range(BASIS_DEGREE + 1)
        for j in range(BASIS_DEGREE + 1)
    )
    leading_principal_minors = [mp.det(mp.matrix([row[:size] for row in gram[:size]])) for size in range(1, BASIS_DEGREE + 2)]
    pivots = cholesky_pivots(gram)
    diagonal = [gram[i][i] for i in range(BASIS_DEGREE + 1)]
    symmetry_residual = max(abs(gram[i][j] - gram[j][i]) for i in range(BASIS_DEGREE + 1) for j in range(BASIS_DEGREE + 1))

    witness_coefficients = [mp.mpf("1e-6"), mp.mpf("-2e-6"), mp.mpf("3e-6"), mp.mpf("-1e-6")]
    witness_energy = mp.fsum(
        witness_coefficients[i] * gram[i][j] * witness_coefficients[j]
        for i in range(BASIS_DEGREE + 1)
        for j in range(BASIS_DEGREE + 1)
    )
    zero_node_residuals = []
    for d in DOMAIN:
        ell_d = mp.log(d)
        value = mp.fsum(witness_coefficients[k] * eval_poly(basis_coeffs[k], ell_d) for k in range(BASIS_DEGREE + 1))
        zero_node_residuals.append(abs(value))

    return {
        "basis_family": "phi_k(ell)=R(ell)*ell^k for k=0..3, R(ell)=ell*prod_{d=2}^{11}(ell-log(d))",
        "basis_dimension": BASIS_DEGREE + 1,
        "integration_interval": ["0", text(mp.log(11), 50)],
        "split_count": len(split_nodes) - 1,
        "primary_integration_method": "closed-form polynomial antiderivative evaluated with 90 decimal digits; split quadrature is retained only as a cross-check",
        "closed_form_vs_split_quadrature_max_abs_residual": text(closed_form_quadrature_residual, 50),
        "gram_matrix_decimal": [[text(entry, 50) for entry in row] for row in gram],
        "gram_diagonal_decimal": [text(entry, 50) for entry in diagonal],
        "cholesky_pivots_decimal": [text(pivot, 50) for pivot in pivots],
        "leading_principal_minors_decimal": [text(minor, 50) for minor in leading_principal_minors],
        "minimum_leading_principal_minor": text(min(leading_principal_minors), 50),
        "all_leading_principal_minors_positive": all(minor > 0 for minor in leading_principal_minors),
        "minimum_cholesky_pivot": text(min(pivots), 50),
        "all_cholesky_pivots_positive": all(pivot > 0 for pivot in pivots),
        "symmetry_residual": text(symmetry_residual, 30),
        "witness_coefficients": [text(value, 20) for value in witness_coefficients],
        "witness_node_max_abs_residual": text(max(zero_node_residuals), 50),
        "witness_energy": text(witness_energy, 70),
        "witness_energy_positive": witness_energy > 0,
        "finite_polynomial_nullspace_coercive": all(pivot > 0 for pivot in pivots) and all(minor > 0 for minor in leading_principal_minors) and witness_energy > 0,
    }


def build_coercivity_candidate(p2506: dict[str, Any]) -> dict[str, Any]:
    symbolic = symbolic_coercivity_argument()
    gram = gram_coercivity_certificate()
    return {
        "frontier_atom_under_attack": "strict_damping_beta_eta_source",
        "p2506_conditional_selector_inherited": p2506.get("constant_flow_selected_if_selector_postulated") is True,
        "selector_type": "conditional minimum-roughness selector coercivity audit on a finite polynomial node-vanishing nullspace",
        "symbolic_coercivity_argument": symbolic,
        "finite_polynomial_gram_audit": gram,
        "finite_polynomial_nullspace_coercive": gram["finite_polynomial_nullspace_coercive"],
        "coercivity_is_selector_support_not_source_derivation": True,
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
## P2507/S1457 strict damping RG roughness-nullspace coercivity certificate

`P2507/S1457` strengthens the conditional P2506 minimum-roughness selector by auditing coercivity on a finite polynomial nullspace family.  For perturbations `p(ell)=R(ell) q(ell)`, `R(ell)=ell prod_{d=2}^{11}(ell-log(d))`, and `deg q<=3`, the roughness quadratic form `int_0^log(11) (p''(ell))^2 d ell` is evaluated by closed-form polynomial antiderivatives, with split quadrature retained only as a cross-check.  The resulting Gram audit has positive Cholesky pivots and positive leading principal minors.  Symbolically, zero roughness would force an affine perturbation, and the node constraints at `d=1,2` eliminate every nonzero affine perturbation.

This supports the conditional selector against a broader finite polynomial nullspace than the single P2505 witness.  It still does not derive the roughness action from nadsoliton dynamics, does not export `strict_damping_beta_eta_source`, and exports no bridge theorem, role-transfer theorem, QW-2191 discharge, physical-value generator, or ToE closure.
"""
    lag_section = """
## P2507/S1457 roughness-nullspace coercivity guard

`P2507/S1457` verifies that the postulated P2506 roughness functional is coercive on an audited finite polynomial node-vanishing nullspace using a closed-form Gram antiderivative, positive Cholesky pivots, positive leading principal minors, and a tight split-quadrature cross-check.  This is selector-support theorem-prep only: it does not turn the postulated action into a strict nonlinear compression-flow source or a role-bearing `L_total` term.
"""
    append_once(DOC_FILES["equation_sheet"], "## P2507/S1457 strict damping RG roughness-nullspace coercivity certificate", eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "## P2507/S1457 roughness-nullspace coercivity guard", lag_section)


def build_payload() -> dict[str, Any]:
    grep = rg_audit()
    sources = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    p2506 = theorem(sources["P2506_ROUGHNESS_SELECTOR"], "strict_damping_rg_minimum_roughness_selector_candidate_certificate")
    cert = build_coercivity_candidate(p2506)
    gram = cert["finite_polynomial_gram_audit"]
    theorem_export = {
        "theorem_name": "P2507_T1_strict_damping_rg_roughness_nullspace_coercivity_certificate",
        "audited_chain": ["P2503/S1453", "P2504/S1454", "P2505/S1455", "P2506/S1456"],
        "strict_damping_rg_roughness_nullspace_coercivity_certificate": cert,
        "frontier_atom_under_attack": cert["frontier_atom_under_attack"],
        "p2506_conditional_selector_inherited": cert["p2506_conditional_selector_inherited"],
        "basis_dimension": gram["basis_dimension"],
        "primary_integration_method": gram["primary_integration_method"],
        "closed_form_vs_split_quadrature_max_abs_residual": gram["closed_form_vs_split_quadrature_max_abs_residual"],
        "all_cholesky_pivots_positive": gram["all_cholesky_pivots_positive"],
        "all_leading_principal_minors_positive": gram["all_leading_principal_minors_positive"],
        "minimum_cholesky_pivot": gram["minimum_cholesky_pivot"],
        "minimum_leading_principal_minor": gram["minimum_leading_principal_minor"],
        "witness_node_max_abs_residual": gram["witness_node_max_abs_residual"],
        "witness_energy": gram["witness_energy"],
        "witness_energy_positive": gram["witness_energy_positive"],
        "finite_polynomial_nullspace_coercive": cert["finite_polynomial_nullspace_coercive"],
        "coercivity_is_selector_support_not_source_derivation": cert["coercivity_is_selector_support_not_source_derivation"],
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
            "P2507 proves coercivity only for the postulated roughness selector on an audited finite polynomial nullspace family.",
            "It does not derive the roughness action, delta, beta0, or strict_damping_beta_eta_source from strict nadsoliton dynamics.",
            "No bridge theorem, role-transfer theorem, QW-2191 discharge, role-bearing L_total, physical-value generator, or ToE closure is exported.",
        ],
        "next_honest_step": "Try to derive the roughness functional/source term from strict nadsoliton dynamics, or expand the coercivity audit toward a genuine functional-analytic theorem with stated admissible spaces and boundary conditions.",
    }
    gatekeepers = {
        "rg_audit_performed": grep["tool"] == "rg",
        "p2506_conditional_selector_inherited": theorem_export["p2506_conditional_selector_inherited"],
        "closed_form_quadrature_crosscheck_tight": mp.mpf(theorem_export["closed_form_vs_split_quadrature_max_abs_residual"]) < mp.mpf("1e-70"),
        "positive_cholesky_pivots": theorem_export["all_cholesky_pivots_positive"],
        "positive_leading_principal_minors": theorem_export["all_leading_principal_minors_positive"],
        "witness_nodes_still_vanish": mp.mpf(theorem_export["witness_node_max_abs_residual"]) < mp.mpf("1e-70"),
        "witness_energy_positive": theorem_export["witness_energy_positive"],
        "finite_polynomial_nullspace_coercive": theorem_export["finite_polynomial_nullspace_coercive"],
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
        "packet_id": "P2507",
        "stage_id": "S1457",
        "status": "STRICT_DAMPING_RG_ROUGHNESS_NULLSPACE_COERCIVITY_CERTIFICATE_CONDITIONAL_SELECTOR_SUPPORT_NO_SOURCE_EXPORT_NO_BRIDGE_THEOREM_NO_ROLE_TRANSFER_NO_QW2191_NO_TOE",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "strict_damping_rg_roughness_nullspace_coercivity_certificate": {
            "theorem_export": theorem_export,
            "source_fingerprints": {name: sha256_json(sources[name]) for name in sorted(sources)},
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["strict_damping_rg_roughness_nullspace_coercivity_certificate"]["theorem_export"]
    lines = [
        "# P2507/S1457 strict damping RG roughness-nullspace coercivity certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Result",
        "",
        f"- Frontier atom under attack: `{t['frontier_atom_under_attack']}`.",
        f"- P2506 conditional selector inherited: `{t['p2506_conditional_selector_inherited']}`.",
        f"- Audited polynomial nullspace basis dimension: `{t['basis_dimension']}`.",
        f"- Primary integration method: `{t['primary_integration_method']}`.",
        f"- Closed-form vs split-quadrature residual: `{t['closed_form_vs_split_quadrature_max_abs_residual']}`.",
        f"- All Gram Cholesky pivots positive: `{t['all_cholesky_pivots_positive']}`.",
        f"- All leading principal minors positive: `{t['all_leading_principal_minors_positive']}`.",
        f"- Minimum Cholesky pivot: `{t['minimum_cholesky_pivot']}`.",
        f"- Minimum leading principal minor: `{t['minimum_leading_principal_minor']}`.",
        f"- Witness node max residual: `{t['witness_node_max_abs_residual']}`.",
        f"- Witness roughness energy: `{t['witness_energy']}`.",
        f"- Finite polynomial nullspace coercive: `{t['finite_polynomial_nullspace_coercive']}`.",
        f"- Roughness action still postulated, not derived: `{t['roughness_action_still_postulated_not_derived']}`.",
        f"- Source theorem exported: `{t['strict_damping_beta_eta_source_exported']}`.",
        "",
        "## Negative controls",
        "",
        "This packet is conditional selector support. It audits roughness coercivity on a finite polynomial node-vanishing nullspace, but does not derive the roughness action from strict dynamics and exports no strict source atom, bridge theorem, role-transfer theorem, selector/QW-2191 closure, role-bearing L_total term, physical-value generator, or ToE closure.",
        "",
        "## Fingerprint",
        "",
        f"`{payload['strict_damping_rg_roughness_nullspace_coercivity_certificate']['theorem_fingerprint_sha256']}`",
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    GEN.mkdir(exist_ok=True)
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    write_markdown(payload)
    append_doc_sections()
    print(json.dumps(payload["strict_damping_rg_roughness_nullspace_coercivity_certificate"]["theorem_export"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
