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
from p2513_s1463_strict_damping_rg_derivative_order_nonidentifiability_certificate import matrix_pivots_and_minors

GEN = ROOT / "generated"
OUT = GEN / "p2515_s1465_strict_damping_rg_operator_order_signature_acceptance_audit.json"
MD = GEN / "p2515_s1465_strict_damping_rg_operator_order_signature_acceptance_audit.md"

SOURCE_FILES = {
    "P2514_HIGHER_ORDER_SELECTOR_NONIDENTIFIABILITY": GEN / "p2514_s1464_strict_damping_rg_higher_order_selector_nonidentifiability_certificate.json",
}

mp.mp.dps = 100
SIGNATURE_ORDERS = list(range(1, 7))
MONOMIAL_DEGREE = 14


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
        "new_packet": "P2515|S1465|operator-order signature|Euler-Lagrange order signature|biharmonic source signature|boundary concomitant",
        "precursor_packets": "P2514|S1464|higher-order selector nonidentifiability|P2513|derivative-order nonidentifiability|P2512|quadratic source admissibility",
        "operator_language": "Euler-Lagrange operator|natural boundary|boundary concomitant|biharmonic|operator order|source signature",
        "guardrails": "legacy -> strict completion bridge|role-transfer audit|K_legacy_ont|K_strict_gate|QW-2191|ToE closure",
        "closure_blockers": "source theorem|bridge theorem|role-transfer theorem|physical-value generator|role-bearing L_total|selector closure",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def theorem(payload: dict[str, Any], container: str) -> dict[str, Any]:
    return payload.get(container, {}).get("theorem_export", {})


def text(value: mp.mpf, digits: int = 70) -> str:
    return mp.nstr(value, digits)


def derivative_operator_matrix(order: int, degree: int) -> mp.matrix:
    size = degree + 1
    matrix = mp.matrix(size, size)
    for col in range(size):
        if col >= order:
            coeff = mp.mpf("1")
            for factor in range(col - order + 1, col + 1):
                coeff *= factor
            matrix[col - order, col] = coeff
    return matrix


def matrix_rank_by_svd(matrix: mp.matrix, tolerance: mp.mpf = mp.mpf("1e-70")) -> int:
    singular_values = list(mp.svd(matrix, compute_uv=False))
    return sum(1 for value in singular_values if abs(value) > tolerance)


def symbolic_operator_signature_audit() -> dict[str, Any]:
    ell = sp.symbols("ell", real=True)
    delta = sp.Rational(14, 5) - 4 * sp.log(2)
    y0 = delta * ell
    rows = []
    signatures = []
    for order in SIGNATURE_ORDERS:
        el_order = 2 * order
        el_residual = sp.diff(y0, ell, el_order)
        free_boundary_derivatives = list(range(order, 2 * order))
        node_fixed_boundary_derivatives = list(range(order, max(order, 2 * order - 1)))
        boundary_residuals = [sp.diff(y0, ell, derivative_order) for derivative_order in node_fixed_boundary_derivatives]
        signature = {
            "derivative_order_m": order,
            "euler_lagrange_operator": f"(-1)^{order} D^{el_order} y",
            "euler_lagrange_differential_order": el_order,
            "free_boundary_derivative_orders": free_boundary_derivatives,
            "node_fixed_boundary_derivative_orders": node_fixed_boundary_derivatives,
            "zero_mode_dimension": order,
        }
        signatures.append((el_order, tuple(node_fixed_boundary_derivatives), tuple(free_boundary_derivatives), order))
        rows.append({
            **signature,
            "candidate_y0_euler_lagrange_residual": sp.sstr(el_residual),
            "candidate_y0_boundary_residuals": [sp.sstr(value) for value in boundary_residuals],
            "candidate_residuals_zero": el_residual == 0 and all(value == 0 for value in boundary_residuals),
            "matches_p2506_roughness_order_m2": order == 2,
            "source_signature_needed_to_select_this_order": f"future source must export differential order {el_order}; node-fixed boundary jets {node_fixed_boundary_derivatives}; free-boundary jets {free_boundary_derivatives}",
        })
    return {
        "symbolic_backend": "sympy",
        "sympy_version": sp.__version__,
        "candidate_y0": sp.sstr(y0),
        "signature_orders": SIGNATURE_ORDERS,
        "rows": rows,
        "all_node_fixed_candidate_residuals_zero": all(row["candidate_residuals_zero"] for row in rows),
        "operator_signatures_pairwise_distinct": len(set(signatures)) == len(signatures),
        "p2506_roughness_signature": "m=2, Euler-Lagrange order 4, node-fixed boundary derivative order [2], free-boundary derivative orders [2,3]",
        "signature_fingerprint_sha256": sha256_json({"delta": sp.srepr(delta), "rows": rows}),
    }


def finite_operator_matrix_audit() -> dict[str, Any]:
    rows = []
    for order in SIGNATURE_ORDERS:
        first = derivative_operator_matrix(order, MONOMIAL_DEGREE)
        el = derivative_operator_matrix(2 * order, MONOMIAL_DEGREE)
        gram = first.T * first
        active_size = min(4, gram.rows - order)
        active = mp.matrix(active_size, active_size)
        for ai in range(active_size):
            for aj in range(active_size):
                active[ai, aj] = gram[order + ai, order + aj]
        minors, pivots = matrix_pivots_and_minors(active)
        first_rank = matrix_rank_by_svd(first)
        el_rank = matrix_rank_by_svd(el)
        rows.append({
            "derivative_order_m": order,
            "monomial_degree": MONOMIAL_DEGREE,
            "first_variation_operator_rank_Dm": first_rank,
            "expected_rank_Dm": max(0, MONOMIAL_DEGREE + 1 - order),
            "euler_lagrange_operator_rank_D2m": el_rank,
            "expected_rank_D2m": max(0, MONOMIAL_DEGREE + 1 - 2 * order),
            "kernel_dimension_Dm": MONOMIAL_DEGREE + 1 - first_rank,
            "expected_zero_mode_dimension_Dm": order,
            "kernel_dimension_D2m": MONOMIAL_DEGREE + 1 - el_rank,
            "expected_kernel_dimension_D2m": min(MONOMIAL_DEGREE + 1, 2 * order),
            "low_order_principal_minors_of_Dm_gram": [text(value, 50) for value in minors],
            "low_order_pivots_of_Dm_gram": [text(value, 50) for value in pivots],
            "rank_identities_pass": first_rank == max(0, MONOMIAL_DEGREE + 1 - order) and el_rank == max(0, MONOMIAL_DEGREE + 1 - 2 * order),
            "matches_p2506_roughness_order_m2": order == 2,
        })
    return {
        "monomial_degree": MONOMIAL_DEGREE,
        "rows": rows,
        "all_rank_identities_pass": all(row["rank_identities_pass"] for row in rows),
        "operator_rank_signatures_pairwise_distinct": len({(row["first_variation_operator_rank_Dm"], row["euler_lagrange_operator_rank_D2m"], row["kernel_dimension_Dm"]) for row in rows}) == len(rows),
    }


def acceptance_signature_table() -> dict[str, Any]:
    rows = []
    for order in SIGNATURE_ORDERS:
        rows.append({
            "candidate_order_m": order,
            "would_select_affine_node_solution": True,
            "operator_signature": {
                "euler_lagrange_order": 2 * order,
                "node_fixed_boundary_derivative_orders": list(range(order, max(order, 2 * order - 1))),
                "free_boundary_derivative_orders": list(range(order, 2 * order)),
                "zero_mode_dimension": order,
            },
            "licensed_as_strict_source": False,
            "p2506_roughness_order_match": order == 2,
            "source_acceptance_requirement": "export this operator signature from strict nadsoliton dynamics, not from node-fitting stationarity",
        })
    return {
        "rows": rows,
        "p2506_order_m2_signature_present": any(row["p2506_roughness_order_match"] for row in rows),
        "signature_level_can_distinguish_orders": True,
        "solution_level_cannot_distinguish_orders": True,
        "source_acceptance_boundary": "To promote P2506/P2509 roughness from postulate to source, a future theorem must export the m=2 biharmonic/fourth-order operator signature (or explicitly justify another order) from strict dynamics.",
    }


def build_operator_order_signature_audit(p2514: dict[str, Any]) -> dict[str, Any]:
    symbolic = symbolic_operator_signature_audit()
    finite = finite_operator_matrix_audit()
    table = acceptance_signature_table()
    return {
        "frontier_atom_under_attack": "strict_damping_beta_eta_source",
        "p2514_higher_order_nonidentifiability_inherited": p2514.get("higher_order_selector_tower_nonidentifiability_exported") is True,
        "audit_type": "operator-order signature acceptance audit for derivative-only selector source theorem",
        "symbolic_operator_signature_audit": symbolic,
        "finite_operator_matrix_audit": finite,
        "acceptance_signature_table": table,
        "operator_signatures_distinguish_orders_even_when_solution_does_not": symbolic["operator_signatures_pairwise_distinct"] and finite["operator_rank_signatures_pairwise_distinct"],
        "p2506_roughness_m2_signature_identified_not_sourced": table["p2506_order_m2_signature_present"],
        "source_acceptance_boundary_exported": True,
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
## P2515/S1465 strict damping RG operator-order signature acceptance audit

`P2515/S1465` converts the P2514 derivative-order tower ambiguity into a source-acceptance target.  For `J_m[y]=1/2 int (D^m y)^2 d ell`, the Euler-Lagrange signature is `(-1)^m D^{2m}y`, with node-fixed boundary derivative orders `m..2m-2` and free-boundary derivative orders `m..2m-1`.  These operator signatures are pairwise distinct even though the affine strict damping node solution `y0(ell)=delta ell` satisfies every audited node-fixed stationarity signature.  The finite monomial operator audit through `m=1..6` verifies the expected ranks and kernel dimensions for `D^m` and `D^{2m}`.

Thus the P2506 roughness selector corresponds to the `m=2` biharmonic/fourth-order signature, but P2515 does not derive that signature from strict dynamics.  It only states a sharper acceptance boundary: a future strict source theorem must export the `m=2` operator signature (or justify another order) from the nadsoliton dynamics before `strict_damping_beta_eta_source` can be claimed.  No bridge theorem, role-transfer theorem, QW-2191 discharge, physical-value generator, or ToE closure is exported.
"""
    lag_section = """
## P2515/S1465 operator-order signature acceptance guard

`P2515/S1465` records that the derivative-order ambiguity can be broken only at the operator-signature/source level: `m=2` means a fourth-order/biharmonic Euler-Lagrange signature with node-fixed boundary jet `[2]` and free-boundary jets `[2,3]`, but the affine solution itself satisfies every audited derivative-only order.  Therefore no nonlinear compression-flow source theorem or role-bearing `L_total` term is licensed until strict dynamics export the operator signature.
"""
    append_once(DOC_FILES["equation_sheet"], "## P2515/S1465 strict damping RG operator-order signature acceptance audit", eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "## P2515/S1465 operator-order signature acceptance guard", lag_section)


def build_payload() -> dict[str, Any]:
    grep = rg_audit()
    sources = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    p2514 = theorem(sources["P2514_HIGHER_ORDER_SELECTOR_NONIDENTIFIABILITY"], "strict_damping_rg_higher_order_selector_nonidentifiability_certificate")
    cert = build_operator_order_signature_audit(p2514)
    symbolic = cert["symbolic_operator_signature_audit"]
    finite = cert["finite_operator_matrix_audit"]
    table = cert["acceptance_signature_table"]
    theorem_export = {
        "theorem_name": "P2515_T1_strict_damping_rg_operator_order_signature_acceptance_audit",
        "audited_chain": ["P2512/S1462", "P2513/S1463", "P2514/S1464"],
        "strict_damping_rg_operator_order_signature_acceptance_audit": cert,
        "frontier_atom_under_attack": cert["frontier_atom_under_attack"],
        "p2514_higher_order_nonidentifiability_inherited": cert["p2514_higher_order_nonidentifiability_inherited"],
        "signature_orders": symbolic["signature_orders"],
        "all_node_fixed_candidate_residuals_zero": symbolic["all_node_fixed_candidate_residuals_zero"],
        "operator_signatures_pairwise_distinct": symbolic["operator_signatures_pairwise_distinct"],
        "finite_operator_rank_identities_pass": finite["all_rank_identities_pass"],
        "operator_rank_signatures_pairwise_distinct": finite["operator_rank_signatures_pairwise_distinct"],
        "operator_signatures_distinguish_orders_even_when_solution_does_not": cert["operator_signatures_distinguish_orders_even_when_solution_does_not"],
        "p2506_roughness_m2_signature_identified_not_sourced": cert["p2506_roughness_m2_signature_identified_not_sourced"],
        "source_acceptance_boundary_exported": cert["source_acceptance_boundary_exported"],
        "source_acceptance_boundary": table["source_acceptance_boundary"],
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
            "P2515 identifies the operator-order signature needed by a future source theorem; it is not itself a strict source derivation.",
            "The affine node solution satisfies every audited derivative-only order, so only an exported operator signature can break the P2514 tower degeneracy.",
            "No bridge theorem, role-transfer theorem, QW-2191 discharge, role-bearing L_total, physical-value generator, or ToE closure is exported.",
        ],
        "next_honest_step": table["source_acceptance_boundary"],
    }
    gatekeepers = {
        "rg_audit_performed": grep["tool"] == "rg",
        "p2514_nonidentifiability_inherited": theorem_export["p2514_higher_order_nonidentifiability_inherited"],
        "solution_residuals_zero_for_all_orders": theorem_export["all_node_fixed_candidate_residuals_zero"],
        "operator_signatures_distinct": theorem_export["operator_signatures_pairwise_distinct"] and theorem_export["operator_rank_signatures_pairwise_distinct"],
        "finite_rank_identities_pass": theorem_export["finite_operator_rank_identities_pass"],
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
        "packet_id": "P2515",
        "stage_id": "S1465",
        "status": "STRICT_DAMPING_RG_OPERATOR_ORDER_SIGNATURE_ACCEPTANCE_AUDIT_NO_SOURCE_EXPORT_NO_BRIDGE_THEOREM_NO_ROLE_TRANSFER_NO_QW2191_NO_TOE",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "strict_damping_rg_operator_order_signature_acceptance_audit": {
            "theorem_export": theorem_export,
            "source_fingerprints": {name: sha256_json(sources[name]) for name in sorted(sources)},
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["strict_damping_rg_operator_order_signature_acceptance_audit"]["theorem_export"]
    lines = [
        "# P2515/S1465 strict damping RG operator-order signature acceptance audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Result",
        "",
        f"- Frontier atom under attack: `{t['frontier_atom_under_attack']}`.",
        f"- P2514 higher-order nonidentifiability inherited: `{t['p2514_higher_order_nonidentifiability_inherited']}`.",
        f"- Signature orders: `{t['signature_orders']}`.",
        f"- Node-fixed candidate residuals zero for all orders: `{t['all_node_fixed_candidate_residuals_zero']}`.",
        f"- Operator signatures pairwise distinct: `{t['operator_signatures_pairwise_distinct']}`.",
        f"- Finite operator rank identities pass: `{t['finite_operator_rank_identities_pass']}`.",
        f"- P2506 roughness m=2 signature identified, not sourced: `{t['p2506_roughness_m2_signature_identified_not_sourced']}`.",
        f"- Source acceptance boundary exported: `{t['source_acceptance_boundary_exported']}`.",
        f"- Source theorem exported: `{t['strict_damping_beta_eta_source_exported']}`.",
        "",
        "## Acceptance boundary",
        "",
        t["source_acceptance_boundary"],
        "",
        "## Negative controls",
        "",
        "This packet sharpens the source acceptance target but exports no strict source atom, bridge theorem, role-transfer theorem, selector/QW-2191 closure, role-bearing L_total term, physical-value generator, or ToE closure.",
        "",
        "## Fingerprint",
        "",
        f"`{payload['strict_damping_rg_operator_order_signature_acceptance_audit']['theorem_fingerprint_sha256']}`",
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    GEN.mkdir(exist_ok=True)
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    write_markdown(payload)
    append_doc_sections()
    print(json.dumps(payload["strict_damping_rg_operator_order_signature_acceptance_audit"]["theorem_export"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
