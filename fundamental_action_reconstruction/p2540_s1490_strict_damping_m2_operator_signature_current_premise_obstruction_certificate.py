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
OUT = GEN / "p2540_s1490_strict_damping_m2_operator_signature_current_premise_obstruction_certificate.json"
MD = GEN / "p2540_s1490_strict_damping_m2_operator_signature_current_premise_obstruction_certificate.md"

SOURCE_FILES = {
    "P2509_VARIATIONAL_WELLPOSEDNESS": GEN / "p2509_s1459_strict_damping_rg_minimum_roughness_variational_wellposedness_certificate.json",
    "P2514_HIGHER_ORDER_NONIDENTIFIABILITY": GEN / "p2514_s1464_strict_damping_rg_higher_order_selector_nonidentifiability_certificate.json",
    "P2515_OPERATOR_SIGNATURE_ACCEPTANCE": GEN / "p2515_s1465_strict_damping_rg_operator_order_signature_acceptance_audit.json",
    "P2516_DUAL_KEY_ACCEPTANCE": GEN / "p2516_s1466_strict_damping_dual_key_source_acceptance_matrix.json",
    "P2530_FOUR_KEY_IRREDUNDANCY": GEN / "p2530_s1480_strict_damping_four_key_irredundancy_witness_certificate.json",
    "P2539_TOE_POTENTIAL_RECOMMENDATION": GEN / "p2539_s1489_strict_damping_toe_potential_recommendation_certificate.json",
}

mp.mp.dps = 110
DOMAIN = list(range(1, 12))
ORDER_RANGE = list(range(1, 11))
BASIS_DEGREES = list(range(7))
COUNTERMODEL_PAIR = (2, 3)


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
    return {"count": len(lines), "samples": lines[:50]}


def rg_audit() -> dict[str, Any]:
    patterns = {
        "new_packet": "P2540|S1490|m2 operator signature current-premise obstruction|current-premise non-entailment|m2 source obstruction",
        "intended_research_nonduplication": "m2_operator_signature_source|operator-signature source|biharmonic source|current-premise obstruction|non-entailment",
        "precursor_packets": "P2509|S1459|P2514|S1464|P2515|S1465|P2516|S1466|P2530|S1480|P2539|S1489",
        "source_key_language": "multiplicative_character_law_source|prime_log_proportionality_source|slope_value_or_prime_anchor_source|m2_operator_signature_source",
        "guardrails": "K_legacy_ont|K_strict_gate|role-transfer audit|QW-2191|ToE closure|selector guardrail",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def theorem(payload: dict[str, Any], container: str) -> dict[str, Any]:
    return payload.get(container, {}).get("theorem_export", {})


def text(value: mp.mpf, digits: int = 70) -> str:
    return mp.nstr(value, digits)


def poly_mul(left: list[mp.mpf], right: list[mp.mpf]) -> list[mp.mpf]:
    out = [mp.mpf("0")] * (len(left) + len(right) - 1)
    for i, left_value in enumerate(left):
        for j, right_value in enumerate(right):
            out[i + j] += left_value * right_value
    return out


def derivative_coefficients(coeffs: list[mp.mpf], order: int) -> list[mp.mpf]:
    current = coeffs[:]
    for _ in range(order):
        if len(current) <= 1:
            return [mp.mpf("0")]
        current = [mp.mpf(index + 1) * current[index + 1] for index in range(len(current) - 1)]
    return current or [mp.mpf("0")]


def integrate_product_closed_form(left: list[mp.mpf], right: list[mp.mpf], upper: mp.mpf) -> mp.mpf:
    product = poly_mul(left, right)
    return mp.fsum(coeff * upper ** (power + 1) / mp.mpf(power + 1) for power, coeff in enumerate(product))


def node_vanishing_r_coefficients() -> list[mp.mpf]:
    coeffs = [mp.mpf("1")]
    for d in DOMAIN:
        coeffs = poly_mul(coeffs, [-mp.log(d), mp.mpf("1")])
    return coeffs


def monomial_times(coeffs: list[mp.mpf], degree: int) -> list[mp.mpf]:
    return [mp.mpf("0")] * degree + coeffs


def tangent_basis() -> list[list[mp.mpf]]:
    r = node_vanishing_r_coefficients()
    return [monomial_times(r, degree) for degree in BASIS_DEGREES]


def derivative_gram_matrix(order: int) -> mp.matrix:
    upper = mp.log(11)
    basis = tangent_basis()
    matrix = mp.matrix(len(basis), len(basis))
    for i, left in enumerate(basis):
        left_d = derivative_coefficients(left, order)
        for j, right in enumerate(basis):
            right_d = derivative_coefficients(right, order)
            matrix[i, j] = integrate_product_closed_form(left_d, right_d, upper)
    return matrix


def cholesky_like_pivots(matrix: mp.matrix) -> list[mp.mpf]:
    work = mp.matrix(matrix)
    pivots: list[mp.mpf] = []
    for k in range(work.rows):
        pivot = work[k, k]
        pivots.append(pivot)
        if abs(pivot) < mp.mpf("1e-90"):
            continue
        for i in range(k + 1, work.rows):
            for j in range(k + 1, work.cols):
                work[i, j] -= work[i, k] * work[k, j] / pivot
    return pivots


def symbolic_order_rows() -> list[dict[str, Any]]:
    ell = sp.symbols("ell", real=True)
    delta = sp.Rational(14, 5) - 4 * sp.log(2)
    y0 = delta * ell
    rows = []
    for order in ORDER_RANGE:
        rows.append({
            "derivative_order_m": order,
            "functional": f"J_{order}[y]=int_0^log(11) (D^{order} y)^2 d ell",
            "euler_lagrange_order": 2 * order,
            "candidate_y0_euler_lagrange_residual": sp.sstr(sp.diff(y0, ell, 2 * order)),
            "candidate_y0_satisfies_euler_lagrange": sp.diff(y0, ell, 2 * order) == 0,
            "zero_mode_degree_bound": order - 1,
            "node_count": len(DOMAIN),
            "node_constraints_kill_zero_modes": len(DOMAIN) > order - 1,
            "matches_target_m2_signature": order == 2,
            "current_premise_model_candidate": True,
        })
    return rows


def extended_finite_gram_audit() -> dict[str, Any]:
    rows = []
    for order in ORDER_RANGE:
        matrix = derivative_gram_matrix(order)
        pivots = cholesky_like_pivots(matrix)
        min_pivot = min(pivots)
        max_pivot = max(pivots)
        rows.append({
            "derivative_order_m": order,
            "basis_family": "R(ell)*ell^k, k=0..6, R(ell)=prod_{d=1}^{11}(ell-log(d))",
            "basis_size": len(BASIS_DEGREES),
            "min_cholesky_like_pivot": text(min_pivot, 70),
            "max_cholesky_like_pivot": text(max_pivot, 70),
            "all_pivots_positive": all(pivot > 0 for pivot in pivots),
            "matches_target_m2_signature": order == 2,
        })
    return {
        "order_rows": rows,
        "basis_degrees": BASIS_DEGREES,
        "order_range": ORDER_RANGE,
        "all_orders_have_positive_tangent_gram_on_extended_basis": all(row["all_pivots_positive"] for row in rows),
        "extension_beyond_p2514_finite_audit": "P2514 finite Gram support was through m=1..6 on k=0..3; P2540 checks m=1..10 on k=0..6 for the current-premise non-entailment countermodel audit.",
    }


def current_premise_countermodels(symbolic_rows: list[dict[str, Any]], finite: dict[str, Any]) -> dict[str, Any]:
    finite_by_order = {row["derivative_order_m"]: row for row in finite["order_rows"]}
    rows = []
    for row in symbolic_rows:
        order = row["derivative_order_m"]
        passes = (
            row["candidate_y0_satisfies_euler_lagrange"]
            and row["node_constraints_kill_zero_modes"]
            and finite_by_order[order]["all_pivots_positive"]
        )
        rows.append({
            "derivative_order_m": order,
            "passes_current_source_free_premises": passes,
            "premises_used": [
                "strict damping node data lie on the affine y0(ell)=delta*ell line",
                "node-fixed derivative-only quadratic action J_m is nonnegative",
                "candidate y0 satisfies the Euler-Lagrange equation D^(2m)y=0",
                "node constraints kill polynomial zero modes for m<=10",
                "finite tangent Gram audit is positive on R(ell)*ell^k, k=0..6",
            ],
            "matches_target_m2_signature": order == 2,
            "m2_operator_signature_source_exported_in_this_model": order == 2 and False,
        })
    passing_orders = [row["derivative_order_m"] for row in rows if row["passes_current_source_free_premises"]]
    pair_rows = [row for row in rows if row["derivative_order_m"] in COUNTERMODEL_PAIR]
    return {
        "countermodel_rows": rows,
        "current_source_free_premise_passing_orders": passing_orders,
        "passing_order_count": len(passing_orders),
        "countermodel_pair_orders": list(COUNTERMODEL_PAIR),
        "countermodel_pair_rows": pair_rows,
        "m2_and_m3_both_pass_current_source_free_premises": all(row["passes_current_source_free_premises"] for row in pair_rows),
        "m2_is_unique_under_current_source_free_premises": passing_orders == [2],
        "current_source_free_premises_entail_m2": passing_orders == [2],
        "model_theoretic_nonentailment_statement": (
            "Because the same exported source-free premise schema has at least two satisfying derivative-order models "
            "(m=2 and m=3), it does not entail the proposition 'the strict operator signature is m=2'."
        ),
    }


def build_obstruction_certificate(sources: dict[str, dict[str, Any]]) -> dict[str, Any]:
    p2509 = theorem(sources["P2509_VARIATIONAL_WELLPOSEDNESS"], "strict_damping_rg_minimum_roughness_variational_wellposedness_certificate")
    p2514 = theorem(sources["P2514_HIGHER_ORDER_NONIDENTIFIABILITY"], "strict_damping_rg_higher_order_selector_nonidentifiability_certificate")
    p2515 = theorem(sources["P2515_OPERATOR_SIGNATURE_ACCEPTANCE"], "strict_damping_rg_operator_order_signature_acceptance_audit")
    p2516 = theorem(sources["P2516_DUAL_KEY_ACCEPTANCE"], "strict_damping_dual_key_source_acceptance_matrix")
    p2530 = theorem(sources["P2530_FOUR_KEY_IRREDUNDANCY"], "strict_damping_four_key_irredundancy_witness_certificate")
    p2539 = theorem(sources["P2539_TOE_POTENTIAL_RECOMMENDATION"], "strict_damping_toe_potential_recommendation_certificate")
    symbolic = symbolic_order_rows()
    finite = extended_finite_gram_audit()
    countermodels = current_premise_countermodels(symbolic, finite)
    return {
        "frontier_source_key_under_attack": "m2_operator_signature_source",
        "p2509_postulated_functional_wellposedness_inherited": p2509.get("minimum_roughness_problem_wellposed_for_postulated_functional") is True,
        "p2514_higher_order_nonidentifiability_inherited": p2514.get("higher_order_selector_tower_nonidentifiability_exported") is True,
        "p2515_m2_signature_identified_but_unsourced_inherited": p2515.get("p2506_roughness_m2_signature_identified_not_sourced") is True,
        "p2516_dual_key_boundary_inherited": p2516.get("dual_key_acceptance_normal_form_exported") is True,
        "p2530_four_key_irredundancy_inherited": p2530.get("four_key_irredundancy_witness_exported") is True,
        "p2539_next_step_recommendation_inherited": p2539.get("recommended_next_honest_step") == "prove_or_refute_one_strict_source_theorem_for_a_single_P2529_source_key_before_more_bookkeeping_layers",
        "symbolic_order_premise_rows": symbolic,
        "extended_finite_gram_audit": finite,
        "current_premise_countermodels": countermodels,
        "current_premise_nonentailment_of_m2_exported": not countermodels["current_source_free_premises_entail_m2"],
        "m2_operator_signature_source_route_refuted_for_current_source_free_premises": countermodels["m2_and_m3_both_pass_current_source_free_premises"],
        "scope_of_refutation": (
            "This refutes only the current source-free derivative-order route to m2. It does not prove that a future "
            "strict nadsoliton theorem cannot source m2 by adding a genuinely new symmetry-breaking/operator-order premise."
        ),
        "required_new_premise_class": "strict nadsoliton dynamics must export an operator-order selection principle that distinguishes m=2 from the admitted derivative-order tower",
    }


def append_doc_sections() -> None:
    eq_section = """
## P2540/S1490 strict damping m2 operator-signature current-premise obstruction certificate

`P2540/S1490` performs the source-key attack recommended by P2539, choosing the `m2_operator_signature_source` key.  The result is negative but proof-relevant: under the current exported source-free premises, the derivative-order choice `m=2` is not entailed.  The audit inherits P2509 well-posedness for the postulated roughness functional, P2514 higher-order nonidentifiability, P2515 identification-but-not-source of the `m=2` signature, P2516 dual-key acceptance, and P2530 four-key irredundancy.  It then checks an extended finite tangent Gram basis `R(ell)*ell^k`, `k=0..6`, for derivative orders `m=1..10`; every checked order passes the same node-fixed nonnegative quadratic-action premise schema.

Thus `m=2` and `m=3` are explicit current-premise countermodels: both satisfy the present source-free derivative-order premises, but only one is the desired biharmonic target.  Therefore the present route cannot export `m2_operator_signature_source`; a future theorem must add a genuinely new strict operator-order selection principle from nadsoliton dynamics.  No strict source key, source-obligation discharge, bridge completion, role-transfer theorem, QW-2191 discharge, role-bearing `L_total`, physical-value generator, or ToE closure is exported.
"""
    lag_section = """
## P2540/S1490 m2 operator-signature source obstruction guard

`P2540/S1490` attacks the `m2_operator_signature_source` key and finds a current-premise obstruction: the existing derivative-order and node-fixed variational premises admit `m=2` and `m=3` (indeed `m=1..10` in the audit) as source-free models.  So the `m=2` biharmonic signature remains a target, not a strict source theorem, until nadsoliton dynamics supply a new operator-order selection principle.
"""
    append_once(DOC_FILES["equation_sheet"], "## P2540/S1490 strict damping m2 operator-signature current-premise obstruction certificate", eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "## P2540/S1490 m2 operator-signature source obstruction guard", lag_section)


def build_payload() -> dict[str, Any]:
    grep = rg_audit()
    sources = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    cert = build_obstruction_certificate(sources)
    countermodels = cert["current_premise_countermodels"]
    finite = cert["extended_finite_gram_audit"]
    theorem_export = {
        "theorem_name": "P2540_T1_strict_damping_m2_operator_signature_current_premise_obstruction_certificate",
        "audited_chain": ["P2509/S1459", "P2514/S1464", "P2515/S1465", "P2516/S1466", "P2530/S1480", "P2539/S1489"],
        "strict_damping_m2_operator_signature_current_premise_obstruction_certificate": cert,
        "frontier_source_key_under_attack": cert["frontier_source_key_under_attack"],
        "p2509_postulated_functional_wellposedness_inherited": cert["p2509_postulated_functional_wellposedness_inherited"],
        "p2514_higher_order_nonidentifiability_inherited": cert["p2514_higher_order_nonidentifiability_inherited"],
        "p2515_m2_signature_identified_but_unsourced_inherited": cert["p2515_m2_signature_identified_but_unsourced_inherited"],
        "p2516_dual_key_boundary_inherited": cert["p2516_dual_key_boundary_inherited"],
        "p2530_four_key_irredundancy_inherited": cert["p2530_four_key_irredundancy_inherited"],
        "p2539_next_step_recommendation_inherited": cert["p2539_next_step_recommendation_inherited"],
        "finite_order_range": finite["order_range"],
        "finite_basis_degrees": finite["basis_degrees"],
        "all_orders_have_positive_tangent_gram_on_extended_basis": finite["all_orders_have_positive_tangent_gram_on_extended_basis"],
        "current_source_free_premise_passing_orders": countermodels["current_source_free_premise_passing_orders"],
        "passing_order_count": countermodels["passing_order_count"],
        "countermodel_pair_orders": countermodels["countermodel_pair_orders"],
        "m2_and_m3_both_pass_current_source_free_premises": countermodels["m2_and_m3_both_pass_current_source_free_premises"],
        "m2_is_unique_under_current_source_free_premises": countermodels["m2_is_unique_under_current_source_free_premises"],
        "current_source_free_premises_entail_m2": countermodels["current_source_free_premises_entail_m2"],
        "current_premise_nonentailment_of_m2_exported": cert["current_premise_nonentailment_of_m2_exported"],
        "m2_operator_signature_source_route_refuted_for_current_source_free_premises": cert["m2_operator_signature_source_route_refuted_for_current_source_free_premises"],
        "required_new_premise_class": cert["required_new_premise_class"],
        "scope_of_refutation": cert["scope_of_refutation"],
        "multiplicative_character_law_source_exported": False,
        "prime_log_proportionality_source_exported": False,
        "slope_value_or_prime_anchor_source_exported": False,
        "beta_eta_numeric_source_exported": False,
        "m2_operator_signature_source_exported": False,
        "strict_damping_beta_eta_source_exported": False,
        "source_obligation_discharge_exported": False,
        "damping_compression_bridge_component_ready": False,
        "full_bridge_theorem_exported": False,
        "role_transfer_theorem_exported": False,
        "selector_closure_exported": False,
        "qw2191_discharged_by_this_certificate": False,
        "role_bearing_ltotal_exported": False,
        "toe_closure_claimed": False,
        "not_licensed": [
            "P2540 refutes only the current source-free derivative-order route to m2; it does not refute future strict nadsoliton source premises.",
            "The m=2 operator signature remains identified as a target but unsourced as a theorem.",
            "No strict damping source key, bridge completion, role-transfer theorem, QW-2191 discharge, role-bearing L_total, physical-value generator, or ToE closure is exported.",
        ],
        "next_honest_step": "Construct a genuinely new strict operator-order selection premise from nadsoliton dynamics, or pivot to a different P2529 source key such as multiplicative_character_law_source.",
    }
    gatekeepers = {
        "rg_audit_performed": grep["tool"] == "rg",
        "all_required_precursors_inherited": all(theorem_export[key] for key in [
            "p2509_postulated_functional_wellposedness_inherited",
            "p2514_higher_order_nonidentifiability_inherited",
            "p2515_m2_signature_identified_but_unsourced_inherited",
            "p2516_dual_key_boundary_inherited",
            "p2530_four_key_irredundancy_inherited",
            "p2539_next_step_recommendation_inherited",
        ]),
        "extended_finite_gram_passes": theorem_export["all_orders_have_positive_tangent_gram_on_extended_basis"],
        "countermodel_pair_passes": theorem_export["m2_and_m3_both_pass_current_source_free_premises"],
        "nonentailment_recorded": theorem_export["current_premise_nonentailment_of_m2_exported"] and not theorem_export["current_source_free_premises_entail_m2"],
        "m2_source_not_exported": not theorem_export["m2_operator_signature_source_exported"],
        "negative_controls_preserved": not any(theorem_export[key] for key in [
            "multiplicative_character_law_source_exported",
            "prime_log_proportionality_source_exported",
            "slope_value_or_prime_anchor_source_exported",
            "beta_eta_numeric_source_exported",
            "m2_operator_signature_source_exported",
            "strict_damping_beta_eta_source_exported",
            "source_obligation_discharge_exported",
            "damping_compression_bridge_component_ready",
            "full_bridge_theorem_exported",
            "role_transfer_theorem_exported",
            "selector_closure_exported",
            "qw2191_discharged_by_this_certificate",
            "role_bearing_ltotal_exported",
            "toe_closure_claimed",
        ]),
    }
    return {
        "packet_id": "P2540",
        "stage_id": "S1490",
        "status": "STRICT_DAMPING_M2_OPERATOR_SIGNATURE_CURRENT_PREMISE_OBSTRUCTION_CERTIFICATE_CURRENT_ROUTE_REFUTED_NO_SOURCE_EXPORT_NO_BRIDGE_NO_ROLE_TRANSFER_NO_QW2191_NO_TOE",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "strict_damping_m2_operator_signature_current_premise_obstruction_certificate": {
            "theorem_export": theorem_export,
            "source_fingerprints": {name: sha256_json(sources[name]) for name in sorted(sources)},
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["strict_damping_m2_operator_signature_current_premise_obstruction_certificate"]["theorem_export"]
    lines = [
        "# P2540/S1490 strict damping m2 operator-signature current-premise obstruction certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Result",
        "",
        f"- Frontier source key under attack: `{t['frontier_source_key_under_attack']}`.",
        f"- P2539 source-key recommendation inherited: `{t['p2539_next_step_recommendation_inherited']}`.",
        f"- Extended finite order range: `{t['finite_order_range']}`.",
        f"- Extended finite basis degrees: `{t['finite_basis_degrees']}`.",
        f"- All audited orders have positive tangent Gram: `{t['all_orders_have_positive_tangent_gram_on_extended_basis']}`.",
        f"- Current source-free premise passing orders: `{t['current_source_free_premise_passing_orders']}`.",
        f"- Countermodel pair: `{t['countermodel_pair_orders']}`.",
        f"- m=2 and m=3 both pass current premises: `{t['m2_and_m3_both_pass_current_source_free_premises']}`.",
        f"- Current premises entail m=2: `{t['current_source_free_premises_entail_m2']}`.",
        f"- Current-premise non-entailment exported: `{t['current_premise_nonentailment_of_m2_exported']}`.",
        f"- m2 operator signature source exported: `{t['m2_operator_signature_source_exported']}`.",
        "",
        "## Interpretation",
        "",
        t["scope_of_refutation"],
        "",
        "## Negative controls",
        "",
        "This packet exports a current-premise obstruction only. It does not source any strict damping key, discharge source_obligation_discharge, complete the bridge, export role transfer, discharge QW-2191, produce role-bearing L_total, or claim ToE closure.",
        "",
        "## Fingerprint",
        "",
        f"`{payload['strict_damping_m2_operator_signature_current_premise_obstruction_certificate']['theorem_fingerprint_sha256']}`",
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    GEN.mkdir(exist_ok=True)
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    write_markdown(payload)
    append_doc_sections()
    print(json.dumps(payload["strict_damping_m2_operator_signature_current_premise_obstruction_certificate"]["theorem_export"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
