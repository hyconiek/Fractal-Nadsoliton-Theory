#!/usr/bin/env python3
"""P2838/S1788: formal variational derivative obstruction audit.

P2837 leaves a formal variational derivative theorem or localization map as the
next admissible single-premise attack.  P2838 attacks the variational-derivative
premise directly.  It distinguishes finite graph differences on adjacency bits
from an Euler-Lagrange/functional derivative into K or L_total and records the
missing chain-map premises.
"""
from __future__ import annotations

import json
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2803_s1753_meringer_scd_import_decode_validation_certificate import SCD, sha
from p2812_s1762_two_wl_refined_collision_audit import read_json

GEN = ROOT / "generated"
P2835 = GEN / "p2835_s1785_combined_witness_source_law_theorem_obligation_audit.json"
P2837 = GEN / "p2837_s1787_combined_graph_functional_typed_domain_codomain_audit.json"
OUT = GEN / "p2838_s1788_combined_graph_functional_variational_derivative_obstruction_audit.json"
MD = GEN / "p2838_s1788_combined_graph_functional_variational_derivative_obstruction_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

VERTEX_COUNT = 16
PAIR_COUNT = VERTEX_COUNT * (VERTEX_COUNT - 1) // 2
TWO_EDGE_PAIR_COUNT = PAIR_COUNT * (PAIR_COUNT - 1) // 2
FULL_CARRIER_COUNT = 16828

REQUIRED_VARIATIONAL_PREMISES = (
    "finite_difference_operator_on_graph_bits",
    "continuous_or_formal_field_variable",
    "localization_or_pullback_Aij_to_fields",
    "action_integral_or_density_embedding",
    "boundary_condition_and_integration_by_parts_rule",
    "target_independent_units_or_coupling_coefficient",
    "chain_rule_from_graph_bits_to_K_or_Ltotal",
)


def finite_difference_inventory() -> dict[str, Any]:
    return {
        "adjacency_bit_variable_count": PAIR_COUNT,
        "first_difference_slots_per_graph": PAIR_COUNT,
        "second_difference_slots_per_graph": TWO_EDGE_PAIR_COUNT,
        "p2833_first_variation_witness_available": True,
        "p2834_second_variation_witness_available_on_residuals": True,
        "finite_difference_layer": "discrete adjacency-bit hypercube",
    }


def candidate_variational_derivative_routes() -> list[dict[str, Any]]:
    routes = [
        {
            "candidate": "adjacency_bit_discrete_gradient",
            "proposed_derivative": "ΔF/ΔA_ij on the finite graph hypercube",
            "premises": {
                "finite_difference_operator_on_graph_bits": True,
                "continuous_or_formal_field_variable": False,
                "localization_or_pullback_Aij_to_fields": False,
                "action_integral_or_density_embedding": False,
                "boundary_condition_and_integration_by_parts_rule": False,
                "target_independent_units_or_coupling_coefficient": False,
                "chain_rule_from_graph_bits_to_K_or_Ltotal": False,
            },
            "obstruction": "This is a valid finite-difference object, but not an Euler-Lagrange derivative into K or L_total.",
        },
        {
            "candidate": "field_variation_delta_F_delta_phi",
            "proposed_derivative": "δF/δφ(x)",
            "premises": {
                "finite_difference_operator_on_graph_bits": True,
                "continuous_or_formal_field_variable": True,
                "localization_or_pullback_Aij_to_fields": False,
                "action_integral_or_density_embedding": False,
                "boundary_condition_and_integration_by_parts_rule": False,
                "target_independent_units_or_coupling_coefficient": False,
                "chain_rule_from_graph_bits_to_K_or_Ltotal": False,
            },
            "obstruction": "No exported A_ij -> φ(x) pullback or local density turns graph differences into field variations.",
        },
        {
            "candidate": "metric_variation_delta_F_delta_g",
            "proposed_derivative": "δF/δg_{μν}(x)",
            "premises": {
                "finite_difference_operator_on_graph_bits": True,
                "continuous_or_formal_field_variable": True,
                "localization_or_pullback_Aij_to_fields": False,
                "action_integral_or_density_embedding": False,
                "boundary_condition_and_integration_by_parts_rule": False,
                "target_independent_units_or_coupling_coefficient": False,
                "chain_rule_from_graph_bits_to_K_or_Ltotal": False,
            },
            "obstruction": "No graph-to-metric embedding or stress-energy variation theorem is exported for the combined witness.",
        },
        {
            "candidate": "kernel_parameter_variation_delta_F_delta_K",
            "proposed_derivative": "δF/δK(d) or ∂F/∂kernel-parameter",
            "premises": {
                "finite_difference_operator_on_graph_bits": True,
                "continuous_or_formal_field_variable": True,
                "localization_or_pullback_Aij_to_fields": False,
                "action_integral_or_density_embedding": False,
                "boundary_condition_and_integration_by_parts_rule": False,
                "target_independent_units_or_coupling_coefficient": False,
                "chain_rule_from_graph_bits_to_K_or_Ltotal": False,
            },
            "obstruction": "P2837 found no typed map from graph rows into K(d), so a graph-bit-to-kernel chain rule is unavailable.",
        },
    ]
    for route in routes:
        route["satisfied_premise_count"] = sum(bool(route["premises"][name]) for name in REQUIRED_VARIATIONAL_PREMISES)
        route["missing_premises"] = [name for name in REQUIRED_VARIATIONAL_PREMISES if not route["premises"][name]]
        route["accepted_as_formal_variational_derivative"] = not route["missing_premises"]
    return routes


def obstruction_summary(routes: list[dict[str, Any]]) -> dict[str, Any]:
    accepted = [route["candidate"] for route in routes if route["accepted_as_formal_variational_derivative"]]
    blocker_histogram: dict[str, int] = {name: 0 for name in REQUIRED_VARIATIONAL_PREMISES}
    for route in routes:
        for missing in route["missing_premises"]:
            blocker_histogram[missing] += 1
    return {
        "candidate_route_count": len(routes),
        "accepted_formal_variational_derivative_count": len(accepted),
        "accepted_routes": accepted,
        "blocker_histogram": blocker_histogram,
        "common_hard_blockers": [name for name, count in blocker_histogram.items() if count == len(routes)],
    }


def build_audit(p2835: dict[str, Any], p2837: dict[str, Any]) -> dict[str, Any]:
    p2835_sep = p2835["combined_witness_source_law_theorem_obligation_audit"]["combined_separator"]
    routes = candidate_variational_derivative_routes()
    return {
        "input_statuses_rechecked": {"P2835": p2835.get("status"), "P2837": p2837.get("status")},
        "combined_separator_rechecked": {
            "combined_class_count": p2835_sep["combined_class_count"],
            "combined_collision_class_count": p2835_sep["combined_collision_class_count"],
        },
        "finite_difference_inventory": finite_difference_inventory(),
        "required_variational_premises": list(REQUIRED_VARIATIONAL_PREMISES),
        "candidate_variational_derivative_routes": routes,
        "variational_obstruction_summary": obstruction_summary(routes),
    }


def acceptance_matrix(audit: dict[str, Any]) -> dict[str, Any]:
    sep = audit["combined_separator_rechecked"]
    summary = audit["variational_obstruction_summary"]
    inventory = audit["finite_difference_inventory"]
    facts = {
        "p2835_combined_separator_rechecked": sep["combined_class_count"] == FULL_CARRIER_COUNT and sep["combined_collision_class_count"] == 0,
        "finite_graph_difference_operators_available": inventory["p2833_first_variation_witness_available"] and inventory["p2834_second_variation_witness_available_on_residuals"],
        "accepted_formal_variational_derivative_exported": summary["accepted_formal_variational_derivative_count"] > 0,
        "typed_domain_codomain_available_from_p2837": False,
        "selector_bridge_or_role_transfer_imported": False,
    }
    accepted_as_formal_derivative = all([
        facts["p2835_combined_separator_rechecked"],
        facts["finite_graph_difference_operators_available"],
        facts["accepted_formal_variational_derivative_exported"],
        facts["typed_domain_codomain_available_from_p2837"],
        not facts["selector_bridge_or_role_transfer_imported"],
    ])
    return {
        "facts": facts,
        "accepted_as_finite_difference_inventory": facts["p2835_combined_separator_rechecked"] and facts["finite_graph_difference_operators_available"],
        "accepted_as_formal_variational_derivative": accepted_as_formal_derivative,
        "accepted_as_variational_derivative_no_go": not accepted_as_formal_derivative,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    audit = payload["combined_graph_functional_variational_derivative_obstruction_audit"]
    summary = audit["variational_obstruction_summary"]
    inventory = audit["finite_difference_inventory"]
    lines = [
        "# P2838/S1788 combined graph functional variational derivative obstruction audit", "", f"Status: `{payload['status']}`", "",
        "## Finite difference inventory",
        f"- adjacency_bit_variable_count={inventory['adjacency_bit_variable_count']}",
        f"- first_difference_slots_per_graph={inventory['first_difference_slots_per_graph']}",
        f"- second_difference_slots_per_graph={inventory['second_difference_slots_per_graph']}", "",
        "## Variational obstruction result",
        f"- candidate_route_count={summary['candidate_route_count']}",
        f"- accepted_formal_variational_derivative_count={summary['accepted_formal_variational_derivative_count']}",
        f"- common_hard_blockers={summary['common_hard_blockers']}", "",
        "## Acceptance",
        f"- accepted_as_finite_difference_inventory={payload['acceptance_matrix']['accepted_as_finite_difference_inventory']}",
        f"- accepted_as_formal_variational_derivative={payload['acceptance_matrix']['accepted_as_formal_variational_derivative']}",
        f"- accepted_as_variational_derivative_no_go={payload['acceptance_matrix']['accepted_as_variational_derivative_no_go']}", "",
        "## Boundary", payload["decision"]["reason"], "", "## Recommendation", payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    p2835 = read_json(P2835)
    p2837 = read_json(P2837)
    audit = build_audit(p2835, p2837)
    payload: dict[str, Any] = {
        "status": "P2838_VARIATIONAL_DERIVATIVE_OBSTRUCTION_NO_GO_NO_COUPLING_NO_CLOSURE",
        "input_hashes": {"P2835": sha(P2835), "P2837": sha(P2837), "16_4_4.scd": sha(SCD)},
        "combined_graph_functional_variational_derivative_obstruction_audit": audit,
        "decision": {
            "negative_export_flags": {
                "formal_variational_derivative_exported": False,
                "typed_K_or_Ltotal_map_exported": False,
                "strict_graph_source_law_exported": False,
                "typed_coupling_to_K_or_Ltotal_exported": False,
                "role_bearing_ltotal_promoted": False,
                "bridge_closure_exported": False,
                "role_transfer_started": False,
                "selector_closure_exported": False,
                "toe_closure_exported": False,
            },
            "reason": "P2838 attacks exactly one remaining premise: a formal variational derivative theorem.  Finite graph difference operators exist on 120 adjacency-bit variables, with 120 first-difference and 7140 second-difference slots per graph, but all four derivative routes fail to become Euler-Lagrange/functional derivatives into K or L_total.  Missing chain maps include localization/pullback from A_ij to fields, an action-density embedding, boundary/integration-by-parts rules, units/coupling, and graph-bit-to-K/L_total chain rules.",
            "next_honest_step": "Do not replay finite differences as variational calculus.  The next admissible proof-grade move should attack the remaining explicit evaluation/pullback/localization map from graph bits to field variables.  If no such localization map is exported, preserve the P2831-P2838 finite-difference/no-variational-derivative/no-coupling boundary and pivot away from graph-source promotion.",
        },
    }
    payload["acceptance_matrix"] = acceptance_matrix(audit)
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2838/S1788 combined graph functional variational derivative obstruction audit", "## P2838/S1788 combined graph functional variational derivative obstruction audit\n\n`P2838/S1788` attacks exactly one remaining theorem premise: a formal variational derivative for the combined graph functional.  Finite graph differences exist on `120` adjacency-bit variables (`120` first-difference slots and `7140` second-difference slots per graph), but candidate routes to `δF/δφ(x)`, `δF/δg_{μν}(x)`, `δF/δK(d)`, or an adjacency-bit gradient do not export an Euler-Lagrange/functional derivative into `K` or `L_total`.  Missing premises remain localization/pullback from graph bits to fields, action-density embedding, boundary/integration-by-parts rules, units/coupling, and graph-bit-to-`K`/`L_total` chain rules.  No strict graph-source law, role-bearing `L_total`, bridge closure, role transfer, selector closure, or ToE closure follows.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2838/S1788 variational derivative Ltotal guard", "## P2838/S1788 variational derivative Ltotal guard\n\n`P2838/S1788` adds no term and no Euler-Lagrange equation to `L_total`.  It confirms finite graph difference operators for the combined witness, but no localization/pullback, action-density embedding, boundary rule, units/coupling, or chain rule turns them into a formal variational derivative into `K` or `L_total`.\n")
    append_once(AGENTS, "Current variational-derivative obstruction guardrail (P2838/S1788, 2026-06-17)", "## Current variational-derivative obstruction guardrail (P2838/S1788, 2026-06-17)\n\n- P2838 attacks one remaining premise after P2837: a formal variational derivative theorem for the combined graph functional.\n- Finite graph differences on `120` adjacency bits exist, but no candidate route exports an Euler-Lagrange/functional derivative into `K` or `L_total`; localization/pullback, action-density embedding, boundary rules, units/coupling, and chain rules remain missing.\n- Do not promote P2838 to a strict graph-source law, `L_total`, bridge closure, role transfer, selector closure, or ToE closure.  A next admissible move must attack the explicit evaluation/pullback/localization map from graph bits to field variables, or preserve the no-variational-coupling boundary.\n")
    return payload


if __name__ == "__main__":
    main()
