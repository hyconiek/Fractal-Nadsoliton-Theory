#!/usr/bin/env python3
"""P2839/S1789: graph-bit localization/pullback obstruction audit.

P2838 leaves one honest next premise: an explicit evaluation/pullback/localization
map from graph adjacency bits to field variables.  P2839 attacks exactly that
premise and separates finite labeled-array conveniences from a strict, label-safe,
unit-bearing localization theorem into K or L_total.
"""
from __future__ import annotations

import math
import json
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2803_s1753_meringer_scd_import_decode_validation_certificate import SCD, sha
from p2812_s1762_two_wl_refined_collision_audit import read_json

GEN = ROOT / "generated"
P2835 = GEN / "p2835_s1785_combined_witness_source_law_theorem_obligation_audit.json"
P2838 = GEN / "p2838_s1788_combined_graph_functional_variational_derivative_obstruction_audit.json"
OUT = GEN / "p2839_s1789_graph_bit_localization_map_obstruction_audit.json"
MD = GEN / "p2839_s1789_graph_bit_localization_map_obstruction_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

VERTEX_COUNT = 16
EDGE_BIT_COUNT = VERTEX_COUNT * (VERTEX_COUNT - 1) // 2
FULL_CARRIER_COUNT = 16828

REQUIRED_LOCALIZATION_PREMISES = (
    "graph_bit_domain",
    "label_gauge_invariant_vertex_or_edge_coordinates",
    "target_field_or_spacetime_support",
    "pullback_formula_Aij_to_field_data",
    "locality_covariance_rule",
    "target_independent_units_or_measure",
    "compatibility_with_variational_chain_rule",
)


def finite_label_gauge_inventory() -> dict[str, Any]:
    return {
        "vertex_count": VERTEX_COUNT,
        "edge_bit_count": EDGE_BIT_COUNT,
        "label_permutation_count_exact": math.factorial(VERTEX_COUNT),
        "label_permutation_count_scientific": f"{math.factorial(VERTEX_COUNT):.6e}",
        "decoded_graphs_are_finite_unlabeled_isomorphism_class_representatives": True,
        "canonical_vertex_coordinate_source_exported": False,
    }


def candidate_localization_maps() -> list[dict[str, Any]]:
    candidates = [
        {
            "candidate": "label_index_to_fixed_points",
            "proposed_map": "vertex label i -> fixed point x_i; edge bit A_ij -> support between x_i,x_j",
            "premises": {
                "graph_bit_domain": True,
                "label_gauge_invariant_vertex_or_edge_coordinates": False,
                "target_field_or_spacetime_support": True,
                "pullback_formula_Aij_to_field_data": True,
                "locality_covariance_rule": False,
                "target_independent_units_or_measure": False,
                "compatibility_with_variational_chain_rule": False,
            },
            "obstruction": "Uses arbitrary decoded labels as coordinates; no strict label-gauge selector or coordinate source is exported.",
        },
        {
            "candidate": "edge_orbit_to_local_source_bins",
            "proposed_map": "edge orbits/features -> unlabeled source bins",
            "premises": {
                "graph_bit_domain": True,
                "label_gauge_invariant_vertex_or_edge_coordinates": True,
                "target_field_or_spacetime_support": False,
                "pullback_formula_Aij_to_field_data": False,
                "locality_covariance_rule": False,
                "target_independent_units_or_measure": False,
                "compatibility_with_variational_chain_rule": False,
            },
            "obstruction": "Can be label-safe as finite bins, but no spacetime/field support or pullback formula is exported.",
        },
        {
            "candidate": "graphon_step_function",
            "proposed_map": "adjacency matrix -> step function W_G(u,v) on [0,1]^2",
            "premises": {
                "graph_bit_domain": True,
                "label_gauge_invariant_vertex_or_edge_coordinates": False,
                "target_field_or_spacetime_support": True,
                "pullback_formula_Aij_to_field_data": True,
                "locality_covariance_rule": False,
                "target_independent_units_or_measure": False,
                "compatibility_with_variational_chain_rule": False,
            },
            "obstruction": "Graphon coordinates depend on vertex ordering unless quotient data or a canonical ordering source is exported; also lacks units/coupling.",
        },
        {
            "candidate": "spectral_embedding_localization",
            "proposed_map": "graph spectrum/eigenvectors -> coordinates/support",
            "premises": {
                "graph_bit_domain": True,
                "label_gauge_invariant_vertex_or_edge_coordinates": False,
                "target_field_or_spacetime_support": True,
                "pullback_formula_Aij_to_field_data": False,
                "locality_covariance_rule": False,
                "target_independent_units_or_measure": False,
                "compatibility_with_variational_chain_rule": False,
            },
            "obstruction": "Eigenvector sign/degeneracy/order choices are not a strict localization theorem and no field pullback is exported.",
        },
    ]
    for candidate in candidates:
        candidate["satisfied_premise_count"] = sum(bool(candidate["premises"][name]) for name in REQUIRED_LOCALIZATION_PREMISES)
        candidate["missing_premises"] = [name for name in REQUIRED_LOCALIZATION_PREMISES if not candidate["premises"][name]]
        candidate["accepted_as_localization_map"] = not candidate["missing_premises"]
    return candidates


def obstruction_summary(candidates: list[dict[str, Any]]) -> dict[str, Any]:
    accepted = [candidate["candidate"] for candidate in candidates if candidate["accepted_as_localization_map"]]
    blocker_histogram = {name: 0 for name in REQUIRED_LOCALIZATION_PREMISES}
    for candidate in candidates:
        for missing in candidate["missing_premises"]:
            blocker_histogram[missing] += 1
    return {
        "candidate_count": len(candidates),
        "accepted_localization_map_count": len(accepted),
        "accepted_localization_maps": accepted,
        "blocker_histogram": blocker_histogram,
        "common_hard_blockers": [name for name, count in blocker_histogram.items() if count == len(candidates)],
    }


def build_audit(p2835: dict[str, Any], p2838: dict[str, Any]) -> dict[str, Any]:
    p2835_sep = p2835["combined_witness_source_law_theorem_obligation_audit"]["combined_separator"]
    candidates = candidate_localization_maps()
    return {
        "input_statuses_rechecked": {"P2835": p2835.get("status"), "P2838": p2838.get("status")},
        "combined_separator_rechecked": {
            "combined_class_count": p2835_sep["combined_class_count"],
            "combined_collision_class_count": p2835_sep["combined_collision_class_count"],
        },
        "finite_label_gauge_inventory": finite_label_gauge_inventory(),
        "required_localization_premises": list(REQUIRED_LOCALIZATION_PREMISES),
        "candidate_localization_maps": candidates,
        "localization_obstruction_summary": obstruction_summary(candidates),
    }


def acceptance_matrix(audit: dict[str, Any]) -> dict[str, Any]:
    sep = audit["combined_separator_rechecked"]
    summary = audit["localization_obstruction_summary"]
    facts = {
        "p2835_combined_separator_rechecked": sep["combined_class_count"] == FULL_CARRIER_COUNT and sep["combined_collision_class_count"] == 0,
        "at_least_one_localization_candidate_tested": summary["candidate_count"] > 0,
        "accepted_localization_map_exported": summary["accepted_localization_map_count"] > 0,
        "canonical_vertex_coordinate_source_exported": audit["finite_label_gauge_inventory"]["canonical_vertex_coordinate_source_exported"],
        "selector_bridge_or_role_transfer_imported": False,
    }
    accepted = all([
        facts["p2835_combined_separator_rechecked"],
        facts["at_least_one_localization_candidate_tested"],
        facts["accepted_localization_map_exported"],
        facts["canonical_vertex_coordinate_source_exported"],
        not facts["selector_bridge_or_role_transfer_imported"],
    ])
    return {
        "facts": facts,
        "accepted_as_localization_obstruction_audit": facts["p2835_combined_separator_rechecked"] and facts["at_least_one_localization_candidate_tested"],
        "accepted_as_evaluation_pullback_localization_map": accepted,
        "accepted_as_localization_no_go": not accepted,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    audit = payload["graph_bit_localization_map_obstruction_audit"]
    inventory = audit["finite_label_gauge_inventory"]
    summary = audit["localization_obstruction_summary"]
    lines = [
        "# P2839/S1789 graph-bit localization map obstruction audit", "", f"Status: `{payload['status']}`", "",
        "## Finite label-gauge inventory",
        f"- vertex_count={inventory['vertex_count']}",
        f"- edge_bit_count={inventory['edge_bit_count']}",
        f"- label_permutation_count_exact={inventory['label_permutation_count_exact']}", "",
        "## Localization obstruction result",
        f"- candidate_count={summary['candidate_count']}",
        f"- accepted_localization_map_count={summary['accepted_localization_map_count']}",
        f"- common_hard_blockers={summary['common_hard_blockers']}", "",
        "## Acceptance",
        f"- accepted_as_localization_obstruction_audit={payload['acceptance_matrix']['accepted_as_localization_obstruction_audit']}",
        f"- accepted_as_evaluation_pullback_localization_map={payload['acceptance_matrix']['accepted_as_evaluation_pullback_localization_map']}",
        f"- accepted_as_localization_no_go={payload['acceptance_matrix']['accepted_as_localization_no_go']}", "",
        "## Boundary", payload["decision"]["reason"], "", "## Recommendation", payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    p2835 = read_json(P2835)
    p2838 = read_json(P2838)
    audit = build_audit(p2835, p2838)
    payload: dict[str, Any] = {
        "status": "P2839_LOCALIZATION_PULLBACK_OBSTRUCTION_NO_GO_NO_COUPLING_NO_CLOSURE",
        "input_hashes": {"P2835": sha(P2835), "P2838": sha(P2838), "16_4_4.scd": sha(SCD)},
        "graph_bit_localization_map_obstruction_audit": audit,
        "decision": {
            "negative_export_flags": {
                "evaluation_pullback_localization_map_exported": False,
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
            "reason": "P2839 attacks exactly one remaining P2838 premise: an evaluation/pullback/localization map from graph bits to field variables.  The finite label-gauge inventory has 16 vertices, 120 edge bits, and 16! possible labelings; current artifacts do not export a canonical vertex-coordinate source.  Four localization candidates are audited, and none satisfies all premises.  Label-index and graphon maps depend on arbitrary ordering, orbit bins lack field support, and spectral embeddings lack strict sign/degeneracy handling, units, and pullback/coupling rules.",
            "next_honest_step": "Do not replay graph-bit localization names without a new strict coordinate/localization source.  The next honest move is a post-graph-source state-map reconciliation/no-new-live-frontier certificate for this lane, unless a genuinely new strict localization object or coupling theorem is supplied.",
        },
    }
    payload["acceptance_matrix"] = acceptance_matrix(audit)
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2839/S1789 graph-bit localization map obstruction audit", "## P2839/S1789 graph-bit localization map obstruction audit\n\n`P2839/S1789` attacks exactly one remaining P2838 premise: an evaluation/pullback/localization map from graph bits to field variables.  The finite label-gauge inventory has `16` vertices, `120` edge bits, and `16!` possible labelings; no canonical vertex-coordinate source is exported.  Four localization candidates are audited (`label index -> fixed points`, edge-orbit bins, graphon step function, and spectral embedding), and none satisfies all premises.  Label-index and graphon maps depend on arbitrary ordering, orbit bins lack field support, and spectral embeddings lack strict sign/degeneracy handling, units, and pullback/coupling rules.  No strict graph-source law, role-bearing `L_total`, bridge closure, role transfer, selector closure, or ToE closure follows.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2839/S1789 localization map Ltotal guard", "## P2839/S1789 localization map Ltotal guard\n\n`P2839/S1789` adds no term and no localization map to `L_total`.  Candidate graph-bit-to-field localization routes remain blocked by label gauge, missing field support, missing pullback formula, missing locality/covariance, missing units/measure, and missing variational chain-rule compatibility.\n")
    append_once(AGENTS, "Current graph-bit localization map obstruction guardrail (P2839/S1789, 2026-06-17)", "## Current graph-bit localization map obstruction guardrail (P2839/S1789, 2026-06-17)\n\n- P2839 attacks one remaining P2838 premise: an evaluation/pullback/localization map from graph bits to field variables.\n- Four localization candidates are audited and all fail: current artifacts do not export a canonical vertex-coordinate source, field support, pullback formula, locality/covariance rule, units/measure, or variational chain-rule compatibility.\n- Do not promote P2839 to a strict graph-source law, `L_total`, bridge closure, role transfer, selector closure, or ToE closure.  Unless a new strict localization object or coupling theorem is supplied, the next move should be a post-graph-source state-map/no-new-live-frontier reconciliation for this lane.\n")
    return payload


if __name__ == "__main__":
    main()
