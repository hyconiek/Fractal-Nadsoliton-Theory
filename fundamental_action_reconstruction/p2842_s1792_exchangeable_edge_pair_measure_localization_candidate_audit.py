#!/usr/bin/env python3
"""P2842/S1792: exchangeable edge-pair measure localization candidate audit.

P2841 allows reopening only with one genuinely new strict typed object/source or
provider.  P2842 introduces exactly one candidate object: the S_16-exchangeable
edge-pair counting/probability measure on the fixed 16-node 4-regular graph
carrier.  The audit checks whether this label-safe finite measure can serve as a
localization/coupling object for the combined graph functional.
"""
from __future__ import annotations

import json
from fractions import Fraction
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2803_s1753_meringer_scd_import_decode_validation_certificate import SCD, decode_scd_bytes, sha
from p2812_s1762_two_wl_refined_collision_audit import read_json

GEN = ROOT / "generated"
P2841 = GEN / "p2841_s1791_post_p2840_broad_state_map_intake_gate.json"
OUT = GEN / "p2842_s1792_exchangeable_edge_pair_measure_localization_candidate_audit.json"
MD = GEN / "p2842_s1792_exchangeable_edge_pair_measure_localization_candidate_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

VERTEX_COUNT = 16
DEGREE = 4
EDGE_COUNT = VERTEX_COUNT * DEGREE // 2
PAIR_COUNT = VERTEX_COUNT * (VERTEX_COUNT - 1) // 2
EXPECTED_GRAPH_COUNT = 16828

REQUIRED_LOCALIZATION_OBJECT_PREMISES = (
    "new_typed_object_declared",
    "label_gauge_invariant",
    "nonconstant_on_full_carrier",
    "separates_or_refines_combined_functional",
    "field_or_spacetime_support_exported",
    "target_independent_units_or_measure_exported",
    "coupling_or_variational_chain_rule_exported",
)


def fraction_payload(value: Fraction) -> dict[str, int | str]:
    return {"numerator": value.numerator, "denominator": value.denominator, "decimal": f"{float(value):.12g}"}


def edge_count(graph: tuple[tuple[int, ...], ...]) -> int:
    return sum(len(neighbors) for neighbors in graph) // 2


def build_measure_candidate() -> dict[str, Any]:
    edge_density = Fraction(EDGE_COUNT, PAIR_COUNT)
    nonedge_density = Fraction(PAIR_COUNT - EDGE_COUNT, PAIR_COUNT)
    return {
        "candidate_object": "mu_edge^exch(G)",
        "definition": "The S_16-exchangeable probability measure on unordered vertex pairs with mass split by edge/non-edge occupancy on the fixed 16-node 4-regular carrier.",
        "vertex_count": VERTEX_COUNT,
        "regular_degree": DEGREE,
        "unordered_pair_count": PAIR_COUNT,
        "edge_count_per_graph": EDGE_COUNT,
        "edge_density": fraction_payload(edge_density),
        "nonedge_density": fraction_payload(nonedge_density),
        "label_gauge_invariant_by_construction": True,
        "finite_probability_measure_exported": True,
    }


def finite_carrier_check() -> dict[str, Any]:
    graphs, _ = decode_scd_bytes(SCD.read_bytes())
    counts = [edge_count(graph) for graph in graphs]
    distinct_counts = sorted(set(counts))
    return {
        "decoded_graph_count": len(graphs),
        "expected_graph_count": EXPECTED_GRAPH_COUNT,
        "distinct_edge_counts": distinct_counts,
        "edge_count_histogram": {str(count): counts.count(count) for count in distinct_counts},
        "measure_constant_on_full_carrier": len(distinct_counts) == 1,
        "constant_edge_count": distinct_counts[0] if len(distinct_counts) == 1 else None,
    }


def candidate_premise_matrix(carrier: dict[str, Any]) -> dict[str, Any]:
    premises = {
        "new_typed_object_declared": True,
        "label_gauge_invariant": True,
        "nonconstant_on_full_carrier": not carrier["measure_constant_on_full_carrier"],
        "separates_or_refines_combined_functional": False,
        "field_or_spacetime_support_exported": False,
        "target_independent_units_or_measure_exported": True,
        "coupling_or_variational_chain_rule_exported": False,
    }
    return {
        "premises": premises,
        "missing_premises": [name for name in REQUIRED_LOCALIZATION_OBJECT_PREMISES if not premises[name]],
        "accepted_as_strict_localization_or_coupling_object": all(premises.values()),
    }


def build_audit(p2841: dict[str, Any]) -> dict[str, Any]:
    candidate = build_measure_candidate()
    carrier = finite_carrier_check()
    matrix = candidate_premise_matrix(carrier)
    return {
        "input_statuses_rechecked": {"P2841": p2841.get("status")},
        "candidate_measure": candidate,
        "finite_carrier_check": carrier,
        "required_localization_object_premises": list(REQUIRED_LOCALIZATION_OBJECT_PREMISES),
        "candidate_premise_matrix": matrix,
    }


def acceptance_matrix(audit: dict[str, Any]) -> dict[str, Any]:
    premises = audit["candidate_premise_matrix"]["premises"]
    facts = {
        "new_single_candidate_object_tested": premises["new_typed_object_declared"],
        "candidate_is_label_gauge_invariant": premises["label_gauge_invariant"],
        "finite_probability_measure_exported": audit["candidate_measure"]["finite_probability_measure_exported"],
        "candidate_nonconstant_on_full_carrier": premises["nonconstant_on_full_carrier"],
        "candidate_separates_or_refines_combined_functional": premises["separates_or_refines_combined_functional"],
        "field_or_spacetime_support_exported": premises["field_or_spacetime_support_exported"],
        "coupling_or_variational_chain_rule_exported": premises["coupling_or_variational_chain_rule_exported"],
        "selector_bridge_or_role_transfer_imported": False,
    }
    accepted = all([
        facts["new_single_candidate_object_tested"],
        facts["candidate_is_label_gauge_invariant"],
        facts["candidate_nonconstant_on_full_carrier"],
        facts["candidate_separates_or_refines_combined_functional"],
        facts["field_or_spacetime_support_exported"],
        facts["coupling_or_variational_chain_rule_exported"],
        not facts["selector_bridge_or_role_transfer_imported"],
    ])
    return {
        "facts": facts,
        "accepted_as_finite_label_safe_measure_candidate": facts["new_single_candidate_object_tested"] and facts["candidate_is_label_gauge_invariant"] and facts["finite_probability_measure_exported"],
        "accepted_as_strict_localization_or_coupling_object": accepted,
        "accepted_as_bounded_candidate_no_go": not accepted,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    audit = payload["exchangeable_edge_pair_measure_localization_candidate_audit"]
    carrier = audit["finite_carrier_check"]
    matrix = audit["candidate_premise_matrix"]
    lines = [
        "# P2842/S1792 exchangeable edge-pair measure localization candidate audit", "", f"Status: `{payload['status']}`", "",
        "## Candidate", audit["candidate_measure"]["definition"], "",
        "## Finite carrier check",
        f"- decoded_graph_count={carrier['decoded_graph_count']}",
        f"- distinct_edge_counts={carrier['distinct_edge_counts']}",
        f"- measure_constant_on_full_carrier={carrier['measure_constant_on_full_carrier']}", "",
        "## Premise result",
        f"- accepted_as_strict_localization_or_coupling_object={matrix['accepted_as_strict_localization_or_coupling_object']}",
        f"- missing_premises={matrix['missing_premises']}", "",
        "## Acceptance",
        f"- accepted_as_finite_label_safe_measure_candidate={payload['acceptance_matrix']['accepted_as_finite_label_safe_measure_candidate']}",
        f"- accepted_as_strict_localization_or_coupling_object={payload['acceptance_matrix']['accepted_as_strict_localization_or_coupling_object']}",
        f"- accepted_as_bounded_candidate_no_go={payload['acceptance_matrix']['accepted_as_bounded_candidate_no_go']}", "",
        "## Boundary", payload["decision"]["reason"], "", "## Recommendation", payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    p2841 = read_json(P2841)
    audit = build_audit(p2841)
    payload: dict[str, Any] = {
        "status": "P2842_EXCHANGEABLE_EDGE_PAIR_MEASURE_LOCALIZATION_CANDIDATE_NO_GO_NO_COUPLING_NO_CLOSURE",
        "input_hashes": {"P2841": sha(P2841), "16_4_4.scd": sha(SCD)},
        "exchangeable_edge_pair_measure_localization_candidate_audit": audit,
        "decision": {
            "negative_export_flags": {
                "strict_localization_object_exported": False,
                "coupling_or_variational_chain_rule_exported": False,
                "strict_graph_source_law_exported": False,
                "typed_coupling_to_K_or_Ltotal_exported": False,
                "role_bearing_ltotal_promoted": False,
                "bridge_closure_exported": False,
                "role_transfer_started": False,
                "selector_closure_exported": False,
                "toe_closure_exported": False,
            },
            "reason": "P2842 introduces exactly one new candidate after P2841: the S_16-exchangeable edge-pair probability measure on the fixed 16-node 4-regular carrier.  It is label-gauge invariant and finite, with edge density 32/120=4/15.  However every decoded graph has the same edge count 32, so the measure is constant on the full 16,828-graph carrier; it does not refine the combined witness and exports no field/spacetime support or coupling/variational chain rule.  It is accepted only as a finite label-safe measure candidate, not as a strict localization/coupling object.",
            "next_honest_step": "Do not replay exchangeable edge-count measures or other carrier-constant summaries.  A next admissible move needs a genuinely nonconstant, label-safe strict localization object with field support and a coupling/chain-rule theorem; otherwise preserve the P2841-P2842 no-new-live-frontier/no-localization boundary.",
        },
    }
    payload["acceptance_matrix"] = acceptance_matrix(audit)
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2842/S1792 exchangeable edge-pair measure localization candidate audit", "## P2842/S1792 exchangeable edge-pair measure localization candidate audit\n\n`P2842/S1792` introduces exactly one candidate after P2841: the `S_16`-exchangeable edge-pair probability measure on the fixed `16`-node `4`-regular carrier.  It is label-gauge invariant and finite, with edge density `32/120 = 4/15`, but every decoded graph has the same edge count `32`; therefore the measure is constant on the full `16,828`-graph carrier, does not refine the combined witness, and exports no field/spacetime support or coupling/variational chain rule.  No strict localization object, graph-source law, role-bearing `L_total`, bridge closure, role transfer, selector closure, or ToE closure follows.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2842/S1792 exchangeable edge measure Ltotal guard", "## P2842/S1792 exchangeable edge measure Ltotal guard\n\n`P2842/S1792` adds no term to `L_total`.  The exchangeable edge-pair measure is label-safe but carrier-constant on 16-node 4-regular graphs and lacks field support, coupling, and variational chain-rule data.\n")
    append_once(AGENTS, "Current exchangeable edge-pair measure guardrail (P2842/S1792, 2026-06-18)", "## Current exchangeable edge-pair measure guardrail (P2842/S1792, 2026-06-18)\n\n- P2842 tests exactly one new candidate after P2841: the `S_16`-exchangeable edge-pair probability measure on the fixed `16`-node `4`-regular carrier.\n- The measure is label-gauge invariant and finite, but all decoded graphs have edge count `32`, so it is constant on the full `16,828`-graph carrier; it does not refine the combined graph witness and exports no field support or coupling/variational chain rule.\n- Do not promote P2842 to a strict localization object, graph-source law, `L_total`, bridge closure, role transfer, selector closure, or ToE closure.  A next admissible move requires a genuinely nonconstant, label-safe strict localization object with field support and a coupling/chain-rule theorem, or preservation of the no-live-frontier boundary.\n")
    return payload


if __name__ == "__main__":
    main()
