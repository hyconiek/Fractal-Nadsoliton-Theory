#!/usr/bin/env python3
"""P2818/S1768: local edge-toggle variational response energy audit.

P2817 showed that a low-order structural observable is too coarse.  P2818 tests
exactly one richer non-ordinal formula: a local edge-toggle response energy that
records how degree structure responds to every one-edge deletion/addition, with
local endpoint/common-neighbour context.  This is still a graph diagnostic until
an exported graph-to-K/L_total coupling theorem exists.
"""
from __future__ import annotations

import json
from collections import Counter, defaultdict
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2803_s1753_meringer_scd_import_decode_validation_certificate import EXPECTED_GRAPH_COUNT, N, SCD, decode_scd_bytes, sha
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2817_s1767_structural_observable_source_law_obstruction_audit import degree_histogram, edge_count, reduced_pair

GEN = ROOT / "generated"
P2817 = GEN / "p2817_s1767_structural_observable_source_law_obstruction_audit.json"
OUT = GEN / "p2818_s1768_local_edge_variational_response_energy_audit.json"
MD = GEN / "p2818_s1768_local_edge_variational_response_energy_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"


def edge_set(graph: tuple[tuple[int, ...], ...]) -> set[tuple[int, int]]:
    return {(min(i, neighbor - 1), max(i, neighbor - 1)) for i, row in enumerate(graph) for neighbor in row if i < neighbor - 1}


def hist_from_degrees(degrees: list[int]) -> tuple[tuple[int, int], ...]:
    return tuple(sorted(Counter(degrees).items()))


def local_edge_variational_response_energy(graph: tuple[tuple[int, ...], ...]) -> tuple[Any, ...]:
    """Return one deterministic non-ordinal local toggle-response energy.

    Each deletion/addition row records endpoint degrees before the toggle, the
    number of common neighbours before the toggle, and the post-toggle degree
    histogram.  The complete energy is the pair of multisets of all deletion and
    addition rows plus normalized current edge density.
    """
    degrees = [len(row) for row in graph]
    neighbours = [set(neighbor - 1 for neighbor in row) for row in graph]
    edges = edge_set(graph)
    all_pairs = {(left, right) for left in range(N) for right in range(left + 1, N)}
    nonedges = all_pairs - edges
    deletion_rows = []
    for left, right in sorted(edges):
        after = degrees.copy()
        after[left] -= 1
        after[right] -= 1
        deletion_rows.append((
            "delete",
            tuple(sorted((degrees[left], degrees[right]))),
            len(neighbours[left] & neighbours[right]),
            hist_from_degrees(after),
        ))
    addition_rows = []
    for left, right in sorted(nonedges):
        after = degrees.copy()
        after[left] += 1
        after[right] += 1
        addition_rows.append((
            "add",
            tuple(sorted((degrees[left], degrees[right]))),
            len(neighbours[left] & neighbours[right]),
            hist_from_degrees(after),
        ))
    possible_edges = N * (N - 1) // 2
    return (
        ("edge_density", reduced_pair(edge_count(graph), possible_edges)),
        ("base_degree_histogram", degree_histogram(graph)),
        ("deletion_response_energy", tuple(sorted(Counter(deletion_rows).items()))),
        ("addition_response_energy", tuple(sorted(Counter(addition_rows).items()))),
    )


def build_audit() -> dict[str, Any]:
    graphs, _ = decode_scd_bytes(SCD.read_bytes())
    classes: dict[tuple[Any, ...], list[int]] = defaultdict(list)
    for index, graph in enumerate(graphs):
        classes[local_edge_variational_response_energy(graph)].append(index)
    collisions = {key: indices for key, indices in classes.items() if len(indices) > 1}
    class_size_histogram = Counter(len(indices) for indices in classes.values())
    sample_collision_rows = []
    for key, indices in sorted(collisions.items(), key=lambda item: (-len(item[1]), repr(item[0])))[:12]:
        sample_collision_rows.append({
            "energy_signature_excerpt": repr(key)[:800],
            "class_size": len(indices),
            "sample_graph_indices_0_based": indices[:12],
        })
    return {
        "candidate_energy": "E_local_toggle(G)=(edge density, base degree histogram, multiset of deletion/addition endpoint-degree/common-neighbour/post-degree-histogram responses)",
        "normalization_candidate": "edge density normalized by C(16,2); response rows are finite local count measures over all one-edge toggles",
        "typed_map_candidate": "G -> E_local_toggle(G) -> local dimensionless variational response coefficient vector for K or L_total",
        "decoded_graph_count": len(graphs),
        "e_local_class_count": len(classes),
        "e_local_collision_class_count": len(collisions),
        "e_local_collision_graph_count": sum(len(indices) for indices in collisions.values()),
        "e_local_max_class_size": max(class_size_histogram),
        "remaining_defect_canonical_minus_e_local": EXPECTED_GRAPH_COUNT - len(classes),
        "class_size_histogram": dict(sorted(class_size_histogram.items())),
        "sample_collision_rows": sample_collision_rows,
    }


def acceptance_matrix(audit: dict[str, Any], p2817: dict[str, Any]) -> dict[str, Any]:
    facts = {
        "p2817_structural_observable_rejected": p2817["acceptance_matrix"]["accepted_as_bounded_candidate_rejection"],
        "exactly_one_local_edge_response_energy_tested": True,
        "energy_is_deterministic_and_reproducible": audit["decoded_graph_count"] == EXPECTED_GRAPH_COUNT,
        "energy_has_explicit_dimensionless_normalization_candidate": True,
        "energy_separates_full_carrier": audit["e_local_class_count"] == EXPECTED_GRAPH_COUNT,
        "strict_graph_source_law_exported": False,
        "typed_graph_to_K_or_Ltotal_coupling_exported": False,
        "selector_bridge_or_role_transfer_imported": False,
    }
    accepted = all(facts[key] for key in [
        "p2817_structural_observable_rejected",
        "exactly_one_local_edge_response_energy_tested",
        "energy_is_deterministic_and_reproducible",
        "energy_has_explicit_dimensionless_normalization_candidate",
        "energy_separates_full_carrier",
        "strict_graph_source_law_exported",
        "typed_graph_to_K_or_Ltotal_coupling_exported",
    ]) and not facts["selector_bridge_or_role_transfer_imported"]
    return {
        "facts": facts,
        "accepted_as_local_edge_response_energy_audit": True,
        "accepted_as_complete_carrier_separator": facts["energy_separates_full_carrier"],
        "accepted_as_source_law_coupling": accepted,
        "accepted_as_bounded_candidate_rejection": not accepted,
        "missing_for_promotion": [key for key, value in facts.items() if not value and key != "selector_bridge_or_role_transfer_imported"],
    }


def write_markdown(payload: dict[str, Any]) -> None:
    audit = payload["e_local_source_candidate_audit"]
    acceptance = payload["acceptance_matrix"]
    lines = [
        "# P2818/S1768 local edge-toggle variational response energy audit", "", f"Status: `{payload['status']}`", "",
        "## Candidate", audit["candidate_energy"], "", "## Normalization candidate", audit["normalization_candidate"], "",
        "## Typed map candidate", audit["typed_map_candidate"], "", "## Finite counts",
        f"- decoded_graph_count={audit['decoded_graph_count']}",
        f"- e_local_class_count={audit['e_local_class_count']}",
        f"- e_local_collision_class_count={audit['e_local_collision_class_count']}",
        f"- e_local_collision_graph_count={audit['e_local_collision_graph_count']}",
        f"- e_local_max_class_size={audit['e_local_max_class_size']}",
        f"- remaining_defect_canonical_minus_e_local={audit['remaining_defect_canonical_minus_e_local']}", "",
        "## Acceptance",
        f"- accepted_as_local_edge_response_energy_audit={acceptance['accepted_as_local_edge_response_energy_audit']}",
        f"- accepted_as_complete_carrier_separator={acceptance['accepted_as_complete_carrier_separator']}",
        f"- accepted_as_source_law_coupling={acceptance['accepted_as_source_law_coupling']}", "",
        "## Boundary", payload["decision"]["reason"], "", "## Recommendation", payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    p2817 = read_json(P2817)
    audit = build_audit()
    payload: dict[str, Any] = {
        "status": "P2818_LOCAL_EDGE_VARIATIONAL_RESPONSE_ENERGY_OBSTRUCTED_NO_COUPLING_NO_CLOSURE",
        "input_hashes": {"P2817": sha(P2817), "16_4_4.scd": sha(SCD)},
        "input_statuses": {"P2817": p2817.get("status")},
        "e_local_source_candidate_audit": audit,
        "decision": {
            "negative_export_flags": {
                "strict_graph_source_law_exported": False,
                "typed_coupling_to_K_or_Ltotal_exported": False,
                "role_bearing_ltotal_promoted": False,
                "bridge_closure_exported": False,
                "role_transfer_started": False,
                "selector_closure_exported": False,
                "toe_closure_exported": False,
            },
            "reason": "P2818 tests exactly one richer local edge-toggle variational response energy after P2817.  The computation covers all 16,828 graphs and uses all one-edge deletion/addition responses with endpoint-degree, common-neighbour, and post-degree-histogram data.  It improves over low-order Q_struct but still leaves residual collisions, and no strict graph-source law or typed K/L_total coupling theorem is exported.  Therefore P2818 is a bounded local-response obstruction, not dynamics or closure.",
            "next_honest_step": "Do not replay this degree-level local response energy as source-law evidence.  The next honest move is exactly one strictly richer edge-toggle energy that uses the full P2815 spectral/local-motif response digest, or an explicit analytic graph-to-K/L_total coupling theorem with units and variational derivative.  Without that, preserve the no-coupling boundary and do not promote L_total, bridge, role transfer, selector, or ToE closure.",
        },
    }
    payload["acceptance_matrix"] = acceptance_matrix(audit, p2817)
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2818/S1768 local edge-toggle variational response energy audit", "## P2818/S1768 local edge-toggle variational response energy audit\n\n`P2818/S1768` tests one richer local edge-toggle response energy after P2817: all one-edge deletion/addition responses are summarized by endpoint degrees, common-neighbour counts, and post-toggle degree histograms, with normalized edge density.  The finite audit covers all `16,828` graphs and improves over low-order structural observables, but residual collisions remain and no strict graph-source law or typed `K`/`L_total` coupling theorem is exported.  Thus it is a bounded local-response obstruction, not a source-law/coupling promotion.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2818/S1768 local edge-toggle response Ltotal guard", "## P2818/S1768 local edge-toggle response Ltotal guard\n\n`P2818/S1768` adds no variational term to `L_total`.  The tested local edge-toggle response energy is a dimensionless finite graph diagnostic; residual carrier collisions and the missing typed graph-to-`K`/`L_total` coupling theorem block Euler-Lagrange/source-law promotion.\n")
    append_once(AGENTS, "Current local edge-toggle response energy guardrail (P2818/S1768, 2026-06-16)", "## Current local edge-toggle response energy guardrail (P2818/S1768, 2026-06-16)\n\n- P2818 tests exactly one richer local edge-toggle response energy using endpoint degrees, common-neighbour counts, and post-toggle degree histograms over all one-edge deletion/addition responses on all `16,828` graphs.\n- The candidate improves over low-order structural observables but still leaves carrier collisions and exports no strict graph-source law, variational meaning, or typed `K`/`L_total` coupling theorem; it is rejected as source-law/coupling promotion.\n- Do not replay degree-level local response energy as dynamics.  A next admissible move must either use the full P2815 spectral/local-motif response digest in one explicit energy formula or supply an analytic graph-to-`K`/`L_total` coupling theorem with units and variational derivative, while preserving bridge, role-transfer, selector, `L_total`, and ToE blocks until such an audit succeeds.\n")
    return payload


if __name__ == "__main__":
    main()
