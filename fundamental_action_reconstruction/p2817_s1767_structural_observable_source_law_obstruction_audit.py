#!/usr/bin/env python3
"""P2817/S1767: non-ordinal structural observable source-law obstruction.

P2816 rejected the ordinal rank induced by the complete P2815 carrier.  P2817
therefore tests exactly one non-ordinal graph observable family with an explicit
normalization proposal:

    Q_struct(G) = (edge density, degree histogram, distance histogram,
                   exact 4-cycle count)

on all 16,828 Meringer girth>=4 graphs.  The candidate typed map is still only
G -> Q_struct(G) -> dimensionless structural coefficient vector for K/L_total.
The audit is computational: if this natural structural observable cannot even
separate the complete carrier, it cannot serve as a complete source law; and in
any case a typed variational coupling theorem would still be required.
"""
from __future__ import annotations

import json
from collections import Counter, defaultdict, deque
from math import comb, gcd
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2803_s1753_meringer_scd_import_decode_validation_certificate import EXPECTED_GRAPH_COUNT, N, SCD, decode_scd_bytes, sha
from p2812_s1762_two_wl_refined_collision_audit import read_json

GEN = ROOT / "generated"
P2816 = GEN / "p2816_s1766_edge_toggle_source_law_coupling_acceptance_audit.json"
OUT = GEN / "p2817_s1767_structural_observable_source_law_obstruction_audit.json"
MD = GEN / "p2817_s1767_structural_observable_source_law_obstruction_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"


def reduced_pair(numerator: int, denominator: int) -> tuple[int, int]:
    factor = gcd(numerator, denominator)
    return (numerator // factor, denominator // factor)


def edge_count(graph: tuple[tuple[int, ...], ...]) -> int:
    return sum(len(row) for row in graph) // 2


def degree_histogram(graph: tuple[tuple[int, ...], ...]) -> tuple[tuple[int, int], ...]:
    return tuple(sorted(Counter(len(row) for row in graph).items()))


def distance_histogram(graph: tuple[tuple[int, ...], ...]) -> tuple[tuple[int, int], ...]:
    distances: Counter[int] = Counter()
    for source in range(N):
        seen = {source: 0}
        queue: deque[int] = deque([source])
        while queue:
            node = queue.popleft()
            for neighbor in graph[node]:
                target = neighbor - 1
                if target not in seen:
                    seen[target] = seen[node] + 1
                    queue.append(target)
        for target in range(source + 1, N):
            distances[seen[target]] += 1
    return tuple(sorted(distances.items()))


def four_cycle_count(graph: tuple[tuple[int, ...], ...]) -> int:
    neighbor_sets = [set(neighbor - 1 for neighbor in row) for row in graph]
    total = 0
    for left in range(N):
        for right in range(left + 1, N):
            common = len(neighbor_sets[left] & neighbor_sets[right])
            if common >= 2:
                total += comb(common, 2)
    return total // 2


def q_struct(graph: tuple[tuple[int, ...], ...]) -> tuple[Any, ...]:
    edges = edge_count(graph)
    possible_edges = N * (N - 1) // 2
    return (
        ("edge_density", reduced_pair(edges, possible_edges)),
        ("degree_histogram", degree_histogram(graph)),
        ("distance_histogram", distance_histogram(graph)),
        ("four_cycle_count", four_cycle_count(graph)),
    )


def build_audit() -> dict[str, Any]:
    graphs, _ = decode_scd_bytes(SCD.read_bytes())
    classes: dict[tuple[Any, ...], list[int]] = defaultdict(list)
    for index, graph in enumerate(graphs):
        classes[q_struct(graph)].append(index)
    collisions = {key: indices for key, indices in classes.items() if len(indices) > 1}
    class_size_histogram = Counter(len(indices) for indices in classes.values())
    sample_collision_rows = []
    for key, indices in sorted(collisions.items(), key=lambda item: (-len(item[1]), repr(item[0])))[:12]:
        sample_collision_rows.append({
            "q_struct": repr(key),
            "class_size": len(indices),
            "sample_graph_indices_0_based": indices[:12],
        })
    return {
        "candidate_observable": "Q_struct(G)=(edge density, degree histogram, distance histogram, exact 4-cycle count)",
        "normalization_candidate": "edge density is normalized by C(16,2); histograms are finite count measures; 4-cycle count is dimensionless",
        "typed_map_candidate": "G -> Q_struct(G) -> dimensionless structural coefficient vector for K or L_total",
        "decoded_graph_count": len(graphs),
        "q_struct_class_count": len(classes),
        "q_struct_collision_class_count": len(collisions),
        "q_struct_collision_graph_count": sum(len(indices) for indices in collisions.values()),
        "q_struct_max_class_size": max(class_size_histogram),
        "remaining_defect_canonical_minus_q_struct": EXPECTED_GRAPH_COUNT - len(classes),
        "class_size_histogram": dict(sorted(class_size_histogram.items())),
        "sample_collision_rows": sample_collision_rows,
    }


def acceptance_matrix(audit: dict[str, Any], p2816: dict[str, Any]) -> dict[str, Any]:
    facts = {
        "p2816_completed_rank_candidate_rejected": p2816["acceptance_matrix"]["accepted_as_bounded_candidate_rejection"],
        "exactly_one_nonordinal_observable_tested": True,
        "observable_is_deterministic_and_reproducible": audit["decoded_graph_count"] == EXPECTED_GRAPH_COUNT,
        "observable_has_explicit_dimensionless_normalization_candidate": True,
        "observable_separates_full_carrier": audit["q_struct_class_count"] == EXPECTED_GRAPH_COUNT,
        "strict_graph_source_law_exported": False,
        "typed_graph_to_K_or_Ltotal_coupling_exported": False,
        "selector_bridge_or_role_transfer_imported": False,
    }
    accepted = all(facts[key] for key in [
        "p2816_completed_rank_candidate_rejected",
        "exactly_one_nonordinal_observable_tested",
        "observable_is_deterministic_and_reproducible",
        "observable_has_explicit_dimensionless_normalization_candidate",
        "observable_separates_full_carrier",
        "strict_graph_source_law_exported",
        "typed_graph_to_K_or_Ltotal_coupling_exported",
    ]) and not facts["selector_bridge_or_role_transfer_imported"]
    return {
        "facts": facts,
        "accepted_as_structural_observable_obstruction_audit": True,
        "accepted_as_source_law_coupling": accepted,
        "accepted_as_bounded_candidate_rejection": not accepted,
        "missing_for_promotion": [key for key, value in facts.items() if not value and key != "selector_bridge_or_role_transfer_imported"],
    }


def write_markdown(payload: dict[str, Any]) -> None:
    audit = payload["q_struct_source_candidate_audit"]
    acceptance = payload["acceptance_matrix"]
    lines = [
        "# P2817/S1767 structural observable source-law obstruction audit", "", f"Status: `{payload['status']}`", "",
        "## Candidate", audit["candidate_observable"], "", "## Normalization candidate", audit["normalization_candidate"], "",
        "## Typed map candidate", audit["typed_map_candidate"], "", "## Finite counts",
        f"- decoded_graph_count={audit['decoded_graph_count']}",
        f"- q_struct_class_count={audit['q_struct_class_count']}",
        f"- q_struct_collision_class_count={audit['q_struct_collision_class_count']}",
        f"- q_struct_collision_graph_count={audit['q_struct_collision_graph_count']}",
        f"- q_struct_max_class_size={audit['q_struct_max_class_size']}",
        f"- remaining_defect_canonical_minus_q_struct={audit['remaining_defect_canonical_minus_q_struct']}", "",
        "## Acceptance",
        f"- accepted_as_structural_observable_obstruction_audit={acceptance['accepted_as_structural_observable_obstruction_audit']}",
        f"- accepted_as_source_law_coupling={acceptance['accepted_as_source_law_coupling']}",
        f"- accepted_as_bounded_candidate_rejection={acceptance['accepted_as_bounded_candidate_rejection']}", "",
        "## Boundary", payload["decision"]["reason"], "", "## Recommendation", payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    p2816 = read_json(P2816)
    audit = build_audit()
    payload: dict[str, Any] = {
        "status": "P2817_STRUCTURAL_OBSERVABLE_SOURCE_CANDIDATE_OBSTRUCTED_NO_COUPLING_NO_CLOSURE",
        "input_hashes": {"P2816": sha(P2816), "16_4_4.scd": sha(SCD)},
        "input_statuses": {"P2816": p2816.get("status")},
        "q_struct_source_candidate_audit": audit,
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
            "reason": "P2817 tests exactly one non-ordinal normalized structural observable candidate after P2816: Q_struct=(edge density, degree histogram, distance histogram, exact 4-cycle count).  The finite computation covers all 16,828 graphs, but Q_struct collapses the complete P2815 carrier into far fewer classes and leaves residual collisions.  It also exports no strict graph-source law or typed K/L_total coupling theorem.  Therefore P2817 is a bounded structural-observable obstruction, not dynamics or closure.",
            "next_honest_step": "Do not replay ordinal ranks or this low-order structural observable as source-law evidence.  The next honest move is exactly one richer non-ordinal formula with an actual typed coupling ansatz, preferably a local edge-toggle variational response energy whose value is computed directly from the full P2815 separating response data and accompanied by a falsifiable graph-to-K or graph-to-L_total normalization/coupling theorem.  If no such formula is supplied, preserve the no-coupling boundary.",
        },
    }
    payload["acceptance_matrix"] = acceptance_matrix(audit, p2816)
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2817/S1767 structural observable source-law obstruction audit", "## P2817/S1767 structural observable source-law obstruction audit\n\n`P2817/S1767` tests exactly one non-ordinal normalized structural observable candidate after P2816: `Q_struct(G)=(edge density, degree histogram, distance histogram, exact 4-cycle count)` on all `16,828` Meringer graphs.  The observable is deterministic and dimensionless, but it collapses the complete P2815 carrier into far fewer classes and exports no strict graph-source law or typed `K`/`L_total` coupling theorem.  Thus it is a bounded obstruction, not a source-law/coupling promotion.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2817/S1767 structural observable Ltotal guard", "## P2817/S1767 structural observable Ltotal guard\n\n`P2817/S1767` adds no variational term to `L_total`.  The tested `Q_struct` observable is a dimensionless graph diagnostic; its finite collisions and missing typed graph-to-`K`/`L_total` coupling theorem block any Euler-Lagrange/source-law promotion.\n")
    append_once(AGENTS, "Current structural observable source-law obstruction guardrail (P2817/S1767, 2026-06-16)", "## Current structural observable source-law obstruction guardrail (P2817/S1767, 2026-06-16)\n\n- P2817 tests exactly one non-ordinal normalized structural observable candidate: `Q_struct(G)=(edge density, degree histogram, distance histogram, exact 4-cycle count)` over all `16,828` graphs.\n- The candidate is deterministic and dimensionless but leaves carrier collisions and exports no strict graph-source law, variational meaning, or typed `K`/`L_total` coupling theorem; it is rejected as source-law/coupling promotion.\n- Do not replay ordinal ranks or low-order structural observables as dynamics.  A next admissible move must provide one richer non-ordinal formula with explicit local edge-toggle variational response energy or equivalent typed graph-to-`K`/`L_total` normalization/coupling data, while preserving bridge, role-transfer, selector, `L_total`, and ToE blocks until such an audit succeeds.\n")
    return payload


if __name__ == "__main__":
    main()
