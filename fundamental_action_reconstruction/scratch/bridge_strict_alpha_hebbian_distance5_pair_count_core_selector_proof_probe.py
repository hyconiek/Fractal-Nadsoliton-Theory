#!/usr/bin/env python3
"""Scratch probe: distance-5 pair-count core selector proof for d5.

The combinatorial-score probe showed that the supplied d5/Hebbian margin can be
replayed by an integer score on pair-distance histograms.  This probe tests an
even smaller core: does the single count h_5, the number of active pairs at
folded distance 5 on Z_12, already isolate the d5 teacher orbit?

Finite result:
- Among all 792 five-active-node supports, the maximum h_5 count is 4.
- Exactly 12 supports attain h_5=4, and they are exactly the d5 translates.
- Graph proof: folded-distance-5 edges form one 12-cycle because gcd(5,12)=1;
  any 5 selected vertices in a cycle induce at most 4 such edges, with equality
  exactly on 5-vertex paths, giving the 12 d5 translates.
- Therefore the previous integer Hebbian score contains a simpler core selector:
  maximize distance-5 pair count.

No false pass: this proves a finite conditional selector once the distance-5
pair-count objective is supplied; it does not derive that objective, the d5
teacher trace, or a Hebbian learning law from strict geometry.
"""
from __future__ import annotations

import json
from collections import Counter
from itertools import combinations
from math import gcd
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
OUT_JSON = HERE / "bridge_strict_alpha_hebbian_distance5_pair_count_core_selector_proof_report.json"
OUT_MD = HERE / "bridge_strict_alpha_hebbian_distance5_pair_count_core_selector_proof_report.md"

N = 12
ACTIVE_COUNT = 5
TARGET_STEP = 5
AUT_UNITS = [1, 5, 7, 11]
TARGET_Q_POWER = "256/243"
TARGET_ETA = "9/5"
MINIMAL_REQUIRED_SUBGROUP = [1, 11]

Support = tuple[int, ...]


def canonical_support(nodes: set[int] | tuple[int, ...] | list[int]) -> Support:
    return tuple(sorted(node % N for node in nodes))


def all_supports() -> list[Support]:
    return [tuple(support) for support in combinations(range(N), ACTIVE_COUNT)]


def teacher_support(step: int, translate: int) -> Support:
    return canonical_support({translate + index * step for index in range(ACTIVE_COUNT)})


def teacher_orbit(step: int) -> list[Support]:
    return sorted({teacher_support(step, translate) for translate in range(N)})


def folded_distance(left: int, right: int) -> int:
    distance = (right - left) % N
    return min(distance, N - distance)


def distance_histogram(support: Support) -> tuple[int, int, int, int, int, int]:
    counts = [0] * (N // 2)
    for left, right in combinations(support, 2):
        counts[folded_distance(left, right) - 1] += 1
    return tuple(counts)  # type: ignore[return-value]


def distance5_pair_count(support: Support) -> int:
    return distance_histogram(support)[TARGET_STEP - 1]


def distance5_edges() -> list[list[int]]:
    edges = {tuple(sorted((node, (node + TARGET_STEP) % N))) for node in range(N)}
    return [list(edge) for edge in sorted(edges)]


def step_cycle_order(step: int) -> list[int]:
    order = []
    node = 0
    while node not in order:
        order.append(node)
        node = (node + step) % N
    return order


def path_supports_in_distance5_cycle() -> list[Support]:
    order = step_cycle_order(TARGET_STEP)
    path_supports = set()
    for start in range(N):
        path_supports.add(canonical_support(order[(start + offset) % N] for offset in range(ACTIVE_COUNT)))
    return sorted(path_supports)


def folded_mode(value: int) -> int:
    residue = value % N
    return min(residue, N - residue)


def distance_objective_stabilizer(distance: int) -> list[int]:
    return [unit for unit in AUT_UNITS if folded_mode(unit * distance) == distance]


def main() -> None:
    supports = all_supports()
    teacher = teacher_orbit(TARGET_STEP)
    teacher_set = set(teacher)
    path_supports = path_supports_in_distance5_cycle()
    path_set = set(path_supports)
    count_by_support = {support: distance5_pair_count(support) for support in supports}
    max_count = max(count_by_support.values())
    maximizers = sorted(support for support, count in count_by_support.items() if count == max_count)
    count_distribution = Counter(count_by_support.values())
    closest_competitors = sorted(support for support, count in count_by_support.items() if count == max_count - 1)
    closest_histogram_distribution = Counter(distance_histogram(support) for support in closest_competitors)
    cycle_order = step_cycle_order(TARGET_STEP)
    edges = distance5_edges()

    report: dict[str, Any] = {
        "result_kind": "SCRATCH_STRICT_ALPHA_HEBBIAN_DISTANCE5_PAIR_COUNT_CORE_SELECTOR_PROOF_PROBE__NOT_A_THEOREM",
        "status": "distance5-pair-count-uniquely-selects-d5-conditionally-not-origin-theorem",
        "finite_model": {
            "ring": "Z_12",
            "active_count": ACTIVE_COUNT,
            "support_count": len(supports),
            "target_distance": TARGET_STEP,
            "teacher_orbit_size": len(teacher),
            "automorphism_units": AUT_UNITS,
            "minimal_required_subgroup_from_previous_probe": MINIMAL_REQUIRED_SUBGROUP,
            "target_q_power": TARGET_Q_POWER,
            "target_eta": TARGET_ETA,
        },
        "distance5_graph_certificate": {
            "edge_count": len(edges),
            "edges": edges,
            "gcd_step_ring": gcd(TARGET_STEP, N),
            "is_single_cycle": gcd(TARGET_STEP, N) == 1 and len(edges) == N,
            "cycle_order_by_adding_5": cycle_order,
            "path_support_count": len(path_supports),
            "path_supports_equal_teacher_orbit": path_set == teacher_set,
            "graph_bound": "A 5-vertex induced subgraph of a cycle has at most 4 cycle edges; equality requires a connected 5-vertex path.",
        },
        "pair_count_selector_certificate": {
            "max_distance5_pair_count": max_count,
            "maximizer_count": len(maximizers),
            "maximizers_equal_d5_teacher_orbit": set(maximizers) == teacher_set,
            "maximizers_equal_distance5_cycle_paths": set(maximizers) == path_set,
            "count_distribution": {str(count): value for count, value in sorted(count_distribution.items())},
            "closest_nonmax_count": max_count - 1,
            "closest_nonmax_support_count": len(closest_competitors),
            "closest_nonmax_examples": [list(support) for support in closest_competitors[:8]],
            "closest_nonmax_histogram_distribution": {
                str(list(histogram)): count for histogram, count in sorted(closest_histogram_distribution.items())
            },
            "objective_unit_stabilizer": distance_objective_stabilizer(TARGET_STEP),
            "objective_stabilizer_equals_required_subgroup": distance_objective_stabilizer(TARGET_STEP) == MINIMAL_REQUIRED_SUBGROUP,
        },
        "relation_to_previous_score": {
            "previous_score_formula": "C(h)=3*h_2 + 2*h_3 + h_4 + 4*h_5",
            "core_objective": "h_5 only",
            "core_score_coefficients_d1_to_d6": [0, 0, 0, 0, 1, 0],
            "finite_gain": "The h_5 term alone already uniquely selects the d5 orbit; the larger Hebbian integer score is sufficient but not minimal as a selector objective.",
        },
        "candidate_source_interpretation": {
            "finite_gain": "The d5 selector can be reduced to the graph fact that five vertices in the distance-5 cycle maximize induced distance-5 edges exactly on 5-paths.",
            "conditional_positive_result": "If the nadsoliton self-record supplies distance-5 pair-count maximization, d5 is selected without the extra Hebbian score terms.",
            "honest_limit": "This probe does not derive why distance 5, pair-count maximization, or the distance-5 cycle should be the strict internal source.",
            "relation_to_previous_probe": "The combinatorial-score proof used C(h); this probe shows the h_5 component alone is already a finite selector.",
        },
        "ontology_guardrail": {
            "allowed_reading": "The nadsoliton itself may later be shown to carry a finite distance-5 pair-count self-record.",
            "forbidden_reading": "No separate informational layer underneath the nadsoliton is introduced to choose distance 5 or the pair-count objective.",
            "preferred_order_preserved": "nadsoliton -> light -> matter -> emergent observer",
        },
        "hard_limits": [
            "No identity K_legacy_ont == K_strict_gate is asserted or used.",
            "No legacy physical-role claim is transferred onto K_strict_gate.",
            "No theorem derives the d5 teacher/self-record trace from strict nadsoliton geometry.",
            "No theorem derives a Hebbian learning law as strict-core dynamics.",
            "No theorem derives the distance-5 pair-count objective as a strict source.",
            "The h_5 core selector is conditional on supplying distance 5 as the relevant pair channel.",
            "This is a conditional graph/combinatorial selector proof, not strict-core selector closure.",
            "No endpoint, arrow orientation, ledger selector, positive lambda action, cycle metric source, or fifth-mode source theorem is claimed.",
            "No QW-2191 discharge and no strict-core selector closure are claimed.",
            "No ToE closure is claimed.",
        ],
        "next_honest_step": "Try to derive the distance-5 pair channel or pair-count maximization internally; otherwise keep h_5 as a conditional selector objective.",
    }

    OUT_JSON.write_text(json.dumps(report, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    OUT_MD.write_text(
        "# Scratch strict-alpha Hebbian distance-5 pair-count core selector proof probe\n\n"
        "Status: h_5 pair count uniquely selects d5 conditionally; not an origin theorem.\n\n"
        f"- Supports scanned: `{len(supports)}`.\n"
        f"- Distance-5 graph is single cycle: `{report['distance5_graph_certificate']['is_single_cycle']}` with order `{cycle_order}`.\n"
        f"- Max h_5 pair count: `{max_count}`; maximizer count `{len(maximizers)}`; equals d5 teacher orbit `{set(maximizers) == teacher_set}`.\n"
        f"- h_5 count distribution: `{report['pair_count_selector_certificate']['count_distribution']}`.\n"
        f"- Closest nonmax h_5 count/supports: `{max_count - 1}` / `{len(closest_competitors)}`.\n"
        f"- Objective unit stabilizer: `{distance_objective_stabilizer(TARGET_STEP)}`.\n"
        f"- Target replay: `q^5={TARGET_Q_POWER}`, eta `{TARGET_ETA}`.\n"
        "- Honest read: h_5 alone is a finite selector, but distance-5 and pair-count maximization remain supplied inputs.\n"
        "- No false pass: no strict d5-origin theorem, no Hebbian-law theorem, no QW-2191 discharge, no ToE closure.\n",
        encoding="utf-8",
    )
    print(OUT_JSON)
    print(OUT_MD)


if __name__ == "__main__":
    main()
