#!/usr/bin/env python3
"""Scratch probe: exact-cover certificate for the cyclic exclusion selector.

The previous exclusion-energy audit found that E_excl=d1+d6 uniquely selects the
A5/d5 dihedral orbit among all 5-subsets of Z_12.  This probe recasts that fact
as a more proof-like exact-cover / SAT certificate:

    variables x_0..x_11 in {0,1}
    sum_i x_i = 5
    x_i + x_{i+1} <= 1  for all cyclic nearest-neighbour edges
    x_i + x_{i+6} <= 1  for all antipodal pairs

The satisfying assignments are exactly the E_excl=0 supports.  The certificate
then checks that there are 12 assignments, one dihedral orbit, and one cyclic gap
necklace [2,2,3,2,3], whose folded distance histogram is [0,3,2,1,4,0].

No false pass: this is still a conditional finite selector certificate.  It does
not derive the exact-cover constraints from strict nadsoliton geometry and does
not discharge QW-2191 or close ToE.
"""
from __future__ import annotations

import json
from collections import Counter
from itertools import combinations
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
OUT_JSON = HERE / "bridge_strict_alpha_cyclic_exclusion_exact_cover_certificate_report.json"
OUT_MD = HERE / "bridge_strict_alpha_cyclic_exclusion_exact_cover_certificate_report.md"

N = 12
ACTIVE_COUNT = 5
TARGET_Q_POWER = "256/243"
TARGET_ETA = "9/5"
A5_SUPPORT = (0, 5, 10, 3, 8)
A1_SUPPORT = tuple(range(5))
A5_HISTOGRAM = [0, 3, 2, 1, 4, 0]


def all_supports() -> list[tuple[int, ...]]:
    return [tuple(support) for support in combinations(range(N), ACTIVE_COUNT)]


def folded(value: int) -> int:
    residue = value % N
    return min(residue, N - residue)


def distance_histogram(support: tuple[int, ...]) -> tuple[int, int, int, int, int, int]:
    counts = [0] * (N // 2)
    for left, right in combinations(support, 2):
        counts[folded(right - left) - 1] += 1
    return tuple(counts)  # type: ignore[return-value]


def canonical_support(support: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(sorted(value % N for value in support))


def nearest_edges() -> list[tuple[int, int]]:
    return [(index, (index + 1) % N) for index in range(N)]


def antipodal_edges() -> list[tuple[int, int]]:
    return [(index, index + N // 2) for index in range(N // 2)]


def violates_pair(support: tuple[int, ...], pair: tuple[int, int]) -> bool:
    support_set = set(support)
    return pair[0] in support_set and pair[1] in support_set


def satisfies_exact_cover_constraints(support: tuple[int, ...]) -> bool:
    return (
        len(support) == ACTIVE_COUNT
        and not any(violates_pair(support, pair) for pair in nearest_edges())
        and not any(violates_pair(support, pair) for pair in antipodal_edges())
    )


def dihedral_orbit(support: tuple[int, ...]) -> set[tuple[int, ...]]:
    support_set = set(support)
    orbit = set()
    for shift in range(N):
        orbit.add(tuple(sorted((value + shift) % N for value in support_set)))
        orbit.add(tuple(sorted((-value + shift) % N for value in support_set)))
    return orbit


def orbit_representatives(supports: list[tuple[int, ...]]) -> list[dict[str, Any]]:
    remaining = set(supports)
    rows = []
    while remaining:
        seed = min(remaining)
        orbit = dihedral_orbit(seed)
        members = sorted(remaining & orbit)
        remaining -= orbit
        representative = members[0]
        rows.append(
            {
                "representative": list(representative),
                "orbit_size": len(members),
                "distance_histogram_d1_to_d6": list(distance_histogram(representative)),
                "is_A5_d5_orbit": canonical_support(A5_SUPPORT) in orbit,
                "is_A1_contiguous_orbit": canonical_support(A1_SUPPORT) in orbit,
            }
        )
    return sorted(rows, key=lambda row: row["representative"])


def cyclic_gaps(support: tuple[int, ...]) -> tuple[int, ...]:
    ordered = sorted(support)
    return tuple((ordered[(index + 1) % len(ordered)] - ordered[index]) % N or N for index in range(len(ordered)))


def canonical_gap_necklace(gaps: tuple[int, ...]) -> tuple[int, ...]:
    rotations = [gaps[index:] + gaps[:index] for index in range(len(gaps))]
    reversed_gaps = tuple(reversed(gaps))
    rotations.extend(reversed_gaps[index:] + reversed_gaps[:index] for index in range(len(reversed_gaps)))
    return min(rotations)


def gap_necklace_audit(supports: list[tuple[int, ...]]) -> list[dict[str, Any]]:
    by_necklace: dict[tuple[int, ...], list[tuple[int, ...]]] = {}
    for support in supports:
        necklace = canonical_gap_necklace(cyclic_gaps(support))
        by_necklace.setdefault(necklace, []).append(support)

    rows = []
    for necklace, members in sorted(by_necklace.items()):
        histograms = sorted({distance_histogram(member) for member in members})
        rows.append(
            {
                "canonical_gap_necklace": list(necklace),
                "support_count": len(members),
                "histograms_d1_to_d6": [list(histogram) for histogram in histograms],
                "has_antipodal_contact": any(histogram[5] > 0 for histogram in histograms),
                "is_A5_d5_necklace": any(list(histogram) == A5_HISTOGRAM for histogram in histograms),
            }
        )
    return rows


def exact_cover_solution_audit(supports: list[tuple[int, ...]]) -> dict[str, Any]:
    solutions = [support for support in supports if satisfies_exact_cover_constraints(support)]
    hist_counter = Counter(distance_histogram(solution) for solution in solutions)
    return {
        "variable_count": N,
        "cardinality_constraint": f"sum x_i = {ACTIVE_COUNT}",
        "nearest_edge_constraint_count": len(nearest_edges()),
        "antipodal_edge_constraint_count": len(antipodal_edges()),
        "forbidden_pair_count": len(nearest_edges()) + len(antipodal_edges()),
        "solution_count": len(solutions),
        "solution_histogram_rows": [
            {"distance_histogram_d1_to_d6": list(histogram), "support_count": count}
            for histogram, count in sorted(hist_counter.items())
        ],
        "solution_orbit_count": len(orbit_representatives(solutions)),
        "solution_orbits": orbit_representatives(solutions),
        "solution_gap_necklaces": gap_necklace_audit(solutions),
        "selects_A5_d5_orbit": len(orbit_representatives(solutions)) == 1 and orbit_representatives(solutions)[0]["is_A5_d5_orbit"],
    }


def d1_zero_gap_classification(supports: list[tuple[int, ...]]) -> dict[str, Any]:
    d1_zero = [support for support in supports if distance_histogram(support)[0] == 0]
    rows = gap_necklace_audit(d1_zero)
    return {
        "d1_zero_support_count": len(d1_zero),
        "d1_zero_gap_necklace_count": len(rows),
        "d1_zero_gap_necklaces": rows,
        "exact_cover_surviving_necklaces": [row for row in rows if not row["has_antipodal_contact"]],
    }


def exact_proof_certificate() -> dict[str, str]:
    return {
        "boolean_formulation": "Use x_0..x_11 in {0,1}, sum x_i=5, forbid nearest-neighbour pairs and antipodal pairs.",
        "enumeration_domain": "The certificate checks all C(12,5)=792 cardinality-compatible assignments exactly.",
        "solution_count": "Exactly 12 assignments satisfy the constraints.",
        "orbit_certificate": "The 12 assignments form one dihedral orbit, containing the A5/d5 support class and not the A1 contiguous class.",
        "gap_certificate": "The surviving cyclic gap necklace is [2,2,3,2,3], with folded histogram [0,3,2,1,4,0].",
        "d1_only_boundary": "Without antipodal constraints, d1=0 gives three gap necklaces; the antipodal ban leaves only the A5/d5 necklace.",
        "missing_source": "The exact-cover constraints are not derived here from strict nadsoliton geometry; they are a conditional selector premise.",
    }


def main() -> None:
    supports = all_supports()
    all_histograms = {distance_histogram(support) for support in supports}
    solution_audit = exact_cover_solution_audit(supports)
    d1_gap_audit = d1_zero_gap_classification(supports)

    report: dict[str, Any] = {
        "result_kind": "SCRATCH_STRICT_ALPHA_CYCLIC_EXCLUSION_EXACT_COVER_CERTIFICATE_PROBE__CONDITIONAL_NOT_A_THEOREM",
        "status": "exact-cover-nearest-plus-antipodal-constraints-select-A5-d5-orbit-conditionally-source-not-derived",
        "finite_model": {
            "ring": "Z_12",
            "active_count": ACTIVE_COUNT,
            "support_count": len(supports),
            "histogram_class_count": len(all_histograms),
            "target_q_power": TARGET_Q_POWER,
            "target_eta": TARGET_ETA,
        },
        "constraint_system": {
            "variables": [f"x_{index}" for index in range(N)],
            "cardinality": f"sum x_i = {ACTIVE_COUNT}",
            "nearest_forbidden_pairs": [list(pair) for pair in nearest_edges()],
            "antipodal_forbidden_pairs": [list(pair) for pair in antipodal_edges()],
        },
        "exact_cover_solution_audit": solution_audit,
        "d1_zero_gap_classification": d1_gap_audit,
        "exact_proof_certificate": exact_proof_certificate(),
        "interpretation": {
            "direct_result": "The exact-cover constraints sum=5, no adjacent pair, and no antipodal pair select exactly the A5/d5 dihedral orbit.",
            "proof_style_upgrade": "This is a Boolean constraint certificate for the same d1+d6 exclusion selector, independent of the previous energy-minimization presentation.",
            "d1_boundary": "The nearest-neighbour constraint alone gives three gap necklaces; the antipodal clauses are the extra finite ingredient that remove the two non-A5 necklaces.",
            "honest_limit": "The constraints are not proven strict-core outputs here; the selector remains conditional/non-strict unless a strict source exports them.",
        },
        "ontology_guardrail": {
            "allowed_reading": "The exact-cover constraints may be investigated as candidate structures of the nadsoliton itself.",
            "forbidden_reading": "No separate informational layer underneath the nadsoliton is introduced to provide the exact-cover constraints or unit-axis bit.",
            "preferred_order_preserved": "nadsoliton -> light -> matter -> emergent observer",
        },
        "hard_limits": [
            "No identity K_legacy_ont == K_strict_gate is asserted or used.",
            "No legacy physical-role or matter-generation claim is transferred onto K_strict_gate.",
            "No theorem derives exact-cover constraints, cyclic exclusion energy, chi_11, or the unit-axis bit from strict nadsoliton geometry.",
            "The exact-cover selector is a conditional finite certificate, not a strict-core theorem.",
            "No endpoint, arrow orientation, ledger selector, positive lambda action, cycle metric source, anti-Nyquist source, fifth-mode source, future-probability source, future-value source, future-path source, matter-bit source, existence-bit source, recursive-self-information source, character-source, meta-character source, cyclic-adjacency source theorem, variational tie-break source theorem, exclusion-energy source theorem, or legacy-strict bridge theorem is claimed.",
            "No QW-2191 discharge and no strict-core selector closure are claimed.",
            "No ToE closure is claimed.",
        ],
        "next_honest_step": "Audit existing strict-side nad12/sigma artifacts for a source of the exact-cover clauses; otherwise keep the clauses as an explicit non-strict selector premise.",
    }

    OUT_JSON.write_text(json.dumps(report, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    OUT_MD.write_text(
        "# Scratch strict-alpha cyclic-exclusion exact-cover certificate probe\n\n"
        "Status: exact-cover nearest-plus-antipodal constraints select A5/d5 conditionally; source premise not derived here.\n\n"
        f"- Supports scanned: `{len(supports)}`; histogram classes: `{len(all_histograms)}`.\n"
        f"- Boolean variables: `{solution_audit['variable_count']}`.\n"
        f"- Nearest forbidden pairs: `{solution_audit['nearest_edge_constraint_count']}`; antipodal forbidden pairs: `{solution_audit['antipodal_edge_constraint_count']}`.\n"
        f"- Exact-cover solution count: `{solution_audit['solution_count']}`.\n"
        f"- Solution orbit count: `{solution_audit['solution_orbit_count']}`.\n"
        f"- Selects A5/d5 orbit: `{solution_audit['selects_A5_d5_orbit']}`.\n"
        f"- d1=0 gap necklace count before antipodal clauses: `{d1_gap_audit['d1_zero_gap_necklace_count']}`.\n"
        f"- Surviving gap necklaces: `{d1_gap_audit['exact_cover_surviving_necklaces']}`.\n"
        f"- Target replay kept conditional: `q^5={TARGET_Q_POWER}`, eta `{TARGET_ETA}`.\n"
        "- No false pass: no strict exact-cover source theorem, no QW-2191 discharge, no ToE closure.\n",
        encoding="utf-8",
    )
    print(OUT_JSON)
    print(OUT_MD)


if __name__ == "__main__":
    main()
