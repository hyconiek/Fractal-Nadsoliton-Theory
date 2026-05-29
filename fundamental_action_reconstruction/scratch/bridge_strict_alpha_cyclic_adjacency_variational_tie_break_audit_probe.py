#!/usr/bin/env python3
"""Scratch probe: cyclic adjacency needs a variational tie-break to select d5.

The previous cyclic-adjacency audit proved a conditional source fact: the
nearest-neighbour shell {+1,-1} has unit stabilizer {1,11}, the kernel of the
required chi_11 character.  This probe checks the next finite question:

    If cyclic adjacency is admitted, does the nearest-neighbour datum alone
    select the d5/A5 support orbit among all 5-subsets of Z_12?

No.  Enumerating all C(12,5)=792 supports shows that minimizing nearest-neighbour
contacts d1 leaves 36 supports, grouped into 3 dihedral orbits / 3 distance
histograms.  The d5/A5 orbit is one of them, not uniquely selected.

A conditional positive refinement exists: among the d1-minimizers, maximizing
the fifth-shell count d5 uniquely selects the A5/d5 dihedral orbit.  But this
uses an additional variational polarity/tie-break rule (min d1, then max d5),
so it is not a strict-core theorem and it does not discharge QW-2191.
"""
from __future__ import annotations

import json
from collections import Counter
from itertools import combinations
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
OUT_JSON = HERE / "bridge_strict_alpha_cyclic_adjacency_variational_tie_break_audit_report.json"
OUT_MD = HERE / "bridge_strict_alpha_cyclic_adjacency_variational_tie_break_audit_report.md"

N = 12
ACTIVE_COUNT = 5
TARGET_Q_POWER = "256/243"
TARGET_ETA = "9/5"
A1_SUPPORT = tuple(range(5))
A5_SUPPORT = (0, 5, 10, 3, 8)
A1 = "A1_k1_contiguous"
A5 = "A5_k5_d5"


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
        hist = distance_histogram(representative)
        rows.append(
            {
                "representative": list(representative),
                "orbit_size": len(members),
                "distance_histogram_d1_to_d6": list(hist),
                "d1_nearest_neighbor_count": hist[0],
                "d5_fifth_shell_count": hist[4],
                "is_A5_d5_orbit": canonical_support(A5_SUPPORT) in orbit,
                "is_A1_contiguous_orbit": canonical_support(A1_SUPPORT) in orbit,
            }
        )
    return sorted(rows, key=lambda row: (row["d1_nearest_neighbor_count"], -row["d5_fifth_shell_count"], row["representative"]))


def d1_distribution(supports: list[tuple[int, ...]]) -> list[dict[str, Any]]:
    rows = []
    for d1 in sorted({distance_histogram(support)[0] for support in supports}):
        bucket = [support for support in supports if distance_histogram(support)[0] == d1]
        rows.append(
            {
                "d1_nearest_neighbor_count": d1,
                "support_count": len(bucket),
                "histogram_class_count": len({distance_histogram(support) for support in bucket}),
                "dihedral_orbit_count": len(orbit_representatives(bucket)),
            }
        )
    return rows


def adjacency_minimizer_audit(supports: list[tuple[int, ...]]) -> dict[str, Any]:
    min_d1 = min(distance_histogram(support)[0] for support in supports)
    minimizers = [support for support in supports if distance_histogram(support)[0] == min_d1]
    hist_counter = Counter(distance_histogram(support) for support in minimizers)
    orbit_rows = orbit_representatives(minimizers)
    return {
        "minimum_d1": min_d1,
        "minimum_d1_support_count": len(minimizers),
        "minimum_d1_histogram_class_count": len(hist_counter),
        "minimum_d1_dihedral_orbit_count": len(orbit_rows),
        "minimum_d1_histogram_rows": [
            {
                "distance_histogram_d1_to_d6": list(hist),
                "support_count": count,
                "is_A5_d5_histogram": list(hist) == list(distance_histogram(A5_SUPPORT)),
            }
            for hist, count in sorted(hist_counter.items())
        ],
        "minimum_d1_dihedral_orbits": orbit_rows,
        "nearest_neighbor_minimization_uniquely_selects_A5": len(orbit_rows) == 1 and orbit_rows[0]["is_A5_d5_orbit"],
    }


def lexicographic_tie_break_audit(supports: list[tuple[int, ...]]) -> dict[str, Any]:
    scored = [(distance_histogram(support)[0], -distance_histogram(support)[4], support) for support in supports]
    best_d1, best_neg_d5, _ = min(scored)
    best_d5 = -best_neg_d5
    winners = [support for d1, neg_d5, support in scored if d1 == best_d1 and -neg_d5 == best_d5]
    orbit_rows = orbit_representatives(winners)
    return {
        "rule": "lexicographic_minimize_d1_then_maximize_d5",
        "best_d1": best_d1,
        "best_d5": best_d5,
        "winner_support_count": len(winners),
        "winner_dihedral_orbit_count": len(orbit_rows),
        "winner_histograms": [list(hist) for hist in sorted({distance_histogram(support) for support in winners})],
        "winner_orbits": orbit_rows,
        "selects_A5_d5_orbit": len(orbit_rows) == 1 and orbit_rows[0]["is_A5_d5_orbit"],
        "honest_status": "conditional: requires both adjacency/anti-adjacency polarity and a fifth-shell maximization tie-break",
    }


def pairwise_A1_A5_audit() -> dict[str, Any]:
    a1_hist = distance_histogram(A1_SUPPORT)
    a5_hist = distance_histogram(A5_SUPPORT)
    return {
        "A1_support": list(A1_SUPPORT),
        "A1_histogram_d1_to_d6": list(a1_hist),
        "A5_support": list(canonical_support(A5_SUPPORT)),
        "A5_histogram_d1_to_d6": list(a5_hist),
        "minimize_d1_prefers_A5_over_A1_pairwise": a5_hist[0] < a1_hist[0],
        "global_warning": "Pairwise A1-vs-A5 preference is weaker than global support-orbit uniqueness; d1 minimization has three global dihedral-orbit minimizers.",
    }


def exact_proof_certificate() -> dict[str, str]:
    return {
        "enumeration_domain": "All 5-subsets of Z_12 are enumerated exactly: C(12,5)=792.",
        "adjacency_energy": "The nearest-neighbour adjacency scalar is d1, the first entry of the folded distance histogram.",
        "no_go": "The minimum d1 value is 0, but its minimizer set has 36 supports, 3 histograms, and 3 dihedral orbits; therefore d1 alone does not uniquely select A5/d5.",
        "conditional_positive": "The lexicographic rule min d1 then max d5 has 12 support winners forming one dihedral orbit with histogram [0,3,2,1,4,0], the A5/d5 orbit.",
        "missing_source": "The polarity/tie-break rule is extra variational structure and is not derived here from strict nadsoliton geometry.",
    }


def main() -> None:
    supports = all_supports()
    all_histograms = {distance_histogram(support) for support in supports}
    pairwise = pairwise_A1_A5_audit()
    minimizer = adjacency_minimizer_audit(supports)
    tie_break = lexicographic_tie_break_audit(supports)

    report: dict[str, Any] = {
        "result_kind": "SCRATCH_STRICT_ALPHA_CYCLIC_ADJACENCY_VARIATIONAL_TIE_BREAK_AUDIT_PROBE__CONDITIONAL_NOT_A_THEOREM",
        "status": "nearest-neighbor-adjacency-alone-not-global-d5-selector-lexicographic-d1-d5-rule-conditionally-selects-A5",
        "finite_model": {
            "ring": "Z_12",
            "active_count": ACTIVE_COUNT,
            "support_count": len(supports),
            "histogram_class_count": len(all_histograms),
            "target_q_power": TARGET_Q_POWER,
            "target_eta": TARGET_ETA,
        },
        "pairwise_A1_A5_audit": pairwise,
        "d1_distribution": d1_distribution(supports),
        "adjacency_minimizer_audit": minimizer,
        "lexicographic_tie_break_audit": tie_break,
        "exact_proof_certificate": exact_proof_certificate(),
        "interpretation": {
            "direct_result": "Cyclic adjacency can provide a chi_11-compatible label, but nearest-neighbour minimization alone does not globally select the A5/d5 support orbit.",
            "pairwise_result": "Against the contiguous A1 support alone, minimizing d1 prefers A5 because A5 has d1=0 while A1 has d1=4.",
            "global_negative_result": "Among all 5-subsets of Z_12, d1=0 has three dihedral-orbit minimizers, so adjacency alone leaves a selector ambiguity.",
            "conditional_positive_result": "Adding a second rule, maximize d5 among d1-minimizers, uniquely selects the A5/d5 orbit.",
            "honest_limit": "The extra d5 polarity/tie-break is not derived from strict core here; it must remain an explicit conditional premise unless exported elsewhere.",
        },
        "ontology_guardrail": {
            "allowed_reading": "Adjacency and any future variational rule may be investigated as candidate structures of the nadsoliton itself.",
            "forbidden_reading": "No separate informational layer underneath the nadsoliton is introduced to provide the tie-break or unit-axis bit.",
            "preferred_order_preserved": "nadsoliton -> light -> matter -> emergent observer",
        },
        "hard_limits": [
            "No identity K_legacy_ont == K_strict_gate is asserted or used.",
            "No legacy physical-role or matter-generation claim is transferred onto K_strict_gate.",
            "No theorem derives cyclic adjacency, anti-adjacency polarity, d5 maximization, chi_11, or the unit-axis bit from strict nadsoliton geometry.",
            "The lexicographic d1/d5 rule is a conditional variational tie-break audit, not a strict-core theorem.",
            "No endpoint, arrow orientation, ledger selector, positive lambda action, cycle metric source, anti-Nyquist source, fifth-mode source, future-probability source, future-value source, future-path source, matter-bit source, existence-bit source, recursive-self-information source, character-source, meta-character source, cyclic-adjacency source theorem, or legacy-strict bridge theorem is claimed.",
            "No QW-2191 discharge and no strict-core selector closure are claimed.",
            "No ToE closure is claimed.",
        ],
        "next_honest_step": "Search existing strict-side nad12/sigma artifacts for an exported variational polarity that justifies min d1 and then max d5; absent that, record d1/d5 selection as an explicit non-strict selector premise.",
    }

    OUT_JSON.write_text(json.dumps(report, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    OUT_MD.write_text(
        "# Scratch strict-alpha cyclic-adjacency variational tie-break audit probe\n\n"
        "Status: nearest-neighbour adjacency alone is not a global d5 selector; d1/d5 tie-break is conditional.\n\n"
        f"- Supports scanned: `{len(supports)}`; histogram classes: `{len(all_histograms)}`.\n"
        f"- A1 histogram d1..d6: `{pairwise['A1_histogram_d1_to_d6']}`.\n"
        f"- A5 histogram d1..d6: `{pairwise['A5_histogram_d1_to_d6']}`.\n"
        f"- Minimum d1 support count: `{minimizer['minimum_d1_support_count']}`.\n"
        f"- Minimum d1 dihedral-orbit count: `{minimizer['minimum_d1_dihedral_orbit_count']}`.\n"
        f"- d1 alone uniquely selects A5: `{minimizer['nearest_neighbor_minimization_uniquely_selects_A5']}`.\n"
        f"- Lexicographic min d1 / max d5 selects A5: `{tie_break['selects_A5_d5_orbit']}`.\n"
        f"- Lexicographic winner support count: `{tie_break['winner_support_count']}`; orbit count: `{tie_break['winner_dihedral_orbit_count']}`.\n"
        f"- Target replay kept conditional: `q^5={TARGET_Q_POWER}`, eta `{TARGET_ETA}`.\n"
        "- No false pass: no strict variational-source theorem, no QW-2191 discharge, no ToE closure.\n",
        encoding="utf-8",
    )
    print(OUT_JSON)
    print(OUT_MD)


if __name__ == "__main__":
    main()
