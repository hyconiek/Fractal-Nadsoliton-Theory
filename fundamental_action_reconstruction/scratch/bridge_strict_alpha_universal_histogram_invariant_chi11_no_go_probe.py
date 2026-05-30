#!/usr/bin/env python3
"""Scratch probe: universal histogram-invariant no-go for sourcing chi_11.

The bounded linear and degree<=2 audits showed that low-complexity full-Aut
histogram expressions cannot export the missing chi_11 bit.  This probe removes
the degree bound for the histogram route.  It enumerates all C(12,5)=792 supports,
all distance-histogram classes, and the full-Aut shell action on histograms.

Key fact: the A5/d5 histogram is exactly the d1<->d5 swap of the A1/contiguous
histogram.  Therefore any function of the distance histogram that is invariant
under full Aut(Z_12) must assign the same value to A1 and A5.  This is not just a
linear or quadratic no-go; it is a universal no-go for full-Aut-invariant
histogram-only sources.

No false pass: anti-invariant histogram functions can still transform like
chi_11, but only by choosing an orientation of the d1/d5 shell pair.  This probe
does not exclude non-histogram strict sources, does not derive chi_11, does not
discharge QW-2191, and does not close ToE.
"""
from __future__ import annotations

import json
from collections import Counter
from itertools import combinations
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
OUT_JSON = HERE / "bridge_strict_alpha_universal_histogram_invariant_chi11_no_go_report.json"
OUT_MD = HERE / "bridge_strict_alpha_universal_histogram_invariant_chi11_no_go_report.md"

N = 12
ACTIVE_COUNT = 5
UNITS = [1, 5, 7, 11]
TARGET_Q_POWER = "256/243"
TARGET_ETA = "9/5"
A1_HISTOGRAM = (4, 3, 2, 1, 0, 0)
A5_HISTOGRAM = (0, 3, 2, 1, 4, 0)


def folded(value: int) -> int:
    residue = value % N
    return min(residue, N - residue)


def distance_histogram(support: tuple[int, ...]) -> tuple[int, int, int, int, int, int]:
    counts = [0] * (N // 2)
    for left, right in combinations(support, 2):
        counts[folded(right - left) - 1] += 1
    return tuple(counts)  # type: ignore[return-value]


def all_supports() -> list[tuple[int, ...]]:
    return [tuple(support) for support in combinations(range(N), ACTIVE_COUNT)]


def swap_d1_d5(histogram: tuple[int, int, int, int, int, int]) -> tuple[int, int, int, int, int, int]:
    return (histogram[4], histogram[1], histogram[2], histogram[3], histogram[0], histogram[5])


def histogram_orbits(histograms: list[tuple[int, int, int, int, int, int]]) -> list[list[tuple[int, int, int, int, int, int]]]:
    remaining = set(histograms)
    orbits = []
    while remaining:
        seed = min(remaining)
        orbit = sorted(({seed, swap_d1_d5(seed)}) & remaining)
        orbits.append(orbit)
        remaining -= set(orbit)
    return sorted(orbits, key=lambda orbit: orbit[0])


def orbit_rows(orbits: list[list[tuple[int, int, int, int, int, int]]], counts: Counter[tuple[int, int, int, int, int, int]]) -> list[dict[str, Any]]:
    rows = []
    for index, orbit in enumerate(orbits):
        rows.append(
            {
                "orbit_index": index,
                "histograms_d1_to_d6": [list(histogram) for histogram in orbit],
                "orbit_size": len(orbit),
                "support_count_total": sum(counts[histogram] for histogram in orbit),
                "contains_A1_contiguous_histogram": A1_HISTOGRAM in orbit,
                "contains_A5_d5_histogram": A5_HISTOGRAM in orbit,
                "is_A1_A5_pair_orbit": set(orbit) == {A1_HISTOGRAM, A5_HISTOGRAM},
            }
        )
    return rows


def exact_proof_certificate(histogram_class_count: int, orbit_count: int) -> dict[str, str]:
    return {
        "finite_domain": "All C(12,5)=792 supports are enumerated and collapsed to distance-histogram classes.",
        "histogram_swap_fact": "A1 has histogram [4,3,2,1,0,0] and A5 has [0,3,2,1,4,0], exactly the d1<->d5 swap.",
        "universal_invariant_no_go": "For any histogram-only source F invariant under full Aut, F(h)=F(swap_d1_d5(h)); hence F(A1)=F(A5) and F cannot export chi_11.",
        "classifier_count": f"The {histogram_class_count} histogram classes form {orbit_count} full-Aut histogram orbits; invariant Boolean classifiers are constant on these orbits, so singleton A5-vs-A1 classification count is 0.",
        "anti_invariant_boundary": "A histogram expression can flip between A1 and A5 only if it is anti-invariant under d1<->d5, which imports the missing shell-label orientation.",
        "scope_limit": "This no-go is universal for histogram-only full-Aut invariant sources, not for all possible non-histogram strict nadsoliton sources.",
    }


def build_payload() -> dict[str, Any]:
    supports = all_supports()
    counts = Counter(distance_histogram(support) for support in supports)
    histograms = sorted(counts)
    orbits = histogram_orbits(histograms)
    a1_a5_orbit = next(orbit for orbit in orbits if set(orbit) == {A1_HISTOGRAM, A5_HISTOGRAM})
    invariant_boolean_classifier_count_power = len(orbits)
    invariant_boolean_classifiers_total = 2 ** invariant_boolean_classifier_count_power
    both_or_neither_a1_a5_classifier_count = invariant_boolean_classifiers_total
    singleton_a5_not_a1_classifier_count = 0
    return {
        "result_kind": "SCRATCH_STRICT_ALPHA_UNIVERSAL_HISTOGRAM_INVARIANT_CHI11_NO_GO_PROBE__NO_GO_NOT_A_THEOREM",
        "status": "full-aut-invariant-histogram-only-sources-cannot-distinguish-a1-from-a5-or-export-chi11",
        "finite_model": {
            "ring": "Z_12",
            "active_count": ACTIVE_COUNT,
            "automorphism_units": UNITS,
            "support_count": len(supports),
            "histogram_class_count": len(histograms),
            "full_Aut_histogram_orbit_count": len(orbits),
            "target_q_power": TARGET_Q_POWER,
            "target_eta": TARGET_ETA,
        },
        "a1_a5_histogram_orbit_certificate": {
            "A1_contiguous_histogram": list(A1_HISTOGRAM),
            "A5_d5_histogram": list(A5_HISTOGRAM),
            "swap_d1_d5_A1_histogram": list(swap_d1_d5(A1_HISTOGRAM)),
            "same_full_Aut_histogram_orbit": set(a1_a5_orbit) == {A1_HISTOGRAM, A5_HISTOGRAM},
            "A1_support_count": counts[A1_HISTOGRAM],
            "A5_support_count": counts[A5_HISTOGRAM],
            "orbit_support_count": sum(counts[histogram] for histogram in a1_a5_orbit),
        },
        "histogram_orbit_summary": {
            "histogram_class_count": len(histograms),
            "full_Aut_histogram_orbit_count": len(orbits),
            "fixed_histogram_orbit_count": sum(1 for orbit in orbits if len(orbit) == 1),
            "two_histogram_orbit_count": sum(1 for orbit in orbits if len(orbit) == 2),
            "invariant_boolean_classifier_count_power": invariant_boolean_classifier_count_power,
            "invariant_boolean_classifier_total": invariant_boolean_classifiers_total,
            "singleton_A5_not_A1_invariant_classifier_count": singleton_a5_not_a1_classifier_count,
            "both_or_neither_A1_A5_invariant_classifier_count": both_or_neither_a1_a5_classifier_count,
        },
        "histogram_orbit_rows": orbit_rows(orbits, counts),
        "exact_proof_certificate": exact_proof_certificate(len(histograms), len(orbits)),
        "interpretation": {
            "honest_positive": "The no-go now covers every possible full-Aut invariant function of the distance histogram, not only bounded linear or quadratic formulas.",
            "honest_negative": "The result still leaves non-histogram strict sources open; it only closes the histogram-only full-Aut invariant route.",
            "relation_to_previous_probe": "This subsumes the linear and quadratic histogram invariant no-go audits at the level of arbitrary invariant histogram functions.",
        },
        "ontology_guardrail": {
            "allowed_reading": "The nadsoliton itself remains the only admissible place where a future non-histogram strict chi_11 source could be derived.",
            "forbidden_reading": "No separate informational layer underneath the nadsoliton is introduced to supply chi_11.",
        },
        "hard_limits": [
            "No identity K_legacy_ont == K_strict_gate is used or claimed.",
            "No legacy role transfer to K_strict_gate, alpha_geo, beta_tors, or D_f is made.",
            "No theorem derives the chi_11-kernel, shell-label d1 vs d5, unit-axis bit, exact-cover clauses, or cardinality 5 from strict nadsoliton geometry.",
            "The result is a universal histogram-only full-Aut invariant no-go, not an exhaustive strict-source theorem.",
            "No QW-2191 discharge.",
            "No ToE closure.",
        ],
    }


def write_markdown(payload: dict[str, Any]) -> str:
    summary = payload["histogram_orbit_summary"]
    cert = payload["a1_a5_histogram_orbit_certificate"]
    lines = [
        "# Strict alpha universal histogram-invariant chi_11 no-go",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Histogram orbit summary",
        "",
        f"- Histogram classes: `{summary['histogram_class_count']}`",
        f"- Full-Aut histogram orbits: `{summary['full_Aut_histogram_orbit_count']}`",
        f"- Fixed histogram orbits: `{summary['fixed_histogram_orbit_count']}`",
        f"- Two-histogram orbits: `{summary['two_histogram_orbit_count']}`",
        f"- Invariant Boolean classifier total: `{summary['invariant_boolean_classifier_total']}`",
        f"- Singleton A5-not-A1 invariant classifiers: `{summary['singleton_A5_not_A1_invariant_classifier_count']}`",
        "",
        "## A1/A5 orbit certificate",
        "",
        f"- A1 histogram: `{cert['A1_contiguous_histogram']}`",
        f"- A5 histogram: `{cert['A5_d5_histogram']}`",
        f"- swap_d1_d5(A1): `{cert['swap_d1_d5_A1_histogram']}`",
        f"- Same full-Aut histogram orbit: `{cert['same_full_Aut_histogram_orbit']}`",
        f"- Orbit support count: `{cert['orbit_support_count']}`",
        "",
        "## Proof certificate",
        "",
    ]
    for key, value in payload["exact_proof_certificate"].items():
        lines.append(f"- `{key}`: {value}")
    lines.extend(["", "## Hard limits", ""])
    lines.extend(f"- {item}" for item in payload["hard_limits"])
    lines.append("")
    return "\n".join(lines)


def main() -> None:
    payload = build_payload()
    OUT_JSON.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    OUT_MD.write_text(write_markdown(payload), encoding="utf-8")
    print(json.dumps(payload, indent=2, sort_keys=True))
    print(OUT_JSON)
    print(OUT_MD)


if __name__ == "__main__":
    main()
