#!/usr/bin/env python3
"""Scratch probe: D12 chi_11 generator selected by max shell imbalance.

The D12 character-module certificate found 13 independent chi_11-covariant
basis generators on the D12 quotient.  This probe asks the next finite question:
if the reduced D12 quotient and the shell-labelled d1/d5 axis are already
available, is the branch generator still one arbitrary generator among 13, or is
it selected by a small computable extremal rule?

We enumerate the same 38 D12 orbits and the residual unit-5 two-cycles.  For
each two-cycle {O,tau(O)}, compute the signed shell imbalance

    I(O) = h_d5(O) - h_d1(O),

where h is the distance histogram of a representative.  Unit 5 swaps d1 and d5,
so I(tau(O))=-I(O).  The branch cycle A1/A11 <-> A5/A7 is the unique two-cycle
with maximal |I|=4.  Therefore a conditional D12+shell-label selector can pick
the branch generator by maximal shell imbalance.

No false pass: this is deliberately conditional on the reduced D12 quotient and
the shell-labelled d1/d5 axis.  It does not derive the unit-axis bit, does not
derive shell labels from strict full-Aut data, does not discharge QW-2191, and
does not close ToE.
"""
from __future__ import annotations

import json
from collections import Counter
from itertools import combinations
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
OUT_JSON = HERE / "bridge_strict_alpha_d12_chi11_max_shell_imbalance_selector_certificate_report.json"
OUT_MD = HERE / "bridge_strict_alpha_d12_chi11_max_shell_imbalance_selector_certificate_report.md"

N = 12
ACTIVE_COUNT = 5
UNITS = [1, 5, 7, 11]
DIHEDRAL_UNITS = [1, 11]
BRANCH_MODES = [1, 5, 7, 11]
TARGET_Q_POWER = "256/243"
TARGET_ETA = "9/5"


def unit_support(mode: int) -> tuple[int, ...]:
    return tuple(sorted((mode * index) % N for index in range(ACTIVE_COUNT)))


def affine_image(support: tuple[int, ...], shift: int, unit: int) -> tuple[int, ...]:
    return tuple(sorted((shift + unit * value) % N for value in support))


def all_supports() -> list[tuple[int, ...]]:
    return [tuple(support) for support in combinations(range(N), ACTIVE_COUNT)]


def affine_orbit(support: tuple[int, ...], units: list[int]) -> frozenset[tuple[int, ...]]:
    return frozenset(affine_image(support, shift, unit) for shift in range(N) for unit in units)


def orbit_partition(supports: list[tuple[int, ...]], units: list[int]) -> list[frozenset[tuple[int, ...]]]:
    remaining = set(supports)
    orbits: list[frozenset[tuple[int, ...]]] = []
    while remaining:
        seed = min(remaining)
        orbit = affine_orbit(seed, units)
        orbits.append(orbit)
        remaining -= orbit
    return sorted(orbits, key=lambda orbit: min(orbit))


def index_by_support(orbits: list[frozenset[tuple[int, ...]]]) -> dict[tuple[int, ...], int]:
    return {support: index for index, orbit in enumerate(orbits) for support in orbit}


def folded(value: int) -> int:
    residue = value % N
    return min(residue, N - residue)


def distance_histogram(support: tuple[int, ...]) -> tuple[int, int, int, int, int, int]:
    counts = [0] * (N // 2)
    for left, right in combinations(support, 2):
        counts[folded(right - left) - 1] += 1
    return tuple(counts)  # type: ignore[return-value]


def cyclic_gaps(support: tuple[int, ...]) -> tuple[int, ...]:
    cyclic = sorted(support)
    gaps = [cyclic[index + 1] - cyclic[index] for index in range(len(cyclic) - 1)]
    gaps.append(N + cyclic[0] - cyclic[-1])
    return tuple(gaps)


def gap_necklace(support: tuple[int, ...]) -> tuple[int, ...]:
    gaps = cyclic_gaps(support)
    rotations = [gaps[index:] + gaps[:index] for index in range(len(gaps))]
    reversals = [tuple(reversed(rotation)) for rotation in rotations]
    return min(rotations + reversals)


def branch_name(mode: int) -> str:
    return f"A{mode}_k{mode}"


def residual_unit5_permutation(d_orbits: list[frozenset[tuple[int, ...]]], d_index: dict[tuple[int, ...], int]) -> list[int]:
    return [d_index[affine_image(min(orbit), 0, 5)] for orbit in d_orbits]


def two_cycles(unit5_perm: list[int]) -> list[tuple[int, int]]:
    visited: set[int] = set()
    cycles: list[tuple[int, int]] = []
    for index, image in enumerate(unit5_perm):
        if index in visited:
            continue
        if image == index:
            visited.add(index)
            continue
        low, high = sorted((index, image))
        cycles.append((low, high))
        visited.update({low, high})
    return cycles


def shell_imbalance(histogram: tuple[int, int, int, int, int, int]) -> int:
    return histogram[4] - histogram[0]


def cycle_rows(cycles: list[tuple[int, int]], d_orbits: list[frozenset[tuple[int, ...]]]) -> list[dict[str, Any]]:
    rows = []
    for cycle_index, (low, high) in enumerate(cycles):
        low_support = min(d_orbits[low])
        high_support = min(d_orbits[high])
        low_hist = distance_histogram(low_support)
        high_hist = distance_histogram(high_support)
        low_imbalance = shell_imbalance(low_hist)
        high_imbalance = shell_imbalance(high_hist)
        rows.append(
            {
                "cycle_index": cycle_index,
                "orbit_pair": [low, high],
                "low_representative_support": list(low_support),
                "high_representative_support": list(high_support),
                "low_gap_necklace": list(gap_necklace(low_support)),
                "high_gap_necklace": list(gap_necklace(high_support)),
                "low_histogram_d1_to_d6": list(low_hist),
                "high_histogram_d1_to_d6": list(high_hist),
                "low_d5_minus_d1": low_imbalance,
                "high_d5_minus_d1": high_imbalance,
                "imbalance_flips_under_unit5": high_imbalance == -low_imbalance,
                "absolute_shell_imbalance": abs(low_imbalance),
                "is_branch_cycle": {low, high} == {0, 37},
            }
        )
    return rows


def branch_rows(d_index: dict[tuple[int, ...], int]) -> list[dict[str, Any]]:
    rows = []
    for mode in BRANCH_MODES:
        name = branch_name(mode)
        support = unit_support(mode)
        histogram = distance_histogram(support)
        rows.append(
            {
                "name": name,
                "support": list(support),
                "dihedral_orbit_index": d_index[support],
                "histogram_d1_to_d6": list(histogram),
                "d5_minus_d1": shell_imbalance(histogram),
                "required_chi11_value": 1 if name in {"A5_k5", "A7_k7"} else -1,
            }
        )
    return rows


def exact_proof_certificate(max_abs: int, max_count: int, distribution: dict[int, int]) -> dict[str, str]:
    return {
        "finite_domain": "All C(12,5)=792 supports are enumerated, quotiented by D12, and reduced to the unit-5 two-cycle basis.",
        "covariance_fact": "The shell-labelled imbalance I=h_d5-h_d1 satisfies I(unit5*O)=-I(O) on every D12 two-cycle.",
        "maximality_certificate": f"Across the 13 chi_11 two-cycles, max |I| is {max_abs} and is attained by {max_count} cycle(s); the amplitude distribution is {distribution}.",
        "branch_selection": "The unique maximizer is the branch cycle [0,37], i.e. A1/A11 versus A5/A7, so max |h_d5-h_d1| conditionally selects the branch generator.",
        "import_boundary": "The rule is not full-Aut strict provenance because it uses both the D12 quotient and the labelled d1/d5 shell axis.",
        "full_aut_intersection": "If full-Aut invariance is restored, d1 and d5 are exchanged and this oriented imbalance is not an invariant source.",
    }


def build_payload() -> dict[str, Any]:
    supports = all_supports()
    d_orbits = orbit_partition(supports, DIHEDRAL_UNITS)
    d_index = index_by_support(d_orbits)
    cycles = two_cycles(residual_unit5_permutation(d_orbits, d_index))
    rows = cycle_rows(cycles, d_orbits)
    distribution = dict(sorted(Counter(row["absolute_shell_imbalance"] for row in rows).items()))
    max_abs = max(distribution)
    max_rows = [row for row in rows if row["absolute_shell_imbalance"] == max_abs]
    branch_max_rows = [row for row in max_rows if row["is_branch_cycle"]]
    assert len(branch_max_rows) == 1
    return {
        "result_kind": "SCRATCH_STRICT_ALPHA_D12_CHI11_MAX_SHELL_IMBALANCE_SELECTOR_CERTIFICATE__CONDITIONAL_NOT_A_THEOREM",
        "status": "branch-generator-uniquely-maximizes-shell-labelled-d1-d5-imbalance-but-imports-axis",
        "finite_model": {
            "ring": "Z_12",
            "active_count": ACTIVE_COUNT,
            "support_count": len(supports),
            "dihedral_units": DIHEDRAL_UNITS,
            "automorphism_units": UNITS,
            "d12_orbit_count": len(d_orbits),
            "unit5_two_cycle_count": len(cycles),
            "target_q_power": TARGET_Q_POWER,
            "target_eta": TARGET_ETA,
        },
        "selector_summary": {
            "candidate_score": "abs(h_d5-h_d1) on each D12/unit5 two-cycle",
            "score_requires_shell_label": True,
            "score_requires_reduced_d12_quotient": True,
            "max_absolute_shell_imbalance": max_abs,
            "maximizer_count": len(max_rows),
            "unique_maximizer_is_branch_cycle": len(branch_max_rows) == 1 and len(max_rows) == 1,
            "absolute_shell_imbalance_distribution": distribution,
            "full_aut_allowed_strict_source": False,
        },
        "branch_maximizer_certificate": branch_max_rows[0],
        "branch_rows": branch_rows(d_index),
        "cycle_score_rows": rows,
        "exact_proof_certificate": exact_proof_certificate(max_abs, len(max_rows), distribution),
        "interpretation": {
            "honest_positive": "Within the already reduced D12 + shell-labelled setting, the branch generator is not arbitrary: it is the unique maximum shell-imbalance generator.",
            "honest_negative": "The selector score itself imports the d1/d5 shell axis and reduced D12 quotient, so it is not a strict full-Aut source of chi_11.",
            "relation_to_previous_probe": "The previous probe computed the rank-13 module; this one audits a finite extremal criterion for selecting one generator from that module.",
        },
        "ontology_guardrail": {
            "allowed_reading": "The nadsoliton itself remains the primordial information in a solitonic state.",
            "forbidden_reading": "No separate informational layer underneath the nadsoliton is introduced.",
        },
        "hard_limits": [
            "No identity K_legacy_ont == K_strict_gate is asserted.",
            "No legacy physical-role transfer onto K_strict_gate is used.",
            "No theorem derives chi_11, shell-labels, unit-axis bit, exact-cover clauses, or cardinality 5 from strict geometry.",
            "No QW-2191 discharge is claimed.",
            "No ToE closure is claimed.",
            "Result is conditional on D12 quotient plus shell-labelled d1/d5 data; it is not a full-Aut strict-source theorem.",
        ],
    }


def write_markdown(payload: dict[str, Any]) -> str:
    model = payload["finite_model"]
    summary = payload["selector_summary"]
    branch = payload["branch_maximizer_certificate"]
    lines = [
        "# D12 chi_11 max shell-imbalance selector certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Finite model",
        "",
        f"- Ring: `{model['ring']}`",
        f"- Enumerated supports: `{model['support_count']}`",
        f"- D12 orbit count: `{model['d12_orbit_count']}`",
        f"- Unit-5 two-cycle count: `{model['unit5_two_cycle_count']}`",
        "",
        "## Selector summary",
        "",
        f"- Candidate score: `{summary['candidate_score']}`",
        f"- Requires shell label: `{summary['score_requires_shell_label']}`",
        f"- Requires reduced D12 quotient: `{summary['score_requires_reduced_d12_quotient']}`",
        f"- Max |h_d5-h_d1|: `{summary['max_absolute_shell_imbalance']}`",
        f"- Maximizer count: `{summary['maximizer_count']}`",
        f"- Unique maximizer is branch cycle: `{summary['unique_maximizer_is_branch_cycle']}`",
        f"- Amplitude distribution: `{summary['absolute_shell_imbalance_distribution']}`",
        f"- Full-Aut allowed strict source: `{summary['full_aut_allowed_strict_source']}`",
        "",
        "## Branch maximizer",
        "",
        f"- Cycle index: `{branch['cycle_index']}`",
        f"- Orbit pair: `{branch['orbit_pair']}`",
        f"- Low/high histograms: `{branch['low_histogram_d1_to_d6']}` / `{branch['high_histogram_d1_to_d6']}`",
        f"- Low/high d5-d1: `{branch['low_d5_minus_d1']}` / `{branch['high_d5_minus_d1']}`",
        f"- Low/high gap necklaces: `{branch['low_gap_necklace']}` / `{branch['high_gap_necklace']}`",
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
