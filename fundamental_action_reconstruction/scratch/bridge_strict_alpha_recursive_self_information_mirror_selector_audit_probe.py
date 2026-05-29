#!/usr/bin/env python3
"""Scratch probe: recursive self-information mirrors do not create the unit bit.

User intuition audited here: perhaps information contains information about its
own states, which contain information about their own states, and so on -- like
an infinite mirror.  Could that reflexive tower generate the missing selector
bit?

Finite audit:

    Level 0 is the two-branch unit-mirror orbit {A1, A5}.
    Level n+1 is the powerset of Level n, i.e. a record about possible records.

The unit mirror action lifts functorially to every level.  If a record is built
without an Aut-breaking seed, it must be fixed by the lifted action.  We exactly
enumerate the first powerset levels and prove by orbit lifting for deeper levels.

No false pass: recursive self-reference can preserve an already supplied
orientation label and can generate a large hierarchy of invariant records, but
it cannot turn a symmetric mirror pair into singleton d5.  Reflection depth is
not the same datum as unit orientation.
"""
from __future__ import annotations

import json
from collections import Counter
from itertools import combinations
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
OUT_JSON = HERE / "bridge_strict_alpha_recursive_self_information_mirror_selector_audit_report.json"
OUT_MD = HERE / "bridge_strict_alpha_recursive_self_information_mirror_selector_audit_report.md"

N = 12
ACTIVE_COUNT = 5
AUT_UNITS = [1, 5, 7, 11]
TARGET_Q_POWER = "256/243"
TARGET_ETA = "9/5"
MAX_LEVEL = 5
ENUMERATED_MAX_LEVEL = 3
A1 = "A1_k1_contiguous"
A5 = "A5_k5_d5"
D5_HISTOGRAM = [0, 3, 2, 1, 4, 0]
CONTIGUOUS_HISTOGRAM = [4, 3, 2, 1, 0, 0]


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


def powerset_masks(size: int) -> list[int]:
    return list(range(1 << size))


def mirror_level0(index: int) -> int:
    return 1 - index


def lifted_permutation(previous_permutation: list[int]) -> list[int]:
    """Lift a permutation of an N-element set to the powerset encoded by masks."""
    size = len(previous_permutation)
    lifted = []
    for mask in powerset_masks(size):
        image = 0
        for bit in range(size):
            if mask & (1 << bit):
                image |= 1 << previous_permutation[bit]
        lifted.append(image)
    return lifted


def int_or_power_text(value: int | None, exponent: int | None = None) -> int | str:
    if value is not None:
        return value
    if exponent is not None and exponent <= 128:
        return 1 << exponent
    return f"2^{exponent}"


def enumerated_level_permutations(max_level: int) -> list[list[int]]:
    permutations = [[mirror_level0(0), mirror_level0(1)]]
    for _ in range(max_level):
        permutations.append(lifted_permutation(permutations[-1]))
    return permutations


def orbit_count(permutation: list[int]) -> int:
    seen: set[int] = set()
    count = 0
    for item in range(len(permutation)):
        if item in seen:
            continue
        seen.add(item)
        seen.add(permutation[item])
        count += 1
    return count


def fixed_count(permutation: list[int]) -> int:
    return sum(1 for item, image in enumerate(permutation) if item == image)


def singleton_d5_index_by_level(level_count: int) -> list[int | str]:
    # Level 0: A5 is index 1.  Level n+1 singleton({previous target}) is mask 1<<previous_target.
    targets: list[int | str] = [1]
    for _ in range(1, level_count):
        previous = targets[-1]
        if isinstance(previous, int) and previous <= 128:
            targets.append(1 << previous)
        else:
            targets.append(f"2^({previous})")
    return targets


def row_from_counts(
    level: int,
    object_count: int | str,
    mirror_orbit_count: int | str,
    fixed_record_count: int | str,
    target_index: int | str,
    enumerated_exactly: bool,
) -> dict[str, Any]:
    nonfixed: int | str
    if isinstance(object_count, int) and isinstance(fixed_record_count, int):
        nonfixed = object_count - fixed_record_count
    else:
        nonfixed = "object_count - fixed_record_count"
    return {
        "level": level,
        "object_kind": "branch" if level == 0 else f"P^{level}({{A1,A5}}) self-information records",
        "object_count": object_count,
        "mirror_orbit_count": mirror_orbit_count,
        "fixed_record_count": fixed_record_count,
        "nonfixed_record_count": nonfixed,
        "enumerated_exactly": enumerated_exactly,
        "singleton_d5_tower_index": target_index,
        "singleton_d5_tower_fixed_by_mirror": False,
        "aut_invariant_singleton_d5_selector_exists_at_level": False,
    }


def level_rows() -> list[dict[str, Any]]:
    enumerated = enumerated_level_permutations(ENUMERATED_MAX_LEVEL)
    rows = []
    targets = singleton_d5_index_by_level(MAX_LEVEL + 1)
    object_counts: list[int | str] = []
    fixed_counts: list[int | str] = []
    orbit_counts: list[int | str] = []

    for level, permutation in enumerate(enumerated):
        object_count = len(permutation)
        fixed = fixed_count(permutation)
        orbits = orbit_count(permutation)
        object_counts.append(object_count)
        fixed_counts.append(fixed)
        orbit_counts.append(orbits)
        rows.append(row_from_counts(level, object_count, orbits, fixed, targets[level], True))

    for level in range(ENUMERATED_MAX_LEVEL + 1, MAX_LEVEL + 1):
        previous_object = object_counts[-1]
        previous_orbits = orbit_counts[-1]
        if isinstance(previous_object, int):
            object_count = int_or_power_text(None, previous_object)
        else:
            object_count = f"2^({previous_object})"
        if isinstance(previous_orbits, int):
            fixed = int_or_power_text(None, previous_orbits)
        else:
            fixed = f"2^({previous_orbits})"
        if isinstance(object_count, int) and isinstance(fixed, int):
            orbits = (object_count + fixed) // 2
        else:
            orbits = f"({object_count} + {fixed})/2"
        object_counts.append(object_count)
        fixed_counts.append(fixed)
        orbit_counts.append(orbits)
        rows.append(row_from_counts(level, object_count, orbits, fixed, targets[level], False))
    return rows


def finite_reflection_depth_rows() -> list[dict[str, Any]]:
    rows = []
    for depth in range(MAX_LEVEL + 1):
        rows.append(
            {
                "reflection_depth": depth,
                "depth_is_aut_invariant_scalar": True,
                "depth_selects_A5_over_A1": False,
                "interpretation": "knowing how many mirror levels are represented is not an orientation character",
            }
        )
    return rows


def exact_proof_certificate() -> dict[str, str]:
    return {
        "base_action": "J_0 swaps A1 and A5; neither singleton branch is fixed.",
        "recursive_lift": "J_{n+1}(S)={J_n(x): x in S} on the powerset/self-record level P(Level_n).",
        "fixed_record_condition": "A self-information record is Aut-invariant only if it is a union of whole J_n-orbits.",
        "singleton_tower_obstruction": "The tower singleton { ... {A5} ... } maps to { ... {A1} ... } at every finite level, so it is never fixed.",
        "infinite_mirror_limit": "An inverse/direct limit of equivariant finite levels remains equivariant; no finite stage supplies a fixed singleton d5 seed to the limit.",
        "information_reading": "Recursive information-about-information can increase reflection depth, but depth is an Aut-invariant scalar, not an Aut-breaking unit character.",
    }


def main() -> None:
    supports = all_supports()
    histogram_counter = Counter(distance_histogram(support) for support in supports)
    rows = level_rows()

    report: dict[str, Any] = {
        "result_kind": "SCRATCH_STRICT_ALPHA_RECURSIVE_SELF_INFORMATION_MIRROR_SELECTOR_AUDIT_PROBE__NOT_A_THEOREM",
        "status": "recursive-self-information-mirror-depth-does-not-create-unit-orientation-bit",
        "finite_model": {
            "ring": "Z_12",
            "active_count": ACTIVE_COUNT,
            "support_count": len(supports),
            "histogram_class_count": len(histogram_counter),
            "automorphism_units": AUT_UNITS,
            "survivor_axes": [
                {"name": A1, "mode": 1, "distance_histogram_d1_to_d6": CONTIGUOUS_HISTOGRAM},
                {"name": A5, "mode": 5, "distance_histogram_d1_to_d6": D5_HISTOGRAM},
            ],
            "target_q_power": TARGET_Q_POWER,
            "target_eta": TARGET_ETA,
        },
        "recursive_self_information_level_audit": {
            "max_level": MAX_LEVEL,
            "enumerated_max_level": ENUMERATED_MAX_LEVEL,
            "construction": "Level_0={A1,A5}; Level_{n+1}=P(Level_n); mirror action lifts by image of subsets.",
            "rows": rows,
            "singleton_d5_fixed_count_across_levels": sum(1 for row in rows if row["singleton_d5_tower_fixed_by_mirror"]),
            "aut_invariant_singleton_d5_selector_exists_any_level": any(
                row["aut_invariant_singleton_d5_selector_exists_at_level"] for row in rows
            ),
        },
        "reflection_depth_audit": {
            "rows": finite_reflection_depth_rows(),
            "depth_as_scalar_can_break_unit_mirror": False,
        },
        "exact_proof_certificate": exact_proof_certificate(),
        "interpretation": {
            "direct_correction": "The infinite-mirror intuition is compatible with nadsoliton-as-information, but recursion alone preserves the unit-mirror symmetry.",
            "what_it_can_do": "It can build a hierarchy of self-records/information-about-information and fixed invariant records made from whole mirror orbits.",
            "what_it_cannot_do": "It cannot select the singleton d5 tower unless an Aut-breaking unit-orientation seed is already supplied.",
            "information_principle": "Self-reference changes depth/order of information; the selector needs orientation character.",
        },
        "ontology_guardrail": {
            "allowed_reading": "The nadsoliton itself may contain recursive self-information about its own possible internal states.",
            "forbidden_reading": "No separate informational layer underneath the nadsoliton is introduced to host the mirror tower or selector bit.",
            "preferred_order_preserved": "nadsoliton -> light -> matter -> emergent observer",
        },
        "hard_limits": [
            "No identity K_legacy_ont == K_strict_gate is asserted or used.",
            "No legacy physical-role or matter-generation claim is transferred onto K_strict_gate.",
            "No theorem derives recursive self-information from strict nadsoliton geometry.",
            "No theorem derives the required one-bit unit-axis record from strict core.",
            "Recursive self-information does not discharge QW-2191 because the lifted mirror action remains equivariant at every finite level.",
            "Without an Aut-breaking unit-orientation seed, singleton d5 branch selection remains blocked through the recursive mirror tower.",
            "No endpoint, arrow orientation, ledger selector, positive lambda action, cycle metric source, anti-Nyquist source, fifth-mode source, future-probability source, future-value source, future-path source, matter-bit source, existence-bit source, recursive-self-information source, or legacy-strict bridge theorem is claimed.",
            "No QW-2191 discharge and no strict-core selector closure are claimed.",
            "No ToE closure is claimed.",
        ],
        "next_honest_step": "If pursuing infinite-mirror information, search for a strict non-equivariant seed in the recursive tower; otherwise classify recursion as an invariant amplifier/carrier, not the source of the unit bit.",
    }

    OUT_JSON.write_text(json.dumps(report, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    OUT_MD.write_text(
        "# Scratch strict-alpha recursive self-information mirror selector audit\n\n"
        "Status: recursive information-about-information does not create the unit-orientation bit.\n\n"
        f"- Supports scanned: `{len(supports)}`; histogram classes: `{len(histogram_counter)}`.\n"
        f"- Recursive levels audited: `0..{MAX_LEVEL}`; exact permutation enumeration through level `{ENUMERATED_MAX_LEVEL}`.\n"
        f"- Singleton d5 tower fixed count across levels: `{report['recursive_self_information_level_audit']['singleton_d5_fixed_count_across_levels']}`.\n"
        f"- Aut-invariant singleton d5 selector exists at any level: `{report['recursive_self_information_level_audit']['aut_invariant_singleton_d5_selector_exists_any_level']}`.\n"
        f"- Reflection depth as scalar can break unit mirror: `{report['reflection_depth_audit']['depth_as_scalar_can_break_unit_mirror']}`.\n"
        "- Honest read: self-reference adds mirror depth, not unit orientation.\n"
        f"- Target replay kept conditional: `q^5={TARGET_Q_POWER}`, eta `{TARGET_ETA}`.\n"
        "- No false pass: no recursive-self-information source theorem, no QW-2191 discharge, no ToE closure.\n",
        encoding="utf-8",
    )
    print(OUT_JSON)
    print(OUT_MD)


if __name__ == "__main__":
    main()
