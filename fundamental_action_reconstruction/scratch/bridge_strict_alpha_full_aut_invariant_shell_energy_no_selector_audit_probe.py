#!/usr/bin/env python3
"""Scratch probe: full-Aut invariant shell energies cannot select A5/d5.

The previous dual-symmetry audit showed that the successful d1+d6 selector has
stabilizer {1,11}, i.e. the chi_11 kernel, because full Aut(Z_12) swaps shell d1
with shell d5.  This probe closes the complementary no-go:

    If a linear folded-shell energy is required to be invariant under the full
    unit action Aut(Z_12), can it uniquely select A5/d5?

No.  Full-Aut invariance imposes w1=w5 on shell weights.  Since
A1/contiguous has histogram [4,3,2,1,0,0] and A5/d5 has [0,3,2,1,4,0], every
full-Aut-invariant linear shell energy assigns them the same score.  Therefore
no such scalar energy can distinguish A5 from A1, let alone uniquely select A5.

The probe also enumerates all 31 nonzero binary full-Aut-invariant shell masks
and confirms that none uniquely selects the A5/d5 orbit.  This is a no-go audit,
not a selector theorem and not a QW-2191 discharge.
"""
from __future__ import annotations

import json
from collections import Counter
from itertools import combinations, product
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
OUT_JSON = HERE / "bridge_strict_alpha_full_aut_invariant_shell_energy_no_selector_audit_report.json"
OUT_MD = HERE / "bridge_strict_alpha_full_aut_invariant_shell_energy_no_selector_audit_report.md"

N = 12
ACTIVE_COUNT = 5
UNITS = [1, 5, 7, 11]
TARGET_Q_POWER = "256/243"
TARGET_ETA = "9/5"
A1_SUPPORT = tuple(range(5))
A5_SUPPORT = (0, 5, 10, 3, 8)
A1_HISTOGRAM = [4, 3, 2, 1, 0, 0]
A5_HISTOGRAM = [0, 3, 2, 1, 4, 0]
SHELL_NAMES = ["d1_nearest", "d2", "d3", "d4", "d5", "d6_antipodal"]


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


def shell_image(unit: int, shell_index_1_based: int) -> int:
    return folded(unit * shell_index_1_based)


def shell_permutation(unit: int) -> tuple[int, int, int, int, int, int]:
    return tuple(shell_image(unit, shell) for shell in range(1, N // 2 + 1))  # type: ignore[return-value]


def transform_weight(weight: tuple[int, ...], unit: int) -> tuple[int, ...]:
    transformed = [0] * len(weight)
    for source_index, coefficient in enumerate(weight, start=1):
        target_shell = shell_image(unit, source_index)
        transformed[target_shell - 1] = coefficient
    return tuple(transformed)


def is_full_aut_invariant_weight(weight: tuple[int, ...]) -> bool:
    return all(transform_weight(weight, unit) == weight for unit in UNITS)


def weight_formula(weight: tuple[int, ...]) -> str:
    terms = [SHELL_NAMES[index] for index, coefficient in enumerate(weight) if coefficient]
    return " + ".join(terms) if terms else "0"


def score(weight: tuple[int, ...], support: tuple[int, ...]) -> int:
    histogram = distance_histogram(support)
    return sum(coefficient * histogram[index] for index, coefficient in enumerate(weight))


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
        histogram = distance_histogram(representative)
        rows.append(
            {
                "representative": list(representative),
                "orbit_size": len(members),
                "distance_histogram_d1_to_d6": list(histogram),
                "is_A1_contiguous_orbit": canonical_support(A1_SUPPORT) in orbit,
                "is_A5_d5_orbit": canonical_support(A5_SUPPORT) in orbit,
            }
        )
    return sorted(rows, key=lambda row: row["representative"])


def invariant_weight_structure_audit() -> dict[str, Any]:
    unit_rows = [
        {"unit": unit, "shell_permutation_source_to_target": list(shell_permutation(unit))}
        for unit in UNITS
    ]
    symbolic_difference = {
        "A1_minus_A5_histogram": [A1_HISTOGRAM[index] - A5_HISTOGRAM[index] for index in range(6)],
        "dot_with_general_weight": "4*w1 - 4*w5",
        "full_aut_invariance_condition": "w1 = w5",
        "therefore_A1_A5_scores_equal": True,
    }
    return {
        "unit_shell_action_rows": unit_rows,
        "shell_orbits_under_full_Aut": [[1, 5], [2], [3], [4], [6]],
        "linear_invariance_condition": "w1=w5; w2,w3,w4,w6 free",
        "A1_A5_symbolic_score_difference": symbolic_difference,
    }


def binary_invariant_masks() -> list[tuple[int, ...]]:
    masks = []
    for d15, d2, d3, d4, d6 in product([0, 1], repeat=5):
        weight = (d15, d2, d3, d4, d15, d6)
        if any(weight):
            masks.append(weight)
    return masks


def binary_mask_minimizer_audit(supports: list[tuple[int, ...]]) -> dict[str, Any]:
    rows = []
    unique_a5_count = 0
    masks_selecting_a1_or_a5 = []
    for weight in binary_invariant_masks():
        assert is_full_aut_invariant_weight(weight)
        scored = [(score(weight, support), support) for support in supports]
        minimum = min(value for value, _support in scored)
        winners = [support for value, support in scored if value == minimum]
        orbit_rows = orbit_representatives(winners)
        hist_counter = Counter(distance_histogram(support) for support in winners)
        selects_a5 = len(orbit_rows) == 1 and orbit_rows[0]["is_A5_d5_orbit"]
        selects_a1 = len(orbit_rows) == 1 and orbit_rows[0]["is_A1_contiguous_orbit"]
        if selects_a5:
            unique_a5_count += 1
        if any(row["is_A1_contiguous_orbit"] or row["is_A5_d5_orbit"] for row in orbit_rows):
            masks_selecting_a1_or_a5.append(weight_formula(weight))
        rows.append(
            {
                "weight": list(weight),
                "formula": weight_formula(weight),
                "minimum_score": minimum,
                "winner_support_count": len(winners),
                "winner_dihedral_orbit_count": len(orbit_rows),
                "winner_histogram_count": len(hist_counter),
                "contains_A1_contiguous_orbit": any(row["is_A1_contiguous_orbit"] for row in orbit_rows),
                "contains_A5_d5_orbit": any(row["is_A5_d5_orbit"] for row in orbit_rows),
                "selects_unique_A1_contiguous_orbit": selects_a1,
                "selects_unique_A5_d5_orbit": selects_a5,
            }
        )
    return {
        "binary_full_aut_invariant_mask_count": len(rows),
        "unique_A5_selecting_mask_count": unique_a5_count,
        "mask_formulas_touching_A1_or_A5_winners": masks_selecting_a1_or_a5,
        "rows": rows,
    }


def exact_proof_certificate() -> dict[str, str]:
    return {
        "full_aut_shell_action": "The full unit action fixes d2,d3,d4,d6 and swaps d1 with d5 via units 5 and 7.",
        "invariant_weight_condition": "A linear shell energy is full-Aut invariant only if w1=w5.",
        "A1_A5_pair_obstruction": "A1-A5 histogram difference is [4,0,0,0,-4,0], so every full-Aut invariant linear shell energy has zero A1/A5 score difference.",
        "binary_enumeration": "All 31 nonzero binary full-Aut invariant shell masks were enumerated; none uniquely selects A5/d5.",
        "source_obstruction": "Any shell-linear selector that distinguishes A5 from A1 must break full Aut at least to the {1,11} chi_11 kernel or add a non-shell-linear strict source not audited here.",
    }


def main() -> None:
    supports = all_supports()
    histograms = {distance_histogram(support) for support in supports}
    structure = invariant_weight_structure_audit()
    binary_audit = binary_mask_minimizer_audit(supports)

    report: dict[str, Any] = {
        "result_kind": "SCRATCH_STRICT_ALPHA_FULL_AUT_INVARIANT_SHELL_ENERGY_NO_SELECTOR_AUDIT_PROBE__NO_GO_NOT_A_THEOREM",
        "status": "full-aut-invariant-linear-shell-energies-cannot-distinguish-A5-from-A1",
        "finite_model": {
            "ring": "Z_12",
            "active_count": ACTIVE_COUNT,
            "support_count": len(supports),
            "histogram_class_count": len(histograms),
            "automorphism_units": UNITS,
            "target_q_power": TARGET_Q_POWER,
            "target_eta": TARGET_ETA,
            "A1_histogram_d1_to_d6": A1_HISTOGRAM,
            "A5_histogram_d1_to_d6": A5_HISTOGRAM,
        },
        "invariant_weight_structure_audit": structure,
        "binary_full_aut_invariant_mask_minimizer_audit": binary_audit,
        "exact_proof_certificate": exact_proof_certificate(),
        "interpretation": {
            "direct_no_go": "Full-Aut invariant shell-linear scalar energies cannot distinguish A5/d5 from A1/contiguous because full Aut forces w1=w5.",
            "computational_confirmation": "The 31 nonzero binary full-Aut invariant shell masks contain zero unique A5 selectors.",
            "relation_to_previous_probe": "The previous d1+d6 selector works precisely because it is not full-Aut invariant; it has chi_11-kernel stabilizer.",
            "honest_limit": "This no-go covers shell-linear full-Aut invariant energies only; it does not rule out a new non-shell-linear strict source theorem.",
        },
        "ontology_guardrail": {
            "allowed_reading": "The nadsoliton itself may later export a non-shell-linear or symmetry-breaking source.",
            "forbidden_reading": "No separate informational layer underneath the nadsoliton is introduced to bypass the full-Aut obstruction.",
            "preferred_order_preserved": "nadsoliton -> light -> matter -> emergent observer",
        },
        "hard_limits": [
            "No identity K_legacy_ont == K_strict_gate is asserted or used.",
            "No legacy physical-role or matter-generation claim is transferred onto K_strict_gate.",
            "No theorem derives a full-Aut invariant selector, chi_11, or the unit-axis bit from strict nadsoliton geometry.",
            "This is a shell-linear no-go audit, not a global impossibility theorem for all conceivable strict sources.",
            "No endpoint, arrow orientation, ledger selector, positive lambda action, cycle metric source, anti-Nyquist source, fifth-mode source, future-probability source, future-value source, future-path source, matter-bit source, existence-bit source, recursive-self-information source, character-source, meta-character source, cyclic-adjacency source theorem, variational tie-break source theorem, exclusion-energy source theorem, exact-cover source theorem, dual-shell-label source theorem, or legacy-strict bridge theorem is claimed.",
            "No QW-2191 discharge and no strict-core selector closure are claimed.",
            "No ToE closure is claimed.",
        ],
        "next_honest_step": "Either find a strict source that breaks full Aut to the chi_11 kernel, or leave A5/d5 selection explicitly conditional; full-Aut invariant shell-linear scalar data cannot do the job.",
    }

    OUT_JSON.write_text(json.dumps(report, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    OUT_MD.write_text(
        "# Scratch strict-alpha full-Aut invariant shell-energy no-selector audit probe\n\n"
        "Status: full-Aut invariant shell-linear energies cannot distinguish A5/d5 from A1/contiguous.\n\n"
        f"- Supports scanned: `{len(supports)}`; histogram classes: `{len(histograms)}`.\n"
        f"- Shell orbits under full Aut: `{structure['shell_orbits_under_full_Aut']}`.\n"
        f"- Linear invariance condition: `{structure['linear_invariance_condition']}`.\n"
        f"- A1-A5 score difference under full Aut invariant weights: `0`.\n"
        f"- Binary full-Aut invariant masks enumerated: `{binary_audit['binary_full_aut_invariant_mask_count']}`.\n"
        f"- Unique A5-selecting full-Aut invariant masks: `{binary_audit['unique_A5_selecting_mask_count']}`.\n"
        f"- Target replay kept conditional: `q^5={TARGET_Q_POWER}`, eta `{TARGET_ETA}`.\n"
        "- No false pass: no full-Aut shell-linear selector, no QW-2191 discharge, no ToE closure.\n",
        encoding="utf-8",
    )
    print(OUT_JSON)
    print(OUT_MD)


if __name__ == "__main__":
    main()
