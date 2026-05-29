#!/usr/bin/env python3
"""Scratch probe: the d1+d6 selector carries a chi_11 shell-label premise.

The previous exact-cover audit gave a strong finite certificate: forbidding
nearest-neighbour and antipodal contacts selects the A5/d5 orbit.  This probe
checks the dual symmetry cost of that success.

Aut(Z_12)={1,5,7,11} acts on folded shell coordinates d1..d6.  Units 5 and 7
swap shell d1 with shell d5, while units 1 and 11 preserve them.  Therefore the
weight vector for E_A5=d1+d6 has stabilizer {1,11}: exactly the chi_11 kernel.
Its full Aut orbit has two weights:

    d1+d6 -> selects A5/d5,
    d5+d6 -> selects A1/contiguous.

So E_A5 is a useful conditional selector, but choosing E_A5 rather than its
Aut-conjugate E_A1 already imports the same unit-label/character datum that the
selector problem was missing.  This is not a failure of the finite certificate;
it is the honest source accounting for why the certificate is conditional.
"""
from __future__ import annotations

import json
from collections import Counter
from itertools import combinations
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
OUT_JSON = HERE / "bridge_strict_alpha_exclusion_energy_dual_symmetry_source_obstruction_report.json"
OUT_MD = HERE / "bridge_strict_alpha_exclusion_energy_dual_symmetry_source_obstruction_report.md"

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
E_A5_WEIGHT = (1, 0, 0, 0, 0, 1)


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
    # Entry i says source shell i+1 maps to target shell entry.
    return tuple(shell_image(unit, shell) for shell in range(1, N // 2 + 1))  # type: ignore[return-value]


def transform_weight(weight: tuple[int, ...], unit: int) -> tuple[int, ...]:
    transformed = [0] * len(weight)
    for source_index, coefficient in enumerate(weight, start=1):
        target_shell = shell_image(unit, source_index)
        transformed[target_shell - 1] = coefficient
    return tuple(transformed)


def weight_name(weight: tuple[int, ...]) -> str:
    if weight == (1, 0, 0, 0, 0, 1):
        return "E_A5_candidate_d1_plus_d6"
    if weight == (0, 0, 0, 0, 1, 1):
        return "E_A1_conjugate_d5_plus_d6"
    return "weight_" + "_".join(str(value) for value in weight)


def weight_formula(weight: tuple[int, ...]) -> str:
    terms = [SHELL_NAMES[index] for index, coefficient in enumerate(weight) if coefficient]
    return " + ".join(terms) if terms else "0"


def score(weight: tuple[int, ...], support: tuple[int, ...]) -> int:
    hist = distance_histogram(support)
    return sum(coefficient * hist[index] for index, coefficient in enumerate(weight))


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
                "is_A1_contiguous_orbit": canonical_support(A1_SUPPORT) in orbit,
                "is_A5_d5_orbit": canonical_support(A5_SUPPORT) in orbit,
            }
        )
    return sorted(rows, key=lambda row: row["representative"])


def weight_minimizer_audit(weight: tuple[int, ...], supports: list[tuple[int, ...]]) -> dict[str, Any]:
    scored = [(score(weight, support), support) for support in supports]
    minimum = min(value for value, _support in scored)
    winners = [support for value, support in scored if value == minimum]
    hist_counter = Counter(distance_histogram(support) for support in winners)
    orbit_rows = orbit_representatives(winners)
    return {
        "weight": list(weight),
        "name": weight_name(weight),
        "formula": weight_formula(weight),
        "minimum_score": minimum,
        "winner_support_count": len(winners),
        "winner_histogram_rows": [
            {"distance_histogram_d1_to_d6": list(hist), "support_count": count}
            for hist, count in sorted(hist_counter.items())
        ],
        "winner_dihedral_orbit_count": len(orbit_rows),
        "winner_orbits": orbit_rows,
        "selects_A1_contiguous_orbit": len(orbit_rows) == 1 and orbit_rows[0]["is_A1_contiguous_orbit"],
        "selects_A5_d5_orbit": len(orbit_rows) == 1 and orbit_rows[0]["is_A5_d5_orbit"],
    }


def shell_action_audit() -> dict[str, Any]:
    rows = []
    for unit in UNITS:
        perm = shell_permutation(unit)
        transformed = transform_weight(E_A5_WEIGHT, unit)
        rows.append(
            {
                "unit": unit,
                "shell_permutation_source_to_target": list(perm),
                "E_A5_maps_to_weight": list(transformed),
                "E_A5_maps_to_name": weight_name(transformed),
                "preserves_E_A5": transformed == E_A5_WEIGHT,
            }
        )
    stabilizer = [row["unit"] for row in rows if row["preserves_E_A5"]]
    orbit = sorted({tuple(row["E_A5_maps_to_weight"]) for row in rows})
    return {
        "rows": rows,
        "E_A5_weight": list(E_A5_WEIGHT),
        "E_A5_formula": weight_formula(E_A5_WEIGHT),
        "stabilizer_units": stabilizer,
        "stabilizer_is_chi_11_kernel": stabilizer == [1, 11],
        "Aut_orbit_size": len(orbit),
        "Aut_orbit_weights": [
            {"weight": list(weight), "name": weight_name(weight), "formula": weight_formula(weight)}
            for weight in orbit
        ],
    }


def exact_proof_certificate() -> dict[str, str]:
    return {
        "shell_action": "Units 5 and 7 send folded shell d1 to d5 and d5 to d1; units 1 and 11 preserve d1 and d5.",
        "stabilizer": "The weight E_A5=d1+d6 is fixed exactly by {1,11}, the chi_11 kernel.",
        "orbit_fact": "The full Aut orbit of E_A5 has two weights: d1+d6 and d5+d6.",
        "selector_conjugates": "Minimizing d1+d6 selects A5/d5; minimizing the Aut-conjugate d5+d6 selects A1/contiguous.",
        "source_obstruction": "Choosing d1+d6 rather than d5+d6 is therefore already a unit-label/chi_11 premise, not a strict-core consequence proven here.",
    }


def main() -> None:
    supports = all_supports()
    histograms = {distance_histogram(support) for support in supports}
    action = shell_action_audit()
    minimizer_rows = [
        weight_minimizer_audit(tuple(row["weight"]), supports)
        for row in action["Aut_orbit_weights"]
    ]

    report: dict[str, Any] = {
        "result_kind": "SCRATCH_STRICT_ALPHA_EXCLUSION_ENERGY_DUAL_SYMMETRY_SOURCE_OBSTRUCTION_PROBE__NOT_A_THEOREM",
        "status": "d1_plus_d6-selector-success-imports-chi_11-shell-label-premise",
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
        "shell_action_audit": action,
        "aut_conjugate_selector_audit": minimizer_rows,
        "exact_proof_certificate": exact_proof_certificate(),
        "interpretation": {
            "direct_result": "The successful d1+d6 selector is not full-Aut invariant; its stabilizer is exactly {1,11}, the chi_11 kernel.",
            "conjugate_warning": "The Aut-conjugate selector d5+d6 is equally natural under abstract unit symmetry and selects A1/contiguous instead of A5/d5.",
            "honest_limit": "The finite selector remains conditional unless strict nadsoliton geometry exports why the d1 shell, not the d5 shell, is the penalized neighbour shell.",
            "research_value": "This identifies the exact symmetry cost of the exclusion-energy certificate and prevents mistaking a conditional selector for a strict-core discharge.",
        },
        "ontology_guardrail": {
            "allowed_reading": "The nadsoliton itself may later export the shell-label premise that distinguishes d1 from d5.",
            "forbidden_reading": "No separate informational layer underneath the nadsoliton is introduced to choose d1+d6 over d5+d6.",
            "preferred_order_preserved": "nadsoliton -> light -> matter -> emergent observer",
        },
        "hard_limits": [
            "No identity K_legacy_ont == K_strict_gate is asserted or used.",
            "No legacy physical-role or matter-generation claim is transferred onto K_strict_gate.",
            "No theorem derives the d1-vs-d5 shell label, chi_11, or the unit-axis bit from strict nadsoliton geometry.",
            "The d1+d6 exclusion selector is classified as chi_11-kernel-stabilized conditional data, not a strict-core theorem.",
            "No endpoint, arrow orientation, ledger selector, positive lambda action, cycle metric source, anti-Nyquist source, fifth-mode source, future-probability source, future-value source, future-path source, matter-bit source, existence-bit source, recursive-self-information source, character-source, meta-character source, cyclic-adjacency source theorem, variational tie-break source theorem, exclusion-energy source theorem, exact-cover source theorem, or legacy-strict bridge theorem is claimed.",
            "No QW-2191 discharge and no strict-core selector closure are claimed.",
            "No ToE closure is claimed.",
        ],
        "next_honest_step": "Search for a strict source that distinguishes shell d1 from its Aut-conjugate d5; absent that, keep d1+d6 as an explicit chi_11-shell-label premise.",
    }

    OUT_JSON.write_text(json.dumps(report, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    OUT_MD.write_text(
        "# Scratch strict-alpha exclusion-energy dual-symmetry source obstruction probe\n\n"
        "Status: d1+d6 selector success imports a chi_11 shell-label premise.\n\n"
        f"- Supports scanned: `{len(supports)}`; histogram classes: `{len(histograms)}`.\n"
        f"- E_A5 formula: `{action['E_A5_formula']}`.\n"
        f"- E_A5 stabilizer units: `{action['stabilizer_units']}`.\n"
        f"- Stabilizer is chi_11 kernel: `{action['stabilizer_is_chi_11_kernel']}`.\n"
        f"- Full Aut orbit size of E_A5: `{action['Aut_orbit_size']}`.\n"
        f"- Aut orbit weights: `{action['Aut_orbit_weights']}`.\n"
        f"- Conjugate selector outcomes: `{[(row['name'], row['selects_A1_contiguous_orbit'], row['selects_A5_d5_orbit']) for row in minimizer_rows]}`.\n"
        f"- Target replay kept conditional: `q^5={TARGET_Q_POWER}`, eta `{TARGET_ETA}`.\n"
        "- No false pass: no strict d1-vs-d5 shell-label theorem, no QW-2191 discharge, no ToE closure.\n",
        encoding="utf-8",
    )
    print(OUT_JSON)
    print(OUT_MD)


if __name__ == "__main__":
    main()
