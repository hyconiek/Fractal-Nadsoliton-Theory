#!/usr/bin/env python3
"""Scratch probe: affine subgroup lattice obstruction for chi_11 sources.

Previous audits showed that D12 data can carry the chi_11 branch sign, while
full-Aut data glues A1/A11 and A5/A7 into one block and erases the sign.  This
probe makes the group-theoretic boundary exact by enumerating every affine unit
subgroup between translations and full Aut(Z_12):

    T⋊{1}, T⋊{1,11}, T⋊{1,5}, T⋊{1,7}, T⋊{1,5,7,11}.

For each subgroup we enumerate all C(12,5)=792 supports, compute orbit counts,
and test whether the branch-level chi_11 assignment is well-defined on the
quotient.  A quotient supports chi_11 iff no orbit contains both required signs
(A1/A11=-1 and A5/A7=+1).

Result: only the cyclic translation quotient and the D12 quotient support the
branch chi_11 sign.  Adjoining either non-D12 unit 5 or unit 7 immediately mixes
opposite signs in the branch orbit.  Thus the exact obstruction is not "lack of
a clever invariant" inside full Aut; it is the presence of either unit-5 or
unit-7 symmetry in the admitted quotient.

No false pass: this is a subgroup-lattice obstruction certificate.  It does not
derive a unit-axis premise, does not discharge QW-2191, and does not close ToE.
"""
from __future__ import annotations

import json
from collections import Counter
from itertools import combinations
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
OUT_JSON = HERE / "bridge_strict_alpha_affine_subgroup_lattice_chi11_source_obstruction_report.json"
OUT_MD = HERE / "bridge_strict_alpha_affine_subgroup_lattice_chi11_source_obstruction_report.md"

N = 12
ACTIVE_COUNT = 5
ALL_UNITS = [1, 5, 7, 11]
UNIT_SUBGROUPS = [
    ("T_semidirect_{1}_cyclic_translation", [1]),
    ("T_semidirect_{1,11}_D12", [1, 11]),
    ("T_semidirect_{1,5}_unit5_axis", [1, 5]),
    ("T_semidirect_{1,7}_unit7_axis", [1, 7]),
    ("T_semidirect_{1,5,7,11}_full_affine_Aut", [1, 5, 7, 11]),
]
BRANCH_MODES = [1, 5, 7, 11]
TARGET_Q_POWER = "256/243"
TARGET_ETA = "9/5"


def unit_support(mode: int) -> tuple[int, ...]:
    return tuple(sorted((mode * index) % N for index in range(ACTIVE_COUNT)))


def branch_name(mode: int) -> str:
    return f"A{mode}_k{mode}"


def branch_pair(name: str) -> str:
    return "d5_pair_A5_A7" if name in {"A5_k5", "A7_k7"} else "contiguous_pair_A1_A11"


def required_chi11_value(name: str) -> int:
    return 1 if branch_pair(name) == "d5_pair_A5_A7" else -1


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


def branch_rows_for_partition(orbits: list[frozenset[tuple[int, ...]]]) -> list[dict[str, Any]]:
    index_by_support = {support: index for index, orbit in enumerate(orbits) for support in orbit}
    rows = []
    for mode in BRANCH_MODES:
        name = branch_name(mode)
        support = unit_support(mode)
        rows.append(
            {
                "name": name,
                "support": list(support),
                "gap_necklace": list(gap_necklace(support)),
                "branch_pair": branch_pair(name),
                "required_chi11_value": required_chi11_value(name),
                "quotient_orbit_index": index_by_support[support],
            }
        )
    return rows


def subgroup_row(name: str, units: list[int], supports: list[tuple[int, ...]]) -> dict[str, Any]:
    orbits = orbit_partition(supports, units)
    branch_rows = branch_rows_for_partition(orbits)
    orbit_to_branch_values: dict[int, list[int]] = {}
    orbit_to_branch_names: dict[int, list[str]] = {}
    for row in branch_rows:
        orbit_index = row["quotient_orbit_index"]
        orbit_to_branch_values.setdefault(orbit_index, []).append(row["required_chi11_value"])
        orbit_to_branch_names.setdefault(orbit_index, []).append(row["name"])
    mixed_orbits = [
        {
            "orbit_index": orbit_index,
            "branch_names": orbit_to_branch_names[orbit_index],
            "required_chi11_values": values,
        }
        for orbit_index, values in sorted(orbit_to_branch_values.items())
        if len(set(values)) > 1
    ]
    supports_chi11 = not mixed_orbits
    orbit_size_counts = dict(sorted(Counter(len(orbit) for orbit in orbits).items()))
    return {
        "subgroup_name": name,
        "unit_subgroup": units,
        "contains_reversal_unit11": 11 in units,
        "contains_unit5": 5 in units,
        "contains_unit7": 7 in units,
        "orbit_count_on_all_supports": len(orbits),
        "orbit_size_counts": orbit_size_counts,
        "branch_rows": branch_rows,
        "branch_orbit_partition": {
            str(orbit_index): orbit_to_branch_names[orbit_index]
            for orbit_index in sorted(orbit_to_branch_names)
        },
        "mixed_chi11_sign_branch_orbits": mixed_orbits,
        "branch_chi11_well_defined_on_quotient": supports_chi11,
        "classification": (
            "admits_branch_chi11_quotient"
            if supports_chi11
            else "kills_branch_chi11_by_mixing_opposite_signs"
        ),
    }


def exact_proof_certificate(rows: list[dict[str, Any]]) -> dict[str, str]:
    admitting = [row["subgroup_name"] for row in rows if row["branch_chi11_well_defined_on_quotient"]]
    killing = [row["subgroup_name"] for row in rows if not row["branch_chi11_well_defined_on_quotient"]]
    return {
        "finite_domain": "All C(12,5)=792 supports are enumerated for every affine unit subgroup between translations and full Aut.",
        "well_definedness_test": "A quotient can host branch chi_11 exactly when no quotient orbit contains both signs A1/A11=-1 and A5/A7=+1.",
        "admitting_subgroups": f"The only admitting subgroups are {admitting}.",
        "killing_subgroups": f"The subgroups {killing} all mix opposite branch signs in one quotient orbit.",
        "unit_obstruction": "The obstruction appears exactly when the admitted unit subgroup contains unit 5 or unit 7; unit 11 alone is harmless because it preserves the chi_11 sign pairs.",
        "full_aut_consequence": "Full affine Aut contains both unit 5 and unit 7, so any full-Aut invariant support source is constant on the mixed branch orbit and cannot export chi_11 polarity.",
        "strict_limit": "This subgroup-lattice fact does not derive the missing unit-axis premise and does not discharge QW-2191.",
    }


def build_payload() -> dict[str, Any]:
    supports = all_supports()
    rows = [subgroup_row(name, units, supports) for name, units in UNIT_SUBGROUPS]
    admitting = [row for row in rows if row["branch_chi11_well_defined_on_quotient"]]
    killing = [row for row in rows if not row["branch_chi11_well_defined_on_quotient"]]
    return {
        "result_kind": "SCRATCH_STRICT_ALPHA_AFFINE_SUBGROUP_LATTICE_CHI11_SOURCE_OBSTRUCTION__NO_FALSE_PASS",
        "status": "only-cyclic-and-D12-quotients-admit-chi11-unit5-or-unit7-kill-it",
        "finite_model": {
            "ring": "Z_12",
            "active_count": ACTIVE_COUNT,
            "support_count": len(supports),
            "all_automorphism_units": ALL_UNITS,
            "affine_unit_subgroup_count": len(UNIT_SUBGROUPS),
            "target_q_power": TARGET_Q_POWER,
            "target_eta": TARGET_ETA,
        },
        "subgroup_rows": rows,
        "lattice_summary": {
            "admitting_subgroup_count": len(admitting),
            "killing_subgroup_count": len(killing),
            "admitting_subgroups": [row["subgroup_name"] for row in admitting],
            "killing_subgroups": [row["subgroup_name"] for row in killing],
            "minimal_nontrivial_admitting_reflection_subgroup": "T_semidirect_{1,11}_D12",
            "minimal_killing_unit_additions": ["unit5", "unit7"],
            "full_aut_admits_chi11": False,
        },
        "exact_proof_certificate": exact_proof_certificate(rows),
        "interpretation": {
            "honest_positive": "The subgroup lattice gives an exact finite boundary: D12 is the largest affine unit subgroup preserving the branch chi_11 quotient.",
            "honest_negative": "Any quotient admitting unit 5 or unit 7 mixes A1/A11 with A5/A7 and therefore erases the branch sign.",
            "relation_to_previous_probe": "The previous probe showed full-Aut amplitude can locate the block; this probe proves which subgroup symmetries erase the polarity.",
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
            "Result classifies affine subgroup quotients; it does not supply an internal strict source for selecting the D12/unit-axis premise.",
        ],
    }


def write_markdown(payload: dict[str, Any]) -> str:
    model = payload["finite_model"]
    summary = payload["lattice_summary"]
    lines = [
        "# Affine subgroup lattice chi_11 source obstruction",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Finite model",
        "",
        f"- Ring: `{model['ring']}`",
        f"- Enumerated supports: `{model['support_count']}`",
        f"- Affine unit subgroup count: `{model['affine_unit_subgroup_count']}`",
        f"- All automorphism units: `{model['all_automorphism_units']}`",
        "",
        "## Lattice summary",
        "",
        f"- Admitting subgroup count: `{summary['admitting_subgroup_count']}`",
        f"- Killing subgroup count: `{summary['killing_subgroup_count']}`",
        f"- Admitting subgroups: `{summary['admitting_subgroups']}`",
        f"- Killing subgroups: `{summary['killing_subgroups']}`",
        f"- Minimal nontrivial admitting reflection subgroup: `{summary['minimal_nontrivial_admitting_reflection_subgroup']}`",
        f"- Minimal killing unit additions: `{summary['minimal_killing_unit_additions']}`",
        f"- Full Aut admits chi_11: `{summary['full_aut_admits_chi11']}`",
        "",
        "## Subgroup rows",
        "",
    ]
    for row in payload["subgroup_rows"]:
        lines.append(
            f"- `{row['subgroup_name']}` units `{row['unit_subgroup']}`: "
            f"orbits `{row['orbit_count_on_all_supports']}`, "
            f"branch chi_11 well-defined `{row['branch_chi11_well_defined_on_quotient']}`, "
            f"classification `{row['classification']}`"
        )
    lines.extend(["", "## Proof certificate", ""])
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
