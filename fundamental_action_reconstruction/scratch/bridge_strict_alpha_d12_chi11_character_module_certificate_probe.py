#!/usr/bin/env python3
"""Scratch probe: exact D12-quotient chi_11 character module certificate.

The previous non-histogram quotient audit showed that A1/A11 and A5/A7 are two
D12 components inside one full-Aut support orbit.  This probe turns that boundary
into a small algebra certificate: enumerate the D12 quotient of all C(12,5)=792
supports, compute the residual unit-5 involution on the 38 D12 orbits, and solve
exactly for every D12-invariant support function f that transforms as chi_11:

    f(unit5 * O) = -f(O), while unit11 is already inside D12.

Over characteristic zero, each two-cycle contributes one free integer generator
and each fixed orbit is forced to value 0.  Therefore the full chi_11-covariant
D12-quotient module has rank 13, with no globally nonzero ±1 character because
12 fixed D12 orbits must vanish.  Full-Aut invariance is the opposite equation
f(unit5 * O)=f(O), so its intersection with chi_11 covariance is zero.

No false pass: the branch block has a valid D12 quotient character generator,
but using it still requires the reduced D12 quotient / unit-axis premise.  This
is not a strict provenance theorem for chi_11, does not discharge QW-2191, and
does not close ToE.
"""
from __future__ import annotations

import json
from itertools import combinations
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
OUT_JSON = HERE / "bridge_strict_alpha_d12_chi11_character_module_certificate_report.json"
OUT_MD = HERE / "bridge_strict_alpha_d12_chi11_character_module_certificate_report.md"

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


def required_chi11_value(name: str) -> int:
    return 1 if name in {"A5_k5", "A7_k7"} else -1


def residual_unit5_permutation(d_orbits: list[frozenset[tuple[int, ...]]], d_index: dict[tuple[int, ...], int]) -> list[int]:
    return [d_index[affine_image(min(orbit), 0, 5)] for orbit in d_orbits]


def cycle_decomposition(unit5_perm: list[int]) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    visited: set[int] = set()
    two_cycles: list[dict[str, Any]] = []
    fixed: list[dict[str, Any]] = []
    for index, image in enumerate(unit5_perm):
        if index in visited:
            continue
        if image == index:
            fixed.append({"orbit_index": index, "forced_chi11_value": 0})
            visited.add(index)
            continue
        low, high = sorted((index, image))
        two_cycles.append(
            {
                "cycle_index": len(two_cycles),
                "orbit_pair": [low, high],
                "basis_values_low_high": [-1, 1],
                "generator_rule": f"e_{low}_{high}({low})=-1, e_{low}_{high}({high})=+1, zero elsewhere",
            }
        )
        visited.update({low, high})
    return two_cycles, fixed


def enrich_cycle_rows(two_cycles: list[dict[str, Any]], d_orbits: list[frozenset[tuple[int, ...]]]) -> list[dict[str, Any]]:
    branch_indices = {0, 37}
    rows = []
    for row in two_cycles:
        low, high = row["orbit_pair"]
        rows.append(
            {
                **row,
                "low_representative_support": list(min(d_orbits[low])),
                "high_representative_support": list(min(d_orbits[high])),
                "low_gap_necklace": list(gap_necklace(min(d_orbits[low]))),
                "high_gap_necklace": list(gap_necklace(min(d_orbits[high]))),
                "is_branch_chi11_generator": set(row["orbit_pair"]) == branch_indices,
            }
        )
    return rows


def branch_rows(d_index: dict[tuple[int, ...], int]) -> list[dict[str, Any]]:
    rows = []
    for mode in BRANCH_MODES:
        name = branch_name(mode)
        support = unit_support(mode)
        orbit_index = d_index[support]
        rows.append(
            {
                "name": name,
                "support": list(support),
                "dihedral_orbit_index": orbit_index,
                "gap_necklace": list(gap_necklace(support)),
                "required_chi11_value": required_chi11_value(name),
            }
        )
    return rows


def exact_proof_certificate(two_cycle_count: int, fixed_count: int) -> dict[str, str]:
    return {
        "finite_domain": "All C(12,5)=792 supports are quotient-enumerated by D12.",
        "residual_action": "Modulo D12, units 1 and 11 are trivial while units 5 and 7 induce the same involution on D12 orbits.",
        "module_equations": "A D12-invariant chi_11-covariant support function is exactly a function f on D12 orbits satisfying f(tau(O))=-f(O), tau=unit5.",
        "rank_certificate": f"The tau action has {two_cycle_count} two-cycles and {fixed_count} fixed points; hence the integer chi_11 module has rank {two_cycle_count} and fixed-point values are forced to 0.",
        "full_aut_intersection": "Full-Aut invariance requires f(tau(O))=f(O); together with chi_11 covariance this gives f(O)=0 on every orbit over characteristic zero.",
        "branch_generator": "The A1/A11 orbit and A5/A7 orbit form one tau two-cycle, so there is a branch-local generator with values -1 and +1 on that pair.",
        "not_strict_source": "Selecting that generator requires the D12 quotient and unit-axis orientation; it is not exported by full-Aut strict support data alone.",
    }


def build_payload() -> dict[str, Any]:
    supports = all_supports()
    d_orbits = orbit_partition(supports, DIHEDRAL_UNITS)
    d_index = index_by_support(d_orbits)
    unit5_perm = residual_unit5_permutation(d_orbits, d_index)
    two_cycles, fixed = cycle_decomposition(unit5_perm)
    two_cycle_rows = enrich_cycle_rows(two_cycles, d_orbits)
    branches = branch_rows(d_index)
    branch_cycle = next(row for row in two_cycle_rows if row["is_branch_chi11_generator"])
    return {
        "result_kind": "SCRATCH_STRICT_ALPHA_D12_CHI11_CHARACTER_MODULE_CERTIFICATE__BOUNDARY_NOT_A_THEOREM",
        "status": "exact-d12-quotient-chi11-module-rank-13-but-requires-unit-axis-premise",
        "finite_model": {
            "ring": "Z_12",
            "active_count": ACTIVE_COUNT,
            "support_count": len(supports),
            "dihedral_units": DIHEDRAL_UNITS,
            "automorphism_units": UNITS,
            "d12_orbit_count": len(d_orbits),
            "residual_unit5_permutation_size": len(unit5_perm),
            "target_q_power": TARGET_Q_POWER,
            "target_eta": TARGET_ETA,
        },
        "module_summary": {
            "tau_two_cycle_count": len(two_cycles),
            "tau_fixed_orbit_count": len(fixed),
            "integer_chi11_covariant_module_rank": len(two_cycles),
            "full_aut_invariant_intersection_rank": 0,
            "global_nonzero_pm1_character_count": 0,
            "ternary_minus_zero_plus_chi11_covariant_count": 3 ** len(two_cycles),
            "branch_normalized_ternary_chi11_covariant_count": 3 ** (len(two_cycles) - 1),
            "fixed_orbit_forced_value": 0,
        },
        "branch_rows": branches,
        "branch_generator_certificate": {
            "cycle_index": branch_cycle["cycle_index"],
            "orbit_pair": branch_cycle["orbit_pair"],
            "basis_values_low_high": branch_cycle["basis_values_low_high"],
            "low_gap_necklace": branch_cycle["low_gap_necklace"],
            "high_gap_necklace": branch_cycle["high_gap_necklace"],
            "values_by_branch": {row["name"]: row["required_chi11_value"] for row in branches},
            "requires_imported_premise": "unit-axis orientation selecting the D12 quotient component pair; not full-Aut strict-core provenance",
        },
        "two_cycle_basis_rows": two_cycle_rows,
        "fixed_orbit_rows": fixed,
        "exact_proof_certificate": exact_proof_certificate(len(two_cycles), len(fixed)),
        "interpretation": {
            "honest_positive": "The complete D12-quotient chi_11-covariant module is explicitly computed and has a branch generator.",
            "honest_negative": "The generator is a reduced-symmetry object; full-Aut invariant support data intersects this module only at zero.",
            "relation_to_previous_probe": "The previous quotient audit identified the boundary; this probe solves the full character module on that quotient.",
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
            "Result is a complete D12-quotient character-module certificate, not a complete strict-source provenance theorem.",
        ],
    }


def write_markdown(payload: dict[str, Any]) -> str:
    model = payload["finite_model"]
    summary = payload["module_summary"]
    branch = payload["branch_generator_certificate"]
    lines = [
        "# D12 chi_11 character module certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Finite model",
        "",
        f"- Ring: `{model['ring']}`",
        f"- Enumerated supports: `{model['support_count']}`",
        f"- D12 orbit count: `{model['d12_orbit_count']}`",
        f"- Residual unit-5 permutation size: `{model['residual_unit5_permutation_size']}`",
        "",
        "## Module summary",
        "",
        f"- Tau two-cycles: `{summary['tau_two_cycle_count']}`",
        f"- Tau fixed D12 orbits: `{summary['tau_fixed_orbit_count']}`",
        f"- Integer chi_11 module rank: `{summary['integer_chi11_covariant_module_rank']}`",
        f"- Full-Aut invariant intersection rank: `{summary['full_aut_invariant_intersection_rank']}`",
        f"- Global nonzero ±1 character count: `{summary['global_nonzero_pm1_character_count']}`",
        f"- Ternary chi_11-covariant count: `{summary['ternary_minus_zero_plus_chi11_covariant_count']}`",
        f"- Branch-normalized ternary count: `{summary['branch_normalized_ternary_chi11_covariant_count']}`",
        "",
        "## Branch generator",
        "",
        f"- Cycle index: `{branch['cycle_index']}`",
        f"- Orbit pair: `{branch['orbit_pair']}`",
        f"- Basis values low/high: `{branch['basis_values_low_high']}`",
        f"- Low/high gap necklaces: `{branch['low_gap_necklace']}` / `{branch['high_gap_necklace']}`",
        f"- Requires premise: {branch['requires_imported_premise']}",
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
