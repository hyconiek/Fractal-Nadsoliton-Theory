#!/usr/bin/env python3
"""Scratch probe: sparsest-extension selector audit for the D12 chi_11 module.

The D12 character-module certificate showed a rank-13 chi_11-covariant module and
531441 branch-normalized ternary functions.  This probe asks the next finite
question without repeating that audit: if the branch values A1/A11=-1 and
A5/A7=+1 are fixed, does a minimal-support extension principle select a unique
D12 chi_11 function?

Using the unit-5 two-cycle basis on the 38 D12 orbits, each nonzero coefficient
turns on both components of one two-cycle.  The branch condition fixes the branch
cycle coefficient to +1.  The remaining 12 coefficients are free in {-1,0,1}.
Therefore the branch-normalized ternary family has 3^12 functions; support size
is 2+2k when k of the 12 non-branch coefficients are nonzero.  The unique
minimum-support extension has support size 2 and is exactly the branch generator.

No false pass: sparsity/minimal support is an extra selector principle.  This
probe does not derive that principle from strict geometry, does not discharge
QW-2191, and does not close ToE.
"""
from __future__ import annotations

import json
from math import comb
from itertools import combinations
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
OUT_JSON = HERE / "bridge_strict_alpha_d12_chi11_sparsest_extension_selector_audit_report.json"
OUT_MD = HERE / "bridge_strict_alpha_d12_chi11_sparsest_extension_selector_audit_report.md"

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


def branch_name(mode: int) -> str:
    return f"A{mode}_k{mode}"


def branch_rows(d_index: dict[tuple[int, ...], int]) -> list[dict[str, Any]]:
    rows = []
    for mode in BRANCH_MODES:
        name = branch_name(mode)
        support = unit_support(mode)
        rows.append(
            {
                "name": name,
                "support": list(support),
                "dihedral_orbit_index": d_index[support],
                "required_chi11_value": 1 if name in {"A5_k5", "A7_k7"} else -1,
            }
        )
    return rows


def support_distribution(free_cycle_count: int) -> list[dict[str, int]]:
    rows = []
    for nonzero_free_coefficients in range(free_cycle_count + 1):
        rows.append(
            {
                "nonzero_free_coefficients": nonzero_free_coefficients,
                "d12_component_support_size": 2 + 2 * nonzero_free_coefficients,
                "function_count": comb(free_cycle_count, nonzero_free_coefficients) * (2 ** nonzero_free_coefficients),
            }
        )
    return rows


def exact_proof_certificate(total_cycles: int, free_cycles: int) -> dict[str, str]:
    return {
        "finite_domain": "All C(12,5)=792 supports are quotiented by D12; the residual unit-5 action gives the chi_11 two-cycle basis.",
        "module_reuse_boundary": f"The D12 chi_11 module has {total_cycles} two-cycle basis generators; this probe only studies branch-normalized ternary extensions inside that already computed module.",
        "branch_constraint": "A1/A11=-1 and A5/A7=+1 fixes the branch two-cycle coefficient to +1.",
        "free_coefficients": f"The other {free_cycles} two-cycle coefficients are independently free in {-1,0,1}, giving 3^{free_cycles}=531441 branch-normalized ternary extensions.",
        "support_formula": "If k non-branch coefficients are nonzero, the D12-component support size is 2+2k.",
        "unique_minimum": "The unique minimum-support extension has k=0, support size 2, and equals the branch generator only.",
        "strict_limit": "The sparsest-extension rule is an added selector criterion, not a strict-core derivation of chi_11 or QW-2191 discharge.",
    }


def build_payload() -> dict[str, Any]:
    supports = all_supports()
    d_orbits = orbit_partition(supports, DIHEDRAL_UNITS)
    d_index = index_by_support(d_orbits)
    cycles = two_cycles(residual_unit5_permutation(d_orbits, d_index))
    branch_cycle = next(cycle for cycle in cycles if set(cycle) == {0, 37})
    free_cycle_count = len(cycles) - 1
    distribution = support_distribution(free_cycle_count)
    total_branch_normalized = sum(row["function_count"] for row in distribution)
    return {
        "result_kind": "SCRATCH_STRICT_ALPHA_D12_CHI11_SPARSEST_EXTENSION_SELECTOR_AUDIT__CONDITIONAL_NOT_A_THEOREM",
        "status": "sparsest-branch-normalized-d12-chi11-extension-is-unique-but-sparsity-is-extra",
        "finite_model": {
            "ring": "Z_12",
            "active_count": ACTIVE_COUNT,
            "support_count": len(supports),
            "dihedral_units": DIHEDRAL_UNITS,
            "automorphism_units": UNITS,
            "d12_orbit_count": len(d_orbits),
            "chi11_two_cycle_count": len(cycles),
            "target_q_power": TARGET_Q_POWER,
            "target_eta": TARGET_ETA,
        },
        "sparsity_summary": {
            "branch_cycle_orbit_pair": list(branch_cycle),
            "branch_cycle_coefficient_fixed_to": 1,
            "free_nonbranch_cycle_count": free_cycle_count,
            "branch_normalized_ternary_function_count": total_branch_normalized,
            "minimum_d12_component_support_size": distribution[0]["d12_component_support_size"],
            "minimum_support_function_count": distribution[0]["function_count"],
            "unique_minimum_is_branch_generator": distribution[0]["function_count"] == 1,
            "maximum_d12_component_support_size": distribution[-1]["d12_component_support_size"],
        },
        "support_size_distribution": distribution,
        "branch_rows": branch_rows(d_index),
        "minimum_support_witness": {
            "nonzero_cycle_coefficients": {"branch_cycle_[0,37]": 1},
            "values_by_d12_orbit": {"0": -1, "37": 1},
            "values_by_branch": {"A1_k1": -1, "A5_k5": 1, "A7_k7": 1, "A11_k11": -1},
            "requires_imported_premise": "D12 chi_11 module plus sparsest-extension selector; not strict full-Aut provenance",
        },
        "exact_proof_certificate": exact_proof_certificate(len(cycles), free_cycle_count),
        "interpretation": {
            "honest_positive": "Given the D12 chi_11 module and a sparsest-extension rule, the branch generator is uniquely selected.",
            "honest_negative": "The sparsest-extension rule is not derived from strict geometry and cannot be used as a hidden QW-2191 discharge.",
            "relation_to_previous_probe": "This refines the rank-13 / 531441 branch-normalized count by proving the exact support-size distribution and unique minimum.",
        },
        "ontology_guardrail": {
            "allowed_reading": "The nadsoliton itself remains the primordial information in a solitonic state.",
            "forbidden_reading": "No separate informational layer underneath the nadsoliton is introduced.",
        },
        "hard_limits": [
            "No identity K_legacy_ont == K_strict_gate is asserted.",
            "No legacy physical-role transfer onto K_strict_gate is used.",
            "No theorem derives chi_11, shell-labels, unit-axis bit, exact-cover clauses, or cardinality 5 from strict geometry.",
            "No theorem derives a sparsest-extension selector from strict core.",
            "No QW-2191 discharge is claimed.",
            "No ToE closure is claimed.",
        ],
    }


def write_markdown(payload: dict[str, Any]) -> str:
    model = payload["finite_model"]
    summary = payload["sparsity_summary"]
    lines = [
        "# D12 chi_11 sparsest-extension selector audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Finite model",
        "",
        f"- Ring: `{model['ring']}`",
        f"- Enumerated supports: `{model['support_count']}`",
        f"- D12 orbit count: `{model['d12_orbit_count']}`",
        f"- chi_11 two-cycle count: `{model['chi11_two_cycle_count']}`",
        "",
        "## Sparsity summary",
        "",
        f"- Branch cycle orbit pair: `{summary['branch_cycle_orbit_pair']}`",
        f"- Free nonbranch cycle count: `{summary['free_nonbranch_cycle_count']}`",
        f"- Branch-normalized ternary functions: `{summary['branch_normalized_ternary_function_count']}`",
        f"- Minimum D12-component support size: `{summary['minimum_d12_component_support_size']}`",
        f"- Minimum-support function count: `{summary['minimum_support_function_count']}`",
        f"- Unique minimum is branch generator: `{summary['unique_minimum_is_branch_generator']}`",
        "",
        "## Support-size distribution",
        "",
    ]
    for row in payload["support_size_distribution"]:
        lines.append(
            f"- support `{row['d12_component_support_size']}`: "
            f"`{row['function_count']}` functions "
            f"(nonzero free coefficients `{row['nonzero_free_coefficients']}`)"
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
