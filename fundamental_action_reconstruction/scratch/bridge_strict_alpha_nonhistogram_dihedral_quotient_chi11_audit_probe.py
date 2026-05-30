#!/usr/bin/env python3
"""Scratch probe: non-histogram dihedral-quotient audit for chi_11 sources.

The universal histogram no-go left open non-histogram support features.  This
probe makes the next finite source audit precise: enumerate all C(12,5)=792
supports, quotient them by the dihedral group D12 (translations plus reversal),
and then compare that quotient with the larger affine/full-Aut action
x -> a + u*x, u in {1,5,7,11}.

Result: the A1/A11 contiguous branch and the A5/A7 d5 branch are distinct on the
D12 quotient, but they are the two D12 components of one full-Aut affine orbit,
swapped by unit 5 (and unit 7).  Therefore D12-invariant non-histogram records
can carry the chi_11 sign only after the full Aut symmetry has already been
reduced to a cyclic/dihedral order.  Full-Aut invariant support-only records,
however non-histogram, still cannot distinguish A1 from A5.

No false pass: this is a support-combinatorics quotient audit, not a strict
nadsoliton provenance theorem.  It does not derive chi_11, does not discharge
QW-2191, and does not close ToE.
"""
from __future__ import annotations

import json
from collections import Counter
from itertools import combinations
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
OUT_JSON = HERE / "bridge_strict_alpha_nonhistogram_dihedral_quotient_chi11_audit_report.json"
OUT_MD = HERE / "bridge_strict_alpha_nonhistogram_dihedral_quotient_chi11_audit_report.md"

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
    return tuple(sorted(((shift + unit * value) % N) for value in support))


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


def branch_name(mode: int) -> str:
    return f"A{mode}_k{mode}"


def branch_pair(name: str) -> str:
    if name in {"A1_k1", "A11_k11"}:
        return "contiguous_pair_A1_A11"
    if name in {"A5_k5", "A7_k7"}:
        return "d5_pair_A5_A7"
    raise ValueError(name)


def required_chi11_value(name: str) -> int:
    return 1 if branch_pair(name) == "d5_pair_A5_A7" else -1


def index_by_support(orbits: list[frozenset[tuple[int, ...]]]) -> dict[tuple[int, ...], int]:
    return {support: index for index, orbit in enumerate(orbits) for support in orbit}


def d_orbit_action_rows(d_orbits: list[frozenset[tuple[int, ...]]], d_index: dict[tuple[int, ...], int]) -> list[dict[str, Any]]:
    rows = []
    for index, orbit in enumerate(d_orbits):
        representative = min(orbit)
        image_by_unit = {str(unit): d_index[affine_image(representative, 0, unit)] for unit in UNITS}
        rows.append(
            {
                "dihedral_orbit_index": index,
                "representative_support": list(representative),
                "gap_necklace": list(gap_necklace(representative)),
                "orbit_size": len(orbit),
                "unit_image_dihedral_orbit_index": image_by_unit,
                "fixed_by_unit5_mod_dihedral": image_by_unit["5"] == index,
            }
        )
    return rows


def full_orbit_decomposition_rows(
    full_orbits: list[frozenset[tuple[int, ...]]],
    d_orbits: list[frozenset[tuple[int, ...]]],
    d_index: dict[tuple[int, ...], int],
) -> list[dict[str, Any]]:
    rows = []
    branch_support_set = {unit_support(mode) for mode in BRANCH_MODES}
    for full_index, full_orbit in enumerate(full_orbits):
        d_indices = sorted({d_index[support] for support in full_orbit})
        rows.append(
            {
                "full_affine_orbit_index": full_index,
                "full_orbit_size": len(full_orbit),
                "dihedral_component_indices": d_indices,
                "dihedral_component_count": len(d_indices),
                "dihedral_component_gap_necklaces": [list(gap_necklace(min(d_orbits[index]))) for index in d_indices],
                "contains_branch_orbit": bool(full_orbit & branch_support_set),
                "is_chi11_capacity_block": len(d_indices) == 2,
            }
        )
    return rows


def branch_rows(d_index: dict[tuple[int, ...], int], full_index: dict[tuple[int, ...], int]) -> list[dict[str, Any]]:
    rows = []
    for mode in BRANCH_MODES:
        name = branch_name(mode)
        support = unit_support(mode)
        rows.append(
            {
                "name": name,
                "mode": mode,
                "support": list(support),
                "gap_necklace": list(gap_necklace(support)),
                "branch_pair": branch_pair(name),
                "required_chi11_value": required_chi11_value(name),
                "dihedral_orbit_index": d_index[support],
                "full_affine_orbit_index": full_index[support],
            }
        )
    return rows


def exact_proof_certificate(d_count: int, full_count: int, paired_count: int, fixed_count: int) -> dict[str, str]:
    return {
        "finite_domain": "All C(12,5)=792 supports are enumerated exactly.",
        "quotient_split": f"The D12 quotient has {d_count} orbits; the full affine Aut quotient has {full_count} orbits.",
        "branch_block": "A1/A11 and A5/A7 are two distinct D12 orbits but one full-Aut affine orbit; unit 5 swaps them.",
        "full_aut_no_go": "Any full-Aut invariant support-only source is constant on the full affine orbit, so it has F(A1)=F(A5) and cannot export nonzero chi_11.",
        "dihedral_positive_boundary": "A D12-invariant sign can be assigned as - on the A1/A11 component and + on the A5/A7 component, matching chi_11 on the branch block.",
        "import_diagnosis": "That D12 sign is not strict full-Aut data: it requires a reduced cyclic/dihedral order, equivalently a unit-axis premise separating {±1} from {±5}.",
        "global_character_space": f"Under the unit-5 action on D12 orbits there are {paired_count} two-component blocks and {fixed_count} fixed blocks; integer chi_11 anti-invariants have dimension {paired_count} and must vanish on fixed blocks.",
        "scope_limit": "This audits support-combinatorics features, not all possible strict nadsoliton provenance mechanisms.",
    }


def build_payload() -> dict[str, Any]:
    supports = all_supports()
    d_orbits = orbit_partition(supports, DIHEDRAL_UNITS)
    full_orbits = orbit_partition(supports, UNITS)
    d_index = index_by_support(d_orbits)
    full_index = index_by_support(full_orbits)
    action_rows = d_orbit_action_rows(d_orbits, d_index)
    decomposition_rows = full_orbit_decomposition_rows(full_orbits, d_orbits, d_index)
    branch = branch_rows(d_index, full_index)
    d_size_counts = Counter(len(orbit) for orbit in d_orbits)
    full_size_counts = Counter(len(orbit) for orbit in full_orbits)
    paired_blocks = [row for row in decomposition_rows if row["dihedral_component_count"] == 2]
    fixed_blocks = [row for row in decomposition_rows if row["dihedral_component_count"] == 1]
    branch_full_rows = [row for row in decomposition_rows if row["contains_branch_orbit"]]
    assert len(branch_full_rows) == 1
    return {
        "result_kind": "SCRATCH_STRICT_ALPHA_NONHISTOGRAM_DIHEDRAL_QUOTIENT_CHI11_AUDIT__BOUNDARY_NOT_A_THEOREM",
        "status": "dihedral-nonhistogram-features-can-carry-chi11-only-with-reduced-symmetry-premise",
        "finite_model": {
            "ring": "Z_12",
            "active_count": ACTIVE_COUNT,
            "support_count": len(supports),
            "automorphism_units": UNITS,
            "dihedral_units": DIHEDRAL_UNITS,
            "dihedral_orbit_count": len(d_orbits),
            "full_affine_aut_orbit_count": len(full_orbits),
            "dihedral_orbit_size_counts": dict(sorted(d_size_counts.items())),
            "full_affine_orbit_size_counts": dict(sorted(full_size_counts.items())),
            "target_q_power": TARGET_Q_POWER,
            "target_eta": TARGET_ETA,
        },
        "branch_orbit_rows": branch,
        "branch_full_orbit_certificate": {
            "full_affine_orbit_index": branch_full_rows[0]["full_affine_orbit_index"],
            "dihedral_component_indices": branch_full_rows[0]["dihedral_component_indices"],
            "component_gap_necklaces": branch_full_rows[0]["dihedral_component_gap_necklaces"],
            "unit5_maps_A1_dihedral_component_to_A5_component": d_index[affine_image(unit_support(1), 0, 5)] == d_index[unit_support(5)],
            "unit7_maps_A1_dihedral_component_to_A7_component": d_index[affine_image(unit_support(1), 0, 7)] == d_index[unit_support(7)],
            "unit11_maps_A1_to_contiguous_reversal_component": d_index[affine_image(unit_support(1), 0, 11)] == d_index[unit_support(11)],
            "full_aut_invariant_singleton_A5_not_A1_classifier_count": 0,
        },
        "quotient_summary": {
            "dihedral_boolean_classifier_count_power": len(d_orbits),
            "dihedral_boolean_classifier_total_exact": str(2 ** len(d_orbits)),
            "dihedral_branch_chi11_classifier_count_power": len(d_orbits) - 2,
            "dihedral_branch_chi11_classifier_total_exact": str(2 ** (len(d_orbits) - 2)),
            "full_aut_boolean_classifier_count_power": len(full_orbits),
            "full_aut_boolean_classifier_total_exact": str(2 ** len(full_orbits)),
            "full_aut_singleton_A5_not_A1_classifier_count": 0,
            "full_aut_blocks_with_two_dihedral_components": len(paired_blocks),
            "full_aut_blocks_fixed_on_dihedral_quotient": len(fixed_blocks),
            "global_integer_chi11_antiinvariant_dimension_on_dihedral_quotient": len(paired_blocks),
            "global_nonzero_pm1_chi11_character_count": 0,
        },
        "branch_compatible_dihedral_sign_witness": {
            "description": "Set the contiguous D12 component to -1 and the d5 D12 component to +1; extend arbitrarily or by zero off selected two-component blocks.",
            "values_by_branch": {row["name"]: row["required_chi11_value"] for row in branch},
            "requires_imported_premise": "reduced cyclic/dihedral order or unit-axis bit; not full-Aut strict-core data",
        },
        "full_affine_orbit_decomposition_rows": decomposition_rows,
        "dihedral_orbit_action_rows": action_rows,
        "exact_proof_certificate": exact_proof_certificate(len(d_orbits), len(full_orbits), len(paired_blocks), len(fixed_blocks)),
        "interpretation": {
            "honest_positive": "The finite D12 quotient contains a branch-level chi_11 carrier: the split between the contiguous and d5 D12 components.",
            "honest_negative": "That carrier is not full-Aut provenance; it appears only after reducing symmetry from full affine Aut to D12.",
            "relation_to_previous_probe": "This goes beyond histogram-only data, but preserves the same no-false-pass diagnosis: distinguishing A5 from A1 imports a missing orientation/unit-axis premise.",
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
            "Result is a finite support-combinatorics boundary audit, not a complete strict-source no-go theorem.",
        ],
    }


def write_markdown(payload: dict[str, Any]) -> str:
    model = payload["finite_model"]
    summary = payload["quotient_summary"]
    cert = payload["branch_full_orbit_certificate"]
    lines = [
        "# Non-histogram dihedral-quotient chi_11 audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Finite model",
        "",
        f"- Ring: `{model['ring']}`",
        f"- Active support count: `{model['active_count']}`",
        f"- Enumerated supports: `{model['support_count']}`",
        f"- D12 orbit count: `{model['dihedral_orbit_count']}`",
        f"- Full affine Aut orbit count: `{model['full_affine_aut_orbit_count']}`",
        f"- Full-Aut two-D12-component blocks: `{summary['full_aut_blocks_with_two_dihedral_components']}`",
        f"- Full-Aut fixed D12 blocks: `{summary['full_aut_blocks_fixed_on_dihedral_quotient']}`",
        "",
        "## Branch certificate",
        "",
        f"- Branch full-Aut orbit index: `{cert['full_affine_orbit_index']}`",
        f"- D12 component indices: `{cert['dihedral_component_indices']}`",
        f"- Component gap necklaces: `{cert['component_gap_necklaces']}`",
        f"- Unit 5 maps A1 component to A5 component: `{cert['unit5_maps_A1_dihedral_component_to_A5_component']}`",
        f"- Full-Aut invariant singleton A5-not-A1 classifier count: `{cert['full_aut_invariant_singleton_A5_not_A1_classifier_count']}`",
        "",
        "## Classifier counts",
        "",
        f"- D12 Boolean classifiers: `2^{summary['dihedral_boolean_classifier_count_power']} = {summary['dihedral_boolean_classifier_total_exact']}`",
        f"- D12 branch-chi_11-compatible classifiers: `2^{summary['dihedral_branch_chi11_classifier_count_power']} = {summary['dihedral_branch_chi11_classifier_total_exact']}`",
        f"- Full-Aut Boolean classifiers: `2^{summary['full_aut_boolean_classifier_count_power']} = {summary['full_aut_boolean_classifier_total_exact']}`",
        f"- Full-Aut singleton A5-not-A1 classifiers: `{summary['full_aut_singleton_A5_not_A1_classifier_count']}`",
        f"- Global integer chi_11 anti-invariant dimension on D12 quotient: `{summary['global_integer_chi11_antiinvariant_dimension_on_dihedral_quotient']}`",
        f"- Global nonzero ±1 chi_11 character count: `{summary['global_nonzero_pm1_chi11_character_count']}`",
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
