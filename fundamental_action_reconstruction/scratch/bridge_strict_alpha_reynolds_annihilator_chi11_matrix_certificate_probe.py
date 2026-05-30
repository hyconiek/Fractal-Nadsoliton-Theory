#!/usr/bin/env python3
"""Scratch probe: exact Reynolds-annihilator matrix certificate for chi_11.

This is deliberately not another orbit-count audit.  The previous cyclic quotient
probe computed character dimensions by traces.  Here we build the actual integer
matrix numerators for the U12 Reynolds projector and the chi_11 character
projector on the 66-dimensional translation quotient of the C(12,5)=792 support
set.  The proof object is the exact annihilator identity

    (sum_u P_u) * (sum_v chi_11(v) P_v) = 0,

plus a branch-generator witness showing that the A1/A11 versus A5/A7 polarity
lives in the chi_11 sector and is killed by full-Aut Reynolds averaging.

No false pass: this is a matrix-level obstruction/certificate.  It does not
export a strict source for the missing unit-axis or chi_11 character premise, it
does not discharge QW-2191, and it does not close ToE.
"""
from __future__ import annotations

import json
from collections import Counter
from fractions import Fraction
from itertools import combinations
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
OUT_JSON = HERE / "bridge_strict_alpha_reynolds_annihilator_chi11_matrix_certificate_report.json"
OUT_MD = HERE / "bridge_strict_alpha_reynolds_annihilator_chi11_matrix_certificate_report.md"

N = 12
ACTIVE_COUNT = 5
UNITS = [1, 5, 7, 11]
TRANSLATION_UNITS = [1]
CHI11 = {1: 1, 5: -1, 7: -1, 11: 1}
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


def translation_index_by_support(translation_orbits: list[frozenset[tuple[int, ...]]]) -> dict[tuple[int, ...], int]:
    return {support: index for index, orbit in enumerate(translation_orbits) for support in orbit}


def unit_permutations(
    translation_orbits: list[frozenset[tuple[int, ...]]],
    index_by_support: dict[tuple[int, ...], int],
) -> dict[int, list[int]]:
    return {
        unit: [index_by_support[affine_image(min(orbit), 0, unit)] for orbit in translation_orbits]
        for unit in UNITS
    }


def unit_orbits(permutations: dict[int, list[int]]) -> list[list[int]]:
    remaining = set(range(len(permutations[1])))
    rows: list[list[int]] = []
    while remaining:
        seed = min(remaining)
        orbit = sorted({permutations[unit][seed] for unit in UNITS})
        rows.append(orbit)
        remaining -= set(orbit)
    return rows


def permutation_matrix(permutation: list[int]) -> list[list[int]]:
    size = len(permutation)
    matrix = [[0 for _ in range(size)] for _ in range(size)]
    for source, target in enumerate(permutation):
        matrix[target][source] = 1
    return matrix


def zero_matrix(size: int) -> list[list[int]]:
    return [[0 for _ in range(size)] for _ in range(size)]


def add_scaled(target: list[list[int]], matrix: list[list[int]], scale: int) -> None:
    for row_index, row in enumerate(matrix):
        target_row = target[row_index]
        for col_index, value in enumerate(row):
            target_row[col_index] += scale * value


def matmul(left: list[list[int]], right: list[list[int]]) -> list[list[int]]:
    size = len(left)
    product = zero_matrix(size)
    for row in range(size):
        for mid in range(size):
            left_value = left[row][mid]
            if left_value == 0:
                continue
            for col in range(size):
                product[row][col] += left_value * right[mid][col]
    return product


def matvec(matrix: list[list[int]], vector: list[int]) -> list[int]:
    return [sum(value * vector[col] for col, value in enumerate(row)) for row in matrix]


def matrix_rank_integer(matrix: list[list[int]]) -> int:
    rows = [[Fraction(value) for value in row] for row in matrix if any(value != 0 for value in row)]
    if not rows:
        return 0
    row_count = len(rows)
    col_count = len(rows[0])
    rank = 0
    for col in range(col_count):
        pivot = next((row for row in range(rank, row_count) if rows[row][col] != 0), None)
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        pivot_value = rows[rank][col]
        rows[rank] = [value / pivot_value for value in rows[rank]]
        for row in range(row_count):
            if row == rank or rows[row][col] == 0:
                continue
            factor = rows[row][col]
            rows[row] = [value - factor * rows[rank][idx] for idx, value in enumerate(rows[row])]
        rank += 1
        if rank == row_count:
            break
    return rank


def support_label(index: int, translation_orbits: list[frozenset[tuple[int, ...]]]) -> list[int]:
    return list(min(translation_orbits[index]))


def build_projector_numerators(permutations: dict[int, list[int]]) -> tuple[list[list[int]], list[list[int]]]:
    size = len(permutations[1])
    reynolds = zero_matrix(size)
    chi11 = zero_matrix(size)
    for unit in UNITS:
        matrix = permutation_matrix(permutations[unit])
        add_scaled(reynolds, matrix, 1)
        add_scaled(chi11, matrix, CHI11[unit])
    return reynolds, chi11


def branch_generator(index_by_support: dict[tuple[int, ...], int]) -> tuple[list[int], dict[str, Any]]:
    a1_index = index_by_support[unit_support(1)]
    a5_index = index_by_support[unit_support(5)]
    a7_index = index_by_support[unit_support(7)]
    a11_index = index_by_support[unit_support(11)]
    vector = [0 for _ in range(66)]
    vector[a1_index] = -1
    vector[a11_index] = -1
    vector[a5_index] = 1
    vector[a7_index] = 1
    # A1/A11 and A5/A7 are translation-identical pairs on the 66 quotient.
    compact_vector = [0 for _ in range(66)]
    compact_vector[a1_index] = -1
    compact_vector[a5_index] = 1
    return compact_vector, {
        "A1_translation_orbit_index": a1_index,
        "A11_translation_orbit_index": a11_index,
        "A5_translation_orbit_index": a5_index,
        "A7_translation_orbit_index": a7_index,
        "quotient_pair_indices": sorted({a1_index, a5_index}),
        "normalized_values_by_branch": {"A1_k1": -1, "A11_k11": -1, "A5_k5": 1, "A7_k7": 1},
        "compact_translation_vector_nonzero_entries": {str(a1_index): -1, str(a5_index): 1},
        "translation_pair_collapse_note": "A1 and A11 share one translation orbit; A5 and A7 share one translation orbit.",
    }


def basis_rows_for_chi11(
    full_unit_orbits: list[list[int]],
    permutations: dict[int, list[int]],
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for orbit in full_unit_orbits:
        seed = orbit[0]
        stabilizer = [unit for unit in UNITS if permutations[unit][seed] == seed]
        if not all(CHI11[unit] == 1 for unit in stabilizer):
            continue
        values: dict[int, int] = {}
        for unit in UNITS:
            values[permutations[unit][seed]] = CHI11[unit]
        rows.append(
            {
                "basis_index": len(rows),
                "translation_orbit_indices": orbit,
                "stabilizer_units": stabilizer,
                "values_by_translation_orbit_index": {str(index): values[index] for index in sorted(values)},
                "unit_orbit_signed_sum": sum(values.values()),
            }
        )
    return rows


def nonzero_entries(vector: list[int]) -> dict[str, int]:
    return {str(index): value for index, value in enumerate(vector) if value != 0}


def build_payload() -> dict[str, Any]:
    supports = all_supports()
    translation_orbits = orbit_partition(supports, TRANSLATION_UNITS)
    index_by_support = translation_index_by_support(translation_orbits)
    permutations = unit_permutations(translation_orbits, index_by_support)
    full_unit_orbits = unit_orbits(permutations)
    reynolds_num, chi11_num = build_projector_numerators(permutations)
    annihilator_product = matmul(reynolds_num, chi11_num)
    reverse_product = matmul(chi11_num, reynolds_num)
    branch_vec, branch_info = branch_generator(index_by_support)
    reynolds_branch = matvec(reynolds_num, branch_vec)
    chi11_branch = matvec(chi11_num, branch_vec)
    basis_rows = basis_rows_for_chi11(full_unit_orbits, permutations)
    unit_orbit_sums = []
    for orbit in full_unit_orbits:
        unit_orbit_sums.append({"translation_orbit_indices": orbit, "branch_sum": sum(branch_vec[index] for index in orbit)})

    annihilator_zero = all(value == 0 for row in annihilator_product for value in row)
    reverse_zero = all(value == 0 for row in reverse_product for value in row)
    branch_reynolds_zero = all(value == 0 for value in reynolds_branch)
    branch_chi11_eigen = chi11_branch == [4 * value for value in branch_vec]
    all_chi11_basis_rows_reynolds_zero = all(row["unit_orbit_signed_sum"] == 0 for row in basis_rows)

    return {
        "result_kind": "SCRATCH_STRICT_ALPHA_REYNOLDS_ANNIHILATOR_CHI11_MATRIX_CERTIFICATE__NO_FALSE_PASS",
        "status": "exact-integer-matrix-annihilator-computed-on-translation-quotient",
        "finite_model": {
            "ring": "Z_12",
            "active_count": ACTIVE_COUNT,
            "support_count": len(supports),
            "translation_orbit_count": len(translation_orbits),
            "residual_unit_group": UNITS,
            "residual_unit_orbit_count": len(full_unit_orbits),
            "residual_unit_orbit_size_counts": dict(sorted(Counter(len(orbit) for orbit in full_unit_orbits).items())),
            "matrix_size": len(translation_orbits),
            "projector_denominator": len(UNITS),
            "target_q_power": TARGET_Q_POWER,
            "target_eta": TARGET_ETA,
        },
        "matrix_certificate": {
            "reynolds_numerator_rank": matrix_rank_integer(reynolds_num),
            "chi11_numerator_rank": matrix_rank_integer(chi11_num),
            "reynolds_times_chi11_numerator_is_zero_matrix": annihilator_zero,
            "chi11_times_reynolds_numerator_is_zero_matrix": reverse_zero,
            "annihilator_zero_entry_count": sum(1 for row in annihilator_product for value in row if value == 0),
            "annihilator_total_entry_count": len(annihilator_product) ** 2,
            "identity": "R_num*C_chi11_num = C_chi11_num*R_num = 0 over integers before dividing by 4.",
        },
        "branch_generator_certificate": {
            **branch_info,
            "A1_representative_support": support_label(branch_info["A1_translation_orbit_index"], translation_orbits),
            "A5_representative_support": support_label(branch_info["A5_translation_orbit_index"], translation_orbits),
            "reynolds_numerator_on_branch_is_zero": branch_reynolds_zero,
            "chi11_numerator_on_branch_equals_four_times_branch": branch_chi11_eigen,
            "reynolds_numerator_image_nonzero_entries": nonzero_entries(reynolds_branch),
            "chi11_numerator_image_nonzero_entries": nonzero_entries(chi11_branch),
            "full_unit_orbit_branch_sums": unit_orbit_sums,
        },
        "chi11_basis_annihilation_summary": {
            "chi11_basis_row_count": len(basis_rows),
            "all_chi11_basis_rows_have_zero_reynolds_sum": all_chi11_basis_rows_reynolds_zero,
            "basis_rows": basis_rows,
        },
        "exact_proof_certificate": {
            "finite_domain": "All C(12,5)=792 supports are first quotiented by translations into a 66-point exact finite representation.",
            "matrix_construction": "For each unit u in U12, P_u is the exact 66x66 permutation matrix on translation orbits; R_num=sum_u P_u and C_chi11_num=sum_u chi_11(u)P_u.",
            "annihilator_identity": "The computed integer products R_num*C_chi11_num and C_chi11_num*R_num are zero matrices, so full-Aut Reynolds data annihilates the chi_11 sector exactly.",
            "branch_witness": "The compact A1/A11=-1 and A5/A7=+1 translation vector is a chi_11 eigenvector with C_chi11_num*v=4v and R_num*v=0.",
            "logical_boundary": "This strengthens the obstruction to full-Aut export of chi_11 polarity, but the chi_11/unit-axis choice remains an extra premise, not a strict-core theorem.",
        },
        "interpretation": {
            "honest_positive": "The previous trace-rank statement is upgraded to an explicit integer matrix annihilator and branch-vector eigen-certificate.",
            "honest_negative": "Any full-Aut Reynolds/invariant support datum loses the branch polarity exactly; no full-Aut invariant source exports chi_11.",
            "relation_to_previous_probe": "This does not repeat subgroup or premise enumeration; it verifies the projector-product obstruction on the 66-dimensional quotient representation.",
        },
        "ontology_guardrail": {
            "allowed_reading": "The nadsoliton itself remains the primordial information in a solitonic state.",
            "forbidden_reading": "No separate informational layer underneath the nadsoliton is introduced.",
        },
        "hard_limits": [
            "No identity K_legacy_ont == K_strict_gate is asserted.",
            "No legacy physical-role transfer onto K_strict_gate is used.",
            "No theorem derives chi_11, shell-labels, unit-axis bit, exact-cover clauses, or cardinality 5 from strict geometry.",
            "No theorem derives a full-Aut internal chi_11 polarity source.",
            "No QW-2191 discharge is claimed.",
            "No ToE closure is claimed.",
        ],
    }


def write_markdown(payload: dict[str, Any]) -> str:
    model = payload["finite_model"]
    matrix = payload["matrix_certificate"]
    branch = payload["branch_generator_certificate"]
    basis = payload["chi11_basis_annihilation_summary"]
    lines = [
        "# Reynolds annihilator chi_11 matrix certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Finite model",
        "",
        f"- Ring: `{model['ring']}`",
        f"- Enumerated supports: `{model['support_count']}`",
        f"- Translation orbit count / matrix size: `{model['translation_orbit_count']}` / `{model['matrix_size']}`",
        f"- Residual unit orbit count: `{model['residual_unit_orbit_count']}`",
        f"- Residual unit orbit size counts: `{model['residual_unit_orbit_size_counts']}`",
        "",
        "## Exact matrix certificate",
        "",
        f"- Reynolds numerator rank: `{matrix['reynolds_numerator_rank']}`",
        f"- chi_11 numerator rank: `{matrix['chi11_numerator_rank']}`",
        f"- R_num*C_chi11_num zero: `{matrix['reynolds_times_chi11_numerator_is_zero_matrix']}`",
        f"- C_chi11_num*R_num zero: `{matrix['chi11_times_reynolds_numerator_is_zero_matrix']}`",
        f"- Zero entries in annihilator product: `{matrix['annihilator_zero_entry_count']}/{matrix['annihilator_total_entry_count']}`",
        "",
        "## Branch generator witness",
        "",
        f"- Quotient pair indices: `{branch['quotient_pair_indices']}`",
        f"- Normalized branch values: `{branch['normalized_values_by_branch']}`",
        f"- A1 representative support: `{branch['A1_representative_support']}`",
        f"- A5 representative support: `{branch['A5_representative_support']}`",
        f"- Reynolds numerator kills branch: `{branch['reynolds_numerator_on_branch_is_zero']}`",
        f"- chi_11 numerator returns 4*branch: `{branch['chi11_numerator_on_branch_equals_four_times_branch']}`",
        "",
        "## chi_11 basis annihilation",
        "",
        f"- chi_11 basis rows: `{basis['chi11_basis_row_count']}`",
        f"- Every row has zero Reynolds sum: `{basis['all_chi11_basis_rows_have_zero_reynolds_sum']}`",
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
