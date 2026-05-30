#!/usr/bin/env python3
"""Scratch probe: candidate-source transform audit for the missing chi_11 bit.

The previous audits established a conditional positive result: once the
chi_11/shell-label premise is supplied, d1+d6 is the unique and locally necessary
pure shell-exclusion selector for A5/d5.  This probe asks the next source
question in a finite, machine-checkable way: do common branch-level candidate
records already export the required chi_11 bit, or do they only do so after a
shell-label / non-full-Aut coordinate has been imported?

We audit the four unit-generated branch supports A1,A5,A7,A11.  The required
selector bit is constant on reversal pairs {A1,A11} and {A5,A7}, and flips under
units 5 and 7; equivalently it transforms as the character with kernel {1,11}.

Result:
- full-Aut invariant records (orbit-symmetrized distance histogram, support
  size, antipodal count) cannot distinguish the two branch pairs;
- raw shell-labelled records such as d5_count-d1_count do transform like chi_11,
  but only because the d1-vs-d5 shell coordinate has already been named;
- dihedral-only records such as gap necklace distinguish A1 from A5 but are not
  full-Aut invariant and therefore also require a reduced-symmetry premise.

No false pass: this is a finite source-candidate classification.  It does not
derive chi_11 from strict nadsoliton geometry, does not discharge QW-2191, and
does not close ToE.
"""
from __future__ import annotations

import json
from itertools import combinations
from pathlib import Path
from typing import Any, Callable

HERE = Path(__file__).resolve().parent
OUT_JSON = HERE / "bridge_strict_alpha_chi11_candidate_source_transform_audit_report.json"
OUT_MD = HERE / "bridge_strict_alpha_chi11_candidate_source_transform_audit_report.md"

N = 12
ACTIVE_COUNT = 5
UNITS = [1, 5, 7, 11]
BRANCH_MODES = [1, 5, 7, 11]
TARGET_Q_POWER = "256/243"
TARGET_ETA = "9/5"
CONTIGUOUS_PAIR = "contiguous_pair_A1_A11"
D5_PAIR = "d5_pair_A5_A7"


def folded(value: int) -> int:
    residue = value % N
    return min(residue, N - residue)


def unit_support(mode: int) -> tuple[int, ...]:
    return tuple(sorted((mode * index) % N for index in range(ACTIVE_COUNT)))


def multiply_support(unit: int, support: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(sorted((unit * value) % N for value in support))


def branch_name(mode: int) -> str:
    return f"A{mode}_k{mode}"


def branch_pair(name: str) -> str:
    if name in {"A1_k1", "A11_k11"}:
        return CONTIGUOUS_PAIR
    if name in {"A5_k5", "A7_k7"}:
        return D5_PAIR
    raise ValueError(name)


def required_chi11_value(name: str) -> int:
    return 1 if branch_pair(name) == D5_PAIR else -1


def distance_histogram(support: tuple[int, ...]) -> tuple[int, int, int, int, int, int]:
    counts = [0] * (N // 2)
    for left, right in combinations(support, 2):
        counts[folded(right - left) - 1] += 1
    return tuple(counts)  # type: ignore[return-value]


def cyclic_gap_tuple(support: tuple[int, ...]) -> tuple[int, ...]:
    cyclic = sorted(support)
    gaps = [cyclic[index + 1] - cyclic[index] for index in range(len(cyclic) - 1)]
    gaps.append(N + cyclic[0] - cyclic[-1])
    return tuple(gaps)


def canonical_necklace(gaps: tuple[int, ...]) -> tuple[int, ...]:
    rotations = [gaps[index:] + gaps[:index] for index in range(len(gaps))]
    reversals = [tuple(reversed(rotation)) for rotation in rotations]
    return min(rotations + reversals)


def full_aut_orbit_histogram(histogram: tuple[int, int, int, int, int, int]) -> tuple[int, int, int, int, int]:
    # Full Aut shell orbits are {d1,d5}, {d2}, {d3}, {d4}, {d6}.
    return (histogram[0] + histogram[4], histogram[1], histogram[2], histogram[3], histogram[5])


def shell_labelled_d5_minus_d1(histogram: tuple[int, int, int, int, int, int]) -> int:
    return histogram[4] - histogram[0]


def shell_labelled_d1_plus_d6_energy(histogram: tuple[int, int, int, int, int, int]) -> int:
    return histogram[0] + histogram[5]


def shell_labelled_d5_plus_d6_energy(histogram: tuple[int, int, int, int, int, int]) -> int:
    return histogram[4] + histogram[5]


def branch_rows() -> list[dict[str, Any]]:
    rows = []
    for mode in BRANCH_MODES:
        name = branch_name(mode)
        support = unit_support(mode)
        histogram = distance_histogram(support)
        rows.append(
            {
                "name": name,
                "mode": mode,
                "support": list(support),
                "branch_pair": branch_pair(name),
                "required_chi11_value": required_chi11_value(name),
                "distance_histogram_d1_to_d6": list(histogram),
                "full_aut_orbit_histogram_O15_O2_O3_O4_O6": list(full_aut_orbit_histogram(histogram)),
                "gap_necklace": list(canonical_necklace(cyclic_gap_tuple(support))),
            }
        )
    return rows


def support_name_map() -> dict[tuple[int, ...], str]:
    return {unit_support(mode): branch_name(mode) for mode in BRANCH_MODES}


def unit_action_rows() -> list[dict[str, Any]]:
    names = support_name_map()
    rows = []
    for unit in UNITS:
        image = {name: names[multiply_support(unit, support)] for support, name in names.items()}
        rows.append(
            {
                "unit": unit,
                "support_action": image,
                "required_character_value": 1 if unit in {1, 11} else -1,
                "preserves_d5_pair": branch_pair(image["A5_k5"]) == D5_PAIR,
            }
        )
    return rows


def candidate_specs() -> list[dict[str, Any]]:
    return [
        {
            "name": "support_size",
            "kind": "full_aut_invariant_scalar",
            "requires_shell_label": False,
            "value_fn": lambda _name, support, _histogram: len(support),
        },
        {
            "name": "full_aut_orbit_histogram_O15_O2_O3_O4_O6",
            "kind": "full_aut_invariant_scalar",
            "requires_shell_label": False,
            "value_fn": lambda _name, _support, histogram: full_aut_orbit_histogram(histogram),
        },
        {
            "name": "antipodal_pair_count_d6",
            "kind": "full_aut_invariant_scalar",
            "requires_shell_label": False,
            "value_fn": lambda _name, _support, histogram: histogram[5],
        },
        {
            "name": "dihedral_gap_necklace",
            "kind": "dihedral_invariant_not_full_aut",
            "requires_shell_label": "reduced cyclic-order premise",
            "value_fn": lambda _name, support, _histogram: canonical_necklace(cyclic_gap_tuple(support)),
        },
        {
            "name": "raw_distance_histogram_d1_to_d6",
            "kind": "shell_labelled_record",
            "requires_shell_label": True,
            "value_fn": lambda _name, _support, histogram: histogram,
        },
        {
            "name": "shell_labelled_d5_minus_d1_count",
            "kind": "chi11_covariant_but_imports_shell_label",
            "requires_shell_label": True,
            "value_fn": lambda _name, _support, histogram: shell_labelled_d5_minus_d1(histogram),
        },
        {
            "name": "shell_labelled_energy_difference_(d5+d6)-(d1+d6)",
            "kind": "chi11_covariant_but_imports_shell_label",
            "requires_shell_label": True,
            "value_fn": lambda _name, _support, histogram: shell_labelled_d5_plus_d6_energy(histogram) - shell_labelled_d1_plus_d6_energy(histogram),
        },
    ]


def normalize_value(value: Any) -> Any:
    if isinstance(value, tuple):
        return list(value)
    return value


def candidate_row(spec: dict[str, Any]) -> dict[str, Any]:
    value_fn: Callable[[str, tuple[int, ...], tuple[int, int, int, int, int, int]], Any] = spec["value_fn"]
    values = {}
    for mode in BRANCH_MODES:
        name = branch_name(mode)
        support = unit_support(mode)
        histogram = distance_histogram(support)
        values[name] = value_fn(name, support, histogram)

    values_by_branch = {name: normalize_value(value) for name, value in values.items()}
    value_packets_by_pair = {
        CONTIGUOUS_PAIR: sorted({json.dumps(normalize_value(values[name]), sort_keys=True) for name in values if branch_pair(name) == CONTIGUOUS_PAIR}),
        D5_PAIR: sorted({json.dumps(normalize_value(values[name]), sort_keys=True) for name in values if branch_pair(name) == D5_PAIR}),
    }
    constant_on_pairs = all(len(packets) == 1 for packets in value_packets_by_pair.values())
    distinguishes_pairs = constant_on_pairs and value_packets_by_pair[CONTIGUOUS_PAIR] != value_packets_by_pair[D5_PAIR]

    numeric_values = [values[name] for name in values]
    numeric_chi11_covariant = all(isinstance(value, int) for value in numeric_values) and all(
        values[name] == required_chi11_value(name) * abs(values["A5_k5"])
        for name in values
    ) and abs(values["A5_k5"]) > 0

    full_aut_invariant_on_branch_orbit = len({json.dumps(normalize_value(value), sort_keys=True) for value in values.values()}) == 1
    exports_allowed_strict_chi11_source = bool(
        numeric_chi11_covariant and not spec["requires_shell_label"] and spec["kind"] != "dihedral_invariant_not_full_aut"
    )
    return {
        "candidate": spec["name"],
        "kind": spec["kind"],
        "requires_shell_label_or_reduced_symmetry": spec["requires_shell_label"],
        "values_by_branch": values_by_branch,
        "constant_on_reversal_pairs": constant_on_pairs,
        "distinguishes_contiguous_pair_from_d5_pair": distinguishes_pairs,
        "full_aut_invariant_on_four_branch_orbit": full_aut_invariant_on_branch_orbit,
        "numeric_chi11_covariant_on_branches": numeric_chi11_covariant,
        "exports_allowed_strict_chi11_source": exports_allowed_strict_chi11_source,
        "interpretation": (
            "candidate would be an allowed strict chi_11 source"
            if exports_allowed_strict_chi11_source
            else "does not supply chi_11 without importing shell-label/reduced symmetry"
        ),
    }


def candidate_rows() -> list[dict[str, Any]]:
    return [candidate_row(spec) for spec in candidate_specs()]


def exact_proof_certificate(rows: list[dict[str, Any]]) -> dict[str, str]:
    allowed_sources = [row["candidate"] for row in rows if row["exports_allowed_strict_chi11_source"]]
    covariant_imports = [row["candidate"] for row in rows if row["numeric_chi11_covariant_on_branches"] and row["requires_shell_label_or_reduced_symmetry"]]
    invariant_distinguishers = [
        row["candidate"]
        for row in rows
        if row["full_aut_invariant_on_four_branch_orbit"] and row["distinguishes_contiguous_pair_from_d5_pair"]
    ]
    return {
        "finite_domain": "The audit is over the four unit branches A1,A5,A7,A11 and the exact unit action of Aut(Z_12)={1,5,7,11}.",
        "required_transform": "A chi_11 source must be constant on reversal pairs, separate contiguous_pair from d5_pair, and flip under units 5 and 7 with kernel {1,11}.",
        "invariant_no_go": f"Full-Aut invariant candidates that distinguish the branch pairs: {invariant_distinguishers}; hence invariant scalar records do not provide the bit.",
        "covariant_imports": f"Candidates with chi_11 covariance are {covariant_imports}, but each requires a d1-vs-d5 shell-label or reduced-symmetry premise.",
        "allowed_source_count": f"Allowed strict chi_11 source candidates found without imported shell labels: {len(allowed_sources)}.",
        "boundary": "The computation classifies candidate records; it does not prove that the listed candidate family exhausts all possible strict nadsoliton sources.",
    }


def build_payload() -> dict[str, Any]:
    rows = candidate_rows()
    return {
        "result_kind": "SCRATCH_STRICT_ALPHA_CHI11_CANDIDATE_SOURCE_TRANSFORM_AUDIT_PROBE__NO_GO_NOT_A_THEOREM",
        "status": "candidate-records-do-not-export-chi11-without-shell-label-or-reduced-symmetry-import",
        "finite_model": {
            "ring": "Z_12",
            "active_count": ACTIVE_COUNT,
            "automorphism_units": UNITS,
            "branch_modes": BRANCH_MODES,
            "required_chi11_kernel": [1, 11],
            "required_chi11_nonzero_coset": [5, 7],
            "target_q_power": TARGET_Q_POWER,
            "target_eta": TARGET_ETA,
        },
        "branch_rows": branch_rows(),
        "unit_action_rows": unit_action_rows(),
        "candidate_transform_rows": rows,
        "candidate_summary": {
            "candidate_count": len(rows),
            "full_aut_invariant_candidate_count": sum(row["full_aut_invariant_on_four_branch_orbit"] for row in rows),
            "pair_distinguishing_candidate_count": sum(row["distinguishes_contiguous_pair_from_d5_pair"] for row in rows),
            "numeric_chi11_covariant_candidate_count": sum(row["numeric_chi11_covariant_on_branches"] for row in rows),
            "allowed_strict_chi11_source_candidate_count": sum(row["exports_allowed_strict_chi11_source"] for row in rows),
            "chi11_covariant_but_imported_candidates": [
                row["candidate"]
                for row in rows
                if row["numeric_chi11_covariant_on_branches"] and row["requires_shell_label_or_reduced_symmetry"]
            ],
        },
        "exact_proof_certificate": exact_proof_certificate(rows),
        "interpretation": {
            "honest_positive": "The audit pinpoints records that have the right chi_11 covariance once shell labels are supplied.",
            "honest_negative": "No audited full-Aut invariant branch record exports chi_11 without importing a shell-label or reduced-symmetry premise.",
            "relation_to_previous_probe": "This attacks the source question left open by the conditional uniqueness and local-necessity certificates for d1+d6.",
        },
        "ontology_guardrail": {
            "allowed_reading": "The nadsoliton itself remains the only admissible place where a future strict chi_11 source could be derived.",
            "forbidden_reading": "No separate informational layer underneath the nadsoliton is introduced to supply chi_11.",
        },
        "hard_limits": [
            "No identity K_legacy_ont == K_strict_gate is used or claimed.",
            "No legacy role transfer to K_strict_gate, alpha_geo, beta_tors, or D_f is made.",
            "No theorem derives the chi_11-kernel, shell-label d1 vs d5, unit-axis bit, exact-cover clauses, or cardinality 5 from strict nadsoliton geometry.",
            "The result is a finite candidate-source transform audit over listed branch records, not an exhaustive strict-source theorem.",
            "No QW-2191 discharge.",
            "No ToE closure.",
        ],
    }


def write_markdown(payload: dict[str, Any]) -> str:
    summary = payload["candidate_summary"]
    proof = payload["exact_proof_certificate"]
    lines = [
        "# Strict alpha chi_11 candidate-source transform audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Candidate summary",
        "",
        f"- Candidate count: `{summary['candidate_count']}`",
        f"- Full-Aut invariant candidate count: `{summary['full_aut_invariant_candidate_count']}`",
        f"- Pair-distinguishing candidate count: `{summary['pair_distinguishing_candidate_count']}`",
        f"- Numeric chi_11-covariant candidate count: `{summary['numeric_chi11_covariant_candidate_count']}`",
        f"- Allowed strict chi_11 source candidate count: `{summary['allowed_strict_chi11_source_candidate_count']}`",
        f"- chi_11-covariant but imported candidates: `{summary['chi11_covariant_but_imported_candidates']}`",
        "",
        "## Proof certificate",
        "",
    ]
    for key, value in proof.items():
        lines.append(f"- `{key}`: {value}")
    lines.extend(["", "## Candidate rows", ""])
    for row in payload["candidate_transform_rows"]:
        lines.append(
            f"- `{row['candidate']}`: kind=`{row['kind']}`, "
            f"distinguishes_pairs=`{row['distinguishes_contiguous_pair_from_d5_pair']}`, "
            f"chi11_covariant=`{row['numeric_chi11_covariant_on_branches']}`, "
            f"allowed_source=`{row['exports_allowed_strict_chi11_source']}`"
        )
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
