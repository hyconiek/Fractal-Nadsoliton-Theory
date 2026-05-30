#!/usr/bin/env python3
"""Scratch probe: quadratic histogram invariants cannot source chi_11.

The bounded linear audit proved that every chi_11-covariant linear histogram
functional imports the d1-vs-d5 shell label.  This probe takes the next finite
step: degree <= 2 histogram features.  Instead of enumerating a huge coefficient
cube, it enumerates the full-Aut monomial orbits under the only nontrivial shell
identification relevant here, d1 <-> d5, and evaluates their orbit sums on the
four unit branches A1,A5,A7,A11.

Result: all 21 full-Aut invariant monomial orbit sums of degree <= 2 have the
same value on the contiguous pair and the d5 pair.  Therefore any quadratic
polynomial built from these invariant orbit sums is unable to distinguish A1
from A5 and cannot export chi_11.  The anti-invariant orbit differences do flip
with chi_11, but exactly because they choose an orientation of the d1/d5 shell
pair, i.e. they import the missing shell-label premise.

No false pass: this is a finite degree<=2 histogram invariant no-go, not an
exhaustive strict-source theorem over all non-histogram or higher-order data.  It
does not derive chi_11, discharge QW-2191, or close ToE.
"""
from __future__ import annotations

import json
from itertools import combinations, combinations_with_replacement
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
OUT_JSON = HERE / "bridge_strict_alpha_quadratic_histogram_invariant_chi11_no_go_report.json"
OUT_MD = HERE / "bridge_strict_alpha_quadratic_histogram_invariant_chi11_no_go_report.md"

N = 12
ACTIVE_COUNT = 5
UNITS = [1, 5, 7, 11]
BRANCH_MODES = [1, 5, 7, 11]
TARGET_Q_POWER = "256/243"
TARGET_ETA = "9/5"
CONTIGUOUS_PAIR = "contiguous_pair_A1_A11"
D5_PAIR = "d5_pair_A5_A7"
SHELL_LABELS = ["d1", "d2", "d3", "d4", "d5", "d6"]
A1_HISTOGRAM = (4, 3, 2, 1, 0, 0)
A5_HISTOGRAM = (0, 3, 2, 1, 4, 0)


def folded(value: int) -> int:
    residue = value % N
    return min(residue, N - residue)


def unit_support(mode: int) -> tuple[int, ...]:
    return tuple(sorted((mode * index) % N for index in range(ACTIVE_COUNT)))


def distance_histogram(support: tuple[int, ...]) -> tuple[int, int, int, int, int, int]:
    counts = [0] * (N // 2)
    for left, right in combinations(support, 2):
        counts[folded(right - left) - 1] += 1
    return tuple(counts)  # type: ignore[return-value]


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


def swap_d1_d5_index(index: int) -> int:
    if index == 0:
        return 4
    if index == 4:
        return 0
    return index


def swap_monomial(monomial: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(sorted(swap_d1_d5_index(index) for index in monomial))


def monomial_label(monomial: tuple[int, ...]) -> str:
    return "*".join(SHELL_LABELS[index] for index in monomial)


def monomial_value(monomial: tuple[int, ...], histogram: tuple[int, ...]) -> int:
    result = 1
    for index in monomial:
        result *= histogram[index]
    return result


def monomial_list_degree_le_2() -> list[tuple[int, ...]]:
    linear = [(index,) for index in range(len(SHELL_LABELS))]
    quadratic = list(combinations_with_replacement(range(len(SHELL_LABELS)), 2))
    return linear + quadratic


def monomial_orbits_degree_le_2() -> list[list[tuple[int, ...]]]:
    seen: set[tuple[int, ...]] = set()
    orbits = []
    for monomial in monomial_list_degree_le_2():
        if monomial in seen:
            continue
        orbit = sorted({monomial, swap_monomial(monomial)})
        seen.update(orbit)
        orbits.append(orbit)
    return orbits


def branch_rows() -> list[dict[str, Any]]:
    return [
        {
            "name": branch_name(mode),
            "mode": mode,
            "support": list(unit_support(mode)),
            "branch_pair": branch_pair(branch_name(mode)),
            "required_chi11_value": required_chi11_value(branch_name(mode)),
            "distance_histogram_d1_to_d6": list(distance_histogram(unit_support(mode))),
        }
        for mode in BRANCH_MODES
    ]


def invariant_orbit_rows() -> list[dict[str, Any]]:
    rows = []
    for orbit in monomial_orbits_degree_le_2():
        value_a1 = sum(monomial_value(monomial, A1_HISTOGRAM) for monomial in orbit)
        value_a5 = sum(monomial_value(monomial, A5_HISTOGRAM) for monomial in orbit)
        rows.append(
            {
                "orbit_monomials": [monomial_label(monomial) for monomial in orbit],
                "degree": len(orbit[0]),
                "orbit_size": len(orbit),
                "invariant_orbit_sum_value_A1_A11": value_a1,
                "invariant_orbit_sum_value_A5_A7": value_a5,
                "distinguishes_branch_pairs": value_a1 != value_a5,
            }
        )
    return rows


def anti_invariant_rows() -> list[dict[str, Any]]:
    rows = []
    for orbit in monomial_orbits_degree_le_2():
        if len(orbit) != 2:
            continue
        left, right = orbit
        value_a1 = monomial_value(right, A1_HISTOGRAM) - monomial_value(left, A1_HISTOGRAM)
        value_a5 = monomial_value(right, A5_HISTOGRAM) - monomial_value(left, A5_HISTOGRAM)
        rows.append(
            {
                "anti_orbit_difference": f"{monomial_label(right)} - {monomial_label(left)}",
                "degree": len(left),
                "value_A1_A11": value_a1,
                "value_A5_A7": value_a5,
                "numeric_chi11_covariant": value_a5 == -value_a1 and value_a5 != 0,
                "imports_d1_vs_d5_shell_orientation": True,
            }
        )
    return rows


def symbolic_certificate() -> dict[str, str]:
    return {
        "histogram_relation": "A1/A11 has h=(4,3,2,1,0,0) and A5/A7 has h'=(0,3,2,1,4,0), exactly the d1<->d5 swap of h.",
        "invariant_principle": "Any full-Aut invariant histogram polynomial P satisfies P(h)=P(swap_d1_d5(h)), so P(A1)=P(A5).",
        "degree2_enumeration": "The degree<=2 monomial basis has 6 linear + 21 quadratic monomials; quotienting by d1<->d5 gives 21 invariant orbit sums.",
        "anti_invariant_boundary": "Orbit differences can flip sign under d1<->d5, but choosing an orbit difference chooses an orientation of the d1/d5 pair and therefore imports the missing shell-label bit.",
    }


def build_payload() -> dict[str, Any]:
    invariant_rows = invariant_orbit_rows()
    anti_rows = anti_invariant_rows()
    return {
        "result_kind": "SCRATCH_STRICT_ALPHA_QUADRATIC_HISTOGRAM_INVARIANT_CHI11_NO_GO_PROBE__NO_GO_NOT_A_THEOREM",
        "status": "degree2-full-aut-histogram-invariants-cannot-distinguish-a1-from-a5-or-export-chi11",
        "finite_model": {
            "ring": "Z_12",
            "active_count": ACTIVE_COUNT,
            "automorphism_units": UNITS,
            "branch_modes": BRANCH_MODES,
            "shell_labels": SHELL_LABELS,
            "monomial_count_degree_le_2": len(monomial_list_degree_le_2()),
            "invariant_orbit_sum_count_degree_le_2": len(invariant_rows),
            "target_q_power": TARGET_Q_POWER,
            "target_eta": TARGET_ETA,
        },
        "branch_rows": branch_rows(),
        "symbolic_certificate": symbolic_certificate(),
        "invariant_orbit_rows": invariant_rows,
        "anti_invariant_rows": anti_rows,
        "quadratic_invariant_summary": {
            "invariant_orbit_sum_count": len(invariant_rows),
            "invariant_orbit_sum_pair_distinguishing_count": sum(row["distinguishes_branch_pairs"] for row in invariant_rows),
            "anti_invariant_orbit_difference_count": len(anti_rows),
            "anti_invariant_chi11_covariant_count": sum(row["numeric_chi11_covariant"] for row in anti_rows),
            "allowed_strict_chi11_source_count": 0,
        },
        "exact_proof_certificate": {
            "finite_domain": "All degree<=2 histogram monomials and their d1<->d5 full-Aut orbits are enumerated exactly.",
            "invariant_no_go": "Every one of the 21 invariant orbit sums has equal value on A1/A11 and A5/A7, so their linear span cannot export chi_11.",
            "symbolic_reason": symbolic_certificate()["invariant_principle"],
            "anti_invariant_boundary": "Five anti-invariant orbit differences are nonzero chi_11-covariant, but each imports an orientation of the d1/d5 shell pair.",
            "scope_limit": "This is exhaustive only for degree<=2 histogram invariants, not for all possible strict nadsoliton source mechanisms.",
        },
        "interpretation": {
            "honest_positive": "The computation identifies low-degree anti-invariant expressions that behave like chi_11 once d1/d5 orientation is supplied.",
            "honest_negative": "No full-Aut invariant degree<=2 histogram expression can distinguish the branch pairs or supply chi_11.",
            "relation_to_previous_probe": "This extends the bounded linear no-go to the complete quadratic invariant orbit-sum basis.",
        },
        "ontology_guardrail": {
            "allowed_reading": "The nadsoliton itself remains the only admissible place where a future non-histogram or higher-order strict chi_11 source could be derived.",
            "forbidden_reading": "No separate informational layer underneath the nadsoliton is introduced to supply chi_11.",
        },
        "hard_limits": [
            "No identity K_legacy_ont == K_strict_gate is used or claimed.",
            "No legacy role transfer to K_strict_gate, alpha_geo, beta_tors, or D_f is made.",
            "No theorem derives the chi_11-kernel, shell-label d1 vs d5, unit-axis bit, exact-cover clauses, or cardinality 5 from strict nadsoliton geometry.",
            "The result is a finite degree<=2 histogram invariant no-go, not an exhaustive strict-source theorem.",
            "No QW-2191 discharge.",
            "No ToE closure.",
        ],
    }


def write_markdown(payload: dict[str, Any]) -> str:
    summary = payload["quadratic_invariant_summary"]
    lines = [
        "# Strict alpha quadratic histogram invariant chi_11 no-go",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Summary",
        "",
        f"- Invariant orbit sums: `{summary['invariant_orbit_sum_count']}`",
        f"- Invariant pair-distinguishing orbit sums: `{summary['invariant_orbit_sum_pair_distinguishing_count']}`",
        f"- Anti-invariant orbit differences: `{summary['anti_invariant_orbit_difference_count']}`",
        f"- Anti-invariant chi_11-covariant differences: `{summary['anti_invariant_chi11_covariant_count']}`",
        f"- Allowed strict chi_11 source count: `{summary['allowed_strict_chi11_source_count']}`",
        "",
        "## Symbolic certificate",
        "",
    ]
    for key, value in payload["symbolic_certificate"].items():
        lines.append(f"- `{key}`: {value}")
    lines.extend(["", "## Proof certificate", ""])
    for key, value in payload["exact_proof_certificate"].items():
        lines.append(f"- `{key}`: {value}")
    lines.extend(["", "## Anti-invariant chi_11-covariant rows", ""])
    for row in payload["anti_invariant_rows"]:
        if row["numeric_chi11_covariant"]:
            lines.append(
                f"- `{row['anti_orbit_difference']}`: A1=`{row['value_A1_A11']}`, "
                f"A5=`{row['value_A5_A7']}`, imports_label=`{row['imports_d1_vs_d5_shell_orientation']}`"
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
