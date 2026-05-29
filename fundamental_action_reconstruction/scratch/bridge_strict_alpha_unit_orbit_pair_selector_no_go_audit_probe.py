#!/usr/bin/env python3
"""Scratch probe: full unit symmetry cannot select the d5 branch pair.

The previous full-Aut shell-energy audit showed a shell-linear no-go: if a
scalar shell energy is invariant under all units of Z_12, then it cannot
distinguish A5/d5 from A1/contiguous.  This probe removes the shell-linearity
assumption and audits the branch set directly.

Use the anchored unit-generated 5-supports

    A1  = {0,1,2,3,4}
    A5  = {0,5,10,3,8}
    A7  = {0,7,2,9,4}
    A11 = {0,11,10,9,8}.

Orientation reversal pairs them into two unoriented branch pairs:

    contiguous_pair = {A1, A11}
    d5_pair         = {A5, A7}.

Full Aut(Z_12)={1,5,7,11} swaps these two pairs; the subgroup {1,11} preserves
the d5 pair.  Therefore a selector on branch pairs that is full-Aut invariant
can select neither singleton pair.  It can only select none or both.  Choosing
just d5_pair is possible only after reducing symmetry to the chi_11 kernel
{1,11}, or after adding an equivalent selector premise.

No false pass: this is a finite no-go/conditionality audit, not a derivation of
chi_11 from strict core, not a QW-2191 discharge, and not ToE closure.
"""
from __future__ import annotations

import json
from itertools import combinations
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
OUT_JSON = HERE / "bridge_strict_alpha_unit_orbit_pair_selector_no_go_audit_report.json"
OUT_MD = HERE / "bridge_strict_alpha_unit_orbit_pair_selector_no_go_audit_report.md"

N = 12
ACTIVE_COUNT = 5
UNITS = [1, 5, 7, 11]
TARGET_Q_POWER = "256/243"
TARGET_ETA = "9/5"
BRANCH_MODES = [1, 5, 7, 11]
CONTIGUOUS_PAIR = "contiguous_pair_A1_A11"
D5_PAIR = "d5_pair_A5_A7"


def unit_support(mode: int) -> tuple[int, ...]:
    return tuple(sorted((mode * index) % N for index in range(ACTIVE_COUNT)))


def folded(value: int) -> int:
    residue = value % N
    return min(residue, N - residue)


def distance_histogram(support: tuple[int, ...]) -> tuple[int, int, int, int, int, int]:
    counts = [0] * (N // 2)
    for left, right in combinations(support, 2):
        counts[folded(right - left) - 1] += 1
    return tuple(counts)  # type: ignore[return-value]


def multiply_support(unit: int, support: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(sorted((unit * value) % N for value in support))


def support_name_by_tuple() -> dict[tuple[int, ...], str]:
    return {unit_support(mode): f"A{mode}_k{mode}" for mode in BRANCH_MODES}


def pair_for_support_name(name: str) -> str:
    if name in {"A1_k1", "A11_k11"}:
        return CONTIGUOUS_PAIR
    if name in {"A5_k5", "A7_k7"}:
        return D5_PAIR
    raise ValueError(f"unknown support name {name}")


def branch_support_rows() -> list[dict[str, Any]]:
    return [
        {
            "name": f"A{mode}_k{mode}",
            "mode": mode,
            "support": list(unit_support(mode)),
            "distance_histogram_d1_to_d6": list(distance_histogram(unit_support(mode))),
            "orientation_pair": pair_for_support_name(f"A{mode}_k{mode}"),
        }
        for mode in BRANCH_MODES
    ]


def unit_action_on_supports() -> list[dict[str, Any]]:
    names = support_name_by_tuple()
    rows = []
    for unit in UNITS:
        image = {}
        for support, name in names.items():
            image[name] = names[multiply_support(unit, support)]
        rows.append(
            {
                "unit": unit,
                "support_action": image,
                "pair_action": {
                    CONTIGUOUS_PAIR: pair_for_support_name(image["A1_k1"]),
                    D5_PAIR: pair_for_support_name(image["A5_k5"]),
                },
                "preserves_d5_pair": pair_for_support_name(image["A5_k5"]) == D5_PAIR,
            }
        )
    return rows


def subset_name(subset: set[str]) -> str:
    if not subset:
        return "empty_selector"
    if subset == {CONTIGUOUS_PAIR, D5_PAIR}:
        return "both_pairs_unoriented_orbit"
    if subset == {CONTIGUOUS_PAIR}:
        return "singleton_contiguous_pair"
    if subset == {D5_PAIR}:
        return "singleton_d5_pair"
    return "unknown_subset"


def transform_pair_subset(subset: set[str], pair_action: dict[str, str]) -> set[str]:
    return {pair_action[pair] for pair in subset}


def selector_subset_audit(action_rows: list[dict[str, Any]]) -> dict[str, Any]:
    pair_universe = [CONTIGUOUS_PAIR, D5_PAIR]
    subsets = [set(), {CONTIGUOUS_PAIR}, {D5_PAIR}, {CONTIGUOUS_PAIR, D5_PAIR}]
    full_rows = []
    reduced_rows = []
    reduced_units = [1, 11]
    for subset in subsets:
        full_invariant = all(
            transform_pair_subset(subset, row["pair_action"]) == subset
            for row in action_rows
        )
        reduced_invariant = all(
            transform_pair_subset(subset, row["pair_action"]) == subset
            for row in action_rows
            if row["unit"] in reduced_units
        )
        row = {
            "selector": subset_name(subset),
            "selected_pairs": sorted(subset),
            "selects_singleton_d5_pair": subset == {D5_PAIR},
        }
        full_rows.append({**row, "invariant_under": "full_Aut_Z12", "is_invariant": full_invariant})
        reduced_rows.append({**row, "invariant_under": "chi_11_kernel_units_{1,11}", "is_invariant": reduced_invariant})
    return {
        "pair_universe": pair_universe,
        "full_Aut_invariant_selectors": [row for row in full_rows if row["is_invariant"]],
        "chi_11_kernel_invariant_selectors": [row for row in reduced_rows if row["is_invariant"]],
        "all_selector_rows": full_rows + reduced_rows,
        "full_Aut_singleton_d5_selector_count": sum(
            1 for row in full_rows if row["is_invariant"] and row["selects_singleton_d5_pair"]
        ),
        "chi_11_kernel_singleton_d5_selector_count": sum(
            1 for row in reduced_rows if row["is_invariant"] and row["selects_singleton_d5_pair"]
        ),
    }


def exact_proof_certificate() -> dict[str, str]:
    return {
        "unit_orbit": "The four anchored unit supports A1,A5,A7,A11 form one full Aut(Z_12) orbit.",
        "orientation_pairs": "Orientation reversal groups them into contiguous_pair={A1,A11} and d5_pair={A5,A7}.",
        "pair_action": "Units 5 and 7 swap contiguous_pair with d5_pair; units 1 and 11 preserve each pair.",
        "selector_no_go": "A full-Aut invariant subset of the two-pair universe must be empty or both pairs; singleton d5_pair is not invariant.",
        "conditional_positive": "After reducing symmetry to {1,11}, singleton d5_pair is invariant, exactly matching the chi_11-kernel condition.",
        "source_obstruction": "The reduction from full Aut to {1,11} is not derived here from strict core; it is the missing selector premise.",
    }


def main() -> None:
    support_rows = branch_support_rows()
    action_rows = unit_action_on_supports()
    selector_audit = selector_subset_audit(action_rows)

    report: dict[str, Any] = {
        "result_kind": "SCRATCH_STRICT_ALPHA_UNIT_ORBIT_PAIR_SELECTOR_NO_GO_AUDIT_PROBE__NO_GO_NOT_A_THEOREM",
        "status": "full-aut-on-branch-pairs-forbids-singleton-d5-selector-chi_11-kernel-allows-it-conditionally",
        "finite_model": {
            "ring": "Z_12",
            "active_count": ACTIVE_COUNT,
            "automorphism_units": UNITS,
            "branch_modes": BRANCH_MODES,
            "target_q_power": TARGET_Q_POWER,
            "target_eta": TARGET_ETA,
        },
        "branch_supports": support_rows,
        "unit_action_on_supports_and_pairs": action_rows,
        "selector_subset_audit": selector_audit,
        "exact_proof_certificate": exact_proof_certificate(),
        "interpretation": {
            "direct_no_go": "On the two unoriented branch pairs, full Aut symmetry swaps contiguous_pair and d5_pair, so it cannot select singleton d5_pair.",
            "conditional_positive": "Reducing symmetry to the chi_11 kernel {1,11} makes singleton d5_pair invariant and therefore selectable.",
            "relation_to_previous_probe": "This removes the shell-linear restriction: the obstruction is already present at the finite branch-pair action level.",
            "honest_limit": "The symmetry reduction is not proven from strict nadsoliton geometry here; QW-2191 remains open.",
        },
        "ontology_guardrail": {
            "allowed_reading": "The nadsoliton itself may later export the symmetry reduction to the chi_11 kernel.",
            "forbidden_reading": "No separate informational layer underneath the nadsoliton is introduced to choose d5_pair.",
            "preferred_order_preserved": "nadsoliton -> light -> matter -> emergent observer",
        },
        "hard_limits": [
            "No identity K_legacy_ont == K_strict_gate is asserted or used.",
            "No legacy physical-role or matter-generation claim is transferred onto K_strict_gate.",
            "No theorem derives the chi_11-kernel symmetry reduction or the unit-axis bit from strict nadsoliton geometry.",
            "This is a finite branch-pair no-go audit, not a global impossibility theorem for all conceivable strict sources.",
            "No endpoint, arrow orientation, ledger selector, positive lambda action, cycle metric source, anti-Nyquist source, fifth-mode source, future-probability source, future-value source, future-path source, matter-bit source, existence-bit source, recursive-self-information source, character-source, meta-character source, cyclic-adjacency source theorem, variational tie-break source theorem, exclusion-energy source theorem, exact-cover source theorem, dual-shell-label source theorem, full-Aut shell-linear selector theorem, or legacy-strict bridge theorem is claimed.",
            "No QW-2191 discharge and no strict-core selector closure are claimed.",
            "No ToE closure is claimed.",
        ],
        "next_honest_step": "Search for a strict source of the symmetry reduction full Aut -> {1,11}; without it, singleton d5_pair selection remains conditional/non-strict.",
    }

    OUT_JSON.write_text(json.dumps(report, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    OUT_MD.write_text(
        "# Scratch strict-alpha unit-orbit pair selector no-go audit probe\n\n"
        "Status: full Aut forbids singleton d5-pair selection; chi_11-kernel reduction allows it conditionally.\n\n"
        f"- Branch modes audited: `{BRANCH_MODES}`.\n"
        f"- Full-Aut invariant singleton d5 selectors: `{selector_audit['full_Aut_singleton_d5_selector_count']}`.\n"
        f"- chi_11-kernel invariant singleton d5 selectors: `{selector_audit['chi_11_kernel_singleton_d5_selector_count']}`.\n"
        f"- Full-Aut invariant selectors: `{selector_audit['full_Aut_invariant_selectors']}`.\n"
        f"- chi_11-kernel invariant selectors: `{selector_audit['chi_11_kernel_invariant_selectors']}`.\n"
        f"- Target replay kept conditional: `q^5={TARGET_Q_POWER}`, eta `{TARGET_ETA}`.\n"
        "- No false pass: no strict symmetry-reduction source theorem, no QW-2191 discharge, no ToE closure.\n",
        encoding="utf-8",
    )
    print(OUT_JSON)
    print(OUT_MD)


if __name__ == "__main__":
    main()
