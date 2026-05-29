#!/usr/bin/env python3
"""Scratch probe: cyclic adjacency is the minimal unit-label candidate for chi_11.

The previous meta-symmetry audit showed that the abstract unit group
Aut(Z_12) ~= V4 does not uniquely choose the required character chi_11.  This
probe checks the next honest candidate: not another scalar record, but a
specific Z_12 cyclic-adjacency datum.  If the strict-side geometry is allowed to
name the nearest-neighbour edge class {+1,-1}, then the unit subgroup preserving
that edge class is exactly {1,11}; this is the kernel of the required chi_11.

This is a conditional source audit, not a strict-core discharge.  It proves an
exact finite statement:

    Stab_Units(shell_1) = {u in {1,5,7,11}: folded_distance(u*1)=1} = {1,11}.

Therefore a cyclic-adjacency / nearest-neighbour premise would be sufficient to
select chi_11 from the three-character meta-orbit.  But the premise itself is
not derived here from strict nadsoliton geometry, so QW-2191 remains open unless
that premise is exported by an accepted strict source theorem.
"""
from __future__ import annotations

import json
from collections import Counter
from itertools import combinations, product
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
OUT_JSON = HERE / "bridge_strict_alpha_cyclic_adjacency_unit_label_source_audit_report.json"
OUT_MD = HERE / "bridge_strict_alpha_cyclic_adjacency_unit_label_source_audit_report.md"

N = 12
ACTIVE_COUNT = 5
UNITS = [1, 5, 7, 11]
TARGET_Q_POWER = "256/243"
TARGET_ETA = "9/5"
A1 = "A1_k1_contiguous"
A5 = "A5_k5_d5"
D5_HISTOGRAM = [0, 3, 2, 1, 4, 0]
CONTIGUOUS_HISTOGRAM = [4, 3, 2, 1, 0, 0]


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


def multiply_units(left: int, right: int) -> int:
    return (left * right) % N


def shell(unit: int) -> int:
    return folded(unit)


def shell_stabilizer(target_shell: int) -> list[int]:
    return [unit for unit in UNITS if shell(unit) == target_shell]


def enumerate_shell_rows() -> list[dict[str, Any]]:
    rows = []
    for target_shell in range(1, N // 2 + 1):
        preserving_units = shell_stabilizer(target_shell)
        rows.append(
            {
                "cyclic_shell": target_shell,
                "preserving_units": preserving_units,
                "preserving_unit_count": len(preserving_units),
                "is_required_kernel_{1,11}": preserving_units == [1, 11],
                "selected_character": character_name_from_kernel(preserving_units),
            }
        )
    return rows


def is_character(values: dict[int, int]) -> bool:
    return values[1] == 0 and all(
        values[multiply_units(left, right)] == (values[left] + values[right]) % 2
        for left in UNITS
        for right in UNITS
    )


def character_name_from_kernel(kernel: list[int]) -> str:
    if kernel == UNITS:
        return "trivial_character"
    if kernel == [1, 11]:
        return "chi_11_required_d5_unit_axis"
    if kernel == [1, 5]:
        return "chi_5_kernel"
    if kernel == [1, 7]:
        return "chi_7_kernel"
    if not kernel:
        return "no_character_empty_kernel"
    return f"no_Z2_character_for_kernel_{kernel}"


def enumerate_characters() -> list[dict[str, Any]]:
    rows = []
    for bits in product([0, 1], repeat=len(UNITS)):
        values = dict(zip(UNITS, bits))
        if not is_character(values):
            continue
        kernel = sorted(unit for unit in UNITS if values[unit] == 0)
        rows.append(
            {
                "name": character_name_from_kernel(kernel),
                "values_mod2": {str(unit): values[unit] for unit in UNITS},
                "kernel": kernel,
                "nonzero_coset": sorted(unit for unit in UNITS if values[unit] == 1),
                "is_required_chi_11": kernel == [1, 11],
            }
        )
    return sorted(rows, key=lambda row: row["name"])


def adjacency_audit() -> dict[str, Any]:
    unit_shells = {str(unit): shell(unit) for unit in UNITS}
    edge_kernel = shell_stabilizer(1)
    long_step_kernel = shell_stabilizer(5)
    return {
        "unit_image_of_generator_plus_one": {str(unit): unit % N for unit in UNITS},
        "unit_image_folded_shell": unit_shells,
        "nearest_neighbor_shell_1_preserving_units": edge_kernel,
        "nearest_neighbor_shell_1_preserving_character": character_name_from_kernel(edge_kernel),
        "long_step_shell_5_preserving_units": long_step_kernel,
        "long_step_shell_5_preserving_character": character_name_from_kernel(long_step_kernel),
        "shell_rows": enumerate_shell_rows(),
        "exact_positive_if_adjacency_admitted": edge_kernel == [1, 11],
        "why_this_breaks_meta_symmetry": "The cyclic graph/ring adjacency labels ±1 as nearest neighbours, so the abstract V4 meta-automorphisms may no longer freely permute 5,7,11.",
    }


def selection_source_audit() -> list[dict[str, Any]]:
    return [
        {
            "candidate_source": "abstract_unit_group_V4_only",
            "selects_kernel": None,
            "selects_required_chi_11": False,
            "strict_status": "insufficient: three nontrivial characters remain meta-symmetric",
        },
        {
            "candidate_source": "cyclic_adjacency_shell_{+1,-1}",
            "selects_kernel": [1, 11],
            "selects_required_chi_11": True,
            "strict_status": "sufficient if this nearest-neighbour geometry is exported as a strict-side source",
        },
        {
            "candidate_source": "choose_long_step_shell_{+5,-5}",
            "selects_kernel": [5, 7],
            "selects_required_chi_11": False,
            "strict_status": "not a character kernel containing identity; does not select required chi_11",
        },
        {
            "candidate_source": "name_11_by_decree_without_adjacency_or_geometry",
            "selects_kernel": [1, 11],
            "selects_required_chi_11": True,
            "strict_status": "works only as explicit non-strict/unit-label axiom unless derived elsewhere",
        },
    ]


def exact_proof_certificate() -> dict[str, str]:
    return {
        "finite_group": "Units(Z_12)={1,5,7,11}; every nonidentity unit has abstract order 2, so V4 alone does not distinguish 11.",
        "cyclic_adjacency_datum": "The cyclic Z_12 nearest-neighbour edge class is the unoriented shell {+1,-1} = {1,11} modulo 12.",
        "stabilizer_computation": "A unit u preserves the nearest-neighbour shell iff folded(u*1)=1; among units this holds exactly for u=1 and u=11.",
        "character_identification": "The subgroup {1,11} is the kernel of chi_11, whose nonzero coset is {5,7}; this is the required d5 unit-axis character.",
        "conditionality": "The proof is conditional on admitting cyclic adjacency as strict-side geometry; it does not derive that source premise from the strict core here.",
    }


def main() -> None:
    supports = all_supports()
    histogram_counter = Counter(distance_histogram(support) for support in supports)
    adj = adjacency_audit()

    report: dict[str, Any] = {
        "result_kind": "SCRATCH_STRICT_ALPHA_CYCLIC_ADJACENCY_UNIT_LABEL_SOURCE_AUDIT_PROBE__CONDITIONAL_NOT_A_THEOREM",
        "status": "cyclic-adjacency-conditionally-selects-chi_11-kernel-but-source-premise-not-derived",
        "finite_model": {
            "ring": "Z_12",
            "active_count": ACTIVE_COUNT,
            "support_count": len(supports),
            "histogram_class_count": len(histogram_counter),
            "automorphism_units": UNITS,
            "survivor_axes": [
                {"name": A1, "mode": 1, "distance_histogram_d1_to_d6": CONTIGUOUS_HISTOGRAM},
                {"name": A5, "mode": 5, "distance_histogram_d1_to_d6": D5_HISTOGRAM},
            ],
            "target_q_power": TARGET_Q_POWER,
            "target_eta": TARGET_ETA,
        },
        "character_table": enumerate_characters(),
        "cyclic_adjacency_source_audit": adj,
        "selection_source_audit": selection_source_audit(),
        "exact_proof_certificate": exact_proof_certificate(),
        "interpretation": {
            "direct_result": "Cyclic nearest-neighbour adjacency is an exact finite source that would name the {1,11} kernel and hence chi_11.",
            "positive_conditional": "If strict nadsoliton geometry exports the unoriented edge shell {+1,-1}, the previous three-character ambiguity collapses to chi_11.",
            "negative_result": "The abstract V4 unit group alone still cannot provide this label; adjacency is extra structure relative to the meta-symmetry audit.",
            "honest_limit": "This audit does not prove that cyclic adjacency is already an accepted strict-core selector source, so QW-2191 is not discharged.",
        },
        "ontology_guardrail": {
            "allowed_reading": "Cyclic adjacency may be investigated as a candidate structure of the nadsoliton itself.",
            "forbidden_reading": "No separate informational layer underneath the nadsoliton is introduced to provide adjacency or the unit-axis bit.",
            "preferred_order_preserved": "nadsoliton -> light -> matter -> emergent observer",
        },
        "hard_limits": [
            "No identity K_legacy_ont == K_strict_gate is asserted or used.",
            "No legacy physical-role or matter-generation claim is transferred onto K_strict_gate.",
            "No theorem derives cyclic adjacency, chi_11, or the unit-axis bit from strict nadsoliton geometry.",
            "Cyclic adjacency selecting {1,11} is a conditional source audit, not a strict-core theorem.",
            "No endpoint, arrow orientation, ledger selector, positive lambda action, cycle metric source, anti-Nyquist source, fifth-mode source, future-probability source, future-value source, future-path source, matter-bit source, existence-bit source, recursive-self-information source, character-source, meta-character source, or legacy-strict bridge theorem is claimed.",
            "No QW-2191 discharge and no strict-core selector closure are claimed.",
            "No ToE closure is claimed.",
        ],
        "next_honest_step": "Audit whether an existing strict-side nad12/cyclic graph object already exports the nearest-neighbour shell {+1,-1}; if yes, test its compatibility with the d5 branch without importing an orientation bit; if no, mark adjacency as an explicit non-strict selector premise.",
    }

    OUT_JSON.write_text(json.dumps(report, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    OUT_MD.write_text(
        "# Scratch strict-alpha cyclic-adjacency unit-label source audit probe\n\n"
        "Status: cyclic adjacency conditionally selects the chi_11 kernel; source premise not derived here.\n\n"
        f"- Supports scanned: `{len(supports)}`; histogram classes: `{len(histogram_counter)}`.\n"
        f"- Units: `{UNITS}`.\n"
        f"- Unit image folded shells: `{adj['unit_image_folded_shell']}`.\n"
        f"- Nearest-neighbour shell preserving units: `{adj['nearest_neighbor_shell_1_preserving_units']}`.\n"
        f"- Selected character from shell 1: `{adj['nearest_neighbor_shell_1_preserving_character']}`.\n"
        f"- Conditional positive result: `{adj['exact_positive_if_adjacency_admitted']}`.\n"
        f"- Target replay kept conditional: `q^5={TARGET_Q_POWER}`, eta `{TARGET_ETA}`.\n"
        "- No false pass: no strict adjacency-source theorem, no QW-2191 discharge, no ToE closure.\n",
        encoding="utf-8",
    )
    print(OUT_JSON)
    print(OUT_MD)


if __name__ == "__main__":
    main()
