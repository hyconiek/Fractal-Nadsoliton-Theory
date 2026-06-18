#!/usr/bin/env python3
"""P2881/S1831: finite variational/unit-law 9/5 derivation no-go audit.

P2880 showed that an endpoint-11 pin plus a denominator-5 coefficient alphabet
only represents 9/5 by coefficient import.  This packet tests the next honest
candidate requested by the P2880 boundary: a small strict-style variational/unit
law on the same radius-1 stencil should derive the 9/5 coefficient as a unique
minimizer, without placing 9/5 directly in a constraint.

The audit enumerates denominator-5 radius-1 stencils, applies source-neutral
unit constraints (sum, moment, boundary balance, and adjacent contrast targets),
and minimizes a finite family of nonnegative local quadratic objectives.  The
finite result is negative: none of the tested variational/unit-law candidates
has a unique minimizer whose center coefficient is 9/5.  The only way to obtain
9/5 remains to put that value in the coefficient alphabet or constraint by hand.
"""
from __future__ import annotations

import json
from fractions import Fraction
from itertools import product
from typing import Any, Callable

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2803_s1753_meringer_scd_import_decode_validation_certificate import sha
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN

P2880 = GEN / "p2880_s1830_origin_pinned_9_over_5_coefficient_import_no_go_audit.json"
OUT = GEN / "p2881_s1831_variational_unit_law_9_over_5_derivation_no_go_audit.json"
MD = GEN / "p2881_s1831_variational_unit_law_9_over_5_derivation_no_go_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

OFFSETS = (-1, 0, 1)
DENOMINATOR = 5
NUMERATORS = tuple(range(-9, 10))
TARGET = Fraction(9, 5)
UNIT_TARGETS = tuple(Fraction(n, 1) for n in (-1, 0, 1))

Stencil = tuple[Fraction, Fraction, Fraction]
Predicate = Callable[[Stencil], bool]
Objective = Callable[[Stencil], Fraction]


def stencils() -> list[Stencil]:
    return [tuple(Fraction(n, DENOMINATOR) for n in nums) for nums in product(NUMERATORS, repeat=3)]


def sum_value(s: Stencil) -> Fraction:
    return sum(s, Fraction(0))


def moment_value(s: Stencil) -> Fraction:
    return sum(Fraction(offset) * coeff for offset, coeff in zip(OFFSETS, s))


def boundary_balance_value(s: Stencil) -> Fraction:
    return s[0] + s[2]


def adjacent_contrast_value(s: Stencil) -> Fraction:
    return s[1] - s[0]


def quadratic_norm(s: Stencil) -> Fraction:
    return sum(x * x for x in s)


def smoothness(s: Stencil) -> Fraction:
    return (s[1] - s[0]) ** 2 + (s[2] - s[1]) ** 2


def boundary_penalty(s: Stencil) -> Fraction:
    return s[0] * s[0] + s[2] * s[2]


def moment_penalty(s: Stencil) -> Fraction:
    return sum(Fraction(offset * offset) * coeff * coeff for offset, coeff in zip(OFFSETS, s))


def curvature_penalty(s: Stencil) -> Fraction:
    return (s[0] - 2 * s[1] + s[2]) ** 2


def constraint_family() -> list[tuple[str, Predicate]]:
    constraints: list[tuple[str, Predicate]] = []
    invariants: list[tuple[str, Callable[[Stencil], Fraction]]] = [
        ("sum", sum_value),
        ("moment", moment_value),
        ("boundary_balance", boundary_balance_value),
        ("adjacent_contrast", adjacent_contrast_value),
    ]
    for name, fn in invariants:
        for target in UNIT_TARGETS:
            constraints.append((f"{name}={target}", lambda s, fn=fn, target=target: fn(s) == target))
    constraints.extend(
        [
            ("sum=1_and_moment=0", lambda s: sum_value(s) == 1 and moment_value(s) == 0),
            ("sum=1_and_boundary_balance=0", lambda s: sum_value(s) == 1 and boundary_balance_value(s) == 0),
            ("sum=1_and_adjacent_contrast=0", lambda s: sum_value(s) == 1 and adjacent_contrast_value(s) == 0),
        ]
    )
    return constraints


def objective_family() -> list[tuple[str, Objective]]:
    return [
        ("quadratic_norm", quadratic_norm),
        ("smoothness", smoothness),
        ("boundary_penalty", boundary_penalty),
        ("moment_penalty", moment_penalty),
        ("curvature_penalty", curvature_penalty),
    ]


def minimizer_record(name: str, pred: Predicate, obj_name: str, obj: Objective, universe: list[Stencil]) -> dict[str, Any] | None:
    feasible = [s for s in universe if pred(s)]
    if not feasible:
        return None
    best = min(obj(s) for s in feasible)
    mins = [s for s in feasible if obj(s) == best]
    center_values = sorted({m[1] for m in mins})
    return {
        "constraint": name,
        "objective": obj_name,
        "feasible_count": len(feasible),
        "minimum_value": str(best),
        "minimizer_count": len(mins),
        "center_values_at_minimum": [str(v) for v in center_values],
        "unique_minimizer": len(mins) == 1,
        "unique_minimizer_center": str(mins[0][1]) if len(mins) == 1 else None,
        "derives_center_9_over_5": len(mins) == 1 and mins[0][1] == TARGET,
        "sample_minimizers": [[str(x) for x in s] for s in mins[:6]],
    }


def build_payload(p2880: dict[str, Any]) -> dict[str, Any]:
    universe = stencils()
    records = [
        rec
        for constraint_name, pred in constraint_family()
        for objective_name, obj in objective_family()
        if (rec := minimizer_record(constraint_name, pred, objective_name, obj, universe)) is not None
    ]
    derivation_records = [r for r in records if r["derives_center_9_over_5"]]
    unique_records = [r for r in records if r["unique_minimizer"]]
    center_9_seen_nonunique = [r for r in records if "9/5" in r["center_values_at_minimum"] and not r["unique_minimizer"]]
    facts = {
        "p2880_rechecked": p2880.get("status") == "P2880_ORIGIN_PINNED_9_OVER_5_COEFFICIENT_IMPORT_NO_GO_AUDIT_NO_CLOSURE",
        "denominator_5_grid_checked": len(universe) == 19**3,
        "candidate_count_checked": len(records) == len(constraint_family()) * len(objective_family()),
        "no_candidate_derives_unique_center_9_over_5": len(derivation_records) == 0,
        "some_unique_minimizers_exist_but_not_9_over_5": len(unique_records) > 0 and all(r["unique_minimizer_center"] != "9/5" for r in unique_records),
        "center_9_over_5_only_appears_in_nonunique_minima_when_seen": all(not r["unique_minimizer"] for r in center_9_seen_nonunique),
    }
    return {
        "status": "P2881_VARIATIONAL_UNIT_LAW_9_OVER_5_DERIVATION_NO_GO_AUDIT_NO_CLOSURE",
        "input_hashes": {"P2880": sha(P2880)},
        "variational_unit_law_9_over_5_derivation_no_go_audit": {
            "input_status_rechecked": p2880.get("status"),
            "candidate_class": "source-neutral finite local variational/unit laws on denominator-5 radius-1 stencils, without explicit center=9/5 constraint",
            "stencil_count": len(universe),
            "constraint_count": len(constraint_family()),
            "objective_count": len(objective_family()),
            "candidate_count": len(records),
            "unique_minimizer_candidate_count": len(unique_records),
            "derives_center_9_over_5_candidate_count": len(derivation_records),
            "center_9_over_5_nonunique_minimum_candidate_count": len(center_9_seen_nonunique),
            "sample_candidate_records": records[:12],
            "proof_certificate": {
                "finite_test": "For each source-neutral unit constraint and local quadratic objective, enumerate all denominator-5 radius-1 stencils and compute exact rational minimizers.",
                "finite_result": "No candidate has a unique minimizer with center coefficient 9/5; candidates that have unique minimizers select smaller unit-scale coefficients instead.",
                "sourcehood_step": "A variational/unit law must derive 9/5 as a unique output without putting 9/5 into the constraint.  The audited finite family does not do so.",
            },
        },
        "acceptance_matrix": {
            "facts": facts,
            "accepted_as_independent_variational_unit_law_for_9_over_5": False,
            "exports_unit_bearing_9_over_5_coupling_theorem": False,
            "exports_unit_bearing_action_density": False,
            "exports_translation_origin_source_law": False,
            "exports_strict_damping_compression_bridge": False,
        },
        "decision": {
            "negative_export_flags": {
                "independent_variational_unit_law_exported": False,
                "unit_bearing_9_over_5_coupling_theorem_exported": False,
                "unit_bearing_action_density_exported": False,
                "translation_origin_source_law_exported": False,
                "strict_damping_compression_bridge_exported": False,
                "full_kernel_bridge_exported": False,
                "role_transfer_exported": False,
                "ltotal_exported": False,
                "eom_closure_exported": False,
                "hamiltonian_closure_exported": False,
                "toe_closure_exported": False,
            },
            "reason": "P2881 audits a finite family of source-neutral local variational/unit laws after P2880.  Exact enumeration over the denominator-5 radius-1 grid finds no candidate whose unique minimizer has center coefficient 9/5.  Therefore the required unit-bearing 9/5 coupling is still not derived from strict data.",
            "next_honest_step": "Do not replay denominator-5 coefficient boxes, imported center=9/5 constraints, endpoint pins, C12/D12 pinned-defect selectors, circulant/Fourier/irrep constructions, or generic local quadratic minimization as sourcehood.  A next proof-grade move must either introduce a new strict dimensional/unit source outside the local endpoint-coefficient family or provide an exported variational functional whose Euler/minimizer equation analytically forces 9/5 without inserting that value as a premise; otherwise preserve no-new-live-frontier.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    audit = payload["variational_unit_law_9_over_5_derivation_no_go_audit"]
    lines = [
        "# P2881/S1831 finite variational/unit-law 9/5 derivation no-go audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Variational/unit-law audit",
        f"- candidate class: `{audit['candidate_class']}`",
        f"- stencil count: `{audit['stencil_count']}`",
        f"- constraint count: `{audit['constraint_count']}`",
        f"- objective count: `{audit['objective_count']}`",
        f"- candidate count: `{audit['candidate_count']}`",
        f"- unique minimizer candidate count: `{audit['unique_minimizer_candidate_count']}`",
        f"- derives center 9/5 candidate count: `{audit['derives_center_9_over_5_candidate_count']}`",
        "",
        "## Boundary",
        payload["decision"]["reason"],
        "",
        "## Recommendation",
        payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    payload = build_payload(read_json(P2880))
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(
        STRICT_EQUATION_SHEET,
        "P2881/S1831 finite variational/unit-law 9/5 derivation no-go audit",
        "## P2881/S1831 finite variational/unit-law 9/5 derivation no-go audit\n\n"
        "`P2881/S1831` audits source-neutral local variational/unit-law candidates after `P2880`: denominator-5 radius-1 stencils under finite unit constraints, minimized against local nonnegative quadratic objectives.  Exact enumeration finds no candidate whose unique minimizer has center coefficient `9/5`; any `9/5` still must be inserted as coefficient alphabet/constraint data rather than derived.  No independent unit-bearing `9/5` coupling theorem, strict damping bridge, `L_total`, EOM, Hamiltonian, role transfer, or ToE closure follows.\n",
    )
    append_once(
        STRICT_LAGRANGIAN_DRAFT,
        "P2881/S1831 finite variational/unit-law 9/5 derivation `L_total` guard",
        "## P2881/S1831 finite variational/unit-law 9/5 derivation `L_total` guard\n\n"
        "`P2881/S1831` adds no strict action term.  The tested local quadratic variational/unit laws do not export a derived `9/5` coupling coefficient, localized unit-bearing boundary/source density, localization/pullback theorem, variational chain rule, nonproxy `L_total`, EOM, or Hamiltonian source.\n",
    )
    append_once(
        AGENTS,
        "Current finite variational/unit-law 9/5 derivation no-go guardrail (P2881/S1831, 2026-06-18)",
        "## Current finite variational/unit-law 9/5 derivation no-go guardrail (P2881/S1831, 2026-06-18)\n\n"
        "- P2881 tests the post-P2880 requirement for an independent variational/unit law deriving `9/5` from strict data, using exact finite minimization over denominator-5 radius-1 stencils and source-neutral local unit constraints/objectives.\n"
        "- No audited candidate has a unique minimizer with center coefficient `9/5`; local quadratic minimization does not derive the missing unit-bearing coupling without inserting `9/5` as premise data.\n"
        "- Do not promote denominator-5 coefficient boxes, imported center=`9/5` constraints, endpoint pins, `C12`/`D12` pinned-defect selectors, circulant/Fourier/irrep constructions, or generic local quadratic minimization to strict damping/compression bridge, selector closure, role transfer, unit-bearing `L_total`, EOM, Hamiltonian, or ToE closure.\n"
        "- A next admissible proof-grade move must introduce a new strict dimensional/unit source outside the local endpoint-coefficient family, or an exported variational functional whose Euler/minimizer equation analytically forces `9/5` without premise insertion; otherwise preserve no-new-live-frontier.\n",
    )
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
