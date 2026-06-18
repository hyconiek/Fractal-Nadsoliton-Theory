#!/usr/bin/env python3
"""P2884/S1834: finite invariant-ratio 9:5 source-law no-go audit.

P2883 found no already-exported dimensional/unit source fixing the primitive
Euler source/stiffness ratio 9:5.  This packet tests the honest next candidate:
can the small source-neutral invariant vocabulary already present around Z12/C12
/D12 force a primitive numerator/denominator pair 9:5?

The audit enumerates a finite expression family over basic group/count invariants
without installing a preferred endpoint, coefficient slot, Euler source, or
literal 9:5 ratio as a law.  It records that 9 and 5 are representable in many
ways, but no expression pair is selected by the invariant vocabulary itself.
Thus this candidate supplies representability/coordinate algebra, not a strict
unit-bearing source law.
"""
from __future__ import annotations

import json
from math import gcd, lcm
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2803_s1753_meringer_scd_import_decode_validation_certificate import sha
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN

P2883 = GEN / "p2883_s1833_dimensional_unit_source_inventory_no_go_audit.json"
OUT = GEN / "p2884_s1834_finite_invariant_ratio_9_to_5_source_law_no_go_audit.json"
MD = GEN / "p2884_s1834_finite_invariant_ratio_9_to_5_source_law_no_go_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

MAX_VALUE = 72
BASE_INVARIANTS: dict[str, int] = {
    "|Z12|": 12,
    "|Aut(Z12)|": 4,
    "|D12_reflections|": 12,
    "|D12_rotations|": 12,
    "|proper_divisors_12|": 4,
    "|divisors_12|": 6,
    "|units_mod_12|": 4,
    "|nonunits_mod_12|": 8,
    "|order2_elements_Z12|": 1,
    "|solutions_x2_eq_1_mod12|": 4,
    "cycle_rank_Z12": 1,
    "boundary_rank_Z12": 11,
}


def primitive_ratio(num: int, den: int) -> tuple[int, int]:
    divisor = gcd(abs(num), abs(den))
    return (num // divisor, den // divisor)


def expression_values() -> dict[str, int]:
    """Build a bounded, source-neutral expression vocabulary from count data."""
    expressions: dict[str, int] = dict(BASE_INVARIANTS)
    items = list(BASE_INVARIANTS.items())
    for left_name, left in items:
        for right_name, right in items:
            candidates = {
                f"({left_name}+{right_name})": left + right,
                f"abs({left_name}-{right_name})": abs(left - right),
                f"gcd({left_name},{right_name})": gcd(left, right),
                f"lcm({left_name},{right_name})": lcm(left, right),
            }
            product = left * right
            if product <= MAX_VALUE:
                candidates[f"({left_name}*{right_name})"] = product
            for name, value in candidates.items():
                if 0 < value <= MAX_VALUE:
                    expressions.setdefault(name, value)
    return expressions


def ratio_records(expressions: dict[str, int]) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    for num_name, num in expressions.items():
        for den_name, den in expressions.items():
            primitive = primitive_ratio(num, den)
            records.append(
                {
                    "numerator_expression": num_name,
                    "numerator_value": num,
                    "denominator_expression": den_name,
                    "denominator_value": den,
                    "primitive_ratio": list(primitive),
                    "is_9_to_5": primitive == (9, 5),
                }
            )
    return records


def build_payload(p2883: dict[str, Any]) -> dict[str, Any]:
    expressions = expression_values()
    records = ratio_records(expressions)
    target_records = [record for record in records if record["is_9_to_5"]]
    value_histogram: dict[str, int] = {}
    for value in expressions.values():
        value_histogram[str(value)] = value_histogram.get(str(value), 0) + 1
    facts = {
        "p2883_rechecked": p2883.get("status") == "P2883_DIMENSIONAL_UNIT_SOURCE_INVENTORY_NO_GO_AUDIT_NO_CLOSURE",
        "finite_invariant_expression_family_nonempty": len(expressions) > 0,
        "target_ratio_representable_in_invariant_algebra": len(target_records) > 0,
        "target_ratio_not_unique_or_selected": len(target_records) != 1,
        "no_exported_selector_for_target_pair": True,
        "no_unit_dimension_attached_to_ratio": True,
    }
    return {
        "status": "P2884_FINITE_INVARIANT_RATIO_9_TO_5_SOURCE_LAW_NO_GO_AUDIT_NO_CLOSURE",
        "input_hashes": {"P2883": sha(P2883)},
        "finite_invariant_ratio_9_to_5_source_law_no_go_audit": {
            "input_status_rechecked": p2883.get("status"),
            "candidate_class": "bounded source-neutral Z12/C12/D12 count-invariant expressions and ordered primitive ratios",
            "base_invariants": BASE_INVARIANTS,
            "max_expression_value": MAX_VALUE,
            "expression_count": len(expressions),
            "ratio_candidate_count": len(records),
            "value_histogram_sample": dict(sorted(value_histogram.items(), key=lambda kv: int(kv[0]))[:20]),
            "target_ratio": "9:5",
            "target_ratio_record_count": len(target_records),
            "sample_target_records": target_records[:20],
            "proof_certificate": {
                "finite_rule": "Generate count expressions from Z12/C12/D12 source-neutral invariants using +, absolute difference, gcd, lcm, and bounded products; then test all ordered primitive ratios.",
                "finite_result": "The ratio 9:5 is representable in the expression algebra, but it is not unique and the vocabulary exports no selector choosing one numerator expression and one denominator expression as source/stiffness law.",
                "sourcehood_step": "A strict derivation still needs an independent unit/dimensional law selecting the 9:5 source/stiffness pair and coupling it to the action; invariant-count representability alone is not sourcehood.",
            },
        },
        "acceptance_matrix": {
            "facts": facts,
            "accepted_as_strict_invariant_source_law_for_9_to_5": False,
            "exports_source_stiffness_ratio_9_to_5": False,
            "exports_unit_bearing_9_over_5_coupling_theorem": False,
            "exports_unit_bearing_action_density": False,
            "exports_strict_damping_compression_bridge": False,
        },
        "decision": {
            "negative_export_flags": {
                "invariant_count_source_law_exported": False,
                "source_stiffness_ratio_9_to_5_exported": False,
                "unit_bearing_9_over_5_coupling_theorem_exported": False,
                "unit_bearing_action_density_exported": False,
                "strict_damping_compression_bridge_exported": False,
                "full_kernel_bridge_exported": False,
                "role_transfer_exported": False,
                "ltotal_exported": False,
                "eom_closure_exported": False,
                "hamiltonian_closure_exported": False,
                "toe_closure_exported": False,
            },
            "reason": "P2884 tests whether finite source-neutral Z12/C12/D12 count-invariant algebra can supply the missing P2883 strict dimensional/unit source for primitive ratio 9:5.  The ratio is representable, but only as one of many ordered expression-pair choices; no invariant selector, unit dimension, action density, or source/stiffness law is exported.",
            "next_honest_step": "Do not replay invariant-count expression algebra, endpoint pins, denominator-5 boxes, local quadratic minimization, scalar Euler ratio transmission, or generated-artifact inventory as sourcehood.  A next proof-grade move must either provide an explicit unit/dimensional action functional whose Euler equation analytically selects the 9:5 source/stiffness pair, or pivot to a genuinely different typed object outside the endpoint/coefficient/source-ratio/invariant-count family; otherwise preserve no-new-live-frontier.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    audit = payload["finite_invariant_ratio_9_to_5_source_law_no_go_audit"]
    lines = [
        "# P2884/S1834 finite invariant-ratio 9:5 source-law no-go audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Finite invariant-ratio audit",
        f"- candidate class: `{audit['candidate_class']}`",
        f"- expression count: `{audit['expression_count']}`",
        f"- ratio candidate count: `{audit['ratio_candidate_count']}`",
        f"- target ratio record count: `{audit['target_ratio_record_count']}`",
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
    payload = build_payload(read_json(P2883))
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(
        STRICT_EQUATION_SHEET,
        "P2884/S1834 finite invariant-ratio 9:5 source-law no-go audit",
        "## P2884/S1834 finite invariant-ratio 9:5 source-law no-go audit\n\n"
        "`P2884/S1834` tests whether bounded source-neutral `Z12/C12/D12` count-invariant algebra can supply the missing strict dimensional/unit source for primitive ratio `9:5`.  Exact finite enumeration shows `9:5` is representable, but only as one of many expression-pair choices; no invariant selector, unit dimension, action density, or source/stiffness law is exported.  No strict damping bridge, `L_total`, EOM, Hamiltonian, role transfer, or ToE closure follows.\n",
    )
    append_once(
        STRICT_LAGRANGIAN_DRAFT,
        "P2884/S1834 finite invariant-ratio source-law `L_total` guard",
        "## P2884/S1834 finite invariant-ratio source-law `L_total` guard\n\n"
        "`P2884/S1834` adds no strict action term.  Finite invariant-count representability of primitive ratio `9:5` does not export a unit-bearing action density, a selected source/stiffness law, localized boundary/source density, nonproxy `L_total`, EOM, or Hamiltonian source.\n",
    )
    append_once(
        AGENTS,
        "Current finite invariant-ratio 9:5 source-law no-go guardrail (P2884/S1834, 2026-06-18)",
        "## Current finite invariant-ratio 9:5 source-law no-go guardrail (P2884/S1834, 2026-06-18)\n\n"
        "- P2884 tests bounded source-neutral `Z12/C12/D12` count-invariant expression algebra as a candidate source for the primitive source/stiffness ratio `9:5`.\n"
        "- The ratio `9:5` is representable, but only as one of many expression-pair choices; current artifacts export no invariant selector, unit dimension, action density, or source/stiffness law selecting it.\n"
        "- Do not promote invariant-count expression algebra, endpoint pins, denominator-5 coefficient boxes, local quadratic minimization, scalar Euler ratio transmission, or generated-artifact inventory to strict damping/compression bridge, selector closure, role transfer, unit-bearing `L_total`, EOM, Hamiltonian, or ToE closure.\n"
        "- A next admissible proof-grade move must provide an explicit unit/dimensional action functional analytically selecting the `9:5` source/stiffness pair, pivot to a genuinely different typed object outside this family, or preserve no-new-live-frontier.\n",
    )
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
