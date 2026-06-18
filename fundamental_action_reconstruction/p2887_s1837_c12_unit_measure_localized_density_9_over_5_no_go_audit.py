#!/usr/bin/env python3
"""P2887/S1837: C12 unit-measure/localized-density 9/5 no-go audit.

P2886 found no already-exported external unit measure or localized action density.
This packet tests a concrete candidate rather than replaying inventory: can a
C12-neutral unit measure and a finite localized action density on Z12 produce a
nonzero 9/5 coupling without importing an endpoint or coefficient?

The finite calculation enumerates all binary supports on Z12 and all ternary
integer densities {-1,0,1}^12.  Rotation-invariant supports are only empty/full,
and rotation-invariant ternary densities are only constant.  Therefore a
C12-neutral unit measure/density can be global but not localized; singleton
localization exists only by choosing an endpoint representative.  The uniform
unit-measure integral of C12-invariant ternary densities is -1, 0, or 1, never
9/5.  Thus this natural external unit-measure/local-density candidate does not
export the missing unit-bearing 9/5 source.
"""
from __future__ import annotations

import itertools
import json
from fractions import Fraction
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2803_s1753_meringer_scd_import_decode_validation_certificate import sha
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN

P2886 = GEN / "p2886_s1836_external_unit_measure_action_density_inventory_no_go_audit.json"
OUT = GEN / "p2887_s1837_c12_unit_measure_localized_density_9_over_5_no_go_audit.json"
MD = GEN / "p2887_s1837_c12_unit_measure_localized_density_9_over_5_no_go_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

N = 12
TARGET = Fraction(9, 5)
TERNARY_VALUES = (-1, 0, 1)


def rotate_tuple(values: tuple[int, ...], shift: int) -> tuple[int, ...]:
    shift %= len(values)
    return values[shift:] + values[:shift]


def is_c12_invariant(values: tuple[int, ...]) -> bool:
    return all(rotate_tuple(values, shift) == values for shift in range(N))


def support_size(values: tuple[int, ...]) -> int:
    return sum(1 for value in values if value != 0)


def binary_support_records() -> dict[str, Any]:
    support_count = 0
    invariant_supports: list[tuple[int, ...]] = []
    singleton_supports = 0
    invariant_singleton_supports = 0
    for bits in itertools.product((0, 1), repeat=N):
        support_count += 1
        if sum(bits) == 1:
            singleton_supports += 1
        if is_c12_invariant(bits):
            invariant_supports.append(bits)
            if sum(bits) == 1:
                invariant_singleton_supports += 1
    return {
        "binary_support_count": support_count,
        "c12_invariant_support_count": len(invariant_supports),
        "c12_invariant_support_sizes": [sum(bits) for bits in invariant_supports],
        "singleton_support_count": singleton_supports,
        "c12_invariant_singleton_support_count": invariant_singleton_supports,
        "sample_c12_invariant_supports": [list(bits) for bits in invariant_supports],
    }


def ternary_density_records() -> dict[str, Any]:
    density_count = 0
    invariant_densities: list[tuple[int, ...]] = []
    singleton_densities = 0
    invariant_singleton_densities = 0
    invariant_integrals: list[Fraction] = []
    for density in itertools.product(TERNARY_VALUES, repeat=N):
        density_count += 1
        size = support_size(density)
        if size == 1:
            singleton_densities += 1
        if is_c12_invariant(density):
            invariant_densities.append(density)
            invariant_integrals.append(Fraction(sum(density), N))
            if size == 1:
                invariant_singleton_densities += 1
    target_integrals = [value for value in invariant_integrals if value == TARGET]
    return {
        "ternary_density_count": density_count,
        "c12_invariant_density_count": len(invariant_densities),
        "c12_invariant_density_values": [list(density) for density in invariant_densities],
        "singleton_density_count": singleton_densities,
        "c12_invariant_singleton_density_count": invariant_singleton_densities,
        "uniform_unit_measure_integrals_of_invariant_densities": [str(value) for value in invariant_integrals],
        "target_9_over_5_integral_count": len(target_integrals),
    }


def build_payload(p2886: dict[str, Any]) -> dict[str, Any]:
    support_audit = binary_support_records()
    density_audit = ternary_density_records()
    facts = {
        "p2886_rechecked": p2886.get("status") == "P2886_EXTERNAL_UNIT_MEASURE_ACTION_DENSITY_INVENTORY_NO_GO_AUDIT_NO_CLOSURE",
        "all_binary_supports_checked": support_audit["binary_support_count"] == 2**N,
        "all_ternary_densities_checked": density_audit["ternary_density_count"] == 3**N,
        "c12_neutral_supports_are_only_empty_or_full": support_audit["c12_invariant_support_sizes"] == [0, N],
        "no_c12_invariant_singleton_localization": support_audit["c12_invariant_singleton_support_count"] == 0 and density_audit["c12_invariant_singleton_density_count"] == 0,
        "uniform_c12_invariant_density_integral_never_9_over_5": density_audit["target_9_over_5_integral_count"] == 0,
    }
    return {
        "status": "P2887_C12_UNIT_MEASURE_LOCALIZED_DENSITY_9_OVER_5_NO_GO_AUDIT_NO_CLOSURE",
        "input_hashes": {"P2886": sha(P2886)},
        "c12_unit_measure_localized_density_9_over_5_no_go_audit": {
            "input_status_rechecked": p2886.get("status"),
            "candidate_class": "C12-neutral unit measure/localized action density on Z12 with ternary integer density values {-1,0,1}",
            "support_audit": support_audit,
            "density_audit": density_audit,
            "proof_certificate": {
                "support_result": "The only C12-invariant binary supports on Z12 are empty and full; every singleton support imports an endpoint representative.",
                "density_result": "The only C12-invariant ternary densities are constant -1, 0, and +1; none is localized.",
                "coupling_result": "With the C12-neutral uniform unit measure, invariant ternary density integrals are -1, 0, and +1, never 9/5.",
                "sourcehood_step": "A strict unit-bearing 9/5 action density therefore needs a new non-C12-neutral origin/source law, a new coefficient theorem, or a different typed object; this candidate does not supply it.",
            },
        },
        "acceptance_matrix": {
            "facts": facts,
            "accepted_as_strict_unit_measure_localized_density_source": False,
            "exports_nonimported_endpoint_localization": False,
            "exports_unit_bearing_9_over_5_coupling": False,
            "exports_variational_chain_rule_to_ltotal": False,
            "exports_strict_damping_compression_bridge": False,
        },
        "decision": {
            "negative_export_flags": {
                "strict_unit_measure_localized_density_exported": False,
                "nonimported_endpoint_localization_exported": False,
                "unit_bearing_9_over_5_coupling_exported": False,
                "localized_action_density_exported": False,
                "variational_chain_rule_exported": False,
                "strict_damping_compression_bridge_exported": False,
                "full_kernel_bridge_exported": False,
                "role_transfer_exported": False,
                "ltotal_exported": False,
                "eom_closure_exported": False,
                "hamiltonian_closure_exported": False,
                "toe_closure_exported": False,
            },
            "reason": "P2887 tests a concrete post-P2886 candidate: C12-neutral unit measure plus localized action density on Z12.  Exhaustive support and ternary-density enumeration shows C12 neutrality permits only empty/full supports and constant densities; singleton localization imports an endpoint, and the uniform unit-measure integrals are -1, 0, and 1 rather than 9/5.",
            "next_honest_step": "Do not replay C12-neutral unit measures, singleton endpoint imports, invariant-count actions, ratio algebra, or scalar Euler transmission.  A next proof-grade move must provide one new non-C12-neutral strict origin/source law with a computed unit-bearing 9/5 coupling, or pivot to a genuinely different typed object; otherwise preserve no-new-live-frontier.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    audit = payload["c12_unit_measure_localized_density_9_over_5_no_go_audit"]
    support = audit["support_audit"]
    density = audit["density_audit"]
    lines = [
        "# P2887/S1837 C12 unit-measure/localized-density 9/5 no-go audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Finite support/density audit",
        f"- binary support count: `{support['binary_support_count']}`",
        f"- C12-invariant support count: `{support['c12_invariant_support_count']}`",
        f"- C12-invariant support sizes: `{support['c12_invariant_support_sizes']}`",
        f"- ternary density count: `{density['ternary_density_count']}`",
        f"- C12-invariant density count: `{density['c12_invariant_density_count']}`",
        f"- uniform integrals of invariant densities: `{density['uniform_unit_measure_integrals_of_invariant_densities']}`",
        f"- target 9/5 integral count: `{density['target_9_over_5_integral_count']}`",
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
    payload = build_payload(read_json(P2886))
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(
        STRICT_EQUATION_SHEET,
        "P2887/S1837 C12 unit-measure/localized-density 9/5 no-go audit",
        "## P2887/S1837 C12 unit-measure/localized-density 9/5 no-go audit\n\n"
        "`P2887/S1837` tests a concrete post-`P2886` candidate: a `C12`-neutral unit measure and localized action density on `Z12`.  Exhaustive enumeration of all `2^12` supports and all `3^12` ternary densities shows the only `C12`-invariant supports are empty/full and the only invariant ternary densities are constant; singleton localization imports an endpoint representative, and uniform unit-measure integrals are `-1`, `0`, and `1`, never `9/5`.  No localized unit-bearing action density, strict damping bridge, `L_total`, EOM, Hamiltonian, role transfer, or ToE closure follows.\n",
    )
    append_once(
        STRICT_LAGRANGIAN_DRAFT,
        "P2887/S1837 C12 unit-measure/localized-density `L_total` guard",
        "## P2887/S1837 C12 unit-measure/localized-density `L_total` guard\n\n"
        "`P2887/S1837` adds no strict action term.  `C12`-neutral supports/densities are global rather than localized, and their uniform integrals do not produce `9/5`; the packet exports no localized unit-bearing action density, variational chain rule into nonproxy `L_total`, EOM, or Hamiltonian source.\n",
    )
    append_once(
        AGENTS,
        "Current C12 unit-measure/localized-density 9/5 no-go guardrail (P2887/S1837, 2026-06-18)",
        "## Current C12 unit-measure/localized-density 9/5 no-go guardrail (P2887/S1837, 2026-06-18)\n\n"
        "- P2887 exhaustively enumerates all binary supports and ternary `{-1,0,1}` densities on `Z12` as a concrete post-P2886 unit-measure/localized-action-density candidate.\n"
        "- `C12` neutrality leaves only empty/full supports and constant ternary densities; singleton localization imports an endpoint, and uniform unit-measure integrals are `-1`, `0`, and `1`, not `9/5`.\n"
        "- Do not promote `C12`-neutral unit measures, singleton endpoint imports, invariant-count actions, ratio algebra, or scalar Euler transmission to strict damping/compression bridge, selector closure, role transfer, unit-bearing `L_total`, EOM, Hamiltonian, or ToE closure.\n"
        "- A next admissible proof-grade move must provide one new non-`C12`-neutral strict origin/source law with computed unit-bearing `9/5` coupling, pivot to a genuinely different typed object, or preserve no-new-live-frontier.\n",
    )
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
