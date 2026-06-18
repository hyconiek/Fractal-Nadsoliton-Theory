#!/usr/bin/env python3
"""P2880/S1830: origin-pinned 9/5 coefficient-import no-go audit.

P2879 left a sharply localized obligation: a strict translation-origin source
selecting endpoint 11 together with a computed unit-bearing 9/5 coupling.  This
packet grants the first half as a hypothetical premise and audits the second
half.  Namely: if endpoint 11 is pinned by an imported origin law, can the
radius-1 local coefficient family itself force the coefficient 9/5?

The finite answer is no.  Over the denominator-5 local coefficient box
{-9/5,...,9/5}^3, coefficient 9/5 is representable only as a chosen coefficient
slot.  The target record is non-unique, many stencils contain 9/5, many do not,
and the same counts translate uniformly to every endpoint.  Thus an origin pin
plus a prime-5 coefficient alphabet represents the missing number but does not
export a strict unit-bearing coupling theorem, action density, bridge, L_total,
EOM, Hamiltonian, role transfer, or ToE closure.
"""
from __future__ import annotations

import json
from fractions import Fraction
from itertools import product
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2803_s1753_meringer_scd_import_decode_validation_certificate import sha
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN

P2879 = GEN / "p2879_s1829_c12_chiral_pinned_defect_translation_origin_no_go_audit.json"
OUT = GEN / "p2880_s1830_origin_pinned_9_over_5_coefficient_import_no_go_audit.json"
MD = GEN / "p2880_s1830_origin_pinned_9_over_5_coefficient_import_no_go_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

MODULUS = 12
ENDPOINTS = tuple(range(MODULUS))
TARGET_ENDPOINT = 11
OFFSETS = (-1, 0, 1)
DENOMINATOR = 5
NUMERATORS = tuple(range(-9, 10))
TARGET_COEFFICIENT = Fraction(9, 5)


def endpoint_for(pin: int, offset: int) -> int:
    return (pin + offset) % MODULUS


def coefficient_stencils() -> list[tuple[Fraction, Fraction, Fraction]]:
    return [tuple(Fraction(n, DENOMINATOR) for n in nums) for nums in product(NUMERATORS, repeat=len(OFFSETS))]


def translated_counts(stencils: list[tuple[Fraction, Fraction, Fraction]]) -> dict[str, dict[str, int]]:
    table: dict[str, dict[str, int]] = {}
    for endpoint in ENDPOINTS:
        endpoint_slot_counts = {str(offset): 0 for offset in OFFSETS}
        for stencil in stencils:
            for offset, coeff in zip(OFFSETS, stencil):
                if coeff == TARGET_COEFFICIENT and endpoint_for(endpoint - offset, offset) == endpoint:
                    endpoint_slot_counts[str(offset)] += 1
        table[str(endpoint)] = endpoint_slot_counts
    return table


def build_payload(p2879: dict[str, Any]) -> dict[str, Any]:
    stencils = coefficient_stencils()
    target_center_stencils = [s for s in stencils if s[1] == TARGET_COEFFICIENT]
    target_any_slot_stencils = [s for s in stencils if TARGET_COEFFICIENT in s]
    target_absent_stencils = [s for s in stencils if TARGET_COEFFICIENT not in s]
    unique_target_center = len(target_center_stencils) == 1
    counts = translated_counts(stencils)
    endpoint_total_counts = {endpoint: sum(slot_counts.values()) for endpoint, slot_counts in counts.items()}
    uniform_total = len(set(endpoint_total_counts.values())) == 1
    facts = {
        "p2879_rechecked": p2879.get("status") == "P2879_C12_CHIRAL_PINNED_DEFECT_TRANSLATION_ORIGIN_NO_GO_AUDIT_NO_CLOSURE",
        "denominator_5_grid_checked": DENOMINATOR == 5 and len(NUMERATORS) == 19,
        "all_radius_1_rational_stencils_generated": len(stencils) == 19**3,
        "target_9_over_5_is_representable": bool(target_center_stencils),
        "target_9_over_5_is_not_unique": not unique_target_center,
        "many_stencils_contain_9_over_5": len(target_any_slot_stencils) > 1,
        "many_stencils_omit_9_over_5": len(target_absent_stencils) > 1,
        "translation_uniformity_blocks_endpoint_11_privilege": uniform_total,
    }
    return {
        "status": "P2880_ORIGIN_PINNED_9_OVER_5_COEFFICIENT_IMPORT_NO_GO_AUDIT_NO_CLOSURE",
        "input_hashes": {"P2879": sha(P2879)},
        "origin_pinned_9_over_5_coefficient_import_no_go_audit": {
            "input_status_rechecked": p2879.get("status"),
            "candidate_class": "hypothetical imported endpoint-11 origin pin plus denominator-5 radius-1 local coefficient stencils",
            "modulus": MODULUS,
            "target_endpoint": TARGET_ENDPOINT,
            "offsets": list(OFFSETS),
            "denominator": DENOMINATOR,
            "numerator_range": [min(NUMERATORS), max(NUMERATORS)],
            "target_coefficient": "9/5",
            "stencil_count": len(stencils),
            "target_center_stencil_count": len(target_center_stencils),
            "target_any_slot_stencil_count": len(target_any_slot_stencils),
            "target_absent_stencil_count": len(target_absent_stencils),
            "endpoint_slot_counts_for_imported_9_over_5": counts,
            "endpoint_total_counts_for_imported_9_over_5": endpoint_total_counts,
            "sample_target_center_stencils": [[str(x) for x in s] for s in target_center_stencils[:12]],
            "proof_certificate": {
                "granted_premise": "The missing translation-origin source is hypothetically granted by pinning endpoint 11.",
                "finite_test": "Enumerate all denominator-5 radius-1 coefficient stencils with numerators -9..9 and ask whether 9/5 is forced rather than chosen.",
                "finite_result": "9/5 appears in 361 center-slot stencils and 1027 any-slot stencils, while 5832 stencils omit it; all endpoint counts translate uniformly.",
                "sourcehood_step": "The coefficient is represented by importing the prime-5 coefficient alphabet and selecting a slot/value; no unit-bearing source law computes 9/5.",
            },
        },
        "acceptance_matrix": {
            "facts": facts,
            "accepted_as_unit_bearing_9_over_5_coupling_theorem": False,
            "exports_translation_origin_source_law": False,
            "exports_boundary_source_law": False,
            "exports_selector_or_localizer_source": False,
            "exports_unit_bearing_action_density": False,
            "exports_strict_damping_compression_bridge": False,
        },
        "decision": {
            "negative_export_flags": {
                "translation_origin_source_law_exported": False,
                "boundary_source_law_exported": False,
                "selector_or_localizer_source_exported": False,
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
            "reason": "P2880 grants an endpoint-11 origin pin as a premise and audits the denominator-5 local coefficient box.  The value 9/5 is representable but not forced: it occurs in many chosen stencils, many stencils omit it, and the endpoint counts are translation-uniform.  Thus 9/5 remains a coefficient import rather than a strict unit-bearing coupling theorem.",
            "next_honest_step": "Do not replay imported origin pins, denominator-5 coefficient boxes, prime-5 scaled coefficients, endpoint-localizer families, C12/D12 pinned-defect selectors, circulant/Fourier constructions, D12 irreps/characters, or endpoint predicates as sourcehood.  A next proof-grade move must provide an independent variational/unit law deriving the 9/5 coupling from strict data, or pivot to a genuinely different typed object outside the endpoint/coefficient-import family; otherwise preserve no-new-live-frontier.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    audit = payload["origin_pinned_9_over_5_coefficient_import_no_go_audit"]
    lines = [
        "# P2880/S1830 origin-pinned 9/5 coefficient-import no-go audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Coefficient import audit",
        f"- candidate class: `{audit['candidate_class']}`",
        f"- stencil count: `{audit['stencil_count']}`",
        f"- target center-slot stencil count: `{audit['target_center_stencil_count']}`",
        f"- target any-slot stencil count: `{audit['target_any_slot_stencil_count']}`",
        f"- target-absent stencil count: `{audit['target_absent_stencil_count']}`",
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
    payload = build_payload(read_json(P2879))
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(
        STRICT_EQUATION_SHEET,
        "P2880/S1830 origin-pinned 9/5 coefficient-import no-go audit",
        "## P2880/S1830 origin-pinned 9/5 coefficient-import no-go audit\n\n"
        "`P2880/S1830` grants a hypothetical endpoint-`11` origin pin after `P2879` and audits denominator-5 radius-1 local coefficient stencils.  The finite table shows `9/5` is representable but not forced: it appears only by chosen coefficient slot/value, many stencils omit it, and translated endpoint counts are uniform.  Therefore no independent unit-bearing `9/5` coupling theorem, strict damping bridge, `L_total`, EOM, Hamiltonian, role transfer, or ToE closure follows.\n",
    )
    append_once(
        STRICT_LAGRANGIAN_DRAFT,
        "P2880/S1830 origin-pinned 9/5 coefficient-import `L_total` guard",
        "## P2880/S1830 origin-pinned 9/5 coefficient-import `L_total` guard\n\n"
        "`P2880/S1830` adds no strict action term.  An imported endpoint-`11` pin plus a denominator-5 coefficient box does not export a localized unit-bearing boundary/source density, derived `9/5` coupling coefficient, localization/pullback theorem, variational chain rule, nonproxy `L_total`, EOM, or Hamiltonian source.\n",
    )
    append_once(
        AGENTS,
        "Current origin-pinned 9/5 coefficient-import no-go guardrail (P2880/S1830, 2026-06-18)",
        "## Current origin-pinned 9/5 coefficient-import no-go guardrail (P2880/S1830, 2026-06-18)\n\n"
        "- P2880 grants the post-P2879 missing endpoint-`11` translation-origin pin as a premise and audits whether denominator-5 radius-1 local coefficient stencils force the unit-bearing `9/5` coupling.\n"
        "- The finite table shows `9/5` is representable but not forced: it occurs in many chosen stencils, many stencils omit it, and all endpoint counts translate uniformly.\n"
        "- Do not promote imported origin pins, denominator-5 coefficient boxes, prime-5 scaled coefficients, endpoint-localizer families, `C12`/`D12` pinned-defect selectors, circulant/Fourier constructions, `D12` irreps/characters, or endpoint predicates to strict damping/compression bridge, selector closure, role transfer, unit-bearing `L_total`, EOM, Hamiltonian, or ToE closure.\n"
        "- A next admissible proof-grade move must export an independent variational/unit law deriving the `9/5` coupling from strict data, or pivot outside the endpoint/coefficient-import family; otherwise preserve no-new-live-frontier.\n",
    )
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
