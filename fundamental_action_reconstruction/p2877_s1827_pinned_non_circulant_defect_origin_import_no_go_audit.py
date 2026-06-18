#!/usr/bin/env python3
"""P2877/S1827: pinned non-circulant defect origin-import no-go audit.

P2876 showed that translation-compatible local circulant operators cannot compute
singleton endpoint 11 from intrinsic Fourier inputs.  This packet tests the next
candidate class named by that boundary: non-circulant local defects.

The audit enumerates all ternary radius-1 pinned defect densities on Z12: choose
a pin p in Z12 and coefficients on offsets {-1,0,+1}.  These defects can localize
only because p is already supplied.  The singleton endpoint-11 witnesses occur only when the imported pin is 11 or
an imported neighboring pin 0/10 is combined with an offset-only stencil;
translating the imported pin-neighborhood translates the singleton.  Thus the
construction demonstrates representability by an imported pin/offset relation,
not a strict non-premise origin/support law or a unit-bearing 9/5 coupling
theorem.
"""
from __future__ import annotations

import json
from itertools import product
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2803_s1753_meringer_scd_import_decode_validation_certificate import sha
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN

P2876 = GEN / "p2876_s1826_local_circulant_source_operator_endpoint11_no_go_audit.json"
OUT = GEN / "p2877_s1827_pinned_non_circulant_defect_origin_import_no_go_audit.json"
MD = GEN / "p2877_s1827_pinned_non_circulant_defect_origin_import_no_go_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

MODULUS = 12
ENDPOINTS = tuple(range(MODULUS))
TARGET_ENDPOINT = 11
OFFSETS = (-1, 0, 1)
COEFFICIENTS = (-1, 0, 1)


def defect_density(pin: int, stencil: tuple[int, ...]) -> list[int]:
    values = [0 for _ in ENDPOINTS]
    for coefficient, offset in zip(stencil, OFFSETS):
        values[(pin + offset) % MODULUS] += coefficient
    return values


def support(values: list[int]) -> list[int]:
    return [endpoint for endpoint, value in enumerate(values) if value != 0]


def all_candidates() -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    for pin in ENDPOINTS:
        for stencil in product(COEFFICIENTS, repeat=len(OFFSETS)):
            values = defect_density(pin, stencil)
            records.append({"pin": pin, "stencil": list(stencil), "support": support(values), "values": values})
    return records


def classify(records: list[dict[str, Any]]) -> dict[str, Any]:
    singleton_records = [r for r in records if len(r["support"]) == 1]
    singleton_11_records = [r for r in singleton_records if r["support"] == [TARGET_ENDPOINT]]
    by_pin = {str(pin): 0 for pin in ENDPOINTS}
    for r in singleton_records:
        by_pin[str(r["pin"])] += 1
    support_histogram: dict[str, int] = {}
    for r in records:
        key = str(len(r["support"]))
        support_histogram[key] = support_histogram.get(key, 0) + 1
    return {
        "candidate_count": len(records),
        "pin_count": len(ENDPOINTS),
        "stencil_count_per_pin": len(COEFFICIENTS) ** len(OFFSETS),
        "singleton_record_count": len(singleton_records),
        "singleton_11_record_count": len(singleton_11_records),
        "singleton_count_by_pin": by_pin,
        "support_size_histogram": support_histogram,
        "singleton_11_records": singleton_11_records,
        "all_singletons_have_one_nonzero_offset_coefficient": all(sum(1 for value in r["stencil"] if value != 0) == 1 for r in singleton_records),
    }


def build_payload(p2876: dict[str, Any]) -> dict[str, Any]:
    records = all_candidates()
    summary = classify(records)
    facts = {
        "p2876_rechecked": p2876.get("status") == "P2876_LOCAL_CIRCULANT_SOURCE_OPERATOR_ENDPOINT11_NO_GO_AUDIT_NO_CLOSURE",
        "all_radius1_pinned_ternary_defects_enumerated": summary["candidate_count"] == 12 * 3**3,
        "singleton_11_representable_with_imported_pin_neighborhood": summary["singleton_11_record_count"] == 6,
        "singleton_witnesses_are_uniform_across_pins": len(set(summary["singleton_count_by_pin"].values())) == 1,
        "all_singletons_are_imported_pin_neighborhood_translates": summary["singleton_record_count"] == 12 * 6,
        "no_independent_origin_law_exported": True,
    }
    return {
        "status": "P2877_PINNED_NON_CIRCULANT_DEFECT_ORIGIN_IMPORT_NO_GO_AUDIT_NO_CLOSURE",
        "input_hashes": {"P2876": sha(P2876)},
        "pinned_non_circulant_defect_origin_import_no_go_audit": {
            "input_status_rechecked": p2876.get("status"),
            "candidate_class": "ternary radius-1 pinned non-circulant defect densities on Z12",
            "modulus": MODULUS,
            "target_endpoint": TARGET_ENDPOINT,
            "offsets": list(OFFSETS),
            "coefficient_alphabet": list(COEFFICIENTS),
            "summary": summary,
            "proof_certificate": {
                "enumeration_step": "All 12*3^3 pinned radius-1 ternary defect densities are enumerated.",
                "localization_step": "Singleton support at endpoint 11 occurs exactly when an imported pin-neighborhood (pins 0, 10, or 11 with one nonzero offset coefficient) already names endpoint 11.",
                "orbit_step": "The same six singleton-neighborhood witnesses occur for every output endpoint; translating the imported pin-neighborhood translates singleton support.",
                "sourcehood_step": "The family represents endpoint 11 only after importing a non-circulant pin/origin plus offset relation that already names 11.  No strict non-premise origin/support law or unit-bearing 9/5 coefficient theorem is exported.",
            },
        },
        "acceptance_matrix": {
            "facts": facts,
            "accepted_as_strict_non_circulant_local_endpoint_11_defect_source": False,
            "exports_independent_origin_support_law": False,
            "exports_boundary_source_law": False,
            "exports_orientation_source_law": False,
            "exports_selector_or_localizer_source": False,
            "exports_unit_bearing_9_over_5_coupling_theorem": False,
            "exports_strict_damping_compression_bridge": False,
        },
        "decision": {
            "negative_export_flags": {
                "independent_origin_support_law_exported": False,
                "boundary_source_law_exported": False,
                "orientation_source_law_exported": False,
                "selector_or_localizer_source_exported": False,
                "unit_bearing_9_over_5_coupling_theorem_exported": False,
                "strict_damping_compression_bridge_exported": False,
                "full_kernel_bridge_exported": False,
                "role_transfer_exported": False,
                "ltotal_exported": False,
                "eom_closure_exported": False,
                "hamiltonian_closure_exported": False,
                "toe_closure_exported": False,
            },
            "reason": "P2877 enumerates all 324 ternary radius-1 pinned non-circulant defect densities.  Singleton endpoint 11 is representable only when an imported pin-neighborhood already names 11; the same six singleton-neighborhood witnesses occur for every output endpoint.  Thus endpoint 11 is not sourced by a strict law, and no unit-bearing 9/5 coupling theorem is exported.",
            "next_honest_step": "Do not replay pinned non-circulant defects, imported origin pins, center-only delta defects, circulant stencils, Fourier reconstruction, D12 irreps/characters, or endpoint predicates as sourcehood.  A next proof-grade move must provide an independently exported strict origin/support law selecting endpoint 11 plus a unit-bearing 9/5 coefficient theorem, or pivot to a genuinely different typed object; otherwise preserve no-new-live-frontier.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    audit = payload["pinned_non_circulant_defect_origin_import_no_go_audit"]
    summary = audit["summary"]
    lines = [
        "# P2877/S1827 pinned non-circulant defect origin-import no-go audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Pinned defect audit",
        f"- candidate class: `{audit['candidate_class']}`",
        f"- candidate count: `{summary['candidate_count']}`",
        f"- singleton endpoint-11 records: `{summary['singleton_11_record_count']}`",
        f"- singleton count by imported pin: `{summary['singleton_count_by_pin']}`",
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
    payload = build_payload(read_json(P2876))
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(
        STRICT_EQUATION_SHEET,
        "P2877/S1827 pinned non-circulant defect origin-import no-go audit",
        "## P2877/S1827 pinned non-circulant defect origin-import no-go audit\n\n"
        "`P2877/S1827` enumerates all `12*3^3=324` ternary radius-1 pinned non-circulant defect densities on `Z12`.  Singleton endpoint `11` is representable exactly when an imported pin-neighborhood already names `11`; the same six singleton-neighborhood witnesses occur for every output endpoint.  This is imported origin/offset representability, not an independently exported strict origin/support law.  No unit-bearing `9/5` coefficient theorem, boundary/orientation source law, selector/localizer source, strict damping bridge, `L_total`, EOM, Hamiltonian, role transfer, or ToE closure follows.\n",
    )
    append_once(
        STRICT_LAGRANGIAN_DRAFT,
        "P2877/S1827 pinned non-circulant defect `L_total` guard",
        "## P2877/S1827 pinned non-circulant defect `L_total` guard\n\n"
        "`P2877/S1827` adds no strict action term.  Pinned non-circulant defect densities can localize endpoint `11` only by importing a pin/origin plus offset relation that already names `11`; they do not export a localized unit-bearing boundary/source density, non-premise origin/support law, coupling coefficient, localization/pullback theorem, variational chain rule, nonproxy `L_total`, EOM, or Hamiltonian source.\n",
    )
    append_once(
        AGENTS,
        "Current pinned non-circulant defect origin-import no-go guardrail (P2877/S1827, 2026-06-18)",
        "## Current pinned non-circulant defect origin-import no-go guardrail (P2877/S1827, 2026-06-18)\n\n"
        "- P2877 enumerates all ternary radius-1 pinned non-circulant defect densities on `Z12` after P2876.\n"
        "- Singleton endpoint `11` appears only when an imported pin-neighborhood already names `11`; identical singleton-neighborhood witnesses exist for every output endpoint.\n"
        "- Do not promote pinned non-circulant defects, imported origin pins, center-only delta defects, circulant stencils, Fourier reconstruction, `D12` irreps/characters, endpoint predicates, endpoint-label imports, or prime-5 scaled coefficients to strict damping/compression bridge, selector closure, role transfer, unit-bearing `L_total`, EOM, Hamiltonian, or ToE closure.\n"
        "- A next admissible proof-grade move must export an independent strict origin/support law selecting endpoint `11` plus a unit-bearing `9/5` coefficient theorem, or use a genuinely different typed object; otherwise preserve no-new-live-frontier.\n",
    )
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
