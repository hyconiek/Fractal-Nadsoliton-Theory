#!/usr/bin/env python3
"""P2879/S1829: C12 chiral pinned-defect translation-origin no-go audit.

P2878 showed that full D12-equivariant selectors on pinned-defect density orbits
cannot export endpoint 11.  This packet tests the next honest weakening: assume
a non-D12-symmetric chiral/orientation source has broken reflections, so only the
orientation-preserving C12 rotations remain.  The question is whether this
chiral reduction can now strict-source the endpoint-11 origin/support law.

The finite answer is still no.  The audit quotients the P2877 pinned defect
12-vectors by C12, computes rotational stabilizers, and applies the same
stabilizer fixed-point criterion for equivariant endpoint selectors.  No C12
orbit forces endpoint 11 and no orbit has a unique endpoint selector: asymmetric
or singleton records have trivial stabilizer and therefore a 12-way translated
choice, while rotation-invariant records have no fixed endpoint.  Chiral
orientation can remove reflection ambiguity but does not source a translation
origin or the unit-bearing 9/5 coupling theorem.
"""
from __future__ import annotations

import json
from itertools import product
from typing import Any, Iterable

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2803_s1753_meringer_scd_import_decode_validation_certificate import sha
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN

P2878 = GEN / "p2878_s1828_d12_equivariant_pinned_defect_origin_law_no_go_audit.json"
OUT = GEN / "p2879_s1829_c12_chiral_pinned_defect_translation_origin_no_go_audit.json"
MD = GEN / "p2879_s1829_c12_chiral_pinned_defect_translation_origin_no_go_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

MODULUS = 12
ENDPOINTS = tuple(range(MODULUS))
TARGET_ENDPOINT = 11
OFFSETS = (-1, 0, 1)
COEFFICIENTS = (-1, 0, 1)
C12 = tuple(range(MODULUS))  # shift: x -> x + shift mod 12


def rotate_endpoint(shift: int, endpoint: int) -> int:
    return (endpoint + shift) % MODULUS


def rotate_vector(shift: int, vector: tuple[int, ...]) -> tuple[int, ...]:
    out = [0 for _ in ENDPOINTS]
    for endpoint, value in enumerate(vector):
        out[rotate_endpoint(shift, endpoint)] = value
    return tuple(out)


def defect_density(pin: int, stencil: tuple[int, ...]) -> tuple[int, ...]:
    values = [0 for _ in ENDPOINTS]
    for coefficient, offset in zip(stencil, OFFSETS):
        values[(pin + offset) % MODULUS] += coefficient
    return tuple(values)


def generated_vectors() -> list[tuple[int, ...]]:
    return [defect_density(pin, stencil) for pin in ENDPOINTS for stencil in product(COEFFICIENTS, repeat=len(OFFSETS))]


def support(vector: tuple[int, ...]) -> list[int]:
    return [endpoint for endpoint, value in enumerate(vector) if value != 0]


def orbit(vector: tuple[int, ...]) -> set[tuple[int, ...]]:
    return {rotate_vector(shift, vector) for shift in C12}


def stabilizer(vector: tuple[int, ...]) -> list[int]:
    return [shift for shift in C12 if rotate_vector(shift, vector) == vector]


def fixed_endpoints(shifts: Iterable[int]) -> list[int]:
    shifts = list(shifts)
    return [endpoint for endpoint in ENDPOINTS if all(rotate_endpoint(shift, endpoint) == endpoint for shift in shifts)]


def orbit_records(vectors: list[tuple[int, ...]]) -> list[dict[str, Any]]:
    unseen = set(vectors)
    records: list[dict[str, Any]] = []
    while unseen:
        seed = next(iter(unseen))
        orb = orbit(seed)
        rep = min(orb)
        stab = stabilizer(rep)
        fixed = fixed_endpoints(stab)
        records.append(
            {
                "representative": list(rep),
                "support": support(rep),
                "orbit_size": len(orb),
                "rotational_stabilizer": stab,
                "stabilizer_size": len(stab),
                "fixed_endpoint_options_for_c12_equivariant_selector": fixed,
                "unique_fixed_endpoint": fixed[0] if len(fixed) == 1 else None,
                "target_11_forced": fixed == [TARGET_ENDPOINT],
            }
        )
        unseen -= orb
    return sorted(records, key=lambda r: (r["orbit_size"], r["stabilizer_size"], r["support"], r["representative"]))


def build_payload(p2878: dict[str, Any]) -> dict[str, Any]:
    vectors = generated_vectors()
    unique_vectors = sorted(set(vectors))
    records = orbit_records(unique_vectors)
    target_forcing_records = [record for record in records if record["target_11_forced"]]
    unique_endpoint_records = [record for record in records if record["unique_fixed_endpoint"] is not None]
    singleton_support_orbits = [record for record in records if len(record["support"]) == 1]
    trivial_stabilizer_records = [record for record in records if record["stabilizer_size"] == 1]
    no_fixed_endpoint_records = [record for record in records if len(record["fixed_endpoint_options_for_c12_equivariant_selector"]) == 0]
    stabilizer_histogram: dict[str, int] = {}
    fixed_option_histogram: dict[str, int] = {}
    for record in records:
        stabilizer_histogram[str(record["stabilizer_size"])] = stabilizer_histogram.get(str(record["stabilizer_size"]), 0) + 1
        fixed_option_histogram[str(len(record["fixed_endpoint_options_for_c12_equivariant_selector"]))] = fixed_option_histogram.get(str(len(record["fixed_endpoint_options_for_c12_equivariant_selector"])), 0) + 1
    facts = {
        "p2878_rechecked": p2878.get("status") == "P2878_D12_EQUIVARIANT_PINNED_DEFECT_ORIGIN_LAW_NO_GO_AUDIT_NO_CLOSURE",
        "all_pinned_defect_vectors_generated": len(vectors) == 12 * 3**3,
        "c12_group_size_checked": len(C12) == 12,
        "no_c12_orbit_forces_endpoint_11": len(target_forcing_records) == 0,
        "no_c12_orbit_has_unique_endpoint_selector": len(unique_endpoint_records) == 0,
        "singleton_orbits_have_translation_choice_not_endpoint_11_law": all(len(record["fixed_endpoint_options_for_c12_equivariant_selector"]) == 12 for record in singleton_support_orbits),
    }
    return {
        "status": "P2879_C12_CHIRAL_PINNED_DEFECT_TRANSLATION_ORIGIN_NO_GO_AUDIT_NO_CLOSURE",
        "input_hashes": {"P2878": sha(P2878)},
        "c12_chiral_pinned_defect_translation_origin_no_go_audit": {
            "input_status_rechecked": p2878.get("status"),
            "candidate_class": "C12-equivariant endpoint selectors after a hypothetical chiral/reflection-breaking source on pinned defect density orbits",
            "modulus": MODULUS,
            "target_endpoint": TARGET_ENDPOINT,
            "raw_pinned_defect_vector_count": len(vectors),
            "unique_density_vector_count": len(unique_vectors),
            "c12_group_size": len(C12),
            "orbit_count": len(records),
            "stabilizer_size_histogram": stabilizer_histogram,
            "fixed_endpoint_option_count_histogram": fixed_option_histogram,
            "target_11_forcing_orbit_count": len(target_forcing_records),
            "unique_endpoint_selector_orbit_count": len(unique_endpoint_records),
            "singleton_support_orbit_count": len(singleton_support_orbits),
            "trivial_stabilizer_orbit_count": len(trivial_stabilizer_records),
            "no_fixed_endpoint_orbit_count": len(no_fixed_endpoint_records),
            "sample_orbit_records": records[:12],
            "proof_certificate": {
                "chiral_reduction_step": "Reflections are removed and only the orientation-preserving C12 rotation action is retained.",
                "equivariant_map_criterion": "For a C12 orbit representative x, an equivariant endpoint selector can choose only endpoints fixed by the rotational stabilizer Stab_C12(x).",
                "finite_result": "No C12 orbit has fixed endpoint set {11} and no C12 orbit has a unique fixed endpoint selector; singleton/asymmetric orbits have a 12-way translated choice.",
                "sourcehood_step": "A chiral orientation premise can break reflection but not translation origin.  Endpoint 11 and the unit-bearing 9/5 coupling remain unsourced.",
            },
        },
        "acceptance_matrix": {
            "facts": facts,
            "accepted_as_chiral_strict_origin_support_law_for_endpoint_11": False,
            "exports_translation_origin_source_law": False,
            "exports_boundary_source_law": False,
            "exports_orientation_source_law": False,
            "exports_selector_or_localizer_source": False,
            "exports_unit_bearing_9_over_5_coupling_theorem": False,
            "exports_strict_damping_compression_bridge": False,
        },
        "decision": {
            "negative_export_flags": {
                "translation_origin_source_law_exported": False,
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
            "reason": "P2879 repeats the stabilizer fixed-point audit after hypothetically breaking reflection symmetry down from D12 to C12.  No C12 orbit forces endpoint 11 and no C12 orbit has a unique endpoint selector; chiral orientation removes reflection ambiguity but leaves a translated 12-way origin choice.  Therefore no strict translation-origin source or unit-bearing 9/5 theorem is exported.",
            "next_honest_step": "Do not replay C12/chiral pinned-defect selectors, orientation-only symmetry breaking, imported translation origins, D12-equivariant selectors, pinned defects, circulant/Fourier constructions, D12 irreps/characters, or endpoint predicates as sourcehood.  A next proof-grade move must either export a strict translation-origin source with a computed unit-bearing 9/5 coupling, or pivot outside the endpoint-localizer/source family; otherwise preserve no-new-live-frontier.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    audit = payload["c12_chiral_pinned_defect_translation_origin_no_go_audit"]
    lines = [
        "# P2879/S1829 C12 chiral pinned-defect translation-origin no-go audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Chiral C12 selector audit",
        f"- candidate class: `{audit['candidate_class']}`",
        f"- raw pinned defect vector count: `{audit['raw_pinned_defect_vector_count']}`",
        f"- unique density vector count: `{audit['unique_density_vector_count']}`",
        f"- C12 orbit count: `{audit['orbit_count']}`",
        f"- target-11 forcing orbit count: `{audit['target_11_forcing_orbit_count']}`",
        f"- unique endpoint selector orbit count: `{audit['unique_endpoint_selector_orbit_count']}`",
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
    payload = build_payload(read_json(P2878))
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(
        STRICT_EQUATION_SHEET,
        "P2879/S1829 C12 chiral pinned-defect translation-origin no-go audit",
        "## P2879/S1829 C12 chiral pinned-defect translation-origin no-go audit\n\n"
        "`P2879/S1829` reruns the pinned-defect origin-law test after a hypothetical chiral/reflection-breaking source reduces `D12` to orientation-preserving `C12`.  The finite stabilizer fixed-point table still has no orbit forcing endpoint `11` and no orbit with a unique endpoint selector: singleton/asymmetric records retain a translated 12-way origin choice, while rotation-invariant records fix no endpoint.  Chiral orientation alone therefore does not export a strict translation-origin source, unit-bearing `9/5` coefficient theorem, boundary source law, selector/localizer source, strict damping bridge, `L_total`, EOM, Hamiltonian, role transfer, or ToE closure.\n",
    )
    append_once(
        STRICT_LAGRANGIAN_DRAFT,
        "P2879/S1829 C12 chiral pinned-defect translation-origin `L_total` guard",
        "## P2879/S1829 C12 chiral pinned-defect translation-origin `L_total` guard\n\n"
        "`P2879/S1829` adds no strict action term.  Reducing `D12` to chiral `C12` symmetry does not export a localized unit-bearing boundary/source density at `11`, translation-origin source law, coupling coefficient, localization/pullback theorem, variational chain rule, nonproxy `L_total`, EOM, or Hamiltonian source.\n",
    )
    append_once(
        AGENTS,
        "Current C12 chiral pinned-defect translation-origin no-go guardrail (P2879/S1829, 2026-06-18)",
        "## Current C12 chiral pinned-defect translation-origin no-go guardrail (P2879/S1829, 2026-06-18)\n\n"
        "- P2879 tests the orientation-only/chiral relaxation of P2878 by reducing the symmetry from full `D12` to rotations `C12` on the pinned-defect density family.\n"
        "- No `C12` orbit forces endpoint `11`, and no orbit has a unique endpoint selector; chiral orientation removes reflection ambiguity but not the translated 12-way origin choice.\n"
        "- Do not promote `C12`/chiral pinned-defect selectors, orientation-only symmetry breaking, imported translation origins, `D12`-equivariant selectors, pinned defects, circulant/Fourier constructions, `D12` irreps/characters, endpoint predicates, endpoint-label imports, or prime-5 scaled coefficients to strict damping/compression bridge, selector closure, role transfer, unit-bearing `L_total`, EOM, Hamiltonian, or ToE closure.\n"
        "- A next admissible proof-grade move must export a strict translation-origin source with computed unit-bearing `9/5` coupling, or pivot outside the endpoint-localizer/source family; otherwise preserve no-new-live-frontier.\n",
    )
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
