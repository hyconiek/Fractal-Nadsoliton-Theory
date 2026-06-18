#!/usr/bin/env python3
"""P2878/S1828: D12-equivariant pinned-defect origin-law no-go audit.

P2877 showed that pinned non-circulant defects localize endpoint 11 only after
an imported pin-neighborhood already names 11.  This packet tests the missing
premise more directly: can the finite family of pinned radius-1 defect densities
itself export an intrinsic D12-equivariant origin/support law into endpoints?

For every ternary radius-1 pinned defect density on Z12, the audit forms its
12-vector, quotients by the D12 action, computes the stabilizer of each orbit
representative, and applies the finite equivariant-map criterion: an equivariant
selector from an orbit to endpoints can choose only endpoints fixed by the
representative stabilizer.  No orbit has a unique fixed endpoint 11; delta-like
singletons have two reflection-fixed endpoints, while asymmetric records have
trivial stabilizer and therefore arbitrary 12-way choice.  Hence no independent
strict origin/support law or unit-bearing 9/5 theorem is exported.
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

P2877 = GEN / "p2877_s1827_pinned_non_circulant_defect_origin_import_no_go_audit.json"
OUT = GEN / "p2878_s1828_d12_equivariant_pinned_defect_origin_law_no_go_audit.json"
MD = GEN / "p2878_s1828_d12_equivariant_pinned_defect_origin_law_no_go_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

MODULUS = 12
ENDPOINTS = tuple(range(MODULUS))
TARGET_ENDPOINT = 11
OFFSETS = (-1, 0, 1)
COEFFICIENTS = (-1, 0, 1)
GroupElement = tuple[int, int]  # (sign, shift): x -> sign*x + shift mod 12, sign in {1,-1}
D12: tuple[GroupElement, ...] = tuple((sign, shift) for sign in (1, -1) for shift in ENDPOINTS)


def apply_group_to_endpoint(element: GroupElement, endpoint: int) -> int:
    sign, shift = element
    return (sign * endpoint + shift) % MODULUS


def apply_group_to_vector(element: GroupElement, vector: tuple[int, ...]) -> tuple[int, ...]:
    out = [0 for _ in ENDPOINTS]
    for endpoint, value in enumerate(vector):
        out[apply_group_to_endpoint(element, endpoint)] = value
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
    return {apply_group_to_vector(element, vector) for element in D12}


def stabilizer(vector: tuple[int, ...]) -> list[GroupElement]:
    return [element for element in D12 if apply_group_to_vector(element, vector) == vector]


def fixed_endpoints(elements: Iterable[GroupElement]) -> list[int]:
    elements = list(elements)
    return [endpoint for endpoint in ENDPOINTS if all(apply_group_to_endpoint(element, endpoint) == endpoint for element in elements)]


def canonical_representative(vectors: set[tuple[int, ...]]) -> tuple[int, ...]:
    return min(vectors)


def orbit_records(vectors: list[tuple[int, ...]]) -> list[dict[str, Any]]:
    unseen = set(vectors)
    records: list[dict[str, Any]] = []
    while unseen:
        seed = next(iter(unseen))
        orb = orbit(seed)
        rep = canonical_representative(orb)
        stab = stabilizer(rep)
        fixed = fixed_endpoints(stab)
        records.append(
            {
                "representative": list(rep),
                "support": support(rep),
                "orbit_size": len(orb),
                "stabilizer_size": len(stab),
                "fixed_endpoint_options_for_equivariant_selector": fixed,
                "unique_fixed_endpoint": fixed[0] if len(fixed) == 1 else None,
                "target_11_forced": fixed == [TARGET_ENDPOINT],
            }
        )
        unseen -= orb
    return sorted(records, key=lambda r: (r["orbit_size"], r["stabilizer_size"], r["support"], r["representative"]))


def build_payload(p2877: dict[str, Any]) -> dict[str, Any]:
    vectors = generated_vectors()
    unique_vectors = sorted(set(vectors))
    records = orbit_records(unique_vectors)
    target_forcing_records = [record for record in records if record["target_11_forced"]]
    unique_endpoint_records = [record for record in records if record["unique_fixed_endpoint"] is not None]
    singleton_support_orbits = [record for record in records if len(record["support"]) == 1]
    stabilizer_histogram: dict[str, int] = {}
    fixed_option_histogram: dict[str, int] = {}
    for record in records:
        stabilizer_histogram[str(record["stabilizer_size"])] = stabilizer_histogram.get(str(record["stabilizer_size"]), 0) + 1
        fixed_option_histogram[str(len(record["fixed_endpoint_options_for_equivariant_selector"]))] = fixed_option_histogram.get(str(len(record["fixed_endpoint_options_for_equivariant_selector"])), 0) + 1
    facts = {
        "p2877_rechecked": p2877.get("status") == "P2877_PINNED_NON_CIRCULANT_DEFECT_ORIGIN_IMPORT_NO_GO_AUDIT_NO_CLOSURE",
        "all_pinned_defect_vectors_generated": len(vectors) == 12 * 3**3,
        "d12_group_size_checked": len(D12) == 24,
        "no_orbit_forces_endpoint_11": len(target_forcing_records) == 0,
        "no_orbit_has_unique_endpoint_selector": len(unique_endpoint_records) == 0,
        "singleton_orbits_are_not_endpoint_11_source_laws": all(record["fixed_endpoint_options_for_equivariant_selector"] != [TARGET_ENDPOINT] for record in singleton_support_orbits),
    }
    return {
        "status": "P2878_D12_EQUIVARIANT_PINNED_DEFECT_ORIGIN_LAW_NO_GO_AUDIT_NO_CLOSURE",
        "input_hashes": {"P2877": sha(P2877)},
        "d12_equivariant_pinned_defect_origin_law_no_go_audit": {
            "input_status_rechecked": p2877.get("status"),
            "candidate_class": "D12-equivariant endpoint selectors sourced by ternary radius-1 pinned defect density orbits",
            "modulus": MODULUS,
            "target_endpoint": TARGET_ENDPOINT,
            "raw_pinned_defect_vector_count": len(vectors),
            "unique_density_vector_count": len(unique_vectors),
            "d12_group_size": len(D12),
            "orbit_count": len(records),
            "stabilizer_size_histogram": stabilizer_histogram,
            "fixed_endpoint_option_count_histogram": fixed_option_histogram,
            "target_11_forcing_orbit_count": len(target_forcing_records),
            "unique_endpoint_selector_orbit_count": len(unique_endpoint_records),
            "singleton_support_orbit_count": len(singleton_support_orbits),
            "sample_orbit_records": records[:12],
            "proof_certificate": {
                "enumeration_step": "All 12*3^3 pinned radius-1 ternary defect densities are generated as endpoint 12-vectors and quotiented by the 24-element D12 action.",
                "equivariant_map_criterion": "For an orbit representative x, a D12-equivariant endpoint selector can choose only endpoints fixed by Stab(x).",
                "finite_result": "No orbit has fixed endpoint set {11}; in fact no orbit has a unique fixed endpoint at all.",
                "sourcehood_step": "Any endpoint-11 choice therefore imports an orbit representative, pin, orientation, or endpoint label; no strict non-premise origin/support law or unit-bearing 9/5 theorem is exported.",
            },
        },
        "acceptance_matrix": {
            "facts": facts,
            "accepted_as_independent_strict_origin_support_law_for_endpoint_11": False,
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
            "reason": "P2878 quotients the P2877 pinned-defect density family by the full D12 action and applies the stabilizer fixed-point criterion for equivariant endpoint selectors.  No orbit forces endpoint 11, and no orbit has any unique fixed endpoint selector.  Thus an endpoint-11 origin/support law is still imported rather than strict-sourced, and no unit-bearing 9/5 theorem is exported.",
            "next_honest_step": "Do not replay D12-equivariant selectors on pinned-defect orbits, imported orbit representatives, origin pins, center-only delta defects, circulant/Fourier constructions, D12 irreps/characters, or endpoint predicates as sourcehood.  A next proof-grade move must introduce a genuinely new non-D12-symmetric strict origin/support source with a computed coupling to the 9/5 coefficient, or pivot outside the endpoint-localizer/source family; otherwise preserve no-new-live-frontier.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    audit = payload["d12_equivariant_pinned_defect_origin_law_no_go_audit"]
    lines = [
        "# P2878/S1828 D12-equivariant pinned-defect origin-law no-go audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Equivariant selector audit",
        f"- candidate class: `{audit['candidate_class']}`",
        f"- raw pinned defect vector count: `{audit['raw_pinned_defect_vector_count']}`",
        f"- unique density vector count: `{audit['unique_density_vector_count']}`",
        f"- D12 orbit count: `{audit['orbit_count']}`",
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
    payload = build_payload(read_json(P2877))
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(
        STRICT_EQUATION_SHEET,
        "P2878/S1828 D12-equivariant pinned-defect origin-law no-go audit",
        "## P2878/S1828 D12-equivariant pinned-defect origin-law no-go audit\n\n"
        "`P2878/S1828` quotients the `P2877` pinned radius-1 defect density family by the full `D12` action and applies the stabilizer fixed-point criterion for equivariant endpoint selectors.  Across the finite orbit table, no orbit forces endpoint `11` and no orbit has a unique fixed endpoint selector.  Endpoint `11` therefore remains an imported orbit-representative/pin/orientation/label choice, not an independently exported strict origin/support law.  No unit-bearing `9/5` coefficient theorem, boundary/orientation source law, selector/localizer source, strict damping bridge, `L_total`, EOM, Hamiltonian, role transfer, or ToE closure follows.\n",
    )
    append_once(
        STRICT_LAGRANGIAN_DRAFT,
        "P2878/S1828 D12-equivariant pinned-defect origin-law `L_total` guard",
        "## P2878/S1828 D12-equivariant pinned-defect origin-law `L_total` guard\n\n"
        "`P2878/S1828` adds no strict action term.  D12-equivariant selectors on pinned-defect density orbits do not export a localized unit-bearing boundary/source density at `11`, non-premise origin/support law, coupling coefficient, localization/pullback theorem, variational chain rule, nonproxy `L_total`, EOM, or Hamiltonian source.\n",
    )
    append_once(
        AGENTS,
        "Current D12-equivariant pinned-defect origin-law no-go guardrail (P2878/S1828, 2026-06-18)",
        "## Current D12-equivariant pinned-defect origin-law no-go guardrail (P2878/S1828, 2026-06-18)\n\n"
        "- P2878 quotients the P2877 pinned radius-1 defect density family by the full `D12` action and applies the stabilizer fixed-point criterion for endpoint selectors.\n"
        "- No orbit forces endpoint `11`, and no orbit has a unique fixed endpoint selector; endpoint `11` remains an imported representative/pin/orientation/label choice.\n"
        "- Do not promote `D12`-equivariant pinned-defect selectors, imported orbit representatives, origin pins, center-only delta defects, circulant/Fourier constructions, `D12` irreps/characters, endpoint predicates, endpoint-label imports, or prime-5 scaled coefficients to strict damping/compression bridge, selector closure, role transfer, unit-bearing `L_total`, EOM, Hamiltonian, or ToE closure.\n"
        "- A next admissible proof-grade move must introduce a genuinely new non-`D12`-symmetric strict origin/support source with computed coupling to the `9/5` coefficient, or pivot outside the endpoint-localizer/source family; otherwise preserve no-new-live-frontier.\n",
    )
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
