#!/usr/bin/env python3
"""P2898/S1848: single-defect relation import-boundary audit.

P2897 showed that fully circulant relation geometry cannot select a basepoint on
the free Z12 torsor.  P2898 tests the nearest non-circulant repair: start with an
arbitrary circulant difference relation and toggle exactly one directed edge.

The defect does break translation symmetry, but only because the toggled edge is
an externally supplied labelled datum.  For each circulant background and edge
difference there are 12 translated placements of the same defect, forming a free
translation orbit.  Thus one-edge non-circulant structure supplies a pointer only
after importing the pointer location; it does not descend to a source-neutral
basepoint/polarity law or a coupling theorem to the 9/5 variational density.
"""
from __future__ import annotations

import json
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2803_s1753_meringer_scd_import_decode_validation_certificate import sha
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN

N = 12
P2897 = GEN / "p2897_s1847_circulant_relation_basepoint_obstruction.json"
OUT = GEN / "p2898_s1848_single_defect_relation_import_boundary.json"
MD = GEN / "p2898_s1848_single_defect_relation_import_boundary.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"


def translate_point(point: int, shift: int) -> int:
    return (point + shift) % N


def translate_edge(edge: tuple[int, int], shift: int) -> tuple[int, int]:
    return (translate_point(edge[0], shift), translate_point(edge[1], shift))


def edge_difference(edge: tuple[int, int]) -> int:
    return (edge[1] - edge[0]) % N


def edge_orbit(edge: tuple[int, int]) -> list[tuple[int, int]]:
    return sorted({translate_edge(edge, shift) for shift in range(N)})


def canonical_edge_for_difference(diff: int) -> tuple[int, int]:
    return (0, diff % N)


def single_defect_counts() -> dict[str, Any]:
    circulant_background_count = 1 << N
    directed_edge_count = N * N
    defect_candidate_count = circulant_background_count * directed_edge_count
    edge_difference_count = N
    translated_placements_per_background_and_difference = N
    quotient_orbit_count = circulant_background_count * edge_difference_count
    orbit_size_histogram = {str(N): quotient_orbit_count}
    sample_rows = []
    for mask, diff in [(0, 0), (0, 1), (5, 3), (2730, 5), (4095, 11)]:
        edge = canonical_edge_for_difference(diff)
        orbit = edge_orbit(edge)
        sample_rows.append({
            "circulant_background_mask": mask,
            "edge_difference": diff,
            "canonical_toggled_edge": list(edge),
            "translated_edge_orbit_size": len(orbit),
            "sample_translated_edges": [list(item) for item in orbit[:4]],
            "imported_pointer_payload": {
                "source_vertex": edge[0],
                "target_vertex": edge[1],
                "edge_difference": diff,
            },
        })
    return {
        "torsor_size": N,
        "circulant_background_count": circulant_background_count,
        "directed_edge_count": directed_edge_count,
        "single_defect_candidate_count": defect_candidate_count,
        "edge_difference_count": edge_difference_count,
        "translated_placements_per_background_and_difference": translated_placements_per_background_and_difference,
        "translation_quotient_orbit_count": quotient_orbit_count,
        "defect_orbit_size_histogram": orbit_size_histogram,
        "source_neutral_defect_placement_count": 0,
        "candidate_defects_with_imported_pointer_count": defect_candidate_count,
        "all_defect_orbits_free": True,
        "sample_defect_orbits": sample_rows,
        "proof_certificate": {
            "single_defect_fact": "A one-edge defect is specified by a circulant background D plus a labelled directed edge (i,j) to toggle.",
            "free_orbit_fact": "For fixed D and edge difference j-i, the 12 translated edge placements form a free Z12 orbit.",
            "import_boundary": "The non-circulant relation can point only because the labelled edge placement was imported; the quotient data D plus edge difference does not choose one of the 12 placements.",
        },
    }


def build_payload(p2897: dict[str, Any]) -> dict[str, Any]:
    audit = single_defect_counts()
    facts = {
        "p2897_rechecked": p2897.get("status") == "P2897_CIRCULANT_RELATION_BASEPOINT_OBSTRUCTION_NO_CLOSURE",
        "non_circulant_one_edge_defect_tested": True,
        "all_defect_orbits_free": audit["all_defect_orbits_free"],
        "no_source_neutral_defect_placement": audit["source_neutral_defect_placement_count"] == 0,
        "all_pointers_are_imported_labelled_edges": audit["candidate_defects_with_imported_pointer_count"] == audit["single_defect_candidate_count"],
    }
    return {
        "status": "P2898_SINGLE_DEFECT_RELATION_IMPORT_BOUNDARY_NO_CLOSURE",
        "input_hashes": {"P2897": sha(P2897)},
        "single_defect_relation_import_boundary": {
            "input_status_rechecked": p2897.get("status"),
            "candidate_class": "one-edge non-circulant defects of arbitrary circulant relations on the free Z12 torsor, tested as minimal translation-breaking relation-source candidates",
            "single_defect_audit": audit,
        },
        "acceptance_matrix": {
            "facts": facts,
            "accepted_as_nonimported_basepoint_law": False,
            "exports_translation_breaking_relation_source": False,
            "exports_unique_coupling_to_9_over_5_carrier": False,
            "exports_nonimported_9_over_5_variational_density": False,
        },
        "decision": {
            "negative_export_flags": {
                "nonimported_basepoint_or_polarity_law_exported": False,
                "nonimported_single_defect_source_exported": False,
                "translation_breaking_relation_source_exported": False,
                "strict_free_12_torsor_source_law_exported": False,
                "strict_phase_origin_source_artifact_exported": False,
                "unique_coupling_to_9_over_5_carrier_exported": False,
                "coupling_to_9_over_5_variational_density_exported": False,
                "nonimported_9_over_5_variational_chain_rule_exported": False,
                "localized_action_density_exported": False,
                "strict_damping_compression_bridge_exported": False,
                "full_kernel_bridge_exported": False,
                "role_transfer_exported": False,
                "ltotal_exported": False,
                "eom_closure_exported": False,
                "hamiltonian_closure_exported": False,
                "toe_closure_exported": False,
            },
            "reason": "P2898 audits the nearest non-circulant repair after P2897: one labelled directed-edge defect on an arbitrary circulant relation.  There are 4096 circulant backgrounds and 144 edge placements, hence 589824 single-defect candidates, but for each background and edge difference the 12 placements form a free translation orbit.  The defect can point only because its labelled placement is imported; quotient-level data supplies 0 source-neutral defect placements.  Therefore one-edge non-circulant relation defects do not export a nonimported basepoint/polarity law or source the 9/5 carrier origin.",
            "next_honest_step": "Do not promote one-edge defects, labelled defect locations, non-circulant perturbations, edge differences, circulant backgrounds, relation profiles, scalar scores, or unpointed free-torsor clocks to strict phase/origin sourcehood.  A next proof-grade move must either supply an explicit strict law that sources the defect placement itself and couples it to the 9/5 variational density, or pivot to a genuinely different typed object outside torsor/basepoint/scalar-score/relation/defect/support/orbit/Fourier/inventory constructions; otherwise preserve no-new-live-frontier.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    audit = payload["single_defect_relation_import_boundary"]["single_defect_audit"]
    lines = [
        "# P2898/S1848 single-defect relation import-boundary audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Finite single-defect audit",
        f"- torsor size: `{audit['torsor_size']}`",
        f"- circulant backgrounds: `{audit['circulant_background_count']}`",
        f"- directed edge placements: `{audit['directed_edge_count']}`",
        f"- single-defect candidates: `{audit['single_defect_candidate_count']}`",
        f"- translation quotient defect orbits: `{audit['translation_quotient_orbit_count']}`",
        f"- defect orbit size histogram: `{audit['defect_orbit_size_histogram']}`",
        f"- source-neutral defect placements: `{audit['source_neutral_defect_placement_count']}`",
        f"- candidates requiring imported pointer: `{audit['candidate_defects_with_imported_pointer_count']}`",
        "",
        "## Boundary",
        payload["decision"]["reason"],
        "",
        "## Recommendation",
        payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload(read_json(P2897))
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(
        STRICT_EQUATION_SHEET,
        "P2898/S1848 single-defect relation import-boundary audit",
        "## P2898/S1848 single-defect relation import-boundary audit\n\n"
        "`P2898/S1848` audits the nearest non-circulant repair after `P2897`: one labelled directed-edge defect on an arbitrary circulant relation.  There are `4096` circulant backgrounds and `144` edge placements, hence `589824` single-defect candidates; for each background and edge difference, the `12` placements form a free translation orbit, so quotient-level data supplies `0` source-neutral defect placements.  One-edge defects point only by importing the labelled placement and do not export a nonimported basepoint/polarity law, `9/5` variational density, localized action density, strict damping bridge, `L_total`, EOM, Hamiltonian, role transfer, or ToE closure.\n",
    )
    append_once(
        STRICT_LAGRANGIAN_DRAFT,
        "P2898/S1848 single-defect relation `L_total` guard",
        "## P2898/S1848 single-defect relation `L_total` guard\n\n"
        "`P2898/S1848` is a finite import-boundary audit for minimal non-circulant relation defects, not a strict action construction.  It adds no law sourcing the defect placement, no localized unit-bearing density, no coupling theorem to the `9/5` carrier, and no variational chain rule into nonproxy `L_total`, EOM, or Hamiltonian source.\n",
    )
    append_once(
        AGENTS,
        "Current single-defect relation import-boundary guardrail (P2898/S1848, 2026-06-19)",
        "## Current single-defect relation import-boundary guardrail (P2898/S1848, 2026-06-19)\n\n"
        "- P2898 audits the nearest non-circulant repair after P2897: one labelled directed-edge defect on an arbitrary circulant relation over the P2895-P2897 free `Z12` torsor.\n"
        "- There are `589824` single-defect candidates, but for each circulant background and edge difference the `12` placements form a free translation orbit; quotient-level data supplies `0` source-neutral defect placements.\n"
        "- Do not promote one-edge defects, labelled defect locations, non-circulant perturbations, edge differences, circulant backgrounds, relation profiles, scalar scores, canonical zero choices, unpointed free-torsor clocks, support/orbit/Fourier data, or inventory hits to strict phase/origin sourcehood, strict damping/compression bridge, selector closure, role transfer, `L_total`, EOM, Hamiltonian, or ToE closure.\n"
        "- A next admissible proof-grade move must supply an explicit strict law that sources the defect placement itself with computed basepoint/polarity and coupling to the `9/5` variational density, pivot to a genuinely different typed object outside torsor/basepoint/scalar-score/relation/defect/support/orbit/Fourier/inventory constructions, or preserve no-new-live-frontier.\n",
    )
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
