#!/usr/bin/env python3
"""P2899/S1849: post-defect potential/readiness matrix.

P2898 closed the nearest non-circulant one-edge defect repair unless a law sources
the defect placement itself.  P2899 is an honest state-map/potential audit rather
than another replay: it reads P2895-P2898 and separates positive structural
symptoms from remaining closure blockers.

The result is deliberately mixed.  There is real conditional potential: a free-12
carrier capacity exists, equivariant maps exist after offset choice, and a
non-circulant defect would break translation if its placement were strictly
sourced.  But every positive symptom is premise/placement/offset conditional;
no current artifact exports the missing source law, coupling theorem, action
density, L_total, EOM, Hamiltonian, role transfer, or ToE closure.
"""
from __future__ import annotations

import json
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2803_s1753_meringer_scd_import_decode_validation_certificate import sha
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN

P2895 = GEN / "p2895_s1845_free_12_torsor_basepoint_polarity_law_no_go.json"
P2896 = GEN / "p2896_s1846_invariant_scalar_basepoint_law_exhaustion.json"
P2897 = GEN / "p2897_s1847_circulant_relation_basepoint_obstruction.json"
P2898 = GEN / "p2898_s1848_single_defect_relation_import_boundary.json"
OUT = GEN / "p2899_s1849_post_defect_potential_readiness_matrix.json"
MD = GEN / "p2899_s1849_post_defect_potential_readiness_matrix.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

INPUTS = {"P2895": P2895, "P2896": P2896, "P2897": P2897, "P2898": P2898}
EXPECTED_STATUSES = {
    "P2895": "P2895_FREE_12_TORSOR_BASEPOINT_POLARITY_LAW_NO_GO_NO_CLOSURE",
    "P2896": "P2896_INVARIANT_SCALAR_BASEPOINT_LAW_EXHAUSTION_NO_CLOSURE",
    "P2897": "P2897_CIRCULANT_RELATION_BASEPOINT_OBSTRUCTION_NO_CLOSURE",
    "P2898": "P2898_SINGLE_DEFECT_RELATION_IMPORT_BOUNDARY_NO_CLOSURE",
}


def nested_get(payload: dict[str, Any], path: list[str], default: Any = None) -> Any:
    value: Any = payload
    for key in path:
        if not isinstance(value, dict) or key not in value:
            return default
        value = value[key]
    return value


def build_symptom_matrix(loaded: dict[str, dict[str, Any]]) -> list[dict[str, Any]]:
    p2895_audit = nested_get(loaded["P2895"], ["free_12_torsor_basepoint_polarity_law_no_go", "torsor_audit"], {})
    p2898_audit = nested_get(loaded["P2898"], ["single_defect_relation_import_boundary", "single_defect_audit"], {})
    return [
        {
            "symptom": "minimal_free_12_capacity_exists",
            "artifact": "P2895",
            "positive": p2895_audit.get("equivariant_map_count") == 12,
            "computed_witness": {"equivariant_map_count": p2895_audit.get("equivariant_map_count")},
            "boundary": "capacity only; offset/basepoint remains unsourced",
        },
        {
            "symptom": "non_circulant_defect_can_break_translation_if_placed",
            "artifact": "P2898",
            "positive": p2898_audit.get("single_defect_candidate_count") == 589824,
            "computed_witness": {"single_defect_candidate_count": p2898_audit.get("single_defect_candidate_count")},
            "boundary": "breaks translation only after importing labelled defect placement",
        },
        {
            "symptom": "finite_obstruction_chain_is_executable_and_reproducible",
            "artifact": "P2895-P2898",
            "positive": all(loaded[key].get("status") == EXPECTED_STATUSES[key] for key in EXPECTED_STATUSES),
            "computed_witness": {key: loaded[key].get("status") for key in EXPECTED_STATUSES},
            "boundary": "proof infrastructure exists, but infrastructure is not sourcehood",
        },
        {
            "symptom": "next_missing_object_is_sharp",
            "artifact": "P2898",
            "positive": p2898_audit.get("source_neutral_defect_placement_count") == 0,
            "computed_witness": {"source_neutral_defect_placement_count": p2898_audit.get("source_neutral_defect_placement_count")},
            "boundary": "the sharp object is a strict law sourcing defect placement and coupling it to 9/5 density; it is not present",
        },
    ]


def build_blocker_matrix(loaded: dict[str, dict[str, Any]]) -> list[dict[str, Any]]:
    p2895_audit = nested_get(loaded["P2895"], ["free_12_torsor_basepoint_polarity_law_no_go", "torsor_audit"], {})
    p2896_audit = nested_get(loaded["P2896"], ["invariant_scalar_basepoint_law_exhaustion", "scalar_exhaustion"], {})
    p2897_audit = nested_get(loaded["P2897"], ["circulant_relation_basepoint_obstruction", "relation_exhaustion"], {})
    p2898_audit = nested_get(loaded["P2898"], ["single_defect_relation_import_boundary", "single_defect_audit"], {})
    return [
        {
            "blocker": "unpointed_free_torsor_has_no_invariant_basepoint",
            "artifact": "P2895",
            "blocked": p2895_audit.get("invariant_basepoint_count") == 0,
            "computed_value": p2895_audit.get("invariant_basepoint_count"),
        },
        {
            "blocker": "invariant_scalar_scores_select_no_unique_level_or_extremum",
            "artifact": "P2896",
            "blocked": p2896_audit.get("total_unique_marked_level_law_count") == 0 and p2896_audit.get("total_unique_argmin_law_count") == 0 and p2896_audit.get("total_unique_argmax_law_count") == 0,
            "computed_value": {
                "unique_marked_levels": p2896_audit.get("total_unique_marked_level_law_count"),
                "unique_argmin": p2896_audit.get("total_unique_argmin_law_count"),
                "unique_argmax": p2896_audit.get("total_unique_argmax_law_count"),
            },
        },
        {
            "blocker": "circulant_relations_select_no_unique_vertex",
            "artifact": "P2897",
            "blocked": p2897_audit.get("relations_selecting_unique_vertex_count") == 0,
            "computed_value": p2897_audit.get("relations_selecting_unique_vertex_count"),
        },
        {
            "blocker": "single_defect_placement_is_imported_not_sourced",
            "artifact": "P2898",
            "blocked": p2898_audit.get("source_neutral_defect_placement_count") == 0,
            "computed_value": p2898_audit.get("source_neutral_defect_placement_count"),
        },
        {
            "blocker": "no_coupling_theorem_to_9_over_5_variational_density",
            "artifact": "P2895-P2898",
            "blocked": all(not nested_get(loaded[key], ["acceptance_matrix", "exports_unique_coupling_to_9_over_5_carrier"], False) for key in loaded),
            "computed_value": {key: nested_get(loaded[key], ["acceptance_matrix", "exports_unique_coupling_to_9_over_5_carrier"], None) for key in loaded},
        },
        {
            "blocker": "no_ltotal_eom_hamiltonian_or_toe_closure",
            "artifact": "P2895-P2898",
            "blocked": all(not any(nested_get(loaded[key], ["decision", "negative_export_flags", flag], False) for flag in ("ltotal_exported", "eom_closure_exported", "hamiltonian_closure_exported", "toe_closure_exported")) for key in loaded),
            "computed_value": {key: {flag: nested_get(loaded[key], ["decision", "negative_export_flags", flag], None) for flag in ("ltotal_exported", "eom_closure_exported", "hamiltonian_closure_exported", "toe_closure_exported")} for key in loaded},
        },
    ]


def build_payload() -> dict[str, Any]:
    loaded = {key: read_json(path) for key, path in INPUTS.items()}
    symptoms = build_symptom_matrix(loaded)
    blockers = build_blocker_matrix(loaded)
    positive_count = sum(1 for row in symptoms if row["positive"])
    blocked_count = sum(1 for row in blockers if row["blocked"])
    closure_ready = positive_count == len(symptoms) and blocked_count == 0
    return {
        "status": "P2899_POST_DEFECT_POTENTIAL_READINESS_MATRIX_NO_CLOSURE",
        "input_hashes": {key: sha(path) for key, path in INPUTS.items()},
        "input_statuses": {key: loaded[key].get("status") for key in loaded},
        "potential_readiness_matrix": {
            "positive_symptom_count": positive_count,
            "blocker_count": blocked_count,
            "closure_ready": closure_ready,
            "toe_potential_class": "conditional_structural_potential_but_no_current_ToE_closure",
            "positive_symptoms": symptoms,
            "blockers": blockers,
        },
        "acceptance_matrix": {
            "all_inputs_rechecked": all(loaded[key].get("status") == EXPECTED_STATUSES[key] for key in EXPECTED_STATUSES),
            "positive_symptoms_exist": positive_count > 0,
            "all_positive_symptoms_are_conditional": True,
            "all_named_blockers_still_blocked": blocked_count == len(blockers),
            "accepted_as_toe_potential_evidence": True,
            "accepted_as_toe_closure": False,
        },
        "decision": {
            "negative_export_flags": {
                "strict_defect_placement_source_law_exported": False,
                "nonimported_basepoint_or_polarity_law_exported": False,
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
            "reason": "P2899 separates potential from closure after P2895-P2898.  Positive symptoms are real but conditional: free-12 capacity exists, equivariant maps exist after offset choice, one-edge defects would break translation after labelled placement, and the missing object is now sharply specified.  The blockers remain decisive: zero invariant basepoints, zero scalar unique selectors, zero circulant unique vertices, zero source-neutral defect placements, no coupling theorem to the 9/5 variational density, and no L_total/EOM/Hamiltonian/ToE export.",
            "next_honest_step": "Do not market the current chain as ToE closure.  It is legitimate positive potential evidence only in the weak/conditional sense that the finite obstruction chain has isolated the exact missing object: a strict law sourcing defect placement/basepoint/polarity and coupling it to the 9/5 variational density.  The next proof-grade move must either construct that law with computed nonconventional sign/phase and a coupling theorem, or pivot to a genuinely different typed object outside torsor/basepoint/scalar-score/relation/defect/support/orbit/Fourier/inventory constructions; otherwise preserve no-new-live-frontier.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    matrix = payload["potential_readiness_matrix"]
    lines = [
        "# P2899/S1849 post-defect potential/readiness matrix",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Readiness summary",
        f"- positive symptom count: `{matrix['positive_symptom_count']}`",
        f"- blocker count: `{matrix['blocker_count']}`",
        f"- closure ready: `{matrix['closure_ready']}`",
        f"- ToE potential class: `{matrix['toe_potential_class']}`",
        "",
        "## Positive symptoms",
    ]
    for row in matrix["positive_symptoms"]:
        lines.append(f"- `{row['symptom']}` from `{row['artifact']}`: positive=`{row['positive']}`; boundary: {row['boundary']}")
    lines.append("")
    lines.append("## Blockers")
    for row in matrix["blockers"]:
        lines.append(f"- `{row['blocker']}` from `{row['artifact']}`: blocked=`{row['blocked']}`; value=`{row['computed_value']}`")
    lines.extend(["", "## Boundary", payload["decision"]["reason"], "", "## Recommendation", payload["decision"]["next_honest_step"]])
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(
        STRICT_EQUATION_SHEET,
        "P2899/S1849 post-defect potential/readiness matrix",
        "## P2899/S1849 post-defect potential/readiness matrix\n\n"
        "`P2899/S1849` separates conditional ToE potential from closure after `P2895-P2898`.  It records `4` positive symptoms (free-`12` capacity, offset-conditioned equivariant maps, placement-conditioned one-edge translation breaking, and a sharply specified missing object) but also `6` still-active blockers: `0` invariant basepoints, `0` scalar unique selectors, `0` circulant unique vertices, `0` source-neutral defect placements, no coupling theorem to the `9/5` variational density, and no `L_total`/EOM/Hamiltonian/ToE export.  This is weak conditional potential evidence, not action-density or ToE closure.\n",
    )
    append_once(
        STRICT_LAGRANGIAN_DRAFT,
        "P2899/S1849 post-defect potential/readiness `L_total` guard",
        "## P2899/S1849 post-defect potential/readiness `L_total` guard\n\n"
        "`P2899/S1849` is a readiness/state-map matrix.  It identifies conditional structural potential but adds no strict defect-placement source law, no localized unit-bearing density, no coupling theorem to the `9/5` carrier, and no variational chain rule into nonproxy `L_total`, EOM, or Hamiltonian source.\n",
    )
    append_once(
        AGENTS,
        "Current post-defect potential/readiness guardrail (P2899/S1849, 2026-06-19)",
        "## Current post-defect potential/readiness guardrail (P2899/S1849, 2026-06-19)\n\n"
        "- P2899 separates conditional potential from closure after P2895-P2898: free-`12` capacity and one-edge defect translation-breaking are real positive symptoms only after offset/placement/source premises.\n"
        "- The active blockers remain: `0` invariant basepoints, `0` scalar unique selectors, `0` circulant unique vertices, `0` source-neutral defect placements, no coupling theorem to the `9/5` variational density, and no `L_total`/EOM/Hamiltonian/ToE export.\n"
        "- Treat the theory's ToE potential here as conditional structural potential, not closure evidence.  Do not promote capacity, equivariant maps, one-edge defects, labelled placements, relation profiles, scalar scores, support/orbit/Fourier data, or inventory hits to strict phase/origin sourcehood, strict damping/compression bridge, role transfer, `L_total`, EOM, Hamiltonian, or ToE closure.\n"
        "- A next admissible proof-grade move must construct a strict law sourcing defect placement/basepoint/polarity with computed nonconventional sign/phase and coupling to the `9/5` variational density, pivot to a genuinely different typed object outside torsor/basepoint/scalar-score/relation/defect/support/orbit/Fourier/inventory constructions, or preserve no-new-live-frontier.\n",
    )
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
