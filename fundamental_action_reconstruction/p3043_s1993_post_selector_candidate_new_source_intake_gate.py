#!/usr/bin/env python3
"""P3043/S1993: post-selector-candidate new-source intake gate.

P3042 closed the P3038-P3041 integrated selector-candidate route as a
no-selector-export certificate.  P3043 follows the guardrail recommendation:
perform a broad state-map/new-object intake rather than replaying the exhausted
receiver classes.  Because no external new strict source law is supplied in the
current artifacts, the script constructs the missing theoretical object as an
acceptance gate: an explicit complement signature and source-law predicate matrix
for any future selector mechanism that may not match familiar human schemas.

This is not closure.  It preserves the computational hint that branch separation
is possible, but requires a genuinely new strict source law outside sine/chiral
phase receivers, rho/path tuning, unit normalization/imports, and nonproxy
placeholders before any further selector promotion is admissible.
"""
from __future__ import annotations

import hashlib, json
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p3037_s1987_selector_mechanism_hint_sheaf_acceptance_matrix import OUT as P3037
from p3042_s1992_integrated_selector_candidate_reconciliation_certificate import OUT as P3042

OUT = GEN / "p3043_s1993_post_selector_candidate_new_source_intake_gate.json"
MD = GEN / "p3043_s1993_post_selector_candidate_new_source_intake_gate.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

EXHAUSTED_CLASSES = [
    "sine_or_chiral_phase_receiver",
    "rho_or_retarded_path_tuning_receiver",
    "unit_normalization_or_import_receiver",
    "nonproxy_variational_placeholder",
    "feature_hint_coverage_without_source_row",
]
REQUIRED_NEW_SOURCE_PREDICATES = [
    "strict_nadsoliton_provenance",
    "nonpremise_source_law",
    "outside_exhausted_receiver_classes",
    "computable_nonzero_signed_or_branch_value",
    "chart_or_gauge_independent_localizer",
    "explicit_coupling_to_selector_torsor_or_readout",
    "unit_or_nonproxy_installation_when_physical_export_is_claimed",
]


def load_inputs() -> dict[str, Any]:
    return {"P3037": read_json(P3037), "P3042": read_json(P3042)}


def candidate_rows() -> list[dict[str, Any]]:
    return [
        {
            "candidate_class": "P3037_hint_sheaf_feature_coverage",
            "novel_outside_exhausted_classes": False,
            "computable_positive": True,
            "accepted_new_source_law": False,
            "failure": "feature coverage is a hint object, not a source row exporting a selector mechanism",
        },
        {
            "candidate_class": "P3038_integrated_branch_separating_receiver",
            "novel_outside_exhausted_classes": False,
            "computable_positive": True,
            "accepted_new_source_law": False,
            "failure": "branch separation is real but P3042 classifies the source premises as open",
        },
        {
            "candidate_class": "P3039_P3041_exhausted_receiver_family",
            "novel_outside_exhausted_classes": False,
            "computable_positive": True,
            "accepted_new_source_law": False,
            "failure": "sine/chiral phase, rho/path tuning, unit normalization/import, and placeholder classes are repetition-gated",
        },
        {
            "candidate_class": "unsupplied_genuinely_new_strict_source_law_slot",
            "novel_outside_exhausted_classes": True,
            "computable_positive": False,
            "accepted_new_source_law": False,
            "failure": "the slot states the required kind of object, but no concrete formula/artifact is supplied in current repo contents",
        },
    ]


def predicate_matrix(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    matrix = []
    for row in rows:
        supplied = row["candidate_class"] != "unsupplied_genuinely_new_strict_source_law_slot"
        outside = row["novel_outside_exhausted_classes"]
        accepted = row["accepted_new_source_law"]
        matrix.append({
            "candidate_class": row["candidate_class"],
            "predicates": {
                "strict_nadsoliton_provenance": supplied and accepted,
                "nonpremise_source_law": accepted,
                "outside_exhausted_receiver_classes": outside,
                "computable_nonzero_signed_or_branch_value": row["computable_positive"],
                "chart_or_gauge_independent_localizer": accepted,
                "explicit_coupling_to_selector_torsor_or_readout": accepted,
                "unit_or_nonproxy_installation_when_physical_export_is_claimed": accepted,
            },
            "accepted_new_source_law": accepted,
            "failure": row["failure"],
        })
    return matrix


def build_matrix() -> dict[str, Any]:
    inputs = load_inputs()
    rows = candidate_rows()
    predicates = predicate_matrix(rows)
    state_lanes = [
        {"lane": "selector_integrated_candidate", "status": "bounded_no_export", "basis": "P3042"},
        {"lane": "unknown_selector_hint_sheaf", "status": "hints_real_but_no_source_row", "basis": "P3037"},
        {"lane": "new_strict_source_law_intake", "status": "no_concrete_new_formula_supplied", "basis": "P3043"},
    ]
    obligations = [
        {"obligation": "broad_state_map_consulted", "satisfied": True, "detail": "selector integrated candidate, hint sheaf, and new-source intake lanes are reconciled"},
        {"obligation": "exhausted_receiver_classes_declared", "satisfied": True, "detail": ", ".join(EXHAUSTED_CLASSES)},
        {"obligation": "new_source_predicate_matrix_constructed", "satisfied": True, "detail": f"{len(REQUIRED_NEW_SOURCE_PREDICATES)} predicates are explicit"},
        {"obligation": "concrete_new_strict_source_law_supplied", "satisfied": any(row["accepted_new_source_law"] for row in predicates), "detail": "no row exports a concrete new strict source law"},
        {"obligation": "selector_lane_reopened", "satisfied": False, "detail": "P3042 remains the active no-selector-export certificate"},
    ]
    return {
        "object": "PostSelectorCandidate_NewSourceIntakeGate",
        "exhausted_receiver_classes": EXHAUSTED_CLASSES,
        "required_new_source_predicates": REQUIRED_NEW_SOURCE_PREDICATES,
        "state_lanes": state_lanes,
        "candidate_rows": rows,
        "predicate_matrix": predicates,
        "proof_obligations": obligations,
        "finite_certificate": {
            "state_lanes": len(state_lanes),
            "exhausted_receiver_classes": len(EXHAUSTED_CLASSES),
            "candidate_rows": len(rows),
            "predicate_count": len(REQUIRED_NEW_SOURCE_PREDICATES),
            "accepted_new_source_law_rows": sum(1 for row in predicates if row["accepted_new_source_law"]),
            "replay_gated_rows": sum(1 for row in rows if not row["novel_outside_exhausted_classes"]),
            "unsupplied_new_source_slots": sum(1 for row in rows if row["candidate_class"] == "unsupplied_genuinely_new_strict_source_law_slot"),
            "proof_obligations": len(obligations),
            "satisfied_proof_obligations": sum(1 for row in obligations if row["satisfied"]),
            "new_live_frontier_unlocked": any(row["accepted_new_source_law"] for row in predicates),
        },
        "input_statuses": {"P3037": inputs["P3037"]["status"], "P3042": inputs["P3042"]["status"]},
    }


def build_payload() -> dict[str, Any]:
    matrix = build_matrix()
    return {
        "status": "P3043_POST_SELECTOR_CANDIDATE_NEW_SOURCE_INTAKE_NO_NEW_LIVE_FRONTIER",
        "input_hashes": {
            "P3037": hashlib.sha256(P3037.read_bytes()).hexdigest() if P3037.exists() else None,
            "P3042": hashlib.sha256(P3042.read_bytes()).hexdigest() if P3042.exists() else None,
        },
        "constructed_theoretical_objects": matrix,
        "finite_certificate": matrix["finite_certificate"],
        "decision": {
            "bounded_no_go": "P3043 builds the post-P3042 intake gate rather than replaying exhausted selector receivers.  Current hints and branch-separating receivers remain real, but no concrete new strict source law is supplied outside the exhausted classes.  Therefore no new live selector frontier is unlocked on current artifacts.",
            "negative_export_flags": {k: False for k in ["new_live_frontier_unlocked", "strict_selector_mechanism_exported", "qw2191_discharged", "strict_selector_closure_exported", "observed_physics_exported", "unit_bearing_action_eom_hamiltonian_exported", "ltotal_exported", "bridge_closure_exported", "role_transfer_exported", "toe_closure_exported"]},
            "next_honest_step": "The next admissible move must provide one concrete formula/artifact satisfying the P3043 new-source predicates, or pivot to a different broad state-map lane with a genuinely new typed object.  Without that, preserve the P3042-P3043 no-new-live-frontier boundary rather than manufacturing selector closure.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    c = payload["finite_certificate"]
    lines = [
        "# P3043/S1993 post-selector-candidate new-source intake gate", "",
        f"Status: `{payload['status']}`", "", "## Finite certificate",
        f"- state lanes: `{c['state_lanes']}`",
        f"- exhausted receiver classes: `{c['exhausted_receiver_classes']}`",
        f"- candidate rows: `{c['candidate_rows']}`",
        f"- predicate count: `{c['predicate_count']}`",
        f"- accepted new source-law rows: `{c['accepted_new_source_law_rows']}`",
        f"- replay-gated rows: `{c['replay_gated_rows']}`",
        f"- unsupplied new-source slots: `{c['unsupplied_new_source_slots']}`",
        f"- satisfied proof obligations: `{c['satisfied_proof_obligations']}/{c['proof_obligations']}`",
        f"- new live frontier unlocked: `{c['new_live_frontier_unlocked']}`", "",
        "## Decision", payload["decision"]["bounded_no_go"], "", "## Recommendation", payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3043/S1993 post-selector-candidate new-source intake gate", "## P3043/S1993 post-selector-candidate new-source intake gate\n\n`P3043/S1993` follows the P3042 recommendation by building a broad state-map/new-source intake gate instead of replaying P3038 receiver classes.  It declares the exhausted receiver complement (sine/chiral phase receivers, rho/path tuning, unit normalization/imports, nonproxy placeholders, and hint coverage without source rows) and a seven-predicate acceptance matrix for any future source law.  Current hints and branch-separating receivers remain real, but no concrete new strict source law is supplied outside the exhausted classes, so no new live selector frontier, `QW-2191` discharge, selector closure, observed physics, unit-bearing action/EOM, `L_total`, bridge/role transfer, or ToE closure follows.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3043/S1993 new-source intake `L_total` guard", "## P3043/S1993 new-source intake `L_total` guard\n\n`P3043/S1993` adds no physical `L_total` term.  It constructs an intake gate for future strict source laws, but current rows are either replay-gated receiver classes or an unsupplied new-source slot with no unit-bearing variational/action/EOM installation.\n")
    append_once(AGENTS, "Current post-selector-candidate new-source intake guardrail (P3043/S1993, 2026-06-23)", "## Current post-selector-candidate new-source intake guardrail (P3043/S1993, 2026-06-23)\n\n- P3043 performs the broad state-map/new-object intake requested after P3042 instead of replaying the integrated selector candidate.\n- Exhausted receiver classes are explicit: sine/chiral phase receivers, rho/path tuning, unit normalization/imports, nonproxy placeholders, and hint coverage without a source row; current rows export zero accepted new strict source laws.\n- Do not promote the P3043 intake gate, unsupplied source slots, or any exhausted receiver class to `QW-2191` discharge, selector closure, observed physics, unit-bearing action/EOM, `L_total`, bridge/role-transfer, or ToE closure.\n- A next admissible move must supply one concrete formula/artifact satisfying the new-source predicates, or pivot to a different broad state-map lane with a genuinely new typed object.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
