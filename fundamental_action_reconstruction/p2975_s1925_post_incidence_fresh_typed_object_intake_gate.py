#!/usr/bin/env python3
"""P2975/S1925: post-incidence fresh typed-object intake gate.

P2974 closed the P2971 incidence theorem triad on current artifacts.  The next
honest move is not another incidence replay, but a finite gate for any genuinely
new strict typed object/theorem outside the incidence lane.  This module builds
that gate and audits the currently visible follow-up routes by content rather
than by number: incidence derivative replay, incidence coupling/localizer replay,
ratio-package scalar replay, Gamma/Jacobian import replay, selector replay,
bridge replay, and an unsupplied fresh-object placeholder.

The result is intentionally conservative: a formal intake object is constructed,
but no current candidate supplies all acceptance atoms needed to reopen a live
frontier or promote L_total/bridge/role-transfer/ToE closure.
"""
from __future__ import annotations

import hashlib, json
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p2974_s1924_incidence_nonproxy_variational_chain_rule_obstruction import OUT as P2974

OUT = GEN / "p2975_s1925_post_incidence_fresh_typed_object_intake_gate.json"
MD = GEN / "p2975_s1925_post_incidence_fresh_typed_object_intake_gate.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

ACCEPTANCE_ATOMS = [
    "genuinely_new_outside_incidence_lane",
    "not_replay_of_closed_ratio_or_gamma_or_selector_or_bridge_lane",
    "strict_nadsoliton_provenance_or_explicit_premise_scope",
    "targets_exactly_one_named_missing_object",
    "finite_witness_or_acceptance_test",
    "nonpromotion_boundaries_preserved_until_pass",
]


def candidate_rows() -> list[dict[str, Any]]:
    rows = [
        {
            "candidate": "incidence_formal_derivative_table_replay",
            "genuinely_new_outside_incidence_lane": False,
            "not_replay_of_closed_ratio_or_gamma_or_selector_or_bridge_lane": True,
            "strict_nadsoliton_provenance_or_explicit_premise_scope": False,
            "targets_exactly_one_named_missing_object": False,
            "finite_witness_or_acceptance_test": True,
            "nonpromotion_boundaries_preserved_until_pass": True,
            "failure": "P2974 already reduced this to formal derivative readiness, not nonproxy variational closure",
        },
        {
            "candidate": "incidence_source_or_coupling_replay",
            "genuinely_new_outside_incidence_lane": False,
            "not_replay_of_closed_ratio_or_gamma_or_selector_or_bridge_lane": True,
            "strict_nadsoliton_provenance_or_explicit_premise_scope": False,
            "targets_exactly_one_named_missing_object": False,
            "finite_witness_or_acceptance_test": True,
            "nonpromotion_boundaries_preserved_until_pass": True,
            "failure": "P2972/P2973 already bounded source-localizer and unit-bearing coupling attempts",
        },
        {
            "candidate": "ratio_package_scalar_9_5_or_unit_replay",
            "genuinely_new_outside_incidence_lane": False,
            "not_replay_of_closed_ratio_or_gamma_or_selector_or_bridge_lane": False,
            "strict_nadsoliton_provenance_or_explicit_premise_scope": False,
            "targets_exactly_one_named_missing_object": False,
            "finite_witness_or_acceptance_test": True,
            "nonpromotion_boundaries_preserved_until_pass": True,
            "failure": "replays primitive-mean, unit, k-selection, or scalar Euler lanes closed by P2965-P2970",
        },
        {
            "candidate": "Gamma_or_P2912_Jacobian_import_replay",
            "genuinely_new_outside_incidence_lane": True,
            "not_replay_of_closed_ratio_or_gamma_or_selector_or_bridge_lane": False,
            "strict_nadsoliton_provenance_or_explicit_premise_scope": False,
            "targets_exactly_one_named_missing_object": False,
            "finite_witness_or_acceptance_test": True,
            "nonpromotion_boundaries_preserved_until_pass": True,
            "failure": "imports an existing readiness skeleton without its missing Gamma/source/continuum theorems",
        },
        {
            "candidate": "selector_or_QW2191_replay_without_new_source",
            "genuinely_new_outside_incidence_lane": True,
            "not_replay_of_closed_ratio_or_gamma_or_selector_or_bridge_lane": False,
            "strict_nadsoliton_provenance_or_explicit_premise_scope": False,
            "targets_exactly_one_named_missing_object": True,
            "finite_witness_or_acceptance_test": True,
            "nonpromotion_boundaries_preserved_until_pass": True,
            "failure": "no new non-premise strict selector/orientation source is supplied",
        },
        {
            "candidate": "generic_legacy_to_strict_bridge_replay",
            "genuinely_new_outside_incidence_lane": True,
            "not_replay_of_closed_ratio_or_gamma_or_selector_or_bridge_lane": False,
            "strict_nadsoliton_provenance_or_explicit_premise_scope": False,
            "targets_exactly_one_named_missing_object": False,
            "finite_witness_or_acceptance_test": False,
            "nonpromotion_boundaries_preserved_until_pass": True,
            "failure": "generic bridge repetition lacks one new typed bridge/source atom or completion-map theorem",
        },
        {
            "candidate": "fresh_strict_typed_object_placeholder",
            "genuinely_new_outside_incidence_lane": False,
            "not_replay_of_closed_ratio_or_gamma_or_selector_or_bridge_lane": True,
            "strict_nadsoliton_provenance_or_explicit_premise_scope": False,
            "targets_exactly_one_named_missing_object": False,
            "finite_witness_or_acceptance_test": False,
            "nonpromotion_boundaries_preserved_until_pass": True,
            "failure": "placeholder names the required shape but supplies no object/formula/provider",
        },
    ]
    for row in rows:
        row["passed_atoms"] = [atom for atom in ACCEPTANCE_ATOMS if row[atom]]
        row["missing_atoms"] = [atom for atom in ACCEPTANCE_ATOMS if not row[atom]]
        row["accepted_current_new_live_frontier"] = all(row[atom] for atom in ACCEPTANCE_ATOMS)
    return rows


def lane_rows() -> list[dict[str, Any]]:
    return [
        {"lane": "incidence theorem triad", "latest_boundary": "P2972-P2974", "blocked_by": "source-localizer, unit-bearing coupling, and nonproxy variational chain rule all bounded no-go", "reopen_requires": "new non-incidence typed object/theorem", "live_frontier_unlocked": False},
        {"lane": "ratio-package scalar/unit/k lane", "latest_boundary": "P2965-P2970", "blocked_by": "unit/source/k installation and coefficient source law are closed on current artifacts", "reopen_requires": "new strict typed object outside scalar ratio arithmetic", "live_frontier_unlocked": False},
        {"lane": "Gamma/Jacobian variational readiness lane", "latest_boundary": "P2912/P2914 and later nonpromotion guards", "blocked_by": "missing strict field-variable provenance, Gamma source, and continuum measure theorem", "reopen_requires": "one new theorem fixing a named provenance atom", "live_frontier_unlocked": False},
        {"lane": "selector/QW-2191 lane", "latest_boundary": "P2699-P2721 and later guardrails", "blocked_by": "missing non-premise strict selector/orientation source", "reopen_requires": "new strict symmetry-breaking provider", "live_frontier_unlocked": False},
        {"lane": "legacy-to-strict bridge and role-transfer lane", "latest_boundary": "kernel split and P2680 source atom closures", "blocked_by": "missing explicit completion map and downstream role-transfer theorem", "reopen_requires": "one new typed bridge/source atom or completion theorem", "live_frontier_unlocked": False},
        {"lane": "L_total/ToE promotion lane", "latest_boundary": "all nonpromotion guards through P2974", "blocked_by": "missing strict source packets and nonproxy variational/coupling theorems", "reopen_requires": "all prerequisite strict exports, not a replay", "live_frontier_unlocked": False},
    ]


def acceptance_matrix() -> list[dict[str, Any]]:
    full = (1 << len(ACCEPTANCE_ATOMS)) - 1
    return [{"mask": m, "present": {atom: bool(m & (1 << i)) for i, atom in enumerate(ACCEPTANCE_ATOMS)}, "accepts_fresh_strict_typed_object": m == full} for m in range(1 << len(ACCEPTANCE_ATOMS))]


def build_payload(p2974_path: Any) -> dict[str, Any]:
    candidates = candidate_rows()
    lanes = lane_rows()
    matrix = acceptance_matrix()
    return {
        "status": "P2975_POST_INCIDENCE_FRESH_TYPED_OBJECT_INTAKE_GATE_NO_NEW_LIVE_FRONTIER",
        "input_hashes": {"P2974": hashlib.sha256(p2974_path.read_bytes()).hexdigest() if p2974_path.exists() else None},
        "constructed_theoretical_objects": {
            "candidate_object": "PostIncidenceFreshStrictTypedObjectIntakeGate",
            "acceptance_atoms": ACCEPTANCE_ATOMS,
            "candidate_intake_rows": candidates,
            "state_map_lane_rows": lanes,
            "finite_acceptance_matrix": matrix,
        },
        "intake_certificate": {
            "candidate_count": len(candidates),
            "accepted_candidate_count": sum(1 for row in candidates if row["accepted_current_new_live_frontier"]),
            "lane_count": len(lanes),
            "unlocked_lane_count": sum(1 for row in lanes if row["live_frontier_unlocked"]),
            "acceptance_matrix_rows": len(matrix),
            "accepted_matrix_rows": sum(1 for row in matrix if row["accepts_fresh_strict_typed_object"]),
            "no_new_live_frontier_certificate_exported": True,
        },
        "decision": {
            "positive_progress": "P2975 constructs the post-incidence fresh typed-object intake gate and audits seven content-defined follow-up routes rather than replaying P2971-P2974.",
            "breakthrough": "No new live frontier is unlocked: every current candidate is either incidence replay, ratio/Gamma/selector/bridge replay, or an unsupplied placeholder.",
            "negative_export_flags": {k: False for k in ["fresh_strict_typed_object_accepted", "new_live_frontier_unlocked", "strict_source_theorem_exported", "nonproxy_ltotal_exported", "bridge_closure_exported", "role_transfer_exported", "toe_closure_exported"]},
            "next_honest_step": "Do not replay the incidence theorem triad, ratio-package scalar/unit/k arithmetic, Gamma/Jacobian readiness, selector/QW-2191 replay, generic bridge replay, or L_total/ToE promotion.  The next proof-grade move must supply exactly one concrete new strict typed object/theorem/provider outside these closed lanes, with a formula or artifact and a finite witness test; otherwise preserve the P2929-P2975 no-new-live-frontier boundary.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    cert = payload["intake_certificate"]
    lines = [
        "# P2975/S1925 post-incidence fresh typed-object intake gate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Intake certificate",
        f"- candidates / accepted: `{cert['candidate_count']}` / `{cert['accepted_candidate_count']}`",
        f"- lanes / unlocked: `{cert['lane_count']}` / `{cert['unlocked_lane_count']}`",
        f"- acceptance matrix rows/accepted: `{cert['acceptance_matrix_rows']}/{cert['accepted_matrix_rows']}`",
        f"- no-new-live-frontier certificate exported: `{cert['no_new_live_frontier_certificate_exported']}`",
        "",
        "## Lay summary",
        payload["decision"]["positive_progress"],
        payload["decision"]["breakthrough"],
        "",
        "## Recommendation",
        payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    read_json(P2974)
    payload = build_payload(P2974)
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2975/S1925 post-incidence fresh typed-object intake gate", "## P2975/S1925 post-incidence fresh typed-object intake gate\n\n`P2975/S1925` constructs the post-incidence fresh typed-object intake gate after P2974.  Seven content-defined candidate routes are tested: incidence derivative replay, incidence source/coupling replay, ratio-package scalar/unit/k replay, Gamma/P2912 Jacobian import replay, selector/`QW-2191` replay, generic legacy-to-strict bridge replay, and an unsupplied fresh-object placeholder.  The finite gate has `64` acceptance profiles and only the all-atoms profile accepts; current accepted candidates are `0`, and unlocked lanes are `0`.  Thus no fresh strict typed object, nonproxy `L_total`, bridge closure, role transfer, or ToE follows.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2975/S1925 post-incidence intake `L_total` guard", "## P2975/S1925 post-incidence intake `L_total` guard\n\n`P2975/S1925` adds no variational or sourced term to `L_total`.  It is an intake/no-new-live-frontier gate showing that incidence replay, ratio-package replay, Gamma/Jacobian import, selector replay, generic bridge replay, and unsupplied placeholders do not provide a new strict typed object/theorem/provider.  Therefore no EOM, Hamiltonian, bridge closure, role transfer, or ToE is exported.\n")
    append_once(AGENTS, "Current post-incidence fresh typed-object intake guardrail (P2975/S1925, 2026-06-20)", "## Current post-incidence fresh typed-object intake guardrail (P2975/S1925, 2026-06-20)\n\n- P2975 constructs the post-P2974 fresh typed-object intake gate across content-defined follow-up routes: incidence replay, ratio-package replay, Gamma/Jacobian import, selector/`QW-2191` replay, generic bridge replay, and unsupplied fresh-object placeholders.\n- The finite gate has `7` candidates, `6` state-map lanes, `64` acceptance profiles, `0` accepted current candidates, and `0` unlocked lanes.\n- Do not promote incidence replay, ratio-package scalar/unit/k arithmetic, Gamma/Jacobian readiness, selector replay, generic bridge replay, or placeholders to strict sourcehood, nonproxy `L_total`, bridge closure, role transfer, or ToE.\n- A next admissible move must supply exactly one concrete new strict typed object/theorem/provider outside these closed lanes, with a formula or artifact and a finite witness test; otherwise preserve the P2929-P2975 no-new-live-frontier boundary.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
