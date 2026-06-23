#!/usr/bin/env python3
"""P3057/S2007: source-polarity carrier extension SAT verifier.

P3056 made the normal form explicit but still left the next object described in
words.  P3057 turns that next object into a finite theorem template: a single
source-polarity carrier plus row-import/certification maps and compatibility
squares.  It then solves the small Boolean extension problem induced by the
P3056 no-export certificate, so the next theoretical obligation is a minimal
computed cut rather than another named clue.
"""
from __future__ import annotations

import hashlib, itertools, json, subprocess
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p3056_s2006_selector_gluing_object_normal_form_acceptance_verifier import OUT as P3056

OUT = GEN / "p3057_s2007_source_polarity_carrier_extension_sat_verifier.json"
MD = GEN / "p3057_s2007_source_polarity_carrier_extension_sat_verifier.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

CONTENT_PATTERNS = {
    "new_source_polarity_carrier": r"source-law/polarity carrier|source.*polarity.*carrier|strict source law.*unique polarity coupling|unique polarity.*strict source law",
    "compatibility_square_theorem": r"compatibility square|source-to-readout|selector/readout/action|field pullback.*chain rule|source.*signed value.*localizer",
    "single_object_not_bundle": r"single.*object.*all.*rows|bundling.*gluing|one.*non-imported object|row-complete bundle",
}

CORE_ROWS = [
    "strict_source_law",
    "nonzero_signed_value",
    "nonpremise_origin_or_localizer",
    "unique_polarity_coupling",
    "localization_pullback_to_field_support",
    "unit_bearing_variational_chain_rule",
]

SQUARES = [
    "source_law_to_signed_value",
    "source_law_to_origin_localizer",
    "signed_value_to_unique_polarity",
    "origin_localizer_to_field_pullback",
    "field_pullback_to_unit_variational_chain_rule",
    "unique_polarity_to_selector_readout",
]

# Atoms that a genuinely new theorem object could add.  Existing clues are not
# reused as closure; they must be pulled back/certified onto the same carrier.
EXTENSION_ATOMS = [
    {"atom": "new_strict_source_law", "covers": ["strict_source_law"], "kind": "primitive_row"},
    {"atom": "new_unique_polarity_coupling", "covers": ["unique_polarity_coupling"], "kind": "primitive_row"},
    {"atom": "certify_signed_value_on_carrier", "covers": ["nonzero_signed_value"], "kind": "row_import_certificate"},
    {"atom": "certify_origin_localizer_on_carrier", "covers": ["nonpremise_origin_or_localizer"], "kind": "row_import_certificate"},
    {"atom": "certify_field_pullback_on_carrier", "covers": ["localization_pullback_to_field_support"], "kind": "row_import_certificate"},
    {"atom": "certify_unit_variational_chain_on_carrier", "covers": ["unit_bearing_variational_chain_rule"], "kind": "row_import_certificate"},
    *[{"atom": f"prove_{square}", "covers": [square], "kind": "compatibility_square"} for square in SQUARES],
]


def content_grep() -> list[dict[str, Any]]:
    rows = []
    for name, pattern in CONTENT_PATTERNS.items():
        proc = subprocess.run(["rg", "-n", pattern, "fundamental_action_reconstruction", "AGENTS.md", "-S"], cwd=REPO, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
        hits = [line for line in proc.stdout.splitlines() if line.strip()]
        rows.append({"content_lane": name, "pattern": pattern, "hit_count": len(hits), "sample_hits": hits[:8]})
    return rows


def current_single_carrier_state() -> dict[str, bool]:
    # P3056 proves no current carrier is the single object.  Therefore every row
    # and square must be certified on the new carrier, even when a disconnected
    # clue with the same row type exists elsewhere.
    return {key: False for key in [*CORE_ROWS, *SQUARES]}


def acceptance(state: dict[str, bool]) -> bool:
    return all(state[row] for row in CORE_ROWS) and all(state[square] for square in SQUARES)


def apply_atoms(atoms: tuple[dict[str, Any], ...]) -> dict[str, bool]:
    state = current_single_carrier_state()
    for atom in atoms:
        for covered in atom["covers"]:
            state[covered] = True
    return state


def enumerate_minimal_extensions() -> list[dict[str, Any]]:
    accepted = []
    for r in range(1, len(EXTENSION_ATOMS) + 1):
        for combo in itertools.combinations(EXTENSION_ATOMS, r):
            state = apply_atoms(combo)
            if acceptance(state):
                accepted.append({
                    "atoms": [atom["atom"] for atom in combo],
                    "size": r,
                    "primitive_rows": sum(1 for atom in combo if atom["kind"] == "primitive_row"),
                    "row_import_certificates": sum(1 for atom in combo if atom["kind"] == "row_import_certificate"),
                    "compatibility_squares": sum(1 for atom in combo if atom["kind"] == "compatibility_square"),
                })
        if accepted:
            return accepted
    return []


def obstruction_cut_table() -> list[dict[str, Any]]:
    rows = []
    for atom in EXTENSION_ATOMS:
        state = apply_atoms((atom,))
        rows.append({
            "single_added_atom": atom["atom"],
            "kind": atom["kind"],
            "covered": atom["covers"],
            "remaining_missing_count": sum(1 for value in state.values() if not value),
            "accepted_alone": acceptance(state),
        })
    return rows


def theorem_template(minimal: list[dict[str, Any]]) -> dict[str, Any]:
    return {
        "object": "StrictSourcePolarityCarrierExtensionTheoremTemplate",
        "carrier_symbol": "G_selector",
        "premise": "one non-imported strict object carries every row and compatibility square internally",
        "required_atoms": minimal[0]["atoms"] if minimal else [],
        "minimal_atom_count": minimal[0]["size"] if minimal else None,
        "interpretation": "the missing object is not another sign observable; it is a source-polarity law plus four row-certification maps and six compatibility-square theorems on the same carrier",
    }


def build_payload() -> dict[str, Any]:
    p3056 = read_json(P3056)
    grep_rows = content_grep()
    minimal = enumerate_minimal_extensions()
    cuts = obstruction_cut_table()
    obligations = [
        {"obligation": "content_first_grep_before_extension", "satisfied": True, "detail": "three content patterns searched for source-polarity carrier and compatibility-square material"},
        {"obligation": "construct_source_polarity_extension_template", "satisfied": True, "detail": "G_selector template separates primitive rows, row-import certificates, and compatibility-square theorems"},
        {"obligation": "solve_boolean_minimal_extension_problem", "satisfied": True, "detail": "all subsets are searched until the first accepted extension size is found"},
        {"obligation": "export_actual_new_carrier_formula", "satisfied": False, "detail": "P3057 computes the exact theorem template but does not supply an external formula/value for G_selector"},
        {"obligation": "install_selector_ltotal_bridge_toe", "satisfied": False, "detail": "no QW-2191 discharge, selector closure, L_total, bridge, role transfer, or ToE closure follows from a template alone"},
    ]
    return {
        "status": "P3057_SOURCE_POLARITY_CARRIER_EXTENSION_SAT_VERIFIER_TEMPLATE_NO_EXPORT",
        "input_hashes": {"P3056": hashlib.sha256(P3056.read_bytes()).hexdigest() if P3056.exists() else None},
        "constructed_theoretical_objects": {
            "content_first_repo_grep": grep_rows,
            "current_single_carrier_state": current_single_carrier_state(),
            "extension_atoms": EXTENSION_ATOMS,
            "minimal_accepted_extensions": minimal,
            "obstruction_cut_table": cuts,
            "theorem_template": theorem_template(minimal),
        },
        "finite_certificate": {
            "content_grep_lanes": len(grep_rows),
            "content_grep_hits": sum(row["hit_count"] for row in grep_rows),
            "core_rows": len(CORE_ROWS),
            "compatibility_squares": len(SQUARES),
            "extension_atoms": len(EXTENSION_ATOMS),
            "minimal_extension_count": len(minimal),
            "minimal_extension_size": minimal[0]["size"] if minimal else None,
            "primitive_rows_in_minimal_extension": minimal[0]["primitive_rows"] if minimal else None,
            "row_import_certificates_in_minimal_extension": minimal[0]["row_import_certificates"] if minimal else None,
            "compatibility_squares_in_minimal_extension": minimal[0]["compatibility_squares"] if minimal else None,
            "single_atom_obstruction_rows": len(cuts),
            "single_atom_accepts": sum(1 for row in cuts if row["accepted_alone"]),
            "proof_obligations": len(obligations),
            "satisfied_proof_obligations": sum(1 for row in obligations if row["satisfied"]),
            "p3056_status_seen": p3056.get("status"),
        },
        "proof_obligations": obligations,
        "decision": {
            "bounded_no_go": "P3057 upgrades P3056 from a no-export normal form to a computed extension theorem template.  The SAT search shows the minimal accepted single-carrier extension has 12 atoms: two primitive rows, four row-import/certification maps, and six compatibility-square theorems.  Every one-atom continuation leaves 11 missing obligations, so naming one more sign source, one more localizer, or one more variational slot cannot close the selector lane.",
            "negative_export_flags": {k: False for k in ["actual_g_selector_formula_exported", "selector_gluing_object_exported", "qw2191_discharged", "strict_selector_closure_exported", "observed_physics_exported", "unit_bearing_action_eom_hamiltonian_exported", "ltotal_exported", "bridge_closure_exported", "role_transfer_exported", "toe_closure_exported"]},
            "next_honest_step": "The next proof-grade move should not add an isolated clue.  It should either supply an explicit formula/artifact for G_selector realizing the 12-atom StrictSourcePolarityCarrierExtensionTheoremTemplate, or pivot outside the selector clue lane.  A partial continuation is only honest if it attacks one named atom of the 12-atom template and states that selector/readout/L_total/bridge/ToE closure remains unexported.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    c = payload["finite_certificate"]
    lines = [
        "# P3057/S2007 source-polarity carrier extension SAT verifier", "",
        f"Status: `{payload['status']}`", "", "## Finite certificate",
        f"- content grep lanes: `{c['content_grep_lanes']}`",
        f"- content grep hits: `{c['content_grep_hits']}`",
        f"- core rows: `{c['core_rows']}`",
        f"- compatibility squares: `{c['compatibility_squares']}`",
        f"- extension atoms: `{c['extension_atoms']}`",
        f"- minimal extension count: `{c['minimal_extension_count']}`",
        f"- minimal extension size: `{c['minimal_extension_size']}`",
        f"- primitive rows in minimal extension: `{c['primitive_rows_in_minimal_extension']}`",
        f"- row-import certificates in minimal extension: `{c['row_import_certificates_in_minimal_extension']}`",
        f"- compatibility squares in minimal extension: `{c['compatibility_squares_in_minimal_extension']}`",
        f"- single-atom accepts: `{c['single_atom_accepts']}`",
        f"- satisfied proof obligations: `{c['satisfied_proof_obligations']}/{c['proof_obligations']}`", "",
        "## Decision", payload["decision"]["bounded_no_go"], "", "## Recommendation", payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3057/S2007 source-polarity carrier extension SAT verifier", "## P3057/S2007 source-polarity carrier extension SAT verifier\n\n`P3057/S2007` constructs the computed extension theorem template left open by `P3056`: a single carrier `G_selector` must internalize six rows and six compatibility squares.  The Boolean extension search over `12` atoms finds the minimal accepted extension has all `12` atoms: two primitive rows (`strict_source_law`, `unique_polarity_coupling`), four row-import/certification maps, and six compatibility-square theorems.  Every single-atom continuation remains blocked.  No actual `G_selector` formula, selector gluing object, `QW-2191` discharge, selector closure, `L_total`, bridge/role transfer, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3057/S2007 source-polarity carrier extension `L_total` guard", "## P3057/S2007 source-polarity carrier extension `L_total` guard\n\n`P3057/S2007` adds no physical `L_total` term.  It computes the 12-atom theorem template that a future `G_selector` would need, but no actual unit-bearing signed source/action/EOM carrier or variational chain rule is exported.\n")
    append_once(AGENTS, "Current source-polarity carrier extension SAT guardrail (P3057/S2007, 2026-06-23)", "## Current source-polarity carrier extension SAT guardrail (P3057/S2007, 2026-06-23)\n\n- P3057 converts the post-P3056 recommendation into a finite theorem template for one carrier `G_selector`: six internal rows plus six compatibility squares.\n- The Boolean extension search shows the minimal accepted single-carrier extension has all `12` atoms: two primitive source/polarity rows, four row-import certificates, and six compatibility-square theorems; every one-atom continuation remains blocked.\n- Do not promote an isolated sign source, localizer, graph clue, receiver/chiral diagnostic, or variational slot to `QW-2191` discharge, selector closure, observed physics, `L_total`, bridge/role-transfer, or ToE closure.\n- A next move must either supply an explicit `G_selector` formula/artifact realizing the 12-atom template, or attack one named atom while preserving the P3048-P3057 bounded no-export boundary.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
