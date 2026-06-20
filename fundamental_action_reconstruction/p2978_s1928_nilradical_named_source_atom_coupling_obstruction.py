#!/usr/bin/env python3
"""P2978/S1928: nilradical named source-atom coupling obstruction.

P2977 left one concrete nilradical route: couple the Z12 nilradical object to a
named missing source atom.  This audit attacks that route only.  It does not
replay nilradical provenance, ratio/Gamma/incidence scans, selector closure,
generic bridge completion, role transfer, or formal L_total placeholders.

The finite calculation treats the nonzero nilpotent 6 as an invariant scalar
anchor and tests whether it supplies a coupling witness for four named atoms:
selector/orientation sign, positive damping beta/Z_beta, legacy-to-strict bridge
source, and action-density/variational source.  The scalar is unit-fixed and
nilpotent, but it is orientation-blind, nonpositive-scale-unsourced, bridge-law
silent, and action-measure silent.
"""
from __future__ import annotations

import hashlib, json
from itertools import product
from math import gcd
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p2977_s1927_nilradical_strict_nadsoliton_provenance_obstruction import OUT as P2977

OUT = GEN / "p2978_s1928_nilradical_named_source_atom_coupling_obstruction.json"
MD = GEN / "p2978_s1928_nilradical_named_source_atom_coupling_obstruction.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

MODULUS = 12
NILPOTENT = 6
UNITS = [x for x in range(MODULUS) if gcd(x, MODULUS) == 1]
NAMED_SOURCE_ATOMS = [
    "selector_orientation_sign_atom",
    "positive_beta_Z_beta_damping_atom",
    "legacy_to_strict_bridge_source_atom",
    "action_density_variational_source_atom",
]


def nilpotent_scalar_witness() -> dict[str, Any]:
    unit_images = {str(u): (u * NILPOTENT) % MODULUS for u in UNITS}
    signed_orientation_scores = {"+omega": NILPOTENT, "-omega": NILPOTENT}
    beta_ratios = {"n_over_12": NILPOTENT / MODULUS, "square_over_12": (NILPOTENT * NILPOTENT % MODULUS) / MODULUS}
    return {
        "modulus": MODULUS,
        "nilpotent": NILPOTENT,
        "nilpotent_square_mod_12": (NILPOTENT * NILPOTENT) % MODULUS,
        "units": UNITS,
        "unit_images": unit_images,
        "unit_fixed": all(v == NILPOTENT for v in unit_images.values()),
        "signed_orientation_scores": signed_orientation_scores,
        "orientation_score_gap": signed_orientation_scores["+omega"] - signed_orientation_scores["-omega"],
        "beta_ratio_samples": beta_ratios,
        "exports_positive_beta_scale": False,
        "exports_bridge_completion_law": False,
        "exports_action_measure": False,
    }


def source_atom_rows(witness: dict[str, Any]) -> list[dict[str, Any]]:
    return [
        {
            "atom": "selector_orientation_sign_atom",
            "finite_nilradical_witness": True,
            "named_atom_targeted": True,
            "strict_provenance_available": False,
            "orientation_sensitive_coupling": witness["orientation_score_gap"] != 0,
            "positive_scale_or_measure_source": False,
            "bridge_completion_law": False,
            "action_variational_installation": False,
            "accepted_current_coupling": False,
            "obstruction": "unit-fixed nilpotent gives identical scores to +omega and -omega",
        },
        {
            "atom": "positive_beta_Z_beta_damping_atom",
            "finite_nilradical_witness": True,
            "named_atom_targeted": True,
            "strict_provenance_available": False,
            "orientation_sensitive_coupling": False,
            "positive_scale_or_measure_source": False,
            "bridge_completion_law": False,
            "action_variational_installation": False,
            "accepted_current_coupling": False,
            "obstruction": "n/12 is a formal ratio sample and n^2=0; no target-independent positive beta/Z_beta law is exported",
        },
        {
            "atom": "legacy_to_strict_bridge_source_atom",
            "finite_nilradical_witness": True,
            "named_atom_targeted": True,
            "strict_provenance_available": False,
            "orientation_sensitive_coupling": False,
            "positive_scale_or_measure_source": False,
            "bridge_completion_law": False,
            "action_variational_installation": False,
            "accepted_current_coupling": False,
            "obstruction": "nilpotent algebra does not supply amplitude, phase/topological, or nonlinear damping completion map",
        },
        {
            "atom": "action_density_variational_source_atom",
            "finite_nilradical_witness": True,
            "named_atom_targeted": True,
            "strict_provenance_available": False,
            "orientation_sensitive_coupling": False,
            "positive_scale_or_measure_source": False,
            "bridge_completion_law": False,
            "action_variational_installation": False,
            "accepted_current_coupling": False,
            "obstruction": "no unit-bearing density, measure, field variable, or variational chain is attached",
        },
        {
            "atom": "completed_strict_nilradical_named_atom_coupling_schema",
            "finite_nilradical_witness": True,
            "named_atom_targeted": True,
            "strict_provenance_available": True,
            "orientation_sensitive_coupling": True,
            "positive_scale_or_measure_source": True,
            "bridge_completion_law": True,
            "action_variational_installation": True,
            "accepted_current_coupling": False,
            "obstruction": "schema row only; no current artifact exports the required theorem package",
        },
    ]


def coupling_obligations(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    current = [r for r in rows if r["atom"] != "completed_strict_nilradical_named_atom_coupling_schema"]
    return [
        {"obligation": "finite_nilradical_witness", "satisfied": all(r["finite_nilradical_witness"] for r in current), "evidence": "the same computed nilpotent 6 and 6^2=0 witness is used for every atom row"},
        {"obligation": "exactly_named_source_atom_targeted", "satisfied": all(r["named_atom_targeted"] for r in current) and len(current) == len(NAMED_SOURCE_ATOMS), "evidence": "selector/orientation, beta/Z_beta, bridge-source, and action-density atoms are enumerated explicitly"},
        {"obligation": "strict_provenance_available", "satisfied": any(r["strict_provenance_available"] for r in current), "evidence": "P2977 found no strict nilradical provenance"},
        {"obligation": "orientation_sensitive_coupling", "satisfied": any(r["orientation_sensitive_coupling"] for r in current), "evidence": "nilpotent scalar scores +omega and -omega equally"},
        {"obligation": "positive_scale_or_measure_source", "satisfied": any(r["positive_scale_or_measure_source"] for r in current), "evidence": "formal n/12 does not export a target-independent beta/Z_beta or measure source"},
        {"obligation": "bridge_completion_law", "satisfied": any(r["bridge_completion_law"] for r in current), "evidence": "no amplitude/phase/damping completion map is supplied"},
        {"obligation": "action_variational_installation", "satisfied": any(r["action_variational_installation"] for r in current), "evidence": "no unit-bearing named density or variational chain is installed"},
        {"obligation": "accepted_current_coupling", "satisfied": any(r["accepted_current_coupling"] for r in current), "evidence": "all current rows fail at least one strict coupling premise"},
    ]


def acceptance_matrix() -> list[dict[str, Any]]:
    names = ["finite_witness", "named_atom", "strict_provenance", "orientation_sensitive", "positive_scale_measure", "bridge_law", "action_variational"]
    rows = []
    for bits in product([False, True], repeat=len(names)):
        present = dict(zip(names, bits))
        rows.append({"present": present, "accepts_named_source_atom_coupling": all(bits)})
    return rows


def build_payload(p2977_path: Any) -> dict[str, Any]:
    witness = nilpotent_scalar_witness()
    rows = source_atom_rows(witness)
    obligations = coupling_obligations(rows)
    matrix = acceptance_matrix()
    return {
        "status": "P2978_NILRADICAL_NAMED_SOURCE_ATOM_COUPLING_OBSTRUCTION_BOUNDED_NO_GO",
        "input_hashes": {"P2977": hashlib.sha256(p2977_path.read_bytes()).hexdigest() if p2977_path.exists() else None},
        "constructed_theoretical_objects": {
            "object": "NilradicalNamedSourceAtomCoupling_ObstructionMatrix",
            "nilpotent_scalar_witness": witness,
            "named_source_atom_rows": rows,
            "proof_obligation_rows": obligations,
            "finite_acceptance_matrix": matrix,
        },
        "coupling_certificate": {
            "named_source_atom_count": len(NAMED_SOURCE_ATOMS),
            "candidate_row_count": len(rows),
            "unit_fixed": witness["unit_fixed"],
            "nilpotent_square_mod_12": witness["nilpotent_square_mod_12"],
            "orientation_score_gap": witness["orientation_score_gap"],
            "accepted_current_couplings": [r["atom"] for r in rows if r["accepted_current_coupling"]],
            "acceptance_matrix_rows": len(matrix),
            "accepted_rows": sum(1 for r in matrix if r["accepts_named_source_atom_coupling"]),
        },
        "decision": {
            "positive_progress": "P2978 attacks exactly one remaining nilradical route: coupling the finite nilpotent anchor 6 to named missing source atoms.",
            "breakthrough": "Bounded no-go: the nilpotent anchor is unit-fixed and finite, but it is orientation-blind, does not export target-independent positive beta/Z_beta or measure, supplies no legacy-to-strict completion law, and installs no unit-bearing action/variational density.",
            "negative_export_flags": {k: False for k in ["named_source_atom_coupling_exported", "strict_nilradical_source_exported", "selector_or_orientation_exported", "damping_source_exported", "bridge_closure_exported", "action_density_exported", "nonproxy_ltotal_exported", "role_transfer_exported", "toe_closure_exported"]},
            "next_honest_step": "Do not replay nilradical source-atom couplings, unit-fixity, n/12 ratio samples, orientation-blind scalar scoring, generic bridge maps, or formal action slots.  The only remaining nilradical-lane proof-grade route is action/variational installation with a genuinely unit-bearing named density and nonproxy chain; if that cannot be supplied, pivot to a new strict typed object outside the nilradical lane and preserve the P2929-P2978 no-strict-export boundary.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    cert = payload["coupling_certificate"]
    lines = [
        "# P2978/S1928 nilradical named source-atom coupling obstruction",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Coupling certificate",
        f"- named source atoms: `{cert['named_source_atom_count']}`",
        f"- candidate rows: `{cert['candidate_row_count']}`",
        f"- unit-fixed nilpotent: `{cert['unit_fixed']}`",
        f"- nilpotent square mod 12: `{cert['nilpotent_square_mod_12']}`",
        f"- orientation score gap: `{cert['orientation_score_gap']}`",
        f"- accepted current couplings: `{cert['accepted_current_couplings']}`",
        f"- acceptance matrix rows/accepted: `{cert['acceptance_matrix_rows']}/{cert['accepted_rows']}`",
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
    read_json(P2977)
    payload = build_payload(P2977)
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2978/S1928 nilradical named source-atom coupling obstruction", "## P2978/S1928 nilradical named source-atom coupling obstruction\n\n`P2978/S1928` attacks exactly one P2977 remaining route: coupling the `Z/12Z` nilradical anchor `6` to named missing source atoms.  The finite matrix audits selector/orientation sign, positive `beta/Z_beta` damping, legacy-to-strict bridge-source, and action-density/variational atoms.  The result is bounded no-go: `6` is unit-fixed and `6^2=0`, but its scalar score is orientation-blind (`+omega` and `-omega` gap `0`), formal `6/12` is not a target-independent positive scale/measure source, no amplitude/phase/damping completion law is supplied, and no unit-bearing named action density or nonproxy variational chain is installed.  No named source-atom coupling, damping source, selector/orientation source, bridge closure, action-density export, nonproxy `L_total`, role transfer, or ToE follows.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2978/S1928 nilradical source-atom coupling `L_total` guard", "## P2978/S1928 nilradical source-atom coupling `L_total` guard\n\n`P2978/S1928` does not add a nilradical coupling term to `L_total`.  The named source-atom audit finds no strict provenance, positive scale/measure source, bridge law, unit-bearing action density, or nonproxy variational chain, so EOM, Hamiltonian, role transfer, bridge closure, and ToE remain unexported.\n")
    append_once(AGENTS, "Current nilradical named source-atom coupling obstruction guardrail (P2978/S1928, 2026-06-20)", "## Current nilradical named source-atom coupling obstruction guardrail (P2978/S1928, 2026-06-20)\n\n- P2978 attacks exactly one remaining nilradical route: named source-atom coupling for the unit-fixed nilpotent `6`.\n- The finite matrix audits `4` named atoms: selector/orientation sign, positive `beta/Z_beta` damping, legacy-to-strict bridge-source, and action-density/variational source; the acceptance matrix has `128` profiles with only the full profile accepting.\n- The current route is bounded no-go: the nilpotent scalar is orientation-blind, formal `6/12` is not a target-independent scale/measure source, no bridge-completion law is supplied, and no unit-bearing named action density or nonproxy variational chain is installed.\n- Do not promote nilradical source-atom rows, unit-fixity, ratio samples, orientation-blind scalar scores, generic bridge maps, or formal action slots to selector closure, damping source, bridge completion, role transfer, nonproxy `L_total`, or ToE.  The only remaining nilradical-lane route is action/variational installation with a genuinely unit-bearing named density and nonproxy chain; otherwise pivot to a new strict typed object while preserving the P2929-P2978 no-strict-export boundary.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
