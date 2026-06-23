#!/usr/bin/env python3
"""P3058/S2008: unique-polarity coupling constraint verifier.

P3057 showed that a future G_selector needs a named atom
`new_unique_polarity_coupling`.  P3058 attacks that one atom directly.  It
models the two already-known sign torsors as a finite two-polarity coupling
problem and asks whether the current compatibility constraints force exactly
one polarity without importing a premise selector.
"""
from __future__ import annotations

import hashlib, itertools, json, subprocess
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p3057_s2007_source_polarity_carrier_extension_sat_verifier import OUT as P3057

OUT = GEN / "p3058_s2008_unique_polarity_coupling_constraint_verifier.json"
MD = GEN / "p3058_s2008_unique_polarity_coupling_constraint_verifier.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

CONTENT_PATTERNS = {
    "unique_polarity_coupling_constraints": r"unique polarity coupling|polarity.*coupling|coupling polarity|source.*polarity",
    "two_polarity_torsor_blindness": r"opposite polarity|two equivariant couplings|orientation-blind|lambda.*unfixed|polarity.*choice",
    "nonpremise_selector_forbidden": r"non-premise.*selector|QW-2191.*open|selector closure.*unexported|premise selector",
}

POLARITIES = ["plus_polarity", "minus_polarity"]

# Current constraints that may be reused without adding a new source law.  Each
# is intentionally expressed as a predicate on the coupling polarity.  A real
# unique-polarity atom would need at least one strict exported polarity-odd law.
CONSTRAINTS = [
    {"name": "aut_equivariance_of_sign_torsor", "kind": "structural", "accepts": POLARITIES, "polarity_odd": False},
    {"name": "translation_orbit_neutrality", "kind": "localizer_absence", "accepts": POLARITIES, "polarity_odd": False},
    {"name": "bispectrum_sign_readout_compatibility", "kind": "readout", "accepts": POLARITIES, "polarity_odd": False},
    {"name": "field_pullback_support_typing", "kind": "installation", "accepts": POLARITIES, "polarity_odd": False},
    {"name": "unit_positive_variational_slot", "kind": "variational", "accepts": POLARITIES, "polarity_odd": False},
]

MISSING_SOURCE_ATOM = {
    "name": "strict_polarity_odd_source_law_boundary_condition",
    "kind": "not_exported_current_artifact",
    "would_accept": ["plus_polarity"],
    "reason_not_used": "using it would be exactly the new_unique_polarity_coupling atom, not a consequence of current artifacts",
}


def content_grep() -> list[dict[str, Any]]:
    rows = []
    for lane, pattern in CONTENT_PATTERNS.items():
        proc = subprocess.run(["rg", "-n", pattern, "fundamental_action_reconstruction", "AGENTS.md", "-S"], cwd=REPO, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
        hits = [line for line in proc.stdout.splitlines() if line.strip()]
        rows.append({"content_lane": lane, "pattern": pattern, "hit_count": len(hits), "sample_hits": hits[:8]})
    return rows


def accepted_by(constraint_subset: tuple[dict[str, Any], ...]) -> list[str]:
    accepted = set(POLARITIES)
    for constraint in constraint_subset:
        accepted &= set(constraint["accepts"])
    return sorted(accepted)


def enumerate_constraint_intersections() -> list[dict[str, Any]]:
    rows = []
    for r in range(1, len(CONSTRAINTS) + 1):
        for combo in itertools.combinations(CONSTRAINTS, r):
            accepted = accepted_by(combo)
            rows.append({
                "constraints": [c["name"] for c in combo],
                "size": r,
                "accepted_polarities": accepted,
                "accepted_count": len(accepted),
                "unique": len(accepted) == 1,
                "polarity_odd_constraint_count": sum(1 for c in combo if c["polarity_odd"]),
            })
    return rows


def construct_missing_object_normal_form() -> dict[str, Any]:
    return {
        "object": "UniquePolarityCouplingAtomNormalForm",
        "target_atom": "new_unique_polarity_coupling",
        "carrier_symbol": "G_selector",
        "source_torsor": "chiral/readout sign torsor {-2,+2}",
        "target_torsor": "orientation selector torsor {-omega,+omega}",
        "current_candidate_couplings": ["plus_polarity", "minus_polarity"],
        "acceptance_rule": "accept only if current strict constraints leave exactly one coupling polarity and include a non-premise polarity-odd source law",
        "missing_source_atom": MISSING_SOURCE_ATOM,
    }


def build_payload() -> dict[str, Any]:
    p3057 = read_json(P3057)
    grep_rows = content_grep()
    intersections = enumerate_constraint_intersections()
    all_current = tuple(CONSTRAINTS)
    all_current_accepted = accepted_by(all_current)
    unique_rows = [row for row in intersections if row["unique"]]
    obligations = [
        {"obligation": "content_first_grep_before_polarity_test", "satisfied": True, "detail": "three content patterns searched for unique polarity coupling and two-polarity obstruction material"},
        {"obligation": "construct_unique_polarity_atom_normal_form", "satisfied": True, "detail": "the missing atom is represented as a two-torsor coupling with an explicit uniqueness acceptance rule"},
        {"obligation": "exhaust_current_constraint_intersections", "satisfied": True, "detail": "all nonempty intersections of five current compatibility constraints are enumerated"},
        {"obligation": "export_strict_polarity_odd_source_law", "satisfied": False, "detail": "current constraints are polarity-even and leave both plus and minus coupling polarities"},
        {"obligation": "discharge_selector_or_ltotal", "satisfied": False, "detail": "no QW-2191 discharge, selector closure, L_total, bridge, role transfer, or ToE closure follows"},
    ]
    return {
        "status": "P3058_UNIQUE_POLARITY_COUPLING_CONSTRAINT_VERIFIER_BOUNDED_NO_EXPORT",
        "input_hashes": {"P3057": hashlib.sha256(P3057.read_bytes()).hexdigest() if P3057.exists() else None},
        "constructed_theoretical_objects": {
            "content_first_repo_grep": grep_rows,
            "unique_polarity_atom_normal_form": construct_missing_object_normal_form(),
            "current_constraints": CONSTRAINTS,
            "constraint_intersections": intersections,
        },
        "finite_certificate": {
            "content_grep_lanes": len(grep_rows),
            "content_grep_hits": sum(row["hit_count"] for row in grep_rows),
            "current_constraints": len(CONSTRAINTS),
            "constraint_intersections": len(intersections),
            "unique_polarity_intersections": len(unique_rows),
            "all_current_constraints_accepted_count": len(all_current_accepted),
            "all_current_constraints_accepted_polarities": all_current_accepted,
            "polarity_odd_current_constraints": sum(1 for c in CONSTRAINTS if c["polarity_odd"]),
            "missing_source_atoms_named": 1,
            "proof_obligations": len(obligations),
            "satisfied_proof_obligations": sum(1 for row in obligations if row["satisfied"]),
            "p3057_status_seen": p3057.get("status"),
        },
        "proof_obligations": obligations,
        "decision": {
            "bounded_no_go": "P3058 attacks exactly one P3057 atom, new_unique_polarity_coupling.  The finite verifier constructs the missing atom as a two-torsor coupling and enumerates all 31 nonempty intersections of five current compatibility constraints.  Every current constraint is polarity-even, every intersection leaves both plus_polarity and minus_polarity, and zero intersections select a unique polarity.  The needed object is therefore sharper: a strict polarity-odd source-law boundary condition coupled to G_selector, not another polarity-blind compatibility check.",
            "negative_export_flags": {k: False for k in ["unique_polarity_coupling_exported", "actual_g_selector_formula_exported", "selector_gluing_object_exported", "qw2191_discharged", "strict_selector_closure_exported", "observed_physics_exported", "unit_bearing_action_eom_hamiltonian_exported", "ltotal_exported", "bridge_closure_exported", "role_transfer_exported", "toe_closure_exported"]},
            "next_honest_step": "Do not repeat polarity-even compatibility constraints.  The next proof-grade move should construct one explicit strict polarity-odd source-law boundary condition for G_selector, with a signed value that selects plus_polarity or minus_polarity non-premise and then proves compatibility with the P3057 row-import and square obligations; otherwise pivot to a different named P3057 atom while preserving no selector/readout/L_total/bridge/ToE export.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    c = payload["finite_certificate"]
    lines = [
        "# P3058/S2008 unique-polarity coupling constraint verifier", "",
        f"Status: `{payload['status']}`", "", "## Finite certificate",
        f"- content grep lanes: `{c['content_grep_lanes']}`",
        f"- content grep hits: `{c['content_grep_hits']}`",
        f"- current constraints: `{c['current_constraints']}`",
        f"- constraint intersections: `{c['constraint_intersections']}`",
        f"- unique-polarity intersections: `{c['unique_polarity_intersections']}`",
        f"- all-current accepted polarities: `{c['all_current_constraints_accepted_polarities']}`",
        f"- polarity-odd current constraints: `{c['polarity_odd_current_constraints']}`",
        f"- satisfied proof obligations: `{c['satisfied_proof_obligations']}/{c['proof_obligations']}`", "",
        "## Decision", payload["decision"]["bounded_no_go"], "", "## Recommendation", payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3058/S2008 unique-polarity coupling constraint verifier", "## P3058/S2008 unique-polarity coupling constraint verifier\n\n`P3058/S2008` attacks exactly one atom from the `P3057` template: `new_unique_polarity_coupling`.  It constructs the missing object as a two-torsor coupling normal form between the chiral/readout sign torsor and the orientation selector torsor, then enumerates all `31` nonempty intersections of five current compatibility constraints.  All current constraints are polarity-even; `0` intersections select a unique polarity, and the full current constraint set still accepts both `plus_polarity` and `minus_polarity`.  The next required object is a strict polarity-odd source-law boundary condition on `G_selector`.  No actual `G_selector`, `QW-2191` discharge, selector closure, `L_total`, bridge/role transfer, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3058/S2008 unique-polarity coupling `L_total` guard", "## P3058/S2008 unique-polarity coupling `L_total` guard\n\n`P3058/S2008` adds no physical `L_total` term.  It tests the `new_unique_polarity_coupling` atom as a finite two-polarity constraint problem and finds that current constraints remain polarity-even; no unit-bearing signed source/action/EOM carrier or nonproxy variational chain rule is exported.\n")
    append_once(AGENTS, "Current unique-polarity coupling constraint guardrail (P3058/S2008, 2026-06-23)", "## Current unique-polarity coupling constraint guardrail (P3058/S2008, 2026-06-23)\n\n- P3058 attacks one named P3057 atom, `new_unique_polarity_coupling`, by constructing a two-torsor coupling normal form for `G_selector`.\n- The finite check enumerates all `31` nonempty intersections of five current compatibility constraints; every constraint is polarity-even, `0` intersections select a unique polarity, and both `plus_polarity` and `minus_polarity` remain accepted under all current constraints.\n- Do not promote polarity-even compatibility, chiral/readout sign agreement, field-pullback typing, or unit-positive variational slots to `QW-2191` discharge, selector closure, observed physics, `L_total`, bridge/role-transfer, or ToE closure.\n- A next move in this lane must construct an explicit strict polarity-odd source-law boundary condition coupled to `G_selector`; otherwise attack a different named P3057 atom while preserving the P3048-P3058 bounded no-export boundary.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
