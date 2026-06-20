#!/usr/bin/env python3
"""P2976/S1926: Z12 nilradical-filtration source candidate.

P2975 requires exactly one concrete new strict typed object/theorem/provider
outside incidence replay, ratio-package scalar/unit/k replay, Gamma/Jacobian
readiness, selector replay, generic bridge replay, and L_total/ToE promotion.
This audit supplies a new finite typed object: the nilradical filtration of the
ring Z/12Z, computed directly from multiplication powers.

The object is real finite algebraic structure: the nilradical is {0, 6}, every
unit fixes 6, and 6^2 = 0.  It is not a strict source theorem: current artifacts
do not export nadsoliton provenance for selecting this filtration, a coupling to
a named missing object, an orientation/selector source, or an action-density/
variational installation.
"""
from __future__ import annotations

import hashlib, json
from math import gcd
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p2975_s1925_post_incidence_fresh_typed_object_intake_gate import OUT as P2975

OUT = GEN / "p2976_s1926_z12_nilradical_filtration_source_candidate.json"
MD = GEN / "p2976_s1926_z12_nilradical_filtration_source_candidate.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

MODULUS = 12
UNITS = [x for x in range(MODULUS) if gcd(x, MODULUS) == 1]


def power_orbit(x: int, limit: int = 8) -> list[int]:
    out = []
    value = x % MODULUS
    for _ in range(1, limit + 1):
        out.append(value)
        value = (value * x) % MODULUS
    return out


def nilpotent_witness_rows() -> list[dict[str, Any]]:
    rows = []
    for x in range(MODULUS):
        orbit = power_orbit(x)
        nilpotent = 0 in orbit
        first_zero_power = next((i + 1 for i, value in enumerate(orbit) if value == 0), None)
        rows.append({"element": x, "power_orbit_prefix": orbit, "nilpotent": nilpotent, "first_zero_power": first_zero_power})
    return rows


def unit_action_on_nilradical(nilradical: list[int]) -> list[dict[str, Any]]:
    return [{"unit": u, "action": {str(x): (u * x) % MODULUS for x in nilradical}} for u in UNITS]


def filtration_certificate() -> dict[str, Any]:
    rows = nilpotent_witness_rows()
    nilradical = [r["element"] for r in rows if r["nilpotent"]]
    nonzero_nilpotents = [x for x in nilradical if x != 0]
    ideals = {"zero_ideal": [0], "nilradical": nilradical, "whole_ring": list(range(MODULUS))}
    action = unit_action_on_nilradical(nilradical)
    return {
        "modulus": MODULUS,
        "units": UNITS,
        "nilpotent_witness_rows": rows,
        "nilradical": nilradical,
        "nonzero_nilpotents": nonzero_nilpotents,
        "nilradical_size": len(nilradical),
        "nonzero_nilpotent_count": len(nonzero_nilpotents),
        "max_first_zero_power": max(r["first_zero_power"] or 0 for r in rows if r["nilpotent"]),
        "ideals_in_filtration": ideals,
        "unit_action_on_nilradical": action,
        "all_units_fix_nonzero_nilpotent": all((u * x) % MODULUS == x for u in UNITS for x in nonzero_nilpotents),
    }


def candidate_rows(cert: dict[str, Any]) -> list[dict[str, Any]]:
    rows = [
        {
            "candidate": "Z12_nilradical_filtration_object",
            "genuinely_new_outside_incidence_lane": True,
            "not_ratio_gamma_selector_bridge_replay": True,
            "finite_witness_computed": True,
            "strict_nadsoliton_provenance_exported": False,
            "couples_to_named_missing_object": False,
            "orientation_or_selector_source_exported": False,
            "action_density_or_variational_installation": False,
            "accepted_current_strict_source_object": False,
            "witness": f"nilradical={cert['nilradical']} with nonzero nilpotent {cert['nonzero_nilpotents']} and max zero power {cert['max_first_zero_power']}",
        },
        {
            "candidate": "unit_fixed_nilpotent_anchor",
            "genuinely_new_outside_incidence_lane": True,
            "not_ratio_gamma_selector_bridge_replay": True,
            "finite_witness_computed": True,
            "strict_nadsoliton_provenance_exported": False,
            "couples_to_named_missing_object": False,
            "orientation_or_selector_source_exported": False,
            "action_density_or_variational_installation": False,
            "accepted_current_strict_source_object": False,
            "witness": "all units fix 6, so the anchor is unit-invariant but orientation/source-blind",
        },
        {
            "candidate": "completed_strict_nilradical_source_schema",
            "genuinely_new_outside_incidence_lane": True,
            "not_ratio_gamma_selector_bridge_replay": True,
            "finite_witness_computed": True,
            "strict_nadsoliton_provenance_exported": True,
            "couples_to_named_missing_object": True,
            "orientation_or_selector_source_exported": True,
            "action_density_or_variational_installation": True,
            "accepted_current_strict_source_object": False,
            "witness": "completed theorem schema only; no current artifact exports provenance/coupling/installation",
        },
    ]
    return rows


def obligation_rows(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    current = [r for r in rows if r["candidate"] != "completed_strict_nilradical_source_schema"]
    return [
        {"obligation": "new_outside_incidence_lane", "satisfied": any(r["genuinely_new_outside_incidence_lane"] for r in current), "evidence": "nilradical filtration is a Z12 ring-theoretic object, not the P2971 incidence object"},
        {"obligation": "not_closed_replay_lane", "satisfied": any(r["not_ratio_gamma_selector_bridge_replay"] for r in current), "evidence": "does not use ratio 9/5, Gamma Jacobian import, selector premise, or bridge replay"},
        {"obligation": "finite_witness_computed", "satisfied": any(r["finite_witness_computed"] for r in current), "evidence": "all 12 element power orbits are enumerated"},
        {"obligation": "strict_nadsoliton_provenance_exported", "satisfied": any(r["strict_nadsoliton_provenance_exported"] for r in current), "evidence": "no theorem states the nadsoliton sources this nilradical filtration"},
        {"obligation": "couples_to_named_missing_object", "satisfied": any(r["couples_to_named_missing_object"] for r in current), "evidence": "no coupling to selector, damping, bridge, or L_total source atom is exported"},
        {"obligation": "orientation_or_selector_source_exported", "satisfied": any(r["orientation_or_selector_source_exported"] for r in current), "evidence": "unit action fixes 6 and gives no directed orientation"},
        {"obligation": "action_density_or_variational_installation", "satisfied": any(r["action_density_or_variational_installation"] for r in current), "evidence": "no unit-bearing density or variational chain is attached"},
        {"obligation": "accepted_current_strict_source_object", "satisfied": any(r["accepted_current_strict_source_object"] for r in current), "evidence": "completed schema remains unavailable"},
    ]


def acceptance_matrix() -> list[dict[str, Any]]:
    names = ["new_outside_incidence", "not_replay", "finite_witness", "strict_provenance", "named_coupling", "orientation_source", "action_variational_installation"]
    full = (1 << len(names)) - 1
    return [{"mask": m, "present": {name: bool(m & (1 << i)) for i, name in enumerate(names)}, "accepts_strict_nilradical_source_object": m == full} for m in range(1 << len(names))]


def build_payload(p2975_path: Any) -> dict[str, Any]:
    cert = filtration_certificate()
    rows = candidate_rows(cert)
    obligations = obligation_rows(rows)
    matrix = acceptance_matrix()
    return {
        "status": "P2976_Z12_NILRADICAL_FILTRATION_SOURCE_CANDIDATE_DEVELOPMENTAL_NO_STRICT_EXPORT",
        "input_hashes": {"P2975": hashlib.sha256(p2975_path.read_bytes()).hexdigest() if p2975_path.exists() else None},
        "constructed_theoretical_objects": {
            "candidate_object": "Z12NilradicalFiltration_SourceCandidate",
            "filtration_certificate": cert,
            "candidate_rows": rows,
            "proof_obligation_rows": obligations,
            "finite_acceptance_matrix": matrix,
        },
        "nilradical_certificate": {
            "modulus": MODULUS,
            "unit_count": len(UNITS),
            "nilradical": cert["nilradical"],
            "nilradical_size": cert["nilradical_size"],
            "nonzero_nilpotents": cert["nonzero_nilpotents"],
            "max_first_zero_power": cert["max_first_zero_power"],
            "all_units_fix_nonzero_nilpotent": cert["all_units_fix_nonzero_nilpotent"],
            "candidate_count": len(rows),
            "accepted_current_strict_source_objects": [r["candidate"] for r in rows if r["accepted_current_strict_source_object"]],
            "acceptance_matrix_rows": len(matrix),
            "accepted_rows": sum(1 for r in matrix if r["accepts_strict_nilradical_source_object"]),
        },
        "decision": {
            "positive_progress": "P2976 supplies a concrete new finite typed object after the P2975 intake gate: the Z12 nilradical filtration, with all power-orbit witnesses computed.",
            "breakthrough": "This is developmental structural progress, not strict closure: the nilradical anchor is finite and unit-invariant, but current artifacts export no nadsoliton provenance, named coupling, orientation source, action-density installation, or variational theorem.",
            "negative_export_flags": {k: False for k in ["strict_nilradical_source_exported", "strict_source_theorem_exported", "selector_or_orientation_exported", "nonproxy_ltotal_exported", "bridge_closure_exported", "role_transfer_exported", "toe_closure_exported"]},
            "next_honest_step": "Do not promote the Z12 nilradical filtration, the unit-fixed element 6, or nilpotent power-orbit witnesses to selector closure, damping source, bridge completion, role transfer, nonproxy L_total, or ToE.  The next proof-grade move must attack exactly one missing theorem for this object: strict nadsoliton provenance, coupling to one named missing source atom, or action/variational installation; otherwise preserve the P2929-P2976 developmental/no-strict-export boundary.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    cert = payload["nilradical_certificate"]
    lines = [
        "# P2976/S1926 Z12 nilradical-filtration source candidate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Nilradical certificate",
        f"- modulus / unit count: `{cert['modulus']}` / `{cert['unit_count']}`",
        f"- nilradical / size: `{cert['nilradical']}` / `{cert['nilradical_size']}`",
        f"- nonzero nilpotents: `{cert['nonzero_nilpotents']}`",
        f"- max first zero power: `{cert['max_first_zero_power']}`",
        f"- all units fix nonzero nilpotent: `{cert['all_units_fix_nonzero_nilpotent']}`",
        f"- accepted current strict source objects: `{cert['accepted_current_strict_source_objects']}`",
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
    read_json(P2975)
    payload = build_payload(P2975)
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2976/S1926 Z12 nilradical-filtration source candidate", "## P2976/S1926 Z12 nilradical-filtration source candidate\n\n`P2976/S1926` supplies one concrete typed object after the P2975 intake gate: the nilradical filtration of `Z/12Z`.  Exhaustive power-orbit enumeration over all `12` elements gives nilradical `{0,6}`, nonzero nilpotent `[6]`, and `6^2=0`; all `4` units fix the nonzero nilpotent.  This is finite algebraic structure outside the incidence lane, but it is developmental only: no strict nadsoliton provenance, named coupling, orientation/selector source, action-density installation, nonproxy `L_total`, bridge closure, role transfer, or ToE follows.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2976/S1926 nilradical-filtration `L_total` guard", "## P2976/S1926 nilradical-filtration `L_total` guard\n\n`P2976/S1926` adds no sourced term to `L_total`.  The Z12 nilradical filtration and the unit-fixed nilpotent `6` are finite algebraic witnesses only; without strict nadsoliton provenance, a named coupling theorem, orientation/source theorem, and action-density or variational installation, they cannot enter EOM, Hamiltonian, bridge closure, role transfer, or ToE.\n")
    append_once(AGENTS, "Current Z12 nilradical-filtration source-candidate guardrail (P2976/S1926, 2026-06-20)", "## Current Z12 nilradical-filtration source-candidate guardrail (P2976/S1926, 2026-06-20)\n\n- P2976 supplies one concrete new finite typed object after P2975: the nilradical filtration of `Z/12Z`, computed by exhaustive power-orbit witnesses.\n- The nilradical is `{0,6}`, the nonzero nilpotent is `[6]`, `6^2=0`, and all `4` units fix `6`; this is developmental algebraic structure, not a strict source theorem.\n- Do not promote the nilradical filtration, the unit-fixed element `6`, or nilpotent power-orbit witnesses to selector closure, damping source, bridge completion, role transfer, nonproxy `L_total`, or ToE.\n- A next admissible move must attack exactly one missing theorem for this object: strict nadsoliton provenance, coupling to one named missing source atom, or action/variational installation; otherwise preserve the P2929-P2976 developmental/no-strict-export boundary.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
