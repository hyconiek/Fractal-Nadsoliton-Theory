#!/usr/bin/env python3
"""P2972/S1922: source-localizer obstruction for the P2971 incidence complex.

P2971 constructed a typed support/provenance incidence complex and left three
missing theorem routes.  This audit attacks exactly one of them: a strict
nadsoliton source-localizer selecting the incidence object without replaying
ratio arithmetic, K/C exchange, unit conventions, k-selection, or formal Euler
placeholders.

The finite computation enumerates the component/weight automorphism orbits of
the P2971 slots and audits candidate localizers.  Internal predicates can select
slots or label fibers, but they do not export a strict source theorem for the
whole incidence object; provenance and aggregate-coordinate choices remain
bookkeeping labels unless a new source law is supplied.
"""
from __future__ import annotations

import hashlib, itertools, json
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p2971_s1921_typed_support_provenance_incidence_complex import OUT as P2971, SLOTS, automorphism_count

OUT = GEN / "p2972_s1922_incidence_source_localizer_obstruction.json"
MD = GEN / "p2972_s1922_incidence_source_localizer_obstruction.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"


def automorphism_permutations(slots: list[dict[str, Any]], preserve_provenance: bool) -> list[tuple[int, ...]]:
    perms = []
    for perm in itertools.permutations(range(len(slots))):
        ok = True
        for i, j in enumerate(perm):
            a, b = slots[i], slots[j]
            if a["component"] != b["component"] or a["weight"] != b["weight"]:
                ok = False
                break
            if preserve_provenance and a["provenance"] != b["provenance"]:
                ok = False
                break
        if ok:
            perms.append(perm)
    return perms


def orbit_partition(perms: list[tuple[int, ...]], n: int) -> list[list[int]]:
    seen: set[int] = set()
    orbits = []
    for i in range(n):
        if i in seen:
            continue
        orbit = sorted({perm[i] for perm in perms})
        seen.update(orbit)
        orbits.append(orbit)
    return orbits


def invariant_subsets(orbits: list[list[int]]) -> list[list[int]]:
    subsets = []
    for mask in range(1 << len(orbits)):
        chosen: list[int] = []
        for i, orbit in enumerate(orbits):
            if mask & (1 << i):
                chosen.extend(orbit)
        subsets.append(sorted(chosen))
    return subsets


def localizer_rows(orbits: list[list[int]]) -> list[dict[str, Any]]:
    rows = [
        {
            "candidate": "whole_incidence_identity_localizer",
            "selected_slots": [s["slot"] for s in SLOTS],
            "automorphism_invariant": True,
            "selects_whole_incidence_object": True,
            "uses_bookkeeping_label": False,
            "strict_source_theorem_exported": False,
            "unit_coupling_exported": False,
            "accepted_current_strict_localizer": False,
            "witness": "tautologically returns the constructed object, but does not source it",
        },
        {
            "candidate": "unique_weight_one_slot_localizer",
            "selected_slots": ["K0"],
            "automorphism_invariant": True,
            "selects_whole_incidence_object": False,
            "uses_bookkeeping_label": False,
            "strict_source_theorem_exported": False,
            "unit_coupling_exported": False,
            "accepted_current_strict_localizer": False,
            "witness": "selects a slot by weight=1; this is an anchor, not a source for the whole incidence object",
        },
        {
            "candidate": "K_component_slot_pair_localizer",
            "selected_slots": ["K0", "K1"],
            "automorphism_invariant": True,
            "selects_whole_incidence_object": False,
            "uses_bookkeeping_label": True,
            "strict_source_theorem_exported": False,
            "unit_coupling_exported": False,
            "accepted_current_strict_localizer": False,
            "witness": "reuses component label K rather than deriving sourcehood",
        },
        {
            "candidate": "C_equal_summand_orbit_localizer",
            "selected_slots": ["C0", "C1", "C2"],
            "automorphism_invariant": True,
            "selects_whole_incidence_object": False,
            "uses_bookkeeping_label": True,
            "strict_source_theorem_exported": False,
            "unit_coupling_exported": False,
            "accepted_current_strict_localizer": False,
            "witness": "selects the size-3 C orbit; this preserves mismatch but is not a strict source-localizer",
        },
        {
            "candidate": "aggregate_coordinate_order_localizer",
            "selected_slots": ["K0"],
            "automorphism_invariant": False,
            "selects_whole_incidence_object": False,
            "uses_bookkeeping_label": True,
            "strict_source_theorem_exported": False,
            "unit_coupling_exported": False,
            "accepted_current_strict_localizer": False,
            "witness": "uses aggregate coordinate order, a bookkeeping section",
        },
        {
            "candidate": "completed_strict_incidence_source_localizer_schema",
            "selected_slots": [s["slot"] for s in SLOTS],
            "automorphism_invariant": True,
            "selects_whole_incidence_object": True,
            "uses_bookkeeping_label": False,
            "strict_source_theorem_exported": True,
            "unit_coupling_exported": False,
            "accepted_current_strict_localizer": False,
            "witness": "schema only; no current theorem exports the source law",
        },
    ]
    for row in rows:
        row["orbit_compatible"] = sorted(row["selected_slots"]) in [[SLOTS[i]["slot"] for i in subset] for subset in invariant_subsets(orbits)]
    return rows


def obligation_rows(rows: list[dict[str, Any]], orbits: list[list[int]]) -> list[dict[str, Any]]:
    return [
        {"obligation": "slot_orbits_computed", "satisfied": True, "evidence": f"component/weight orbits: {orbits}"},
        {"obligation": "invariant_localizer_candidate_exists", "satisfied": any(r["automorphism_invariant"] for r in rows), "evidence": "whole-object, K0, K-pair, and C-orbit candidates are invariant"},
        {"obligation": "whole_object_candidate_exists", "satisfied": any(r["selects_whole_incidence_object"] for r in rows), "evidence": "identity localizer returns the incidence complex but is tautological"},
        {"obligation": "strict_source_theorem_exported", "satisfied": any(r["strict_source_theorem_exported"] and r["candidate"] != "completed_strict_incidence_source_localizer_schema" for r in rows), "evidence": "no available row sources the incidence complex"},
        {"obligation": "unit_coupling_exported", "satisfied": any(r["unit_coupling_exported"] for r in rows), "evidence": "localizer audit does not install unit-bearing action coupling"},
        {"obligation": "accepted_current_strict_localizer", "satisfied": any(r["accepted_current_strict_localizer"] for r in rows), "evidence": "completed schema remains unavailable"},
    ]


def acceptance_matrix() -> list[dict[str, Any]]:
    names = ["incidence_object", "automorphism_invariant", "whole_object", "not_bookkeeping", "strict_source_theorem", "unit_coupling"]
    full = (1 << len(names)) - 1
    return [{"mask": m, "present": {n: bool(m & (1 << i)) for i, n in enumerate(names)}, "accepts_strict_source_localizer": m == full} for m in range(1 << len(names))]


def build_payload(p2971_path: Any) -> dict[str, Any]:
    perms = automorphism_permutations(SLOTS, preserve_provenance=False)
    orbits = orbit_partition(perms, len(SLOTS))
    rows = localizer_rows(orbits)
    obligations = obligation_rows(rows, orbits)
    matrix = acceptance_matrix()
    return {
        "status": "P2972_INCIDENCE_SOURCE_LOCALIZER_OBSTRUCTION_NO_STRICT_EXPORT",
        "input_hashes": {"P2971": hashlib.sha256(p2971_path.read_bytes()).hexdigest() if p2971_path.exists() else None},
        "constructed_theoretical_objects": {"candidate_object": "IncidenceSourceLocalizer_ObstructionMatrix", "component_weight_automorphism_count": automorphism_count(SLOTS, preserve_provenance=False), "component_weight_orbits": orbits, "invariant_subsets": invariant_subsets(orbits), "localizer_rows": rows, "proof_obligation_rows": obligations, "finite_acceptance_matrix": matrix},
        "localizer_certificate": {"automorphism_count": len(perms), "orbit_count": len(orbits), "invariant_subset_count": len(invariant_subsets(orbits)), "candidate_count": len(rows), "accepted_current_strict_localizers": [r["candidate"] for r in rows if r["accepted_current_strict_localizer"]], "acceptance_matrix_rows": len(matrix), "accepted_rows": sum(1 for r in matrix if r["accepts_strict_source_localizer"])},
        "decision": {
            "positive_progress": "P2972 turns the P2971 missing source-localizer theorem into a finite orbit/localizer audit for the typed incidence complex.",
            "breakthrough": "No strict incidence source-localizer is exported: available candidates are tautological, slot/component anchors, or bookkeeping-label sections, and no unit-bearing coupling theorem accompanies them.",
            "negative_export_flags": {k: False for k in ["strict_source_localizer_exported", "unit_bearing_coupling_exported", "nonproxy_variational_chain_rule_exported", "nonproxy_ltotal_exported", "bridge_closure_exported", "role_transfer_exported", "toe_closure_exported"]},
            "next_honest_step": "Do not replay incidence identity, unique-slot anchors, K/C component labels, aggregate-coordinate order, ratio arithmetic, unit conventions, k-selection predicates, or formal Euler placeholders.  The next proof-grade move must add exactly one of the remaining P2971 theorems: a unit-bearing coupling into a named action density or a nonproxy variational chain rule for the incidence complex; otherwise preserve the P2929-P2972 no-strict-export boundary.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    cert = payload["localizer_certificate"]
    lines = ["# P2972/S1922 incidence source-localizer obstruction", "", f"Status: `{payload['status']}`", "", "## Localizer certificate", f"- automorphism count: `{cert['automorphism_count']}`", f"- orbit count / invariant subsets: `{cert['orbit_count']}` / `{cert['invariant_subset_count']}`", f"- candidate count: `{cert['candidate_count']}`", f"- accepted current strict localizers: `{cert['accepted_current_strict_localizers']}`", f"- acceptance matrix rows/accepted: `{cert['acceptance_matrix_rows']}/{cert['accepted_rows']}`", "", "## Lay summary", payload["decision"]["positive_progress"], payload["decision"]["breakthrough"], "", "## Recommendation", payload["decision"]["next_honest_step"]]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    read_json(P2971)
    payload = build_payload(P2971)
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2972/S1922 incidence source-localizer obstruction", "## P2972/S1922 incidence source-localizer obstruction\n\n`P2972/S1922` attacks one missing theorem left by P2971: a strict source-localizer for the typed support/provenance incidence complex.  The finite audit computes `6` component/weight automorphisms, `3` slot orbits, and `8` invariant subsets, then tests whole-object identity, unique weight-one slot, K component pair, C equal-summand orbit, aggregate-coordinate order, and a completed schema.  Current candidates are tautological, slot/component anchors, or bookkeeping sections; no strict source theorem or unit-bearing coupling is exported.  Thus no strict source-localizer, nonproxy `L_total`, bridge closure, role transfer, or ToE follows.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2972/S1922 incidence-localizer `L_total` guard", "## P2972/S1922 incidence-localizer `L_total` guard\n\n`P2972/S1922` adds no sourced term to `L_total`: invariant source-localizer candidates for the P2971 incidence complex either return the object tautologically, select a slot/component anchor, or use bookkeeping labels.  Since no strict source theorem or unit-bearing action-density coupling is installed, no EOM, Hamiltonian, bridge closure, role transfer, or ToE follows.\n")
    append_once(AGENTS, "Current incidence source-localizer obstruction guardrail (P2972/S1922, 2026-06-20)", "## Current incidence source-localizer obstruction guardrail (P2972/S1922, 2026-06-20)\n\n- P2972 audits the strict source-localizer theorem missing from P2971 for the typed support/provenance incidence complex.\n- The finite orbit/localizer audit finds invariant candidates, but current rows are tautological, slot/component anchors, or bookkeeping-label sections; no strict nadsoliton source theorem or unit-bearing coupling is exported.\n- Do not promote incidence identity, unique-slot anchors, K/C component labels, aggregate-coordinate order, ratio arithmetic, unit conventions, k-selection predicates, or formal Euler placeholders to strict sourcehood, nonproxy `L_total`, bridge closure, role transfer, or ToE.\n- A next admissible move must add exactly one remaining theorem to the incidence object: a unit-bearing coupling into a named action density or a nonproxy variational chain rule; otherwise preserve the P2929-P2972 boundary.\n")
    return payload

if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
