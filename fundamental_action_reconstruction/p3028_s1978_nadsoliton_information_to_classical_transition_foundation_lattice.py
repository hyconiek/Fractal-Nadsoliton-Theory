#!/usr/bin/env python3
"""P3028/S1978: nadsoliton-information to classical-transition foundation lattice.

Pivot outside the closed dissipation time-order lane.  Construct the missing
foundation lattice for promoting primordial nadsoliton information to classical
physics readout.  This is an obligation lattice, not a closure claim.
"""
from __future__ import annotations

import hashlib, itertools, json
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p3013_s1963_observer_physics_readout_strict_kernel_selector_obligation_matrix import OUT as P3013
from p3027_s1977_dissipation_external_unit_source_acceptance_gate import OUT as P3027

OUT = GEN / "p3028_s1978_nadsoliton_information_to_classical_transition_foundation_lattice.json"
MD = GEN / "p3028_s1978_nadsoliton_information_to_classical_transition_foundation_lattice.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

FOUNDATION_ATOMS = [
    "strict_selector_or_branch_source",
    "external_physical_unit_source",
    "unit_bearing_action_eom_hamiltonian",
    "observer_independent_observable_generator",
    "classical_coarse_graining_limit",
]

READOUT_ROWS = ["spacetime_geometry", "time", "matter_fields", "energy_hamiltonian", "observer_readout"]


def row_matrix() -> list[dict[str, Any]]:
    blockers = {
        "spacetime_geometry": ["metric/causal chart source", "unit-bearing geometric action", "nonproxy residual closure"],
        "time": ["nonpremise arrow selector", "physical clock unit", "Hamiltonian evolution theorem"],
        "matter_fields": ["field representation/localizer", "mass/coupling provenance", "stable excitation theorem"],
        "energy_hamiltonian": ["action quantum", "physical tick", "stress-energy/Hamiltonian coupling"],
        "observer_readout": ["observer-independent observable generator", "downstream coarse-graining", "actual observer closure without retrocausal sourcing"],
    }
    rows = []
    for row in READOUT_ROWS:
        obligations = {atom: False for atom in FOUNDATION_ATOMS}
        rows.append({
            "classical_row": row,
            "ontology_order": "nadsoliton -> light -> matter -> emergent observer",
            "foundation_atoms": obligations,
            "missing_blockers": blockers[row],
            "accepted_as_classical_export": all(obligations.values()),
        })
    return rows


def closure_lattice() -> list[dict[str, Any]]:
    rows = []
    for bits in itertools.product([False, True], repeat=len(FOUNDATION_ATOMS)):
        profile = dict(zip(FOUNDATION_ATOMS, bits))
        rows.append({
            "profile": profile,
            "closed_atom_count": sum(bits),
            "accepts_classical_transition": all(bits),
        })
    return rows


def build_payload(paths: dict[str, Any]) -> dict[str, Any]:
    payloads = {name: read_json(path) for name, path in paths.items()}
    matrix = row_matrix()
    lattice = closure_lattice()
    return {
        "status": "P3028_NADSOLITON_INFORMATION_TO_CLASSICAL_TRANSITION_FOUNDATION_LATTICE_NO_CLOSURE",
        "input_hashes": {name: hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else None for name, path in paths.items()},
        "constructed_theoretical_objects": {
            "object": "NadsolitonInformationToClassicalTransition_FoundationObligationLattice",
            "input_statuses": {name: payload.get("status") for name, payload in payloads.items()},
            "foundation_atoms": FOUNDATION_ATOMS,
            "classical_readout_rows": matrix,
            "closure_lattice": lattice,
        },
        "finite_certificate": {
            "foundation_atom_count": len(FOUNDATION_ATOMS),
            "classical_readout_row_count": len(matrix),
            "accepted_readout_rows": sum(1 for row in matrix if row["accepted_as_classical_export"]),
            "closure_profile_count": len(lattice),
            "accepting_profile_count": sum(1 for row in lattice if row["accepts_classical_transition"]),
            "current_closed_foundation_atoms": 0,
        },
        "decision": {
            "breakthrough": "The information-to-classical transition is made explicit as a five-atom foundation lattice over spacetime, time, matter, energy/Hamiltonian, and observer-readout rows.  Current artifacts supply important scaffolds, but zero readout rows have all required foundations; only the all-five-atoms-closed profile in the 32-profile lattice would license a classical-transition export.",
            "negative_export_flags": {k: False for k in ["classical_transition_exported", "observed_spacetime_exported", "time_arrow_exported", "matter_sector_exported", "energy_hamiltonian_exported", "observer_readout_exported", "ltotal_exported", "selector_source_exported", "bridge_closure_exported", "role_transfer_exported", "toe_closure_exported"]},
            "next_honest_step": "Attack exactly one P3028 foundation atom with a new strict source object.  The sharpest non-replay target is an observer-independent observable generator for one classical row with explicit domain/codomain and a finite acceptance test; do not promote the full classical transition until all five atoms close.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    c = payload["finite_certificate"]
    lines = [
        "# P3028/S1978 nadsoliton information-to-classical transition foundation lattice", "",
        f"Status: `{payload['status']}`", "", "## Finite certificate",
        f"- foundation atoms: `{c['foundation_atom_count']}`",
        f"- classical readout rows: `{c['classical_readout_row_count']}`",
        f"- accepted readout rows: `{c['accepted_readout_rows']}`",
        f"- closure profiles: `{c['closure_profile_count']}`",
        f"- accepting profiles: `{c['accepting_profile_count']}`",
        f"- current closed foundation atoms: `{c['current_closed_foundation_atoms']}`", "",
        "## Decision", payload["decision"]["breakthrough"], "", "## Recommendation", payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload({"P3013": P3013, "P3027": P3027})
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3028/S1978 information-to-classical transition foundation lattice", "## P3028/S1978 information-to-classical transition foundation lattice\n\n`P3028/S1978` pivots outside the closed dissipation time-order lane and constructs the missing foundation lattice for promoting primordial nadsoliton information to classical physics readout.  The lattice has five required atoms: strict selector/branch source, external physical unit source, unit-bearing action/EOM/Hamiltonian, observer-independent observable generator, and classical coarse-graining limit.  Across spacetime, time, matter, energy/Hamiltonian, and observer-readout rows, zero rows currently close; the Boolean lattice has `32` profiles with only the all-five-atoms-closed profile accepting a classical-transition export.  No observed spacetime, time arrow, matter sector, energy/Hamiltonian, observer readout, `L_total`, bridge/role transfer, or ToE closure follows.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3028/S1978 information-to-classical transition `L_total` guard", "## P3028/S1978 information-to-classical transition `L_total` guard\n\n`P3028/S1978` adds no physical `L_total` term.  It is a foundation-obligation lattice for the nadsoliton-information to classical-readout transition; zero rows close because selector, external unit source, unit-bearing action/EOM/Hamiltonian, observable generator, and coarse-graining atoms are not all supplied.\n")
    append_once(AGENTS, "Current information-to-classical transition foundation guardrail (P3028/S1978, 2026-06-22)", "## Current information-to-classical transition foundation guardrail (P3028/S1978, 2026-06-22)\n\n- P3028 constructs a five-atom foundation lattice for promoting primordial nadsoliton information to classical physics readout.\n- Required atoms are strict selector/branch source, external physical unit source, unit-bearing action/EOM/Hamiltonian, observer-independent observable generator, and classical coarse-graining limit.\n- Current artifacts have zero accepted classical readout rows; only the all-five-atoms-closed profile in the `32`-profile lattice would license classical-transition export.\n- Do not promote this lattice to observed spacetime, time arrow, matter sector, energy/Hamiltonian, observer readout, `L_total`, bridge/role-transfer, or ToE closure; the next move must attack exactly one foundation atom with a new strict source object.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
