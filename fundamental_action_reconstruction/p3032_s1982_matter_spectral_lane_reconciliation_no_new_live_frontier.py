#!/usr/bin/env python3
"""P3032/S1982: matter spectral lane reconciliation/no-new-live-frontier.

Reconcile the P3029-P3031 matter-spectral lane instead of replaying spectral
magnitude localizers or mass/coupling proxies.  The lane has a real positive
observer-independent observable generator, but no field/sector localizer and no
mass/coupling provenance; therefore it does not export matter physics on current
artifacts.
"""
from __future__ import annotations

import hashlib, itertools, json
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p3029_s1979_matter_spectral_observable_generator_obstruction import OUT as P3029
from p3030_s1980_matter_spectral_field_localizer_obstruction import OUT as P3030
from p3031_s1981_matter_spectral_mass_coupling_provenance_obstruction import OUT as P3031

OUT = GEN / "p3032_s1982_matter_spectral_lane_reconciliation_no_new_live_frontier.json"
MD = GEN / "p3032_s1982_matter_spectral_lane_reconciliation_no_new_live_frontier.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

ATOM_ORDER = [
    "observer_independent_observable_generator",
    "field_representation_localizer",
    "mass_coupling_provenance",
    "selector_sector_unit_action_insertion",
]


def load_payloads() -> tuple[dict[str, Any], dict[str, Any], dict[str, Any]]:
    return read_json(P3029), read_json(P3030), read_json(P3031)


def build_atom_ledger(p3029: dict[str, Any], p3030: dict[str, Any], p3031: dict[str, Any]) -> list[dict[str, Any]]:
    p3029_cert = p3029["finite_certificate"]
    p3030_cert = p3030["finite_certificate"]
    p3031_cert = p3031["finite_certificate"]
    return [
        {
            "atom": "observer_independent_observable_generator",
            "source": "P3029/S1979 sorted DFT magnitude signature",
            "constructed_object_present": True,
            "strict_source_closed": bool(p3029_cert["observer_independent_generator_accepted"]),
            "blocker": "positive pre-physical observable only; not a matter-sector export by itself",
        },
        {
            "atom": "field_representation_localizer",
            "source": "P3030/S1980 phase-retrieval/dihedral orbit matrix",
            "constructed_object_present": True,
            "strict_source_closed": bool(p3030_cert["field_representation_localizer_exported"]),
            "blocker": "24/24 translation/reflection rows preserve magnitude signature; 0/3 localizers accepted",
        },
        {
            "atom": "mass_coupling_provenance",
            "source": "P3031/S1981 rescaling/normalization matrix",
            "constructed_object_present": True,
            "strict_source_closed": bool(p3031_cert["mass_coupling_provenance_exported"]),
            "blocker": "0/4 mass/coupling provenance candidates accepted under K -> 3K test",
        },
        {
            "atom": "selector_sector_unit_action_insertion",
            "source": "P3029/P3030/P3031 negative export flags",
            "constructed_object_present": False,
            "strict_source_closed": False,
            "blocker": "no selector/sector source, no unit-bearing action/EOM insertion, and no Hamiltonian theorem exported",
        },
    ]


def closure_lattice() -> list[dict[str, Any]]:
    rows = []
    for bits in itertools.product([False, True], repeat=len(ATOM_ORDER)):
        profile = dict(zip(ATOM_ORDER, bits))
        rows.append({"profile": profile, "matter_sector_export_licensed": all(bits)})
    return rows


def live_frontier_intake() -> list[dict[str, Any]]:
    return [
        {
            "admissible_input": "new nonpremise field/sector localizer theorem for the P3029 magnitude observable",
            "currently_supplied": False,
            "replay_if_attempted_without_new_object": "P3030 phase-retrieval/orbit-localizer replay",
        },
        {
            "admissible_input": "new target-independent unit-bearing mass/coupling provenance theorem for the P3029 generator",
            "currently_supplied": False,
            "replay_if_attempted_without_new_object": "P3031 spectral proxy/rescaling replay",
        },
        {
            "admissible_input": "different P3028 foundation atom with a genuinely new strict typed object",
            "currently_supplied": False,
            "replay_if_attempted_without_new_object": "spectral magnitude matter-lane replay",
        },
    ]


def build_matrix() -> dict[str, Any]:
    p3029, p3030, p3031 = load_payloads()
    ledger = build_atom_ledger(p3029, p3030, p3031)
    lattice = closure_lattice()
    intake = live_frontier_intake()
    current_profile = {row["atom"]: row["strict_source_closed"] for row in ledger}
    return {
        "object": "MatterSpectralLaneReconciliation_NoNewLiveFrontierCertificate",
        "scope": "P3029-P3031 matter spectral observable lane",
        "atom_ledger": ledger,
        "closure_lattice": lattice,
        "current_profile": current_profile,
        "live_frontier_intake": intake,
        "finite_certificate": {
            "ledger_atoms": len(ledger),
            "closed_atoms": sum(1 for row in ledger if row["strict_source_closed"]),
            "closure_profiles": len(lattice),
            "accepting_profiles": sum(1 for row in lattice if row["matter_sector_export_licensed"]),
            "current_profile_accepts_matter_sector": all(current_profile.values()),
            "new_live_frontier_count": sum(1 for row in intake if row["currently_supplied"]),
        },
    }


def build_payload() -> dict[str, Any]:
    matrix = build_matrix()
    return {
        "status": "P3032_MATTER_SPECTRAL_LANE_RECONCILIATION_NO_NEW_LIVE_FRONTIER",
        "input_hashes": {
            "P3029": hashlib.sha256(P3029.read_bytes()).hexdigest() if P3029.exists() else None,
            "P3030": hashlib.sha256(P3030.read_bytes()).hexdigest() if P3030.exists() else None,
            "P3031": hashlib.sha256(P3031.read_bytes()).hexdigest() if P3031.exists() else None,
        },
        "constructed_theoretical_objects": matrix,
        "finite_certificate": matrix["finite_certificate"],
        "decision": {
            "bounded_no_go": "The P3029-P3031 matter spectral lane has one real positive atom: an observer-independent DFT magnitude observable generator.  The field/localizer, mass/coupling provenance, and selector/sector/unit-action insertion atoms remain open, so the finite closure lattice licenses matter-sector export only in the all-four-atoms-closed profile, not in the current one-closed-atom profile.",
            "negative_export_flags": {k: False for k in ["matter_sector_exported", "field_representation_localizer_exported", "mass_coupling_provenance_exported", "selector_source_exported", "unit_bearing_action_eom_source_exported", "energy_hamiltonian_exported", "observed_physics_exported", "ltotal_exported", "bridge_closure_exported", "role_transfer_exported", "toe_closure_exported"]},
            "next_honest_step": "Preserve the P3029-P3032 matter spectral no-new-live-frontier certificate unless a new nonpremise field/sector localizer theorem, a new target-independent unit-bearing mass/coupling theorem, or a different genuinely new P3028 strict typed object is supplied.  Otherwise pivot to another P3028 foundation atom rather than replaying spectral magnitude matter proxies.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    c = payload["finite_certificate"]
    lines = [
        "# P3032/S1982 matter spectral lane reconciliation", "",
        f"Status: `{payload['status']}`", "", "## Finite certificate",
        f"- ledger atoms: `{c['ledger_atoms']}`",
        f"- closed atoms: `{c['closed_atoms']}`",
        f"- closure profiles: `{c['closure_profiles']}`",
        f"- accepting profiles: `{c['accepting_profiles']}`",
        f"- current profile accepts matter sector: `{c['current_profile_accepts_matter_sector']}`",
        f"- new live frontier count: `{c['new_live_frontier_count']}`", "",
        "## Decision", payload["decision"]["bounded_no_go"], "", "## Recommendation", payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3032/S1982 matter spectral lane reconciliation", "## P3032/S1982 matter spectral lane reconciliation\n\n`P3032/S1982` reconciles the P3029-P3031 matter spectral lane instead of replaying magnitude localizers or mass/coupling proxies.  The finite atom ledger has four atoms: observer-independent observable generator, field-representation/localizer, mass/coupling provenance, and selector/sector/unit-action insertion.  Only the first atom is closed; the `16`-profile closure lattice licenses matter-sector export only in the all-four-atoms-closed profile, while the current profile has `1/4` closed atoms and `0` newly supplied live-frontier inputs.  No matter sector, observed physics, `L_total`, bridge/role transfer, or ToE closure follows.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3032/S1982 matter spectral lane `L_total` guard", "## P3032/S1982 matter spectral lane `L_total` guard\n\n`P3032/S1982` adds no physical `L_total` term.  The matter spectral lane currently closes only the observer-independent observable-generator atom; field/localizer, mass/coupling provenance, selector/sector source, and unit-bearing action/EOM insertion remain unexported.\n")
    append_once(AGENTS, "Current matter spectral lane reconciliation guardrail (P3032/S1982, 2026-06-23)", "## Current matter spectral lane reconciliation guardrail (P3032/S1982, 2026-06-23)\n\n- P3032 reconciles P3029-P3031 instead of replaying spectral magnitude localizers or mass/coupling proxies.\n- The finite ledger has four atoms and only the observer-independent observable-generator atom is closed; the `16`-profile lattice licenses matter-sector export only in the all-four-atoms-closed profile.\n- Do not promote the P3029 spectral magnitude generator, P3030 localizer anchors, or P3031 mass/coupling proxies to matter sector, observed physics, unit-bearing action/EOM, `L_total`, selector, bridge/role-transfer, or ToE closure.\n- A next move must supply a genuinely new nonpremise field/sector localizer theorem, a target-independent unit-bearing mass/coupling theorem, or a different P3028 strict typed object; otherwise preserve this no-new-live-frontier certificate.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
