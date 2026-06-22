#!/usr/bin/env python3
"""P3022/S1972: post-P3017-P3021 unit-atom reconciliation.

This is not a new normalization replay.  It reconciles the typed time-observable
unit sequence after P3017-P3021 and asks whether any unit atom now licenses a
unit-bearing EOM/Hamiltonian, time arrow, observed-physics readout, or ToE move.

The finite ledger has one row for each attacked unit atom: formal EOM/action,
lambda/action normalization, observable-unit readout, clock-unit theorem, and
action-quantum/reference-cell source.  Every row has a constructed finite object,
but every row also has a named obstruction.  The acceptance lattice therefore
has exactly one accepting profile (all atoms closed) and the current profile
closes zero source atoms.
"""
from __future__ import annotations

import hashlib, itertools, json
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p3017_s1967_time_observable_action_eom_source_obstruction import OUT as P3017
from p3018_s1968_lambda_action_unit_normalization_candidate_obstruction import OUT as P3018
from p3019_s1969_observable_unit_readout_source_obstruction import OUT as P3019
from p3020_s1970_clock_unit_theorem_candidate_obstruction import OUT as P3020
from p3021_s1971_action_quantum_reference_cell_source_obstruction import OUT as P3021

OUT = GEN / "p3022_s1972_unit_atom_reconciliation_no_new_live_frontier.json"
MD = GEN / "p3022_s1972_unit_atom_reconciliation_no_new_live_frontier.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

INPUTS = {
    "P3017": P3017,
    "P3018": P3018,
    "P3019": P3019,
    "P3020": P3020,
    "P3021": P3021,
}

UNIT_ATOMS = [
    {
        "atom": "formal_time_observable_action_eom",
        "constructed_object": "cyclic Laplacian formal EOM for T_K",
        "finite_positive": "rank/nullity 11/1 formal EOM scaffold",
        "source_closed": False,
        "blocker": "T_K has nonzero residuals and no strict lambda/action/clock/observable/Hamiltonian units",
    },
    {
        "atom": "lambda_action_unit_normalization",
        "constructed_object": "lambda_* action-to-one normalization",
        "finite_positive": "positive lambda_* normalizes formal action to one",
        "source_closed": False,
        "blocker": "lambda_* rescales as c^-2 under T_K -> c T_K and imports target action quantum",
    },
    {
        "atom": "observable_unit_readout",
        "constructed_object": "RMS/L1/L∞/total-variation observable units",
        "finite_positive": "four positive observer-independent internal units",
        "source_closed": False,
        "blocker": "all units rescale with T_K -> c T_K and no detector/readout unit theorem is exported",
    },
    {
        "atom": "clock_unit_theorem",
        "constructed_object": "label tick, Z12 cycle period, and dominant DFT period",
        "finite_positive": "three dimensionless clock scaffolds",
        "source_closed": False,
        "blocker": "+1 tick is not U(12)-invariant; cycle/DFT periods have no physical frequency/Hamiltonian coupling",
    },
    {
        "atom": "action_quantum_reference_cell_source",
        "constructed_object": "edge/support/rank/DFT-mode reference-cell partitions",
        "finite_positive": "four positive action-per-cell candidates",
        "source_closed": False,
        "blocker": "action-per-cell candidates rescale as c^2 and no independent physical reference cell/action unit is exported",
    },
]


def closure_lattice(atom_names: list[str]) -> dict[str, Any]:
    profiles = []
    for bits in itertools.product([False, True], repeat=len(atom_names)):
        closed = dict(zip(atom_names, bits))
        accepts = all(bits)
        profiles.append({"closed_atoms": closed, "accepts_unit_bearing_eom_hamiltonian": accepts})
    return {
        "profile_count": len(profiles),
        "accepting_profile_count": sum(1 for row in profiles if row["accepts_unit_bearing_eom_hamiltonian"]),
        "current_closed_atom_count": 0,
        "current_profile_accepts": False,
    }


def build_reconciliation() -> dict[str, Any]:
    atom_names = [row["atom"] for row in UNIT_ATOMS]
    lattice = closure_lattice(atom_names)
    live_frontier_candidates = [
        {
            "candidate": "new_strict_time_order_object_with_directed_successor_and_physical_unit_theorem",
            "currently_supplied": False,
            "reason": "no new object beyond the P3017-P3021 internal T_K normalizations is present in this reconciliation",
        }
    ]
    return {
        "object": "TimeObservableUnitAtomReconciliation_NoNewLiveFrontierCertificate",
        "typed_observable": "T_K(d)=K_strict_gate(d+1)-K_strict_gate(d)",
        "unit_atom_rows": UNIT_ATOMS,
        "closure_lattice": lattice,
        "live_frontier_candidates": live_frontier_candidates,
        "new_live_frontier_count": sum(1 for row in live_frontier_candidates if row["currently_supplied"]),
        "accepted_unit_bearing_eom_hamiltonian": lattice["current_profile_accepts"],
    }


def build_payload() -> dict[str, Any]:
    for path in INPUTS.values():
        read_json(path)
    matrix = build_reconciliation()
    return {
        "status": "P3022_UNIT_ATOM_RECONCILIATION_NO_NEW_LIVE_FRONTIER_NO_CLOSURE",
        "input_hashes": {name: hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else None for name, path in INPUTS.items()},
        "constructed_theoretical_objects": matrix,
        "finite_certificate": {
            "unit_atom_count": len(matrix["unit_atom_rows"]),
            "closed_source_atom_count": sum(1 for row in matrix["unit_atom_rows"] if row["source_closed"]),
            "closure_profile_count": matrix["closure_lattice"]["profile_count"],
            "accepting_profile_count": matrix["closure_lattice"]["accepting_profile_count"],
            "current_profile_accepts": matrix["closure_lattice"]["current_profile_accepts"],
            "new_live_frontier_count": matrix["new_live_frontier_count"],
        },
        "decision": {
            "breakthrough": "The P3017-P3021 unit sequence is now reconciled as a finite five-atom ledger.  Each atom has a real constructed object, but none closes as a strict source; the 32-profile closure lattice accepts only the all-atoms-closed profile, while the current profile closes zero atoms and supplies zero new strict time-order candidates.",
            "negative_export_flags": {k: False for k in ["unit_bearing_action_eom_source_exported", "hamiltonian_exported", "time_arrow_exported", "observed_physics_exported", "ltotal_exported", "selector_source_exported", "bridge_closure_exported", "role_transfer_exported", "toe_closure_exported"]},
            "next_honest_step": "Preserve the P3017-P3022 no-new-live-frontier certificate for the T_K unit lane.  The next proof-grade move must introduce a genuinely new strict typed object outside these internal unit normalizations, preferably a strict time-order object with both directed successor and physical unit theorem, before any EOM/Hamiltonian/ToE promotion is reconsidered.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    c = payload["finite_certificate"]
    lines = [
        "# P3022/S1972 unit-atom reconciliation no-new-live-frontier certificate", "",
        f"Status: `{payload['status']}`", "", "## Finite certificate",
        f"- closed source atoms / total: `{c['closed_source_atom_count']}/{c['unit_atom_count']}`",
        f"- accepting closure profiles / total: `{c['accepting_profile_count']}/{c['closure_profile_count']}`",
        f"- current profile accepts: `{c['current_profile_accepts']}`",
        f"- new live frontier count: `{c['new_live_frontier_count']}`", "",
        "## Decision", payload["decision"]["breakthrough"], "", "## Recommendation", payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3022/S1972 unit-atom reconciliation no-new-live-frontier certificate", "## P3022/S1972 unit-atom reconciliation no-new-live-frontier certificate\n\n`P3022/S1972` reconciles the P3017-P3021 typed `T_K` unit lane instead of replaying another internal normalization.  The finite ledger has five attacked atoms: formal action/EOM scaffold, `lambda_*` action normalization, observable-unit readout, clock-unit theorem, and action-quantum/reference-cell source.  Each atom has a constructed finite object, but zero atoms close as strict sources; the Boolean closure lattice has `32` profiles with only the all-atoms-closed profile accepting unit-bearing EOM/Hamiltonian, while the current profile closes zero atoms and supplies no new strict time-order object.  No unit-bearing EOM/Hamiltonian, time arrow, observed-physics export, `L_total`, selector, bridge/role transfer, or ToE closure follows.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3022/S1972 unit-atom reconciliation `L_total` guard", "## P3022/S1972 unit-atom reconciliation `L_total` guard\n\n`P3022/S1972` adds no physical `L_total` term.  It is a no-new-live-frontier certificate for the typed `T_K` unit lane: P3017-P3021 supply finite scaffolds and internal normalizations, but no strict unit source, Hamiltonian normalization, directed time-order theorem, selector source, bridge, or role-transfer theorem.\n")
    append_once(AGENTS, "Current T_K unit-atom reconciliation guardrail (P3022/S1972, 2026-06-22)", "## Current T_K unit-atom reconciliation guardrail (P3022/S1972, 2026-06-22)\n\n- P3022 reconciles the P3017-P3021 typed `T_K` unit lane as a five-atom no-new-live-frontier certificate.\n- Formal action/EOM, `lambda_*`, observable units, clock candidates, and action-quantum/reference-cell partitions all have finite scaffolds, but zero atoms close as strict unit sources; only the all-atoms-closed profile in the `32`-profile lattice would license unit-bearing EOM/Hamiltonian.\n- Do not replay internal `T_K` unit normalizations as unit-bearing EOM/Hamiltonian, time arrow, observed-physics, `L_total`, selector, bridge/role-transfer, or ToE closure.\n- The next honest move must introduce a genuinely new strict typed object outside these internal unit normalizations, preferably a strict time-order object carrying both a directed successor and physical unit theorem, before any EOM/Hamiltonian/ToE promotion is reconsidered.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
