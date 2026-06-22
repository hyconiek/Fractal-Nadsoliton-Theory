#!/usr/bin/env python3
"""P3013/S1963: observer-physics readout of strict kernel vs selector.

This is a bounded answer to the user's physical-perspective question: if the
program is a ToE candidate and we see universe/time/matter/energy, what does the
strict kernel describe and what would the selector be in that perspective?

The constructed object is an observer-physics readout obligation matrix.  It
keeps the ontology/order intact: nadsoliton first, then light/matter/observer as
emergent readouts.  In this matrix K_strict_gate is the pre-physical/strict
working correlation-compression law; a selector is the missing physical-sector
choice/source that turns symmetric/projective carrier data into one directed
branch/readout sector.  The matrix checks what extra atoms are needed before this
becomes our observed spacetime, time, matter, energy, and observer physics.
"""
from __future__ import annotations

import hashlib, json
from itertools import product
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p3012_s1962_strict_kernel_phase_gradient_selector_source_obstruction import OUT as P3012
from p3011_s1961_toe_strict_kernel_selector_role_separation_matrix import OUT as P3011

OUT = GEN / "p3013_s1963_observer_physics_readout_strict_kernel_selector_obligation_matrix.json"
MD = GEN / "p3013_s1963_observer_physics_readout_strict_kernel_selector_obligation_matrix.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

PHYSICS_ROWS = [
    {
        "row": "spacetime_geometry",
        "what_we_see": "3+1 dimensional spacetime/geometry and causal locality",
        "kernel_readout": "strict separation/correlation/compression profile over d, not a metric manifold by itself",
        "selector_readout": "choice of oriented/causal/geometric sector or chart when quotient data leave equivalent branches",
        "missing_atoms": ["canonical metric/causal chart", "unit-bearing geometric action", "nonproxy GR residual closure"],
    },
    {
        "row": "time",
        "what_we_see": "directed time, clocks, before/after ordering, thermodynamic arrow",
        "kernel_readout": "ordered damping/compression profile can parameterize scale or separation but does not select a time arrow",
        "selector_readout": "orientation/arrow source selecting one directed temporal branch rather than the reversed/projective partner",
        "missing_atoms": ["nonpremise arrow source", "clock-unit map", "causal evolution theorem"],
    },
    {
        "row": "matter",
        "what_we_see": "particles/fields/masses/flavour sectors and stable excitations",
        "kernel_readout": "correlation spectrum/shape candidate constraining sectors, not field representations or masses alone",
        "selector_readout": "sector/localizer choice assigning carrier branches to physical field slots",
        "missing_atoms": ["field representation map", "mass/coupling provenance", "sector-localizer theorem"],
    },
    {
        "row": "energy",
        "what_we_see": "energy, Hamiltonian, stress-energy and conserved densities",
        "kernel_readout": "compression/tail profile can feed formal coefficient candidates but is not a unit-bearing energy density",
        "selector_readout": "branch/sign/source choice for positive densities and Hamiltonian sector",
        "missing_atoms": ["unit-bearing action density", "sign convention theorem", "Hamiltonian/stress-energy export"],
    },
    {
        "row": "observer_readout",
        "what_we_see": "emergent observer measuring a single experienced physical branch",
        "kernel_readout": "pre-observer nadsoliton correlation law constraining possible readouts",
        "selector_readout": "physical-sector/readout selection that makes one branch available to emergent observers without putting observer underneath nadsoliton",
        "missing_atoms": ["observer-independent observable map", "measurement/readout functor", "no observer-as-foundation premise"],
    },
]

REQUIRED_ATOMS = ["kernel_profile", "selector_source", "unit_action", "eom_hamiltonian", "observable_generator"]


def row_obligation(row: dict[str, Any]) -> dict[str, Any]:
    current = {
        "kernel_profile": True,
        "selector_source": False,
        "unit_action": False,
        "eom_hamiltonian": False,
        "observable_generator": False,
    }
    return {
        **row,
        "current_atoms": current,
        "accepted_as_observed_physics_export": all(current.values()),
        "blocked_atoms": [k for k, v in current.items() if not v],
        "minimal_acceptance_profiles": [dict(zip(REQUIRED_ATOMS, bits)) for bits in product([False, True], repeat=len(REQUIRED_ATOMS)) if all(bits)],
    }


def readout_matrix() -> dict[str, Any]:
    rows = [row_obligation(r) for r in PHYSICS_ROWS]
    all_profiles = [dict(zip(REQUIRED_ATOMS, bits)) for bits in product([False, True], repeat=len(REQUIRED_ATOMS))]
    return {
        "object": "ObserverPhysicsReadout_StrictKernelSelector_ObligationMatrix",
        "perspective_answer": {
            "strict_kernel_in_our_physics": "the pre-physical strict working law of nadsoliton correlations/compression that must be read out through units, fields, EOM, and observables before it becomes spacetime/matter/energy physics",
            "selector_in_our_physics": "the missing physical-sector/branch/orientation source that chooses one directed readout sector from symmetric/projective carrier data, not an observer underneath the nadsoliton",
            "toe_requirement": "a ToE must add readout maps from the strict kernel plus selector into spacetime, time, matter, energy, and observer observables; current artifacts do not yet close those maps",
        },
        "rows": rows,
        "row_count": len(rows),
        "accepted_row_count": sum(1 for r in rows if r["accepted_as_observed_physics_export"]),
        "global_acceptance_profiles": all_profiles,
        "global_acceptance_profile_count": len(all_profiles),
        "global_accepting_profile_count": sum(1 for p in all_profiles if all(p.values())),
    }


def build_payload(p3011_path: Any, p3012_path: Any) -> dict[str, Any]:
    matrix = readout_matrix()
    return {
        "status": "P3013_OBSERVER_PHYSICS_READOUT_STRICT_KERNEL_SELECTOR_OBLIGATION_MATRIX_NO_CLOSURE",
        "input_hashes": {
            "P3011": hashlib.sha256(p3011_path.read_bytes()).hexdigest() if p3011_path.exists() else None,
            "P3012": hashlib.sha256(p3012_path.read_bytes()).hexdigest() if p3012_path.exists() else None,
        },
        "constructed_theoretical_objects": matrix,
        "readout_certificate": {
            "physics_rows": [r["row"] for r in matrix["rows"]],
            "row_count": matrix["row_count"],
            "accepted_row_count": matrix["accepted_row_count"],
            "global_acceptance_profile_count": matrix["global_acceptance_profile_count"],
            "global_accepting_profile_count": matrix["global_accepting_profile_count"],
            "blocked_atom_union": sorted({a for r in matrix["rows"] for a in r["blocked_atoms"]}),
        },
        "decision": {
            "plain_answer": "In our physics perspective, the strict kernel would describe the pre-physical nadsoliton correlation/compression law whose readouts become spacetime, time, matter, and energy only after unit-bearing action, EOM/Hamiltonian, observable, and selector/source maps are supplied.  The selector would be the physical-sector/branch/orientation source that selects one directed experienced branch, not the kernel itself and not an observer under the nadsoliton.",
            "breakthrough": "Bounded progress: five observer-physics rows are made explicit, and every row is blocked by the same missing selector/unit-action/EOM/observable atoms; zero rows are accepted as observed-physics export on current artifacts.",
            "negative_export_flags": {k: False for k in ["observed_spacetime_exported", "time_arrow_exported", "matter_sector_exported", "energy_hamiltonian_exported", "observer_readout_exported", "strict_selector_source_exported", "unit_bearing_ltotal_exported", "toe_closure_exported"]},
            "next_honest_step": "Attack exactly one readout atom, preferably an observer-independent observable generator for one row (time, matter, energy, spacetime, or observer_readout) with explicit input/output types and a finite acceptance test; do not claim ToE until selector source, unit action, EOM/Hamiltonian, and observable generator are all present.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    cert = payload["readout_certificate"]
    lines = [
        "# P3013/S1963 observer-physics readout strict-kernel/selector obligation matrix", "",
        f"Status: `{payload['status']}`", "", "## Direct answer",
        payload["decision"]["plain_answer"], "", "## Readout certificate",
        f"- rows: `{cert['physics_rows']}`",
        f"- row count / accepted: `{cert['row_count']}/{cert['accepted_row_count']}`",
        f"- global acceptance profiles / accepting: `{cert['global_acceptance_profile_count']}/{cert['global_accepting_profile_count']}`",
        f"- blocked atom union: `{cert['blocked_atom_union']}`", "", "## Decision", payload["decision"]["breakthrough"], "", "## Recommendation", payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    read_json(P3011)
    read_json(P3012)
    payload = build_payload(P3011, P3012)
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3013/S1963 observer-physics readout strict-kernel/selector obligation matrix", "## P3013/S1963 observer-physics readout strict-kernel/selector obligation matrix\n\n`P3013/S1963` answers the observer-physics version of the strict-kernel/selector question.  In the current honest readout, `K_strict_gate` describes the pre-physical strict working law of nadsoliton correlations/compression; it is not yet spacetime, time, matter, energy, or observer physics by itself.  The selector would be the physical-sector/branch/orientation source that chooses one directed readout sector from symmetric/projective carrier data, not the kernel itself and not an observer placed underneath the nadsoliton.  The five-row matrix for spacetime geometry, time, matter, energy, and observer readout has zero accepted rows because selector source, unit-bearing action, EOM/Hamiltonian, and observable-generator atoms remain missing.  No observed-spacetime export, time-arrow export, matter-sector export, energy/Hamiltonian export, observer-readout export, `L_total`, role transfer, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3013/S1963 observer-physics readout `L_total` guard", "## P3013/S1963 observer-physics readout `L_total` guard\n\n`P3013/S1963` adds no `L_total` term.  It is a readout-obligation matrix translating the strict kernel/selector distinction into our-physics rows: spacetime, time, matter, energy, and observer readout.  Because selector source, unit-bearing action, EOM/Hamiltonian, and observable-generator atoms remain missing, no physical sector, stress-energy/Hamiltonian, observer-readout functor, bridge/role transfer, or ToE closure may be installed.\n")
    append_once(AGENTS, "Current observer-physics readout strict-kernel/selector guardrail (P3013/S1963, 2026-06-22)", "## Current observer-physics readout strict-kernel/selector guardrail (P3013/S1963, 2026-06-22)\n\n- P3013 translates the strict-kernel/selector distinction into our-physics rows: spacetime geometry, time, matter, energy, and observer readout.\n- In this perspective `K_strict_gate` is the pre-physical strict working law of nadsoliton correlations/compression; the selector would be a physical-sector/branch/orientation source choosing one directed readout sector from symmetric/projective carrier data.\n- The five-row readout matrix has zero accepted observed-physics exports because selector source, unit-bearing action, EOM/Hamiltonian, and observable-generator atoms remain missing.\n- Do not promote P3013 to observed spacetime, time arrow, matter sector, energy/Hamiltonian, observer-readout, `L_total`, bridge/role transfer, or ToE closure.  The next proof-grade move should attack exactly one readout atom, preferably an observer-independent observable generator for one row with explicit input/output types and a finite acceptance test.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
