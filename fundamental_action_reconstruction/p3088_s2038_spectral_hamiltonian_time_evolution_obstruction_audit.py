#!/usr/bin/env python3
"""P3088/S2038: spectral-to-Hamiltonian/time-evolution obstruction audit.

P3087 left the Z12 Dirichlet/Laplacian branch with finite dimensionless
thermodynamic algebra but no non-imported thermodynamic source.  P3088 attacks
exactly one adjacent standard-physics interface atom: whether the internal
finite Z12 spectrum sources a physical Hamiltonian/time-evolution sector with a
self-adjoint Hamiltonian, sourced time parameter, unitary evolution law, energy
units, and observable dynamics without importing quantum mechanics, measurement
units, spacetime EOM, observed radiation/light, selector closure, L_total,
bridge/role-transfer, or ToE.
"""
from __future__ import annotations

import hashlib, json, math, subprocess
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p3087_s2037_thermodynamic_statistical_ensemble_obstruction_audit import OUT as P3087

OUT = GEN / "p3088_s2038_spectral_hamiltonian_time_evolution_obstruction_audit.json"
MD = GEN / "p3088_s2038_spectral_hamiltonian_time_evolution_obstruction_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"
Z12_SIZE = 12
TIME_GRID = (0.0, 0.125, 0.25, 0.5, 1.0, 2.0)
TIME_SCALES = (0.5, 1.0, 2.0)

CONTENT_PATTERNS = {
    "hamiltonian_atom": r"Hamiltonian|time-evolution|unitary|self-adjoint|observable dynamics",
    "predecessor_thermo_atom": r"thermodynamic|statistical-ensemble|partition function|temperature|Boltzmann",
    "blocked_promotions": r"observed photons|observed light|observed radiation|spacetime EOM|L_total|ToE|selector closure|bridge/role-transfer",
}

HAMILTONIAN_CANDIDATES = (
    {
        "id": "z12_laplacian_spectral_multiplier",
        "description": "the real symmetric Z12 Laplacian diagonalized as a formal spectral multiplier",
        "internal_spectrum_exported": True,
        "self_adjoint_operator_exported": True,
        "sourced_time_parameter_exported": False,
        "unitary_evolution_law_exported": False,
        "energy_unit_source_exported": False,
        "observable_dynamics_exported": False,
        "blocker": "self-adjoint finite linear algebra is present, but time, units, and physical dynamics are not sourced",
    },
    {
        "id": "formal_schrodinger_lift_exp_minus_i_lambda_t",
        "description": "formal modal phases exp(-i lambda_m t) built from the dimensionless spectrum",
        "internal_spectrum_exported": True,
        "self_adjoint_operator_exported": True,
        "sourced_time_parameter_exported": False,
        "unitary_evolution_law_exported": True,
        "energy_unit_source_exported": False,
        "observable_dynamics_exported": False,
        "blocker": "unitarity is only a formal imported Schrödinger-style template with dimensionless t",
    },
    {
        "id": "p3087_dimensionless_partition_generator",
        "description": "reuse of P3087 dimensionless energy labels as a would-be generator",
        "internal_spectrum_exported": True,
        "self_adjoint_operator_exported": True,
        "sourced_time_parameter_exported": False,
        "unitary_evolution_law_exported": False,
        "energy_unit_source_exported": False,
        "observable_dynamics_exported": False,
        "blocker": "thermal beta algebra does not source Hamiltonian time evolution or energy units",
    },
    {
        "id": "imported_quantum_mechanics_template",
        "description": "external Hilbert-space quantum mechanics template with hbar and Schrödinger time",
        "internal_spectrum_exported": False,
        "self_adjoint_operator_exported": True,
        "sourced_time_parameter_exported": True,
        "unitary_evolution_law_exported": True,
        "energy_unit_source_exported": True,
        "observable_dynamics_exported": True,
        "blocker": "passes only by importing quantum mechanics and hbar/time/energy units",
    },
    {
        "id": "imported_spacetime_field_eom_template",
        "description": "external spacetime field-equation template with Hamiltonian density and observables",
        "internal_spectrum_exported": False,
        "self_adjoint_operator_exported": True,
        "sourced_time_parameter_exported": True,
        "unitary_evolution_law_exported": False,
        "energy_unit_source_exported": True,
        "observable_dynamics_exported": True,
        "blocker": "passes only as imported spacetime/EOM physics, not as a strict Z12 source",
    },
)
REQUIRED_GATES = ("internal_spectrum_exported", "self_adjoint_operator_exported", "sourced_time_parameter_exported", "unitary_evolution_law_exported", "energy_unit_source_exported", "observable_dynamics_exported")


def content_grep() -> list[dict[str, Any]]:
    rows = []
    for lane, pattern in CONTENT_PATTERNS.items():
        proc = subprocess.run(["rg", "-n", pattern, "fundamental_action_reconstruction", "AGENTS.md", "-S"], cwd=REPO, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
        hits = [line for line in proc.stdout.splitlines() if line.strip()]
        rows.append({"content_lane": lane, "pattern": pattern, "hit_count": len(hits), "sample_hits": hits[:8]})
    return rows


def z12_spectrum() -> list[float]:
    return [2.0 - 2.0 * math.cos(2.0 * math.pi * m / Z12_SIZE) for m in range(Z12_SIZE)]


def spectral_hamiltonian_rows(spectrum: list[float]) -> list[dict[str, Any]]:
    return [{
        "mode": m,
        "lambda_m": round(lam, 12),
        "formal_hamiltonian_eigenvalue": round(lam, 12),
        "finite_self_adjoint_real_symmetric_witness": True,
        "energy_unit_attached": False,
        "physical_hamiltonian_exported": False,
    } for m, lam in enumerate(spectrum)]


def formal_unitary_phase_rows(spectrum: list[float]) -> list[dict[str, Any]]:
    rows = []
    for t in TIME_GRID:
        max_modulus_error = max(abs(abs(complex(math.cos(-lam * t), math.sin(-lam * t))) - 1.0) for lam in spectrum)
        rows.append({
            "dimensionless_time_parameter": t,
            "modal_phase_count": len(spectrum),
            "max_unitarity_modulus_error": round(max_modulus_error, 15),
            "formal_unitary_identity_holds": max_modulus_error <= 1e-12,
            "time_unit_attached": False,
            "hbar_or_action_unit_attached": False,
            "observable_dynamics_exported": False,
        })
    return rows


def time_energy_scale_orbit_rows(spectrum: list[float]) -> list[dict[str, Any]]:
    rows = []
    for scale in TIME_SCALES:
        for t in TIME_GRID:
            max_phase_error = max(abs((scale * lam) * t - lam * (scale * t)) for lam in spectrum)
            rows.append({
                "energy_scale": scale,
                "dimensionless_time_parameter": t,
                "scaled_energy_phase_matches_rescaled_time_phase": max_phase_error <= 1e-12,
                "max_phase_error": round(max_phase_error, 15),
                "canonical_time_or_energy_scale_selected": scale == 1.0,
                "selection_is_conventional_without_unit_source": True,
            })
    return rows


def gate_rows() -> list[dict[str, Any]]:
    return [{"candidate": c["id"], "required_gate": gate, "gate_passed": bool(c[gate]), "blocking_if_failed": not bool(c[gate]), "detail": "passed" if c[gate] else c["blocker"]} for c in HAMILTONIAN_CANDIDATES for gate in REQUIRED_GATES]


def aggregate(gates: list[dict[str, Any]]) -> list[dict[str, Any]]:
    out = []
    for c in HAMILTONIAN_CANDIDATES:
        subset = [r for r in gates if r["candidate"] == c["id"]]
        out.append({"candidate": c["id"], "passed_gates": sum(r["gate_passed"] for r in subset), "failed_gates": sum(not r["gate_passed"] for r in subset), "accepted_nonimported_hamiltonian_time_evolution_source": all(r["gate_passed"] for r in subset) and bool(c["internal_spectrum_exported"])})
    return out


def build_payload() -> dict[str, Any]:
    p3087 = read_json(P3087)
    greps = content_grep(); spectrum = z12_spectrum(); hrows = spectral_hamiltonian_rows(spectrum); urows = formal_unitary_phase_rows(spectrum); scale_rows = time_energy_scale_orbit_rows(spectrum); gates = gate_rows(); aggs = aggregate(gates)
    accepted = [r for r in aggs if r["accepted_nonimported_hamiltonian_time_evolution_source"]]
    obligations = [
        {"obligation": "read_p3087_next_atom", "satisfied": True, "detail": "P3087 selected spectral-to-Hamiltonian/time-evolution audit as the next interface atom"},
        {"obligation": "construct_self_adjoint_spectral_rows", "satisfied": all(r["finite_self_adjoint_real_symmetric_witness"] for r in hrows), "detail": "the finite Z12 Laplacian is represented as a real symmetric/self-adjoint spectral multiplier"},
        {"obligation": "compute_formal_unitary_phase_grid", "satisfied": all(r["formal_unitary_identity_holds"] for r in urows), "detail": "formal exp(-i lambda t) phases have unit modulus on the dimensionless time grid"},
        {"obligation": "test_time_energy_scale_orbit", "satisfied": all(r["scaled_energy_phase_matches_rescaled_time_phase"] for r in scale_rows), "detail": "energy scaling can be absorbed into dimensionless time scaling, so no canonical time/energy unit is selected"},
        {"obligation": "export_nonimported_observable_dynamics", "satisfied": False, "detail": "0 candidates pass as internal non-imported Hamiltonian/time-evolution observable sources"},
    ]
    return {
        "status": "P3088_SPECTRAL_HAMILTONIAN_TIME_EVOLUTION_OBSTRUCTION_BOUNDED_NO_GO",
        "input_hashes": {"P3087": hashlib.sha256(P3087.read_bytes()).hexdigest() if P3087.exists() else None},
        "constructed_theoretical_objects": {
            "content_first_repo_grep": greps,
            "hamiltonian_time_evolution_audit_object": {"object": "Z12DirichletSpectralHamiltonianTimeEvolutionObstructionAudit", "source_reused": "P3087 recommendation: bounded spectral-to-Hamiltonian/time-evolution obstruction audit", "required_gates": list(REQUIRED_GATES), "candidate_hamiltonian_sources": [c["id"] for c in HAMILTONIAN_CANDIDATES], "acceptance_predicate": "internal spectrum plus self-adjoint operator, sourced time parameter, unitary evolution law, energy/action unit source, and observable dynamics export"},
            "z12_laplacian_spectrum": [round(v, 12) for v in spectrum],
            "spectral_hamiltonian_rows": hrows,
            "formal_unitary_phase_rows": urows,
            "time_energy_scale_orbit_rows": scale_rows,
            "candidate_gate_rows": gates,
            "candidate_aggregate_certificate": aggs,
        },
        "finite_certificate": {
            "content_grep_lanes": len(greps), "content_grep_hits": sum(r["hit_count"] for r in greps),
            "p3087_accepted_nonimported_thermodynamic_ensemble_sources": p3087.get("finite_certificate", {}).get("accepted_nonimported_thermodynamic_ensemble_sources"),
            "spectrum_rows": len(spectrum), "spectral_hamiltonian_rows": len(hrows), "spectral_rows_with_energy_units": sum(r["energy_unit_attached"] for r in hrows),
            "formal_unitary_phase_rows": len(urows), "unitarity_modulus_failures": sum(not r["formal_unitary_identity_holds"] for r in urows), "phase_rows_with_time_units": sum(r["time_unit_attached"] for r in urows), "phase_rows_with_action_units": sum(r["hbar_or_action_unit_attached"] for r in urows),
            "time_energy_scale_orbit_rows": len(scale_rows), "time_energy_scale_identity_failures": sum(not r["scaled_energy_phase_matches_rescaled_time_phase"] for r in scale_rows), "canonical_time_or_energy_sources_exported": 0,
            "hamiltonian_candidates": len(HAMILTONIAN_CANDIDATES), "required_gates": len(REQUIRED_GATES), "candidate_gate_rows": len(gates), "accepted_nonimported_hamiltonian_time_evolution_sources": len(accepted),
            "proof_obligations": len(obligations), "satisfied_proof_obligations": sum(r["satisfied"] for r in obligations),
        },
        "proof_obligations": obligations,
        "decision": {
            "bounded_result": "P3088 constructs the requested spectral-to-Hamiltonian/time-evolution obstruction audit.  The finite Z12 Laplacian supplies a real symmetric/self-adjoint spectral multiplier and formal exp(-i lambda t) phase rows are exactly unit-modulus on a dimensionless grid.  However, the time parameter, hbar/action normalization, energy units, and observable-dynamics readout remain unsourced; energy scaling can be absorbed by rescaling dimensionless time.  Imported quantum-mechanics and spacetime-EOM templates pass only as imported templates.  Therefore no non-imported physical Hamiltonian/time-evolution observable source is exported.",
            "negative_export_flags": {key: False for key in ["sourced_time_parameter_exported", "energy_unit_source_exported", "hbar_or_action_unit_source_exported", "observable_dynamics_exported", "physical_hamiltonian_exported", "measurement_unit_source_exported", "empirical_observable_exported", "observed_radiation_exported", "observed_photons_exported", "observed_light_exported", "spacetime_eom_exported", "thermodynamic_source_exported", "QW_2191_discharged", "strict_selector_closure_exported", "ltotal_exported", "bridge_role_transfer_exported", "toe_closure_exported"]},
            "positive_scoped_flags": {"finite_self_adjoint_spectral_multiplier_computed": True, "formal_unitary_phase_grid_computed": True, "time_energy_scale_orbit_verified": True},
            "next_honest_step": "Pivot to exactly one adjacent standard-physics interface atom outside selector replay: construct a bounded spectral-observable/Born-rule probability-readout obstruction audit for the Z12 Dirichlet/Laplacian branch, testing whether internal eigenmodes and formal unitary phases supply a Hilbert-state normalization, positive probability measure, Born-rule map, measurement basis/source, and empirical probability readout without importing quantum measurement theory, apparatus units, observed radiation/light, spacetime EOM, L_total, bridge/role-transfer, or ToE.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    c = payload["finite_certificate"]
    lines = ["# P3088/S2038 spectral-to-Hamiltonian/time-evolution obstruction audit", "", f"Status: `{payload['status']}`", "", "## Finite certificate",
        f"- content grep lanes: `{c['content_grep_lanes']}`", f"- content grep hits: `{c['content_grep_hits']}`", f"- P3087 accepted non-imported thermodynamic ensemble sources: `{c['p3087_accepted_nonimported_thermodynamic_ensemble_sources']}`", f"- spectrum rows: `{c['spectrum_rows']}`", f"- spectral Hamiltonian rows: `{c['spectral_hamiltonian_rows']}`", f"- spectral rows with energy units: `{c['spectral_rows_with_energy_units']}`", f"- formal unitary phase rows: `{c['formal_unitary_phase_rows']}`", f"- unitarity modulus failures: `{c['unitarity_modulus_failures']}`", f"- phase rows with time units: `{c['phase_rows_with_time_units']}`", f"- phase rows with action units: `{c['phase_rows_with_action_units']}`", f"- time-energy scale orbit rows: `{c['time_energy_scale_orbit_rows']}`", f"- time-energy scale identity failures: `{c['time_energy_scale_identity_failures']}`", f"- canonical time or energy sources exported: `{c['canonical_time_or_energy_sources_exported']}`", f"- Hamiltonian candidates: `{c['hamiltonian_candidates']}`", f"- candidate gate rows: `{c['candidate_gate_rows']}`", f"- accepted non-imported Hamiltonian/time-evolution sources: `{c['accepted_nonimported_hamiltonian_time_evolution_sources']}`", f"- satisfied proof obligations: `{c['satisfied_proof_obligations']}/{c['proof_obligations']}`", "", "## Decision", payload["decision"]["bounded_result"], "", "## Recommendation", payload["decision"]["next_honest_step"]]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3088/S2038 spectral-to-Hamiltonian/time-evolution obstruction audit", "## P3088/S2038 spectral-to-Hamiltonian/time-evolution obstruction audit\n\n`P3088/S2038` attacks exactly one post-P3087 interface atom: a non-imported spectral-to-Hamiltonian/time-evolution source for the Z12 Dirichlet/Laplacian branch.  It enumerates `12` spectral Hamiltonian rows, computes `6` formal unitary phase rows over a dimensionless time grid, verifies `18` time-energy scale compensation rows, and builds a `5 x 6 = 30` candidate gate matrix.  The finite self-adjoint and formal-unitary algebra remains dimensionless; no sourced time parameter, action/energy unit, observable dynamics, physical Hamiltonian, observed radiation/light, spacetime EOM, `L_total`, bridge/role-transfer, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3088/S2038 Hamiltonian/time evolution remains unsourced", "## P3088/S2038 Hamiltonian/time evolution remains unsourced\n\n`P3088/S2038` confirms that a finite Z12 Laplacian can be read as a real symmetric/self-adjoint spectral multiplier and can support formal unit-modulus phase factors.  A Lagrangian/EOM, quantum, or empirical-physics reading still needs strict sources for time, action/energy units, Hamiltonian dynamics, and observable readout; imported quantum mechanics or spacetime field-equation templates do not supply those strict sources.\n")
    append_once(AGENTS, "Current spectral-to-Hamiltonian/time-evolution obstruction guardrail (P3088/S2038, 2026-06-25)", "## Current spectral-to-Hamiltonian/time-evolution obstruction guardrail (P3088/S2038, 2026-06-25)\n\n- P3088 follows the P3087 recommendation and audits one standard-physics interface atom: a non-imported spectral-to-Hamiltonian/time-evolution source for the Z12 Dirichlet/Laplacian branch.\n- The finite audit computes `12` spectral Hamiltonian rows, `6` formal unitary phase rows, `18` time-energy scale orbit rows, and `30` candidate gate rows; `0` candidates export an internal non-imported Hamiltonian/time-evolution observable source.\n- Do not promote finite self-adjoint spectral multipliers, formal `exp(-i lambda t)` phases, dimensionless time parameters, time-energy-normalized phase orbits, imported quantum-mechanics templates, or imported spacetime/EOM templates to observed photons/light, spacetime EOM, physical Hamiltonian, empirical observable, selector closure, `L_total`, bridge/role-transfer, or ToE closure.\n- The next honest move is exactly one spectral-observable/Born-rule probability-readout obstruction audit for the Z12 Dirichlet/Laplacian branch, unless a genuinely new Hamiltonian/time source theorem is introduced.\n")
    return payload

if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
