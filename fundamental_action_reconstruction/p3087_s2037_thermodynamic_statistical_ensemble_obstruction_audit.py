#!/usr/bin/env python3
"""P3087/S2037: thermodynamic/statistical-ensemble obstruction audit.

P3086 left the Z12 Dirichlet/Laplacian branch with finite dimensionless
spectral/scalar witnesses but no unit-calibrated empirical observable.  P3087
attacks exactly one adjacent standard-physics interface atom: whether the
internal finite spectrum sources a thermodynamic/statistical ensemble with a
canonical temperature, partition function, entropy/energy units, and equilibrium
observable without importing Boltzmann constants, measurement units, observed
radiation, spacetime EOM, selector closure, L_total, bridge/role-transfer, or
ToE.
"""
from __future__ import annotations

import hashlib, json, math, subprocess
from collections import Counter
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p3086_s2036_empirical_readout_observable_calibration_obstruction_audit import OUT as P3086

OUT = GEN / "p3087_s2037_thermodynamic_statistical_ensemble_obstruction_audit.json"
MD = GEN / "p3087_s2037_thermodynamic_statistical_ensemble_obstruction_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"
Z12_SIZE = 12
BETA_GRID = (0.0, 0.125, 0.25, 0.5, 1.0, 2.0, 4.0)
ENERGY_SCALES = (0.5, 1.0, 2.0)

CONTENT_PATTERNS = {
    "thermodynamic_atom": r"thermodynamic|statistical-ensemble|partition function|temperature|Boltzmann|entropy/energy",
    "predecessor_observable_atom": r"empirical observable|observable-calibration|measurement units|unit-calibrated|readout",
    "blocked_promotions": r"observed photons|observed light|observed radiation|spacetime EOM|L_total|ToE|selector closure|bridge/role-transfer",
}

ENSEMBLE_CANDIDATES = (
    {
        "id": "finite_z12_laplacian_boltzmann_weights",
        "description": "formal canonical weights exp(-beta lambda_m) on the Z12 Laplacian spectrum",
        "internal_spectrum_exported": True,
        "canonical_temperature_source_exported": False,
        "boltzmann_constant_or_unit_source_exported": False,
        "partition_function_computed": True,
        "entropy_energy_units_exported": False,
        "equilibrium_observable_exported": False,
        "blocker": "formal beta is dimensionless and not sourced as physical inverse temperature",
    },
    {
        "id": "z12_microcanonical_degeneracy_counts",
        "description": "finite degeneracy counts of distinct Z12 Laplacian energy labels",
        "internal_spectrum_exported": True,
        "canonical_temperature_source_exported": False,
        "boltzmann_constant_or_unit_source_exported": False,
        "partition_function_computed": False,
        "entropy_energy_units_exported": False,
        "equilibrium_observable_exported": False,
        "blocker": "degeneracy counting gives dimensionless combinatorics, not thermodynamic units",
    },
    {
        "id": "p3086_scale_orbit_normalized_energy_family",
        "description": "energy-scale orbit inherited from P3086 observable-calibration controls",
        "internal_spectrum_exported": True,
        "canonical_temperature_source_exported": False,
        "boltzmann_constant_or_unit_source_exported": False,
        "partition_function_computed": True,
        "entropy_energy_units_exported": False,
        "equilibrium_observable_exported": False,
        "blocker": "energy rescalings change beta-energy products unless an external unit convention is chosen",
    },
    {
        "id": "imported_boltzmann_gibbs_template",
        "description": "external continuum/statistical mechanics Boltzmann-Gibbs ensemble template",
        "internal_spectrum_exported": False,
        "canonical_temperature_source_exported": True,
        "boltzmann_constant_or_unit_source_exported": True,
        "partition_function_computed": True,
        "entropy_energy_units_exported": True,
        "equilibrium_observable_exported": True,
        "blocker": "passes only by importing thermodynamic units and ensemble machinery",
    },
    {
        "id": "imported_blackbody_radiation_template",
        "description": "external observed-radiation thermodynamic readout template",
        "internal_spectrum_exported": False,
        "canonical_temperature_source_exported": True,
        "boltzmann_constant_or_unit_source_exported": True,
        "partition_function_computed": True,
        "entropy_energy_units_exported": True,
        "equilibrium_observable_exported": True,
        "blocker": "passes only by importing observed radiation/temperature physics",
    },
)
REQUIRED_GATES = ("internal_spectrum_exported", "canonical_temperature_source_exported", "boltzmann_constant_or_unit_source_exported", "partition_function_computed", "entropy_energy_units_exported", "equilibrium_observable_exported")


def content_grep() -> list[dict[str, Any]]:
    rows = []
    for lane, pattern in CONTENT_PATTERNS.items():
        proc = subprocess.run(["rg", "-n", pattern, "fundamental_action_reconstruction", "AGENTS.md", "-S"], cwd=REPO, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
        hits = [line for line in proc.stdout.splitlines() if line.strip()]
        rows.append({"content_lane": lane, "pattern": pattern, "hit_count": len(hits), "sample_hits": hits[:8]})
    return rows


def z12_spectrum() -> list[float]:
    return [2.0 - 2.0 * math.cos(2.0 * math.pi * m / Z12_SIZE) for m in range(Z12_SIZE)]


def degeneracy_rows(spectrum: list[float]) -> list[dict[str, Any]]:
    counts = Counter(round(v, 12) for v in spectrum)
    return [{"energy_label": energy, "degeneracy": counts[energy], "dimensionless_count": True, "entropy_unit_attached": False} for energy in sorted(counts)]


def partition_rows(spectrum: list[float]) -> list[dict[str, Any]]:
    rows = []
    for beta in BETA_GRID:
        weights = [math.exp(-beta * energy) for energy in spectrum]
        z = sum(weights)
        mean_e = sum(energy * weight for energy, weight in zip(spectrum, weights)) / z
        mean_e2 = sum((energy ** 2) * weight for energy, weight in zip(spectrum, weights)) / z
        variance = mean_e2 - mean_e ** 2
        entropy_like = math.log(z) + beta * mean_e
        rows.append({
            "dimensionless_beta": beta,
            "partition_function_Z": round(z, 12),
            "mean_dimensionless_energy": round(mean_e, 12),
            "energy_variance": round(variance, 12),
            "entropy_like_logZ_plus_betaE": round(entropy_like, 12),
            "temperature_unit_attached": False,
            "energy_unit_attached": False,
            "physical_equilibrium_observable": False,
        })
    return rows


def scale_beta_orbit_rows(spectrum: list[float]) -> list[dict[str, Any]]:
    rows = []
    for scale in ENERGY_SCALES:
        scaled = [scale * energy for energy in spectrum]
        for beta in BETA_GRID:
            z_scaled = sum(math.exp(-beta * energy) for energy in scaled)
            z_compensated = sum(math.exp(-(beta * scale) * energy) for energy in spectrum)
            rows.append({
                "energy_scale": scale,
                "dimensionless_beta": beta,
                "scaled_partition_function": round(z_scaled, 12),
                "compensated_original_partition_function": round(z_compensated, 12),
                "scale_beta_compensation_identity_holds": abs(z_scaled - z_compensated) <= 1e-10,
                "canonical_energy_scale_selected": scale == 1.0,
                "selection_is_conventional_without_unit_source": True,
            })
    return rows


def gate_rows() -> list[dict[str, Any]]:
    return [{"candidate": c["id"], "required_gate": gate, "gate_passed": bool(c[gate]), "blocking_if_failed": not bool(c[gate]), "detail": "passed" if c[gate] else c["blocker"]} for c in ENSEMBLE_CANDIDATES for gate in REQUIRED_GATES]


def aggregate(gates: list[dict[str, Any]]) -> list[dict[str, Any]]:
    out = []
    for c in ENSEMBLE_CANDIDATES:
        subset = [r for r in gates if r["candidate"] == c["id"]]
        out.append({"candidate": c["id"], "passed_gates": sum(r["gate_passed"] for r in subset), "failed_gates": sum(not r["gate_passed"] for r in subset), "accepted_nonimported_thermodynamic_ensemble_source": all(r["gate_passed"] for r in subset) and bool(c["internal_spectrum_exported"])})
    return out


def build_payload() -> dict[str, Any]:
    p3086 = read_json(P3086)
    greps = content_grep(); spectrum = z12_spectrum(); degeneracies = degeneracy_rows(spectrum); partitions = partition_rows(spectrum); scale_rows = scale_beta_orbit_rows(spectrum); gates = gate_rows(); aggs = aggregate(gates)
    accepted = [r for r in aggs if r["accepted_nonimported_thermodynamic_ensemble_source"]]
    obligations = [
        {"obligation": "read_p3086_next_atom", "satisfied": True, "detail": "P3086 selected thermodynamic/statistical-ensemble audit as the next interface atom"},
        {"obligation": "construct_finite_degeneracy_table", "satisfied": True, "detail": "distinct Z12 Laplacian energy labels and degeneracies are enumerated"},
        {"obligation": "compute_formal_partition_grid", "satisfied": True, "detail": "dimensionless partition/mean-energy/variance/entropy-like rows are computed on seven beta values"},
        {"obligation": "test_energy_scale_beta_orbit", "satisfied": all(r["scale_beta_compensation_identity_holds"] for r in scale_rows), "detail": "energy-scale changes can be compensated by beta rescaling, so no canonical physical temperature is selected"},
        {"obligation": "export_nonimported_thermodynamic_observable", "satisfied": False, "detail": "0 candidates pass as internal non-imported thermodynamic/equilibrium observable sources"},
    ]
    return {
        "status": "P3087_THERMODYNAMIC_STATISTICAL_ENSEMBLE_OBSTRUCTION_BOUNDED_NO_GO",
        "input_hashes": {"P3086": hashlib.sha256(P3086.read_bytes()).hexdigest() if P3086.exists() else None},
        "constructed_theoretical_objects": {
            "content_first_repo_grep": greps,
            "thermodynamic_ensemble_audit_object": {"object": "Z12DirichletThermodynamicStatisticalEnsembleObstructionAudit", "source_reused": "P3086 recommendation: bounded thermodynamic/statistical-ensemble obstruction audit", "required_gates": list(REQUIRED_GATES), "candidate_ensemble_sources": [c["id"] for c in ENSEMBLE_CANDIDATES], "acceptance_predicate": "internal spectrum plus canonical temperature source, Boltzmann/unit source, partition function, entropy/energy units, and equilibrium observable export"},
            "z12_laplacian_spectrum": [round(v, 12) for v in spectrum],
            "microcanonical_degeneracy_rows": degeneracies,
            "formal_partition_function_rows": partitions,
            "energy_scale_beta_orbit_rows": scale_rows,
            "candidate_gate_rows": gates,
            "candidate_aggregate_certificate": aggs,
        },
        "finite_certificate": {
            "content_grep_lanes": len(greps), "content_grep_hits": sum(r["hit_count"] for r in greps),
            "p3086_accepted_nonimported_empirical_observable_sources": p3086.get("finite_certificate", {}).get("accepted_nonimported_empirical_observable_sources"),
            "spectrum_rows": len(spectrum), "microcanonical_degeneracy_rows": len(degeneracies), "formal_partition_function_rows": len(partitions),
            "partition_rows_with_temperature_units": sum(r["temperature_unit_attached"] for r in partitions), "partition_rows_with_energy_units": sum(r["energy_unit_attached"] for r in partitions),
            "energy_scale_beta_orbit_rows": len(scale_rows), "scale_beta_identity_failures": sum(not r["scale_beta_compensation_identity_holds"] for r in scale_rows), "canonical_temperature_sources_exported": 0,
            "ensemble_candidates": len(ENSEMBLE_CANDIDATES), "required_gates": len(REQUIRED_GATES), "candidate_gate_rows": len(gates), "accepted_nonimported_thermodynamic_ensemble_sources": len(accepted),
            "proof_obligations": len(obligations), "satisfied_proof_obligations": sum(r["satisfied"] for r in obligations),
        },
        "proof_obligations": obligations,
        "decision": {
            "bounded_result": "P3087 constructs the requested thermodynamic/statistical-ensemble obstruction audit.  The finite Z12 Laplacian spectrum supports a formal microcanonical degeneracy table and a dimensionless canonical partition grid, but beta remains a formal parameter.  Energy-scale/beta orbit rows show that rescaling the energy labels can be exactly compensated by rescaling beta, so no canonical temperature or energy unit is selected internally.  Boltzmann-Gibbs and blackbody rows pass only as imported templates.  Therefore no non-imported thermodynamic/equilibrium observable source is exported.",
            "negative_export_flags": {key: False for key in ["canonical_temperature_source_exported", "boltzmann_constant_or_unit_source_exported", "entropy_energy_units_exported", "equilibrium_observable_exported", "measurement_unit_source_exported", "empirical_observable_exported", "observed_radiation_exported", "observed_photons_exported", "observed_light_exported", "spacetime_eom_exported", "physical_hamiltonian_exported", "QW_2191_discharged", "strict_selector_closure_exported", "ltotal_exported", "bridge_role_transfer_exported", "toe_closure_exported"]},
            "positive_scoped_flags": {"finite_degeneracy_table_computed": True, "formal_partition_grid_computed": True, "energy_scale_beta_orbit_verified": True},
            "next_honest_step": "Pivot to exactly one adjacent standard-physics interface atom outside selector replay: construct a bounded spectral-to-Hamiltonian/time-evolution obstruction audit for the Z12 Dirichlet/Laplacian branch, testing whether the internal finite spectrum supplies a self-adjoint Hamiltonian with a sourced time parameter, unitary evolution, energy units, and observable dynamics without importing quantum mechanics, measurement units, spacetime EOM, observed radiation/light, L_total, bridge/role-transfer, or ToE.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    c = payload["finite_certificate"]
    lines = ["# P3087/S2037 thermodynamic/statistical-ensemble obstruction audit", "", f"Status: `{payload['status']}`", "", "## Finite certificate",
        f"- content grep lanes: `{c['content_grep_lanes']}`", f"- content grep hits: `{c['content_grep_hits']}`", f"- P3086 accepted non-imported empirical observable sources: `{c['p3086_accepted_nonimported_empirical_observable_sources']}`", f"- spectrum rows: `{c['spectrum_rows']}`", f"- microcanonical degeneracy rows: `{c['microcanonical_degeneracy_rows']}`", f"- formal partition function rows: `{c['formal_partition_function_rows']}`", f"- partition rows with temperature units: `{c['partition_rows_with_temperature_units']}`", f"- partition rows with energy units: `{c['partition_rows_with_energy_units']}`", f"- energy-scale/beta orbit rows: `{c['energy_scale_beta_orbit_rows']}`", f"- scale-beta identity failures: `{c['scale_beta_identity_failures']}`", f"- canonical temperature sources exported: `{c['canonical_temperature_sources_exported']}`", f"- ensemble candidates: `{c['ensemble_candidates']}`", f"- candidate gate rows: `{c['candidate_gate_rows']}`", f"- accepted non-imported thermodynamic ensemble sources: `{c['accepted_nonimported_thermodynamic_ensemble_sources']}`", f"- satisfied proof obligations: `{c['satisfied_proof_obligations']}/{c['proof_obligations']}`", "", "## Decision", payload["decision"]["bounded_result"], "", "## Recommendation", payload["decision"]["next_honest_step"]]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3087/S2037 thermodynamic/statistical-ensemble obstruction audit", "## P3087/S2037 thermodynamic/statistical-ensemble obstruction audit\n\n`P3087/S2037` attacks exactly one post-P3086 interface atom: a non-imported thermodynamic/statistical-ensemble source for the Z12 Dirichlet/Laplacian branch.  It enumerates `7` microcanonical degeneracy rows, computes `7` formal partition-function rows over a dimensionless beta grid, verifies `21` energy-scale/beta compensation rows, and builds a `5 x 6 = 30` candidate gate matrix.  The finite ensemble algebra remains dimensionless; no canonical temperature source, Boltzmann/unit source, entropy/energy units, equilibrium observable, observed radiation/light, spacetime EOM, `L_total`, bridge/role-transfer, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3087/S2037 thermodynamic ensemble remains unsourced", "## P3087/S2037 thermodynamic ensemble remains unsourced\n\n`P3087/S2037` confirms that a finite Z12 Laplacian spectrum can support formal degeneracy and partition-function algebra, but not yet a physical thermodynamic sector.  A Lagrangian/EOM or empirical-physics reading still needs a strict source for temperature, Boltzmann/energy units, entropy units, and equilibrium observables; imported statistical mechanics or blackbody radiation templates do not supply those strict sources.\n")
    append_once(AGENTS, "Current thermodynamic/statistical-ensemble obstruction guardrail (P3087/S2037, 2026-06-25)", "## Current thermodynamic/statistical-ensemble obstruction guardrail (P3087/S2037, 2026-06-25)\n\n- P3087 follows the P3086 recommendation and audits one standard-physics interface atom: a non-imported thermodynamic/statistical-ensemble source for the Z12 Dirichlet/Laplacian branch.\n- The finite audit computes `7` microcanonical degeneracy rows, `7` formal partition-function rows, `21` energy-scale/beta orbit rows, and `30` candidate gate rows; `0` candidates export an internal non-imported thermodynamic/equilibrium observable source.\n- Do not promote finite degeneracy counts, dimensionless partition functions, formal beta parameters, energy-scale-normalized ensembles, imported Boltzmann-Gibbs templates, or imported blackbody/observed-radiation templates to observed photons/light, spacetime EOM, physical Hamiltonian, empirical observable, selector closure, `L_total`, bridge/role-transfer, or ToE closure.\n- The next honest move is exactly one spectral-to-Hamiltonian/time-evolution obstruction audit for the Z12 Dirichlet/Laplacian branch, unless a genuinely new thermodynamic source theorem is introduced.\n")
    return payload

if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
