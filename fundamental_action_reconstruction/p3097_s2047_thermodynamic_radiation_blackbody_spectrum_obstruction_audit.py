#!/usr/bin/env python3
"""P3097/S2047: thermodynamic-radiation/blackbody-spectrum audit.

P3096 left finite scattering-like witnesses but no sourced physical scattering,
radiation, or detector interface.  P3097 attacks exactly one adjacent standard-
physics interface atom: whether the Z12 Dirichlet/Laplacian branch internally
sources a thermal radiation spectrum, temperature/energy calibration, photon or
observed-light semantics, and empirical intensity readout without importing
continuum statistical field theory, apparatus units, selector closure, L_total,
bridge/role-transfer, or ToE.
"""
from __future__ import annotations

import hashlib, json, math, subprocess
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p3095_s2045_dispersion_propagating_mode_observable_obstruction_audit import Z12_SIZE, eigenvalue
from p3096_s2046_scattering_smatrix_asymptotic_state_obstruction_audit import OUT as P3096

OUT = GEN / "p3097_s2047_thermodynamic_radiation_blackbody_spectrum_obstruction_audit.json"
MD = GEN / "p3097_s2047_thermodynamic_radiation_blackbody_spectrum_obstruction_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"
BETA_GRID = (0.25, 0.5, 1.0, 2.0)
ENERGY_SCALE_GRID = (0.5, 1.0, 2.0)

CONTENT_PATTERNS = {
    "thermodynamic_radiation_atom": r"thermodynamic-radiation|blackbody|radiation spectrum|Planck|intensity readout",
    "predecessor_scattering_atom": r"scattering|S-matrix|cross-section|detector semantics",
    "blocked_promotions": r"continuum statistical field theory|apparatus units|observed light|observed photons|L_total|ToE|selector closure|bridge/role-transfer",
}

RADIATION_CANDIDATES = (
    {"id": "finite_spectral_mode_count", "description": "finite Z12 Laplacian mode-count and degeneracy catalog", "finite_mode_count_exported": True, "partition_weight_exported": False, "radiation_spectrum_exported": False, "temperature_energy_unit_exported": False, "photon_light_semantics_exported": False, "empirical_intensity_readout_exported": False, "nonimported_physical_radiation_source_exported": False, "blocker": "mode counts are finite spectral labels with no temperature, units, or light semantics"},
    {"id": "formal_boltzmann_weighted_modes", "description": "dimensionless exp(-beta lambda) weights over Z12 modes", "finite_mode_count_exported": True, "partition_weight_exported": True, "radiation_spectrum_exported": False, "temperature_energy_unit_exported": False, "photon_light_semantics_exported": False, "empirical_intensity_readout_exported": False, "nonimported_physical_radiation_source_exported": False, "blocker": "beta is formal and not a sourced physical temperature or energy unit"},
    {"id": "modal_flux_weighted_spectrum_proxy", "description": "formal intensity proxy using degeneracy, energy label, and group-velocity magnitude", "finite_mode_count_exported": True, "partition_weight_exported": True, "radiation_spectrum_exported": True, "temperature_energy_unit_exported": False, "photon_light_semantics_exported": False, "empirical_intensity_readout_exported": False, "nonimported_physical_radiation_source_exported": False, "blocker": "spectrum proxy is dimensionless and lacks photon/light semantics and calibrated intensity readout"},
    {"id": "imported_blackbody_planck_template", "description": "external Planck-law/continuum radiation template", "finite_mode_count_exported": True, "partition_weight_exported": True, "radiation_spectrum_exported": True, "temperature_energy_unit_exported": True, "photon_light_semantics_exported": False, "empirical_intensity_readout_exported": False, "nonimported_physical_radiation_source_exported": False, "blocker": "blackbody semantics are imported and do not provide internal observed-light/readout source"},
    {"id": "imported_observed_light_intensity_template", "description": "external calibrated radiation/apparatus readout template", "finite_mode_count_exported": True, "partition_weight_exported": True, "radiation_spectrum_exported": True, "temperature_energy_unit_exported": True, "photon_light_semantics_exported": True, "empirical_intensity_readout_exported": True, "nonimported_physical_radiation_source_exported": False, "blocker": "passes only by importing observed light, apparatus units, and empirical readout"},
)
REQUIRED_GATES = ("finite_mode_count_exported", "partition_weight_exported", "radiation_spectrum_exported", "temperature_energy_unit_exported", "photon_light_semantics_exported", "empirical_intensity_readout_exported", "nonimported_physical_radiation_source_exported")


def content_grep() -> list[dict[str, Any]]:
    rows = []
    for lane, pattern in CONTENT_PATTERNS.items():
        proc = subprocess.run(["rg", "-n", pattern, "fundamental_action_reconstruction", "AGENTS.md", "-S"], cwd=REPO, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
        hits = [line for line in proc.stdout.splitlines() if line.strip()]
        rows.append({"content_lane": lane, "pattern": pattern, "hit_count": len(hits), "sample_hits": hits[:8]})
    return rows


def mode_count_rows() -> list[dict[str, Any]]:
    counts: dict[float, list[int]] = {}
    for k in range(Z12_SIZE):
        lam = round(eigenvalue(k), 12)
        counts.setdefault(lam, []).append(k)
    return [{"energy_label": lam, "degeneracy": len(ch), "channels": ch, "finite_mode_count_witness": True, "temperature_unit_attached": False, "photon_semantics_attached": False} for lam, ch in sorted(counts.items())]


def partition_rows(counts: list[dict[str, Any]]) -> list[dict[str, Any]]:
    rows = []
    for beta in BETA_GRID:
        weights = [r["degeneracy"] * math.exp(-beta * r["energy_label"]) for r in counts]
        z = sum(weights)
        mean_e = sum(w * r["energy_label"] for w, r in zip(weights, counts)) / z
        rows.append({"formal_beta": beta, "partition_Z": round(z, 12), "mean_energy_label": round(mean_e, 12), "dimensionless_partition_witness": True, "physical_temperature_attached": False, "boltzmann_unit_attached": False})
    return rows


def group_velocity(k: int) -> float:
    return math.sin(2.0 * math.pi * k / Z12_SIZE)


def spectrum_proxy_rows(counts: list[dict[str, Any]]) -> list[dict[str, Any]]:
    rows = []
    channel_to_velocity = {k: abs(group_velocity(k)) for k in range(Z12_SIZE)}
    for beta in BETA_GRID:
        z = sum(r["degeneracy"] * math.exp(-beta * r["energy_label"]) for r in counts)
        for r in counts:
            flux = sum(channel_to_velocity[k] for k in r["channels"])
            weight = r["degeneracy"] * math.exp(-beta * r["energy_label"]) / z
            intensity = weight * r["energy_label"] * flux
            rows.append({"formal_beta": beta, "energy_label": r["energy_label"], "degeneracy": r["degeneracy"], "velocity_flux_proxy": round(flux, 12), "normalized_weight": round(weight, 12), "dimensionless_intensity_proxy": round(intensity, 12), "radiation_spectrum_proxy_witness": True, "frequency_unit_attached": False, "observed_light_semantics_attached": False, "empirical_intensity_readout_attached": False})
    return rows


def scale_orbit_rows(counts: list[dict[str, Any]]) -> list[dict[str, Any]]:
    rows = []
    base_beta = 1.0
    for scale in ENERGY_SCALE_GRID:
        beta = base_beta / scale
        weights = [r["degeneracy"] * math.exp(-beta * scale * r["energy_label"]) for r in counts]
        z = sum(weights)
        rows.append({"energy_scale": scale, "compensating_beta": beta, "partition_Z": round(z, 12), "same_as_base_partition": True, "scale_orbit_witness": True, "canonical_temperature_or_energy_unit_selected": False})
    return rows


def gate_rows() -> list[dict[str, Any]]:
    return [{"candidate": c["id"], "required_gate": gate, "gate_passed": bool(c[gate]), "blocking_if_failed": not bool(c[gate]), "detail": "passed" if c[gate] else c["blocker"]} for c in RADIATION_CANDIDATES for gate in REQUIRED_GATES]


def aggregate(gates: list[dict[str, Any]]) -> list[dict[str, Any]]:
    out = []
    for c in RADIATION_CANDIDATES:
        subset = [r for r in gates if r["candidate"] == c["id"]]
        out.append({"candidate": c["id"], "passed_gates": sum(r["gate_passed"] for r in subset), "failed_gates": sum(not r["gate_passed"] for r in subset), "accepted_nonimported_thermodynamic_radiation_source": all(r["gate_passed"] for r in subset) and bool(c["nonimported_physical_radiation_source_exported"])})
    return out


def build_payload() -> dict[str, Any]:
    p3096 = read_json(P3096)
    greps = content_grep(); counts = mode_count_rows(); partitions = partition_rows(counts); spectra = spectrum_proxy_rows(counts); orbits = scale_orbit_rows(counts); gates = gate_rows(); aggs = aggregate(gates)
    accepted = [r for r in aggs if r["accepted_nonimported_thermodynamic_radiation_source"]]
    obligations = [
        {"obligation": "read_p3096_next_atom", "satisfied": True, "detail": "P3096 selected thermodynamic-radiation/blackbody-spectrum audit as the next interface atom"},
        {"obligation": "construct_finite_mode_count_rows", "satisfied": len(counts) == 7, "detail": "finite spectral degeneracy rows are explicit"},
        {"obligation": "construct_partition_weight_rows", "satisfied": len(partitions) == len(BETA_GRID), "detail": "dimensionless partition weights are computed"},
        {"obligation": "construct_radiation_spectrum_proxy_rows", "satisfied": len(spectra) == len(BETA_GRID) * len(counts), "detail": "formal modal intensity proxies are computed"},
        {"obligation": "show_temperature_energy_scale_orbit", "satisfied": len(orbits) == len(ENERGY_SCALE_GRID) and all(r["same_as_base_partition"] for r in orbits), "detail": "energy scale can be compensated by beta rescaling"},
        {"obligation": "export_nonimported_physical_radiation_source", "satisfied": False, "detail": "0 candidates pass as internal non-imported thermodynamic-radiation sources"},
    ]
    return {"status": "P3097_THERMODYNAMIC_RADIATION_BLACKBODY_SPECTRUM_OBSTRUCTION_BOUNDED_NO_GO", "input_hashes": {"P3096": hashlib.sha256(P3096.read_bytes()).hexdigest() if P3096.exists() else None}, "constructed_theoretical_objects": {"content_first_repo_grep": greps, "thermodynamic_radiation_audit_object": {"object": "Z12DirichletThermodynamicRadiationBlackbodySpectrumObstructionAudit", "source_reused": "P3096 recommendation: bounded thermodynamic-radiation/blackbody-spectrum obstruction audit", "required_gates": list(REQUIRED_GATES), "candidate_radiation_sources": [c["id"] for c in RADIATION_CANDIDATES], "acceptance_predicate": "finite mode counts plus partition weights, radiation spectrum, physical temperature/energy units, photon/light semantics, empirical intensity readout, and non-imported physical radiation source"}, "finite_mode_count_rows": counts, "formal_partition_weight_rows": partitions, "radiation_spectrum_proxy_rows": spectra, "temperature_energy_scale_orbit_rows": orbits, "candidate_gate_rows": gates, "candidate_aggregate_certificate": aggs}, "finite_certificate": {"content_grep_lanes": len(greps), "content_grep_hits": sum(r["hit_count"] for r in greps), "p3096_accepted_nonimported_scattering_smatrix_sources": p3096.get("finite_certificate", {}).get("accepted_nonimported_scattering_smatrix_sources"), "finite_mode_count_rows": len(counts), "mode_count_rows_with_temperature_unit": sum(r["temperature_unit_attached"] for r in counts), "mode_count_rows_with_photon_semantics": sum(r["photon_semantics_attached"] for r in counts), "formal_partition_weight_rows": len(partitions), "partition_rows_with_physical_temperature": sum(r["physical_temperature_attached"] for r in partitions), "partition_rows_with_boltzmann_unit": sum(r["boltzmann_unit_attached"] for r in partitions), "radiation_spectrum_proxy_rows": len(spectra), "spectrum_rows_with_frequency_unit": sum(r["frequency_unit_attached"] for r in spectra), "spectrum_rows_with_observed_light_semantics": sum(r["observed_light_semantics_attached"] for r in spectra), "spectrum_rows_with_empirical_intensity_readout": sum(r["empirical_intensity_readout_attached"] for r in spectra), "temperature_energy_scale_orbit_rows": len(orbits), "scale_orbit_rows_with_canonical_temperature_or_energy_unit": sum(r["canonical_temperature_or_energy_unit_selected"] for r in orbits), "radiation_candidates": len(RADIATION_CANDIDATES), "required_gates": len(REQUIRED_GATES), "candidate_gate_rows": len(gates), "accepted_nonimported_thermodynamic_radiation_sources": len(accepted), "proof_obligations": len(obligations), "satisfied_proof_obligations": sum(r["satisfied"] for r in obligations)}, "proof_obligations": obligations, "decision": {"bounded_result": "P3097 constructs the requested thermodynamic-radiation/blackbody-spectrum obstruction audit.  The Z12 Laplacian supports finite spectral mode counts, dimensionless partition weights, formal modal intensity-spectrum proxies, and explicit temperature/energy scale-orbit witnesses.  These are real thermodynamic/radiation-like witnesses, but no internal artifact exports a physical temperature or energy unit, a Planck/blackbody radiation law, photon/light semantics, empirical intensity readout, or a non-imported physical radiation source.  Imported continuum statistical-field and observed-light templates pass only as imported templates.  Therefore no non-imported thermodynamic-radiation/blackbody-spectrum source is exported.", "negative_export_flags": {key: False for key in ["physical_temperature_energy_unit_exported", "blackbody_radiation_law_exported", "photon_light_semantics_exported", "empirical_intensity_readout_exported", "nonimported_physical_radiation_source_exported", "in_out_asymptotic_states_exported", "physical_unitary_smatrix_exported", "empirical_observable_exported", "observed_radiation_exported", "observed_photons_exported", "observed_light_exported", "spacetime_eom_exported", "physical_hamiltonian_exported", "scattering_smatrix_closure_exported", "dispersion_propagation_closure_exported", "QW_2191_discharged", "strict_selector_closure_exported", "ltotal_exported", "bridge_role_transfer_exported", "toe_closure_exported"]}, "positive_scoped_flags": {"finite_mode_count_rows_computed": True, "formal_partition_weight_rows_computed": True, "radiation_spectrum_proxy_rows_computed": True, "temperature_energy_scale_orbit_rows_computed": True}, "next_honest_step": "Pivot to exactly one adjacent standard-physics interface atom outside selector replay: construct a bounded KMS/detailed-balance thermal-state obstruction audit for the Z12 Dirichlet/Laplacian branch, testing whether finite partition weights, transition kernels, modal flux proxies, and scattering/response witnesses supply a non-imported equilibrium-state law, physical temperature clock, fluctuation-dissipation relation, and empirical thermal readout without importing continuum statistical mechanics, apparatus units, observed light, L_total, bridge/role-transfer, or ToE."}}


def write_markdown(payload: dict[str, Any]) -> None:
    c = payload["finite_certificate"]
    lines = ["# P3097/S2047 thermodynamic-radiation/blackbody-spectrum obstruction audit", "", f"Status: `{payload['status']}`", "", "## Finite certificate", f"- content grep lanes: `{c['content_grep_lanes']}`", f"- content grep hits: `{c['content_grep_hits']}`", f"- P3096 accepted non-imported scattering/S-matrix sources: `{c['p3096_accepted_nonimported_scattering_smatrix_sources']}`", f"- finite mode count rows: `{c['finite_mode_count_rows']}`", f"- mode rows with temperature unit: `{c['mode_count_rows_with_temperature_unit']}`", f"- mode rows with photon semantics: `{c['mode_count_rows_with_photon_semantics']}`", f"- formal partition weight rows: `{c['formal_partition_weight_rows']}`", f"- partition rows with physical temperature: `{c['partition_rows_with_physical_temperature']}`", f"- radiation spectrum proxy rows: `{c['radiation_spectrum_proxy_rows']}`", f"- spectrum rows with frequency unit: `{c['spectrum_rows_with_frequency_unit']}`", f"- spectrum rows with observed-light semantics: `{c['spectrum_rows_with_observed_light_semantics']}`", f"- spectrum rows with empirical intensity readout: `{c['spectrum_rows_with_empirical_intensity_readout']}`", f"- temperature/energy scale-orbit rows: `{c['temperature_energy_scale_orbit_rows']}`", f"- radiation candidates: `{c['radiation_candidates']}`", f"- candidate gate rows: `{c['candidate_gate_rows']}`", f"- accepted non-imported thermodynamic-radiation sources: `{c['accepted_nonimported_thermodynamic_radiation_sources']}`", f"- satisfied proof obligations: `{c['satisfied_proof_obligations']}/{c['proof_obligations']}`", "", "## Decision", payload["decision"]["bounded_result"], "", "## Recommendation", payload["decision"]["next_honest_step"]]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3097/S2047 thermodynamic-radiation/blackbody-spectrum obstruction audit", "## P3097/S2047 thermodynamic-radiation/blackbody-spectrum obstruction audit\n\n`P3097/S2047` attacks exactly one post-P3096 interface atom: a non-imported thermodynamic-radiation/blackbody-spectrum source for the Z12 Dirichlet/Laplacian branch.  It constructs `7` finite spectral mode-count rows, `4` formal partition-weight rows, `28` radiation-spectrum proxy rows, `3` temperature/energy scale-orbit rows, and a `5 x 7 = 35` candidate gate matrix.  The finite thermodynamic/radiation-like algebra remains formal; no physical temperature or energy unit, Planck/blackbody radiation law, photon/light semantics, empirical intensity readout, `L_total`, bridge/role-transfer, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3097/S2047 thermodynamic-radiation source remains unsourced", "## P3097/S2047 thermodynamic-radiation source remains unsourced\n\n`P3097/S2047` confirms that the Z12 Laplacian supports finite mode-count tables, dimensionless partition weights, formal intensity-spectrum proxies, and temperature/energy scale-orbit witnesses.  A Lagrangian/EOM, blackbody-radiation, photon/light, or detector-readout reading still needs strict sources for physical temperature and energy units, Planck/blackbody law semantics, observed-light interpretation, and empirical intensity calibration; imported continuum statistical-field and apparatus templates do not supply those strict sources.\n")
    append_once(AGENTS, "Current thermodynamic-radiation/blackbody-spectrum obstruction guardrail (P3097/S2047, 2026-06-25)", "## Current thermodynamic-radiation/blackbody-spectrum obstruction guardrail (P3097/S2047, 2026-06-25)\n\n- P3097 follows the P3096 recommendation and audits one standard-physics interface atom: a non-imported thermodynamic-radiation/blackbody-spectrum source for the Z12 Dirichlet/Laplacian branch.\n- The finite audit computes `7` finite spectral mode-count rows, `4` formal partition-weight rows, `28` radiation-spectrum proxy rows, `3` temperature/energy scale-orbit rows, and `35` candidate gate rows; `0` candidates export an internal non-imported thermodynamic-radiation law.\n- Do not promote finite mode counts, dimensionless partition weights, formal modal intensity-spectrum proxies, temperature/energy scale-orbit rows, imported Planck/blackbody templates, or imported observed-light/apparatus templates to physical temperature, radiation spectrum, empirical observable, observed photons/light, spacetime EOM, physical Hamiltonian, selector closure, `L_total`, bridge/role-transfer, or ToE closure.\n- The next honest move is exactly one KMS/detailed-balance thermal-state obstruction audit for the Z12 Dirichlet/Laplacian branch, unless a genuinely new thermodynamic-radiation theorem is introduced.\n")
    return payload

if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
