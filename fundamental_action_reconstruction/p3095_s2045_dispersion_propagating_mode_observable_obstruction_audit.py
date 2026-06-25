#!/usr/bin/env python3
"""P3095/S2045: dispersion/propagating-mode observable audit.

P3094 left finite stress/metric-response-like witnesses but no sourced physical
stress-energy or metric-response law.  P3095 attacks exactly one adjacent
standard-physics interface atom: whether the Z12 Dirichlet/Laplacian branch
internally sources a propagating field mode, dispersion relation, detector-
independent observable, or observed-light/radiation interface without importing
continuum waves, spacetime metric, apparatus units, selector closure, L_total,
bridge/role-transfer, or ToE.
"""
from __future__ import annotations

import hashlib, json, math, subprocess
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p3094_s2044_stress_energy_metric_response_obstruction_audit import OUT as P3094, Z12_SIZE

OUT = GEN / "p3095_s2045_dispersion_propagating_mode_observable_obstruction_audit.json"
MD = GEN / "p3095_s2045_dispersion_propagating_mode_observable_obstruction_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"
TIME_GRID = (0.0, 0.5, 1.0, 2.0)
PACKET_CENTERS = (1, 2, 3)
SIGMA = 1.25
MASS2_GRID = (0.0, 0.25, 1.0, 4.0)

CONTENT_PATTERNS = {
    "dispersion_propagating_mode_atom": r"dispersion|propagating-mode|propagating field|group-velocity|observed-light|radiation interface",
    "predecessor_stress_metric_atom": r"stress-energy|stress tensor|metric-response|metric coupling|gravitational/field-response",
    "blocked_promotions": r"continuum waves|spacetime metric|apparatus units|observed photons|observed light|L_total|ToE|selector closure|bridge/role-transfer",
}

PROPAGATION_CANDIDATES = (
    {"id": "finite_laplacian_dispersion_curve", "description": "finite Z12 Laplacian eigenvalue curve lambda(k)=2-2cos(2pi k/12)", "finite_dispersion_relation_exported": True, "group_velocity_proxy_exported": True, "mode_packet_evolution_exported": False, "green_pole_structure_exported": False, "propagating_field_mode_exported": False, "detector_independent_observable_exported": False, "observed_light_radiation_interface_exported": False, "nonimported_physical_propagation_source_exported": False, "blocker": "finite eigenvalue curve is a spectral graph relation, not a physical propagation law or observed-light interface"},
    {"id": "formal_unitary_mode_packet_evolution", "description": "dimensionless exp(-i lambda_k t) packet phases on finite modes", "finite_dispersion_relation_exported": True, "group_velocity_proxy_exported": True, "mode_packet_evolution_exported": True, "green_pole_structure_exported": False, "propagating_field_mode_exported": False, "detector_independent_observable_exported": False, "observed_light_radiation_interface_exported": False, "nonimported_physical_propagation_source_exported": False, "blocker": "packet phases use a formal dimensionless time and lack spacetime speed, detector observable, or radiation semantics"},
    {"id": "mass_regularized_green_pole_catalog", "description": "finite spectral pole catalog for resolvents 1/(lambda_k+m^2)", "finite_dispersion_relation_exported": True, "group_velocity_proxy_exported": False, "mode_packet_evolution_exported": False, "green_pole_structure_exported": True, "propagating_field_mode_exported": False, "detector_independent_observable_exported": False, "observed_light_radiation_interface_exported": False, "nonimported_physical_propagation_source_exported": False, "blocker": "poles are finite spectral denominators and not asymptotic propagation/readout states"},
    {"id": "imported_continuum_wave_dispersion_template", "description": "external continuum wave/field dispersion template", "finite_dispersion_relation_exported": True, "group_velocity_proxy_exported": True, "mode_packet_evolution_exported": True, "green_pole_structure_exported": True, "propagating_field_mode_exported": True, "detector_independent_observable_exported": False, "observed_light_radiation_interface_exported": False, "nonimported_physical_propagation_source_exported": False, "blocker": "propagating-field semantics are imported and lack detector/readout and observed radiation"},
    {"id": "imported_observed_light_radiation_template", "description": "external observed photon/light/radiation template", "finite_dispersion_relation_exported": True, "group_velocity_proxy_exported": True, "mode_packet_evolution_exported": True, "green_pole_structure_exported": True, "propagating_field_mode_exported": True, "detector_independent_observable_exported": True, "observed_light_radiation_interface_exported": True, "nonimported_physical_propagation_source_exported": False, "blocker": "passes only by importing observed light, detector units, and radiation semantics"},
)
REQUIRED_GATES = ("finite_dispersion_relation_exported", "group_velocity_proxy_exported", "mode_packet_evolution_exported", "green_pole_structure_exported", "propagating_field_mode_exported", "detector_independent_observable_exported", "observed_light_radiation_interface_exported", "nonimported_physical_propagation_source_exported")


def content_grep() -> list[dict[str, Any]]:
    rows = []
    for lane, pattern in CONTENT_PATTERNS.items():
        proc = subprocess.run(["rg", "-n", pattern, "fundamental_action_reconstruction", "AGENTS.md", "-S"], cwd=REPO, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
        hits = [line for line in proc.stdout.splitlines() if line.strip()]
        rows.append({"content_lane": lane, "pattern": pattern, "hit_count": len(hits), "sample_hits": hits[:8]})
    return rows


def wave_number(k: int) -> float:
    return 2.0 * math.pi * k / Z12_SIZE


def eigenvalue(k: int) -> float:
    return 2.0 - 2.0 * math.cos(wave_number(k))


def group_velocity_proxy(k: int) -> float:
    return 2.0 * math.sin(wave_number(k))


def dispersion_rows() -> list[dict[str, Any]]:
    rows = []
    for k in range(Z12_SIZE):
        rows.append({"mode": k, "wave_number": round(wave_number(k), 12), "laplacian_eigenvalue": round(eigenvalue(k), 12), "group_velocity_proxy_dlambda_dq": round(group_velocity_proxy(k), 12), "opposite_mode": (-k) % Z12_SIZE, "finite_dispersion_witness": True, "spacetime_speed_attached": False, "observed_light_semantics_attached": False})
    return rows


def packet_weight(k: int, center: int) -> float:
    circular_distance = min((k - center) % Z12_SIZE, (center - k) % Z12_SIZE)
    return math.exp(-(circular_distance ** 2) / (2.0 * SIGMA ** 2))


def mode_packet_rows() -> list[dict[str, Any]]:
    rows = []
    for center in PACKET_CENTERS:
        norm = sum(packet_weight(k, center) for k in range(Z12_SIZE))
        mean_velocity = sum(packet_weight(k, center) * group_velocity_proxy(k) for k in range(Z12_SIZE)) / norm
        for t in TIME_GRID:
            phase_cos = sum(packet_weight(k, center) * math.cos(eigenvalue(k) * t) for k in range(Z12_SIZE)) / norm
            phase_sin = -sum(packet_weight(k, center) * math.sin(eigenvalue(k) * t) for k in range(Z12_SIZE)) / norm
            rows.append({"packet_center_mode": center, "dimensionless_time": t, "mean_group_velocity_proxy": round(mean_velocity, 12), "normalized_phase_real": round(phase_cos, 12), "normalized_phase_imag": round(phase_sin, 12), "formal_mode_packet_evolution_witness": True, "physical_time_unit_attached": False, "detector_observable_attached": False})
    return rows


def green_pole_rows() -> list[dict[str, Any]]:
    rows = []
    for m2 in MASS2_GRID:
        for k in range(Z12_SIZE):
            denom = eigenvalue(k) + m2
            rows.append({"dimensionless_mass_squared": m2, "mode": k, "pole_location_negative_lambda_minus_m2": round(-denom, 12), "resolvent_denominator": round(denom, 12), "singular_denominator": abs(denom) < 1e-12, "green_pole_catalog_witness": True, "asymptotic_state_attached": False, "scattering_readout_attached": False})
    return rows


def energy_flux_proxy_rows() -> list[dict[str, Any]]:
    rows = []
    for k in range(Z12_SIZE):
        lam = eigenvalue(k)
        vel = group_velocity_proxy(k)
        rows.append({"mode": k, "modal_energy_proxy": round(lam, 12), "group_velocity_proxy": round(vel, 12), "energy_flux_proxy": round(lam * vel, 12), "stress_energy_link_witness": True, "physical_flux_unit_attached": False, "radiation_interface_attached": False})
    return rows


def gate_rows() -> list[dict[str, Any]]:
    return [{"candidate": c["id"], "required_gate": gate, "gate_passed": bool(c[gate]), "blocking_if_failed": not bool(c[gate]), "detail": "passed" if c[gate] else c["blocker"]} for c in PROPAGATION_CANDIDATES for gate in REQUIRED_GATES]


def aggregate(gates: list[dict[str, Any]]) -> list[dict[str, Any]]:
    out = []
    for c in PROPAGATION_CANDIDATES:
        subset = [r for r in gates if r["candidate"] == c["id"]]
        out.append({"candidate": c["id"], "passed_gates": sum(r["gate_passed"] for r in subset), "failed_gates": sum(not r["gate_passed"] for r in subset), "accepted_nonimported_dispersion_propagating_observable_source": all(r["gate_passed"] for r in subset) and bool(c["nonimported_physical_propagation_source_exported"])})
    return out


def build_payload() -> dict[str, Any]:
    p3094 = read_json(P3094)
    greps = content_grep(); dispersion = dispersion_rows(); packets = mode_packet_rows(); poles = green_pole_rows(); flux = energy_flux_proxy_rows(); gates = gate_rows(); aggs = aggregate(gates)
    accepted = [r for r in aggs if r["accepted_nonimported_dispersion_propagating_observable_source"]]
    obligations = [
        {"obligation": "read_p3094_next_atom", "satisfied": True, "detail": "P3094 selected dispersion/propagating-mode observable audit as the next interface atom"},
        {"obligation": "construct_finite_dispersion_and_group_velocity_rows", "satisfied": len(dispersion) == Z12_SIZE and all(r["finite_dispersion_witness"] for r in dispersion), "detail": "finite Z12 Laplacian dispersion and group-velocity proxy rows are explicit"},
        {"obligation": "construct_mode_packet_and_green_pole_rows", "satisfied": len(packets) == len(PACKET_CENTERS) * len(TIME_GRID) and len(poles) == len(MASS2_GRID) * Z12_SIZE, "detail": "formal packet evolution and Green-pole catalog rows are computed"},
        {"obligation": "construct_energy_flux_proxy_rows", "satisfied": len(flux) == Z12_SIZE and all(r["stress_energy_link_witness"] for r in flux), "detail": "modal energy-flux proxy rows are explicit"},
        {"obligation": "export_nonimported_physical_propagating_observable_source", "satisfied": False, "detail": "0 candidates pass as internal non-imported dispersion/propagating observable sources"},
    ]
    return {
        "status": "P3095_DISPERSION_PROPAGATING_MODE_OBSERVABLE_OBSTRUCTION_BOUNDED_NO_GO",
        "input_hashes": {"P3094": hashlib.sha256(P3094.read_bytes()).hexdigest() if P3094.exists() else None},
        "constructed_theoretical_objects": {"content_first_repo_grep": greps, "dispersion_propagating_observable_audit_object": {"object": "Z12DirichletDispersionPropagatingModeObservableObstructionAudit", "source_reused": "P3094 recommendation: bounded dispersion/propagating-mode observable obstruction audit", "required_gates": list(REQUIRED_GATES), "candidate_propagation_sources": [c["id"] for c in PROPAGATION_CANDIDATES], "acceptance_predicate": "finite dispersion relation plus group-velocity proxy, mode-packet evolution, Green-pole structure, propagating field mode, detector-independent observable, observed-light/radiation interface, and non-imported physical propagation source"}, "dispersion_group_velocity_rows": dispersion, "mode_packet_evolution_rows": packets, "green_pole_catalog_rows": poles, "energy_flux_proxy_rows": flux, "candidate_gate_rows": gates, "candidate_aggregate_certificate": aggs},
        "finite_certificate": {"content_grep_lanes": len(greps), "content_grep_hits": sum(r["hit_count"] for r in greps), "p3094_accepted_nonimported_stress_energy_metric_response_sources": p3094.get("finite_certificate", {}).get("accepted_nonimported_stress_energy_metric_response_sources"), "dispersion_group_velocity_rows": len(dispersion), "dispersion_rows_with_spacetime_speed": sum(r["spacetime_speed_attached"] for r in dispersion), "dispersion_rows_with_observed_light_semantics": sum(r["observed_light_semantics_attached"] for r in dispersion), "mode_packet_evolution_rows": len(packets), "packet_rows_with_physical_time_unit": sum(r["physical_time_unit_attached"] for r in packets), "packet_rows_with_detector_observable": sum(r["detector_observable_attached"] for r in packets), "green_pole_catalog_rows": len(poles), "green_pole_rows_with_asymptotic_state": sum(r["asymptotic_state_attached"] for r in poles), "green_pole_rows_with_scattering_readout": sum(r["scattering_readout_attached"] for r in poles), "energy_flux_proxy_rows": len(flux), "energy_flux_rows_with_physical_flux_unit": sum(r["physical_flux_unit_attached"] for r in flux), "energy_flux_rows_with_radiation_interface": sum(r["radiation_interface_attached"] for r in flux), "propagation_candidates": len(PROPAGATION_CANDIDATES), "required_gates": len(REQUIRED_GATES), "candidate_gate_rows": len(gates), "accepted_nonimported_dispersion_propagating_observable_sources": len(accepted), "proof_obligations": len(obligations), "satisfied_proof_obligations": sum(r["satisfied"] for r in obligations)},
        "proof_obligations": obligations,
        "decision": {"bounded_result": "P3095 constructs the requested dispersion/propagating-mode observable obstruction audit.  The Z12 Laplacian supports a finite dispersion curve, group-velocity proxies, formal dimensionless mode-packet evolution, mass-regularized Green-pole catalogs, and modal energy-flux proxies.  These are real propagation-like witnesses, but no internal artifact exports a physical propagating field mode, spacetime speed/light cone, detector-independent observable, observed-light/radiation interface, physical time/apparatus units, or a non-imported physical propagation source.  Imported continuum wave and observed-radiation templates pass only as imported templates.  Therefore no non-imported dispersion/propagating observable source is exported.", "negative_export_flags": {key: False for key in ["physical_propagating_field_mode_exported", "physical_dispersion_relation_exported", "spacetime_speed_or_lightcone_exported", "detector_independent_observable_exported", "observed_light_radiation_interface_exported", "nonimported_physical_propagation_source_exported", "physical_stress_energy_tensor_exported", "metric_coupling_exported", "empirical_observable_exported", "observed_radiation_exported", "observed_photons_exported", "observed_light_exported", "spacetime_eom_exported", "physical_hamiltonian_exported", "green_response_closure_exported", "effective_action_closure_exported", "rg_scale_flow_closure_exported", "ward_current_charge_closure_exported", "stress_metric_response_closure_exported", "born_rule_readout_exported", "QW_2191_discharged", "strict_selector_closure_exported", "ltotal_exported", "bridge_role_transfer_exported", "toe_closure_exported"]}, "positive_scoped_flags": {"finite_dispersion_group_velocity_rows_computed": True, "formal_mode_packet_evolution_rows_computed": True, "green_pole_catalog_rows_computed": True, "energy_flux_proxy_rows_computed": True}, "next_honest_step": "Pivot to exactly one adjacent standard-physics interface atom outside selector replay: construct a bounded scattering/S-matrix asymptotic-state obstruction audit for the Z12 Dirichlet/Laplacian branch, testing whether Green-pole catalogs, finite mode packets, source/response kernels, and propagation proxies supply non-imported in/out states, an S-matrix or cross-section, unit-normalized scattering amplitudes, and empirical detector semantics without importing continuum scattering theory, spacetime asymptotics, apparatus units, L_total, bridge/role-transfer, or ToE."},
    }


def write_markdown(payload: dict[str, Any]) -> None:
    c = payload["finite_certificate"]
    lines = ["# P3095/S2045 dispersion/propagating-mode observable obstruction audit", "", f"Status: `{payload['status']}`", "", "## Finite certificate", f"- content grep lanes: `{c['content_grep_lanes']}`", f"- content grep hits: `{c['content_grep_hits']}`", f"- P3094 accepted non-imported stress-energy/metric-response sources: `{c['p3094_accepted_nonimported_stress_energy_metric_response_sources']}`", f"- dispersion/group-velocity rows: `{c['dispersion_group_velocity_rows']}`", f"- dispersion rows with spacetime speed: `{c['dispersion_rows_with_spacetime_speed']}`", f"- dispersion rows with observed-light semantics: `{c['dispersion_rows_with_observed_light_semantics']}`", f"- mode-packet evolution rows: `{c['mode_packet_evolution_rows']}`", f"- packet rows with physical time unit: `{c['packet_rows_with_physical_time_unit']}`", f"- packet rows with detector observable: `{c['packet_rows_with_detector_observable']}`", f"- Green-pole catalog rows: `{c['green_pole_catalog_rows']}`", f"- Green-pole rows with asymptotic state: `{c['green_pole_rows_with_asymptotic_state']}`", f"- Green-pole rows with scattering readout: `{c['green_pole_rows_with_scattering_readout']}`", f"- energy-flux proxy rows: `{c['energy_flux_proxy_rows']}`", f"- energy-flux rows with physical flux unit: `{c['energy_flux_rows_with_physical_flux_unit']}`", f"- energy-flux rows with radiation interface: `{c['energy_flux_rows_with_radiation_interface']}`", f"- propagation candidates: `{c['propagation_candidates']}`", f"- candidate gate rows: `{c['candidate_gate_rows']}`", f"- accepted non-imported dispersion/propagating observable sources: `{c['accepted_nonimported_dispersion_propagating_observable_sources']}`", f"- satisfied proof obligations: `{c['satisfied_proof_obligations']}/{c['proof_obligations']}`", "", "## Decision", payload["decision"]["bounded_result"], "", "## Recommendation", payload["decision"]["next_honest_step"]]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3095/S2045 dispersion/propagating-mode observable obstruction audit", "## P3095/S2045 dispersion/propagating-mode observable obstruction audit\n\n`P3095/S2045` attacks exactly one post-P3094 interface atom: a non-imported dispersion/propagating-mode observable source for the Z12 Dirichlet/Laplacian branch.  It constructs `12` dispersion/group-velocity rows, `12` formal mode-packet evolution rows, `48` Green-pole catalog rows, `12` energy-flux proxy rows, and a `5 x 8 = 40` candidate gate matrix.  The finite propagation-like algebra remains formal/dimensionless; no physical propagating field mode, spacetime speed/light cone, detector-independent observable, observed-light/radiation interface, apparatus units, `L_total`, bridge/role-transfer, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3095/S2045 dispersion/propagating observable source remains unsourced", "## P3095/S2045 dispersion/propagating observable source remains unsourced\n\n`P3095/S2045` confirms that the Z12 Laplacian supports a finite dispersion curve, group-velocity proxies, formal mode-packet phases, Green-pole catalogs, and modal energy-flux proxies.  A Lagrangian/EOM, propagating field, or observed-radiation reading still needs strict sources for physical time and distance units, spacetime speed/light-cone semantics, detector-independent observables, empirical radiation/readout, and observed-light interpretation; imported continuum wave templates do not supply those strict sources.\n")
    append_once(AGENTS, "Current dispersion/propagating-mode observable obstruction guardrail (P3095/S2045, 2026-06-25)", "## Current dispersion/propagating-mode observable obstruction guardrail (P3095/S2045, 2026-06-25)\n\n- P3095 follows the P3094 recommendation and audits one standard-physics interface atom: a non-imported dispersion/propagating-mode observable source for the Z12 Dirichlet/Laplacian branch.\n- The finite audit computes `12` dispersion/group-velocity rows, `12` formal mode-packet evolution rows, `48` Green-pole catalog rows, `12` energy-flux proxy rows, and `40` candidate gate rows; `0` candidates export an internal non-imported propagating observable/radiation law.\n- Do not promote finite dispersion curves, group-velocity proxies, formal mode-packet phases, Green-pole catalogs, modal energy-flux proxies, imported continuum wave templates, or imported observed-radiation templates to physical propagating field mode, detector-independent observable, observed photons/light, spacetime EOM, physical Hamiltonian, selector closure, `L_total`, bridge/role-transfer, or ToE closure.\n- The next honest move is exactly one scattering/S-matrix asymptotic-state obstruction audit for the Z12 Dirichlet/Laplacian branch, unless a genuinely new propagation/observable theorem is introduced.\n")
    return payload

if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
