#!/usr/bin/env python3
"""P3096/S2046: scattering/S-matrix asymptotic-state audit.

P3095 left finite propagation-like witnesses but no sourced physical propagating
observable or observed-radiation interface.  P3096 attacks exactly one adjacent
standard-physics interface atom: whether the Z12 Dirichlet/Laplacian branch
internally sources in/out asymptotic states, an S-matrix/cross-section, or
empirical detector semantics without importing continuum scattering theory,
spacetime asymptotics, apparatus units, selector closure, L_total,
bridge/role-transfer, or ToE.
"""
from __future__ import annotations

import hashlib, json, math, subprocess
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p3095_s2045_dispersion_propagating_mode_observable_obstruction_audit import OUT as P3095, Z12_SIZE, eigenvalue

OUT = GEN / "p3096_s2046_scattering_smatrix_asymptotic_state_obstruction_audit.json"
MD = GEN / "p3096_s2046_scattering_smatrix_asymptotic_state_obstruction_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"
POTENTIAL_MASK = 0b001001001001
COUPLING_GRID = (0.1, 0.25, 0.5)
CHANNELS = tuple(range(Z12_SIZE))

CONTENT_PATTERNS = {
    "scattering_smatrix_atom": r"scattering|S-matrix|asymptotic state|cross-section|detector semantics",
    "predecessor_propagation_atom": r"dispersion|propagating-mode|propagating field|group-velocity|observed-light|radiation interface",
    "blocked_promotions": r"continuum scattering theory|spacetime asymptotics|apparatus units|observed photons|observed light|L_total|ToE|selector closure|bridge/role-transfer",
}

SCATTERING_CANDIDATES = (
    {"id": "finite_fourier_channel_catalog", "description": "finite Z12 Fourier channels with degenerate energy labels", "finite_channel_basis_exported": True, "in_out_asymptotic_states_exported": False, "transition_amplitude_exported": False, "s_matrix_exported": False, "cross_section_exported": False, "unit_normalized_scattering_amplitude_exported": False, "empirical_detector_semantics_exported": False, "nonimported_physical_scattering_source_exported": False, "blocker": "finite channels are not asymptotic in/out states and carry no detector/readout semantics"},
    {"id": "finite_born_transition_matrix", "description": "finite source-coupled Born-like transition amplitudes <k'|V|k>", "finite_channel_basis_exported": True, "in_out_asymptotic_states_exported": False, "transition_amplitude_exported": True, "s_matrix_exported": False, "cross_section_exported": False, "unit_normalized_scattering_amplitude_exported": False, "empirical_detector_semantics_exported": False, "nonimported_physical_scattering_source_exported": False, "blocker": "transition matrix is finite/formal and lacks asymptotic states, units, and detector semantics"},
    {"id": "formal_unitarity_defect_s_matrix_proxy", "description": "formal S=I+i g T proxy with finite unitarity-defect rows", "finite_channel_basis_exported": True, "in_out_asymptotic_states_exported": False, "transition_amplitude_exported": True, "s_matrix_exported": True, "cross_section_exported": False, "unit_normalized_scattering_amplitude_exported": False, "empirical_detector_semantics_exported": False, "nonimported_physical_scattering_source_exported": False, "blocker": "S proxy is algebraic and not a physical unitary scattering operator"},
    {"id": "imported_continuum_scattering_template", "description": "external continuum scattering/asymptotic-state template", "finite_channel_basis_exported": True, "in_out_asymptotic_states_exported": True, "transition_amplitude_exported": True, "s_matrix_exported": True, "cross_section_exported": True, "unit_normalized_scattering_amplitude_exported": False, "empirical_detector_semantics_exported": False, "nonimported_physical_scattering_source_exported": False, "blocker": "asymptotic and S-matrix semantics are imported and lack internal units/readout"},
    {"id": "imported_empirical_detector_scattering_template", "description": "external empirically calibrated detector/cross-section template", "finite_channel_basis_exported": True, "in_out_asymptotic_states_exported": True, "transition_amplitude_exported": True, "s_matrix_exported": True, "cross_section_exported": True, "unit_normalized_scattering_amplitude_exported": True, "empirical_detector_semantics_exported": True, "nonimported_physical_scattering_source_exported": False, "blocker": "passes only by importing detector units and empirical scattering semantics"},
)
REQUIRED_GATES = ("finite_channel_basis_exported", "in_out_asymptotic_states_exported", "transition_amplitude_exported", "s_matrix_exported", "cross_section_exported", "unit_normalized_scattering_amplitude_exported", "empirical_detector_semantics_exported", "nonimported_physical_scattering_source_exported")


def content_grep() -> list[dict[str, Any]]:
    rows = []
    for lane, pattern in CONTENT_PATTERNS.items():
        proc = subprocess.run(["rg", "-n", pattern, "fundamental_action_reconstruction", "AGENTS.md", "-S"], cwd=REPO, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
        hits = [line for line in proc.stdout.splitlines() if line.strip()]
        rows.append({"content_lane": lane, "pattern": pattern, "hit_count": len(hits), "sample_hits": hits[:8]})
    return rows


def potential(site: int) -> int:
    return int((POTENTIAL_MASK >> site) & 1)


def born_amplitude(k_in: int, k_out: int) -> complex:
    total = 0j
    for x in range(Z12_SIZE):
        phase = 2.0 * math.pi * (k_in - k_out) * x / Z12_SIZE
        total += potential(x) * complex(math.cos(phase), math.sin(phase))
    return total / Z12_SIZE


def channel_rows() -> list[dict[str, Any]]:
    rows = []
    for k in CHANNELS:
        rows.append({"channel": k, "laplacian_energy_label": round(eigenvalue(k), 12), "opposite_channel": (-k) % Z12_SIZE, "finite_channel_basis_witness": True, "in_state_attached": False, "out_state_attached": False, "asymptotic_region_attached": False})
    return rows


def transition_rows() -> list[dict[str, Any]]:
    rows = []
    for kin in CHANNELS:
        for kout in CHANNELS:
            amp = born_amplitude(kin, kout)
            rows.append({"in_channel": kin, "out_channel": kout, "amplitude_real": round(amp.real, 12), "amplitude_imag": round(amp.imag, 12), "amplitude_abs_squared": round((amp.real ** 2 + amp.imag ** 2), 12), "energy_label_difference": round(eigenvalue(kout) - eigenvalue(kin), 12), "finite_transition_amplitude_witness": True, "unit_normalized_amplitude_attached": False, "cross_section_semantics_attached": False})
    return rows


def unitarity_proxy_rows() -> list[dict[str, Any]]:
    rows = []
    diagonal_trace = sum(born_amplitude(k, k).real for k in CHANNELS)
    frob_t2 = sum(abs(born_amplitude(kin, kout)) ** 2 for kin in CHANNELS for kout in CHANNELS)
    for g in COUPLING_GRID:
        # For S = I + i g T with Hermitian T, S^dagger S - I = g^2 T^2.
        rows.append({"dimensionless_coupling": g, "formal_s_proxy": "I + i*g*T_born", "trace_T": round(diagonal_trace, 12), "frobenius_norm_unitarity_defect_proxy": round((g ** 2) * frob_t2, 12), "finite_s_matrix_proxy_witness": True, "exact_unitarity_attached": False, "physical_scattering_operator_attached": False})
    return rows


def cross_section_proxy_rows() -> list[dict[str, Any]]:
    rows = []
    for kin in CHANNELS:
        total = sum(abs(born_amplitude(kin, kout)) ** 2 for kout in CHANNELS)
        elastic = abs(born_amplitude(kin, kin)) ** 2
        rows.append({"in_channel": kin, "total_formal_cross_section_proxy": round(total, 12), "elastic_proxy": round(elastic, 12), "inelastic_proxy": round(total - elastic, 12), "finite_cross_section_proxy_witness": True, "area_unit_attached": False, "detector_semantics_attached": False})
    return rows


def gate_rows() -> list[dict[str, Any]]:
    return [{"candidate": c["id"], "required_gate": gate, "gate_passed": bool(c[gate]), "blocking_if_failed": not bool(c[gate]), "detail": "passed" if c[gate] else c["blocker"]} for c in SCATTERING_CANDIDATES for gate in REQUIRED_GATES]


def aggregate(gates: list[dict[str, Any]]) -> list[dict[str, Any]]:
    out = []
    for c in SCATTERING_CANDIDATES:
        subset = [r for r in gates if r["candidate"] == c["id"]]
        out.append({"candidate": c["id"], "passed_gates": sum(r["gate_passed"] for r in subset), "failed_gates": sum(not r["gate_passed"] for r in subset), "accepted_nonimported_scattering_smatrix_source": all(r["gate_passed"] for r in subset) and bool(c["nonimported_physical_scattering_source_exported"])})
    return out


def build_payload() -> dict[str, Any]:
    p3095 = read_json(P3095)
    greps = content_grep(); channels = channel_rows(); transitions = transition_rows(); unitarity = unitarity_proxy_rows(); xs = cross_section_proxy_rows(); gates = gate_rows(); aggs = aggregate(gates)
    accepted = [r for r in aggs if r["accepted_nonimported_scattering_smatrix_source"]]
    obligations = [
        {"obligation": "read_p3095_next_atom", "satisfied": True, "detail": "P3095 selected scattering/S-matrix asymptotic-state audit as the next interface atom"},
        {"obligation": "construct_finite_channel_and_transition_rows", "satisfied": len(channels) == Z12_SIZE and len(transitions) == Z12_SIZE * Z12_SIZE, "detail": "finite channels and Born-like transition amplitudes are explicit"},
        {"obligation": "construct_s_matrix_unitarity_proxy_rows", "satisfied": len(unitarity) == len(COUPLING_GRID) and all(r["finite_s_matrix_proxy_witness"] for r in unitarity), "detail": "formal S-proxy unitarity-defect rows are computed"},
        {"obligation": "construct_cross_section_proxy_rows", "satisfied": len(xs) == Z12_SIZE and all(r["finite_cross_section_proxy_witness"] for r in xs), "detail": "formal cross-section proxy rows are explicit"},
        {"obligation": "export_nonimported_physical_scattering_source", "satisfied": False, "detail": "0 candidates pass as internal non-imported scattering/S-matrix sources"},
    ]
    return {
        "status": "P3096_SCATTERING_SMATRIX_ASYMPTOTIC_STATE_OBSTRUCTION_BOUNDED_NO_GO",
        "input_hashes": {"P3095": hashlib.sha256(P3095.read_bytes()).hexdigest() if P3095.exists() else None},
        "constructed_theoretical_objects": {"content_first_repo_grep": greps, "scattering_smatrix_audit_object": {"object": "Z12DirichletScatteringSMatrixAsymptoticStateObstructionAudit", "source_reused": "P3095 recommendation: bounded scattering/S-matrix asymptotic-state obstruction audit", "required_gates": list(REQUIRED_GATES), "candidate_scattering_sources": [c["id"] for c in SCATTERING_CANDIDATES], "acceptance_predicate": "finite channel basis plus in/out asymptotic states, transition amplitudes, S-matrix, cross-section, unit-normalized amplitudes, empirical detector semantics, and non-imported physical scattering source"}, "finite_channel_rows": channels, "born_transition_amplitude_rows": transitions, "s_matrix_unitarity_proxy_rows": unitarity, "cross_section_proxy_rows": xs, "candidate_gate_rows": gates, "candidate_aggregate_certificate": aggs},
        "finite_certificate": {"content_grep_lanes": len(greps), "content_grep_hits": sum(r["hit_count"] for r in greps), "p3095_accepted_nonimported_dispersion_propagating_observable_sources": p3095.get("finite_certificate", {}).get("accepted_nonimported_dispersion_propagating_observable_sources"), "finite_channel_rows": len(channels), "channel_rows_with_in_state": sum(r["in_state_attached"] for r in channels), "channel_rows_with_out_state": sum(r["out_state_attached"] for r in channels), "channel_rows_with_asymptotic_region": sum(r["asymptotic_region_attached"] for r in channels), "born_transition_amplitude_rows": len(transitions), "transition_rows_with_unit_normalized_amplitude": sum(r["unit_normalized_amplitude_attached"] for r in transitions), "transition_rows_with_cross_section_semantics": sum(r["cross_section_semantics_attached"] for r in transitions), "s_matrix_unitarity_proxy_rows": len(unitarity), "s_proxy_rows_with_exact_unitarity": sum(r["exact_unitarity_attached"] for r in unitarity), "s_proxy_rows_with_physical_scattering_operator": sum(r["physical_scattering_operator_attached"] for r in unitarity), "cross_section_proxy_rows": len(xs), "cross_section_rows_with_area_unit": sum(r["area_unit_attached"] for r in xs), "cross_section_rows_with_detector_semantics": sum(r["detector_semantics_attached"] for r in xs), "scattering_candidates": len(SCATTERING_CANDIDATES), "required_gates": len(REQUIRED_GATES), "candidate_gate_rows": len(gates), "accepted_nonimported_scattering_smatrix_sources": len(accepted), "proof_obligations": len(obligations), "satisfied_proof_obligations": sum(r["satisfied"] for r in obligations)},
        "proof_obligations": obligations,
        "decision": {"bounded_result": "P3096 constructs the requested scattering/S-matrix asymptotic-state obstruction audit.  The Z12 Laplacian supports finite Fourier channel labels, Born-like transition amplitudes for a finite source potential, formal S=I+i g T unitarity-defect proxies, and formal cross-section proxies.  These are real scattering-like witnesses, but no internal artifact exports in/out asymptotic states, a physical unitary S-matrix, unit-normalized scattering amplitudes, empirical detector semantics, spacetime asymptotics, or a non-imported physical scattering source.  Imported continuum scattering and empirical detector templates pass only as imported templates.  Therefore no non-imported scattering/S-matrix source is exported.", "negative_export_flags": {key: False for key in ["in_out_asymptotic_states_exported", "physical_unitary_smatrix_exported", "physical_cross_section_exported", "unit_normalized_scattering_amplitude_exported", "empirical_detector_semantics_exported", "spacetime_asymptotics_exported", "nonimported_physical_scattering_source_exported", "physical_propagating_field_mode_exported", "detector_independent_observable_exported", "observed_light_radiation_interface_exported", "empirical_observable_exported", "observed_radiation_exported", "observed_photons_exported", "observed_light_exported", "spacetime_eom_exported", "physical_hamiltonian_exported", "green_response_closure_exported", "effective_action_closure_exported", "rg_scale_flow_closure_exported", "ward_current_charge_closure_exported", "stress_metric_response_closure_exported", "dispersion_propagation_closure_exported", "QW_2191_discharged", "strict_selector_closure_exported", "ltotal_exported", "bridge_role_transfer_exported", "toe_closure_exported"]}, "positive_scoped_flags": {"finite_channel_rows_computed": True, "born_transition_amplitude_rows_computed": True, "s_matrix_unitarity_proxy_rows_computed": True, "cross_section_proxy_rows_computed": True}, "next_honest_step": "Pivot to exactly one adjacent standard-physics interface atom outside selector replay: construct a bounded thermodynamic-radiation/blackbody-spectrum obstruction audit for the Z12 Dirichlet/Laplacian branch, testing whether finite spectral mode counts, formal propagation/scattering proxies, partition weights, and energy-flux rows supply a non-imported radiation spectrum, temperature/energy unit calibration, photon/light interpretation, and empirical intensity readout without importing continuum statistical field theory, apparatus units, observed light, L_total, bridge/role-transfer, or ToE."},
    }


def write_markdown(payload: dict[str, Any]) -> None:
    c = payload["finite_certificate"]
    lines = ["# P3096/S2046 scattering/S-matrix asymptotic-state obstruction audit", "", f"Status: `{payload['status']}`", "", "## Finite certificate", f"- content grep lanes: `{c['content_grep_lanes']}`", f"- content grep hits: `{c['content_grep_hits']}`", f"- P3095 accepted non-imported dispersion/propagating observable sources: `{c['p3095_accepted_nonimported_dispersion_propagating_observable_sources']}`", f"- finite channel rows: `{c['finite_channel_rows']}`", f"- channel rows with in-state: `{c['channel_rows_with_in_state']}`", f"- channel rows with out-state: `{c['channel_rows_with_out_state']}`", f"- channel rows with asymptotic region: `{c['channel_rows_with_asymptotic_region']}`", f"- Born transition amplitude rows: `{c['born_transition_amplitude_rows']}`", f"- transition rows with unit-normalized amplitude: `{c['transition_rows_with_unit_normalized_amplitude']}`", f"- transition rows with cross-section semantics: `{c['transition_rows_with_cross_section_semantics']}`", f"- S-matrix unitarity proxy rows: `{c['s_matrix_unitarity_proxy_rows']}`", f"- S-proxy rows with exact unitarity: `{c['s_proxy_rows_with_exact_unitarity']}`", f"- S-proxy rows with physical scattering operator: `{c['s_proxy_rows_with_physical_scattering_operator']}`", f"- cross-section proxy rows: `{c['cross_section_proxy_rows']}`", f"- cross-section rows with area unit: `{c['cross_section_rows_with_area_unit']}`", f"- cross-section rows with detector semantics: `{c['cross_section_rows_with_detector_semantics']}`", f"- scattering candidates: `{c['scattering_candidates']}`", f"- candidate gate rows: `{c['candidate_gate_rows']}`", f"- accepted non-imported scattering/S-matrix sources: `{c['accepted_nonimported_scattering_smatrix_sources']}`", f"- satisfied proof obligations: `{c['satisfied_proof_obligations']}/{c['proof_obligations']}`", "", "## Decision", payload["decision"]["bounded_result"], "", "## Recommendation", payload["decision"]["next_honest_step"]]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3096/S2046 scattering/S-matrix asymptotic-state obstruction audit", "## P3096/S2046 scattering/S-matrix asymptotic-state obstruction audit\n\n`P3096/S2046` attacks exactly one post-P3095 interface atom: a non-imported scattering/S-matrix asymptotic-state source for the Z12 Dirichlet/Laplacian branch.  It constructs `12` finite channel rows, `144` Born-like transition-amplitude rows, `3` S-matrix unitarity-proxy rows, `12` cross-section proxy rows, and a `5 x 8 = 40` candidate gate matrix.  The finite scattering-like algebra remains formal; no in/out asymptotic states, physical unitary S-matrix, unit-normalized scattering amplitudes, empirical detector semantics, spacetime asymptotics, `L_total`, bridge/role-transfer, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3096/S2046 scattering/S-matrix source remains unsourced", "## P3096/S2046 scattering/S-matrix source remains unsourced\n\n`P3096/S2046` confirms that the Z12 Laplacian supports finite Fourier channel labels, Born-like transition amplitudes, formal S-proxy unitarity-defect rows, and cross-section proxies.  A Lagrangian/EOM, scattering operator, or detector-readout reading still needs strict sources for asymptotic in/out states, unit-normalized amplitudes, spacetime asymptotics, empirical detector semantics, and observed scattering interpretation; imported continuum scattering templates do not supply those strict sources.\n")
    append_once(AGENTS, "Current scattering/S-matrix asymptotic-state obstruction guardrail (P3096/S2046, 2026-06-25)", "## Current scattering/S-matrix asymptotic-state obstruction guardrail (P3096/S2046, 2026-06-25)\n\n- P3096 follows the P3095 recommendation and audits one standard-physics interface atom: a non-imported scattering/S-matrix asymptotic-state source for the Z12 Dirichlet/Laplacian branch.\n- The finite audit computes `12` finite channel rows, `144` Born-like transition-amplitude rows, `3` S-matrix unitarity-proxy rows, `12` cross-section proxy rows, and `40` candidate gate rows; `0` candidates export an internal non-imported scattering/S-matrix law.\n- Do not promote finite Fourier channel labels, Born-like transition amplitudes, formal S-proxy unitarity-defect rows, cross-section proxies, imported continuum scattering templates, or imported detector templates to in/out asymptotic states, physical unitary S-matrix, empirical observable, observed photons/light, spacetime EOM, physical Hamiltonian, selector closure, `L_total`, bridge/role-transfer, or ToE closure.\n- The next honest move is exactly one thermodynamic-radiation/blackbody-spectrum obstruction audit for the Z12 Dirichlet/Laplacian branch, unless a genuinely new scattering/S-matrix theorem is introduced.\n")
    return payload

if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
