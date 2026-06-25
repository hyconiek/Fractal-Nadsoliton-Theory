#!/usr/bin/env python3
"""P3094/S2044: stress-energy/metric-response obstruction audit.

P3093 left finite Ward/current/effective-charge-like witnesses but no sourced
physical current or charge law.  P3094 attacks exactly one adjacent
standard-physics interface atom: whether the Z12 Dirichlet/Laplacian branch
internally sources a stress-energy tensor, metric coupling, or metric-response
law without importing spacetime geometry, continuum EOM, observed
radiation/light, apparatus units, selector closure, L_total, bridge/role-
transfer, or ToE.
"""
from __future__ import annotations

import hashlib, json, math, subprocess
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p3093_s2043_ward_identity_symmetry_current_effective_charge_obstruction_audit import OUT as P3093, Z12_SIZE, profile

OUT = GEN / "p3094_s2044_stress_energy_metric_response_obstruction_audit.json"
MD = GEN / "p3094_s2044_stress_energy_metric_response_obstruction_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"
SOURCE_MASKS = (1, 3, 5, 15, 31, 63, 127, 255, 511, 1023, 2047, 4095)
MASS2_GRID = (0.125, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0)
REFERENCE_MASK = 0b001111001111

CONTENT_PATTERNS = {
    "stress_metric_response_atom": r"stress-energy|stress tensor|metric-response|metric coupling|gravitational/field-response",
    "predecessor_ward_current_atom": r"Ward|symmetry-current|conserved current|effective charge|gauge-charge",
    "blocked_promotions": r"spacetime geometry|continuum EOM|observed photons|observed light|L_total|ToE|selector closure|bridge/role-transfer",
}

STRESS_CANDIDATES = (
    {"id": "edge_weight_metric_variation_energy", "description": "finite derivative of graph quadratic energy with respect to edge weights", "finite_metric_variation_exported": True, "graph_energy_quadratic_exported": True, "stress_tensor_candidate_exported": True, "metric_coupling_exported": False, "conservation_law_exported": False, "empirical_gravitational_response_exported": False, "spacetime_geometry_exported": False, "nonimported_physical_stress_source_exported": False, "blocker": "edge-weight variation is a graph-energy derivative, not a physical metric coupling or stress tensor"},
    {"id": "spectral_logdet_pressure_derivative", "description": "mass-regularized log-det derivative under uniform edge/metric scaling", "finite_metric_variation_exported": True, "graph_energy_quadratic_exported": True, "stress_tensor_candidate_exported": True, "metric_coupling_exported": False, "conservation_law_exported": False, "empirical_gravitational_response_exported": False, "spacetime_geometry_exported": False, "nonimported_physical_stress_source_exported": False, "blocker": "spectral pressure-like derivative is formal and unitless"},
    {"id": "formal_edge_stress_divergence_balance", "description": "finite divergence rows for edge stress candidates on one profile", "finite_metric_variation_exported": True, "graph_energy_quadratic_exported": True, "stress_tensor_candidate_exported": True, "metric_coupling_exported": False, "conservation_law_exported": False, "empirical_gravitational_response_exported": False, "spacetime_geometry_exported": False, "nonimported_physical_stress_source_exported": False, "blocker": "divergence rows are not a sourced covariant conservation law"},
    {"id": "imported_continuum_stress_tensor_template", "description": "external continuum variational stress tensor template", "finite_metric_variation_exported": True, "graph_energy_quadratic_exported": True, "stress_tensor_candidate_exported": True, "metric_coupling_exported": True, "conservation_law_exported": True, "empirical_gravitational_response_exported": False, "spacetime_geometry_exported": True, "nonimported_physical_stress_source_exported": False, "blocker": "stress/metric semantics are imported and lack internal empirical response"},
    {"id": "imported_empirical_gravity_response_template", "description": "external empirically calibrated gravitational/field-response template", "finite_metric_variation_exported": True, "graph_energy_quadratic_exported": True, "stress_tensor_candidate_exported": True, "metric_coupling_exported": True, "conservation_law_exported": True, "empirical_gravitational_response_exported": True, "spacetime_geometry_exported": True, "nonimported_physical_stress_source_exported": False, "blocker": "passes only by importing spacetime geometry and empirical field-response semantics"},
)
REQUIRED_GATES = ("finite_metric_variation_exported", "graph_energy_quadratic_exported", "stress_tensor_candidate_exported", "metric_coupling_exported", "conservation_law_exported", "empirical_gravitational_response_exported", "spacetime_geometry_exported", "nonimported_physical_stress_source_exported")


def content_grep() -> list[dict[str, Any]]:
    rows = []
    for lane, pattern in CONTENT_PATTERNS.items():
        proc = subprocess.run(["rg", "-n", pattern, "fundamental_action_reconstruction", "AGENTS.md", "-S"], cwd=REPO, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
        hits = [line for line in proc.stdout.splitlines() if line.strip()]
        rows.append({"content_lane": lane, "pattern": pattern, "hit_count": len(hits), "sample_hits": hits[:8]})
    return rows


def edge_gradient(p: list[int], edge: int) -> int:
    return p[(edge + 1) % Z12_SIZE] - p[edge]


def graph_energy(p: list[int]) -> float:
    return 0.5 * sum(edge_gradient(p, e) ** 2 for e in range(Z12_SIZE))


def metric_variation_rows() -> list[dict[str, Any]]:
    rows = []
    for mask in SOURCE_MASKS:
        p = profile(mask)
        for edge in range(Z12_SIZE):
            grad = edge_gradient(p, edge)
            rows.append({"source_mask": mask, "edge": edge, "edge_gradient": grad, "graph_energy_density": 0.5 * grad * grad, "d_energy_d_edge_weight": 0.5 * grad * grad, "finite_metric_variation_witness": True, "metric_coupling_attached": False, "physical_stress_tensor_attached": False})
    return rows


def energy_quadratic_rows() -> list[dict[str, Any]]:
    rows = []
    for mask in SOURCE_MASKS:
        p = profile(mask)
        rows.append({"source_mask": mask, "support_size": sum(p), "graph_quadratic_energy": graph_energy(p), "laplacian_quadratic_form_witness": True, "action_measure_attached": False, "spacetime_metric_attached": False})
    return rows


def green_value(distance: int, m2: float) -> float:
    return sum(math.cos(2.0 * math.pi * k * distance / Z12_SIZE) / ((2.0 - 2.0 * math.cos(2.0 * math.pi * k / Z12_SIZE)) + m2) for k in range(Z12_SIZE)) / Z12_SIZE


def spectral_pressure_rows() -> list[dict[str, Any]]:
    rows = []
    for m2 in MASS2_GRID:
        edge_logdet_derivative = 2.0 * (green_value(0, m2) - green_value(1, m2))
        rows.append({"dimensionless_mass_squared": m2, "d_logdet_d_uniform_edge_weight_per_edge": round(edge_logdet_derivative, 12), "total_uniform_edge_pressure_like_derivative": round(Z12_SIZE * edge_logdet_derivative, 12), "spectral_pressure_like_witness": True, "physical_pressure_unit_attached": False, "metric_response_law_attached": False})
    return rows


def stress_divergence_rows() -> list[dict[str, Any]]:
    p = profile(REFERENCE_MASK)
    edge_stress = [0.5 * edge_gradient(p, e) ** 2 for e in range(Z12_SIZE)]
    rows = []
    for site in range(Z12_SIZE):
        divergence = edge_stress[site] - edge_stress[(site - 1) % Z12_SIZE]
        rows.append({"site": site, "outgoing_edge_stress": edge_stress[site], "incoming_edge_stress": edge_stress[(site - 1) % Z12_SIZE], "formal_stress_divergence": divergence, "formal_divergence_row_computed": True, "covariant_conservation_law_attached": False, "empirical_field_response_attached": False})
    return rows


def gate_rows() -> list[dict[str, Any]]:
    return [{"candidate": c["id"], "required_gate": gate, "gate_passed": bool(c[gate]), "blocking_if_failed": not bool(c[gate]), "detail": "passed" if c[gate] else c["blocker"]} for c in STRESS_CANDIDATES for gate in REQUIRED_GATES]


def aggregate(gates: list[dict[str, Any]]) -> list[dict[str, Any]]:
    out = []
    for c in STRESS_CANDIDATES:
        subset = [r for r in gates if r["candidate"] == c["id"]]
        out.append({"candidate": c["id"], "passed_gates": sum(r["gate_passed"] for r in subset), "failed_gates": sum(not r["gate_passed"] for r in subset), "accepted_nonimported_stress_energy_metric_response_source": all(r["gate_passed"] for r in subset) and bool(c["nonimported_physical_stress_source_exported"])})
    return out


def build_payload() -> dict[str, Any]:
    p3093 = read_json(P3093)
    greps = content_grep(); metric = metric_variation_rows(); energies = energy_quadratic_rows(); pressure = spectral_pressure_rows(); divergence = stress_divergence_rows(); gates = gate_rows(); aggs = aggregate(gates)
    accepted = [r for r in aggs if r["accepted_nonimported_stress_energy_metric_response_source"]]
    obligations = [
        {"obligation": "read_p3093_next_atom", "satisfied": True, "detail": "P3093 selected stress-energy/metric-response audit as the next interface atom"},
        {"obligation": "construct_metric_variation_rows", "satisfied": len(metric) == len(SOURCE_MASKS) * Z12_SIZE and all(r["finite_metric_variation_witness"] for r in metric), "detail": "edge-weight metric-variation derivatives of graph energy are explicit"},
        {"obligation": "construct_graph_energy_and_spectral_pressure_rows", "satisfied": len(energies) == len(SOURCE_MASKS) and len(pressure) == len(MASS2_GRID), "detail": "quadratic graph energies and mass-regularized pressure-like log-det derivatives are computed"},
        {"obligation": "construct_formal_stress_divergence_rows", "satisfied": len(divergence) == Z12_SIZE and all(r["formal_divergence_row_computed"] for r in divergence), "detail": "formal edge-stress divergence rows are explicit"},
        {"obligation": "export_nonimported_physical_stress_metric_response_source", "satisfied": False, "detail": "0 candidates pass as internal non-imported stress-energy/metric-response sources"},
    ]
    return {
        "status": "P3094_STRESS_ENERGY_METRIC_RESPONSE_OBSTRUCTION_BOUNDED_NO_GO",
        "input_hashes": {"P3093": hashlib.sha256(P3093.read_bytes()).hexdigest() if P3093.exists() else None},
        "constructed_theoretical_objects": {"content_first_repo_grep": greps, "stress_energy_metric_response_audit_object": {"object": "Z12DirichletStressEnergyMetricResponseObstructionAudit", "source_reused": "P3093 recommendation: bounded stress-energy/metric-response obstruction audit", "required_gates": list(REQUIRED_GATES), "candidate_stress_sources": [c["id"] for c in STRESS_CANDIDATES], "acceptance_predicate": "finite metric variation plus graph energy quadratic, stress tensor candidate, metric coupling, conservation law, empirical gravitational/field response, spacetime geometry, and non-imported physical stress source"}, "metric_variation_rows": metric, "graph_energy_quadratic_rows": energies, "spectral_pressure_like_rows": pressure, "formal_stress_divergence_rows": divergence, "candidate_gate_rows": gates, "candidate_aggregate_certificate": aggs},
        "finite_certificate": {"content_grep_lanes": len(greps), "content_grep_hits": sum(r["hit_count"] for r in greps), "p3093_accepted_nonimported_ward_current_effective_charge_sources": p3093.get("finite_certificate", {}).get("accepted_nonimported_ward_current_effective_charge_sources"), "metric_variation_rows": len(metric), "metric_variation_rows_with_metric_coupling": sum(r["metric_coupling_attached"] for r in metric), "metric_variation_rows_with_physical_stress_tensor": sum(r["physical_stress_tensor_attached"] for r in metric), "graph_energy_quadratic_rows": len(energies), "graph_energy_rows_with_action_measure": sum(r["action_measure_attached"] for r in energies), "graph_energy_rows_with_spacetime_metric": sum(r["spacetime_metric_attached"] for r in energies), "spectral_pressure_like_rows": len(pressure), "spectral_pressure_rows_with_physical_pressure_unit": sum(r["physical_pressure_unit_attached"] for r in pressure), "spectral_pressure_rows_with_metric_response_law": sum(r["metric_response_law_attached"] for r in pressure), "formal_stress_divergence_rows": len(divergence), "formal_stress_rows_with_covariant_conservation_law": sum(r["covariant_conservation_law_attached"] for r in divergence), "formal_stress_rows_with_empirical_field_response": sum(r["empirical_field_response_attached"] for r in divergence), "stress_candidates": len(STRESS_CANDIDATES), "required_gates": len(REQUIRED_GATES), "candidate_gate_rows": len(gates), "accepted_nonimported_stress_energy_metric_response_sources": len(accepted), "proof_obligations": len(obligations), "satisfied_proof_obligations": sum(r["satisfied"] for r in obligations)},
        "proof_obligations": obligations,
        "decision": {"bounded_result": "P3094 constructs the requested stress-energy/metric-response obstruction audit.  The Z12 Laplacian supports finite edge-weight metric-variation derivatives, graph-energy quadratic forms, mass-regularized spectral pressure-like log-det derivatives, and formal stress-divergence rows.  These are real stress/metric-like witnesses, but no internal artifact exports a physical stress-energy tensor, metric coupling, covariant conservation law, empirical gravitational/field-response semantics, spacetime geometry, observed radiation/light, or a non-imported physical stress source.  Imported continuum stress-tensor and empirical gravity-response templates pass only as imported templates.  Therefore no non-imported stress-energy/metric-response source is exported.", "negative_export_flags": {key: False for key in ["physical_stress_energy_tensor_exported", "metric_coupling_exported", "covariant_conservation_law_exported", "empirical_gravitational_response_exported", "spacetime_geometry_exported", "nonimported_physical_stress_source_exported", "physical_action_functional_exported", "empirical_observable_exported", "observed_radiation_exported", "observed_photons_exported", "observed_light_exported", "spacetime_eom_exported", "physical_hamiltonian_exported", "green_response_closure_exported", "effective_action_closure_exported", "rg_scale_flow_closure_exported", "ward_current_charge_closure_exported", "born_rule_readout_exported", "QW_2191_discharged", "strict_selector_closure_exported", "ltotal_exported", "bridge_role_transfer_exported", "toe_closure_exported"]}, "positive_scoped_flags": {"finite_metric_variation_rows_computed": True, "graph_energy_quadratic_rows_computed": True, "spectral_pressure_like_rows_computed": True, "formal_stress_divergence_rows_computed": True}, "next_honest_step": "Pivot to exactly one adjacent standard-physics interface atom outside selector replay: construct a bounded dispersion/propagating-mode observable obstruction audit for the Z12 Dirichlet/Laplacian branch, testing whether finite spectral group-velocity proxies, mode packet evolution, Green-response poles, and stress/energy witnesses supply a non-imported propagating field mode, dispersion relation, detector-independent observable, and observed-light/radiation interface without importing continuum waves, spacetime metric, apparatus units, L_total, bridge/role-transfer, or ToE."},
    }


def write_markdown(payload: dict[str, Any]) -> None:
    c = payload["finite_certificate"]
    lines = ["# P3094/S2044 stress-energy/metric-response obstruction audit", "", f"Status: `{payload['status']}`", "", "## Finite certificate", f"- content grep lanes: `{c['content_grep_lanes']}`", f"- content grep hits: `{c['content_grep_hits']}`", f"- P3093 accepted non-imported Ward/current/effective-charge sources: `{c['p3093_accepted_nonimported_ward_current_effective_charge_sources']}`", f"- metric variation rows: `{c['metric_variation_rows']}`", f"- metric variation rows with metric coupling: `{c['metric_variation_rows_with_metric_coupling']}`", f"- metric variation rows with physical stress tensor: `{c['metric_variation_rows_with_physical_stress_tensor']}`", f"- graph energy quadratic rows: `{c['graph_energy_quadratic_rows']}`", f"- graph energy rows with action measure: `{c['graph_energy_rows_with_action_measure']}`", f"- graph energy rows with spacetime metric: `{c['graph_energy_rows_with_spacetime_metric']}`", f"- spectral pressure-like rows: `{c['spectral_pressure_like_rows']}`", f"- spectral pressure rows with physical pressure unit: `{c['spectral_pressure_rows_with_physical_pressure_unit']}`", f"- spectral pressure rows with metric response law: `{c['spectral_pressure_rows_with_metric_response_law']}`", f"- formal stress-divergence rows: `{c['formal_stress_divergence_rows']}`", f"- formal stress rows with covariant conservation law: `{c['formal_stress_rows_with_covariant_conservation_law']}`", f"- formal stress rows with empirical field response: `{c['formal_stress_rows_with_empirical_field_response']}`", f"- stress candidates: `{c['stress_candidates']}`", f"- candidate gate rows: `{c['candidate_gate_rows']}`", f"- accepted non-imported stress-energy/metric-response sources: `{c['accepted_nonimported_stress_energy_metric_response_sources']}`", f"- satisfied proof obligations: `{c['satisfied_proof_obligations']}/{c['proof_obligations']}`", "", "## Decision", payload["decision"]["bounded_result"], "", "## Recommendation", payload["decision"]["next_honest_step"]]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3094/S2044 stress-energy/metric-response obstruction audit", "## P3094/S2044 stress-energy/metric-response obstruction audit\n\n`P3094/S2044` attacks exactly one post-P3093 interface atom: a non-imported stress-energy/metric-response source for the Z12 Dirichlet/Laplacian branch.  It constructs `144` edge-weight metric-variation rows, `12` graph-energy quadratic rows, `7` spectral pressure-like log-det derivative rows, `12` formal stress-divergence rows, and a `5 x 8 = 40` candidate gate matrix.  The finite stress/metric-like algebra remains formal; no physical stress-energy tensor, metric coupling, covariant conservation law, empirical gravitational/field-response interface, spacetime geometry, `L_total`, bridge/role-transfer, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3094/S2044 stress-energy/metric-response source remains unsourced", "## P3094/S2044 stress-energy/metric-response source remains unsourced\n\n`P3094/S2044` confirms that the Z12 Laplacian supports finite edge-weight metric variations, graph-energy quadratics, spectral pressure-like log-det derivatives, and formal stress-divergence rows.  A Lagrangian/EOM, stress tensor, or metric-response reading still needs strict sources for physical metric coupling, covariant conservation, spacetime geometry, empirical gravitational/field-response semantics, and observed-field interpretation; imported continuum stress/gravity templates do not supply those strict sources.\n")
    append_once(AGENTS, "Current stress-energy/metric-response obstruction guardrail (P3094/S2044, 2026-06-25)", "## Current stress-energy/metric-response obstruction guardrail (P3094/S2044, 2026-06-25)\n\n- P3094 follows the P3093 recommendation and audits one standard-physics interface atom: a non-imported stress-energy/metric-response source for the Z12 Dirichlet/Laplacian branch.\n- The finite audit computes `144` edge-weight metric-variation rows, `12` graph-energy quadratic rows, `7` spectral pressure-like log-det derivative rows, `12` formal stress-divergence rows, and `40` candidate gate rows; `0` candidates export an internal non-imported stress-energy/metric-response law.\n- Do not promote edge-weight metric variations, graph-energy quadratics, spectral pressure-like derivatives, formal stress-divergence rows, imported continuum stress-tensor templates, or imported empirical gravity-response templates to physical stress-energy tensor, metric coupling, empirical observable, observed photons/light, spacetime EOM, physical Hamiltonian, selector closure, `L_total`, bridge/role-transfer, or ToE closure.\n- The next honest move is exactly one dispersion/propagating-mode observable obstruction audit for the Z12 Dirichlet/Laplacian branch, unless a genuinely new stress-energy/metric-response theorem is introduced.\n")
    return payload

if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
