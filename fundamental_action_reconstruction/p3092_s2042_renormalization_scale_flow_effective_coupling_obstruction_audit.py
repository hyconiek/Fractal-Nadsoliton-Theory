#!/usr/bin/env python3
"""P3092/S2042: renormalization/scale-flow effective-coupling audit.

P3091 left finite spectral-action/effective-action-like witnesses but no sourced
physical action/generating-functional law.  P3092 attacks exactly one adjacent
standard-physics interface atom: whether the Z12 Dirichlet/Laplacian branch
internally sources a renormalization/scale-flow effective-coupling structure
without importing continuum QFT RG, spacetime EOM, observed radiation/light,
apparatus units, selector closure, L_total, bridge/role-transfer, or ToE.
"""
from __future__ import annotations

import hashlib, json, math, subprocess
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p3091_s2041_spectral_action_effective_action_obstruction_audit import OUT as P3091, Z12_SIZE, lambdas

OUT = GEN / "p3092_s2042_renormalization_scale_flow_effective_coupling_obstruction_audit.json"
MD = GEN / "p3092_s2042_renormalization_scale_flow_effective_coupling_obstruction_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"
MASS2_GRID = (0.125, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0)
COUPLING_GRID = (0.25, 0.5, 1.0, 2.0)
REFERENCE_DISTANCE = 1

CONTENT_PATTERNS = {
    "renormalization_scale_flow_atom": r"renormalization|scale-flow|running coupling|beta function|effective-coupling",
    "predecessor_action_atom": r"spectral-action|effective-action|generating-functional|log-det|1/2 JGJ",
    "blocked_promotions": r"continuum QFT RG|observed photons|observed light|spacetime EOM|L_total|ToE|selector closure|bridge/role-transfer",
}

FLOW_CANDIDATES = (
    {"id": "finite_logdet_mass_scale_derivative", "description": "dimensionless derivative of log det(L+m^2 I) along the mass parameter", "logdet_scale_dependence_exported": True, "mass_parameter_flow_exported": True, "sourced_beta_function_exported": False, "physical_renormalization_scale_exported": False, "unit_normalized_running_coupling_exported": False, "empirical_matching_condition_exported": False, "nonimported_physical_rg_source_exported": False, "blocker": "finite mass derivatives do not source a physical RG scale, beta function, units, or matching condition"},
    {"id": "green_ratio_effective_coupling_orbit", "description": "effective coupling inferred from ratios of mass-regularized Green kernels", "logdet_scale_dependence_exported": False, "mass_parameter_flow_exported": True, "sourced_beta_function_exported": False, "physical_renormalization_scale_exported": False, "unit_normalized_running_coupling_exported": False, "empirical_matching_condition_exported": False, "nonimported_physical_rg_source_exported": False, "blocker": "ratio flow is formal and target-dependent; it lacks a sourced RG law and physical normalization"},
    {"id": "source_coupling_rescaling_profile", "description": "dimensionless rescaling orbit for source coupling g over finite 1/2 JGJ witnesses", "logdet_scale_dependence_exported": False, "mass_parameter_flow_exported": True, "sourced_beta_function_exported": False, "physical_renormalization_scale_exported": False, "unit_normalized_running_coupling_exported": False, "empirical_matching_condition_exported": False, "nonimported_physical_rg_source_exported": False, "blocker": "coupling rescaling is a chosen orbit rather than an internally sourced running law"},
    {"id": "imported_continuum_qft_rg_template", "description": "external Callan-Symanzik/continuum QFT RG template", "logdet_scale_dependence_exported": True, "mass_parameter_flow_exported": True, "sourced_beta_function_exported": True, "physical_renormalization_scale_exported": True, "unit_normalized_running_coupling_exported": False, "empirical_matching_condition_exported": False, "nonimported_physical_rg_source_exported": False, "blocker": "RG semantics are imported and still lack internal units/matching"},
    {"id": "imported_empirical_running_coupling_template", "description": "external experimentally matched running-coupling template", "logdet_scale_dependence_exported": True, "mass_parameter_flow_exported": True, "sourced_beta_function_exported": True, "physical_renormalization_scale_exported": True, "unit_normalized_running_coupling_exported": True, "empirical_matching_condition_exported": True, "nonimported_physical_rg_source_exported": False, "blocker": "passes only by importing physical scale units and empirical matching semantics"},
)
REQUIRED_GATES = ("logdet_scale_dependence_exported", "mass_parameter_flow_exported", "sourced_beta_function_exported", "physical_renormalization_scale_exported", "unit_normalized_running_coupling_exported", "empirical_matching_condition_exported", "nonimported_physical_rg_source_exported")


def content_grep() -> list[dict[str, Any]]:
    rows = []
    for lane, pattern in CONTENT_PATTERNS.items():
        proc = subprocess.run(["rg", "-n", pattern, "fundamental_action_reconstruction", "AGENTS.md", "-S"], cwd=REPO, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
        hits = [line for line in proc.stdout.splitlines() if line.strip()]
        rows.append({"content_lane": lane, "pattern": pattern, "hit_count": len(hits), "sample_hits": hits[:8]})
    return rows


def logdet(m2: float) -> float:
    return sum(math.log(lam + m2) for lam in lambdas())


def trace_resolvent(m2: float) -> float:
    return sum(1.0 / (lam + m2) for lam in lambdas())


def green_value(distance: int, m2: float) -> float:
    return sum(math.cos(2.0 * math.pi * k * distance / Z12_SIZE) / (lam + m2) for k, lam in enumerate(lambdas())) / Z12_SIZE


def scale_rows() -> list[dict[str, Any]]:
    rows = []
    for m2 in MASS2_GRID:
        rows.append({"dimensionless_mass_squared": m2, "dimensionless_log_mass": round(math.log(m2), 12), "logdet_L_plus_m2I": round(logdet(m2), 12), "d_logdet_d_log_m2_trace_witness": round(m2 * trace_resolvent(m2), 12), "finite_scale_dependence_witness": True, "physical_rg_scale_attached": False, "action_unit_attached": False})
    return rows


def beta_estimate_rows() -> list[dict[str, Any]]:
    rows = []
    for lo, hi in zip(MASS2_GRID, MASS2_GRID[1:]):
        g_lo = green_value(REFERENCE_DISTANCE, lo)
        g_hi = green_value(REFERENCE_DISTANCE, hi)
        beta = (g_hi - g_lo) / (math.log(hi) - math.log(lo))
        rows.append({"from_mass_squared": lo, "to_mass_squared": hi, "reference_distance": REFERENCE_DISTANCE, "green_effective_coupling_lo": round(g_lo, 12), "green_effective_coupling_hi": round(g_hi, 12), "finite_difference_beta_like_slope": round(beta, 12), "beta_like_witness_computed": True, "sourced_beta_function_attached": False, "empirical_matching_condition_attached": False})
    return rows


def coupling_rescaling_rows() -> list[dict[str, Any]]:
    base = green_value(REFERENCE_DISTANCE, 1.0)
    rows = []
    for m2 in MASS2_GRID:
        for g0 in COUPLING_GRID:
            formal_running = g0 / math.sqrt(1.0 + m2)
            rows.append({"dimensionless_mass_squared": m2, "bare_dimensionless_coupling": g0, "formal_rescaled_coupling": round(formal_running, 12), "formal_response_weight": round((formal_running ** 2) * base, 12), "coupling_rescaling_witness": True, "unit_normalized_running_coupling_attached": False, "physical_normalization_fixed": False})
    return rows


def gate_rows() -> list[dict[str, Any]]:
    return [{"candidate": c["id"], "required_gate": gate, "gate_passed": bool(c[gate]), "blocking_if_failed": not bool(c[gate]), "detail": "passed" if c[gate] else c["blocker"]} for c in FLOW_CANDIDATES for gate in REQUIRED_GATES]


def aggregate(gates: list[dict[str, Any]]) -> list[dict[str, Any]]:
    out = []
    for c in FLOW_CANDIDATES:
        subset = [r for r in gates if r["candidate"] == c["id"]]
        out.append({"candidate": c["id"], "passed_gates": sum(r["gate_passed"] for r in subset), "failed_gates": sum(not r["gate_passed"] for r in subset), "accepted_nonimported_renormalization_scale_flow_source": all(r["gate_passed"] for r in subset) and bool(c["nonimported_physical_rg_source_exported"])})
    return out


def build_payload() -> dict[str, Any]:
    p3091 = read_json(P3091)
    greps = content_grep(); scales = scale_rows(); betas = beta_estimate_rows(); rescalings = coupling_rescaling_rows(); gates = gate_rows(); aggs = aggregate(gates)
    accepted = [r for r in aggs if r["accepted_nonimported_renormalization_scale_flow_source"]]
    obligations = [
        {"obligation": "read_p3091_next_atom", "satisfied": True, "detail": "P3091 selected renormalization/scale-flow effective-coupling audit as the next interface atom"},
        {"obligation": "construct_logdet_scale_dependence", "satisfied": len(scales) == len(MASS2_GRID) and all(r["finite_scale_dependence_witness"] for r in scales), "detail": "finite log-det scale rows and trace derivative witnesses are explicit"},
        {"obligation": "construct_beta_like_mass_flow", "satisfied": len(betas) == len(MASS2_GRID) - 1 and all(r["beta_like_witness_computed"] for r in betas), "detail": "finite-difference beta-like slopes over Green effective couplings are computed"},
        {"obligation": "construct_coupling_rescaling_orbit", "satisfied": len(rescalings) == len(MASS2_GRID) * len(COUPLING_GRID), "detail": "dimensionless source-coupling rescaling rows are explicit"},
        {"obligation": "export_nonimported_physical_rg_source", "satisfied": False, "detail": "0 candidates pass as internal non-imported renormalization/scale-flow effective-coupling sources"},
    ]
    return {
        "status": "P3092_RENORMALIZATION_SCALE_FLOW_EFFECTIVE_COUPLING_OBSTRUCTION_BOUNDED_NO_GO",
        "input_hashes": {"P3091": hashlib.sha256(P3091.read_bytes()).hexdigest() if P3091.exists() else None},
        "constructed_theoretical_objects": {"content_first_repo_grep": greps, "renormalization_scale_flow_audit_object": {"object": "Z12DirichletRenormalizationScaleFlowEffectiveCouplingObstructionAudit", "source_reused": "P3091 recommendation: bounded renormalization/scale-flow effective-coupling obstruction audit", "required_gates": list(REQUIRED_GATES), "candidate_rg_sources": [c["id"] for c in FLOW_CANDIDATES], "acceptance_predicate": "finite scale dependence plus mass flow, sourced beta function, physical renormalization scale, unit-normalized running coupling, empirical matching condition, and non-imported physical RG source"}, "logdet_scale_dependence_rows": scales, "green_effective_coupling_beta_like_rows": betas, "coupling_rescaling_orbit_rows": rescalings, "candidate_gate_rows": gates, "candidate_aggregate_certificate": aggs},
        "finite_certificate": {"content_grep_lanes": len(greps), "content_grep_hits": sum(r["hit_count"] for r in greps), "p3091_accepted_nonimported_spectral_action_effective_action_sources": p3091.get("finite_certificate", {}).get("accepted_nonimported_spectral_action_effective_action_sources"), "logdet_scale_dependence_rows": len(scales), "logdet_rows_with_physical_rg_scale": sum(r["physical_rg_scale_attached"] for r in scales), "logdet_rows_with_action_unit": sum(r["action_unit_attached"] for r in scales), "green_beta_like_rows": len(betas), "green_beta_like_rows_with_sourced_beta_function": sum(r["sourced_beta_function_attached"] for r in betas), "green_beta_like_rows_with_empirical_matching": sum(r["empirical_matching_condition_attached"] for r in betas), "coupling_rescaling_orbit_rows": len(rescalings), "coupling_rows_with_unit_normalized_running_coupling": sum(r["unit_normalized_running_coupling_attached"] for r in rescalings), "coupling_rows_with_physical_normalization": sum(r["physical_normalization_fixed"] for r in rescalings), "rg_candidates": len(FLOW_CANDIDATES), "required_gates": len(REQUIRED_GATES), "candidate_gate_rows": len(gates), "accepted_nonimported_renormalization_scale_flow_sources": len(accepted), "proof_obligations": len(obligations), "satisfied_proof_obligations": sum(r["satisfied"] for r in obligations)},
        "proof_obligations": obligations,
        "decision": {"bounded_result": "P3092 constructs the requested renormalization/scale-flow effective-coupling obstruction audit.  The Z12 Laplacian supports finite log-det scale-dependence rows, trace derivative witnesses, Green-kernel beta-like finite differences, and dimensionless source-coupling rescaling orbits.  These are real flow-like witnesses, but no internal artifact exports a sourced beta function, physical renormalization scale, action/coupling units, empirical matching condition, spacetime EOM, observed radiation/light, or a non-imported physical RG source.  Imported continuum QFT RG and empirical running-coupling templates pass only as imported templates.  Therefore no non-imported renormalization/scale-flow effective-coupling source is exported.", "negative_export_flags": {key: False for key in ["sourced_beta_function_exported", "physical_renormalization_scale_exported", "unit_normalized_running_coupling_exported", "empirical_matching_condition_exported", "nonimported_physical_rg_source_exported", "physical_action_functional_exported", "empirical_observable_exported", "observed_radiation_exported", "observed_photons_exported", "observed_light_exported", "spacetime_eom_exported", "physical_hamiltonian_exported", "green_response_closure_exported", "effective_action_closure_exported", "born_rule_readout_exported", "QW_2191_discharged", "strict_selector_closure_exported", "ltotal_exported", "bridge_role_transfer_exported", "toe_closure_exported"]}, "positive_scoped_flags": {"finite_logdet_scale_dependence_computed": True, "trace_derivative_witness_computed": True, "green_beta_like_finite_differences_computed": True, "dimensionless_coupling_rescaling_orbit_computed": True}, "next_honest_step": "Pivot to exactly one adjacent standard-physics interface atom outside selector replay: construct a bounded Ward-identity/symmetry-current effective-charge obstruction audit for the Z12 Dirichlet/Laplacian branch, testing whether finite graph symmetries, Laplacian commutators, source variations, and spectral charges supply a non-imported conserved current, Ward identity, gauge-charge normalization, and empirical charge/readout interface without importing continuum gauge theory, spacetime EOM, observed photons/light, apparatus units, L_total, bridge/role-transfer, or ToE."},
    }


def write_markdown(payload: dict[str, Any]) -> None:
    c = payload["finite_certificate"]
    lines = ["# P3092/S2042 renormalization/scale-flow effective-coupling obstruction audit", "", f"Status: `{payload['status']}`", "", "## Finite certificate", f"- content grep lanes: `{c['content_grep_lanes']}`", f"- content grep hits: `{c['content_grep_hits']}`", f"- P3091 accepted non-imported spectral-action/effective-action sources: `{c['p3091_accepted_nonimported_spectral_action_effective_action_sources']}`", f"- log-det scale-dependence rows: `{c['logdet_scale_dependence_rows']}`", f"- log-det rows with physical RG scale: `{c['logdet_rows_with_physical_rg_scale']}`", f"- log-det rows with action unit: `{c['logdet_rows_with_action_unit']}`", f"- Green beta-like rows: `{c['green_beta_like_rows']}`", f"- Green beta-like rows with sourced beta function: `{c['green_beta_like_rows_with_sourced_beta_function']}`", f"- Green beta-like rows with empirical matching: `{c['green_beta_like_rows_with_empirical_matching']}`", f"- coupling rescaling orbit rows: `{c['coupling_rescaling_orbit_rows']}`", f"- coupling rows with unit-normalized running coupling: `{c['coupling_rows_with_unit_normalized_running_coupling']}`", f"- coupling rows with physical normalization: `{c['coupling_rows_with_physical_normalization']}`", f"- RG candidates: `{c['rg_candidates']}`", f"- candidate gate rows: `{c['candidate_gate_rows']}`", f"- accepted non-imported renormalization/scale-flow sources: `{c['accepted_nonimported_renormalization_scale_flow_sources']}`", f"- satisfied proof obligations: `{c['satisfied_proof_obligations']}/{c['proof_obligations']}`", "", "## Decision", payload["decision"]["bounded_result"], "", "## Recommendation", payload["decision"]["next_honest_step"]]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3092/S2042 renormalization/scale-flow effective-coupling obstruction audit", "## P3092/S2042 renormalization/scale-flow effective-coupling obstruction audit\n\n`P3092/S2042` attacks exactly one post-P3091 interface atom: a non-imported renormalization/scale-flow effective-coupling source for the Z12 Dirichlet/Laplacian branch.  It constructs `7` log-det scale-dependence rows, `6` Green-kernel beta-like finite-difference rows, `28` coupling-rescaling orbit rows, and a `5 x 7 = 35` candidate gate matrix.  The finite flow-like algebra remains formal/dimensionless; no sourced beta function, physical renormalization scale, unit-normalized running coupling, empirical matching condition, spacetime EOM, `L_total`, bridge/role-transfer, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3092/S2042 RG scale-flow source remains unsourced", "## P3092/S2042 RG scale-flow source remains unsourced\n\n`P3092/S2042` confirms that the Z12 Laplacian supports finite log-det scale dependence, trace derivative witnesses, Green-kernel beta-like slopes, and source-coupling rescaling orbits.  A Lagrangian/EOM, effective coupling, or empirical running-coupling reading still needs strict sources for a physical RG scale, sourced beta function, unit-normalized coupling, empirical matching, and spacetime/EOM interpretation; imported continuum RG templates do not supply those strict sources.\n")
    append_once(AGENTS, "Current renormalization/scale-flow effective-coupling obstruction guardrail (P3092/S2042, 2026-06-25)", "## Current renormalization/scale-flow effective-coupling obstruction guardrail (P3092/S2042, 2026-06-25)\n\n- P3092 follows the P3091 recommendation and audits one standard-physics interface atom: a non-imported renormalization/scale-flow effective-coupling source for the Z12 Dirichlet/Laplacian branch.\n- The finite audit computes `7` log-det scale-dependence rows, `6` Green-kernel beta-like finite-difference rows, `28` coupling-rescaling orbit rows, and `35` candidate gate rows; `0` candidates export an internal non-imported physical RG/source-flow law.\n- Do not promote finite log-det scale derivatives, Green-kernel beta-like slopes, dimensionless source-coupling rescaling orbits, imported continuum QFT RG templates, or imported empirical running-coupling templates to sourced beta function, physical renormalization scale, unit-normalized running coupling, empirical observable, observed photons/light, spacetime EOM, physical Hamiltonian, selector closure, `L_total`, bridge/role-transfer, or ToE closure.\n- The next honest move is exactly one Ward-identity/symmetry-current effective-charge obstruction audit for the Z12 Dirichlet/Laplacian branch, unless a genuinely new RG/source-flow theorem is introduced.\n")
    return payload

if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
