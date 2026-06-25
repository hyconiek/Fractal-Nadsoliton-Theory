#!/usr/bin/env python3
"""P3091/S2041: spectral-action/effective-action obstruction audit.

P3090 left finite Green/resolvent/correlation witnesses but no sourced physical
Green/response law.  P3091 attacks exactly one adjacent standard-physics
interface atom: whether the Z12 Dirichlet/Laplacian branch internally sources a
spectral action or effective-action generating functional without importing QFT
path integrals, spacetime EOM, observed radiation/light, apparatus units,
selector closure, L_total, bridge/role-transfer, or ToE.
"""
from __future__ import annotations

import hashlib, json, math, subprocess
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p3090_s2040_spectral_correlation_green_function_response_obstruction_audit import OUT as P3090, Z12_SIZE

OUT = GEN / "p3091_s2041_spectral_action_effective_action_obstruction_audit.json"
MD = GEN / "p3091_s2041_spectral_action_effective_action_obstruction_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"
MASS2_GRID = (0.25, 1.0, 4.0)
SOURCE_MASKS = (1, 3, 7, 15, 31, 63, 127, 255, 511, 1023, 2047, 4095)
COUPLING_GRID = (0.0, 0.25, 1.0, 4.0)

CONTENT_PATTERNS = {
    "spectral_action_atom": r"spectral-action|effective-action|generating-functional|log-det|quadratic form",
    "predecessor_green_response_atom": r"Green-function|response kernel|resolvent|i-epsilon|two-point",
    "blocked_promotions": r"observed photons|observed light|observed radiation|spacetime EOM|L_total|ToE|selector closure|bridge/role-transfer",
}

ACTION_CANDIDATES = (
    {"id": "finite_z12_logdet_spectral_action", "description": "dimensionless log det(L + m^2 I) spectral functional", "spectral_action_functional_exported": True, "source_coupled_generating_functional_exported": False, "variation_rule_exported": False, "unit_normalized_coupling_exported": False, "empirical_response_generator_exported": False, "nonimported_physical_action_source_exported": False, "blocker": "log-det is finite but not a sourced action, variation principle, unit-coupled generator, or empirical response law"},
    {"id": "source_coupled_quadratic_form_JGJ", "description": "finite source-coupled quadratic form 1/2 J G_m J", "spectral_action_functional_exported": True, "source_coupled_generating_functional_exported": True, "variation_rule_exported": False, "unit_normalized_coupling_exported": False, "empirical_response_generator_exported": False, "nonimported_physical_action_source_exported": False, "blocker": "quadratic generator is formal/static and lacks a strict variational principle, units, and empirical interface"},
    {"id": "p3090_modal_correlation_effective_generator", "description": "formal generator fitted to P3090 modal correlations", "spectral_action_functional_exported": True, "source_coupled_generating_functional_exported": True, "variation_rule_exported": False, "unit_normalized_coupling_exported": False, "empirical_response_generator_exported": False, "nonimported_physical_action_source_exported": False, "blocker": "matching finite correlators does not source the action semantics, coupling normalization, or readout"},
    {"id": "imported_qft_gaussian_path_integral_template", "description": "external Gaussian path-integral effective-action template", "spectral_action_functional_exported": True, "source_coupled_generating_functional_exported": True, "variation_rule_exported": True, "unit_normalized_coupling_exported": False, "empirical_response_generator_exported": False, "nonimported_physical_action_source_exported": False, "blocker": "variation/generator semantics are imported from QFT and remain uncalibrated"},
    {"id": "imported_empirical_effective_field_theory_template", "description": "external unit-bearing EFT/source-response template", "spectral_action_functional_exported": True, "source_coupled_generating_functional_exported": True, "variation_rule_exported": True, "unit_normalized_coupling_exported": True, "empirical_response_generator_exported": True, "nonimported_physical_action_source_exported": False, "blocker": "passes only by importing EFT units, apparatus/readout, and empirical semantics"},
)
REQUIRED_GATES = ("spectral_action_functional_exported", "source_coupled_generating_functional_exported", "variation_rule_exported", "unit_normalized_coupling_exported", "empirical_response_generator_exported", "nonimported_physical_action_source_exported")


def content_grep() -> list[dict[str, Any]]:
    rows = []
    for lane, pattern in CONTENT_PATTERNS.items():
        proc = subprocess.run(["rg", "-n", pattern, "fundamental_action_reconstruction", "AGENTS.md", "-S"], cwd=REPO, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
        hits = [line for line in proc.stdout.splitlines() if line.strip()]
        rows.append({"content_lane": lane, "pattern": pattern, "hit_count": len(hits), "sample_hits": hits[:8]})
    return rows


def lambdas() -> list[float]:
    return [2.0 - 2.0 * math.cos(2.0 * math.pi * k / Z12_SIZE) for k in range(Z12_SIZE)]


def logdet_rows() -> list[dict[str, Any]]:
    vals = lambdas(); rows = []
    for m2 in MASS2_GRID:
        eigs = [lam + m2 for lam in vals]
        rows.append({"dimensionless_mass_squared": m2, "positive_eigenvalues": all(e > 0 for e in eigs), "logdet_L_plus_m2I": round(sum(math.log(e) for e in eigs), 12), "spectral_action_value": round(0.5 * sum(math.log(e) for e in eigs), 12), "action_unit_attached": False, "variation_principle_attached": False})
    return rows


def profile(mask: int) -> list[float]:
    return [float((mask >> i) & 1) for i in range(Z12_SIZE)]


def resolvent_value(distance: int, m2: float) -> float:
    return sum(math.cos(2.0 * math.pi * k * distance / Z12_SIZE) / (lam + m2) for k, lam in enumerate(lambdas())) / Z12_SIZE


def source_quadratic_rows() -> list[dict[str, Any]]:
    rows = []
    for m2 in MASS2_GRID:
        kernel = [resolvent_value(d, m2) for d in range(Z12_SIZE)]
        for mask in SOURCE_MASKS:
            src = profile(mask)
            q = 0.0
            for i, si in enumerate(src):
                for j, sj in enumerate(src):
                    q += 0.5 * si * kernel[(i - j) % Z12_SIZE] * sj
            rows.append({"dimensionless_mass_squared": m2, "source_mask": mask, "support_size": int(sum(src)), "half_JGJ": round(q, 12), "finite_source_coupled_quadratic_witness": True, "unit_normalized_coupling_attached": False, "empirical_response_generator_attached": False})
    return rows


def coupling_scale_orbit_rows() -> list[dict[str, Any]]:
    base = source_quadratic_rows()[0]["half_JGJ"]
    rows = []
    for g in COUPLING_GRID:
        rows.append({"dimensionless_coupling": g, "scaled_generator_value": round(g * g * base, 12), "quadratic_coupling_scaling_witness": True, "coupling_unit_source_attached": False, "absolute_physical_normalization_fixed": False})
    return rows


def finite_variation_rows() -> list[dict[str, Any]]:
    m2 = 1.0; kernel = [resolvent_value(d, m2) for d in range(Z12_SIZE)]; src = profile(31); rows = []
    for i in range(Z12_SIZE):
        grad = sum(kernel[(i - j) % Z12_SIZE] * src[j] for j in range(Z12_SIZE))
        rows.append({"site": i, "formal_delta_half_JGJ_delta_J": round(grad, 12), "finite_formal_variation_witness": True, "eom_source_attached": False, "physical_variation_rule_attached": False})
    return rows


def gate_rows() -> list[dict[str, Any]]:
    return [{"candidate": c["id"], "required_gate": gate, "gate_passed": bool(c[gate]), "blocking_if_failed": not bool(c[gate]), "detail": "passed" if c[gate] else c["blocker"]} for c in ACTION_CANDIDATES for gate in REQUIRED_GATES]


def aggregate(gates: list[dict[str, Any]]) -> list[dict[str, Any]]:
    out = []
    for c in ACTION_CANDIDATES:
        subset = [r for r in gates if r["candidate"] == c["id"]]
        out.append({"candidate": c["id"], "passed_gates": sum(r["gate_passed"] for r in subset), "failed_gates": sum(not r["gate_passed"] for r in subset), "accepted_nonimported_spectral_action_effective_action_source": all(r["gate_passed"] for r in subset) and bool(c["nonimported_physical_action_source_exported"])})
    return out


def build_payload() -> dict[str, Any]:
    p3090 = read_json(P3090)
    greps = content_grep(); logdets = logdet_rows(); quadratics = source_quadratic_rows(); scales = coupling_scale_orbit_rows(); variations = finite_variation_rows(); gates = gate_rows(); aggs = aggregate(gates)
    accepted = [r for r in aggs if r["accepted_nonimported_spectral_action_effective_action_source"]]
    obligations = [
        {"obligation": "read_p3090_next_atom", "satisfied": True, "detail": "P3090 selected spectral-action/effective-action generating-functional audit as the next interface atom"},
        {"obligation": "construct_finite_logdet_spectral_action", "satisfied": len(logdets) == len(MASS2_GRID) and all(r["positive_eigenvalues"] for r in logdets), "detail": "positive mass-regularized log-det rows are finite"},
        {"obligation": "construct_source_coupled_quadratic_generator", "satisfied": len(quadratics) == len(MASS2_GRID) * len(SOURCE_MASKS), "detail": "finite 1/2 JGJ rows are computed over representative source masks"},
        {"obligation": "construct_formal_variation_and_coupling_orbit", "satisfied": len(scales) == len(COUPLING_GRID) and len(variations) == Z12_SIZE, "detail": "formal derivative rows and dimensionless coupling scale orbit are explicit"},
        {"obligation": "export_nonimported_spectral_action_effective_action_source", "satisfied": False, "detail": "0 candidates pass as internal non-imported physical action/generating-functional sources"},
    ]
    return {
        "status": "P3091_SPECTRAL_ACTION_EFFECTIVE_ACTION_GENERATING_FUNCTIONAL_OBSTRUCTION_BOUNDED_NO_GO",
        "input_hashes": {"P3090": hashlib.sha256(P3090.read_bytes()).hexdigest() if P3090.exists() else None},
        "constructed_theoretical_objects": {"content_first_repo_grep": greps, "spectral_action_audit_object": {"object": "Z12DirichletSpectralActionEffectiveActionGeneratingFunctionalObstructionAudit", "source_reused": "P3090 recommendation: bounded spectral-action/effective-action generating-functional obstruction audit", "required_gates": list(REQUIRED_GATES), "candidate_action_sources": [c["id"] for c in ACTION_CANDIDATES], "acceptance_predicate": "spectral action functional plus source-coupled generator, variation rule, unit-normalized coupling, empirical response generator, and non-imported physical action source"}, "mass_regularized_logdet_spectral_action_rows": logdets, "source_coupled_quadratic_generator_rows": quadratics, "coupling_scale_orbit_rows": scales, "finite_formal_variation_rows": variations, "candidate_gate_rows": gates, "candidate_aggregate_certificate": aggs},
        "finite_certificate": {"content_grep_lanes": len(greps), "content_grep_hits": sum(r["hit_count"] for r in greps), "p3090_accepted_nonimported_spectral_correlation_green_response_sources": p3090.get("finite_certificate", {}).get("accepted_nonimported_spectral_correlation_green_response_sources"), "logdet_spectral_action_rows": len(logdets), "logdet_positive_eigenvalue_failures": sum(not r["positive_eigenvalues"] for r in logdets), "logdet_rows_with_action_unit": sum(r["action_unit_attached"] for r in logdets), "source_coupled_quadratic_generator_rows": len(quadratics), "quadratic_rows_with_unit_normalized_coupling": sum(r["unit_normalized_coupling_attached"] for r in quadratics), "quadratic_rows_with_empirical_response_generator": sum(r["empirical_response_generator_attached"] for r in quadratics), "coupling_scale_orbit_rows": len(scales), "coupling_rows_with_unit_source": sum(r["coupling_unit_source_attached"] for r in scales), "finite_formal_variation_rows": len(variations), "variation_rows_with_eom_source": sum(r["eom_source_attached"] for r in variations), "action_candidates": len(ACTION_CANDIDATES), "required_gates": len(REQUIRED_GATES), "candidate_gate_rows": len(gates), "accepted_nonimported_spectral_action_effective_action_sources": len(accepted), "proof_obligations": len(obligations), "satisfied_proof_obligations": sum(r["satisfied"] for r in obligations)},
        "proof_obligations": obligations,
        "decision": {"bounded_result": "P3091 constructs the requested spectral-action/effective-action generating-functional obstruction audit.  The Z12 Laplacian supports finite positive mass-regularized log-det spectral-action values, finite source-coupled 1/2 JGJ generators, a dimensionless coupling scale orbit, and formal finite variation rows.  These are real action-like witnesses, but no internal artifact exports a physical variation principle, action or coupling units, empirical response-generator semantics, spacetime EOM, observed radiation/light, or a non-imported physical action source.  Imported Gaussian path-integral and EFT templates pass only as imported templates.  Therefore no non-imported spectral-action/effective-action source is exported.", "negative_export_flags": {key: False for key in ["variation_rule_exported", "unit_normalized_coupling_exported", "empirical_response_generator_exported", "nonimported_physical_action_source_exported", "physical_action_functional_exported", "empirical_observable_exported", "observed_radiation_exported", "observed_photons_exported", "observed_light_exported", "spacetime_eom_exported", "physical_hamiltonian_exported", "green_response_closure_exported", "born_rule_readout_exported", "QW_2191_discharged", "strict_selector_closure_exported", "ltotal_exported", "bridge_role_transfer_exported", "toe_closure_exported"]}, "positive_scoped_flags": {"finite_logdet_spectral_action_computed": True, "source_coupled_quadratic_generator_computed": True, "dimensionless_coupling_scale_orbit_computed": True, "finite_formal_variation_rows_computed": True}, "next_honest_step": "Pivot to exactly one adjacent standard-physics interface atom outside selector replay: construct a bounded renormalization/scale-flow effective-coupling obstruction audit for the Z12 Dirichlet/Laplacian branch, testing whether log-det scale dependence, mass-parameter flow, and source-coupling rescaling supply a sourced beta function, physical renormalization scale, unit-normalized running coupling, and empirical matching condition without importing continuum QFT RG, spacetime EOM, observed radiation/light, apparatus units, L_total, bridge/role-transfer, or ToE."},
    }


def write_markdown(payload: dict[str, Any]) -> None:
    c = payload["finite_certificate"]
    lines = ["# P3091/S2041 spectral-action/effective-action generating-functional obstruction audit", "", f"Status: `{payload['status']}`", "", "## Finite certificate", f"- content grep lanes: `{c['content_grep_lanes']}`", f"- content grep hits: `{c['content_grep_hits']}`", f"- P3090 accepted non-imported spectral-correlation/Green response sources: `{c['p3090_accepted_nonimported_spectral_correlation_green_response_sources']}`", f"- log-det spectral action rows: `{c['logdet_spectral_action_rows']}`", f"- log-det positive eigenvalue failures: `{c['logdet_positive_eigenvalue_failures']}`", f"- log-det rows with action unit: `{c['logdet_rows_with_action_unit']}`", f"- source-coupled quadratic generator rows: `{c['source_coupled_quadratic_generator_rows']}`", f"- quadratic rows with unit-normalized coupling: `{c['quadratic_rows_with_unit_normalized_coupling']}`", f"- quadratic rows with empirical response generator: `{c['quadratic_rows_with_empirical_response_generator']}`", f"- coupling scale orbit rows: `{c['coupling_scale_orbit_rows']}`", f"- coupling rows with unit source: `{c['coupling_rows_with_unit_source']}`", f"- finite formal variation rows: `{c['finite_formal_variation_rows']}`", f"- variation rows with EOM source: `{c['variation_rows_with_eom_source']}`", f"- action candidates: `{c['action_candidates']}`", f"- candidate gate rows: `{c['candidate_gate_rows']}`", f"- accepted non-imported spectral-action/effective-action sources: `{c['accepted_nonimported_spectral_action_effective_action_sources']}`", f"- satisfied proof obligations: `{c['satisfied_proof_obligations']}/{c['proof_obligations']}`", "", "## Decision", payload["decision"]["bounded_result"], "", "## Recommendation", payload["decision"]["next_honest_step"]]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3091/S2041 spectral-action/effective-action generating-functional obstruction audit", "## P3091/S2041 spectral-action/effective-action generating-functional obstruction audit\n\n`P3091/S2041` attacks exactly one post-P3090 interface atom: a non-imported spectral-action/effective-action generating-functional source for the Z12 Dirichlet/Laplacian branch.  It constructs `3` positive mass-regularized log-det spectral-action rows, `36` source-coupled quadratic-generator rows, `4` coupling-scale orbit rows, `12` finite formal-variation rows, and a `5 x 6 = 30` candidate gate matrix.  The finite action-like algebra remains formal/dimensionless; no physical variation rule, unit-normalized coupling, empirical response generator, spacetime EOM, `L_total`, bridge/role-transfer, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3091/S2041 effective-action source remains unsourced", "## P3091/S2041 effective-action source remains unsourced\n\n`P3091/S2041` confirms that the Z12 Laplacian supports finite log-det, source-coupled quadratic-generator, coupling-scale, and formal-variation witnesses.  A Lagrangian/EOM, effective-action, or empirical-response reading still needs strict sources for the physical variation principle, action/coupling units, empirical response-generator semantics, and spacetime/EOM interpretation; imported path-integral or EFT templates do not supply those strict sources.\n")
    append_once(AGENTS, "Current spectral-action/effective-action obstruction guardrail (P3091/S2041, 2026-06-25)", "## Current spectral-action/effective-action obstruction guardrail (P3091/S2041, 2026-06-25)\n\n- P3091 follows the P3090 recommendation and audits one standard-physics interface atom: a non-imported spectral-action/effective-action generating-functional source for the Z12 Dirichlet/Laplacian branch.\n- The finite audit computes `3` log-det spectral-action rows, `36` source-coupled quadratic-generator rows, `4` coupling-scale orbit rows, `12` finite formal-variation rows, and `30` candidate gate rows; `0` candidates export an internal non-imported physical action/generating-functional source.\n- Do not promote finite log-dets, source-coupled `1/2 JGJ` generators, formal variations, dimensionless coupling-scale orbits, imported Gaussian path-integral templates, or imported EFT templates to empirical observable, observed photons/light, spacetime EOM, physical Hamiltonian, Green/response closure, selector closure, `L_total`, bridge/role-transfer, or ToE closure.\n- The next honest move is exactly one renormalization/scale-flow effective-coupling obstruction audit for the Z12 Dirichlet/Laplacian branch, unless a genuinely new action/generating-functional source theorem is introduced.\n")
    return payload

if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
