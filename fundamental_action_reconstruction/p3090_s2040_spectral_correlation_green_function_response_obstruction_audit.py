#!/usr/bin/env python3
"""P3090/S2040: spectral-correlation/Green-function response obstruction audit.

P3089 left normalized Fourier-power witnesses but no Born-rule/readout source.
P3090 attacks exactly one adjacent standard-physics interface atom: whether the
Z12 Dirichlet/Laplacian branch internally sources a two-point correlator,
Green/response kernel, causal-retarded prescription, unit-calibrated spectral
density, and empirical scattering/readout interface without importing QFT
measurement theory, spacetime EOM, observed radiation/light, apparatus units,
selector closure, L_total, bridge/role-transfer, or ToE.
"""
from __future__ import annotations

import hashlib, json, math, subprocess
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p3089_s2039_spectral_observable_born_rule_readout_obstruction_audit import OUT as P3089, dft_power_probabilities, SOURCE_PROFILE

OUT = GEN / "p3090_s2040_spectral_correlation_green_function_response_obstruction_audit.json"
MD = GEN / "p3090_s2040_spectral_correlation_green_function_response_obstruction_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"
Z12_SIZE = 12
MASS2_GRID = (0.0, 0.25, 1.0, 4.0)
ETA_GRID = (0.01, 0.1, 1.0)

CONTENT_PATTERNS = {
    "green_response_atom": r"Green-function|two-point correlator|response kernel|spectral density|retarded",
    "predecessor_born_rule_atom": r"Born-rule|probability-readout|Fourier-power|measurement basis",
    "blocked_promotions": r"observed photons|observed light|observed radiation|spacetime EOM|L_total|ToE|selector closure|bridge/role-transfer",
}

RESPONSE_CANDIDATES = (
    {"id": "z12_laplacian_pseudoinverse_green_kernel", "description": "Moore-Penrose Green kernel of the finite Z12 Laplacian", "two_point_correlator_exported": True, "response_kernel_exported": True, "causal_retarded_prescription_exported": False, "unit_calibrated_spectral_density_exported": False, "empirical_scattering_readout_exported": False, "nonimported_physical_response_source_exported": False, "blocker": "finite inverse kernel is Euclidean/static and dimensionless; no causal arrow, units, or readout is sourced"},
    {"id": "mass_regularized_z12_resolvent_family", "description": "finite resolvents (L + m^2 I)^-1 over dimensionless mass parameters", "two_point_correlator_exported": True, "response_kernel_exported": True, "causal_retarded_prescription_exported": False, "unit_calibrated_spectral_density_exported": False, "empirical_scattering_readout_exported": False, "nonimported_physical_response_source_exported": False, "blocker": "mass regularization is a formal parameter family, not an internally sourced physical response law"},
    {"id": "p3089_probability_weighted_modal_correlator", "description": "Fourier-power weighted modal two-point correlation witness", "two_point_correlator_exported": True, "response_kernel_exported": False, "causal_retarded_prescription_exported": False, "unit_calibrated_spectral_density_exported": False, "empirical_scattering_readout_exported": False, "nonimported_physical_response_source_exported": False, "blocker": "probability-like weights do not become a response kernel, causal rule, or empirical cross section"},
    {"id": "imported_retarded_green_template", "description": "external i-epsilon/retarded Green-function template", "two_point_correlator_exported": True, "response_kernel_exported": True, "causal_retarded_prescription_exported": True, "unit_calibrated_spectral_density_exported": False, "empirical_scattering_readout_exported": False, "nonimported_physical_response_source_exported": False, "blocker": "causality is imported from external spacetime/QFT structure and lacks internal units/readout"},
    {"id": "imported_scattering_spectral_density_template", "description": "external unit-calibrated spectral-density/scattering readout template", "two_point_correlator_exported": True, "response_kernel_exported": True, "causal_retarded_prescription_exported": True, "unit_calibrated_spectral_density_exported": True, "empirical_scattering_readout_exported": True, "nonimported_physical_response_source_exported": False, "blocker": "passes only by importing QFT/scattering apparatus, units, and empirical semantics"},
)
REQUIRED_GATES = ("two_point_correlator_exported", "response_kernel_exported", "causal_retarded_prescription_exported", "unit_calibrated_spectral_density_exported", "empirical_scattering_readout_exported", "nonimported_physical_response_source_exported")


def content_grep() -> list[dict[str, Any]]:
    rows = []
    for lane, pattern in CONTENT_PATTERNS.items():
        proc = subprocess.run(["rg", "-n", pattern, "fundamental_action_reconstruction", "AGENTS.md", "-S"], cwd=REPO, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
        hits = [line for line in proc.stdout.splitlines() if line.strip()]
        rows.append({"content_lane": lane, "pattern": pattern, "hit_count": len(hits), "sample_hits": hits[:8]})
    return rows


def laplacian_spectrum_rows() -> list[dict[str, Any]]:
    rows = []
    for k in range(Z12_SIZE):
        lam = 2.0 - 2.0 * math.cos(2.0 * math.pi * k / Z12_SIZE)
        rows.append({"mode": k, "laplacian_eigenvalue": round(lam, 12), "nonnegative_spectrum_witness": lam >= -1e-12, "zero_mode": abs(lam) <= 1e-12, "physical_energy_unit_attached": False})
    return rows


def green_kernel_rows() -> list[dict[str, Any]]:
    lambdas = [2.0 - 2.0 * math.cos(2.0 * math.pi * k / Z12_SIZE) for k in range(Z12_SIZE)]
    rows = []
    for distance in range(Z12_SIZE):
        value = sum((0.0 if abs(lam) <= 1e-12 else math.cos(2.0 * math.pi * k * distance / Z12_SIZE) / lam) for k, lam in enumerate(lambdas)) / Z12_SIZE
        rows.append({"cycle_distance": distance, "pseudoinverse_green_value": round(value, 12), "translation_invariant_kernel_witness": True, "static_euclidean_kernel_only": True, "retarded_prescription_attached": False})
    return rows


def resolvent_rows() -> list[dict[str, Any]]:
    lambdas = [2.0 - 2.0 * math.cos(2.0 * math.pi * k / Z12_SIZE) for k in range(Z12_SIZE)]
    rows = []
    for m2 in MASS2_GRID:
        for distance in range(Z12_SIZE):
            if m2 == 0.0:
                value = sum((0.0 if abs(lam) <= 1e-12 else math.cos(2.0 * math.pi * k * distance / Z12_SIZE) / lam) for k, lam in enumerate(lambdas)) / Z12_SIZE
            else:
                value = sum(math.cos(2.0 * math.pi * k * distance / Z12_SIZE) / (lam + m2) for k, lam in enumerate(lambdas)) / Z12_SIZE
            rows.append({"dimensionless_mass_squared": m2, "cycle_distance": distance, "resolvent_value": round(value, 12), "finite_response_kernel_witness": True, "unit_calibrated_spectral_density_attached": False, "empirical_scattering_readout_attached": False})
    return rows


def iepsilon_rows() -> list[dict[str, Any]]:
    lambdas = [2.0 - 2.0 * math.cos(2.0 * math.pi * k / Z12_SIZE) for k in range(Z12_SIZE)]
    rows = []
    for eta in ETA_GRID:
        for k, lam in enumerate(lambdas):
            denom = lam * lam + eta * eta
            rows.append({"eta_regulator": eta, "mode": k, "real_part": round(lam / denom, 12), "imag_part": round(-eta / denom, 12), "formal_bounded_regulator_witness": True, "causal_retarded_source_attached": False, "spacetime_eom_attached": False})
    return rows


def modal_correlation_rows() -> list[dict[str, Any]]:
    probs = dft_power_probabilities(SOURCE_PROFILE)
    rows = []
    for distance in range(Z12_SIZE):
        corr = sum(probs[k] * math.cos(2.0 * math.pi * k * distance / Z12_SIZE) for k in range(Z12_SIZE))
        rows.append({"cycle_distance": distance, "p3089_weighted_modal_correlation": round(corr, 12), "two_point_correlation_witness": True, "response_law_attached": False, "empirical_readout_attached": False})
    return rows


def gate_rows() -> list[dict[str, Any]]:
    return [{"candidate": c["id"], "required_gate": gate, "gate_passed": bool(c[gate]), "blocking_if_failed": not bool(c[gate]), "detail": "passed" if c[gate] else c["blocker"]} for c in RESPONSE_CANDIDATES for gate in REQUIRED_GATES]


def aggregate(gates: list[dict[str, Any]]) -> list[dict[str, Any]]:
    out = []
    for c in RESPONSE_CANDIDATES:
        subset = [r for r in gates if r["candidate"] == c["id"]]
        out.append({"candidate": c["id"], "passed_gates": sum(r["gate_passed"] for r in subset), "failed_gates": sum(not r["gate_passed"] for r in subset), "accepted_nonimported_spectral_correlation_green_response_source": all(r["gate_passed"] for r in subset) and bool(c["nonimported_physical_response_source_exported"])})
    return out


def build_payload() -> dict[str, Any]:
    p3089 = read_json(P3089)
    greps = content_grep(); spectrum = laplacian_spectrum_rows(); green = green_kernel_rows(); resolvents = resolvent_rows(); ieps = iepsilon_rows(); corr = modal_correlation_rows(); gates = gate_rows(); aggs = aggregate(gates)
    accepted = [r for r in aggs if r["accepted_nonimported_spectral_correlation_green_response_source"]]
    obligations = [
        {"obligation": "read_p3089_next_atom", "satisfied": True, "detail": "P3089 selected spectral-correlation/Green-function response as the next interface atom"},
        {"obligation": "construct_finite_laplacian_spectrum_and_green_kernel", "satisfied": len(spectrum) == Z12_SIZE and len(green) == Z12_SIZE, "detail": "finite spectrum and translation-invariant pseudoinverse Green rows constructed"},
        {"obligation": "construct_resolvent_and_iepsilon_witnesses", "satisfied": len(resolvents) == len(MASS2_GRID) * Z12_SIZE and len(ieps) == len(ETA_GRID) * Z12_SIZE, "detail": "formal resolvent and regulator rows are finite and explicit"},
        {"obligation": "reuse_p3089_probability_weights_as_two_point_witness", "satisfied": len(corr) == Z12_SIZE and all(r["two_point_correlation_witness"] for r in corr), "detail": "P3089 weights define a dimensionless modal correlation witness"},
        {"obligation": "export_nonimported_spectral_correlation_green_response_source", "satisfied": False, "detail": "0 candidates pass as internal non-imported physical Green/response sources"},
    ]
    return {
        "status": "P3090_SPECTRAL_CORRELATION_GREEN_FUNCTION_RESPONSE_OBSTRUCTION_BOUNDED_NO_GO",
        "input_hashes": {"P3089": hashlib.sha256(P3089.read_bytes()).hexdigest() if P3089.exists() else None},
        "constructed_theoretical_objects": {"content_first_repo_grep": greps, "green_response_audit_object": {"object": "Z12DirichletSpectralCorrelationGreenFunctionResponseObstructionAudit", "source_reused": "P3089 recommendation: bounded spectral-correlation/Green-function response obstruction audit", "required_gates": list(REQUIRED_GATES), "candidate_response_sources": [c["id"] for c in RESPONSE_CANDIDATES], "acceptance_predicate": "two-point correlator plus response kernel, causal/retarded prescription, unit-calibrated spectral density, empirical scattering/readout, and non-imported physical response source"}, "laplacian_spectrum_rows": spectrum, "pseudoinverse_green_kernel_rows": green, "mass_regularized_resolvent_rows": resolvents, "formal_iepsilon_regulator_rows": ieps, "p3089_weighted_modal_correlation_rows": corr, "candidate_gate_rows": gates, "candidate_aggregate_certificate": aggs},
        "finite_certificate": {"content_grep_lanes": len(greps), "content_grep_hits": sum(r["hit_count"] for r in greps), "p3089_accepted_nonimported_born_rule_probability_readout_sources": p3089.get("finite_certificate", {}).get("accepted_nonimported_born_rule_probability_readout_sources"), "laplacian_spectrum_rows": len(spectrum), "nonnegative_spectrum_failures": sum(not r["nonnegative_spectrum_witness"] for r in spectrum), "zero_mode_rows": sum(r["zero_mode"] for r in spectrum), "pseudoinverse_green_kernel_rows": len(green), "green_rows_with_retarded_prescription": sum(r["retarded_prescription_attached"] for r in green), "mass_regularized_resolvent_rows": len(resolvents), "resolvent_rows_with_unit_calibrated_spectral_density": sum(r["unit_calibrated_spectral_density_attached"] for r in resolvents), "resolvent_rows_with_empirical_scattering_readout": sum(r["empirical_scattering_readout_attached"] for r in resolvents), "formal_iepsilon_regulator_rows": len(ieps), "iepsilon_rows_with_causal_retarded_source": sum(r["causal_retarded_source_attached"] for r in ieps), "p3089_weighted_modal_correlation_rows": len(corr), "modal_correlation_rows_with_response_law": sum(r["response_law_attached"] for r in corr), "response_candidates": len(RESPONSE_CANDIDATES), "required_gates": len(REQUIRED_GATES), "candidate_gate_rows": len(gates), "accepted_nonimported_spectral_correlation_green_response_sources": len(accepted), "proof_obligations": len(obligations), "satisfied_proof_obligations": sum(r["satisfied"] for r in obligations)},
        "proof_obligations": obligations,
        "decision": {"bounded_result": "P3090 constructs the requested spectral-correlation/Green-function response obstruction audit.  The Z12 Laplacian has a finite nonnegative spectrum, a translation-invariant pseudoinverse Green kernel, formal mass-regularized resolvents, formal i-epsilon regulator rows, and P3089-weighted modal two-point correlations.  These are real finite response-like witnesses, but no internal artifact exports a causal/retarded prescription, hbar/action or spectral-density units, empirical scattering/readout semantics, observed radiation/light, or a physical response law.  Imported retarded/scattering templates pass only as imported templates.  Therefore no non-imported spectral-correlation/Green-function response source is exported.", "negative_export_flags": {key: False for key in ["causal_retarded_prescription_exported", "unit_calibrated_spectral_density_exported", "empirical_scattering_readout_exported", "nonimported_physical_response_source_exported", "empirical_observable_exported", "observed_radiation_exported", "observed_photons_exported", "observed_light_exported", "spacetime_eom_exported", "physical_hamiltonian_exported", "born_rule_readout_exported", "QW_2191_discharged", "strict_selector_closure_exported", "ltotal_exported", "bridge_role_transfer_exported", "toe_closure_exported"]}, "positive_scoped_flags": {"finite_laplacian_spectrum_computed": True, "pseudoinverse_green_kernel_computed": True, "mass_regularized_resolvent_family_computed": True, "formal_iepsilon_regulator_rows_computed": True, "p3089_weighted_modal_correlations_computed": True}, "next_honest_step": "Pivot to exactly one adjacent standard-physics interface atom outside selector replay: construct a bounded spectral-action/effective-action generating-functional obstruction audit for the Z12 Dirichlet/Laplacian branch, testing whether the finite determinant/log-det, source-coupled quadratic form, and modal correlators supply a sourced action functional, variation rule, unit-normalized coupling, and empirical response generator without importing QFT path integrals, spacetime EOM, observed radiation/light, apparatus units, L_total, bridge/role-transfer, or ToE."},
    }


def write_markdown(payload: dict[str, Any]) -> None:
    c = payload["finite_certificate"]
    lines = ["# P3090/S2040 spectral-correlation/Green-function response obstruction audit", "", f"Status: `{payload['status']}`", "", "## Finite certificate", f"- content grep lanes: `{c['content_grep_lanes']}`", f"- content grep hits: `{c['content_grep_hits']}`", f"- P3089 accepted non-imported Born-rule/probability-readout sources: `{c['p3089_accepted_nonimported_born_rule_probability_readout_sources']}`", f"- Laplacian spectrum rows: `{c['laplacian_spectrum_rows']}`", f"- nonnegative spectrum failures: `{c['nonnegative_spectrum_failures']}`", f"- zero-mode rows: `{c['zero_mode_rows']}`", f"- pseudoinverse Green kernel rows: `{c['pseudoinverse_green_kernel_rows']}`", f"- Green rows with retarded prescription: `{c['green_rows_with_retarded_prescription']}`", f"- mass-regularized resolvent rows: `{c['mass_regularized_resolvent_rows']}`", f"- resolvent rows with unit-calibrated spectral density: `{c['resolvent_rows_with_unit_calibrated_spectral_density']}`", f"- resolvent rows with empirical scattering readout: `{c['resolvent_rows_with_empirical_scattering_readout']}`", f"- formal i-epsilon regulator rows: `{c['formal_iepsilon_regulator_rows']}`", f"- i-epsilon rows with causal/retarded source: `{c['iepsilon_rows_with_causal_retarded_source']}`", f"- P3089-weighted modal correlation rows: `{c['p3089_weighted_modal_correlation_rows']}`", f"- modal correlation rows with response law: `{c['modal_correlation_rows_with_response_law']}`", f"- response candidates: `{c['response_candidates']}`", f"- candidate gate rows: `{c['candidate_gate_rows']}`", f"- accepted non-imported spectral-correlation/Green response sources: `{c['accepted_nonimported_spectral_correlation_green_response_sources']}`", f"- satisfied proof obligations: `{c['satisfied_proof_obligations']}/{c['proof_obligations']}`", "", "## Decision", payload["decision"]["bounded_result"], "", "## Recommendation", payload["decision"]["next_honest_step"]]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3090/S2040 spectral-correlation/Green-function response obstruction audit", "## P3090/S2040 spectral-correlation/Green-function response obstruction audit\n\n`P3090/S2040` attacks exactly one post-P3089 interface atom: a non-imported spectral-correlation/Green-function response source for the Z12 Dirichlet/Laplacian branch.  It constructs `12` finite Laplacian spectrum rows, `12` pseudoinverse Green kernel rows, `48` mass-regularized resolvent rows, `36` formal i-epsilon regulator rows, `12` P3089-weighted modal correlation rows, and a `5 x 6 = 30` candidate gate matrix.  The finite response-like algebra remains static/dimensionless; no causal/retarded prescription, unit-calibrated spectral density, empirical scattering/readout, observed radiation/light, spacetime EOM, `L_total`, bridge/role-transfer, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3090/S2040 Green-function/response source remains unsourced", "## P3090/S2040 Green-function/response source remains unsourced\n\n`P3090/S2040` confirms that the Z12 Laplacian supports formal finite Green/resolvent/correlation witnesses.  A Lagrangian/EOM, response-theory, or empirical-scattering reading still needs strict sources for causal/retarded structure, action/energy units, spectral-density calibration, and empirical readout semantics; imported QFT Green-function or scattering templates do not supply those strict sources.\n")
    append_once(AGENTS, "Current spectral-correlation/Green-function response obstruction guardrail (P3090/S2040, 2026-06-25)", "## Current spectral-correlation/Green-function response obstruction guardrail (P3090/S2040, 2026-06-25)\n\n- P3090 follows the P3089 recommendation and audits one standard-physics interface atom: a non-imported spectral-correlation/Green-function response source for the Z12 Dirichlet/Laplacian branch.\n- The finite audit computes `12` spectrum rows, `12` pseudoinverse Green rows, `48` mass-regularized resolvent rows, `36` formal i-epsilon regulator rows, `12` P3089-weighted modal correlation rows, and `30` candidate gate rows; `0` candidates export an internal non-imported physical Green/response source.\n- Do not promote finite Green kernels, formal resolvents, i-epsilon regulators, P3089-weighted modal correlations, imported retarded Green templates, or imported scattering/spectral-density templates to empirical observable, observed photons/light, spacetime EOM, physical Hamiltonian, Born-rule/readout closure, selector closure, `L_total`, bridge/role-transfer, or ToE closure.\n- The next honest move is exactly one spectral-action/effective-action generating-functional obstruction audit for the Z12 Dirichlet/Laplacian branch, unless a genuinely new Green/response source theorem is introduced.\n")
    return payload

if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
