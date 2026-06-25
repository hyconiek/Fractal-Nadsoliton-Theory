#!/usr/bin/env python3
"""P3101/S2051: detector/readout calibration obstruction audit.

P3100 found formal open-system-like rows but no non-imported bath/preparation
source.  P3101 attacks exactly one adjacent standard-physics interface atom:
whether finite bath-coupling candidates, preparation-map proxies, clock-orbit
rows, and thermalization-interface witnesses internally source a detector map,
unit-calibrated empirical readout, threshold/noise semantics, or observed-light
interface without importing apparatus templates, observed photons/light,
selector closure, L_total, bridge/role-transfer, or ToE.
"""
from __future__ import annotations

import hashlib, json, math, subprocess
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p3095_s2045_dispersion_propagating_mode_observable_obstruction_audit import Z12_SIZE, eigenvalue
from p3100_s2050_open_system_bath_preparation_source_obstruction_audit import OUT as P3100

OUT = GEN / "p3101_s2051_detector_readout_calibration_obstruction_audit.json"
MD = GEN / "p3101_s2051_detector_readout_calibration_obstruction_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"
CHANNELS = tuple(range(Z12_SIZE))
BETA_GRID = (0.25, 0.5, 1.0, 2.0)
CALIBRATION_SCALES = (0.5, 1.0, 2.0, math.pi, 4.0)
THRESHOLDS = (0.1, 0.25, 0.5, 0.75)
NOISE_LEVELS = (0.0, 0.05, 0.1)

CONTENT_PATTERNS = {
    "detector_atom": r"detector|readout|calibration|empirical observable|apparatus",
    "predecessor_open_system_atom": r"open-system|bath/preparation|physical bath|thermalization interface",
    "blocked_promotions": r"observed light|observed photons|L_total|ToE|selector closure|bridge/role-transfer",
}

DETECTOR_CANDIDATES = (
    {"id": "finite_modal_response_map", "finite_response_map_exported": True, "calibration_orbit_witness_exported": True, "threshold_semantics_exported": False, "noise_model_exported": False, "unit_calibrated_readout_exported": False, "empirical_detector_map_exported": False, "observed_light_interface_exported": False, "nonimported_physical_detector_source_exported": False, "blocker": "modal response is dimensionless and has no detector threshold, noise, units, or empirical apparatus semantics"},
    {"id": "preparation_to_readout_probability_proxy", "finite_response_map_exported": True, "calibration_orbit_witness_exported": True, "threshold_semantics_exported": True, "noise_model_exported": False, "unit_calibrated_readout_exported": False, "empirical_detector_map_exported": False, "observed_light_interface_exported": False, "nonimported_physical_detector_source_exported": False, "blocker": "probability proxy is formal and lacks physical noise, unit calibration, and detector realization"},
    {"id": "clock_scaled_intensity_readout_proxy", "finite_response_map_exported": True, "calibration_orbit_witness_exported": True, "threshold_semantics_exported": True, "noise_model_exported": True, "unit_calibrated_readout_exported": False, "empirical_detector_map_exported": False, "observed_light_interface_exported": False, "nonimported_physical_detector_source_exported": False, "blocker": "clock/scale rows remain an orbit of conventions, not a sourced unit-calibrated readout"},
    {"id": "imported_apparatus_template", "finite_response_map_exported": True, "calibration_orbit_witness_exported": True, "threshold_semantics_exported": True, "noise_model_exported": True, "unit_calibrated_readout_exported": True, "empirical_detector_map_exported": True, "observed_light_interface_exported": False, "nonimported_physical_detector_source_exported": False, "blocker": "detector and calibration semantics pass only by importing apparatus engineering"},
    {"id": "imported_observed_light_detector_template", "finite_response_map_exported": True, "calibration_orbit_witness_exported": True, "threshold_semantics_exported": True, "noise_model_exported": True, "unit_calibrated_readout_exported": True, "empirical_detector_map_exported": True, "observed_light_interface_exported": True, "nonimported_physical_detector_source_exported": False, "blocker": "observed-light interface is imported and not internally sourced by the Z12 branch"},
)
REQUIRED_GATES = ("finite_response_map_exported", "calibration_orbit_witness_exported", "threshold_semantics_exported", "noise_model_exported", "unit_calibrated_readout_exported", "empirical_detector_map_exported", "observed_light_interface_exported", "nonimported_physical_detector_source_exported")


def content_grep() -> list[dict[str, Any]]:
    rows = []
    for lane, pattern in CONTENT_PATTERNS.items():
        proc = subprocess.run(["rg", "-n", pattern, "fundamental_action_reconstruction", "AGENTS.md", "-S"], cwd=REPO, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
        hits = [line for line in proc.stdout.splitlines() if line.strip()]
        rows.append({"content_lane": lane, "pattern": pattern, "hit_count": len(hits), "sample_hits": hits[:8]})
    return rows


def energies() -> list[float]:
    return [round(eigenvalue(k), 12) for k in CHANNELS]


def gibbs(e: list[float], beta: float) -> list[float]:
    weights = [math.exp(-beta * value) for value in e]
    z = sum(weights)
    return [w / z for w in weights]


def detector_response_rows(e: list[float]) -> list[dict[str, Any]]:
    rows = []
    max_e = max(e) or 1.0
    for beta in BETA_GRID:
        pi = gibbs(e, beta)
        for k, value in enumerate(e):
            response = pi[k] * value / max_e
            rows.append({"formal_beta": beta, "channel": k, "energy_label": value, "formal_response_weight": round(response, 12), "finite_response_map_witness": True, "empirical_detector_attached": False, "unit_calibration_attached": False})
    return rows


def calibration_orbit_rows(e: list[float]) -> list[dict[str, Any]]:
    rows = []
    for scale in CALIBRATION_SCALES:
        for k, value in enumerate(e):
            rows.append({"calibration_scale": round(scale, 12), "channel": k, "scaled_readout_label": round(scale * value, 12), "scale_orbit_witness": True, "canonical_unit_fixed": False, "apparatus_calibration_attached": False})
    return rows


def threshold_rows(e: list[float]) -> list[dict[str, Any]]:
    rows = []
    max_e = max(e) or 1.0
    for threshold in THRESHOLDS:
        for k, value in enumerate(e):
            normalized = value / max_e
            rows.append({"threshold": threshold, "channel": k, "normalized_signal": round(normalized, 12), "click_proxy": normalized >= threshold, "threshold_classifier_witness": True, "physical_threshold_source_attached": False})
    return rows


def noise_stability_rows(e: list[float]) -> list[dict[str, Any]]:
    rows = []
    for beta in BETA_GRID:
        pi = gibbs(e, beta)
        mean = sum(p * value for p, value in zip(pi, e))
        for noise in NOISE_LEVELS:
            for k, value in enumerate(e):
                noisy = (1.0 - noise) * value + noise * mean
                rows.append({"formal_beta": beta, "noise_level": noise, "channel": k, "noise_mixed_signal": round(noisy, 12), "noise_stability_witness": True, "physical_noise_model_attached": False, "empirical_readout_attached": False})
    return rows


def gate_rows() -> list[dict[str, Any]]:
    return [{"candidate": c["id"], "required_gate": gate, "gate_passed": bool(c[gate]), "blocking_if_failed": not bool(c[gate]), "detail": "passed" if c[gate] else c["blocker"]} for c in DETECTOR_CANDIDATES for gate in REQUIRED_GATES]


def aggregate(gates: list[dict[str, Any]]) -> list[dict[str, Any]]:
    rows = []
    for c in DETECTOR_CANDIDATES:
        subset = [r for r in gates if r["candidate"] == c["id"]]
        rows.append({"candidate": c["id"], "passed_gates": sum(r["gate_passed"] for r in subset), "failed_gates": sum(not r["gate_passed"] for r in subset), "accepted_nonimported_detector_source": all(r["gate_passed"] for r in subset) and bool(c["nonimported_physical_detector_source_exported"])})
    return rows


def build_payload() -> dict[str, Any]:
    p3100 = read_json(P3100)
    greps = content_grep(); e = energies(); response = detector_response_rows(e); calib = calibration_orbit_rows(e); thresholds = threshold_rows(e); noise = noise_stability_rows(e); gates = gate_rows(); aggs = aggregate(gates)
    accepted = [r for r in aggs if r["accepted_nonimported_detector_source"]]
    obligations = [
        {"obligation": "read_p3100_next_atom", "satisfied": True, "detail": "P3100 selected detector/readout calibration as the next interface atom"},
        {"obligation": "construct_detector_response_rows", "satisfied": len(response) == len(BETA_GRID) * Z12_SIZE, "detail": "finite modal response weights are explicit"},
        {"obligation": "construct_calibration_orbit_rows", "satisfied": len(calib) == len(CALIBRATION_SCALES) * Z12_SIZE and all(not r["canonical_unit_fixed"] for r in calib), "detail": "scale-calibration orbit remains unfixed"},
        {"obligation": "construct_threshold_classifier_rows", "satisfied": len(thresholds) == len(THRESHOLDS) * Z12_SIZE, "detail": "formal click-threshold classifiers are explicit"},
        {"obligation": "construct_noise_stability_rows", "satisfied": len(noise) == len(BETA_GRID) * len(NOISE_LEVELS) * Z12_SIZE, "detail": "formal noise-mixing witnesses are explicit"},
        {"obligation": "export_nonimported_physical_detector_source", "satisfied": False, "detail": "0 candidates pass as internal non-imported detector/readout sources"},
    ]
    return {"status": "P3101_DETECTOR_READOUT_CALIBRATION_OBSTRUCTION_BOUNDED_NO_GO", "input_hashes": {"P3100": hashlib.sha256(P3100.read_bytes()).hexdigest() if P3100.exists() else None}, "constructed_theoretical_objects": {"content_first_repo_grep": greps, "detector_readout_calibration_audit_object": {"object": "Z12DirichletDetectorReadoutCalibrationObstructionAudit", "source_reused": "P3100 recommendation: bounded detector/readout calibration obstruction audit", "required_gates": list(REQUIRED_GATES), "candidate_detector_sources": [c["id"] for c in DETECTOR_CANDIDATES], "acceptance_predicate": "finite response map plus calibration-orbit witness, threshold semantics, noise model, unit-calibrated readout, empirical detector map, observed-light interface, and non-imported physical detector source"}, "detector_response_rows": response, "calibration_orbit_rows": calib, "threshold_classifier_rows": thresholds, "noise_stability_rows": noise, "candidate_gate_rows": gates, "candidate_aggregate_certificate": aggs}, "finite_certificate": {"content_grep_lanes": len(greps), "content_grep_hits": sum(r["hit_count"] for r in greps), "p3100_accepted_nonimported_bath_preparation_sources": p3100.get("finite_certificate", {}).get("accepted_nonimported_bath_preparation_sources"), "detector_response_rows": len(response), "detector_response_rows_with_empirical_detector": sum(r["empirical_detector_attached"] for r in response), "calibration_orbit_rows": len(calib), "calibration_rows_with_canonical_unit": sum(r["canonical_unit_fixed"] for r in calib), "threshold_classifier_rows": len(thresholds), "threshold_rows_with_physical_threshold_source": sum(r["physical_threshold_source_attached"] for r in thresholds), "noise_stability_rows": len(noise), "noise_rows_with_physical_noise_model": sum(r["physical_noise_model_attached"] for r in noise), "detector_candidates": len(DETECTOR_CANDIDATES), "required_gates": len(REQUIRED_GATES), "candidate_gate_rows": len(gates), "accepted_nonimported_detector_sources": len(accepted), "proof_obligations": len(obligations), "satisfied_proof_obligations": sum(r["satisfied"] for r in obligations)}, "proof_obligations": obligations, "decision": {"bounded_result": "P3101 constructs the requested detector/readout calibration obstruction audit.  The Z12 Laplacian supports finite modal response maps, scale-calibration orbits, formal threshold click classifiers, and noise-mixing stability witnesses.  These are real readout-like witnesses, but no internal artifact exports a physical detector map, canonical unit calibration, physical threshold/noise source, observed-light interface, or a non-imported empirical readout source.  Imported apparatus and observed-light templates pass only as imported templates.  Therefore no non-imported detector/readout calibration source is exported.", "negative_export_flags": {key: False for key in ["physical_detector_map_exported", "canonical_unit_calibration_exported", "physical_threshold_source_exported", "physical_noise_model_exported", "empirical_readout_source_exported", "observed_light_interface_exported", "nonimported_physical_detector_source_exported", "physical_bath_preparation_mechanism_exported", "spacetime_eom_exported", "physical_hamiltonian_exported", "QW_2191_discharged", "strict_selector_closure_exported", "ltotal_exported", "bridge_role_transfer_exported", "toe_closure_exported"]}, "positive_scoped_flags": {"detector_response_rows_computed": True, "calibration_orbit_rows_computed": True, "threshold_classifier_rows_computed": True, "noise_stability_rows_computed": True}, "next_honest_step": "Pivot to exactly one adjacent standard-physics interface atom outside selector replay: construct a bounded Born-rule/probability-measure empirical-readout obstruction audit for the Z12 Dirichlet/Laplacian branch, testing whether finite response weights, threshold classifiers, and noise-stability rows supply a non-imported probability measure, event sigma-algebra, detector calibration, and empirical frequency interface without importing quantum measurement postulates, apparatus templates, observed light, L_total, bridge/role-transfer, or ToE."}}


def write_markdown(payload: dict[str, Any]) -> None:
    c = payload["finite_certificate"]
    lines = ["# P3101/S2051 detector/readout calibration obstruction audit", "", f"Status: `{payload['status']}`", "", "## Finite certificate", f"- content grep lanes: `{c['content_grep_lanes']}`", f"- content grep hits: `{c['content_grep_hits']}`", f"- P3100 accepted non-imported bath/preparation sources: `{c['p3100_accepted_nonimported_bath_preparation_sources']}`", f"- detector response rows: `{c['detector_response_rows']}`", f"- detector response rows with empirical detector: `{c['detector_response_rows_with_empirical_detector']}`", f"- calibration orbit rows: `{c['calibration_orbit_rows']}`", f"- calibration rows with canonical unit: `{c['calibration_rows_with_canonical_unit']}`", f"- threshold classifier rows: `{c['threshold_classifier_rows']}`", f"- threshold rows with physical threshold source: `{c['threshold_rows_with_physical_threshold_source']}`", f"- noise stability rows: `{c['noise_stability_rows']}`", f"- noise rows with physical noise model: `{c['noise_rows_with_physical_noise_model']}`", f"- detector candidates: `{c['detector_candidates']}`", f"- candidate gate rows: `{c['candidate_gate_rows']}`", f"- accepted non-imported detector sources: `{c['accepted_nonimported_detector_sources']}`", f"- satisfied proof obligations: `{c['satisfied_proof_obligations']}/{c['proof_obligations']}`", "", "## Decision", payload["decision"]["bounded_result"], "", "## Recommendation", payload["decision"]["next_honest_step"]]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3101/S2051 detector/readout calibration obstruction audit", "## P3101/S2051 detector/readout calibration obstruction audit\n\n`P3101/S2051` attacks exactly one post-P3100 interface atom: a non-imported detector/readout calibration source for the Z12 Dirichlet/Laplacian branch.  It constructs `48` detector response rows, `60` calibration orbit rows, `48` threshold classifier rows, `144` noise stability rows, and a `5 x 8 = 40` candidate gate matrix.  The finite readout-like algebra remains formal; no physical detector map, canonical unit calibration, physical threshold/noise source, observed-light interface, empirical readout source, `L_total`, bridge/role-transfer, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3101/S2051 detector/readout calibration remains unsourced", "## P3101/S2051 detector/readout calibration remains unsourced\n\n`P3101/S2051` confirms that the Z12 Laplacian supports formal detector-response weights, scale-calibration orbits, threshold click classifiers, and noise-mixing witnesses.  A Lagrangian/EOM, measurement, observed-light, or empirical-readout reading still needs strict sources for detector maps, canonical units, physical thresholds/noise, observed-light coupling, and frequency/readout calibration; imported apparatus and measurement templates do not supply those strict sources.\n")
    append_once(AGENTS, "Current detector/readout calibration obstruction guardrail (P3101/S2051, 2026-06-25)", "## Current detector/readout calibration obstruction guardrail (P3101/S2051, 2026-06-25)\n\n- P3101 follows the P3100 recommendation and audits one standard-physics interface atom: a non-imported detector/readout calibration source for the Z12 Dirichlet/Laplacian branch.\n- The finite audit computes `48` detector response rows, `60` calibration orbit rows, `48` threshold classifier rows, `144` noise stability rows, and `40` candidate gate rows; `0` candidates export an internal non-imported detector/readout law.\n- Do not promote finite response maps, scale-calibration orbits, formal threshold classifiers, noise-mixing witnesses, imported apparatus templates, or imported observed-light templates to empirical detector maps, canonical units, observed photons/light, spacetime EOM, physical Hamiltonian, selector closure, `L_total`, bridge/role-transfer, or ToE closure.\n- The next honest move is exactly one Born-rule/probability-measure empirical-readout obstruction audit for the Z12 Dirichlet/Laplacian branch, unless a genuinely new detector/readout source theorem is introduced.\n")
    return payload

if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
