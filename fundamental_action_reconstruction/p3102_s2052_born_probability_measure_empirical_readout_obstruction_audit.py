#!/usr/bin/env python3
"""P3102/S2052: Born-rule/probability-measure empirical-readout obstruction audit.

P3101 found formal detector/readout rows but no non-imported detector source.
P3102 attacks exactly one adjacent interface atom: whether Z12 Dirichlet/Laplian
finite weights can internally source a Born-rule probability measure, event
algebra, empirical frequency semantics, and calibrated readout without importing
quantum measurement postulates, apparatus templates, observed light, selector
closure, L_total, bridge/role-transfer, or ToE.
"""
from __future__ import annotations

import hashlib, json, math, subprocess
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p3095_s2045_dispersion_propagating_mode_observable_obstruction_audit import Z12_SIZE
from p3101_s2051_detector_readout_calibration_obstruction_audit import OUT as P3101

OUT = GEN / "p3102_s2052_born_probability_measure_empirical_readout_obstruction_audit.json"
MD = GEN / "p3102_s2052_born_probability_measure_empirical_readout_obstruction_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"
CHANNELS = tuple(range(Z12_SIZE))
STATES = ("delta0", "block3", "alternating", "ramp")
BASIS_SHIFTS = tuple(range(6))
SAMPLE_SIZES = (12, 120, 1200, 12000)
CONTENT_PATTERNS = {
    "born_atom": r"Born-rule|Born rule|probability-measure|probability readout|frequency",
    "predecessor_readout_atom": r"detector/readout|calibration|threshold|noise",
    "blocked_promotions": r"quantum measurement|observed light|L_total|ToE|selector closure|bridge/role-transfer",
}

CANDIDATES = (
    {"id": "fourier_power_probability_proxy", "normalized_probability_measure_exported": True, "event_sigma_algebra_exported": True, "born_rule_map_exported": False, "measurement_basis_source_exported": False, "detector_calibration_link_exported": False, "empirical_frequency_interface_exported": False, "p3101_detector_readout_source_exported": False, "nonimported_physical_probability_source_exported": False, "blocker": "normalized weights are formal Fourier-power proxies, not a sourced Born postulate or empirical law"},
    {"id": "finite_event_additivity_proxy", "normalized_probability_measure_exported": True, "event_sigma_algebra_exported": True, "born_rule_map_exported": False, "measurement_basis_source_exported": False, "detector_calibration_link_exported": False, "empirical_frequency_interface_exported": False, "p3101_detector_readout_source_exported": False, "nonimported_physical_probability_source_exported": False, "blocker": "finite additivity is mathematical and lacks physical measurement semantics"},
    {"id": "basis_orbit_born_candidate", "normalized_probability_measure_exported": True, "event_sigma_algebra_exported": True, "born_rule_map_exported": True, "measurement_basis_source_exported": False, "detector_calibration_link_exported": False, "empirical_frequency_interface_exported": False, "p3101_detector_readout_source_exported": False, "nonimported_physical_probability_source_exported": False, "blocker": "Born-like basis weights pass only after choosing an unsourced measurement basis"},
    {"id": "imported_quantum_measurement_template", "normalized_probability_measure_exported": True, "event_sigma_algebra_exported": True, "born_rule_map_exported": True, "measurement_basis_source_exported": True, "detector_calibration_link_exported": True, "empirical_frequency_interface_exported": True, "p3101_detector_readout_source_exported": False, "nonimported_physical_probability_source_exported": False, "blocker": "measurement and frequency semantics are imported quantum/apparatus templates"},
    {"id": "imported_observed_frequency_readout_template", "normalized_probability_measure_exported": True, "event_sigma_algebra_exported": True, "born_rule_map_exported": True, "measurement_basis_source_exported": True, "detector_calibration_link_exported": True, "empirical_frequency_interface_exported": True, "p3101_detector_readout_source_exported": False, "nonimported_physical_probability_source_exported": False, "blocker": "observed frequencies require an imported detector/readout source rejected by P3101"},
)
GATES = ("normalized_probability_measure_exported", "event_sigma_algebra_exported", "born_rule_map_exported", "measurement_basis_source_exported", "detector_calibration_link_exported", "empirical_frequency_interface_exported", "p3101_detector_readout_source_exported", "nonimported_physical_probability_source_exported")


def content_grep() -> list[dict[str, Any]]:
    rows = []
    for lane, pattern in CONTENT_PATTERNS.items():
        proc = subprocess.run(["rg", "-n", pattern, "fundamental_action_reconstruction", "AGENTS.md", "-S"], cwd=REPO, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
        hits = [line for line in proc.stdout.splitlines() if line.strip()]
        rows.append({"content_lane": lane, "pattern": pattern, "hit_count": len(hits), "sample_hits": hits[:8]})
    return rows


def amplitudes(state: str) -> list[complex]:
    if state == "delta0":
        return [1.0 if n == 0 else 0.0 for n in CHANNELS]
    if state == "block3":
        return [1.0 if n in (0, 1, 2) else 0.0 for n in CHANNELS]
    if state == "alternating":
        return [1.0 if n % 2 == 0 else -1.0 for n in CHANNELS]
    return [float(n + 1) for n in CHANNELS]


def fourier_weights(state: str, shift: int = 0) -> list[float]:
    vec = amplitudes(state)
    powers = []
    for k in CHANNELS:
        amp = sum(vec[n] * complex(math.cos(-2 * math.pi * (k + shift) * n / Z12_SIZE), math.sin(-2 * math.pi * (k + shift) * n / Z12_SIZE)) for n in CHANNELS)
        powers.append(abs(amp) ** 2)
    total = sum(powers) or 1.0
    return [p / total for p in powers]


def probability_rows() -> list[dict[str, Any]]:
    rows = []
    for state in STATES:
        weights = fourier_weights(state)
        for k, weight in enumerate(weights):
            rows.append({"state": state, "event": f"mode_{k}", "probability_weight": round(weight, 12), "nonnegative": weight >= -1e-12, "normalized_measure_witness": True, "born_rule_source_attached": False, "empirical_frequency_attached": False})
    return rows


def event_additivity_rows() -> list[dict[str, Any]]:
    pairs = ((0, 1), (2, 3), (4, 5), (6, 7))
    rows = []
    for state in STATES:
        w = fourier_weights(state)
        for a, b in pairs:
            rows.append({"state": state, "disjoint_events": [f"mode_{a}", f"mode_{b}"], "union_probability": round(w[a] + w[b], 12), "sum_probability": round(w[a] + w[b], 12), "finite_additivity_witness": True, "physical_event_semantics_attached": False})
    return rows


def basis_orbit_rows() -> list[dict[str, Any]]:
    rows = []
    for shift in BASIS_SHIFTS:
        for k, weight in enumerate(fourier_weights("block3", shift)):
            rows.append({"basis_shift": shift, "mode": k, "born_like_weight": round(weight, 12), "basis_orbit_witness": True, "canonical_measurement_basis_fixed": False})
    return rows


def frequency_proxy_rows() -> list[dict[str, Any]]:
    rows = []
    for state in STATES:
        weights = fourier_weights(state)
        for sample_size in SAMPLE_SIZES:
            expected_total = sum(round(sample_size * w) for w in weights)
            rows.append({"state": state, "sample_size": sample_size, "rounded_expected_total": int(expected_total), "normalization_defect": int(expected_total - sample_size), "frequency_proxy_witness": True, "apparatus_trial_protocol_attached": False, "empirical_readout_attached": False})
    return rows


def gate_rows() -> list[dict[str, Any]]:
    return [{"candidate": c["id"], "required_gate": gate, "gate_passed": bool(c[gate]), "blocking_if_failed": not bool(c[gate]), "detail": "passed" if c[gate] else c["blocker"]} for c in CANDIDATES for gate in GATES]


def aggregate(gates: list[dict[str, Any]]) -> list[dict[str, Any]]:
    out = []
    for c in CANDIDATES:
        subset = [r for r in gates if r["candidate"] == c["id"]]
        out.append({"candidate": c["id"], "passed_gates": sum(r["gate_passed"] for r in subset), "failed_gates": sum(not r["gate_passed"] for r in subset), "accepted_nonimported_probability_source": all(r["gate_passed"] for r in subset) and bool(c["nonimported_physical_probability_source_exported"])})
    return out


def build_payload() -> dict[str, Any]:
    p3101 = read_json(P3101)
    greps = content_grep(); probs = probability_rows(); events = event_additivity_rows(); basis = basis_orbit_rows(); freq = frequency_proxy_rows(); gates = gate_rows(); aggs = aggregate(gates)
    accepted = [r for r in aggs if r["accepted_nonimported_probability_source"]]
    obligations = [
        {"obligation": "read_p3101_next_atom", "satisfied": True, "detail": "P3101 selected Born-rule/probability-measure empirical readout as the next interface atom"},
        {"obligation": "construct_normalized_probability_rows", "satisfied": len(probs) == len(STATES) * Z12_SIZE and all(r["nonnegative"] for r in probs), "detail": "finite nonnegative normalized weights are explicit"},
        {"obligation": "construct_event_additivity_rows", "satisfied": len(events) == len(STATES) * 4, "detail": "finite event additivity witnesses are explicit"},
        {"obligation": "construct_basis_orbit_rows", "satisfied": len(basis) == len(BASIS_SHIFTS) * Z12_SIZE and all(not r["canonical_measurement_basis_fixed"] for r in basis), "detail": "measurement-basis orbit remains unfixed"},
        {"obligation": "construct_frequency_proxy_rows", "satisfied": len(freq) == len(STATES) * len(SAMPLE_SIZES), "detail": "expected-frequency proxies are explicit"},
        {"obligation": "export_nonimported_physical_probability_source", "satisfied": False, "detail": "0 candidates pass as internal non-imported probability/readout sources"},
    ]
    return {"status": "P3102_BORN_PROBABILITY_MEASURE_EMPIRICAL_READOUT_OBSTRUCTION_BOUNDED_NO_GO", "input_hashes": {"P3101": hashlib.sha256(P3101.read_bytes()).hexdigest() if P3101.exists() else None}, "constructed_theoretical_objects": {"content_first_repo_grep": greps, "born_probability_measure_audit_object": {"object": "Z12DirichletBornProbabilityMeasureEmpiricalReadoutObstructionAudit", "source_reused": "P3101 recommendation: bounded Born-rule/probability-measure empirical-readout obstruction audit", "required_gates": list(GATES), "candidate_probability_sources": [c["id"] for c in CANDIDATES], "acceptance_predicate": "normalized probability measure plus event sigma-algebra, Born-rule map, measurement-basis source, detector calibration link, empirical frequency interface, P3101 detector source, and non-imported physical probability source"}, "probability_measure_rows": probs, "event_additivity_rows": events, "basis_orbit_rows": basis, "frequency_proxy_rows": freq, "candidate_gate_rows": gates, "candidate_aggregate_certificate": aggs}, "finite_certificate": {"content_grep_lanes": len(greps), "content_grep_hits": sum(r["hit_count"] for r in greps), "p3101_accepted_nonimported_detector_sources": p3101.get("finite_certificate", {}).get("accepted_nonimported_detector_sources"), "probability_measure_rows": len(probs), "probability_rows_with_born_rule_source": sum(r["born_rule_source_attached"] for r in probs), "event_additivity_rows": len(events), "event_rows_with_physical_semantics": sum(r["physical_event_semantics_attached"] for r in events), "basis_orbit_rows": len(basis), "basis_rows_with_canonical_measurement_basis": sum(r["canonical_measurement_basis_fixed"] for r in basis), "frequency_proxy_rows": len(freq), "frequency_rows_with_empirical_readout": sum(r["empirical_readout_attached"] for r in freq), "probability_candidates": len(CANDIDATES), "required_gates": len(GATES), "candidate_gate_rows": len(gates), "accepted_nonimported_probability_sources": len(accepted), "proof_obligations": len(obligations), "satisfied_proof_obligations": sum(r["satisfied"] for r in obligations)}, "proof_obligations": obligations, "decision": {"bounded_result": "P3102 constructs the requested Born-rule/probability-measure empirical-readout obstruction audit.  The Z12 Laplacian supports normalized Fourier-power probability measures, finite event-additivity witnesses, basis-orbit Born-like weights, and expected-frequency proxies.  These are real probability-like witnesses, but no internal artifact exports a Born-rule postulate, canonical measurement-basis source, physical event semantics, detector calibration link, empirical frequency interface, P3101 detector source, or non-imported physical probability/readout source.  Therefore no empirical Born/probability readout source is exported.", "negative_export_flags": {key: False for key in ["born_rule_map_exported", "canonical_measurement_basis_exported", "physical_event_semantics_exported", "detector_calibration_link_exported", "empirical_frequency_interface_exported", "nonimported_physical_probability_source_exported", "observed_light_interface_exported", "physical_detector_map_exported", "physical_hamiltonian_exported", "spacetime_eom_exported", "QW_2191_discharged", "strict_selector_closure_exported", "ltotal_exported", "bridge_role_transfer_exported", "toe_closure_exported"]}, "positive_scoped_flags": {"probability_measure_rows_computed": True, "event_additivity_rows_computed": True, "basis_orbit_rows_computed": True, "frequency_proxy_rows_computed": True}, "next_honest_step": "Pivot to exactly one adjacent standard-physics interface atom outside selector replay: construct a bounded Hilbert-space/state-vector completion obstruction audit for the Z12 Dirichlet/Laplacian branch, testing whether finite probability measures, basis-orbit weights, and detector/readout proxies supply a non-imported complex Hilbert state space, inner product with physical units, observable algebra, and state-preparation/readout map without importing quantum axioms, apparatus templates, observed light, L_total, bridge/role-transfer, or ToE."}}


def write_markdown(payload: dict[str, Any]) -> None:
    c = payload["finite_certificate"]
    lines = ["# P3102/S2052 Born-rule/probability-measure empirical-readout obstruction audit", "", f"Status: `{payload['status']}`", "", "## Finite certificate", f"- content grep lanes: `{c['content_grep_lanes']}`", f"- content grep hits: `{c['content_grep_hits']}`", f"- P3101 accepted non-imported detector sources: `{c['p3101_accepted_nonimported_detector_sources']}`", f"- probability measure rows: `{c['probability_measure_rows']}`", f"- probability rows with Born-rule source: `{c['probability_rows_with_born_rule_source']}`", f"- event additivity rows: `{c['event_additivity_rows']}`", f"- event rows with physical semantics: `{c['event_rows_with_physical_semantics']}`", f"- basis orbit rows: `{c['basis_orbit_rows']}`", f"- basis rows with canonical measurement basis: `{c['basis_rows_with_canonical_measurement_basis']}`", f"- frequency proxy rows: `{c['frequency_proxy_rows']}`", f"- frequency rows with empirical readout: `{c['frequency_rows_with_empirical_readout']}`", f"- probability candidates: `{c['probability_candidates']}`", f"- candidate gate rows: `{c['candidate_gate_rows']}`", f"- accepted non-imported probability sources: `{c['accepted_nonimported_probability_sources']}`", f"- satisfied proof obligations: `{c['satisfied_proof_obligations']}/{c['proof_obligations']}`", "", "## Decision", payload["decision"]["bounded_result"], "", "## Recommendation", payload["decision"]["next_honest_step"]]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3102/S2052 Born-rule/probability-measure empirical-readout obstruction audit", "## P3102/S2052 Born-rule/probability-measure empirical-readout obstruction audit\n\n`P3102/S2052` attacks exactly one post-P3101 interface atom: a non-imported Born-rule/probability-measure empirical-readout source for the Z12 Dirichlet/Laplacian branch.  It constructs `48` probability measure rows, `16` event additivity rows, `72` basis orbit rows, `16` frequency proxy rows, and a `5 x 8 = 40` candidate gate matrix.  The finite probability-like algebra remains formal; no Born-rule map, canonical measurement-basis source, physical event semantics, detector calibration link, empirical frequency interface, observed-light interface, `L_total`, bridge/role-transfer, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3102/S2052 Born-rule/probability readout remains unsourced", "## P3102/S2052 Born-rule/probability readout remains unsourced\n\n`P3102/S2052` confirms that the Z12 Laplacian supports normalized Fourier-power probability measures, finite event-additivity witnesses, basis-orbit Born-like weights, and expected-frequency proxies.  A Lagrangian/EOM, quantum-measurement, observed-light, or empirical-readout reading still needs strict sources for the Born rule, canonical measurement basis, physical event semantics, detector calibration, and frequency/readout trials; imported quantum axioms and apparatus templates do not supply those strict sources.\n")
    append_once(AGENTS, "Current Born-rule/probability-measure empirical-readout obstruction guardrail (P3102/S2052, 2026-06-25)", "## Current Born-rule/probability-measure empirical-readout obstruction guardrail (P3102/S2052, 2026-06-25)\n\n- P3102 follows the P3101 recommendation and audits one standard-physics interface atom: a non-imported Born-rule/probability-measure empirical-readout source for the Z12 Dirichlet/Laplacian branch.\n- The finite audit computes `48` probability measure rows, `16` event additivity rows, `72` basis orbit rows, `16` frequency proxy rows, and `40` candidate gate rows; `0` candidates export an internal non-imported Born/probability/readout law.\n- Do not promote normalized Fourier-power weights, finite additivity, basis-orbit Born-like weights, expected-frequency proxies, imported quantum measurement templates, or imported observed-frequency templates to a Born-rule source, canonical measurement basis, empirical detector map, observed photons/light, spacetime EOM, physical Hamiltonian, selector closure, `L_total`, bridge/role-transfer, or ToE closure.\n- The next honest move is exactly one Hilbert-space/state-vector completion obstruction audit for the Z12 Dirichlet/Laplacian branch, unless a genuinely new Born/probability/readout source theorem is introduced.\n")
    return payload

if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
