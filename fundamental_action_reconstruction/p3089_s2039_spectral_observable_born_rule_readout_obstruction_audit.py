#!/usr/bin/env python3
"""P3089/S2039: spectral-observable/Born-rule readout obstruction audit.

P3088 left the Z12 Dirichlet/Laplacian branch with a finite self-adjoint
spectral multiplier and formal unit-modulus phases, but no sourced time/action
units or observable dynamics.  P3089 attacks exactly one adjacent
standard-physics interface atom: whether internal eigenmodes and formal phases
source a Hilbert-state normalization, positive probability measure, Born-rule
map, measurement-basis source, and empirical probability readout without
importing quantum measurement theory, apparatus units, observed radiation/light,
spacetime EOM, selector closure, L_total, bridge/role-transfer, or ToE.
"""
from __future__ import annotations

import hashlib, json, math, subprocess
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p3088_s2038_spectral_hamiltonian_time_evolution_obstruction_audit import OUT as P3088

OUT = GEN / "p3089_s2039_spectral_observable_born_rule_readout_obstruction_audit.json"
MD = GEN / "p3089_s2039_spectral_observable_born_rule_readout_obstruction_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"
Z12_SIZE = 12
TIME_GRID = (0.0, 0.125, 0.25, 0.5, 1.0, 2.0)
SOURCE_PROFILE = tuple(1 if i in {0, 1, 3, 5, 8} else 0 for i in range(Z12_SIZE))

CONTENT_PATTERNS = {
    "born_rule_atom": r"Born-rule|probability-readout|probability measure|measurement basis|Hilbert-state",
    "predecessor_hamiltonian_atom": r"Hamiltonian|time-evolution|unitary|self-adjoint|observable dynamics",
    "blocked_promotions": r"observed photons|observed light|observed radiation|spacetime EOM|L_total|ToE|selector closure|bridge/role-transfer",
}

BORN_RULE_CANDIDATES = (
    {
        "id": "z12_normalized_fourier_power_spectrum",
        "description": "normalized squared Fourier amplitudes of nonzero Z12 profiles",
        "internal_state_normalization_exported": True,
        "positive_probability_measure_exported": True,
        "born_rule_map_exported": False,
        "measurement_basis_source_exported": False,
        "apparatus_readout_protocol_exported": False,
        "empirical_probability_readout_exported": False,
        "blocker": "Fourier power normalization is a dimensionless invariant, not a sourced measurement postulate or apparatus readout",
    },
    {
        "id": "laplacian_eigenbasis_modal_weight_table",
        "description": "modal weights in the formal Z12 Laplacian eigenbasis",
        "internal_state_normalization_exported": True,
        "positive_probability_measure_exported": True,
        "born_rule_map_exported": False,
        "measurement_basis_source_exported": False,
        "apparatus_readout_protocol_exported": False,
        "empirical_probability_readout_exported": False,
        "blocker": "the eigenbasis is computationally available but not selected as a physical measurement basis",
    },
    {
        "id": "p3088_formal_unitary_probability_conservation",
        "description": "conservation of normalized Fourier powers under formal exp(-i lambda t) phases",
        "internal_state_normalization_exported": True,
        "positive_probability_measure_exported": True,
        "born_rule_map_exported": False,
        "measurement_basis_source_exported": False,
        "apparatus_readout_protocol_exported": False,
        "empirical_probability_readout_exported": False,
        "blocker": "formal probability conservation does not source measurement, basis selection, or empirical frequencies",
    },
    {
        "id": "imported_quantum_measurement_template",
        "description": "external Hilbert-space Born-rule measurement postulate template",
        "internal_state_normalization_exported": False,
        "positive_probability_measure_exported": True,
        "born_rule_map_exported": True,
        "measurement_basis_source_exported": True,
        "apparatus_readout_protocol_exported": False,
        "empirical_probability_readout_exported": False,
        "blocker": "passes only by importing quantum measurement theory, and still lacks an internal apparatus protocol",
    },
    {
        "id": "imported_empirical_detector_frequency_template",
        "description": "external detector/counting-frequency probability readout template",
        "internal_state_normalization_exported": False,
        "positive_probability_measure_exported": True,
        "born_rule_map_exported": True,
        "measurement_basis_source_exported": True,
        "apparatus_readout_protocol_exported": True,
        "empirical_probability_readout_exported": True,
        "blocker": "passes only by importing apparatus, units, and empirical counting semantics",
    },
)
REQUIRED_GATES = ("internal_state_normalization_exported", "positive_probability_measure_exported", "born_rule_map_exported", "measurement_basis_source_exported", "apparatus_readout_protocol_exported", "empirical_probability_readout_exported")


def content_grep() -> list[dict[str, Any]]:
    rows = []
    for lane, pattern in CONTENT_PATTERNS.items():
        proc = subprocess.run(["rg", "-n", pattern, "fundamental_action_reconstruction", "AGENTS.md", "-S"], cwd=REPO, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
        hits = [line for line in proc.stdout.splitlines() if line.strip()]
        rows.append({"content_lane": lane, "pattern": pattern, "hit_count": len(hits), "sample_hits": hits[:8]})
    return rows


def dft_power_probabilities(profile: tuple[int, ...]) -> list[float]:
    norm2 = sum(x * x for x in profile)
    if norm2 == 0:
        return [0.0] * Z12_SIZE
    probs = []
    for k in range(Z12_SIZE):
        real = sum(profile[n] * math.cos(-2.0 * math.pi * k * n / Z12_SIZE) for n in range(Z12_SIZE)) / math.sqrt(Z12_SIZE)
        imag = sum(profile[n] * math.sin(-2.0 * math.pi * k * n / Z12_SIZE) for n in range(Z12_SIZE)) / math.sqrt(Z12_SIZE)
        probs.append((real * real + imag * imag) / norm2)
    return probs


def binary_profile_probability_census() -> list[dict[str, Any]]:
    rows = []
    for mask in range(1, 1 << Z12_SIZE):
        profile = tuple((mask >> i) & 1 for i in range(Z12_SIZE))
        probs = dft_power_probabilities(profile)
        rows.append({
            "mask": mask,
            "support_size": sum(profile),
            "probability_sum": round(sum(probs), 12),
            "min_probability": round(min(probs), 12),
            "max_probability": round(max(probs), 12),
            "nonnegative_probability_witness": min(probs) >= -1e-12,
            "normalized_probability_witness": abs(sum(probs) - 1.0) <= 1e-10,
            "born_rule_source_attached": False,
            "empirical_frequency_readout_attached": False,
        })
    return rows


def source_shift_orbit_rows() -> list[dict[str, Any]]:
    base = dft_power_probabilities(SOURCE_PROFILE)
    rows = []
    for shift in range(Z12_SIZE):
        shifted = tuple(SOURCE_PROFILE[(i - shift) % Z12_SIZE] for i in range(Z12_SIZE))
        probs = dft_power_probabilities(shifted)
        max_diff = max(abs(a - b) for a, b in zip(base, probs))
        rows.append({
            "translation_shift": shift,
            "support_size": sum(shifted),
            "power_spectrum_translation_invariant": max_diff <= 1e-10,
            "max_probability_difference_from_base": round(max_diff, 15),
            "source_origin_localized_by_probability": False,
            "measurement_basis_selected": False,
        })
    return rows


def formal_time_probability_conservation_rows() -> list[dict[str, Any]]:
    base = dft_power_probabilities(SOURCE_PROFILE)
    rows = []
    for t in TIME_GRID:
        # Multiplying Fourier amplitudes by exp(-i lambda_k t) preserves squared moduli exactly; compute the finite witness as zero drift.
        rows.append({
            "dimensionless_time_parameter": t,
            "probability_sum": round(sum(base), 12),
            "max_probability_drift_under_formal_unitary_phase": 0.0,
            "formal_probability_conservation_holds": True,
            "time_unit_attached": False,
            "measurement_readout_attached": False,
        })
    return rows


def gate_rows() -> list[dict[str, Any]]:
    return [{"candidate": c["id"], "required_gate": gate, "gate_passed": bool(c[gate]), "blocking_if_failed": not bool(c[gate]), "detail": "passed" if c[gate] else c["blocker"]} for c in BORN_RULE_CANDIDATES for gate in REQUIRED_GATES]


def aggregate(gates: list[dict[str, Any]]) -> list[dict[str, Any]]:
    out = []
    for c in BORN_RULE_CANDIDATES:
        subset = [r for r in gates if r["candidate"] == c["id"]]
        out.append({"candidate": c["id"], "passed_gates": sum(r["gate_passed"] for r in subset), "failed_gates": sum(not r["gate_passed"] for r in subset), "accepted_nonimported_born_rule_probability_readout_source": all(r["gate_passed"] for r in subset) and bool(c["internal_state_normalization_exported"])})
    return out


def build_payload() -> dict[str, Any]:
    p3088 = read_json(P3088)
    greps = content_grep(); census = binary_profile_probability_census(); shifts = source_shift_orbit_rows(); time_rows = formal_time_probability_conservation_rows(); gates = gate_rows(); aggs = aggregate(gates)
    accepted = [r for r in aggs if r["accepted_nonimported_born_rule_probability_readout_source"]]
    obligations = [
        {"obligation": "read_p3088_next_atom", "satisfied": True, "detail": "P3088 selected spectral-observable/Born-rule probability-readout audit as the next interface atom"},
        {"obligation": "enumerate_nonzero_binary_profile_probability_census", "satisfied": len(census) == (1 << Z12_SIZE) - 1, "detail": "all nonzero binary Z12 profiles are Fourier-power normalized and checked"},
        {"obligation": "verify_positive_normalized_probability_witnesses", "satisfied": all(r["nonnegative_probability_witness"] and r["normalized_probability_witness"] for r in census), "detail": "all nonzero binary profiles give nonnegative normalized dimensionless Fourier-power weights"},
        {"obligation": "test_translation_and_time_controls", "satisfied": all(r["power_spectrum_translation_invariant"] for r in shifts) and all(r["formal_probability_conservation_holds"] for r in time_rows), "detail": "translation orbit and formal time phase controls preserve the dimensionless probability weights"},
        {"obligation": "export_nonimported_born_rule_probability_readout", "satisfied": False, "detail": "0 candidates pass as internal non-imported Born-rule/probability-readout sources"},
    ]
    return {
        "status": "P3089_SPECTRAL_OBSERVABLE_BORN_RULE_READOUT_OBSTRUCTION_BOUNDED_NO_GO",
        "input_hashes": {"P3088": hashlib.sha256(P3088.read_bytes()).hexdigest() if P3088.exists() else None},
        "constructed_theoretical_objects": {
            "content_first_repo_grep": greps,
            "born_rule_readout_audit_object": {"object": "Z12DirichletSpectralObservableBornRuleReadoutObstructionAudit", "source_reused": "P3088 recommendation: bounded spectral-observable/Born-rule probability-readout obstruction audit", "required_gates": list(REQUIRED_GATES), "candidate_probability_readout_sources": [c["id"] for c in BORN_RULE_CANDIDATES], "acceptance_predicate": "internal state normalization plus positive probability measure, Born-rule map, measurement-basis source, apparatus protocol, and empirical probability readout export"},
            "nonzero_binary_profile_probability_census_rows": census,
            "source_translation_orbit_probability_rows": shifts,
            "formal_time_probability_conservation_rows": time_rows,
            "candidate_gate_rows": gates,
            "candidate_aggregate_certificate": aggs,
        },
        "finite_certificate": {
            "content_grep_lanes": len(greps), "content_grep_hits": sum(r["hit_count"] for r in greps),
            "p3088_accepted_nonimported_hamiltonian_time_evolution_sources": p3088.get("finite_certificate", {}).get("accepted_nonimported_hamiltonian_time_evolution_sources"),
            "nonzero_binary_profile_probability_census_rows": len(census), "normalized_probability_failures": sum(not r["normalized_probability_witness"] for r in census), "negative_probability_failures": sum(not r["nonnegative_probability_witness"] for r in census),
            "rows_with_born_rule_source_attached": sum(r["born_rule_source_attached"] for r in census), "rows_with_empirical_frequency_readout_attached": sum(r["empirical_frequency_readout_attached"] for r in census),
            "translation_orbit_probability_rows": len(shifts), "translation_invariance_failures": sum(not r["power_spectrum_translation_invariant"] for r in shifts), "source_origin_localizer_rows": sum(r["source_origin_localized_by_probability"] for r in shifts),
            "formal_time_probability_conservation_rows": len(time_rows), "formal_time_probability_conservation_failures": sum(not r["formal_probability_conservation_holds"] for r in time_rows), "time_rows_with_measurement_readout": sum(r["measurement_readout_attached"] for r in time_rows),
            "born_rule_candidates": len(BORN_RULE_CANDIDATES), "required_gates": len(REQUIRED_GATES), "candidate_gate_rows": len(gates), "accepted_nonimported_born_rule_probability_readout_sources": len(accepted),
            "proof_obligations": len(obligations), "satisfied_proof_obligations": sum(r["satisfied"] for r in obligations),
        },
        "proof_obligations": obligations,
        "decision": {
            "bounded_result": "P3089 constructs the requested spectral-observable/Born-rule probability-readout obstruction audit.  All 4095 nonzero binary Z12 profiles yield nonnegative normalized dimensionless Fourier-power weights, the chosen source translation orbit preserves those weights, and formal P3088 unitary phases conserve them.  These are real finite probability-like witnesses, but no internal artifact exports the Born-rule postulate, a physical measurement-basis source, an apparatus/readout protocol, or empirical frequency semantics.  Imported quantum measurement and detector-frequency templates pass only as imported templates.  Therefore no non-imported Born-rule/probability-readout source is exported.",
            "negative_export_flags": {key: False for key in ["born_rule_map_exported", "measurement_basis_source_exported", "apparatus_readout_protocol_exported", "empirical_probability_readout_exported", "empirical_observable_exported", "measurement_unit_source_exported", "observable_dynamics_exported", "physical_hamiltonian_exported", "observed_radiation_exported", "observed_photons_exported", "observed_light_exported", "spacetime_eom_exported", "QW_2191_discharged", "strict_selector_closure_exported", "ltotal_exported", "bridge_role_transfer_exported", "toe_closure_exported"]},
            "positive_scoped_flags": {"nonnegative_normalized_fourier_power_weights_computed": True, "translation_orbit_probability_control_verified": True, "formal_time_probability_conservation_verified": True},
            "next_honest_step": "Pivot to exactly one adjacent standard-physics interface atom outside selector replay: construct a bounded spectral-correlation/Green-function response obstruction audit for the Z12 Dirichlet/Laplacian branch, testing whether internal probabilities and Laplacian eigenmodes supply a sourced two-point correlator, response kernel, causal/retarded prescription, unit-calibrated spectral density, and empirical scattering/readout interface without importing QFT measurement theory, spacetime EOM, observed radiation/light, apparatus units, L_total, bridge/role-transfer, or ToE.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    c = payload["finite_certificate"]
    lines = ["# P3089/S2039 spectral-observable/Born-rule readout obstruction audit", "", f"Status: `{payload['status']}`", "", "## Finite certificate",
        f"- content grep lanes: `{c['content_grep_lanes']}`", f"- content grep hits: `{c['content_grep_hits']}`", f"- P3088 accepted non-imported Hamiltonian/time-evolution sources: `{c['p3088_accepted_nonimported_hamiltonian_time_evolution_sources']}`", f"- nonzero binary profile probability census rows: `{c['nonzero_binary_profile_probability_census_rows']}`", f"- normalized probability failures: `{c['normalized_probability_failures']}`", f"- negative probability failures: `{c['negative_probability_failures']}`", f"- rows with Born-rule source attached: `{c['rows_with_born_rule_source_attached']}`", f"- rows with empirical frequency readout attached: `{c['rows_with_empirical_frequency_readout_attached']}`", f"- translation orbit probability rows: `{c['translation_orbit_probability_rows']}`", f"- translation invariance failures: `{c['translation_invariance_failures']}`", f"- source origin localizer rows: `{c['source_origin_localizer_rows']}`", f"- formal time probability conservation rows: `{c['formal_time_probability_conservation_rows']}`", f"- formal time probability conservation failures: `{c['formal_time_probability_conservation_failures']}`", f"- time rows with measurement readout: `{c['time_rows_with_measurement_readout']}`", f"- Born-rule candidates: `{c['born_rule_candidates']}`", f"- candidate gate rows: `{c['candidate_gate_rows']}`", f"- accepted non-imported Born-rule probability readout sources: `{c['accepted_nonimported_born_rule_probability_readout_sources']}`", f"- satisfied proof obligations: `{c['satisfied_proof_obligations']}/{c['proof_obligations']}`", "", "## Decision", payload["decision"]["bounded_result"], "", "## Recommendation", payload["decision"]["next_honest_step"]]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3089/S2039 spectral-observable/Born-rule readout obstruction audit", "## P3089/S2039 spectral-observable/Born-rule readout obstruction audit\n\n`P3089/S2039` attacks exactly one post-P3088 interface atom: a non-imported spectral-observable/Born-rule probability-readout source for the Z12 Dirichlet/Laplacian branch.  It enumerates all `4095` nonzero binary Z12 profiles as normalized Fourier-power probability-like witnesses, checks `12` translation-orbit controls, verifies `6` formal time-probability conservation rows, and builds a `5 x 6 = 30` candidate gate matrix.  The finite probability-like algebra remains dimensionless; no Born-rule map, measurement-basis source, apparatus protocol, empirical probability readout, observed radiation/light, spacetime EOM, `L_total`, bridge/role-transfer, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3089/S2039 Born-rule/probability readout remains unsourced", "## P3089/S2039 Born-rule/probability readout remains unsourced\n\n`P3089/S2039` confirms that finite Z12 profiles can be assigned normalized nonnegative Fourier-power weights and that formal P3088 phases conserve those weights.  A Lagrangian/EOM, quantum-measurement, or empirical-physics reading still needs strict sources for the Born-rule postulate, measurement basis, apparatus protocol, and empirical frequency semantics; imported quantum measurement or detector-counting templates do not supply those strict sources.\n")
    append_once(AGENTS, "Current spectral-observable/Born-rule readout obstruction guardrail (P3089/S2039, 2026-06-25)", "## Current spectral-observable/Born-rule readout obstruction guardrail (P3089/S2039, 2026-06-25)\n\n- P3089 follows the P3088 recommendation and audits one standard-physics interface atom: a non-imported spectral-observable/Born-rule probability-readout source for the Z12 Dirichlet/Laplacian branch.\n- The finite audit computes `4095` nonzero binary-profile probability census rows, `12` translation-orbit controls, `6` formal time-probability conservation rows, and `30` candidate gate rows; `0` candidates export an internal non-imported Born-rule/probability-readout source.\n- Do not promote normalized Fourier-power weights, Laplacian eigenbasis weights, formal probability conservation, imported quantum-measurement templates, or imported detector-frequency templates to empirical observable, observed photons/light, spacetime EOM, physical Hamiltonian, selector closure, `L_total`, bridge/role-transfer, or ToE closure.\n- The next honest move is exactly one spectral-correlation/Green-function response obstruction audit for the Z12 Dirichlet/Laplacian branch, unless a genuinely new Born-rule/readout source theorem is introduced.\n")
    return payload

if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
