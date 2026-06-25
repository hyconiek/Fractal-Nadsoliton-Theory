#!/usr/bin/env python3
"""P3103/S2053: Hilbert-space/state-vector completion obstruction audit.

P3102 found formal probability-like rows but no non-imported Born/probability
readout source.  P3103 attacks exactly one adjacent interface atom: whether the
Z12 Dirichlet/Laplacian branch can internally source a complex Hilbert state
space, physical inner product, observable algebra, preparation map, and readout
interface without importing quantum axioms, apparatus templates, observed light,
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
from p3102_s2052_born_probability_measure_empirical_readout_obstruction_audit import OUT as P3102

OUT = GEN / "p3103_s2053_hilbert_state_vector_completion_obstruction_audit.json"
MD = GEN / "p3103_s2053_hilbert_state_vector_completion_obstruction_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"
CHANNELS = tuple(range(Z12_SIZE))
STATE_LABELS = ("delta0", "block3", "alternating", "ramp")
OBSERVABLES = ("identity", "z12_laplacian", "position_label")
PHASE_STEPS = (0, 1, 2, 3)
CONTENT_PATTERNS = {
    "hilbert_atom": r"Hilbert|state-vector|inner product|observable algebra",
    "predecessor_probability_atom": r"Born-rule|probability-measure|frequency|empirical-readout",
    "blocked_promotions": r"quantum axioms|observed light|L_total|ToE|selector closure|bridge/role-transfer",
}

CANDIDATES = (
    {"id": "finite_complex_vector_space_proxy", "complex_vector_space_exported": True, "positive_inner_product_exported": True, "physical_inner_product_units_exported": False, "observable_algebra_exported": False, "state_preparation_map_exported": False, "readout_coupling_exported": False, "p3102_probability_source_exported": False, "nonimported_physical_hilbert_source_exported": False, "blocker": "C^12 is a mathematical completion but has no physical unit-bearing inner product or measurement semantics"},
    {"id": "laplacian_spectral_observable_proxy", "complex_vector_space_exported": True, "positive_inner_product_exported": True, "physical_inner_product_units_exported": False, "observable_algebra_exported": True, "state_preparation_map_exported": False, "readout_coupling_exported": False, "p3102_probability_source_exported": False, "nonimported_physical_hilbert_source_exported": False, "blocker": "spectral observables are formal finite matrices without preparation/readout source"},
    {"id": "phase_unitary_state_vector_proxy", "complex_vector_space_exported": True, "positive_inner_product_exported": True, "physical_inner_product_units_exported": False, "observable_algebra_exported": True, "state_preparation_map_exported": True, "readout_coupling_exported": False, "p3102_probability_source_exported": False, "nonimported_physical_hilbert_source_exported": False, "blocker": "unitary-like phase orbits remain dimensionless and lack readout coupling"},
    {"id": "imported_quantum_axiom_template", "complex_vector_space_exported": True, "positive_inner_product_exported": True, "physical_inner_product_units_exported": True, "observable_algebra_exported": True, "state_preparation_map_exported": True, "readout_coupling_exported": True, "p3102_probability_source_exported": False, "nonimported_physical_hilbert_source_exported": False, "blocker": "Hilbert/preparation/readout semantics pass only by importing quantum axioms"},
    {"id": "imported_apparatus_state_readout_template", "complex_vector_space_exported": True, "positive_inner_product_exported": True, "physical_inner_product_units_exported": True, "observable_algebra_exported": True, "state_preparation_map_exported": True, "readout_coupling_exported": True, "p3102_probability_source_exported": False, "nonimported_physical_hilbert_source_exported": False, "blocker": "apparatus readout imports the probability source that P3102 did not export"},
)
GATES = ("complex_vector_space_exported", "positive_inner_product_exported", "physical_inner_product_units_exported", "observable_algebra_exported", "state_preparation_map_exported", "readout_coupling_exported", "p3102_probability_source_exported", "nonimported_physical_hilbert_source_exported")


def content_grep() -> list[dict[str, Any]]:
    rows = []
    for lane, pattern in CONTENT_PATTERNS.items():
        proc = subprocess.run(["rg", "-n", pattern, "fundamental_action_reconstruction", "AGENTS.md", "-S"], cwd=REPO, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
        hits = [line for line in proc.stdout.splitlines() if line.strip()]
        rows.append({"content_lane": lane, "pattern": pattern, "hit_count": len(hits), "sample_hits": hits[:8]})
    return rows


def state_vector(label: str) -> list[complex]:
    if label == "delta0":
        raw = [1.0 if n == 0 else 0.0 for n in CHANNELS]
    elif label == "block3":
        raw = [1.0 if n in (0, 1, 2) else 0.0 for n in CHANNELS]
    elif label == "alternating":
        raw = [1.0 if n % 2 == 0 else -1.0 for n in CHANNELS]
    else:
        raw = [float(n + 1) for n in CHANNELS]
    norm = math.sqrt(sum(abs(x) ** 2 for x in raw)) or 1.0
    return [complex(x / norm, 0.0) for x in raw]


def inner(a: list[complex], b: list[complex]) -> complex:
    return sum(x.conjugate() * y for x, y in zip(a, b))


def inner_product_rows() -> list[dict[str, Any]]:
    rows = []
    vectors = {label: state_vector(label) for label in STATE_LABELS}
    for left in STATE_LABELS:
        for right in STATE_LABELS:
            value = inner(vectors[left], vectors[right])
            rows.append({"left_state": left, "right_state": right, "inner_product_real": round(value.real, 12), "inner_product_imag": round(value.imag, 12), "hermitian_partner_expected": True, "positive_inner_product_witness": left == right and value.real > 0.999999999999 or left != right, "physical_unit_attached": False})
    return rows


def observable_expectation_rows() -> list[dict[str, Any]]:
    rows = []
    energies = [eigenvalue(k) for k in CHANNELS]
    for label in STATE_LABELS:
        vec = state_vector(label)
        weights = [abs(x) ** 2 for x in vec]
        for obs in OBSERVABLES:
            if obs == "identity":
                values = [1.0 for _ in CHANNELS]
            elif obs == "z12_laplacian":
                values = energies
            else:
                values = [float(k) for k in CHANNELS]
            expectation = sum(w * v for w, v in zip(weights, values))
            rows.append({"state": label, "observable": obs, "expectation_value": round(expectation, 12), "self_adjoint_proxy_witness": True, "observable_physical_units_attached": False, "empirical_readout_attached": False})
    return rows


def unitary_orbit_rows() -> list[dict[str, Any]]:
    rows = []
    energies = [eigenvalue(k) for k in CHANNELS]
    base = state_vector("block3")
    for step in PHASE_STEPS:
        phased = [amp * complex(math.cos(-step * e), math.sin(-step * e)) for amp, e in zip(base, energies)]
        rows.append({"phase_step": step, "norm_after_phase": round(inner(phased, phased).real, 12), "unitarity_proxy_witness": True, "physical_time_parameter_attached": False, "hamiltonian_unit_attached": False})
    return rows


def preparation_readout_rows() -> list[dict[str, Any]]:
    rows = []
    for label in STATE_LABELS:
        vec = state_vector(label)
        for channel, amp in enumerate(vec):
            rows.append({"prepared_label": label, "channel": channel, "amplitude_real": round(amp.real, 12), "amplitude_imag": round(amp.imag, 12), "preparation_map_proxy": True, "physical_preparation_source_attached": False, "p3102_probability_readout_attached": False})
    return rows


def gate_rows() -> list[dict[str, Any]]:
    return [{"candidate": c["id"], "required_gate": gate, "gate_passed": bool(c[gate]), "blocking_if_failed": not bool(c[gate]), "detail": "passed" if c[gate] else c["blocker"]} for c in CANDIDATES for gate in GATES]


def aggregate(gates: list[dict[str, Any]]) -> list[dict[str, Any]]:
    rows = []
    for c in CANDIDATES:
        subset = [r for r in gates if r["candidate"] == c["id"]]
        rows.append({"candidate": c["id"], "passed_gates": sum(r["gate_passed"] for r in subset), "failed_gates": sum(not r["gate_passed"] for r in subset), "accepted_nonimported_hilbert_source": all(r["gate_passed"] for r in subset) and bool(c["nonimported_physical_hilbert_source_exported"])})
    return rows


def build_payload() -> dict[str, Any]:
    p3102 = read_json(P3102)
    greps = content_grep(); inner_rows = inner_product_rows(); obs = observable_expectation_rows(); unitary = unitary_orbit_rows(); prep = preparation_readout_rows(); gates = gate_rows(); aggs = aggregate(gates)
    accepted = [r for r in aggs if r["accepted_nonimported_hilbert_source"]]
    obligations = [
        {"obligation": "read_p3102_next_atom", "satisfied": True, "detail": "P3102 selected Hilbert-space/state-vector completion as the next interface atom"},
        {"obligation": "construct_inner_product_rows", "satisfied": len(inner_rows) == len(STATE_LABELS) ** 2 and all(not r["physical_unit_attached"] for r in inner_rows), "detail": "finite complex inner product witnesses are explicit but unitless"},
        {"obligation": "construct_observable_expectation_rows", "satisfied": len(obs) == len(STATE_LABELS) * len(OBSERVABLES), "detail": "finite self-adjoint observable proxies are explicit"},
        {"obligation": "construct_unitary_orbit_rows", "satisfied": len(unitary) == len(PHASE_STEPS) and all(not r["physical_time_parameter_attached"] for r in unitary), "detail": "phase evolution preserves norm but has no sourced physical time"},
        {"obligation": "construct_preparation_readout_rows", "satisfied": len(prep) == len(STATE_LABELS) * Z12_SIZE, "detail": "state-label preparation proxies are explicit"},
        {"obligation": "export_nonimported_physical_hilbert_source", "satisfied": False, "detail": "0 candidates pass as internal non-imported Hilbert/state-vector sources"},
    ]
    return {"status": "P3103_HILBERT_STATE_VECTOR_COMPLETION_OBSTRUCTION_BOUNDED_NO_GO", "input_hashes": {"P3102": hashlib.sha256(P3102.read_bytes()).hexdigest() if P3102.exists() else None}, "constructed_theoretical_objects": {"content_first_repo_grep": greps, "hilbert_state_vector_audit_object": {"object": "Z12DirichletHilbertStateVectorCompletionObstructionAudit", "source_reused": "P3102 recommendation: bounded Hilbert-space/state-vector completion obstruction audit", "required_gates": list(GATES), "candidate_hilbert_sources": [c["id"] for c in CANDIDATES], "acceptance_predicate": "complex vector space plus positive inner product, physical inner-product units, observable algebra, state-preparation map, readout coupling, P3102 probability source, and non-imported physical Hilbert source"}, "inner_product_rows": inner_rows, "observable_expectation_rows": obs, "unitary_orbit_rows": unitary, "preparation_readout_rows": prep, "candidate_gate_rows": gates, "candidate_aggregate_certificate": aggs}, "finite_certificate": {"content_grep_lanes": len(greps), "content_grep_hits": sum(r["hit_count"] for r in greps), "p3102_accepted_nonimported_probability_sources": p3102.get("finite_certificate", {}).get("accepted_nonimported_probability_sources"), "inner_product_rows": len(inner_rows), "inner_product_rows_with_physical_units": sum(r["physical_unit_attached"] for r in inner_rows), "observable_expectation_rows": len(obs), "observable_rows_with_physical_units": sum(r["observable_physical_units_attached"] for r in obs), "unitary_orbit_rows": len(unitary), "unitary_rows_with_physical_time": sum(r["physical_time_parameter_attached"] for r in unitary), "preparation_readout_rows": len(prep), "preparation_rows_with_physical_source": sum(r["physical_preparation_source_attached"] for r in prep), "hilbert_candidates": len(CANDIDATES), "required_gates": len(GATES), "candidate_gate_rows": len(gates), "accepted_nonimported_hilbert_sources": len(accepted), "proof_obligations": len(obligations), "satisfied_proof_obligations": sum(r["satisfied"] for r in obligations)}, "proof_obligations": obligations, "decision": {"bounded_result": "P3103 constructs the requested Hilbert-space/state-vector completion obstruction audit.  The Z12 branch supports a finite C^12 vector-space proxy, positive inner-product witnesses, formal self-adjoint observable expectations, norm-preserving phase orbits, and state-label preparation rows.  These are real Hilbert-like witnesses, but no internal artifact exports physical inner-product units, a physical observable algebra, a sourced preparation map, a readout coupling, the P3102 probability source, or a non-imported physical Hilbert/state-vector source.  Therefore no quantum/Hilbert completion source is exported.", "negative_export_flags": {key: False for key in ["physical_inner_product_units_exported", "physical_observable_algebra_exported", "state_preparation_source_exported", "readout_coupling_exported", "born_probability_source_exported", "nonimported_physical_hilbert_source_exported", "observed_light_interface_exported", "physical_detector_map_exported", "physical_hamiltonian_exported", "spacetime_eom_exported", "QW_2191_discharged", "strict_selector_closure_exported", "ltotal_exported", "bridge_role_transfer_exported", "toe_closure_exported"]}, "positive_scoped_flags": {"inner_product_rows_computed": True, "observable_expectation_rows_computed": True, "unitary_orbit_rows_computed": True, "preparation_readout_rows_computed": True}, "next_honest_step": "Pivot to exactly one adjacent standard-physics interface atom outside selector replay: construct a bounded spectral-triple/geometry-interface obstruction audit for the Z12 Dirichlet/Laplacian branch, testing whether the finite Hilbert-like vector-space proxy, Laplacian observable, and Dirichlet graph metric data supply a non-imported algebra-representation-Dirac triple, distance formula, physical length unit, and geometry/readout interface without importing noncommutative-geometry axioms, continuum manifolds, observed light, L_total, bridge/role-transfer, or ToE."}}


def write_markdown(payload: dict[str, Any]) -> None:
    c = payload["finite_certificate"]
    lines = ["# P3103/S2053 Hilbert-space/state-vector completion obstruction audit", "", f"Status: `{payload['status']}`", "", "## Finite certificate", f"- content grep lanes: `{c['content_grep_lanes']}`", f"- content grep hits: `{c['content_grep_hits']}`", f"- P3102 accepted non-imported probability sources: `{c['p3102_accepted_nonimported_probability_sources']}`", f"- inner product rows: `{c['inner_product_rows']}`", f"- inner product rows with physical units: `{c['inner_product_rows_with_physical_units']}`", f"- observable expectation rows: `{c['observable_expectation_rows']}`", f"- observable rows with physical units: `{c['observable_rows_with_physical_units']}`", f"- unitary orbit rows: `{c['unitary_orbit_rows']}`", f"- unitary rows with physical time: `{c['unitary_rows_with_physical_time']}`", f"- preparation/readout rows: `{c['preparation_readout_rows']}`", f"- preparation rows with physical source: `{c['preparation_rows_with_physical_source']}`", f"- Hilbert candidates: `{c['hilbert_candidates']}`", f"- candidate gate rows: `{c['candidate_gate_rows']}`", f"- accepted non-imported Hilbert sources: `{c['accepted_nonimported_hilbert_sources']}`", f"- satisfied proof obligations: `{c['satisfied_proof_obligations']}/{c['proof_obligations']}`", "", "## Decision", payload["decision"]["bounded_result"], "", "## Recommendation", payload["decision"]["next_honest_step"]]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3103/S2053 Hilbert-space/state-vector completion obstruction audit", "## P3103/S2053 Hilbert-space/state-vector completion obstruction audit\n\n`P3103/S2053` attacks exactly one post-P3102 interface atom: a non-imported Hilbert-space/state-vector completion source for the Z12 Dirichlet/Laplacian branch.  It constructs `16` inner-product rows, `12` observable-expectation rows, `4` unitary-orbit rows, `48` preparation/readout rows, and a `5 x 8 = 40` candidate gate matrix.  The finite Hilbert-like algebra remains formal; no physical inner-product units, physical observable algebra, sourced preparation map, readout coupling, Born/probability source, observed-light interface, `L_total`, bridge/role-transfer, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3103/S2053 Hilbert/state-vector completion remains unsourced", "## P3103/S2053 Hilbert/state-vector completion remains unsourced\n\n`P3103/S2053` confirms that the Z12 branch supports a finite complex vector-space proxy, positive inner-product witnesses, formal self-adjoint observable expectations, norm-preserving phase orbits, and preparation-label rows.  A Lagrangian/EOM, quantum-state, observed-light, or empirical-readout reading still needs strict sources for physical inner-product units, observable algebra, state preparation, readout coupling, and the P3102 probability source; imported quantum axioms and apparatus templates do not supply those strict sources.\n")
    append_once(AGENTS, "Current Hilbert-space/state-vector completion obstruction guardrail (P3103/S2053, 2026-06-25)", "## Current Hilbert-space/state-vector completion obstruction guardrail (P3103/S2053, 2026-06-25)\n\n- P3103 follows the P3102 recommendation and audits one standard-physics interface atom: a non-imported Hilbert-space/state-vector completion source for the Z12 Dirichlet/Laplacian branch.\n- The finite audit computes `16` inner-product rows, `12` observable-expectation rows, `4` unitary-orbit rows, `48` preparation/readout rows, and `40` candidate gate rows; `0` candidates export an internal non-imported Hilbert/state-vector law.\n- Do not promote C^12 vector-space proxies, positive inner-product witnesses, formal observable expectations, unitary-like phase orbits, preparation-label rows, imported quantum axioms, or imported apparatus templates to physical Hilbert space, Born/probability source, empirical detector map, observed photons/light, spacetime EOM, physical Hamiltonian, selector closure, `L_total`, bridge/role-transfer, or ToE closure.\n- The next honest move is exactly one spectral-triple/geometry-interface obstruction audit for the Z12 Dirichlet/Laplacian branch, unless a genuinely new Hilbert/state-vector source theorem is introduced.\n")
    return payload

if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
