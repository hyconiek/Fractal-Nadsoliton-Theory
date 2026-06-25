#!/usr/bin/env python3
"""P3099/S2049: irreversibility/entropy-production thermalization audit.

P3098 found exact finite detailed-balance/KMS-like witnesses but no sourced
physical thermal state.  P3099 attacks exactly one adjacent interface atom:
whether detailed-balance kernels, Gibbs weights, relative-entropy monotonicity
proxies, and modal flux witnesses internally source irreversibility,
entropy-production, a dissipative semigroup, a bath/preparation mechanism, or
empirical thermalization readout without importing continuum nonequilibrium
thermodynamics, apparatus units, observed light, selector closure, L_total,
bridge/role-transfer, or ToE.
"""
from __future__ import annotations

import hashlib, json, math, subprocess
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p3095_s2045_dispersion_propagating_mode_observable_obstruction_audit import Z12_SIZE, eigenvalue
from p3098_s2048_kms_detailed_balance_thermal_state_obstruction_audit import OUT as P3098

OUT = GEN / "p3099_s2049_irreversibility_entropy_production_thermalization_obstruction_audit.json"
MD = GEN / "p3099_s2049_irreversibility_entropy_production_thermalization_obstruction_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"
BETA_GRID = (0.25, 0.5, 1.0, 2.0)
TIME_STEPS = tuple(range(6))
CHANNELS = tuple(range(Z12_SIZE))
INITIAL_DISTRIBUTIONS = ("delta_0", "delta_1", "even_uniform", "high_energy_pair")

CONTENT_PATTERNS = {
    "irreversibility_atom": r"irreversibility|entropy-production|thermalization|time arrow|dissipative semigroup",
    "predecessor_kms_atom": r"KMS|detailed-balance|thermal-state|fluctuation-dissipation",
    "blocked_promotions": r"continuum nonequilibrium thermodynamics|apparatus units|observed light|L_total|ToE|selector closure|bridge/role-transfer",
}

THERMALIZATION_CANDIDATES = (
    {"id": "finite_markov_relaxation_semigroup", "description": "finite stochastic powers of the P3098 detailed-balance kernel", "finite_detailed_balance_kernel_exported": True, "relative_entropy_monotonicity_exported": True, "entropy_production_proxy_exported": False, "physical_time_arrow_exported": False, "dissipative_semigroup_exported": False, "bath_preparation_mechanism_exported": False, "empirical_thermalization_readout_exported": False, "nonimported_physical_thermalization_source_exported": False, "blocker": "finite stochastic relaxation is formal and lacks physical time, bath, and readout semantics"},
    {"id": "schnakenberg_entropy_production_proxy", "description": "edge-current log-ratio entropy-production proxy for the finite kernel", "finite_detailed_balance_kernel_exported": True, "relative_entropy_monotonicity_exported": True, "entropy_production_proxy_exported": True, "physical_time_arrow_exported": False, "dissipative_semigroup_exported": False, "bath_preparation_mechanism_exported": False, "empirical_thermalization_readout_exported": False, "nonimported_physical_thermalization_source_exported": False, "blocker": "entropy-production formula is a dimensionless proxy without a physical clock or bath"},
    {"id": "modal_flux_thermalization_proxy", "description": "formal relaxation of modal energy/flux labels under the finite kernel", "finite_detailed_balance_kernel_exported": True, "relative_entropy_monotonicity_exported": True, "entropy_production_proxy_exported": True, "physical_time_arrow_exported": False, "dissipative_semigroup_exported": False, "bath_preparation_mechanism_exported": False, "empirical_thermalization_readout_exported": False, "nonimported_physical_thermalization_source_exported": False, "blocker": "modal relaxation has no physical dissipative dynamics or apparatus semantics"},
    {"id": "imported_nonequilibrium_thermodynamics_template", "description": "external irreversible thermodynamics/time-arrow template", "finite_detailed_balance_kernel_exported": True, "relative_entropy_monotonicity_exported": True, "entropy_production_proxy_exported": True, "physical_time_arrow_exported": True, "dissipative_semigroup_exported": True, "bath_preparation_mechanism_exported": False, "empirical_thermalization_readout_exported": False, "nonimported_physical_thermalization_source_exported": False, "blocker": "time-arrow and dissipation semantics are imported and lack an internal bath/readout source"},
    {"id": "imported_bath_apparatus_readout_template", "description": "external thermal bath/preparation and apparatus-calibrated readout template", "finite_detailed_balance_kernel_exported": True, "relative_entropy_monotonicity_exported": True, "entropy_production_proxy_exported": True, "physical_time_arrow_exported": True, "dissipative_semigroup_exported": True, "bath_preparation_mechanism_exported": True, "empirical_thermalization_readout_exported": True, "nonimported_physical_thermalization_source_exported": False, "blocker": "passes only by importing bath/preparation and empirical apparatus semantics"},
)
REQUIRED_GATES = ("finite_detailed_balance_kernel_exported", "relative_entropy_monotonicity_exported", "entropy_production_proxy_exported", "physical_time_arrow_exported", "dissipative_semigroup_exported", "bath_preparation_mechanism_exported", "empirical_thermalization_readout_exported", "nonimported_physical_thermalization_source_exported")


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
    z = sum(math.exp(-beta * value) for value in e)
    return [math.exp(-beta * value) / z for value in e]


def transition_matrix(e: list[float], beta: float) -> list[list[float]]:
    proposal = 1.0 / (Z12_SIZE - 1)
    matrix = [[0.0 for _ in CHANNELS] for _ in CHANNELS]
    for i in CHANNELS:
        off_sum = 0.0
        for j in CHANNELS:
            if i == j:
                continue
            rate = proposal * min(1.0, math.exp(-beta * (e[j] - e[i])))
            matrix[i][j] = rate
            off_sum += rate
        matrix[i][i] = 1.0 - off_sum
    return matrix


def initial_distribution(name: str, e: list[float]) -> list[float]:
    dist = [0.0] * Z12_SIZE
    if name == "delta_0":
        dist[0] = 1.0
    elif name == "delta_1":
        dist[1] = 1.0
    elif name == "even_uniform":
        for k in range(0, Z12_SIZE, 2):
            dist[k] = 1.0 / 6.0
    elif name == "high_energy_pair":
        high = sorted(CHANNELS, key=lambda k: e[k], reverse=True)[:2]
        for k in high:
            dist[k] = 0.5
    return dist


def step(dist: list[float], matrix: list[list[float]]) -> list[float]:
    return [sum(dist[i] * matrix[i][j] for i in CHANNELS) for j in CHANNELS]


def relative_entropy(dist: list[float], pi: list[float]) -> float:
    return sum(p * math.log(p / q) for p, q in zip(dist, pi) if p > 0.0 and q > 0.0)


def entropy_monotonicity_rows(e: list[float]) -> list[dict[str, Any]]:
    rows = []
    for beta in BETA_GRID:
        pi = gibbs(e, beta)
        matrix = transition_matrix(e, beta)
        for name in INITIAL_DISTRIBUTIONS:
            dist = initial_distribution(name, e)
            previous = None
            for t in TIME_STEPS:
                d = relative_entropy(dist, pi)
                rows.append({"formal_beta": beta, "initial_distribution": name, "step": t, "relative_entropy_to_gibbs": round(d, 12), "monotone_nonincreasing_from_previous": True if previous is None else d <= previous + 1e-12, "relative_entropy_monotonicity_witness": True, "physical_time_step_attached": False, "thermalization_readout_attached": False})
                previous = d
                dist = step(dist, matrix)
    return rows


def entropy_production_rows(e: list[float]) -> list[dict[str, Any]]:
    rows = []
    for beta in BETA_GRID:
        pi = gibbs(e, beta)
        matrix = transition_matrix(e, beta)
        seed = initial_distribution("delta_1", e)
        dist = step(seed, matrix)
        for i in CHANNELS:
            for j in CHANNELS:
                if i == j:
                    continue
                flow_ij = dist[i] * matrix[i][j]
                flow_ji = dist[j] * matrix[j][i]
                if flow_ij > 0.0 and flow_ji > 0.0:
                    sigma = 0.5 * (flow_ij - flow_ji) * math.log(flow_ij / flow_ji)
                else:
                    sigma = 0.0
                rows.append({"formal_beta": beta, "from_channel": i, "to_channel": j, "flow_ij": round(flow_ij, 12), "flow_ji": round(flow_ji, 12), "entropy_production_proxy": round(sigma, 12), "finite_entropy_production_proxy_witness": True, "physical_bath_attached": False, "dissipative_time_arrow_attached": False})
    return rows


def semigroup_rows(e: list[float]) -> list[dict[str, Any]]:
    rows = []
    for beta in BETA_GRID:
        matrix = transition_matrix(e, beta)
        row_sums = [sum(row) for row in matrix]
        min_entry = min(min(row) for row in matrix)
        pi = gibbs(e, beta)
        stationarity_residual = max(abs(sum(pi[i] * matrix[i][j] for i in CHANNELS) - pi[j]) for j in CHANNELS)
        rows.append({"formal_beta": beta, "min_transition_entry": round(min_entry, 12), "max_row_sum_error": round(max(abs(s - 1.0) for s in row_sums), 12), "gibbs_stationarity_residual": round(stationarity_residual, 12), "finite_stochastic_semigroup_proxy_witness": True, "physical_dissipative_semigroup_attached": False, "continuous_time_generator_attached": False})
    return rows


def modal_flux_rows(e: list[float]) -> list[dict[str, Any]]:
    rows = []
    velocities = [abs(math.sin(2.0 * math.pi * k / Z12_SIZE)) for k in CHANNELS]
    for beta in BETA_GRID:
        matrix = transition_matrix(e, beta)
        dist = initial_distribution("high_energy_pair", e)
        for t in TIME_STEPS:
            energy_mean = sum(dist[k] * e[k] for k in CHANNELS)
            flux_mean = sum(dist[k] * velocities[k] for k in CHANNELS)
            rows.append({"formal_beta": beta, "step": t, "mean_energy_label": round(energy_mean, 12), "mean_modal_flux_proxy": round(flux_mean, 12), "modal_relaxation_witness": True, "physical_energy_unit_attached": False, "empirical_flux_readout_attached": False})
            dist = step(dist, matrix)
    return rows


def gate_rows() -> list[dict[str, Any]]:
    return [{"candidate": c["id"], "required_gate": gate, "gate_passed": bool(c[gate]), "blocking_if_failed": not bool(c[gate]), "detail": "passed" if c[gate] else c["blocker"]} for c in THERMALIZATION_CANDIDATES for gate in REQUIRED_GATES]


def aggregate(gates: list[dict[str, Any]]) -> list[dict[str, Any]]:
    out = []
    for c in THERMALIZATION_CANDIDATES:
        subset = [r for r in gates if r["candidate"] == c["id"]]
        out.append({"candidate": c["id"], "passed_gates": sum(r["gate_passed"] for r in subset), "failed_gates": sum(not r["gate_passed"] for r in subset), "accepted_nonimported_thermalization_source": all(r["gate_passed"] for r in subset) and bool(c["nonimported_physical_thermalization_source_exported"])})
    return out


def build_payload() -> dict[str, Any]:
    p3098 = read_json(P3098)
    greps = content_grep(); e = energies(); entropy = entropy_monotonicity_rows(e); production = entropy_production_rows(e); semigroup = semigroup_rows(e); flux = modal_flux_rows(e); gates = gate_rows(); aggs = aggregate(gates)
    accepted = [r for r in aggs if r["accepted_nonimported_thermalization_source"]]
    obligations = [
        {"obligation": "read_p3098_next_atom", "satisfied": True, "detail": "P3098 selected irreversibility/entropy-production thermalization audit as the next interface atom"},
        {"obligation": "construct_relative_entropy_monotonicity_rows", "satisfied": len(entropy) == len(BETA_GRID) * len(INITIAL_DISTRIBUTIONS) * len(TIME_STEPS) and all(r["monotone_nonincreasing_from_previous"] for r in entropy), "detail": "finite relative-entropy rows are monotone under the formal kernel"},
        {"obligation": "construct_entropy_production_proxy_rows", "satisfied": len(production) == len(BETA_GRID) * Z12_SIZE * (Z12_SIZE - 1), "detail": "finite edge-current entropy-production proxies are explicit"},
        {"obligation": "construct_stochastic_semigroup_rows", "satisfied": len(semigroup) == len(BETA_GRID) and all(r["max_row_sum_error"] == 0.0 for r in semigroup), "detail": "formal kernels are stochastic and Gibbs-stationary"},
        {"obligation": "construct_modal_flux_relaxation_rows", "satisfied": len(flux) == len(BETA_GRID) * len(TIME_STEPS), "detail": "modal energy/flux relaxation proxies are explicit"},
        {"obligation": "export_nonimported_physical_thermalization_source", "satisfied": False, "detail": "0 candidates pass as internal non-imported irreversibility/thermalization sources"},
    ]
    return {"status": "P3099_IRREVERSIBILITY_ENTROPY_PRODUCTION_THERMALIZATION_OBSTRUCTION_BOUNDED_NO_GO", "input_hashes": {"P3098": hashlib.sha256(P3098.read_bytes()).hexdigest() if P3098.exists() else None}, "constructed_theoretical_objects": {"content_first_repo_grep": greps, "irreversibility_entropy_production_audit_object": {"object": "Z12DirichletIrreversibilityEntropyProductionThermalizationObstructionAudit", "source_reused": "P3098 recommendation: bounded irreversibility/entropy-production thermalization obstruction audit", "required_gates": list(REQUIRED_GATES), "candidate_thermalization_sources": [c["id"] for c in THERMALIZATION_CANDIDATES], "acceptance_predicate": "finite detailed-balance kernel plus relative-entropy monotonicity, entropy-production proxy, physical time arrow, dissipative semigroup, bath/preparation mechanism, empirical thermalization readout, and non-imported physical thermalization source"}, "relative_entropy_monotonicity_rows": entropy, "entropy_production_proxy_rows": production, "stochastic_semigroup_proxy_rows": semigroup, "modal_flux_relaxation_rows": flux, "candidate_gate_rows": gates, "candidate_aggregate_certificate": aggs}, "finite_certificate": {"content_grep_lanes": len(greps), "content_grep_hits": sum(r["hit_count"] for r in greps), "p3098_accepted_nonimported_kms_thermal_state_sources": p3098.get("finite_certificate", {}).get("accepted_nonimported_kms_thermal_state_sources"), "relative_entropy_monotonicity_rows": len(entropy), "relative_entropy_rows_monotone": sum(r["monotone_nonincreasing_from_previous"] for r in entropy), "relative_entropy_rows_with_physical_time_step": sum(r["physical_time_step_attached"] for r in entropy), "entropy_production_proxy_rows": len(production), "entropy_production_rows_with_physical_bath": sum(r["physical_bath_attached"] for r in production), "entropy_production_rows_with_dissipative_time_arrow": sum(r["dissipative_time_arrow_attached"] for r in production), "stochastic_semigroup_proxy_rows": len(semigroup), "semigroup_rows_with_physical_dissipative_semigroup": sum(r["physical_dissipative_semigroup_attached"] for r in semigroup), "modal_flux_relaxation_rows": len(flux), "modal_flux_rows_with_empirical_flux_readout": sum(r["empirical_flux_readout_attached"] for r in flux), "thermalization_candidates": len(THERMALIZATION_CANDIDATES), "required_gates": len(REQUIRED_GATES), "candidate_gate_rows": len(gates), "accepted_nonimported_thermalization_sources": len(accepted), "proof_obligations": len(obligations), "satisfied_proof_obligations": sum(r["satisfied"] for r in obligations)}, "proof_obligations": obligations, "decision": {"bounded_result": "P3099 constructs the requested irreversibility/entropy-production thermalization obstruction audit.  The Z12 Laplacian supports formal detailed-balance kernels, relative-entropy-to-Gibbs monotonicity rows, finite edge-current entropy-production proxies, stochastic semigroup/stationarity proxies, and modal energy/flux relaxation witnesses.  These are real thermalization-like witnesses, but no internal artifact exports a physical time arrow, physical bath/preparation mechanism, dissipative semigroup source, empirical thermalization readout, or a non-imported physical irreversibility source.  Imported nonequilibrium thermodynamics and apparatus/bath templates pass only as imported templates.  Therefore no non-imported irreversibility/entropy-production thermalization source is exported.", "negative_export_flags": {key: False for key in ["physical_time_arrow_exported", "physical_bath_preparation_mechanism_exported", "dissipative_semigroup_source_exported", "empirical_thermalization_readout_exported", "nonimported_physical_irreversibility_source_exported", "physical_temperature_clock_exported", "operator_algebra_kms_state_exported", "physical_fluctuation_dissipation_relation_exported", "observed_radiation_exported", "observed_light_exported", "spacetime_eom_exported", "physical_hamiltonian_exported", "kms_thermal_state_closure_exported", "thermodynamic_radiation_closure_exported", "QW_2191_discharged", "strict_selector_closure_exported", "ltotal_exported", "bridge_role_transfer_exported", "toe_closure_exported"]}, "positive_scoped_flags": {"relative_entropy_monotonicity_rows_computed": True, "entropy_production_proxy_rows_computed": True, "stochastic_semigroup_proxy_rows_computed": True, "modal_flux_relaxation_rows_computed": True}, "next_honest_step": "Pivot to exactly one adjacent standard-physics interface atom outside selector replay: construct a bounded open-system bath/preparation source obstruction audit for the Z12 Dirichlet/Laplacian branch, testing whether finite Markov kernels, entropy-production proxies, modal flux relaxation, and Gibbs stationarity supply a non-imported bath coupling, preparation map, physical clock, and empirical thermalization interface without importing apparatus units, continuum open-system dynamics, observed light, L_total, bridge/role-transfer, or ToE."}}


def write_markdown(payload: dict[str, Any]) -> None:
    c = payload["finite_certificate"]
    lines = ["# P3099/S2049 irreversibility/entropy-production thermalization obstruction audit", "", f"Status: `{payload['status']}`", "", "## Finite certificate", f"- content grep lanes: `{c['content_grep_lanes']}`", f"- content grep hits: `{c['content_grep_hits']}`", f"- P3098 accepted non-imported KMS thermal-state sources: `{c['p3098_accepted_nonimported_kms_thermal_state_sources']}`", f"- relative-entropy monotonicity rows: `{c['relative_entropy_monotonicity_rows']}`", f"- relative-entropy rows monotone: `{c['relative_entropy_rows_monotone']}`", f"- relative-entropy rows with physical time step: `{c['relative_entropy_rows_with_physical_time_step']}`", f"- entropy-production proxy rows: `{c['entropy_production_proxy_rows']}`", f"- entropy-production rows with physical bath: `{c['entropy_production_rows_with_physical_bath']}`", f"- stochastic semigroup proxy rows: `{c['stochastic_semigroup_proxy_rows']}`", f"- semigroup rows with physical dissipative semigroup: `{c['semigroup_rows_with_physical_dissipative_semigroup']}`", f"- modal flux relaxation rows: `{c['modal_flux_relaxation_rows']}`", f"- modal flux rows with empirical readout: `{c['modal_flux_rows_with_empirical_flux_readout']}`", f"- thermalization candidates: `{c['thermalization_candidates']}`", f"- candidate gate rows: `{c['candidate_gate_rows']}`", f"- accepted non-imported thermalization sources: `{c['accepted_nonimported_thermalization_sources']}`", f"- satisfied proof obligations: `{c['satisfied_proof_obligations']}/{c['proof_obligations']}`", "", "## Decision", payload["decision"]["bounded_result"], "", "## Recommendation", payload["decision"]["next_honest_step"]]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3099/S2049 irreversibility/entropy-production thermalization obstruction audit", "## P3099/S2049 irreversibility/entropy-production thermalization obstruction audit\n\n`P3099/S2049` attacks exactly one post-P3098 interface atom: a non-imported irreversibility/entropy-production thermalization source for the Z12 Dirichlet/Laplacian branch.  It constructs `96` relative-entropy monotonicity rows, `528` entropy-production proxy rows, `4` stochastic semigroup proxy rows, `24` modal flux relaxation rows, and a `5 x 8 = 40` candidate gate matrix.  The finite thermalization-like algebra remains formal; no physical time arrow, bath/preparation mechanism, dissipative semigroup source, empirical thermalization readout, `L_total`, bridge/role-transfer, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3099/S2049 irreversibility thermalization source remains unsourced", "## P3099/S2049 irreversibility thermalization source remains unsourced\n\n`P3099/S2049` confirms that the Z12 Laplacian supports formal entropy monotonicity, edge-current entropy-production proxies, stochastic stationarity, and modal relaxation witnesses.  A Lagrangian/EOM, nonequilibrium, dissipative, or empirical-thermalization reading still needs strict sources for physical time arrows, bath/preparation mechanisms, dissipative semigroup semantics, physical clocks, and readout calibration; imported nonequilibrium thermodynamics and apparatus templates do not supply those strict sources.\n")
    append_once(AGENTS, "Current irreversibility/entropy-production thermalization obstruction guardrail (P3099/S2049, 2026-06-25)", "## Current irreversibility/entropy-production thermalization obstruction guardrail (P3099/S2049, 2026-06-25)\n\n- P3099 follows the P3098 recommendation and audits one standard-physics interface atom: a non-imported irreversibility/entropy-production thermalization source for the Z12 Dirichlet/Laplacian branch.\n- The finite audit computes `96` relative-entropy monotonicity rows, `528` entropy-production proxy rows, `4` stochastic semigroup proxy rows, `24` modal flux relaxation rows, and `40` candidate gate rows; `0` candidates export an internal non-imported irreversibility/thermalization law.\n- Do not promote formal relative-entropy monotonicity, finite entropy-production proxies, stochastic semigroup/stationarity rows, modal relaxation witnesses, imported nonequilibrium thermodynamics templates, or imported apparatus/bath templates to physical time arrow, physical bath, empirical observable, observed photons/light, spacetime EOM, physical Hamiltonian, selector closure, `L_total`, bridge/role-transfer, or ToE closure.\n- The next honest move is exactly one open-system bath/preparation source obstruction audit for the Z12 Dirichlet/Laplacian branch, unless a genuinely new irreversibility/thermalization theorem is introduced.\n")
    return payload

if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
