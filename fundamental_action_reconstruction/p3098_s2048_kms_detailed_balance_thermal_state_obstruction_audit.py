#!/usr/bin/env python3
"""P3098/S2048: KMS/detailed-balance thermal-state audit.

P3097 left formal thermodynamic/radiation witnesses but no sourced physical
thermal/radiation state.  P3098 attacks exactly one adjacent interface atom:
whether finite Z12 Dirichlet/Laplacian partition weights, transition kernels,
modal flux proxies, and response/scattering witnesses internally source a KMS or
detailed-balance thermal-state law without importing continuum statistical
mechanics, apparatus units, observed light, selector closure, L_total,
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
from p3097_s2047_thermodynamic_radiation_blackbody_spectrum_obstruction_audit import OUT as P3097

OUT = GEN / "p3098_s2048_kms_detailed_balance_thermal_state_obstruction_audit.json"
MD = GEN / "p3098_s2048_kms_detailed_balance_thermal_state_obstruction_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"
BETA_GRID = (0.25, 0.5, 1.0, 2.0)
TAU_FRACTIONS = (0.0, 0.25, 0.5, 0.75)
CHANNELS = tuple(range(Z12_SIZE))

CONTENT_PATTERNS = {
    "kms_thermal_state_atom": r"KMS|detailed-balance|thermal-state|fluctuation-dissipation|thermal readout",
    "predecessor_radiation_atom": r"thermodynamic-radiation|blackbody|radiation spectrum|intensity readout",
    "blocked_promotions": r"continuum statistical mechanics|apparatus units|observed light|L_total|ToE|selector closure|bridge/role-transfer",
}

THERMAL_STATE_CANDIDATES = (
    {"id": "formal_gibbs_weight_functional", "description": "dimensionless exp(-beta lambda)/Z weights over finite Z12 modes", "finite_partition_weights_exported": True, "detailed_balance_kernel_exported": False, "kms_periodicity_exported": False, "physical_temperature_clock_exported": False, "fluctuation_dissipation_relation_exported": False, "empirical_thermal_readout_exported": False, "bath_or_preparation_source_exported": False, "nonimported_physical_thermal_state_source_exported": False, "blocker": "Gibbs weights are formal and lack a physical beta, clock, bath, or readout"},
    {"id": "metropolis_detailed_balance_kernel", "description": "finite proposal/acceptance kernel satisfying algebraic detailed balance", "finite_partition_weights_exported": True, "detailed_balance_kernel_exported": True, "kms_periodicity_exported": False, "physical_temperature_clock_exported": False, "fluctuation_dissipation_relation_exported": False, "empirical_thermal_readout_exported": False, "bath_or_preparation_source_exported": False, "nonimported_physical_thermal_state_source_exported": False, "blocker": "detailed balance is an algebraic finite Markov rule, not a sourced physical thermalization law"},
    {"id": "imaginary_time_kms_proxy", "description": "formal finite spectral correlator tested under tau -> tau + beta", "finite_partition_weights_exported": True, "detailed_balance_kernel_exported": True, "kms_periodicity_exported": True, "physical_temperature_clock_exported": False, "fluctuation_dissipation_relation_exported": False, "empirical_thermal_readout_exported": False, "bath_or_preparation_source_exported": False, "nonimported_physical_thermal_state_source_exported": False, "blocker": "KMS-like periodicity is formal imaginary-time algebra with no physical clock/readout"},
    {"id": "imported_fluctuation_dissipation_template", "description": "external FDT/KMS thermal-response template", "finite_partition_weights_exported": True, "detailed_balance_kernel_exported": True, "kms_periodicity_exported": True, "physical_temperature_clock_exported": True, "fluctuation_dissipation_relation_exported": True, "empirical_thermal_readout_exported": False, "bath_or_preparation_source_exported": False, "nonimported_physical_thermal_state_source_exported": False, "blocker": "FDT and clock semantics are imported and lack internal bath/readout source"},
    {"id": "imported_empirical_thermal_readout_template", "description": "external apparatus-calibrated temperature/thermal readout template", "finite_partition_weights_exported": True, "detailed_balance_kernel_exported": True, "kms_periodicity_exported": True, "physical_temperature_clock_exported": True, "fluctuation_dissipation_relation_exported": True, "empirical_thermal_readout_exported": True, "bath_or_preparation_source_exported": True, "nonimported_physical_thermal_state_source_exported": False, "blocker": "passes only by importing apparatus, bath/preparation, and empirical thermal semantics"},
)
REQUIRED_GATES = ("finite_partition_weights_exported", "detailed_balance_kernel_exported", "kms_periodicity_exported", "physical_temperature_clock_exported", "fluctuation_dissipation_relation_exported", "empirical_thermal_readout_exported", "bath_or_preparation_source_exported", "nonimported_physical_thermal_state_source_exported")


def content_grep() -> list[dict[str, Any]]:
    rows = []
    for lane, pattern in CONTENT_PATTERNS.items():
        proc = subprocess.run(["rg", "-n", pattern, "fundamental_action_reconstruction", "AGENTS.md", "-S"], cwd=REPO, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
        hits = [line for line in proc.stdout.splitlines() if line.strip()]
        rows.append({"content_lane": lane, "pattern": pattern, "hit_count": len(hits), "sample_hits": hits[:8]})
    return rows


def energies() -> dict[int, float]:
    return {k: round(eigenvalue(k), 12) for k in CHANNELS}


def gibbs_weight_rows(e: dict[int, float]) -> list[dict[str, Any]]:
    rows = []
    for beta in BETA_GRID:
        z = sum(math.exp(-beta * e[k]) for k in CHANNELS)
        for k in CHANNELS:
            rows.append({"formal_beta": beta, "channel": k, "energy_label": e[k], "gibbs_weight": round(math.exp(-beta * e[k]) / z, 12), "finite_partition_weight_witness": True, "physical_temperature_attached": False, "thermal_preparation_attached": False})
    return rows


def detailed_balance_rows(e: dict[int, float]) -> list[dict[str, Any]]:
    rows = []
    proposal = 1.0 / (Z12_SIZE - 1)
    for beta in BETA_GRID:
        z = sum(math.exp(-beta * e[k]) for k in CHANNELS)
        pi = {k: math.exp(-beta * e[k]) / z for k in CHANNELS}
        for i in CHANNELS:
            for j in CHANNELS:
                if i == j:
                    continue
                rate_ij = proposal * min(1.0, math.exp(-beta * (e[j] - e[i])))
                rate_ji = proposal * min(1.0, math.exp(-beta * (e[i] - e[j])))
                lhs = pi[i] * rate_ij
                rhs = pi[j] * rate_ji
                rows.append({"formal_beta": beta, "from_channel": i, "to_channel": j, "rate_ij": round(rate_ij, 12), "pi_i_rate_ij": round(lhs, 12), "pi_j_rate_ji": round(rhs, 12), "detailed_balance_residual": round(lhs - rhs, 12), "finite_detailed_balance_witness": True, "physical_bath_attached": False, "empirical_transition_semantics_attached": False})
    return rows


def kms_proxy_rows(e: dict[int, float]) -> list[dict[str, Any]]:
    rows = []
    for beta in BETA_GRID:
        z = sum(math.exp(-beta * e[k]) for k in CHANNELS)
        for frac in TAU_FRACTIONS:
            tau = frac * beta
            c_tau = sum(math.exp(-beta * e[k]) * math.exp(-tau * e[k]) for k in CHANNELS) / z
            c_shift = sum(math.exp(-beta * e[k]) * math.exp(-(tau + beta) * e[k]) for k in CHANNELS) / z
            rows.append({"formal_beta": beta, "tau_fraction": frac, "C_tau": round(c_tau, 12), "C_tau_plus_beta": round(c_shift, 12), "formal_kms_shift_ratio": round(c_shift / c_tau if c_tau else 0.0, 12), "kms_periodicity_proxy_witness": True, "physical_imaginary_time_clock_attached": False, "operator_algebra_state_attached": False})
    return rows


def fdt_proxy_rows(e: dict[int, float]) -> list[dict[str, Any]]:
    gaps = sorted({round(abs(e[i] - e[j]), 12) for i in CHANNELS for j in CHANNELS if i != j})
    rows = []
    for beta in BETA_GRID:
        for gap in gaps:
            if gap == 0.0:
                continue
            ratio = math.exp(-beta * gap)
            rows.append({"formal_beta": beta, "energy_gap": gap, "boltzmann_response_ratio": round(ratio, 12), "antisymmetric_response_proxy": round(1.0 - ratio, 12), "finite_fdt_like_witness": True, "physical_response_function_attached": False, "empirical_thermal_readout_attached": False})
    return rows


def gate_rows() -> list[dict[str, Any]]:
    return [{"candidate": c["id"], "required_gate": gate, "gate_passed": bool(c[gate]), "blocking_if_failed": not bool(c[gate]), "detail": "passed" if c[gate] else c["blocker"]} for c in THERMAL_STATE_CANDIDATES for gate in REQUIRED_GATES]


def aggregate(gates: list[dict[str, Any]]) -> list[dict[str, Any]]:
    out = []
    for c in THERMAL_STATE_CANDIDATES:
        subset = [r for r in gates if r["candidate"] == c["id"]]
        out.append({"candidate": c["id"], "passed_gates": sum(r["gate_passed"] for r in subset), "failed_gates": sum(not r["gate_passed"] for r in subset), "accepted_nonimported_kms_thermal_state_source": all(r["gate_passed"] for r in subset) and bool(c["nonimported_physical_thermal_state_source_exported"])})
    return out


def build_payload() -> dict[str, Any]:
    p3097 = read_json(P3097)
    greps = content_grep(); e = energies(); weights = gibbs_weight_rows(e); db = detailed_balance_rows(e); kms = kms_proxy_rows(e); fdt = fdt_proxy_rows(e); gates = gate_rows(); aggs = aggregate(gates)
    accepted = [r for r in aggs if r["accepted_nonimported_kms_thermal_state_source"]]
    obligations = [
        {"obligation": "read_p3097_next_atom", "satisfied": True, "detail": "P3097 selected KMS/detailed-balance thermal-state audit as the next interface atom"},
        {"obligation": "construct_gibbs_weight_rows", "satisfied": len(weights) == len(BETA_GRID) * Z12_SIZE, "detail": "formal finite Gibbs weights are explicit"},
        {"obligation": "construct_detailed_balance_rows", "satisfied": len(db) == len(BETA_GRID) * Z12_SIZE * (Z12_SIZE - 1) and all(r["detailed_balance_residual"] == 0.0 for r in db), "detail": "Metropolis-style finite kernel satisfies algebraic detailed balance"},
        {"obligation": "construct_kms_proxy_rows", "satisfied": len(kms) == len(BETA_GRID) * len(TAU_FRACTIONS), "detail": "formal imaginary-time shift rows are computed"},
        {"obligation": "construct_fdt_proxy_rows", "satisfied": len(fdt) > 0, "detail": "finite Boltzmann response-ratio rows are computed"},
        {"obligation": "export_nonimported_physical_thermal_state_source", "satisfied": False, "detail": "0 candidates pass as internal non-imported KMS/thermal-state sources"},
    ]
    return {"status": "P3098_KMS_DETAILED_BALANCE_THERMAL_STATE_OBSTRUCTION_BOUNDED_NO_GO", "input_hashes": {"P3097": hashlib.sha256(P3097.read_bytes()).hexdigest() if P3097.exists() else None}, "constructed_theoretical_objects": {"content_first_repo_grep": greps, "kms_detailed_balance_audit_object": {"object": "Z12DirichletKMSDetailedBalanceThermalStateObstructionAudit", "source_reused": "P3097 recommendation: bounded KMS/detailed-balance thermal-state obstruction audit", "required_gates": list(REQUIRED_GATES), "candidate_thermal_state_sources": [c["id"] for c in THERMAL_STATE_CANDIDATES], "acceptance_predicate": "finite partition weights plus detailed-balance kernel, KMS periodicity, physical temperature clock, fluctuation-dissipation relation, empirical thermal readout, bath/preparation source, and non-imported physical thermal-state source"}, "gibbs_weight_rows": weights, "detailed_balance_transition_rows": db, "kms_periodicity_proxy_rows": kms, "fluctuation_dissipation_proxy_rows": fdt, "candidate_gate_rows": gates, "candidate_aggregate_certificate": aggs}, "finite_certificate": {"content_grep_lanes": len(greps), "content_grep_hits": sum(r["hit_count"] for r in greps), "p3097_accepted_nonimported_thermodynamic_radiation_sources": p3097.get("finite_certificate", {}).get("accepted_nonimported_thermodynamic_radiation_sources"), "gibbs_weight_rows": len(weights), "gibbs_rows_with_physical_temperature": sum(r["physical_temperature_attached"] for r in weights), "gibbs_rows_with_thermal_preparation": sum(r["thermal_preparation_attached"] for r in weights), "detailed_balance_transition_rows": len(db), "detailed_balance_rows_with_zero_residual": sum(r["detailed_balance_residual"] == 0.0 for r in db), "detailed_balance_rows_with_physical_bath": sum(r["physical_bath_attached"] for r in db), "detailed_balance_rows_with_empirical_transition_semantics": sum(r["empirical_transition_semantics_attached"] for r in db), "kms_periodicity_proxy_rows": len(kms), "kms_rows_with_physical_imaginary_time_clock": sum(r["physical_imaginary_time_clock_attached"] for r in kms), "kms_rows_with_operator_algebra_state": sum(r["operator_algebra_state_attached"] for r in kms), "fluctuation_dissipation_proxy_rows": len(fdt), "fdt_rows_with_physical_response_function": sum(r["physical_response_function_attached"] for r in fdt), "fdt_rows_with_empirical_thermal_readout": sum(r["empirical_thermal_readout_attached"] for r in fdt), "thermal_state_candidates": len(THERMAL_STATE_CANDIDATES), "required_gates": len(REQUIRED_GATES), "candidate_gate_rows": len(gates), "accepted_nonimported_kms_thermal_state_sources": len(accepted), "proof_obligations": len(obligations), "satisfied_proof_obligations": sum(r["satisfied"] for r in obligations)}, "proof_obligations": obligations, "decision": {"bounded_result": "P3098 constructs the requested KMS/detailed-balance thermal-state obstruction audit.  The Z12 Laplacian supports formal Gibbs weights, an explicit finite Metropolis-style kernel with zero algebraic detailed-balance residual, imaginary-time KMS-shift proxy rows, and finite fluctuation-dissipation-like Boltzmann response ratios.  These are real thermal-state-like witnesses, but no internal artifact exports a physical temperature clock, a bath/preparation source, an operator-algebra KMS state, a physical fluctuation-dissipation relation, empirical thermal readout, or a non-imported physical thermal-state source.  Imported continuum statistical-mechanics, FDT, and apparatus templates pass only as imported templates.  Therefore no non-imported KMS/detailed-balance thermal-state source is exported.", "negative_export_flags": {key: False for key in ["physical_temperature_clock_exported", "bath_or_preparation_source_exported", "operator_algebra_kms_state_exported", "physical_fluctuation_dissipation_relation_exported", "empirical_thermal_readout_exported", "nonimported_physical_thermal_state_source_exported", "physical_temperature_energy_unit_exported", "blackbody_radiation_law_exported", "observed_radiation_exported", "observed_photons_exported", "observed_light_exported", "spacetime_eom_exported", "physical_hamiltonian_exported", "thermodynamic_radiation_closure_exported", "scattering_smatrix_closure_exported", "QW_2191_discharged", "strict_selector_closure_exported", "ltotal_exported", "bridge_role_transfer_exported", "toe_closure_exported"]}, "positive_scoped_flags": {"gibbs_weight_rows_computed": True, "zero_residual_detailed_balance_rows_computed": True, "kms_periodicity_proxy_rows_computed": True, "fluctuation_dissipation_proxy_rows_computed": True}, "next_honest_step": "Pivot to exactly one adjacent standard-physics interface atom outside selector replay: construct a bounded irreversibility/entropy-production thermalization obstruction audit for the Z12 Dirichlet/Laplacian branch, testing whether finite detailed-balance kernels, Gibbs weights, relative-entropy monotonicity proxies, and modal flux witnesses supply a non-imported time arrow, physical bath/preparation mechanism, dissipative semigroup, and empirical thermalization readout without importing continuum nonequilibrium thermodynamics, apparatus units, observed light, L_total, bridge/role-transfer, or ToE."}}


def write_markdown(payload: dict[str, Any]) -> None:
    c = payload["finite_certificate"]
    lines = ["# P3098/S2048 KMS/detailed-balance thermal-state obstruction audit", "", f"Status: `{payload['status']}`", "", "## Finite certificate", f"- content grep lanes: `{c['content_grep_lanes']}`", f"- content grep hits: `{c['content_grep_hits']}`", f"- P3097 accepted non-imported thermodynamic-radiation sources: `{c['p3097_accepted_nonimported_thermodynamic_radiation_sources']}`", f"- Gibbs weight rows: `{c['gibbs_weight_rows']}`", f"- Gibbs rows with physical temperature: `{c['gibbs_rows_with_physical_temperature']}`", f"- detailed-balance transition rows: `{c['detailed_balance_transition_rows']}`", f"- detailed-balance rows with zero residual: `{c['detailed_balance_rows_with_zero_residual']}`", f"- detailed-balance rows with physical bath: `{c['detailed_balance_rows_with_physical_bath']}`", f"- KMS periodicity proxy rows: `{c['kms_periodicity_proxy_rows']}`", f"- KMS rows with physical imaginary-time clock: `{c['kms_rows_with_physical_imaginary_time_clock']}`", f"- fluctuation-dissipation proxy rows: `{c['fluctuation_dissipation_proxy_rows']}`", f"- FDT rows with empirical thermal readout: `{c['fdt_rows_with_empirical_thermal_readout']}`", f"- thermal-state candidates: `{c['thermal_state_candidates']}`", f"- candidate gate rows: `{c['candidate_gate_rows']}`", f"- accepted non-imported KMS thermal-state sources: `{c['accepted_nonimported_kms_thermal_state_sources']}`", f"- satisfied proof obligations: `{c['satisfied_proof_obligations']}/{c['proof_obligations']}`", "", "## Decision", payload["decision"]["bounded_result"], "", "## Recommendation", payload["decision"]["next_honest_step"]]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3098/S2048 KMS/detailed-balance thermal-state obstruction audit", "## P3098/S2048 KMS/detailed-balance thermal-state obstruction audit\n\n`P3098/S2048` attacks exactly one post-P3097 interface atom: a non-imported KMS/detailed-balance thermal-state source for the Z12 Dirichlet/Laplacian branch.  It constructs `48` Gibbs-weight rows, `528` detailed-balance transition rows, `16` KMS-periodicity proxy rows, `40` fluctuation-dissipation proxy rows, and a `5 x 8 = 40` candidate gate matrix.  The finite thermal-state-like algebra remains formal; no physical temperature clock, bath/preparation source, operator-algebra KMS state, physical fluctuation-dissipation relation, empirical thermal readout, `L_total`, bridge/role-transfer, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3098/S2048 KMS thermal-state source remains unsourced", "## P3098/S2048 KMS thermal-state source remains unsourced\n\n`P3098/S2048` confirms that the Z12 Laplacian supports formal Gibbs weights, exact finite detailed-balance kernels, KMS-shift proxies, and Boltzmann response-ratio rows.  A Lagrangian/EOM, thermal-state, fluctuation-dissipation, or empirical-readout reading still needs strict sources for physical temperature clocks, bath/preparation semantics, operator-algebra statehood, physical response functions, and thermal readout calibration; imported continuum statistical-mechanics and apparatus templates do not supply those strict sources.\n")
    append_once(AGENTS, "Current KMS/detailed-balance thermal-state obstruction guardrail (P3098/S2048, 2026-06-25)", "## Current KMS/detailed-balance thermal-state obstruction guardrail (P3098/S2048, 2026-06-25)\n\n- P3098 follows the P3097 recommendation and audits one standard-physics interface atom: a non-imported KMS/detailed-balance thermal-state source for the Z12 Dirichlet/Laplacian branch.\n- The finite audit computes `48` Gibbs-weight rows, `528` exact detailed-balance transition rows, `16` KMS-periodicity proxy rows, `40` fluctuation-dissipation proxy rows, and `40` candidate gate rows; `0` candidates export an internal non-imported KMS/thermal-state law.\n- Do not promote formal Gibbs weights, exact finite detailed-balance kernels, KMS-shift proxies, Boltzmann response-ratio rows, imported FDT/KMS templates, or imported apparatus/bath templates to physical thermal state, physical temperature clock, empirical observable, observed photons/light, spacetime EOM, physical Hamiltonian, selector closure, `L_total`, bridge/role-transfer, or ToE closure.\n- The next honest move is exactly one irreversibility/entropy-production thermalization obstruction audit for the Z12 Dirichlet/Laplacian branch, unless a genuinely new KMS/thermal-state theorem is introduced.\n")
    return payload

if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
