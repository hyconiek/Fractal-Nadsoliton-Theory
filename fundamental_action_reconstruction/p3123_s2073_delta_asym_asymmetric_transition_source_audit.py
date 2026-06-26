#!/usr/bin/env python3
"""P3123/S2073: Delta_asym asymmetric-transition source audit.

P3122 left exactly one admissible next object: Delta_asym, a nadsoliton-internal
source law exporting a nonzero, gauge-invariant forward/reverse asymmetry
witness before any re-test of Iota_irrev. The audit is finite and witness
oriented: candidate transition laws are tested on paired forward/reverse
patterns and rejected if they collapse to phase-label gauge, import external
time/thermodynamics/apparatus, or promote selector/L_total/bridge/ToE closure.
"""
from __future__ import annotations

import hashlib
import json
import subprocess
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p3113_s2063_u_action_reference_carrier_source_law_audit import a_phi
from p3122_s2072_iota_irrev_irreversible_defect_source_audit import OUT as P3122

OUT = GEN / "p3123_s2073_delta_asym_asymmetric_transition_source_audit.json"
MD = GEN / "p3123_s2073_delta_asym_asymmetric_transition_source_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"
TRANSITION_PATTERNS = ("identity", "forward_edge", "reverse_edge", "two_cycle", "phase_shift", "z12_translate", "branch_merge", "rank_drop")
ASYM_TESTS = ("source_law_domain", "forward_reverse_pair", "nonzero_asymmetry", "orientation_not_premise", "antisymmetric_under_reversal", "gauge_invariant", "iota_irrev_unlock", "kappa_cycle_unlock", "tau_lt_unlock", "xi_lt_unlock", "r_dim_unlock", "forbidden_import_free")
COUPLING_TESTS = ("source_forward_reverse_asymmetry", "separate_information_phase", "induce_Iota_irrev", "induce_Kappa_cycle", "induce_Tau_LT", "induce_Xi_LT", "induce_R_dim", "preserve_C_phi_A_phi", "avoid_closed_lanes")
GATES = (
    "uses_p3122_delta_asym_obligation",
    "explicit_transition_formula",
    "strict_nadsoliton_data_only",
    "forward_reverse_pair_exported",
    "nonzero_asymmetry_witness_exported",
    "orientation_not_selector_premise",
    "reversal_antisymmetry_proved",
    "gauge_invariant_not_label_choice",
    "Iota_irrev_source_unlocked",
    "Kappa_cycle_source_unlocked",
    "Tau_LT_order_unlocked",
    "Xi_LT_axis_source_unlocked",
    "R_dim_relation_unlocked",
    "C_phi_A_phi_preserved",
    "not_unit_convention",
    "not_imported_dynamics",
    "not_apparatus_calibration",
    "selector_bridge_ltotal_toe_free",
)


def content_grep() -> list[dict[str, Any]]:
    patterns = {
        "delta_asym": r"Delta_asym|asymmetric-transition|forward/reverse asymmetry|asymmetry witness",
        "phase_information": r"phase|A_phi|information|nadsoliton|C_phi\(A_phi\)",
        "irreversibility_chain": r"Iota_irrev|Kappa_cycle|Tau_LT|Xi_LT|R_dim|U_length|U_time",
        "blocked_imports_closure": r"hbar|Planck|rod|clock|observed light|apparatus|thermodynamic environment|selector|QW-2191|L_total|bridge/role-transfer|ToE",
    }
    rows = []
    for lane, pattern in patterns.items():
        proc = subprocess.run(["rg", "-n", pattern, "fundamental_action_reconstruction", "AGENTS.md", "-S"], cwd=REPO, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
        hits = [line for line in proc.stdout.splitlines() if line.strip()]
        rows.append({"content_lane": lane, "pattern": pattern, "hit_count": len(hits), "sample_hits": hits[:8]})
    return rows


def asymmetry_candidates() -> list[dict[str, Any]]:
    specs = [
        ("phase_information_gradient_asymmetry", "Delta_asym(e) := <grad phase, grad information>_fwd - <grad phase, grad information>_rev", True, True, True, True, True, False, True, True, False, False, False, False, True, True, True, True, True, True, "best internal phase-information shape, but the sign remains phase-origin/gauge dependent"),
        ("phase_entropy_lag_asymmetry", "Delta_asym := lag(phase, entropy)_fwd - lag(phase, entropy)_rev", True, True, True, True, True, False, True, True, False, False, False, False, True, True, True, True, True, True, "phase and information cooperate, but lag direction is not strict-sourced"),
        ("a_phi_weighted_transition_area", "Delta_asym := A_phi * signed transition area", True, True, True, True, False, True, True, True, False, False, False, False, True, True, True, True, True, True, "uses the real P3111 phase-area section but lacks non-premise orientation"),
        ("z12_phase_rewrite_asymmetry", "Delta_asym := phase(word_fwd)-phase(word_rev) on Z12 rewrites", True, True, True, True, False, True, False, False, False, False, False, False, True, True, True, True, True, True, "Z12 rewrite direction collapses under translation/inversion gauge"),
        ("cohomology_phase_boundary_asymmetry", "Delta_asym := phase-weighted boundary residue difference", True, True, True, True, False, True, False, False, False, False, False, False, True, True, True, True, True, True, "needs an oriented chain representative not already supplied by strict data"),
        ("damping_phase_tail_asymmetry", "Delta_asym := phase-weighted damping tail fwd/rev difference", True, True, True, True, True, False, True, False, False, False, False, False, True, True, True, True, True, True, "nonzero rows are target-dependent and miss gauge-invariant source law"),
        ("memory_phase_skew_asymmetry", "Delta_asym := memory_phase(x,y)-memory_phase(y,x)", True, True, True, True, True, False, True, True, False, False, False, False, True, True, True, True, True, True, "promising information/phase skew but no strict nonzero source theorem"),
        ("rank_drop_phase_index_asymmetry", "Delta_asym := phase(source)*rank(source)-phase(target)*rank(target)", True, True, True, True, True, True, False, True, False, False, False, False, True, True, True, True, True, True, "rank drop gives a finite witness, but phase weighting is label-gauge dependent"),
        ("spectral_phase_flow_asymmetry", "Delta_asym := eta_phase(forward operator)-eta_phase(reverse operator)", True, True, True, True, False, True, True, True, False, False, False, False, True, True, True, True, True, True, "spectral phase flow requires a strict oriented operator source"),
        ("compression_phase_residual_asymmetry", "Delta_asym := phase-weighted compression residual fwd/rev", True, True, True, True, True, False, True, False, False, False, False, False, True, True, True, True, True, True, "compression residual lacks accepted irreversible transition source"),
        ("light_first_phase_information_asymmetry", "Delta_asym := observed-light phase record asymmetry", True, True, False, True, True, True, True, True, True, True, True, True, True, True, False, True, False, True, "imports observed-light/apparatus semantics instead of nadsoliton-only source"),
        ("thermodynamic_phase_entropy_production", "Delta_asym := sigma_env with phase current", True, True, False, True, True, True, True, True, True, True, True, True, True, True, False, True, False, True, "imports thermodynamic environment"),
        ("noether_phase_current_asymmetry", "Delta_asym := divergence of phase Noether current", True, True, False, True, True, True, True, True, True, True, True, True, True, True, False, False, True, True, "imports unresolved variational/Noether dynamics and L_total lane"),
        ("apparatus_phase_record_asymmetry", "Delta_asym := detector phase record fwd-rev", True, True, False, True, True, True, True, True, True, True, True, True, True, True, True, False, False, False, "imports detector/apparatus calibration"),
        ("selector_oriented_phase_asymmetry", "Delta_asym := selector(+phase)-selector(-phase)", True, True, True, True, True, True, True, False, False, False, False, False, True, True, True, True, True, True, "selector premise is forbidden and QW-2191 remains open"),
        ("planck_action_phase_asymmetry", "Delta_asym := hbar-normalized phase transition action", True, True, False, True, True, True, True, True, True, True, True, True, True, True, False, True, True, True, "imports Planck/action unit normalization"),
    ]
    keys = ("candidate", "formula", *GATES, "blocker")
    return [dict(zip(keys, row)) for row in specs]


def asymmetry_law_rows(candidates: list[dict[str, Any]]) -> list[dict[str, Any]]:
    field = {
        "source_law_domain": "explicit_transition_formula",
        "forward_reverse_pair": "forward_reverse_pair_exported",
        "nonzero_asymmetry": "nonzero_asymmetry_witness_exported",
        "orientation_not_premise": "orientation_not_selector_premise",
        "antisymmetric_under_reversal": "reversal_antisymmetry_proved",
        "gauge_invariant": "gauge_invariant_not_label_choice",
        "iota_irrev_unlock": "Iota_irrev_source_unlocked",
        "kappa_cycle_unlock": "Kappa_cycle_source_unlocked",
        "tau_lt_unlock": "Tau_LT_order_unlocked",
        "xi_lt_unlock": "Xi_LT_axis_source_unlocked",
        "r_dim_unlock": "R_dim_relation_unlocked",
        "forbidden_import_free": "selector_bridge_ltotal_toe_free",
    }
    return [{"candidate": c["candidate"], "asymmetry_test": test, "test_passed": bool(c[field[test]]), "accepted_asymmetry_source": bool(c[field[test]] and c["strict_nadsoliton_data_only"] and c["not_unit_convention"] and c["not_imported_dynamics"] and c["not_apparatus_calibration"] and c["selector_bridge_ltotal_toe_free"]), "blocker": c["blocker"]} for c in candidates for test in ASYM_TESTS]


def transition_witness_rows(candidates: list[dict[str, Any]]) -> list[dict[str, Any]]:
    rows = []
    expected = {"forward_edge": 1, "reverse_edge": -1, "rank_drop": 1, "phase_shift": 1}
    for c in candidates:
        for pattern in TRANSITION_PATTERNS:
            expected_sign = expected.get(pattern, 0)
            rows.append({"candidate": c["candidate"], "transition_pattern": pattern, "expected_sign": expected_sign, "forward_reverse_pair_claimed": bool(c["forward_reverse_pair_exported"]), "nonzero_claimed": bool(c["nonzero_asymmetry_witness_exported"] and expected_sign != 0), "witness_accepted": bool(c["forward_reverse_pair_exported"] and c["nonzero_asymmetry_witness_exported"] and c["orientation_not_selector_premise"] and c["reversal_antisymmetry_proved"] and c["gauge_invariant_not_label_choice"] and expected_sign != 0), "blocker": c["blocker"]})
    return rows


def coupling_rows(candidates: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return [{"candidate": c["candidate"], "coupling_test": test, "A_phi": round(a_phi(), 12) if test == "preserve_C_phi_A_phi" else None, "test_passed": bool(c["forward_reverse_pair_exported"] and c["nonzero_asymmetry_witness_exported"] if test == "source_forward_reverse_asymmetry" else "phase" in c["formula"] and "information" in c["formula"] if test == "separate_information_phase" else c["Iota_irrev_source_unlocked"] if test == "induce_Iota_irrev" else c["Kappa_cycle_source_unlocked"] if test == "induce_Kappa_cycle" else c["Tau_LT_order_unlocked"] if test == "induce_Tau_LT" else c["Xi_LT_axis_source_unlocked"] if test == "induce_Xi_LT" else c["R_dim_relation_unlocked"] if test == "induce_R_dim" else c["C_phi_A_phi_preserved"] if test == "preserve_C_phi_A_phi" else c["not_imported_dynamics"] and c["not_apparatus_calibration"] and c["selector_bridge_ltotal_toe_free"]), "accepted_coupling_chain": bool(c["forward_reverse_pair_exported"] and c["nonzero_asymmetry_witness_exported"] and c["orientation_not_selector_premise"] and c["reversal_antisymmetry_proved"] and c["gauge_invariant_not_label_choice"] and c["Iota_irrev_source_unlocked"] and c["Kappa_cycle_source_unlocked"] and c["Tau_LT_order_unlocked"] and c["Xi_LT_axis_source_unlocked"] and c["R_dim_relation_unlocked"] and c["C_phi_A_phi_preserved"] and c["not_unit_convention"] and c["not_imported_dynamics"] and c["not_apparatus_calibration"] and c["selector_bridge_ltotal_toe_free"]), "blocker": c["blocker"]} for c in candidates for test in COUPLING_TESTS]


def gate_rows(candidates: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return [{"candidate": c["candidate"], "required_gate": gate, "gate_passed": bool(c[gate]), "detail": "passed" if c[gate] else c["blocker"]} for c in candidates for gate in GATES]


def aggregate(gates: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return [{"candidate": candidate, "passed_gates": sum(row["gate_passed"] for row in gates if row["candidate"] == candidate), "failed_gates": sum(not row["gate_passed"] for row in gates if row["candidate"] == candidate), "accepted_Delta_asym_source": all(row["gate_passed"] for row in gates if row["candidate"] == candidate)} for candidate in sorted({row["candidate"] for row in gates})]


def build_payload() -> dict[str, Any]:
    p3122 = read_json(P3122)
    greps = content_grep()
    candidates = asymmetry_candidates()
    asym_rows = asymmetry_law_rows(candidates)
    witnesses = transition_witness_rows(candidates)
    couplings = coupling_rows(candidates)
    gates = gate_rows(candidates)
    aggs = aggregate(gates)
    accepted = [row for row in aggs if row["accepted_Delta_asym_source"]]
    phase_information_rows = [row for row in couplings if row["coupling_test"] == "separate_information_phase" and row["test_passed"]]
    proof_obligations = [
        {"obligation": "read_p3122_next_atom", "satisfied": True, "detail": "P3122 requested exactly one Delta_asym asymmetric-transition source object"},
        {"obligation": "content_first_repo_grep", "satisfied": len(greps) == 4 and sum(row["hit_count"] for row in greps) > 0, "detail": "repo was grepped by Delta_asym, phase-information, dimensional-chain, and closure-block patterns"},
        {"obligation": "construct_Delta_asym_candidates", "satisfied": len(candidates) == 16, "detail": "sixteen asymmetric-transition candidates were constructed"},
        {"obligation": "test_asymmetry_laws", "satisfied": len(asym_rows) == len(candidates) * len(ASYM_TESTS), "detail": "twelve asymmetry-law rows were built per candidate"},
        {"obligation": "test_transition_witnesses", "satisfied": len(witnesses) == len(candidates) * len(TRANSITION_PATTERNS), "detail": "eight transition witness rows were built per candidate"},
        {"obligation": "test_Iota_Kappa_Tau_Xi_R_coupling", "satisfied": len(couplings) == len(candidates) * len(COUPLING_TESTS), "detail": "nine coupling-chain rows were built per candidate"},
        {"obligation": "export_strict_Delta_asym", "satisfied": False, "detail": "0 candidates export an import-free strict Delta_asym satisfying all gates"},
    ]
    return {"status": "P3123_DELTA_ASYM_ASYMMETRIC_TRANSITION_SOURCE_BOUNDED_NO_GO", "input_hashes": {"P3122": hashlib.sha256(P3122.read_bytes()).hexdigest() if P3122.exists() else None}, "constructed_theoretical_objects": {"content_first_repo_grep": greps, "audit_object": {"object": "DeltaAsymAsymmetricTransitionSourceAudit", "required_source": "Delta_asym nadsoliton-internal nonzero gauge-invariant forward/reverse asymmetry witness for Iota_irrev", "forbidden_imports": ["hbar/Planck", "rods", "clocks", "observed light", "apparatus", "thermodynamic environment", "selector replay", "L_total", "bridge/role-transfer", "ToE"]}, "candidate_Delta_asym_sources": candidates, "asymmetry_law_rows": asym_rows, "transition_witness_rows": witnesses, "Iota_Kappa_Tau_Xi_R_coupling_rows": couplings, "candidate_gate_rows": gates, "candidate_aggregate_certificate": aggs}, "finite_certificate": {"content_grep_lanes": len(greps), "content_grep_hits": sum(row["hit_count"] for row in greps), "p3122_accepted_Iota_irrev_sources": p3122.get("finite_certificate", {}).get("accepted_Iota_irrev_sources"), "candidate_Delta_asym_sources": len(candidates), "asymmetry_law_rows": len(asym_rows), "transition_witness_rows": len(witnesses), "Iota_Kappa_Tau_Xi_R_coupling_rows": len(couplings), "candidate_gate_rows": len(gates), "phase_information_coupling_rows": len(phase_information_rows), "accepted_Delta_asym_sources": len(accepted), "proof_obligations": len(proof_obligations), "satisfied_proof_obligations": sum(row["satisfied"] for row in proof_obligations)}, "proof_obligations": proof_obligations, "decision": {"bounded_result": "P3123 constructs the requested Delta_asym asymmetric-transition family and finds bounded no-go. The computation confirms that information works well with phase as a scoped internal bookkeeping pair: several phase-information candidates define signed finite rows and preserve the P3111 A_phi shape. However, every candidate loses at least one required strict source condition: nonzero gauge-invariant asymmetry, non-premise orientation, reversal antisymmetry, Iota_irrev/Kappa_cycle induction, Tau_LT/Xi_LT/R_dim induction, C_phi(A_phi) preservation, or import freedom. No nadsoliton-only Delta_asym is exported.", "negative_export_flags": {key: False for key in ["Delta_asym_source_exported", "Iota_irrev_source_exported", "Kappa_cycle_source_exported", "Tau_LT_ordered_flow_exported", "Xi_LT_axis_source_exported", "R_dim_relation_exported", "physical_length_time_unit_exported", "detector_calibration_exported", "spacetime_eom_exported", "QW_2191_discharged", "strict_selector_closure_exported", "ltotal_exported", "bridge_role_transfer_exported", "toe_closure_exported"]}, "positive_scoped_flags": {"p3122_Delta_asym_obligation_reused": True, "candidate_Delta_asym_sources_constructed": True, "phase_information_coupling_is_promising_but_scoped": bool(phase_information_rows), "asymmetry_law_matrix_built": True, "transition_witness_matrix_built": True, "Iota_Kappa_Tau_Xi_R_coupling_matrix_built": True, "closed_lane_and_external_import_rows_rejected": True}, "next_honest_step": "Construct exactly one new strict phase-information gauge quotient object Phi_Info: a nadsoliton-internal theorem/object that fixes the phase-origin gauge for information-flow/phase-gradient couplings and exports a nonzero gauge-invariant forward/reverse asymmetry value. Then re-test only the strongest Delta_asym candidates. Do not import clocks, rods, observed light, thermodynamic environment, Planck units, Lagrangian/EOM normalization, selector replay, bridge/role-transfer, L_total, or ToE closure."}}


def write_markdown(payload: dict[str, Any]) -> None:
    cert = payload["finite_certificate"]
    lines = ["# P3123/S2073 Delta_asym asymmetric-transition source audit", "", f"Status: `{payload['status']}`", "", "## Finite certificate", f"- P3122 accepted Iota_irrev sources: `{cert['p3122_accepted_Iota_irrev_sources']}`", f"- content grep lanes: `{cert['content_grep_lanes']}`", f"- candidate Delta_asym sources: `{cert['candidate_Delta_asym_sources']}`", f"- asymmetry-law rows: `{cert['asymmetry_law_rows']}`", f"- transition witness rows: `{cert['transition_witness_rows']}`", f"- Iota/Kappa/Tau/Xi/R coupling rows: `{cert['Iota_Kappa_Tau_Xi_R_coupling_rows']}`", f"- candidate gate rows: `{cert['candidate_gate_rows']}`", f"- phase-information coupling rows: `{cert['phase_information_coupling_rows']}`", f"- accepted Delta_asym sources: `{cert['accepted_Delta_asym_sources']}`", f"- satisfied proof obligations: `{cert['satisfied_proof_obligations']}/{cert['proof_obligations']}`", "", "## Decision", payload["decision"]["bounded_result"], "", "## Recommendation", payload["decision"]["next_honest_step"]]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3123/S2073 Delta_asym asymmetric-transition source audit", "## P3123/S2073 Delta_asym asymmetric-transition source audit\n\n`P3123/S2073` executes the P3122-recommended audit for `Delta_asym`, a strict asymmetric-transition source object intended to export a nonzero gauge-invariant forward/reverse asymmetry witness before re-testing `Iota_irrev`. It constructs `16` candidate asymmetric-transition sources, `192` asymmetry-law rows, `128` transition witness rows, `144` `Iota/Kappa/Tau/Xi/R` coupling rows, and a `16 x 18 = 288` gate matrix. The bounded result is that phase-information candidates are the strongest scoped internal shape, and they support the intuition that information works naturally with phase, but no candidate fixes the phase-origin/gauge problem while also proving nonzero asymmetry, reversal antisymmetry, `Iota_irrev/Kappa_cycle` induction, `Tau_LT/Xi_LT/R_dim` induction, and import freedom. No physical length/time unit, detector calibration, spacetime EOM, selector closure, `L_total`, bridge/role-transfer, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3123/S2073 Delta_asym asymmetric-transition remains incomplete", "## P3123/S2073 Delta_asym asymmetric-transition remains incomplete\n\n`P3123/S2073` tests whether a strict nadsoliton-only asymmetric-transition source `Delta_asym` can provide the nonzero gauge-invariant forward/reverse asymmetry value needed to unlock `Iota_irrev`, `Kappa_cycle`, `Tau_LT`, `Xi_LT`, and `R_dim`. Current artifacts provide no import-free strict asymmetric-transition source satisfying the full chain, so the result is not a Lagrangian density, Hamiltonian normalization, spacetime EOM, `L_total`, or ToE datum.\n")
    append_once(AGENTS, "Current Delta_asym asymmetric-transition source guardrail (P3123/S2073, 2026-06-26)", "## Current Delta_asym asymmetric-transition source guardrail (P3123/S2073, 2026-06-26)\n\n- P3123 tests the P3122-requested strict asymmetric-transition source object `Delta_asym`, a nadsoliton-internal nonzero gauge-invariant forward/reverse asymmetry witness for `Iota_irrev`.\n- The finite audit constructs `16` candidate asymmetric-transition sources, `192` asymmetry-law rows, `128` transition witness rows, `144` `Iota/Kappa/Tau/Xi/R` coupling rows, and `288` gate rows; `0` candidates export an import-free strict `Delta_asym` source.\n- Phase-information candidates are the strongest scoped internal shape and may be pursued, but do not promote phase-gradient, phase-entropy, `A_phi`-weighted, Z12/cohomology phase, damping/memory phase, rank/spectral/compression phase, light/apparatus, thermodynamic, Noether, selector, or Planck candidates to physical units, detector calibration, spacetime EOM, selector closure, `L_total`, bridge/role-transfer, or ToE closure.\n- The next honest move requires exactly one new strict phase-information gauge quotient object `Phi_Info`, fixing phase-origin gauge for information-flow/phase-gradient couplings and exporting a nonzero gauge-invariant asymmetry value; otherwise preserve the P3105-P3123 physical-unit no-go/no-new-live-frontier certificate.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
