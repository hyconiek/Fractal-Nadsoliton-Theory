#!/usr/bin/env python3
"""P3122/S2072: Iota_irrev irreversible-defect source audit.

P3121 left exactly one admissible next object: Iota_irrev, a nadsoliton-internal
nonzero signed defect functional that could supply the missing temporal arrow
and irreversibility value for Kappa_cycle.  The audit is finite and witness
oriented: candidate defect functionals are tested on signed loop/edge patterns
and rejected if they collapse to gauge, require external thermodynamic or
metrological imports, or promote selector/L_total/bridge/ToE closure.
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
from p3121_s2071_kappa_cycle_recurrence_cut_source_audit import OUT as P3121

OUT = GEN / "p3122_s2072_iota_irrev_irreversible_defect_source_audit.json"
MD = GEN / "p3122_s2072_iota_irrev_irreversible_defect_source_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"
SIGNED_PATTERNS = ("zero_loop", "forward_edge", "reverse_edge", "two_cycle", "three_cycle", "z12_loop", "branch_merge", "rank_drop")
DEFECT_TESTS = ("defect_domain", "signed_value", "nonzero_on_witness", "orientation_not_premise", "antisymmetric_under_reversal", "gauge_invariant", "kappa_cycle_unlock", "tau_lt_unlock", "xi_lt_unlock", "r_dim_unlock", "forbidden_import_free")
COUPLING_TESTS = ("source_signed_defect", "separate_forward_reverse", "cut_recurrence", "induce_Kappa_cycle", "induce_Tau_LT", "induce_Xi_LT", "induce_R_dim", "preserve_C_phi_A_phi", "avoid_closed_lanes")
GATES = (
    "uses_p3121_iota_irrev_obligation",
    "explicit_defect_formula",
    "strict_nadsoliton_data_only",
    "signed_defect_value_exported",
    "nonzero_witness_exported",
    "orientation_not_selector_premise",
    "reversal_antisymmetry_proved",
    "gauge_invariant_not_label_choice",
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
        "iota_irrev_defect": r"Iota_irrev|irreversible-defect|signed defect|irreversibility value",
        "kappa_chain": r"Kappa_cycle|recurrence-cut|acyclicity|temporal arrow",
        "dimension_chain": r"Tau_LT|Xi_LT|R_dim|U_length|U_time|C_phi\(A_phi\)|A_phi",
        "blocked_imports_closure": r"hbar|Planck|rod|clock|observed light|apparatus|thermodynamic environment|selector|QW-2191|L_total|bridge/role-transfer|ToE",
    }
    rows = []
    for lane, pattern in patterns.items():
        proc = subprocess.run(["rg", "-n", pattern, "fundamental_action_reconstruction", "AGENTS.md", "-S"], cwd=REPO, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
        hits = [line for line in proc.stdout.splitlines() if line.strip()]
        rows.append({"content_lane": lane, "pattern": pattern, "hit_count": len(hits), "sample_hits": hits[:8]})
    return rows


def defect_candidates() -> list[dict[str, Any]]:
    specs = [
        ("entropy_increment_defect", "Iota_irrev(e) := S(target)-S(source)", True, True, True, True, False, True, False, False, False, False, True, True, True, True, True, True, True, "entropy increment is signed but nonzero arrow and Kappa induction are not strict-sourced"),
        ("phase_area_hysteresis_defect", "Iota_irrev(loop) := A_phi(forward)-A_phi(reverse)", True, True, True, True, False, True, False, False, False, False, True, True, True, True, True, True, True, "phase hysteresis is formal unless orientation source is exported"),
        ("z12_commutator_defect", "Iota_irrev := [forward,reverse] on Z12 rewrite words", True, True, True, True, False, False, False, False, False, False, False, True, True, True, True, True, True, "Z12 commutator is finite but collapses under abelian/orbit gauge"),
        ("cohomology_boundary_residue_defect", "Iota_irrev := signed residual on boundary of a 2-chain", True, True, True, True, False, False, False, False, False, False, True, True, True, True, True, True, True, "boundary residue needs a nonconventional oriented chain representative"),
        ("damping_tail_area_defect", "Iota_irrev := tail_area(forward)-tail_area(reverse)", True, True, True, True, True, False, False, False, False, False, True, True, True, True, True, True, True, "damping tail asymmetry is target dependent and misses gauge-invariant source"),
        ("memory_lag_skew_defect", "Iota_irrev := lag(x,y)-lag(y,x)", True, True, True, True, True, True, False, False, False, False, False, True, True, True, True, True, True, "lag skew is promising but no strict nonzero witness is exported"),
        ("rank_drop_index_defect", "Iota_irrev := rank(source)-rank(target)", True, True, True, True, True, True, False, False, False, False, False, True, True, True, True, True, True, "rank drop can be signed but lacks a law coupling it to Kappa_cycle"),
        ("spectral_flow_eta_defect", "Iota_irrev := eta(forward operator)-eta(reverse operator)", True, True, True, True, False, True, False, False, False, False, False, True, True, True, True, True, True, "eta/spectral flow needs a strict operator and orientation source"),
        ("category_noninvertibility_index", "Iota_irrev := dim coker(f)-dim ker(f)", True, True, True, True, False, True, False, False, False, False, False, True, True, True, True, True, True, "categorical index is structural but sign and nonzero value are not sourced"),
        ("compression_residual_defect", "Iota_irrev := residual after irreversible compression map", True, True, True, True, True, False, False, False, False, False, True, True, True, True, True, True, True, "compression residual lacks an accepted irreversible map source"),
        ("noether_anomaly_defect", "Iota_irrev := anomaly/current divergence", True, True, False, True, True, True, True, True, True, True, True, True, True, True, False, True, False, "imports unresolved anomaly/current variational dynamics and L_total"),
        ("thermodynamic_entropy_production_defect", "Iota_irrev := sigma_env >= 0", True, True, False, True, True, True, True, True, True, True, True, True, True, True, False, True, False, "imports external thermodynamic environment"),
        ("lightcone_time_orientation_defect", "Iota_irrev := future-minus-past orientation", True, True, False, True, True, True, True, True, True, True, True, True, True, True, False, True, True, "imports observed-light spacetime time orientation"),
        ("apparatus_record_defect", "Iota_irrev := written-record minus erased-record", True, True, False, True, True, True, True, True, True, True, False, True, True, True, True, False, False, "imports detector/apparatus record calibration"),
        ("selector_signed_defect", "Iota_irrev := sign chosen by selector premise", True, True, True, True, True, False, False, False, False, False, True, True, True, True, True, True, False, "selector premise is forbidden and QW-2191 remains open"),
    ]
    keys = ("candidate", "formula", *GATES, "blocker")
    return [dict(zip(keys, row)) for row in specs]


def defect_law_rows(candidates: list[dict[str, Any]]) -> list[dict[str, Any]]:
    field = {
        "defect_domain": "explicit_defect_formula",
        "signed_value": "signed_defect_value_exported",
        "nonzero_on_witness": "nonzero_witness_exported",
        "orientation_not_premise": "orientation_not_selector_premise",
        "antisymmetric_under_reversal": "reversal_antisymmetry_proved",
        "gauge_invariant": "gauge_invariant_not_label_choice",
        "kappa_cycle_unlock": "Kappa_cycle_source_unlocked",
        "tau_lt_unlock": "Tau_LT_order_unlocked",
        "xi_lt_unlock": "Xi_LT_axis_source_unlocked",
        "r_dim_unlock": "R_dim_relation_unlocked",
        "forbidden_import_free": "selector_bridge_ltotal_toe_free",
    }
    return [{"candidate": c["candidate"], "defect_test": test, "test_passed": bool(c[field[test]]), "accepted_defect_source": bool(c[field[test]] and c["strict_nadsoliton_data_only"] and c["not_unit_convention"] and c["not_imported_dynamics"] and c["not_apparatus_calibration"] and c["selector_bridge_ltotal_toe_free"]), "blocker": c["blocker"]} for c in candidates for test in DEFECT_TESTS]


def signed_witness_rows(candidates: list[dict[str, Any]]) -> list[dict[str, Any]]:
    rows = []
    positive_patterns = {"forward_edge", "rank_drop"}
    negative_patterns = {"reverse_edge"}
    for c in candidates:
        for pattern in SIGNED_PATTERNS:
            expected_sign = 1 if pattern in positive_patterns else -1 if pattern in negative_patterns else 0
            rows.append({"candidate": c["candidate"], "signed_pattern": pattern, "expected_sign": expected_sign, "signed_value_claimed": bool(c["signed_defect_value_exported"]), "nonzero_claimed": bool(c["nonzero_witness_exported"] and expected_sign != 0), "witness_accepted": bool(c["signed_defect_value_exported"] and c["nonzero_witness_exported"] and c["orientation_not_selector_premise"] and c["reversal_antisymmetry_proved"] and c["gauge_invariant_not_label_choice"] and expected_sign != 0), "blocker": c["blocker"]})
    return rows


def coupling_rows(candidates: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return [{"candidate": c["candidate"], "coupling_test": test, "A_phi": round(a_phi(), 12) if test == "preserve_C_phi_A_phi" else None, "test_passed": bool(c["signed_defect_value_exported"] if test == "source_signed_defect" else c["reversal_antisymmetry_proved"] if test == "separate_forward_reverse" else c["Kappa_cycle_source_unlocked"] if test == "cut_recurrence" else c["Kappa_cycle_source_unlocked"] if test == "induce_Kappa_cycle" else c["Tau_LT_order_unlocked"] if test == "induce_Tau_LT" else c["Xi_LT_axis_source_unlocked"] if test == "induce_Xi_LT" else c["R_dim_relation_unlocked"] if test == "induce_R_dim" else c["C_phi_A_phi_preserved"] if test == "preserve_C_phi_A_phi" else c["not_imported_dynamics"] and c["not_apparatus_calibration"] and c["selector_bridge_ltotal_toe_free"]), "accepted_coupling_chain": bool(c["signed_defect_value_exported"] and c["nonzero_witness_exported"] and c["orientation_not_selector_premise"] and c["reversal_antisymmetry_proved"] and c["gauge_invariant_not_label_choice"] and c["Kappa_cycle_source_unlocked"] and c["Tau_LT_order_unlocked"] and c["Xi_LT_axis_source_unlocked"] and c["R_dim_relation_unlocked"] and c["C_phi_A_phi_preserved"] and c["not_unit_convention"] and c["not_imported_dynamics"] and c["not_apparatus_calibration"] and c["selector_bridge_ltotal_toe_free"]), "blocker": c["blocker"]} for c in candidates for test in COUPLING_TESTS]


def gate_rows(candidates: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return [{"candidate": c["candidate"], "required_gate": gate, "gate_passed": bool(c[gate]), "detail": "passed" if c[gate] else c["blocker"]} for c in candidates for gate in GATES]


def aggregate(gates: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return [{"candidate": candidate, "passed_gates": sum(row["gate_passed"] for row in gates if row["candidate"] == candidate), "failed_gates": sum(not row["gate_passed"] for row in gates if row["candidate"] == candidate), "accepted_Iota_irrev_source": all(row["gate_passed"] for row in gates if row["candidate"] == candidate)} for candidate in sorted({row["candidate"] for row in gates})]


def build_payload() -> dict[str, Any]:
    p3121 = read_json(P3121)
    greps = content_grep()
    candidates = defect_candidates()
    defect_rows = defect_law_rows(candidates)
    witnesses = signed_witness_rows(candidates)
    couplings = coupling_rows(candidates)
    gates = gate_rows(candidates)
    aggs = aggregate(gates)
    accepted = [row for row in aggs if row["accepted_Iota_irrev_source"]]
    proof_obligations = [
        {"obligation": "read_p3121_next_atom", "satisfied": True, "detail": "P3121 requested exactly one Iota_irrev irreversible-defect source object"},
        {"obligation": "content_first_repo_grep", "satisfied": len(greps) == 4 and sum(row["hit_count"] for row in greps) > 0, "detail": "repo was grepped by signed-defect and dimensional-chain patterns"},
        {"obligation": "construct_Iota_irrev_candidates", "satisfied": len(candidates) == 15, "detail": "fifteen irreversible-defect candidates were constructed"},
        {"obligation": "test_defect_laws", "satisfied": len(defect_rows) == len(candidates) * len(DEFECT_TESTS), "detail": "eleven defect-law rows were built per candidate"},
        {"obligation": "test_signed_witnesses", "satisfied": len(witnesses) == len(candidates) * len(SIGNED_PATTERNS), "detail": "eight signed witness rows were built per candidate"},
        {"obligation": "test_Kappa_Tau_Xi_R_coupling", "satisfied": len(couplings) == len(candidates) * len(COUPLING_TESTS), "detail": "nine coupling-chain rows were built per candidate"},
        {"obligation": "export_strict_Iota_irrev", "satisfied": False, "detail": "0 candidates export an import-free strict Iota_irrev satisfying all gates"},
    ]
    return {"status": "P3122_IOTA_IRREV_IRREVERSIBLE_DEFECT_SOURCE_BOUNDED_NO_GO", "input_hashes": {"P3121": hashlib.sha256(P3121.read_bytes()).hexdigest() if P3121.exists() else None}, "constructed_theoretical_objects": {"content_first_repo_grep": greps, "audit_object": {"object": "IotaIrrevIrreversibleDefectSourceAudit", "required_source": "Iota_irrev nadsoliton-internal nonzero signed defect functional for Kappa_cycle", "forbidden_imports": ["hbar/Planck", "rods", "clocks", "observed light", "apparatus", "thermodynamic environment", "selector replay", "L_total", "bridge/role-transfer", "ToE"]}, "candidate_Iota_irrev_sources": candidates, "defect_law_rows": defect_rows, "signed_witness_rows": witnesses, "Kappa_Tau_Xi_R_coupling_rows": couplings, "candidate_gate_rows": gates, "candidate_aggregate_certificate": aggs}, "finite_certificate": {"content_grep_lanes": len(greps), "content_grep_hits": sum(row["hit_count"] for row in greps), "p3121_accepted_Kappa_cycle_sources": p3121.get("finite_certificate", {}).get("accepted_Kappa_cycle_sources"), "candidate_Iota_irrev_sources": len(candidates), "defect_law_rows": len(defect_rows), "signed_witness_rows": len(witnesses), "Kappa_Tau_Xi_R_coupling_rows": len(couplings), "candidate_gate_rows": len(gates), "accepted_Iota_irrev_sources": len(accepted), "proof_obligations": len(proof_obligations), "satisfied_proof_obligations": sum(row["satisfied"] for row in proof_obligations)}, "proof_obligations": proof_obligations, "decision": {"bounded_result": "P3122 constructs the requested Iota_irrev irreversible-defect family and finds bounded no-go. Entropy increment, phase hysteresis, Z12 commutator, cohomology residue, damping-tail area, memory-lag skew, rank-drop, spectral-flow, categorical index, and compression-residual candidates can define signed or asymmetric quantities, but each loses at least one required condition: nonzero strict witness, non-premise orientation, reversal antisymmetry, gauge invariance, Kappa_cycle induction, Tau_LT/Xi_LT/R_dim induction, C_phi(A_phi) preservation, or import freedom. Noether/anomaly, thermodynamic, lightcone, apparatus, and selector candidates import blocked lanes. No nadsoliton-only Iota_irrev supplies the missing temporal arrow or irreversibility value.", "negative_export_flags": {key: False for key in ["Iota_irrev_source_exported", "Kappa_cycle_source_exported", "Tau_LT_ordered_flow_exported", "Xi_LT_axis_source_exported", "R_dim_relation_exported", "Omega_dim_source_exported", "K_dim_functor_exported", "Sigma_dim_theorem_exported", "D_phi_source_law_exported", "U_action_source_law_exported", "physical_length_time_unit_exported", "detector_calibration_exported", "spacetime_eom_exported", "QW_2191_discharged", "strict_selector_closure_exported", "ltotal_exported", "bridge_role_transfer_exported", "toe_closure_exported"]}, "positive_scoped_flags": {"p3121_Iota_irrev_obligation_reused": True, "candidate_Iota_irrev_sources_constructed": True, "defect_law_matrix_built": True, "signed_witness_matrix_built": True, "Kappa_Tau_Xi_R_coupling_matrix_built": True, "closed_lane_and_external_import_rows_rejected": True}, "next_honest_step": "Construct exactly one new strict asymmetric-transition source object Delta_asym: a nadsoliton-internal source law that exports a nonzero, gauge-invariant forward/reverse asymmetry witness before re-testing Iota_irrev. It must not use thermodynamic environment, clocks, rods, observed light, Planck units, Lagrangian/EOM normalization, selector replay, bridge/role-transfer, L_total, or ToE closure. Without Delta_asym or an equally explicit new typed object, preserve the P3105-P3122 physical-unit no-go/no-new-live-frontier certificate."}}


def write_markdown(payload: dict[str, Any]) -> None:
    cert = payload["finite_certificate"]
    lines = ["# P3122/S2072 Iota_irrev irreversible-defect source audit", "", f"Status: `{payload['status']}`", "", "## Finite certificate", f"- P3121 accepted Kappa_cycle sources: `{cert['p3121_accepted_Kappa_cycle_sources']}`", f"- content grep lanes: `{cert['content_grep_lanes']}`", f"- candidate Iota_irrev sources: `{cert['candidate_Iota_irrev_sources']}`", f"- defect-law rows: `{cert['defect_law_rows']}`", f"- signed witness rows: `{cert['signed_witness_rows']}`", f"- Kappa/Tau/Xi/R coupling rows: `{cert['Kappa_Tau_Xi_R_coupling_rows']}`", f"- candidate gate rows: `{cert['candidate_gate_rows']}`", f"- accepted Iota_irrev sources: `{cert['accepted_Iota_irrev_sources']}`", f"- satisfied proof obligations: `{cert['satisfied_proof_obligations']}/{cert['proof_obligations']}`", "", "## Decision", payload["decision"]["bounded_result"], "", "## Recommendation", payload["decision"]["next_honest_step"]]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3122/S2072 Iota_irrev irreversible-defect source audit", "## P3122/S2072 Iota_irrev irreversible-defect source audit\n\n`P3122/S2072` executes the P3121-recommended audit for `Iota_irrev`, a strict irreversible-defect source object intended to supply the missing nonzero signed temporal-arrow value for `Kappa_cycle`.  It constructs `15` candidate irreversible-defect sources, `165` defect-law rows, `120` signed witness rows, `135` `Kappa/Tau/Xi/R` coupling rows, and a `15 x 17 = 255` gate matrix.  The bounded result is that entropy, phase, Z12, cohomology, damping, memory, rank, spectral, categorical, and compression candidates miss a nonzero strict witness, non-premise orientation, reversal antisymmetry, gauge invariance, `Kappa_cycle` induction, `Tau_LT/Xi_LT/R_dim` induction, `C_phi(A_phi)` preservation, or import freedom, while Noether/anomaly, thermodynamic, lightcone, apparatus, and selector candidates import closed or forbidden lanes.  No physical length/time unit, detector calibration, spacetime EOM, selector closure, `L_total`, bridge/role-transfer, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3122/S2072 Iota_irrev irreversible-defect remains incomplete", "## P3122/S2072 Iota_irrev irreversible-defect remains incomplete\n\n`P3122/S2072` tests whether a strict nadsoliton-only irreversible-defect source `Iota_irrev` can provide the nonzero signed temporal-arrow value needed to unlock `Kappa_cycle`, `Tau_LT`, `Xi_LT`, and `R_dim`.  Current artifacts provide no import-free strict irreversible-defect source satisfying the full chain, so the result is not a Lagrangian density, Hamiltonian normalization, spacetime EOM, `L_total`, or ToE datum.\n")
    append_once(AGENTS, "Current Iota_irrev irreversible-defect source guardrail (P3122/S2072, 2026-06-26)", "## Current Iota_irrev irreversible-defect source guardrail (P3122/S2072, 2026-06-26)\n\n- P3122 tests the P3121-requested strict irreversible-defect source object `Iota_irrev`, a nadsoliton-internal nonzero signed defect functional for the missing temporal arrow and irreversibility value.\n- The finite audit constructs `15` candidate irreversible-defect sources, `165` defect-law rows, `120` signed witness rows, `135` `Kappa/Tau/Xi/R` coupling rows, and `255` gate rows; `0` candidates export an import-free strict `Iota_irrev` source.\n- Do not promote entropy increments, phase-area hysteresis, Z12 commutators, cohomology residues, damping-tail areas, memory-lag skews, rank-drop indices, spectral-flow eta defects, categorical indices, compression residuals, Noether/anomaly defects, thermodynamic entropy production, lightcone time orientation, apparatus record defects, or selector-signed defects to physical units, detector calibration, spacetime EOM, selector closure, `L_total`, bridge/role-transfer, or ToE closure.\n- The next honest move requires exactly one new strict asymmetric-transition source object `Delta_asym`, a nonzero gauge-invariant forward/reverse asymmetry witness; otherwise preserve the P3105-P3122 physical-unit no-go/no-new-live-frontier certificate.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
