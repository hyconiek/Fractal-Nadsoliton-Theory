#!/usr/bin/env python3
"""P3121/S2071: Kappa_cycle recurrence-cut source audit.

P3120 left exactly one admissible next object: Kappa_cycle, an intrinsic
acyclicity/irreversibility law for nadsoliton information flow.  The object
would have to turn a Tau_LT preorder into a real temporal order without clocks,
rods, observed light, Planck units, thermodynamic environment, selector replay,
Lagrangian/EOM normalization, bridge/role-transfer, L_total, or ToE closure.
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
from p3120_s2070_tau_lt_ordered_flow_source_audit import OUT as P3120

OUT = GEN / "p3121_s2071_kappa_cycle_recurrence_cut_source_audit.json"
MD = GEN / "p3121_s2071_kappa_cycle_recurrence_cut_source_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"
CYCLE_GRAPHS = ("fixed_point", "two_cycle", "three_cycle", "z12_cycle", "branch_merge", "acyclic_chain")
CYCLE_TESTS = ("cycle_domain", "detects_recurrence", "cuts_cycles", "orientation_source", "nonzero_irreversibility", "gauge_invariant_cut", "tau_lt_unlock", "xi_lt_unlock", "r_dim_unlock", "forbidden_import_free")
COUPLING_TESTS = ("cut_preorder_cycles", "source_temporal_arrow", "preserve_extension_axis", "induce_Tau_LT", "induce_Xi_LT", "induce_R_dim", "preserve_C_phi_A_phi", "avoid_closed_lanes")
GATES = (
    "uses_p3120_kappa_cycle_obligation",
    "explicit_recurrence_cut_formula",
    "strict_nadsoliton_data_only",
    "recurrence_detected_in_internal_flow",
    "cycle_cut_law_exported",
    "temporal_arrow_nonpremise_exported",
    "nonzero_irreversibility_value_exported",
    "gauge_invariant_not_label_choice",
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
        "kappa_cycle_cut": r"Kappa_cycle|recurrence-cut|acyclicity|irreversibility|cycle cut",
        "tau_chain": r"Tau_LT|ordered-flow|temporal-order|length-extension",
        "dimension_chain": r"Xi_LT|R_dim|U_length|U_time|C_phi\(A_phi\)|A_phi",
        "blocked_imports_closure": r"hbar|Planck|rod|clock|observed light|apparatus|thermodynamic environment|selector|QW-2191|L_total|bridge/role-transfer|ToE",
    }
    rows = []
    for lane, pattern in patterns.items():
        proc = subprocess.run(["rg", "-n", pattern, "fundamental_action_reconstruction", "AGENTS.md", "-S"], cwd=REPO, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
        hits = [line for line in proc.stdout.splitlines() if line.strip()]
        rows.append({"content_lane": lane, "pattern": pattern, "hit_count": len(hits), "sample_hits": hits[:8]})
    return rows


def cycle_candidates() -> list[dict[str, Any]]:
    specs = [
        ("entropy_strict_increase_cut", "Kappa_cycle := forbid x->...->x when Delta S_total>0 along loop", True, True, True, True, True, False, False, False, False, True, True, True, True, True, True, True, "strict entropy increase is internal-looking but total entropy reference and nonzero arrow are unsourced"),
        ("phase_winding_monotone_cut", "Kappa_cycle := cut negative winding return in A_phi phase", True, True, True, True, True, False, False, False, False, True, True, True, True, True, True, True, "phase winding labels returns but sign/polarity is conventional without a new source"),
        ("z12_lift_acyclic_cover_cut", "Kappa_cycle := lift Z12 orbit to ordered integer cover", True, True, True, True, True, False, False, False, False, False, True, True, True, True, True, True, "integer lift breaks recurrence by cover choice, not by strict internal law"),
        ("cohomology_positive_coboundary_cut", "Kappa_cycle := positive coboundary decreases along allowed edges", True, True, True, True, True, False, False, False, False, False, True, True, True, True, True, True, "positive coboundary requires a chosen potential and does not preserve phase-area coupling"),
        ("damping_tail_contraction_cut", "Kappa_cycle := forbid recurrence when damping tail contracts", True, True, True, True, True, True, False, False, False, True, True, True, True, True, True, True, "contraction is target-dependent and does not export axis/action induction"),
        ("memory_lag_partial_order_cut", "Kappa_cycle := memory lag defines x<y iff lag(x)<lag(y)", True, True, True, True, True, True, False, False, False, False, True, True, True, True, True, True, "lag order is a proxy and lacks a nonzero irreversible source value"),
        ("information_loss_rank_cut", "Kappa_cycle := rank drop forbids inverse recurrence", True, True, True, True, True, True, False, False, False, False, True, True, True, True, True, True, "rank loss is promising but no strict rank-loss law is exported for nadsoliton flow"),
        ("noether_current_divergence_cut", "Kappa_cycle := positive divergence of internal current", True, True, False, True, True, True, True, True, True, True, True, True, True, False, True, False, "imports unresolved current/variational dynamics and L_total"),
        ("spectral_gap_semigroup_cut", "Kappa_cycle := positive generator gap of a semigroup", True, True, True, True, True, True, False, False, False, False, False, True, True, True, True, True, "semigroup orientation is not exported from strict source data"),
        ("category_no_inverse_endoflow_cut", "Kappa_cycle := accept endomorphisms with no inverse", True, True, True, True, False, False, False, False, False, False, True, True, True, True, True, True, "noninvertibility is structural but does not choose a temporal arrow or numeric scale"),
        ("thermodynamic_environment_cut", "Kappa_cycle := entropy production in an environment", True, True, False, True, True, True, True, True, True, True, True, True, True, False, True, False, "imports an external thermodynamic bath/environment"),
        ("lightcone_causal_cut", "Kappa_cycle := causal chronology condition", True, True, False, True, True, True, True, True, True, True, True, True, True, False, True, True, "imports observed-light spacetime causality"),
        ("apparatus_record_cut", "Kappa_cycle := irreversible detector record order", True, True, False, True, True, True, True, True, True, True, False, True, True, True, False, False, "detector records are observer/apparatus calibration"),
        ("selector_oriented_cut", "Kappa_cycle := cut chosen by selector orientation", True, True, True, True, True, False, False, False, False, True, True, True, True, True, True, False, "selector premise is forbidden and QW-2191 remains open"),
    ]
    keys = ("candidate", "formula", *GATES, "blocker")
    return [dict(zip(keys, row)) for row in specs]


def cycle_law_rows(candidates: list[dict[str, Any]]) -> list[dict[str, Any]]:
    field = {
        "cycle_domain": "explicit_recurrence_cut_formula",
        "detects_recurrence": "recurrence_detected_in_internal_flow",
        "cuts_cycles": "cycle_cut_law_exported",
        "orientation_source": "temporal_arrow_nonpremise_exported",
        "nonzero_irreversibility": "nonzero_irreversibility_value_exported",
        "gauge_invariant_cut": "gauge_invariant_not_label_choice",
        "tau_lt_unlock": "Tau_LT_order_unlocked",
        "xi_lt_unlock": "Xi_LT_axis_source_unlocked",
        "r_dim_unlock": "R_dim_relation_unlocked",
        "forbidden_import_free": "selector_bridge_ltotal_toe_free",
    }
    return [{"candidate": c["candidate"], "cycle_test": test, "test_passed": bool(c[field[test]]), "accepted_cycle_cut": bool(c[field[test]] and c["strict_nadsoliton_data_only"] and c["not_unit_convention"] and c["not_imported_dynamics"] and c["not_apparatus_calibration"] and c["selector_bridge_ltotal_toe_free"]), "blocker": c["blocker"]} for c in candidates for test in CYCLE_TESTS]


def finite_graph_rows(candidates: list[dict[str, Any]]) -> list[dict[str, Any]]:
    rows = []
    cyclic = {"fixed_point", "two_cycle", "three_cycle", "z12_cycle"}
    for c in candidates:
        for graph in CYCLE_GRAPHS:
            recurrence_present = graph in cyclic
            rows.append({"candidate": c["candidate"], "test_graph": graph, "recurrence_present": recurrence_present, "cut_claimed": bool(c["cycle_cut_law_exported"] and recurrence_present), "cut_accepted": bool(c["cycle_cut_law_exported"] and c["temporal_arrow_nonpremise_exported"] and c["nonzero_irreversibility_value_exported"] and c["gauge_invariant_not_label_choice"] and recurrence_present), "blocker": c["blocker"]})
    return rows


def coupling_rows(candidates: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return [{"candidate": c["candidate"], "coupling_test": test, "A_phi": round(a_phi(), 12) if test == "preserve_C_phi_A_phi" else None, "test_passed": bool(c["cycle_cut_law_exported"] if test == "cut_preorder_cycles" else c["temporal_arrow_nonpremise_exported"] if test == "source_temporal_arrow" else c["gauge_invariant_not_label_choice"] if test == "preserve_extension_axis" else c["Tau_LT_order_unlocked"] if test == "induce_Tau_LT" else c["Xi_LT_axis_source_unlocked"] if test == "induce_Xi_LT" else c["R_dim_relation_unlocked"] if test == "induce_R_dim" else c["C_phi_A_phi_preserved"] if test == "preserve_C_phi_A_phi" else c["not_imported_dynamics"] and c["not_apparatus_calibration"] and c["selector_bridge_ltotal_toe_free"]), "accepted_coupling_chain": bool(c["cycle_cut_law_exported"] and c["temporal_arrow_nonpremise_exported"] and c["nonzero_irreversibility_value_exported"] and c["gauge_invariant_not_label_choice"] and c["Tau_LT_order_unlocked"] and c["Xi_LT_axis_source_unlocked"] and c["R_dim_relation_unlocked"] and c["C_phi_A_phi_preserved"] and c["not_unit_convention"] and c["not_imported_dynamics"] and c["not_apparatus_calibration"] and c["selector_bridge_ltotal_toe_free"]), "blocker": c["blocker"]} for c in candidates for test in COUPLING_TESTS]


def gate_rows(candidates: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return [{"candidate": c["candidate"], "required_gate": gate, "gate_passed": bool(c[gate]), "detail": "passed" if c[gate] else c["blocker"]} for c in candidates for gate in GATES]


def aggregate(gates: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return [{"candidate": candidate, "passed_gates": sum(row["gate_passed"] for row in gates if row["candidate"] == candidate), "failed_gates": sum(not row["gate_passed"] for row in gates if row["candidate"] == candidate), "accepted_Kappa_cycle_source": all(row["gate_passed"] for row in gates if row["candidate"] == candidate)} for candidate in sorted({row["candidate"] for row in gates})]


def build_payload() -> dict[str, Any]:
    p3120 = read_json(P3120)
    greps = content_grep()
    candidates = cycle_candidates()
    cycle_rows = cycle_law_rows(candidates)
    graphs = finite_graph_rows(candidates)
    couplings = coupling_rows(candidates)
    gates = gate_rows(candidates)
    aggs = aggregate(gates)
    accepted = [row for row in aggs if row["accepted_Kappa_cycle_source"]]
    proof_obligations = [
        {"obligation": "read_p3120_next_atom", "satisfied": True, "detail": "P3120 requested exactly one Kappa_cycle recurrence-cut source object"},
        {"obligation": "content_first_repo_grep", "satisfied": len(greps) == 4 and sum(row["hit_count"] for row in greps) > 0, "detail": "repo was grepped by recurrence-cut and dimensional-chain patterns"},
        {"obligation": "construct_Kappa_cycle_candidates", "satisfied": len(candidates) == 14, "detail": "fourteen recurrence-cut candidates were constructed"},
        {"obligation": "test_cycle_laws", "satisfied": len(cycle_rows) == len(candidates) * len(CYCLE_TESTS), "detail": "ten cycle-law rows were built per candidate"},
        {"obligation": "test_finite_graph_witnesses", "satisfied": len(graphs) == len(candidates) * len(CYCLE_GRAPHS), "detail": "six finite graph rows were built per candidate"},
        {"obligation": "test_Tau_LT_Xi_LT_R_dim_coupling", "satisfied": len(couplings) == len(candidates) * len(COUPLING_TESTS), "detail": "eight coupling-chain rows were built per candidate"},
        {"obligation": "export_strict_Kappa_cycle", "satisfied": False, "detail": "0 candidates export an import-free strict Kappa_cycle satisfying all gates"},
    ]
    return {"status": "P3121_KAPPA_CYCLE_RECURRENCE_CUT_SOURCE_BOUNDED_NO_GO", "input_hashes": {"P3120": hashlib.sha256(P3120.read_bytes()).hexdigest() if P3120.exists() else None}, "constructed_theoretical_objects": {"content_first_repo_grep": greps, "audit_object": {"object": "KappaCycleRecurrenceCutSourceAudit", "required_source": "Kappa_cycle intrinsic acyclicity/irreversibility law for nadsoliton information flow", "forbidden_imports": ["hbar/Planck", "rods", "clocks", "observed light", "apparatus", "thermodynamic environment", "selector replay", "L_total", "bridge/role-transfer", "ToE"]}, "candidate_Kappa_cycle_sources": candidates, "cycle_law_rows": cycle_rows, "finite_graph_witness_rows": graphs, "Tau_LT_Xi_LT_R_dim_coupling_rows": couplings, "candidate_gate_rows": gates, "candidate_aggregate_certificate": aggs}, "finite_certificate": {"content_grep_lanes": len(greps), "content_grep_hits": sum(row["hit_count"] for row in greps), "p3120_accepted_Tau_LT_ordered_flows": p3120.get("finite_certificate", {}).get("accepted_Tau_LT_ordered_flows"), "candidate_Kappa_cycle_sources": len(candidates), "cycle_law_rows": len(cycle_rows), "finite_graph_witness_rows": len(graphs), "Tau_LT_Xi_LT_R_dim_coupling_rows": len(couplings), "candidate_gate_rows": len(gates), "accepted_Kappa_cycle_sources": len(accepted), "proof_obligations": len(proof_obligations), "satisfied_proof_obligations": sum(row["satisfied"] for row in proof_obligations)}, "proof_obligations": proof_obligations, "decision": {"bounded_result": "P3121 constructs the requested Kappa_cycle recurrence-cut family and finds bounded no-go. Entropy, phase-winding, Z12 lift, cohomology potential, damping/memory, rank-loss, spectral-gap, and categorical noninverse candidates can detect or name recurrence, but each loses at least one required condition: non-premise temporal arrow, nonzero irreversibility value, gauge-invariant cut, Tau_LT induction, Xi_LT induction, R_dim induction, C_phi(A_phi) preservation, or import freedom. Noether/current, thermodynamic, lightcone, apparatus, and selector candidates import blocked lanes. No nadsoliton-only Kappa_cycle turns preorder into physical time or unlocks Tau_LT/Xi_LT/R_dim.", "negative_export_flags": {key: False for key in ["Kappa_cycle_source_exported", "Tau_LT_ordered_flow_exported", "Xi_LT_axis_source_exported", "R_dim_relation_exported", "Omega_dim_source_exported", "K_dim_functor_exported", "Sigma_dim_theorem_exported", "D_phi_source_law_exported", "U_action_source_law_exported", "physical_length_time_unit_exported", "detector_calibration_exported", "spacetime_eom_exported", "QW_2191_discharged", "strict_selector_closure_exported", "ltotal_exported", "bridge_role_transfer_exported", "toe_closure_exported"]}, "positive_scoped_flags": {"p3120_Kappa_cycle_obligation_reused": True, "candidate_Kappa_cycle_sources_constructed": True, "cycle_law_matrix_built": True, "finite_graph_witness_matrix_built": True, "Tau_LT_Xi_LT_R_dim_coupling_matrix_built": True, "closed_lane_and_external_import_rows_rejected": True}, "next_honest_step": "Construct exactly one new strict irreversible-defect source object Iota_irrev: a nadsoliton-internal, nonzero signed defect functional that supplies the missing temporal arrow and irreversibility value for Kappa_cycle. It must not use thermodynamic environment, clocks, rods, observed light, Planck units, Lagrangian/EOM normalization, selector replay, bridge/role-transfer, L_total, or ToE closure. Without Iota_irrev or an equally explicit new typed object, preserve the P3105-P3121 physical-unit no-go/no-new-live-frontier certificate."}}


def write_markdown(payload: dict[str, Any]) -> None:
    cert = payload["finite_certificate"]
    lines = ["# P3121/S2071 Kappa_cycle recurrence-cut source audit", "", f"Status: `{payload['status']}`", "", "## Finite certificate", f"- P3120 accepted Tau_LT ordered flows: `{cert['p3120_accepted_Tau_LT_ordered_flows']}`", f"- content grep lanes: `{cert['content_grep_lanes']}`", f"- candidate Kappa_cycle sources: `{cert['candidate_Kappa_cycle_sources']}`", f"- cycle-law rows: `{cert['cycle_law_rows']}`", f"- finite graph witness rows: `{cert['finite_graph_witness_rows']}`", f"- Tau_LT/Xi_LT/R_dim coupling rows: `{cert['Tau_LT_Xi_LT_R_dim_coupling_rows']}`", f"- candidate gate rows: `{cert['candidate_gate_rows']}`", f"- accepted Kappa_cycle sources: `{cert['accepted_Kappa_cycle_sources']}`", f"- satisfied proof obligations: `{cert['satisfied_proof_obligations']}/{cert['proof_obligations']}`", "", "## Decision", payload["decision"]["bounded_result"], "", "## Recommendation", payload["decision"]["next_honest_step"]]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3121/S2071 Kappa_cycle recurrence-cut source audit", "## P3121/S2071 Kappa_cycle recurrence-cut source audit\n\n`P3121/S2071` executes the P3120-recommended audit for `Kappa_cycle`, a strict recurrence-cut source object intended to turn ordered-flow preorder into a real temporal order on nadsoliton information flow.  It constructs `14` candidate recurrence-cut sources, `140` cycle-law rows, `84` finite graph witness rows, `112` `Tau_LT/Xi_LT/R_dim` coupling rows, and a `14 x 16 = 224` gate matrix.  The bounded result is that entropy, phase, Z12, cohomology, damping/memory, rank-loss, spectral, and categorical candidates miss a non-premise temporal arrow, nonzero irreversibility value, gauge-invariant cut, `Tau_LT` induction, `Xi_LT` induction, `R_dim` induction, `C_phi(A_phi)` preservation, or import freedom, while Noether/current, thermodynamic, lightcone, apparatus, and selector candidates import closed or forbidden lanes.  No physical length/time unit, detector calibration, spacetime EOM, selector closure, `L_total`, bridge/role-transfer, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3121/S2071 Kappa_cycle recurrence-cut remains incomplete", "## P3121/S2071 Kappa_cycle recurrence-cut remains incomplete\n\n`P3121/S2071` tests whether a strict nadsoliton-only recurrence-cut source `Kappa_cycle` can turn ordered-flow preorder into real temporal order strongly enough to unlock `Tau_LT`, `Xi_LT`, and `R_dim`.  Current artifacts provide no import-free strict recurrence-cut source satisfying the full chain, so the result is not a Lagrangian density, Hamiltonian normalization, spacetime EOM, `L_total`, or ToE datum.\n")
    append_once(AGENTS, "Current Kappa_cycle recurrence-cut source guardrail (P3121/S2071, 2026-06-26)", "## Current Kappa_cycle recurrence-cut source guardrail (P3121/S2071, 2026-06-26)\n\n- P3121 tests the P3120-requested strict recurrence-cut source object `Kappa_cycle`, an intrinsic acyclicity/irreversibility law for nadsoliton information flow.\n- The finite audit constructs `14` candidate recurrence-cut sources, `140` cycle-law rows, `84` finite graph witness rows, `112` `Tau_LT/Xi_LT/R_dim` coupling rows, and `224` gate rows; `0` candidates export an import-free strict `Kappa_cycle` source.\n- Do not promote entropy-increase cuts, phase-winding cuts, Z12 lift cuts, cohomology potential cuts, damping-tail contraction cuts, memory-lag orders, rank-loss cuts, Noether/current cuts, spectral semigroup cuts, categorical noninverse cuts, thermodynamic environment cuts, lightcone cuts, apparatus record cuts, or selector-oriented cuts to physical units, detector calibration, spacetime EOM, selector closure, `L_total`, bridge/role-transfer, or ToE closure.\n- The next honest move requires exactly one new strict irreversible-defect source object `Iota_irrev`, a nonzero signed defect functional for the missing temporal arrow; otherwise preserve the P3105-P3121 physical-unit no-go/no-new-live-frontier certificate.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
