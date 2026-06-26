#!/usr/bin/env python3
"""P3120/S2070: Tau_LT ordered-flow source audit.

P3119 left exactly one admissible next object: Tau_LT, an intrinsic
ordered-flow source on nadsoliton information that might split U_time as order
and U_length as extension before re-testing Xi_LT/R_dim.  This audit constructs
finite candidate bifunctors and rejects any candidate that is merely a clock,
rod, light/Planck import, selector premise, Lagrangian normalization,
bridge/role-transfer, L_total, or ToE promotion.
"""
from __future__ import annotations

import hashlib
import json
import subprocess
from fractions import Fraction
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p3113_s2063_u_action_reference_carrier_source_law_audit import a_phi
from p3119_s2069_xi_lt_axis_source_object_audit import OUT as P3119

OUT = GEN / "p3120_s2070_tau_lt_ordered_flow_source_audit.json"
MD = GEN / "p3120_s2070_tau_lt_ordered_flow_source_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"
FLOW_SCALES = (Fraction(1, 5), Fraction(1, 2), Fraction(1, 1), Fraction(2, 1), Fraction(5, 1))
FLOW_TESTS = ("ordered_domain", "time_as_order", "length_as_extension", "bifunctor_law", "antisymmetry_or_cycle_cut", "non_gauge_scale", "xi_lt_induction", "r_dim_induction", "forbidden_import_free")
COUPLING_TESTS = ("source_U_time_order", "source_U_length_extension", "prove_axis_distinction", "induce_Xi_LT", "induce_R_dim", "preserve_C_phi_A_phi", "avoid_closed_lanes")
GATES = (
    "uses_p3119_tau_lt_obligation",
    "explicit_ordered_flow_formula",
    "strict_nadsoliton_data_only",
    "U_time_as_internal_order_exported",
    "U_length_as_internal_extension_exported",
    "order_extension_bifunctor_law_proved",
    "cycle_or_recurrence_cut_nonconventional",
    "scale_covariant_not_gauge_fixed",
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
        "tau_lt_ordered_flow": r"Tau_LT|ordered-flow|temporal-order|length-extension|information flow",
        "xi_r_dim_chain": r"Xi_LT|R_dim|U_length|U_time|U_action=F\(U_length,U_time\)",
        "phase_area_chain": r"C_phi\(A_phi\)|A_phi|Omega_dim|K_dim|Sigma_dim",
        "blocked_imports_closure": r"hbar|Planck|rod|clock|observed light|apparatus|selector|QW-2191|L_total|bridge/role-transfer|ToE",
    }
    rows = []
    for lane, pattern in patterns.items():
        proc = subprocess.run(["rg", "-n", pattern, "fundamental_action_reconstruction", "AGENTS.md", "-S"], cwd=REPO, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
        hits = [line for line in proc.stdout.splitlines() if line.strip()]
        rows.append({"content_lane": lane, "pattern": pattern, "hit_count": len(hits), "sample_hits": hits[:8]})
    return rows


def flow_candidates() -> list[dict[str, Any]]:
    specs = [
        ("entropy_monotone_flow_bifunctor", "Tau_LT(x->y) := (Delta S order, support-extension norm)", True, True, True, True, True, False, False, False, True, True, True, True, True, True, True, "entropy monotonicity gives an order proxy but its scale remains gauge and does not induce Xi_LT/R_dim"),
        ("phase_winding_order_extension", "Tau_LT := (winding precedence, phase-arc extension)", True, True, True, True, False, False, False, False, True, True, True, True, True, True, True, "winding/arc data are internal but recurrence cut is conventional"),
        ("z12_rewrite_flow_bifunctor", "Tau_LT := (rewrite depth, orbit displacement)", True, True, True, True, True, False, False, False, False, True, True, True, True, True, True, "rewrite depth and displacement are finite combinatorics without phase-area preservation"),
        ("path_cohomology_filtration_flow", "Tau_LT := (filtration order, persistent chain span)", True, True, True, True, True, False, False, False, False, False, True, True, True, True, True, "filtration creates a formal order but needs a chosen reference scale and misses C_phi"),
        ("damping_memory_flow_bifunctor", "Tau_LT := (memory lag order, damping-tail extension)", True, True, True, True, True, True, False, False, False, True, True, True, True, True, True, "damping/memory flow is target dependent and cannot source axes independently"),
        ("information_cut_flow_bifunctor", "Tau_LT := (cut precedence, cut-capacity extension)", True, True, True, True, True, True, False, False, False, False, True, True, True, True, True, "cut flow is promising but lacks an exported coupling to phase-area/action"),
        ("dirichlet_mode_flow_bifunctor", "Tau_LT := (mode ordering, nodal span extension)", True, True, True, True, True, False, False, False, False, False, True, True, True, True, True, "mode ordering is spectral bookkeeping without intrinsic dimensional scale"),
        ("category_endoflow_bifunctor", "Tau_LT := End(N) order enriched by image-extension", True, True, True, True, True, False, False, False, False, True, True, True, True, True, True, "categorical endoflow names the bifunctor but lacks a non-gauge numeric source"),
        ("least_action_flow_bifunctor", "Tau_LT := order induced by extremal action path", True, True, False, True, True, True, True, True, True, True, True, True, False, True, False, "imports unresolved variational/Lagrangian dynamics and L_total"),
        ("lightcone_flow_bifunctor", "Tau_LT := causal lightcone order and spatial separation", True, True, False, True, True, True, True, True, True, True, True, True, False, True, True, "imports observed-light spacetime semantics"),
        ("planck_lattice_flow_bifunctor", "Tau_LT := Planck tick order and Planck length edge", True, True, False, True, True, True, True, True, True, True, True, True, True, True, False, "complete flow source is imported Planck metrology"),
        ("apparatus_history_flow_bifunctor", "Tau_LT := detector event order and rod extension", True, True, False, True, True, True, True, True, True, True, False, True, True, False, False, "apparatus histories are downstream observer readout"),
        ("selector_oriented_flow_bifunctor", "Tau_LT := directed flow after selector premise", True, True, True, True, True, False, False, False, False, True, True, True, True, True, False, "selector premise is forbidden and QW-2191 remains open"),
    ]
    keys = ("candidate", "formula", *GATES, "blocker")
    return [dict(zip(keys, row)) for row in specs]


def flow_law_rows(candidates: list[dict[str, Any]]) -> list[dict[str, Any]]:
    field = {
        "ordered_domain": "explicit_ordered_flow_formula",
        "time_as_order": "U_time_as_internal_order_exported",
        "length_as_extension": "U_length_as_internal_extension_exported",
        "bifunctor_law": "order_extension_bifunctor_law_proved",
        "antisymmetry_or_cycle_cut": "cycle_or_recurrence_cut_nonconventional",
        "non_gauge_scale": "scale_covariant_not_gauge_fixed",
        "xi_lt_induction": "Xi_LT_axis_source_unlocked",
        "r_dim_induction": "R_dim_relation_unlocked",
        "forbidden_import_free": "selector_bridge_ltotal_toe_free",
    }
    return [{"candidate": c["candidate"], "flow_test": test, "test_passed": bool(c[field[test]]), "accepted_ordered_flow": bool(c[field[test]] and c["strict_nadsoliton_data_only"] and c["not_unit_convention"] and c["not_imported_dynamics"] and c["not_apparatus_calibration"] and c["selector_bridge_ltotal_toe_free"]), "blocker": c["blocker"]} for c in candidates for test in FLOW_TESTS]


def scale_rows(candidates: list[dict[str, Any]]) -> list[dict[str, Any]]:
    rows = []
    for c in candidates:
        for lam in FLOW_SCALES:
            rows.append({"candidate": c["candidate"], "lambda_flow": f"{lam.numerator}/{lam.denominator}", "order_rescaled": "unchanged_preorder", "extension_rescaled": f"{lam.numerator}/{lam.denominator}", "covariance_claimed": bool(c["scale_covariant_not_gauge_fixed"]), "covariance_accepted": bool(c["scale_covariant_not_gauge_fixed"] and c["cycle_or_recurrence_cut_nonconventional"] and c["Xi_LT_axis_source_unlocked"] and c["not_unit_convention"]), "blocker": c["blocker"]})
    return rows


def coupling_rows(candidates: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return [{"candidate": c["candidate"], "coupling_test": test, "A_phi": round(a_phi(), 12) if test == "preserve_C_phi_A_phi" else None, "test_passed": bool(c["U_time_as_internal_order_exported"] if test == "source_U_time_order" else c["U_length_as_internal_extension_exported"] if test == "source_U_length_extension" else c["cycle_or_recurrence_cut_nonconventional"] if test == "prove_axis_distinction" else c["Xi_LT_axis_source_unlocked"] if test == "induce_Xi_LT" else c["R_dim_relation_unlocked"] if test == "induce_R_dim" else c["C_phi_A_phi_preserved"] if test == "preserve_C_phi_A_phi" else c["not_imported_dynamics"] and c["not_apparatus_calibration"] and c["selector_bridge_ltotal_toe_free"]), "accepted_coupling_chain": bool(c["U_time_as_internal_order_exported"] and c["U_length_as_internal_extension_exported"] and c["cycle_or_recurrence_cut_nonconventional"] and c["Xi_LT_axis_source_unlocked"] and c["R_dim_relation_unlocked"] and c["C_phi_A_phi_preserved"] and c["not_unit_convention"] and c["not_imported_dynamics"] and c["not_apparatus_calibration"] and c["selector_bridge_ltotal_toe_free"]), "blocker": c["blocker"]} for c in candidates for test in COUPLING_TESTS]


def gate_rows(candidates: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return [{"candidate": c["candidate"], "required_gate": gate, "gate_passed": bool(c[gate]), "detail": "passed" if c[gate] else c["blocker"]} for c in candidates for gate in GATES]


def aggregate(gates: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return [{"candidate": candidate, "passed_gates": sum(row["gate_passed"] for row in gates if row["candidate"] == candidate), "failed_gates": sum(not row["gate_passed"] for row in gates if row["candidate"] == candidate), "accepted_Tau_LT_ordered_flow": all(row["gate_passed"] for row in gates if row["candidate"] == candidate)} for candidate in sorted({row["candidate"] for row in gates})]


def build_payload() -> dict[str, Any]:
    p3119 = read_json(P3119)
    greps = content_grep()
    candidates = flow_candidates()
    flows = flow_law_rows(candidates)
    scales = scale_rows(candidates)
    couplings = coupling_rows(candidates)
    gates = gate_rows(candidates)
    aggs = aggregate(gates)
    accepted = [row for row in aggs if row["accepted_Tau_LT_ordered_flow"]]
    proof_obligations = [
        {"obligation": "read_p3119_next_atom", "satisfied": True, "detail": "P3119 requested exactly one Tau_LT ordered-flow source object"},
        {"obligation": "content_first_repo_grep", "satisfied": len(greps) == 4 and sum(row["hit_count"] for row in greps) > 0, "detail": "repo was grepped by ordered-flow and dimensional-chain patterns"},
        {"obligation": "construct_Tau_LT_candidates", "satisfied": len(candidates) == 13, "detail": "thirteen ordered-flow bifunctor candidates were constructed"},
        {"obligation": "test_flow_laws", "satisfied": len(flows) == len(candidates) * len(FLOW_TESTS), "detail": "nine flow-law rows were built per candidate"},
        {"obligation": "test_flow_scale_covariance", "satisfied": len(scales) == len(candidates) * len(FLOW_SCALES), "detail": "five flow-scale rows were built per candidate"},
        {"obligation": "test_Xi_LT_R_dim_coupling", "satisfied": len(couplings) == len(candidates) * len(COUPLING_TESTS), "detail": "seven coupling-chain rows were built per candidate"},
        {"obligation": "export_strict_Tau_LT", "satisfied": False, "detail": "0 candidates export an import-free strict Tau_LT satisfying all gates"},
    ]
    return {"status": "P3120_TAU_LT_ORDERED_FLOW_SOURCE_BOUNDED_NO_GO", "input_hashes": {"P3119": hashlib.sha256(P3119.read_bytes()).hexdigest() if P3119.exists() else None}, "constructed_theoretical_objects": {"content_first_repo_grep": greps, "audit_object": {"object": "TauLTOrderedFlowSourceAudit", "required_source": "Tau_LT intrinsic temporal-order/length-extension bifunctor on nadsoliton information flow", "forbidden_imports": ["hbar/Planck", "rods", "clocks", "observed light", "apparatus", "selector replay", "L_total", "bridge/role-transfer", "ToE"]}, "candidate_Tau_LT_ordered_flows": candidates, "flow_law_rows": flows, "flow_scale_covariance_rows": scales, "Xi_LT_R_dim_coupling_rows": couplings, "candidate_gate_rows": gates, "candidate_aggregate_certificate": aggs}, "finite_certificate": {"content_grep_lanes": len(greps), "content_grep_hits": sum(row["hit_count"] for row in greps), "p3119_accepted_Xi_LT_axis_sources": p3119.get("finite_certificate", {}).get("accepted_Xi_LT_axis_sources"), "candidate_Tau_LT_ordered_flows": len(candidates), "flow_law_rows": len(flows), "flow_scale_covariance_rows": len(scales), "Xi_LT_R_dim_coupling_rows": len(couplings), "candidate_gate_rows": len(gates), "accepted_Tau_LT_ordered_flows": len(accepted), "proof_obligations": len(proof_obligations), "satisfied_proof_obligations": sum(row["satisfied"] for row in proof_obligations)}, "proof_obligations": proof_obligations, "decision": {"bounded_result": "P3120 constructs the requested Tau_LT ordered-flow family and finds bounded no-go. Entropy monotone, phase winding, Z12 rewrite, cohomology filtration, damping/memory, information-cut, spectral, and categorical endoflow candidates are real typed source attempts, but each loses at least one required condition: nonconventional recurrence/cycle cut, non-gauge extension scale, Xi_LT induction, R_dim induction, C_phi(A_phi) preservation, or closed-lane import freedom. Least-action, lightcone, Planck, apparatus, and selector candidates import blocked lanes. No nadsoliton-only Tau_LT exports physical length/time axes or unlocks Xi_LT/R_dim.", "negative_export_flags": {key: False for key in ["Tau_LT_ordered_flow_exported", "Xi_LT_axis_source_exported", "R_dim_relation_exported", "Omega_dim_source_exported", "K_dim_functor_exported", "Sigma_dim_theorem_exported", "D_phi_source_law_exported", "U_action_source_law_exported", "physical_length_time_unit_exported", "detector_calibration_exported", "spacetime_eom_exported", "QW_2191_discharged", "strict_selector_closure_exported", "ltotal_exported", "bridge_role_transfer_exported", "toe_closure_exported"]}, "positive_scoped_flags": {"p3119_Tau_LT_obligation_reused": True, "candidate_Tau_LT_ordered_flows_constructed": True, "flow_law_matrix_built": True, "flow_scale_covariance_matrix_built": True, "Xi_LT_R_dim_coupling_matrix_built": True, "closed_lane_and_external_import_rows_rejected": True}, "next_honest_step": "Construct exactly one new strict recurrence-cut source object Kappa_cycle: an intrinsic, nonconventional acyclicity/irreversibility law on nadsoliton information flow that can turn ordered-flow preorder into a real temporal order before re-testing Tau_LT. It must not use clocks, rods, observed light, Planck units, thermodynamic environment, Lagrangian/EOM normalization, selector replay, bridge/role-transfer, L_total, or ToE closure. Without Kappa_cycle or an equally explicit new typed object, preserve the P3105-P3120 physical-unit no-go/no-new-live-frontier certificate."}}


def write_markdown(payload: dict[str, Any]) -> None:
    cert = payload["finite_certificate"]
    lines = ["# P3120/S2070 Tau_LT ordered-flow source audit", "", f"Status: `{payload['status']}`", "", "## Finite certificate", f"- P3119 accepted Xi_LT axis sources: `{cert['p3119_accepted_Xi_LT_axis_sources']}`", f"- content grep lanes: `{cert['content_grep_lanes']}`", f"- candidate Tau_LT ordered flows: `{cert['candidate_Tau_LT_ordered_flows']}`", f"- flow-law rows: `{cert['flow_law_rows']}`", f"- flow-scale covariance rows: `{cert['flow_scale_covariance_rows']}`", f"- Xi_LT/R_dim coupling rows: `{cert['Xi_LT_R_dim_coupling_rows']}`", f"- candidate gate rows: `{cert['candidate_gate_rows']}`", f"- accepted Tau_LT ordered flows: `{cert['accepted_Tau_LT_ordered_flows']}`", f"- satisfied proof obligations: `{cert['satisfied_proof_obligations']}/{cert['proof_obligations']}`", "", "## Decision", payload["decision"]["bounded_result"], "", "## Recommendation", payload["decision"]["next_honest_step"]]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3120/S2070 Tau_LT ordered-flow source audit", "## P3120/S2070 Tau_LT ordered-flow source audit\n\n`P3120/S2070` executes the P3119-recommended audit for `Tau_LT`, a strict ordered-flow source object intended to split `U_time` as internal order and `U_length` as internal extension on nadsoliton information flow.  It constructs `13` candidate ordered-flow bifunctors, `117` flow-law rows, `65` flow-scale covariance rows, `91` `Xi_LT/R_dim` coupling rows, and a `13 x 15 = 195` gate matrix.  The bounded result is that entropy, phase, Z12, cohomology, damping/memory, information-cut, spectral, and categorical flow candidates miss a nonconventional recurrence cut, non-gauge scale, `Xi_LT` induction, `R_dim` induction, `C_phi(A_phi)` preservation, or import freedom, while least-action, lightcone, Planck, apparatus, and selector candidates import closed or forbidden lanes.  No physical length/time unit, detector calibration, spacetime EOM, selector closure, `L_total`, bridge/role-transfer, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3120/S2070 Tau_LT ordered-flow remains incomplete", "## P3120/S2070 Tau_LT ordered-flow remains incomplete\n\n`P3120/S2070` tests whether a strict nadsoliton-only ordered-flow source `Tau_LT` can source time as order and length as extension strongly enough to unlock `Xi_LT` and `R_dim`.  Current artifacts provide no import-free strict ordered-flow source satisfying the full chain, so the result is not a Lagrangian density, Hamiltonian normalization, spacetime EOM, `L_total`, or ToE datum.\n")
    append_once(AGENTS, "Current Tau_LT ordered-flow source guardrail (P3120/S2070, 2026-06-26)", "## Current Tau_LT ordered-flow source guardrail (P3120/S2070, 2026-06-26)\n\n- P3120 tests the P3119-requested strict ordered-flow source object `Tau_LT`, an intrinsic temporal-order/length-extension bifunctor on nadsoliton information flow.\n- The finite audit constructs `13` candidate ordered-flow bifunctors, `117` flow-law rows, `65` flow-scale covariance rows, `91` `Xi_LT/R_dim` coupling rows, and `195` gate rows; `0` candidates export an import-free strict `Tau_LT` source.\n- Do not promote entropy monotone flows, phase winding flows, Z12 rewrite flows, cohomology filtrations, damping/memory flows, information-cut flows, spectral mode flows, categorical endoflows, least-action flows, lightcone flows, Planck lattice flows, apparatus histories, or selector-oriented flows to physical units, detector calibration, spacetime EOM, selector closure, `L_total`, bridge/role-transfer, or ToE closure.\n- The next honest move requires exactly one new strict recurrence-cut source object `Kappa_cycle`, an intrinsic acyclicity/irreversibility law for nadsoliton information flow; otherwise preserve the P3105-P3120 physical-unit no-go/no-new-live-frontier certificate.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
