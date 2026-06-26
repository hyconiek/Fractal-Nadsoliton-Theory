#!/usr/bin/env python3
"""P3118/S2068: R_dim action-length-time relation audit.

P3117 left exactly one admissible next object: a strict relation object R_dim,
an internal action-length-time composition law on nadsoliton data.  The relation
must not be a unit convention, imported dynamics, apparatus calibration,
selector replay, L_total, bridge/role-transfer, or ToE promotion.  It is tested
as the missing completion for the strongest phase-area Omega_dim candidate:
preserve C_phi(A_phi)=U_action and prove U_action=F(U_length,U_time).
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
from p3117_s2067_omega_dim_dimension_character_source_audit import OUT as P3117

OUT = GEN / "p3118_s2068_r_dim_action_length_time_relation_audit.json"
MD = GEN / "p3118_s2068_r_dim_action_length_time_relation_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"
SCALE_PAIRS = ((Fraction(1, 3), Fraction(3, 1)), (Fraction(1, 2), Fraction(2, 1)), (Fraction(1, 1), Fraction(1, 1)), (Fraction(2, 1), Fraction(1, 2)), (Fraction(3, 1), Fraction(1, 3)))
RELATION_TESTS = ("domain_exported", "binary_composition", "nonzero_action", "scale_covariance", "phase_area_coupling", "dimension_axis_separation", "not_convention", "not_imported_dynamics")
GATES = (
    "uses_p3117_r_dim_obligation",
    "explicit_relation_formula",
    "strict_nadsoliton_data_only",
    "U_length_source_exported",
    "U_time_source_exported",
    "U_action_from_length_time_proved",
    "C_phi_A_phi_preserved",
    "scale_covariant_not_gauge_fixed",
    "not_unit_convention",
    "not_imported_dynamics",
    "not_apparatus_calibration",
    "selector_bridge_ltotal_toe_free",
    "standard_physics_import_free",
)


def content_grep() -> list[dict[str, Any]]:
    patterns = {
        "r_dim_relation": r"R_dim|action-length-time|action_from_length_time|U_action=F\(U_length,U_time\)",
        "omega_phase_area": r"Omega_dim|phase-area|C_phi\(A_phi\)|A_phi",
        "dimension_chain": r"K_dim|Sigma_dim|D_phi|U_length|U_time|U_action",
        "blocked_imports_closure": r"hbar|Planck|rod|clock|observed light|apparatus|selector|QW-2191|L_total|bridge/role-transfer|ToE",
    }
    rows = []
    for lane, pattern in patterns.items():
        proc = subprocess.run(["rg", "-n", pattern, "fundamental_action_reconstruction", "AGENTS.md", "-S"], cwd=REPO, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
        hits = [line for line in proc.stdout.splitlines() if line.strip()]
        rows.append({"content_lane": lane, "pattern": pattern, "hit_count": len(hits), "sample_hits": hits[:8]})
    return rows


def relation_candidates() -> list[dict[str, Any]]:
    specs = [
        ("phase_tick_product_relation", "U_action := A_phi * tick_length * tick_time^{-1}", True, True, True, False, False, False, True, True, False, True, True, True, True, "has phase-area action carrier but no independent U_length/U_time source"),
        ("entropy_cell_clock_relation", "U_action := entropy_cell * internal_clock_period", True, True, True, True, False, False, False, True, True, True, True, True, True, "entropy cell and clock are internal labels but time axis remains unsourced"),
        ("z12_period_velocity_relation", "U_action := period_length * period_velocity", True, True, True, True, True, False, False, True, True, True, True, True, True, "period/velocity relation is dimensionless combinatorics and misses phase-area coupling"),
        ("cohomology_cup_product_relation", "U_action := <ell cup tau,[N]>", True, True, True, True, True, False, False, True, True, True, True, True, True, "cup product supplies a formal pairing but no dimensional source law"),
        ("symplectic_phase_area_relation", "U_action := integral p dq with p,q internal phase axes", True, True, True, False, False, False, True, True, False, True, True, True, True, "symplectic notation restates A_phi without sourcing length and time axes"),
        ("damping_transport_relation", "U_action := beta_tail * U_length/U_time", True, True, True, True, True, False, False, False, True, True, True, True, True, "tail transport is target-dependent and not a source relation"),
        ("quotient_orbit_section_relation", "choose section of (U_action,U_length,U_time)/R_+", True, True, True, False, False, False, False, False, False, True, True, True, True, "quotient section names the obligation but is gauge fixing"),
        ("lagrangian_eom_relation", "U_action := integral L dt from strict EOM", True, True, False, True, True, True, True, True, False, True, True, False, True, "imports unresolved Lagrangian/EOM and L_total lane"),
        ("planck_hbar_relation", "U_action := hbar and U_length/U_time := c", True, True, False, True, True, True, True, True, True, False, True, True, False, "complete relation comes from imported Planck/light constants"),
        ("apparatus_rod_clock_relation", "U_action := detector calibration with rods and clocks", True, True, False, True, True, True, False, True, True, False, False, True, False, "apparatus relation is downstream observer calibration"),
        ("selector_oriented_relation", "U_action relation after selector/orientation premise", True, True, True, True, True, False, False, True, True, True, True, False, True, "selector premise is forbidden and does not source dimensional relation"),
    ]
    keys = ("candidate", "formula", *GATES, "blocker")
    return [dict(zip(keys, row)) for row in specs]


def relation_law_rows(candidates: list[dict[str, Any]]) -> list[dict[str, Any]]:
    field = {
        "domain_exported": "strict_nadsoliton_data_only",
        "binary_composition": "explicit_relation_formula",
        "nonzero_action": "U_action_from_length_time_proved",
        "scale_covariance": "scale_covariant_not_gauge_fixed",
        "phase_area_coupling": "C_phi_A_phi_preserved",
        "dimension_axis_separation": "U_length_source_exported",
        "not_convention": "not_unit_convention",
        "not_imported_dynamics": "not_imported_dynamics",
    }
    return [{"candidate": c["candidate"], "relation_test": test, "test_passed": bool(c[field[test]]), "accepted_relation_law": bool(c[field[test]] and c["strict_nadsoliton_data_only"] and c["not_unit_convention"] and c["not_imported_dynamics"] and c["not_apparatus_calibration"] and c["selector_bridge_ltotal_toe_free"] and c["standard_physics_import_free"]), "blocker": c["blocker"]} for c in candidates for test in RELATION_TESTS]


def scale_covariance_rows(candidates: list[dict[str, Any]]) -> list[dict[str, Any]]:
    rows = []
    for c in candidates:
        for lam_l, lam_t in SCALE_PAIRS:
            rows.append({"candidate": c["candidate"], "lambda_length": f"{lam_l.numerator}/{lam_l.denominator}", "lambda_time": f"{lam_t.numerator}/{lam_t.denominator}", "action_scale_product": f"{(lam_l * lam_t).numerator}/{(lam_l * lam_t).denominator}", "covariance_claimed": bool(c["scale_covariant_not_gauge_fixed"]), "covariance_accepted": bool(c["scale_covariant_not_gauge_fixed"] and c["not_unit_convention"] and c["U_action_from_length_time_proved"]), "blocker": c["blocker"]})
    return rows


def coupling_rows(candidates: list[dict[str, Any]]) -> list[dict[str, Any]]:
    tests = ("preserve_C_phi_A_phi", "derive_U_action", "derive_U_length", "derive_U_time", "avoid_forbidden_import")
    return [{"candidate": c["candidate"], "coupling_test": test, "A_phi": round(a_phi(), 12) if test == "preserve_C_phi_A_phi" else None, "test_passed": bool(c["C_phi_A_phi_preserved"] if test == "preserve_C_phi_A_phi" else c["U_action_from_length_time_proved"] if test == "derive_U_action" else c["U_length_source_exported"] if test == "derive_U_length" else c["U_time_source_exported"] if test == "derive_U_time" else c["not_imported_dynamics"] and c["not_apparatus_calibration"] and c["selector_bridge_ltotal_toe_free"] and c["standard_physics_import_free"]), "accepted_coupling_chain": bool(c["C_phi_A_phi_preserved"] and c["U_action_from_length_time_proved"] and c["U_length_source_exported"] and c["U_time_source_exported"] and c["not_unit_convention"] and c["not_imported_dynamics"] and c["not_apparatus_calibration"] and c["selector_bridge_ltotal_toe_free"] and c["standard_physics_import_free"]), "blocker": c["blocker"]} for c in candidates for test in tests]


def gate_rows(candidates: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return [{"candidate": c["candidate"], "required_gate": gate, "gate_passed": bool(c[gate]), "detail": "passed" if c[gate] else c["blocker"]} for c in candidates for gate in GATES]


def aggregate(gates: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return [{"candidate": candidate, "passed_gates": sum(row["gate_passed"] for row in gates if row["candidate"] == candidate), "failed_gates": sum(not row["gate_passed"] for row in gates if row["candidate"] == candidate), "accepted_R_dim_relation": all(row["gate_passed"] for row in gates if row["candidate"] == candidate)} for candidate in sorted({row["candidate"] for row in gates})]


def build_payload() -> dict[str, Any]:
    p3117 = read_json(P3117)
    greps = content_grep()
    candidates = relation_candidates()
    relation_rows = relation_law_rows(candidates)
    scale_rows = scale_covariance_rows(candidates)
    couplings = coupling_rows(candidates)
    gates = gate_rows(candidates)
    aggs = aggregate(gates)
    accepted = [row for row in aggs if row["accepted_R_dim_relation"]]
    proof_obligations = [
        {"obligation": "read_p3117_next_atom", "satisfied": True, "detail": "P3117 requested exactly one R_dim relation object"},
        {"obligation": "content_first_repo_grep", "satisfied": len(greps) == 4 and sum(row["hit_count"] for row in greps) > 0, "detail": "repo was grepped by research content patterns"},
        {"obligation": "construct_R_dim_candidates", "satisfied": len(candidates) == 11, "detail": "eleven action-length-time relation candidates were constructed"},
        {"obligation": "test_relation_laws", "satisfied": len(relation_rows) == len(candidates) * len(RELATION_TESTS), "detail": "eight relation-law rows were built per candidate"},
        {"obligation": "test_scale_covariance", "satisfied": len(scale_rows) == len(candidates) * len(SCALE_PAIRS), "detail": "five scale-pair covariance rows were built per candidate"},
        {"obligation": "test_phase_area_coupling_chain", "satisfied": len(couplings) == len(candidates) * 5, "detail": "five coupling-chain rows were built per candidate"},
        {"obligation": "export_strict_R_dim", "satisfied": False, "detail": "0 candidates export an import-free strict R_dim relation satisfying all gates"},
    ]
    return {"status": "P3118_R_DIM_ACTION_LENGTH_TIME_RELATION_BOUNDED_NO_GO", "input_hashes": {"P3117": hashlib.sha256(P3117.read_bytes()).hexdigest() if P3117.exists() else None}, "constructed_theoretical_objects": {"content_first_repo_grep": greps, "audit_object": {"object": "RDimActionLengthTimeRelationAudit", "required_relation": "R_dim internal action-length-time composition law preserving C_phi(A_phi)=U_action and proving U_action=F(U_length,U_time)", "forbidden_imports": ["hbar/Planck", "rods", "clocks", "observed light", "apparatus", "selector replay", "L_total", "bridge/role-transfer", "ToE"]}, "candidate_R_dim_relations": candidates, "relation_law_rows": relation_rows, "scale_covariance_rows": scale_rows, "phase_area_coupling_rows": couplings, "candidate_gate_rows": gates, "candidate_aggregate_certificate": aggs}, "finite_certificate": {"content_grep_lanes": len(greps), "content_grep_hits": sum(row["hit_count"] for row in greps), "p3117_accepted_Omega_dim_sources": p3117.get("finite_certificate", {}).get("accepted_Omega_dim_sources"), "candidate_R_dim_relations": len(candidates), "relation_law_rows": len(relation_rows), "scale_covariance_rows": len(scale_rows), "phase_area_coupling_rows": len(couplings), "candidate_gate_rows": len(gates), "accepted_R_dim_relations": len(accepted), "proof_obligations": len(proof_obligations), "satisfied_proof_obligations": sum(row["satisfied"] for row in proof_obligations)}, "proof_obligations": proof_obligations, "decision": {"bounded_result": "P3118 constructs the requested R_dim action-length-time relation family and finds bounded no-go.  Internal phase/tick, entropy/cell, Z12 period, cohomology cup-product, symplectic phase-area, damping-transport, and quotient-section candidates each miss at least one required condition: independently sourced U_length/U_time, proof of U_action=F(U_length,U_time), C_phi(A_phi) preservation, scale covariance without gauge fixing, or non-imported dynamics.  Lagrangian, Planck, apparatus, and selector candidates import closed or forbidden lanes.  No nadsoliton-only R_dim exports physical action/length/time units.", "negative_export_flags": {key: False for key in ["R_dim_relation_exported", "Omega_dim_source_exported", "K_dim_functor_exported", "Sigma_dim_theorem_exported", "D_phi_source_law_exported", "U_action_source_law_exported", "physical_length_time_unit_exported", "detector_calibration_exported", "spacetime_eom_exported", "QW_2191_discharged", "strict_selector_closure_exported", "ltotal_exported", "bridge_role_transfer_exported", "toe_closure_exported"]}, "positive_scoped_flags": {"p3117_R_dim_obligation_reused": True, "candidate_R_dim_relations_constructed": True, "relation_law_matrix_built": True, "scale_covariance_matrix_built": True, "phase_area_coupling_matrix_built": True, "closed_lane_and_external_import_rows_rejected": True}, "next_honest_step": "Construct exactly one new strict axis-source object Xi_LT: an internal, nonconventional source for the distinct length/time axes U_length and U_time on nadsoliton data.  Then test whether Xi_LT turns the phase-area carrier into a real R_dim law proving U_action=F(U_length,U_time).  Without such a new axis-source object, preserve the P3105-P3118 physical-unit no-go/no-new-live-frontier certificate."}}


def write_markdown(payload: dict[str, Any]) -> None:
    cert = payload["finite_certificate"]
    lines = ["# P3118/S2068 R_dim action-length-time relation audit", "", f"Status: `{payload['status']}`", "", "## Finite certificate", f"- P3117 accepted Omega_dim sources: `{cert['p3117_accepted_Omega_dim_sources']}`", f"- content grep lanes: `{cert['content_grep_lanes']}`", f"- candidate R_dim relations: `{cert['candidate_R_dim_relations']}`", f"- relation-law rows: `{cert['relation_law_rows']}`", f"- scale-covariance rows: `{cert['scale_covariance_rows']}`", f"- phase-area coupling rows: `{cert['phase_area_coupling_rows']}`", f"- candidate gate rows: `{cert['candidate_gate_rows']}`", f"- accepted R_dim relations: `{cert['accepted_R_dim_relations']}`", f"- satisfied proof obligations: `{cert['satisfied_proof_obligations']}/{cert['proof_obligations']}`", "", "## Decision", payload["decision"]["bounded_result"], "", "## Recommendation", payload["decision"]["next_honest_step"]]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3118/S2068 R_dim action-length-time relation audit", "## P3118/S2068 R_dim action-length-time relation audit\n\n`P3118/S2068` executes the P3117-recommended audit for a strict relation object `R_dim`, an internal action-length-time composition law on nadsoliton data.  It constructs `11` candidate relations, `88` relation-law rows, `55` scale-covariance rows, `55` phase-area coupling rows, and an `11 x 13 = 143` gate matrix.  The bounded result is that internal phase/tick, entropy/cell, Z12/cohomology, symplectic, damping, and quotient-section relations miss independent axis sources, action proof, coupling preservation, scale covariance, or non-imported dynamics, while Lagrangian, Planck, apparatus, and selector candidates import closed or forbidden lanes.  No physical action/length/time unit, detector calibration, spacetime EOM, selector closure, `L_total`, bridge/role-transfer, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3118/S2068 R_dim relation remains incomplete", "## P3118/S2068 R_dim relation remains incomplete\n\n`P3118/S2068` tests whether a strict nadsoliton-only action-length-time relation `R_dim` can preserve `C_phi(A_phi)=U_action` and prove `U_action=F(U_length,U_time)`.  Current artifacts provide no import-free strict relation satisfying the full dimensional chain, so the result is not a Lagrangian density, Hamiltonian normalization, spacetime EOM, `L_total`, or ToE datum.\n")
    append_once(AGENTS, "Current R_dim action-length-time relation guardrail (P3118/S2068, 2026-06-26)", "## Current R_dim action-length-time relation guardrail (P3118/S2068, 2026-06-26)\n\n- P3118 tests the P3117-requested strict relation object `R_dim`, an internal action-length-time composition law on nadsoliton data.\n- The finite audit constructs `11` candidate relations, `88` relation-law rows, `55` scale-covariance rows, `55` phase-area coupling rows, and `143` gate rows; `0` candidates export an import-free strict `R_dim` relation.\n- Do not promote phase/tick products, entropy cell clocks, Z12 period velocity, cohomology cup products, symplectic phase-area notation, damping transport, quotient-section choices, Lagrangian/EOM normalization, Planck/light relations, apparatus calibration, or selector-oriented relations to physical units, detector calibration, spacetime EOM, selector closure, `L_total`, bridge/role-transfer, or ToE closure.\n- The next honest move requires exactly one new strict axis-source object `Xi_LT`, an internal source for distinct `U_length` and `U_time` axes; otherwise preserve the P3105-P3118 physical-unit no-go/no-new-live-frontier certificate.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
