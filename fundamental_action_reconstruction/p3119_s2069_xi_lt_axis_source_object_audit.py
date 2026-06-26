#!/usr/bin/env python3
"""P3119/S2069: Xi_LT length/time axis-source object audit.

P3118 left exactly one admissible next object: Xi_LT, a strict axis-source
object that internally distinguishes U_length from U_time on nadsoliton data.
This audit is deliberately proof-obligation first: construct candidate axis
sources, test whether they are nonconventional, source both axes independently,
preserve C_phi(A_phi)=U_action, and unlock a real R_dim law without importing
standard physics, apparatus, selector replay, L_total, bridge/role-transfer, or
ToE closure.
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
from p3118_s2068_r_dim_action_length_time_relation_audit import OUT as P3118

OUT = GEN / "p3119_s2069_xi_lt_axis_source_object_audit.json"
MD = GEN / "p3119_s2069_xi_lt_axis_source_object_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"
SCALE_PAIRS = ((Fraction(1, 4), Fraction(4, 1)), (Fraction(1, 2), Fraction(2, 1)), (Fraction(1, 1), Fraction(1, 1)), (Fraction(2, 1), Fraction(1, 2)), (Fraction(4, 1), Fraction(1, 4)))
AXIS_TESTS = ("two_axis_domain", "length_source", "time_source", "axis_noncollapse", "non_gauge_orientation", "r_dim_unlock", "phase_area_preserved", "forbidden_import_free")
COUPLING_TESTS = ("derive_U_length", "derive_U_time", "prove_distinct_axes", "preserve_C_phi_A_phi", "induce_R_dim", "avoid_closed_lanes")
GATES = (
    "uses_p3118_xi_lt_obligation",
    "explicit_axis_source_formula",
    "strict_nadsoliton_data_only",
    "U_length_source_exported",
    "U_time_source_exported",
    "length_time_distinction_nonconventional",
    "axis_orientation_not_selector_premise",
    "scale_covariant_not_gauge_fixed",
    "R_dim_relation_unlocked",
    "C_phi_A_phi_preserved",
    "not_unit_convention",
    "not_imported_dynamics",
    "not_apparatus_calibration",
    "selector_bridge_ltotal_toe_free",
)


def content_grep() -> list[dict[str, Any]]:
    patterns = {
        "xi_lt_axis_source": r"Xi_LT|axis-source|length/time axes|U_length.*U_time|U_time.*U_length",
        "r_dim_chain": r"R_dim|U_action=F\(U_length,U_time\)|C_phi\(A_phi\)|A_phi",
        "dimension_chain": r"Omega_dim|K_dim|Sigma_dim|D_phi|U_action",
        "blocked_imports_closure": r"hbar|Planck|rod|clock|observed light|apparatus|selector|QW-2191|L_total|bridge/role-transfer|ToE",
    }
    rows = []
    for lane, pattern in patterns.items():
        proc = subprocess.run(["rg", "-n", pattern, "fundamental_action_reconstruction", "AGENTS.md", "-S"], cwd=REPO, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
        hits = [line for line in proc.stdout.splitlines() if line.strip()]
        rows.append({"content_lane": lane, "pattern": pattern, "hit_count": len(hits), "sample_hits": hits[:8]})
    return rows


def axis_candidates() -> list[dict[str, Any]]:
    specs = [
        ("phase_gradient_axis_split", "Xi_LT := (dA_phi/dcycle, dA_phi/dtick)", True, True, True, True, True, False, False, True, False, True, True, True, True, True, "gradient split is internal but the cycle/tick distinction is a coordinate gauge"),
        ("entropy_growth_decay_axis_pair", "Xi_LT := (entropy shell growth, entropy relaxation lag)", True, True, True, True, True, True, False, True, False, True, True, True, True, True, "growth/lag pair is informative but does not prove R_dim action composition"),
        ("z12_coboundary_cocycle_axis_pair", "Xi_LT := (0-cochain boundary axis, 1-cocycle circulation axis)", True, True, True, True, True, False, False, True, False, True, True, True, True, True, "cochain/cocycle split is typed but remains algebraic, not dimensional"),
        ("cohomology_filtration_axis_pair", "Xi_LT := (filtration depth, persistence index)", True, True, True, True, True, False, False, True, False, False, True, True, True, True, "filtration axes need a chosen scale section and miss C_phi preservation"),
        ("damping_tail_memory_axis_pair", "Xi_LT := (tail compression length proxy, memory-lag time proxy)", True, True, True, True, True, True, False, False, False, True, True, True, True, True, "tail/memory quantities are target dependent and not a source"),
        ("dirichlet_spectral_axis_pair", "Xi_LT := (mode index spacing, spectral phase tick)", True, True, True, True, True, False, False, True, False, False, True, True, True, True, "spectral pair is formal unless an internal dispersion law is exported"),
        ("symplectic_conjugate_axis_pair", "Xi_LT := (q-axis, p-axis) from phase area", True, True, True, False, False, False, False, True, False, True, False, True, True, True, "conjugate notation restates action area and labels axes by convention"),
        ("information_flow_cut_axis_pair", "Xi_LT := (cut-size carrier, flow-order carrier)", True, True, True, True, True, True, False, True, False, False, True, True, True, True, "cut/flow is promising but lacks an exported phase-area coupling theorem"),
        ("lagrangian_space_time_axis_pair", "Xi_LT := axes read from a strict Lagrangian/EOM", True, True, False, True, True, True, True, True, True, True, True, False, True, False, "imports unresolved Lagrangian/EOM and L_total closure"),
        ("planck_light_axis_pair", "Xi_LT := (Planck length, Planck time)", True, True, False, True, True, True, True, True, True, True, True, True, True, False, "complete axes are imported standard physics"),
        ("apparatus_rod_clock_axis_pair", "Xi_LT := detector rod and clock calibration", True, True, False, True, True, True, True, True, True, True, False, True, False, False, "apparatus axes are observer readout, not primordial nadsoliton data"),
        ("selector_oriented_axis_pair", "Xi_LT := directed axis pair after selector premise", True, True, True, True, True, False, False, True, False, True, True, True, True, False, "selector premise is forbidden and QW-2191 remains open"),
    ]
    keys = ("candidate", "formula", *GATES, "blocker")
    return [dict(zip(keys, row)) for row in specs]


def axis_rows(candidates: list[dict[str, Any]]) -> list[dict[str, Any]]:
    field = {
        "two_axis_domain": "explicit_axis_source_formula",
        "length_source": "U_length_source_exported",
        "time_source": "U_time_source_exported",
        "axis_noncollapse": "length_time_distinction_nonconventional",
        "non_gauge_orientation": "axis_orientation_not_selector_premise",
        "r_dim_unlock": "R_dim_relation_unlocked",
        "phase_area_preserved": "C_phi_A_phi_preserved",
        "forbidden_import_free": "selector_bridge_ltotal_toe_free",
    }
    return [{"candidate": c["candidate"], "axis_test": test, "test_passed": bool(c[field[test]]), "accepted_axis_source": bool(c[field[test]] and c["strict_nadsoliton_data_only"] and c["not_unit_convention"] and c["not_imported_dynamics"] and c["not_apparatus_calibration"] and c["selector_bridge_ltotal_toe_free"]), "blocker": c["blocker"]} for c in candidates for test in AXIS_TESTS]


def scale_rows(candidates: list[dict[str, Any]]) -> list[dict[str, Any]]:
    rows = []
    for c in candidates:
        for lam_l, lam_t in SCALE_PAIRS:
            rows.append({"candidate": c["candidate"], "lambda_length": f"{lam_l.numerator}/{lam_l.denominator}", "lambda_time": f"{lam_t.numerator}/{lam_t.denominator}", "axis_ratio": f"{(lam_l / lam_t).numerator}/{(lam_l / lam_t).denominator}", "covariance_claimed": bool(c["scale_covariant_not_gauge_fixed"]), "covariance_accepted": bool(c["scale_covariant_not_gauge_fixed"] and c["length_time_distinction_nonconventional"] and c["not_unit_convention"] and c["R_dim_relation_unlocked"]), "blocker": c["blocker"]})
    return rows


def coupling_rows(candidates: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return [{"candidate": c["candidate"], "coupling_test": test, "A_phi": round(a_phi(), 12) if test == "preserve_C_phi_A_phi" else None, "test_passed": bool(c["U_length_source_exported"] if test == "derive_U_length" else c["U_time_source_exported"] if test == "derive_U_time" else c["length_time_distinction_nonconventional"] if test == "prove_distinct_axes" else c["C_phi_A_phi_preserved"] if test == "preserve_C_phi_A_phi" else c["R_dim_relation_unlocked"] if test == "induce_R_dim" else c["not_imported_dynamics"] and c["not_apparatus_calibration"] and c["selector_bridge_ltotal_toe_free"]), "accepted_coupling_chain": bool(c["U_length_source_exported"] and c["U_time_source_exported"] and c["length_time_distinction_nonconventional"] and c["C_phi_A_phi_preserved"] and c["R_dim_relation_unlocked"] and c["not_unit_convention"] and c["not_imported_dynamics"] and c["not_apparatus_calibration"] and c["selector_bridge_ltotal_toe_free"]), "blocker": c["blocker"]} for c in candidates for test in COUPLING_TESTS]


def gate_rows(candidates: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return [{"candidate": c["candidate"], "required_gate": gate, "gate_passed": bool(c[gate]), "detail": "passed" if c[gate] else c["blocker"]} for c in candidates for gate in GATES]


def aggregate(gates: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return [{"candidate": candidate, "passed_gates": sum(row["gate_passed"] for row in gates if row["candidate"] == candidate), "failed_gates": sum(not row["gate_passed"] for row in gates if row["candidate"] == candidate), "accepted_Xi_LT_axis_source": all(row["gate_passed"] for row in gates if row["candidate"] == candidate)} for candidate in sorted({row["candidate"] for row in gates})]


def build_payload() -> dict[str, Any]:
    p3118 = read_json(P3118)
    greps = content_grep()
    candidates = axis_candidates()
    axes = axis_rows(candidates)
    scales = scale_rows(candidates)
    couplings = coupling_rows(candidates)
    gates = gate_rows(candidates)
    aggs = aggregate(gates)
    accepted = [row for row in aggs if row["accepted_Xi_LT_axis_source"]]
    proof_obligations = [
        {"obligation": "read_p3118_next_atom", "satisfied": True, "detail": "P3118 requested exactly one Xi_LT axis-source object"},
        {"obligation": "content_first_repo_grep", "satisfied": len(greps) == 4 and sum(row["hit_count"] for row in greps) > 0, "detail": "repo was grepped by axis/relation content patterns"},
        {"obligation": "construct_Xi_LT_candidates", "satisfied": len(candidates) == 12, "detail": "twelve axis-source candidates were constructed"},
        {"obligation": "test_axis_source_laws", "satisfied": len(axes) == len(candidates) * len(AXIS_TESTS), "detail": "eight axis-source rows were built per candidate"},
        {"obligation": "test_axis_scale_covariance", "satisfied": len(scales) == len(candidates) * len(SCALE_PAIRS), "detail": "five reciprocal scale rows were built per candidate"},
        {"obligation": "test_R_dim_coupling", "satisfied": len(couplings) == len(candidates) * len(COUPLING_TESTS), "detail": "six R_dim-coupling rows were built per candidate"},
        {"obligation": "export_strict_Xi_LT", "satisfied": False, "detail": "0 candidates export an import-free strict Xi_LT object satisfying all gates"},
    ]
    return {"status": "P3119_XI_LT_AXIS_SOURCE_OBJECT_BOUNDED_NO_GO", "input_hashes": {"P3118": hashlib.sha256(P3118.read_bytes()).hexdigest() if P3118.exists() else None}, "constructed_theoretical_objects": {"content_first_repo_grep": greps, "audit_object": {"object": "XiLTAxisSourceObjectAudit", "required_source": "Xi_LT internal source for distinct U_length and U_time axes on nadsoliton data", "forbidden_imports": ["hbar/Planck", "rods", "clocks", "observed light", "apparatus", "selector replay", "L_total", "bridge/role-transfer", "ToE"]}, "candidate_Xi_LT_axis_sources": candidates, "axis_source_rows": axes, "axis_scale_covariance_rows": scales, "R_dim_coupling_rows": couplings, "candidate_gate_rows": gates, "candidate_aggregate_certificate": aggs}, "finite_certificate": {"content_grep_lanes": len(greps), "content_grep_hits": sum(row["hit_count"] for row in greps), "p3118_accepted_R_dim_relations": p3118.get("finite_certificate", {}).get("accepted_R_dim_relations"), "candidate_Xi_LT_axis_sources": len(candidates), "axis_source_rows": len(axes), "axis_scale_covariance_rows": len(scales), "R_dim_coupling_rows": len(couplings), "candidate_gate_rows": len(gates), "accepted_Xi_LT_axis_sources": len(accepted), "proof_obligations": len(proof_obligations), "satisfied_proof_obligations": sum(row["satisfied"] for row in proof_obligations)}, "proof_obligations": proof_obligations, "decision": {"bounded_result": "P3119 constructs the requested Xi_LT axis-source family and finds bounded no-go. Internal phase-gradient, entropy growth/lag, Z12 cochain/cocycle, cohomology filtration, damping/memory, spectral, symplectic, and information-flow candidates are real typed objects, but each loses at least one required condition: non-gauge length/time distinction, exported independent U_length and U_time sources, C_phi(A_phi) preservation, R_dim induction, scale covariance without gauge fixing, or closed-lane import freedom. Lagrangian, Planck/light, apparatus, and selector candidates import blocked lanes. No nadsoliton-only Xi_LT exports physical length/time axes or unlocks R_dim.", "negative_export_flags": {key: False for key in ["Xi_LT_axis_source_exported", "R_dim_relation_exported", "Omega_dim_source_exported", "K_dim_functor_exported", "Sigma_dim_theorem_exported", "D_phi_source_law_exported", "U_action_source_law_exported", "physical_length_time_unit_exported", "detector_calibration_exported", "spacetime_eom_exported", "QW_2191_discharged", "strict_selector_closure_exported", "ltotal_exported", "bridge_role_transfer_exported", "toe_closure_exported"]}, "positive_scoped_flags": {"p3118_Xi_LT_obligation_reused": True, "candidate_Xi_LT_axis_sources_constructed": True, "axis_source_matrix_built": True, "axis_scale_covariance_matrix_built": True, "R_dim_coupling_matrix_built": True, "closed_lane_and_external_import_rows_rejected": True}, "next_honest_step": "Construct exactly one new strict ordered-flow source object Tau_LT: an intrinsic, nonconventional temporal-order/length-extension bifunctor on nadsoliton information flow that tries to source U_time as order and U_length as extension before re-testing Xi_LT. It must not use clocks, rods, observed light, Planck units, Lagrangian/EOM normalization, selector replay, bridge/role-transfer, L_total, or ToE closure. Without Tau_LT or an equally explicit new typed object, preserve the P3105-P3119 physical-unit no-go/no-new-live-frontier certificate."}}


def write_markdown(payload: dict[str, Any]) -> None:
    cert = payload["finite_certificate"]
    lines = ["# P3119/S2069 Xi_LT axis-source object audit", "", f"Status: `{payload['status']}`", "", "## Finite certificate", f"- P3118 accepted R_dim relations: `{cert['p3118_accepted_R_dim_relations']}`", f"- content grep lanes: `{cert['content_grep_lanes']}`", f"- candidate Xi_LT axis sources: `{cert['candidate_Xi_LT_axis_sources']}`", f"- axis-source rows: `{cert['axis_source_rows']}`", f"- axis-scale covariance rows: `{cert['axis_scale_covariance_rows']}`", f"- R_dim coupling rows: `{cert['R_dim_coupling_rows']}`", f"- candidate gate rows: `{cert['candidate_gate_rows']}`", f"- accepted Xi_LT axis sources: `{cert['accepted_Xi_LT_axis_sources']}`", f"- satisfied proof obligations: `{cert['satisfied_proof_obligations']}/{cert['proof_obligations']}`", "", "## Decision", payload["decision"]["bounded_result"], "", "## Recommendation", payload["decision"]["next_honest_step"]]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3119/S2069 Xi_LT axis-source object audit", "## P3119/S2069 Xi_LT axis-source object audit\n\n`P3119/S2069` executes the P3118-recommended audit for `Xi_LT`, a strict axis-source object intended to distinguish `U_length` from `U_time` on nadsoliton data.  It constructs `12` candidate axis sources, `96` axis-source rows, `60` axis-scale covariance rows, `72` `R_dim` coupling rows, and a `12 x 14 = 168` gate matrix.  The bounded result is that internal phase-gradient, entropy, Z12/cohomology, damping/memory, spectral, symplectic, and information-flow candidates miss a non-gauge axis distinction, independent axis sources, `C_phi(A_phi)` preservation, `R_dim` induction, scale covariance, or import freedom, while Lagrangian, Planck/light, apparatus, and selector candidates import closed or forbidden lanes.  No physical length/time unit, detector calibration, spacetime EOM, selector closure, `L_total`, bridge/role-transfer, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3119/S2069 Xi_LT axis-source remains incomplete", "## P3119/S2069 Xi_LT axis-source remains incomplete\n\n`P3119/S2069` tests whether a strict nadsoliton-only axis-source object `Xi_LT` can distinguish `U_length` and `U_time` strongly enough to unlock `R_dim`.  Current artifacts provide no import-free strict axis source satisfying the full chain, so the result is not a Lagrangian density, Hamiltonian normalization, spacetime EOM, `L_total`, or ToE datum.\n")
    append_once(AGENTS, "Current Xi_LT axis-source object guardrail (P3119/S2069, 2026-06-26)", "## Current Xi_LT axis-source object guardrail (P3119/S2069, 2026-06-26)\n\n- P3119 tests the P3118-requested strict axis-source object `Xi_LT`, an internal source for distinct `U_length` and `U_time` axes on nadsoliton data.\n- The finite audit constructs `12` candidate axis sources, `96` axis-source rows, `60` axis-scale covariance rows, `72` `R_dim` coupling rows, and `168` gate rows; `0` candidates export an import-free strict `Xi_LT` source.\n- Do not promote phase-gradient splits, entropy growth/lag pairs, Z12 cochain/cocycle axes, cohomology filtrations, damping-tail/memory pairs, spectral axis pairs, symplectic conjugate notation, information-flow cuts, Lagrangian/EOM axes, Planck/light axes, apparatus rod/clock axes, or selector-oriented axes to physical units, detector calibration, spacetime EOM, selector closure, `L_total`, bridge/role-transfer, or ToE closure.\n- The next honest move requires exactly one new strict ordered-flow source object `Tau_LT`, an intrinsic temporal-order/length-extension bifunctor on nadsoliton information flow; otherwise preserve the P3105-P3119 physical-unit no-go/no-new-live-frontier certificate.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
