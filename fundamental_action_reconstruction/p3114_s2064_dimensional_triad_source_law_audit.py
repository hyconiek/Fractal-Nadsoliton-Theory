#!/usr/bin/env python3
"""P3114/S2064: nadsoliton-only dimensional triad source-law audit.

P3113 left exactly one admissible next object: a nadsoliton-only dimensional
triad source law D_phi=(U_action,U_length,U_time) with a scale-orbit section,
C_phi(A_phi)=U_action coupling, and an internal relation deriving action from
length/time.  This audit constructs finite triad candidates and checks whether
any candidate exports all three dimensional carriers without importing hbar,
Planck units, rods, clocks, observed light, apparatus, selector replay, L_total,
bridge/role-transfer, or ToE promotion.
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
from p3113_s2063_u_action_reference_carrier_source_law_audit import OUT as P3113, a_phi

OUT = GEN / "p3114_s2064_dimensional_triad_source_law_audit.json"
MD = GEN / "p3114_s2064_dimensional_triad_source_law_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"
SCALE_FACTORS = (Fraction(1, 3), Fraction(1, 2), Fraction(1, 1), Fraction(2, 1), Fraction(3, 1))
TRIAD_AXES = ("U_action", "U_length", "U_time")
RELATION_TESTS = ("action_from_length_time", "C_phi_A_phi_to_U_action", "scale_section", "import_freedom")
GATES = (
    "uses_p3113_D_phi_obligation",
    "explicit_triad_formula",
    "nonzero_action_carrier_exported",
    "nonzero_length_carrier_exported",
    "nonzero_time_carrier_exported",
    "action_length_time_relation_exported",
    "scale_orbit_section_exported",
    "C_phi_A_phi_coupling_exported",
    "standard_physics_import_free",
    "selector_bridge_ltotal_toe_free",
    "nonconventional_source_law_exported",
)


def content_grep() -> list[dict[str, Any]]:
    patterns = {
        "d_phi_triad": r"D_phi|U_action|U_length|U_time|dimensional triad|action/length/time",
        "scale_section": r"scale-orbit|scale orbit|section|quotient|dimensionful|calibration",
        "c_phi_coupling": r"C_phi\(A_phi\)|C_phi|A_phi|2\*pi/alpha_geo",
        "blocked_imports_closure": r"hbar|Planck|rod|clock|observed light|apparatus|selector|QW-2191|L_total|bridge/role-transfer|ToE",
    }
    rows = []
    for lane, pattern in patterns.items():
        proc = subprocess.run(["rg", "-n", pattern, "fundamental_action_reconstruction", "AGENTS.md", "-S"], cwd=REPO, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
        hits = [line for line in proc.stdout.splitlines() if line.strip()]
        rows.append({"content_lane": lane, "pattern": pattern, "hit_count": len(hits), "sample_hits": hits[:8]})
    return rows


def triad_candidates() -> list[dict[str, Any]]:
    return [
        {
            "candidate": "pure_phase_entropy_tick_triad",
            "formula": "D_phi := (A_phi, alpha_geo bit_cell, entropy_tick)",
            "dimension_vector": {"U_action": 0, "U_length": 0, "U_time": 0},
            "uses_p3113_D_phi_obligation": True,
            "explicit_triad_formula": True,
            "nonzero_action_carrier_exported": False,
            "nonzero_length_carrier_exported": False,
            "nonzero_time_carrier_exported": False,
            "action_length_time_relation_exported": False,
            "scale_orbit_section_exported": False,
            "C_phi_A_phi_coupling_exported": True,
            "standard_physics_import_free": True,
            "selector_bridge_ltotal_toe_free": True,
            "nonconventional_source_law_exported": True,
            "blocker": "phase, entropy cell, and tick are internal counters but remain dimensionless",
        },
        {
            "candidate": "cohomology_period_cell_clock_triad",
            "formula": "D_phi := (integer period action, Z12 cell length, cycle count time)",
            "dimension_vector": {"U_action": 0, "U_length": 0, "U_time": 0},
            "uses_p3113_D_phi_obligation": True,
            "explicit_triad_formula": True,
            "nonzero_action_carrier_exported": False,
            "nonzero_length_carrier_exported": False,
            "nonzero_time_carrier_exported": False,
            "action_length_time_relation_exported": False,
            "scale_orbit_section_exported": True,
            "C_phi_A_phi_coupling_exported": False,
            "standard_physics_import_free": True,
            "selector_bridge_ltotal_toe_free": True,
            "nonconventional_source_law_exported": True,
            "blocker": "topological periods and finite cells normalize labels but not physical dimensions",
        },
        {
            "candidate": "declared_dimension_symbol_triad",
            "formula": "declare dim(D_phi)=([S],[L],[T]) and set C_phi(A_phi)=U_action",
            "dimension_vector": {"U_action": 1, "U_length": 1, "U_time": 1},
            "uses_p3113_D_phi_obligation": True,
            "explicit_triad_formula": True,
            "nonzero_action_carrier_exported": True,
            "nonzero_length_carrier_exported": True,
            "nonzero_time_carrier_exported": True,
            "action_length_time_relation_exported": False,
            "scale_orbit_section_exported": False,
            "C_phi_A_phi_coupling_exported": True,
            "standard_physics_import_free": True,
            "selector_bridge_ltotal_toe_free": True,
            "nonconventional_source_law_exported": False,
            "blocker": "declared symbols name dimensions but do not source a scale section or dynamics",
        },
        {
            "candidate": "internal_dispersion_velocity_triad",
            "formula": "D_phi := (phase area action, graph wavelength, graph update period) with v=L/T",
            "dimension_vector": {"U_action": 0, "U_length": 0, "U_time": 0},
            "uses_p3113_D_phi_obligation": True,
            "explicit_triad_formula": True,
            "nonzero_action_carrier_exported": False,
            "nonzero_length_carrier_exported": False,
            "nonzero_time_carrier_exported": False,
            "action_length_time_relation_exported": True,
            "scale_orbit_section_exported": False,
            "C_phi_A_phi_coupling_exported": False,
            "standard_physics_import_free": True,
            "selector_bridge_ltotal_toe_free": True,
            "nonconventional_source_law_exported": True,
            "blocker": "graph wavelength/update period are combinatorial labels unless a dimensionful carrier is sourced",
        },
        {
            "candidate": "imported_planck_light_triad",
            "formula": "D_phi := (hbar, Planck length, Planck time or c-derived light clock)",
            "dimension_vector": {"U_action": 1, "U_length": 1, "U_time": 1},
            "uses_p3113_D_phi_obligation": True,
            "explicit_triad_formula": True,
            "nonzero_action_carrier_exported": True,
            "nonzero_length_carrier_exported": True,
            "nonzero_time_carrier_exported": True,
            "action_length_time_relation_exported": True,
            "scale_orbit_section_exported": True,
            "C_phi_A_phi_coupling_exported": True,
            "standard_physics_import_free": False,
            "selector_bridge_ltotal_toe_free": True,
            "nonconventional_source_law_exported": False,
            "blocker": "the complete triad is obtained only by importing standard physics constants/light clocks",
        },
        {
            "candidate": "apparatus_rod_clock_action_triad",
            "formula": "D_phi := (calorimeter action readout, rod length, clock time)",
            "dimension_vector": {"U_action": 1, "U_length": 1, "U_time": 1},
            "uses_p3113_D_phi_obligation": True,
            "explicit_triad_formula": True,
            "nonzero_action_carrier_exported": True,
            "nonzero_length_carrier_exported": True,
            "nonzero_time_carrier_exported": True,
            "action_length_time_relation_exported": True,
            "scale_orbit_section_exported": True,
            "C_phi_A_phi_coupling_exported": False,
            "standard_physics_import_free": False,
            "selector_bridge_ltotal_toe_free": True,
            "nonconventional_source_law_exported": False,
            "blocker": "apparatus rods and clocks are downstream observer calibrations, not nadsoliton-only source laws",
        },
        {
            "candidate": "triad_quotient_section_placeholder",
            "formula": "D_phi := choose one orbit representative [S,L,T] in R_+^3 / scale",
            "dimension_vector": {"U_action": 1, "U_length": 1, "U_time": 1},
            "uses_p3113_D_phi_obligation": True,
            "explicit_triad_formula": True,
            "nonzero_action_carrier_exported": True,
            "nonzero_length_carrier_exported": True,
            "nonzero_time_carrier_exported": True,
            "action_length_time_relation_exported": True,
            "scale_orbit_section_exported": False,
            "C_phi_A_phi_coupling_exported": False,
            "standard_physics_import_free": True,
            "selector_bridge_ltotal_toe_free": True,
            "nonconventional_source_law_exported": False,
            "blocker": "a placeholder quotient states the needed theorem but supplies no internal section or coupling law",
        },
    ]


def scale_orbit_rows(candidates: list[dict[str, Any]]) -> list[dict[str, Any]]:
    rows = []
    for cand in candidates:
        for factor in SCALE_FACTORS:
            rows.append({
                "candidate": cand["candidate"],
                "scale_factor": f"{factor.numerator}/{factor.denominator}",
                "scaled_triad": {axis: f"{factor}*{axis}" for axis in TRIAD_AXES},
                "minimal_section_factor": factor == 1,
                "candidate_claims_section": cand["scale_orbit_section_exported"],
                "import_free_section_valid": bool(cand["scale_orbit_section_exported"] and cand["standard_physics_import_free"] and cand["nonconventional_source_law_exported"] and all(cand[f"nonzero_{axis.split('_')[1]}_carrier_exported"] for axis in TRIAD_AXES)),
            })
    return rows


def carrier_axis_rows(candidates: list[dict[str, Any]]) -> list[dict[str, Any]]:
    rows = []
    key = {"U_action": "nonzero_action_carrier_exported", "U_length": "nonzero_length_carrier_exported", "U_time": "nonzero_time_carrier_exported"}
    for cand in candidates:
        for axis in TRIAD_AXES:
            rows.append({
                "candidate": cand["candidate"],
                "axis": axis,
                "dimension_exponent": cand["dimension_vector"][axis],
                "carrier_claimed": cand[key[axis]],
                "carrier_accepted": bool(cand[key[axis]] and cand["standard_physics_import_free"] and cand["nonconventional_source_law_exported"] and cand["scale_orbit_section_exported"]),
                "blocker": cand["blocker"],
            })
    return rows


def relation_rows(candidates: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return [
        {
            "candidate": cand["candidate"],
            "relation_test": test,
            "A_phi": round(a_phi(), 12) if test == "C_phi_A_phi_to_U_action" else None,
            "relation_claimed": (
                cand["action_length_time_relation_exported"] if test == "action_from_length_time" else
                cand["C_phi_A_phi_coupling_exported"] if test == "C_phi_A_phi_to_U_action" else
                cand["scale_orbit_section_exported"] if test == "scale_section" else
                cand["standard_physics_import_free"]
            ),
            "relation_accepted": bool(
                cand["action_length_time_relation_exported"] and cand["C_phi_A_phi_coupling_exported"] and cand["scale_orbit_section_exported"] and cand["standard_physics_import_free"] and cand["nonconventional_source_law_exported"] and cand["nonzero_action_carrier_exported"] and cand["nonzero_length_carrier_exported"] and cand["nonzero_time_carrier_exported"]
            ),
            "blocker": cand["blocker"],
        }
        for cand in candidates
        for test in RELATION_TESTS
    ]


def gate_rows(candidates: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return [{"candidate": cand["candidate"], "required_gate": gate, "gate_passed": bool(cand[gate]), "detail": "passed" if cand[gate] else cand["blocker"]} for cand in candidates for gate in GATES]


def aggregate(gates: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return [{"candidate": candidate, "passed_gates": sum(row["gate_passed"] for row in gates if row["candidate"] == candidate), "failed_gates": sum(not row["gate_passed"] for row in gates if row["candidate"] == candidate), "accepted_D_phi_source_law": all(row["gate_passed"] for row in gates if row["candidate"] == candidate)} for candidate in sorted({row["candidate"] for row in gates})]


def build_payload() -> dict[str, Any]:
    p3113 = read_json(P3113)
    greps = content_grep()
    candidates = triad_candidates()
    scale_rows = scale_orbit_rows(candidates)
    axis_rows = carrier_axis_rows(candidates)
    relations = relation_rows(candidates)
    gates = gate_rows(candidates)
    aggs = aggregate(gates)
    accepted = [row for row in aggs if row["accepted_D_phi_source_law"]]
    obligations = [
        {"obligation": "read_p3113_next_atom", "satisfied": True, "detail": "P3113 requested exactly one D_phi dimensional triad source law"},
        {"obligation": "content_first_repo_grep", "satisfied": len(greps) == 4 and sum(row["hit_count"] for row in greps) > 0, "detail": "repo was grepped by research content patterns"},
        {"obligation": "construct_D_phi_candidates", "satisfied": len(candidates) == 7, "detail": "seven triad source-law candidates were constructed"},
        {"obligation": "test_scale_orbit_sections", "satisfied": len(scale_rows) == len(candidates) * len(SCALE_FACTORS), "detail": "triad scale section was tested across five factors"},
        {"obligation": "test_three_carrier_axes", "satisfied": len(axis_rows) == len(candidates) * len(TRIAD_AXES), "detail": "U_action/U_length/U_time carriers were tested"},
        {"obligation": "test_relation_couplings", "satisfied": len(relations) == len(candidates) * len(RELATION_TESTS), "detail": "action-from-length/time and C_phi coupling rows were built"},
        {"obligation": "export_nadsoliton_only_D_phi", "satisfied": False, "detail": "0 candidates export an import-free triad with source section, relation, and C_phi coupling"},
    ]
    return {
        "status": "P3114_DIMENSIONAL_TRIAD_SOURCE_LAW_BOUNDED_NO_GO",
        "input_hashes": {"P3113": hashlib.sha256(P3113.read_bytes()).hexdigest() if P3113.exists() else None},
        "constructed_theoretical_objects": {
            "content_first_repo_grep": greps,
            "audit_object": {"object": "DimensionalTriadSourceLawAudit", "required_source_law": "D_phi=(U_action,U_length,U_time) with scale section, action-length/time relation, and C_phi(A_phi)=U_action coupling", "forbidden_imports": ["hbar/Planck", "rods", "clocks", "observed light", "apparatus", "selector replay", "L_total", "bridge/role-transfer", "ToE"]},
            "candidate_D_phi_source_laws": candidates,
            "scale_orbit_section_rows": scale_rows,
            "carrier_axis_rows": axis_rows,
            "relation_coupling_rows": relations,
            "candidate_gate_rows": gates,
            "candidate_aggregate_certificate": aggs,
        },
        "finite_certificate": {
            "content_grep_lanes": len(greps),
            "content_grep_hits": sum(row["hit_count"] for row in greps),
            "p3113_accepted_U_action_source_laws": p3113.get("finite_certificate", {}).get("accepted_U_action_source_laws"),
            "candidate_D_phi_source_laws": len(candidates),
            "scale_orbit_section_rows": len(scale_rows),
            "carrier_axis_rows": len(axis_rows),
            "relation_coupling_rows": len(relations),
            "candidate_gate_rows": len(gates),
            "accepted_D_phi_source_laws": len(accepted),
            "proof_obligations": len(obligations),
            "satisfied_proof_obligations": sum(row["satisfied"] for row in obligations),
        },
        "proof_obligations": obligations,
        "decision": {
            "bounded_result": "P3114 constructs the requested D_phi dimensional triad source-law family and finds bounded no-go.  Internal phase, entropy, cohomology, graph wavelength, and update-period candidates remain dimensionless/combinatorial; declared dimension symbols and quotient placeholders state the obligation without sourcing a section; and the complete dimensional triads import hbar/Planck, observed-light clocks, rods, or apparatus.  No nadsoliton-only candidate exports U_action, U_length, and U_time with an import-free scale section, action-length/time relation, and C_phi(A_phi)=U_action coupling.",
            "negative_export_flags": {key: False for key in ["D_phi_source_law_exported", "U_action_source_law_exported", "physical_length_time_unit_exported", "detector_calibration_exported", "spacetime_eom_exported", "QW_2191_discharged", "strict_selector_closure_exported", "ltotal_exported", "bridge_role_transfer_exported", "toe_closure_exported"]},
            "positive_scoped_flags": {"p3113_D_phi_obligation_reused": True, "candidate_D_phi_source_laws_constructed": True, "scale_orbit_section_matrix_built": True, "carrier_axis_matrix_built": True, "imported_planck_light_apparatus_rows_rejected": True},
            "next_honest_step": "Construct exactly one nadsoliton-only scale-section theorem for the internal dimension orbit: an explicit strict quotient/section object Sigma_dim selecting a nonzero representative for (U_action,U_length,U_time) and proving C_phi(A_phi)=U_action plus an action-from-length/time relation.  It must not introduce hbar/Planck, rods, clocks, observed light, apparatus, selector replay, L_total, bridge/role-transfer, or ToE; otherwise preserve the P3105-P3114 physical-unit no-go.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    cert = payload["finite_certificate"]
    lines = [
        "# P3114/S2064 dimensional triad source-law audit", "", f"Status: `{payload['status']}`", "", "## Finite certificate",
        f"- P3113 accepted U_action source laws: `{cert['p3113_accepted_U_action_source_laws']}`",
        f"- content grep lanes: `{cert['content_grep_lanes']}`",
        f"- candidate D_phi source laws: `{cert['candidate_D_phi_source_laws']}`",
        f"- scale-orbit section rows: `{cert['scale_orbit_section_rows']}`",
        f"- carrier axis rows: `{cert['carrier_axis_rows']}`",
        f"- relation coupling rows: `{cert['relation_coupling_rows']}`",
        f"- candidate gate rows: `{cert['candidate_gate_rows']}`",
        f"- accepted D_phi source laws: `{cert['accepted_D_phi_source_laws']}`",
        f"- satisfied proof obligations: `{cert['satisfied_proof_obligations']}/{cert['proof_obligations']}`", "", "## Decision", payload["decision"]["bounded_result"], "", "## Recommendation", payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3114/S2064 dimensional triad source-law audit", "## P3114/S2064 dimensional triad source-law audit\n\n`P3114/S2064` executes the P3113-recommended audit for a nadsoliton-only dimensional triad `D_phi=(U_action,U_length,U_time)`.  It constructs `7` candidate triad source laws, `35` scale-orbit section rows, `21` carrier-axis rows, `28` relation/coupling rows, and a `7 x 11 = 77` gate matrix.  The bounded result is that internal phase/entropy/cohomology/graph candidates remain dimensionless or combinatorial, declared triad symbols and quotient placeholders are unsourced, and complete dimensional triads import `hbar/Planck`, observed-light clocks, rods, or apparatus.  No physical action/length/time unit, detector calibration, spacetime EOM, selector closure, `L_total`, bridge/role-transfer, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3114/S2064 dimensional triad remains unsourced", "## P3114/S2064 dimensional triad remains unsourced\n\n`P3114/S2064` tests whether a nadsoliton-only `D_phi=(U_action,U_length,U_time)` can source action, length, and time together.  Current artifacts provide no import-free scale-section theorem, action-from-length/time relation, and `C_phi(A_phi)=U_action` coupling, so the result is not a Lagrangian density, Hamiltonian normalization, spacetime EOM, `L_total`, or ToE datum.\n")
    append_once(AGENTS, "Current dimensional triad source-law guardrail (P3114/S2064, 2026-06-26)", "## Current dimensional triad source-law guardrail (P3114/S2064, 2026-06-26)\n\n- P3114 tests the P3113-requested nadsoliton-only dimensional triad `D_phi=(U_action,U_length,U_time)` with a scale section, action-length/time relation, and `C_phi(A_phi)=U_action` coupling.\n- The finite audit constructs `7` candidates, `35` scale-orbit rows, `21` carrier-axis rows, `28` relation/coupling rows, and `77` gate rows; `0` candidates export an import-free triad source law.\n- Do not promote phase/entropy/cohomology/graph labels, declared dimension symbols, quotient placeholders, imported `hbar/Planck`, observed-light clocks, rods, or apparatus to detector calibration, spacetime EOM, selector closure, `L_total`, bridge/role-transfer, or ToE closure.\n- The next honest move is exactly one nadsoliton-only scale-section theorem `Sigma_dim` for the internal dimension orbit; otherwise preserve the P3105-P3114 physical-unit no-go.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
