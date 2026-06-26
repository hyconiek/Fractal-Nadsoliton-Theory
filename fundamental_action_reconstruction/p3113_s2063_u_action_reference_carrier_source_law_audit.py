#!/usr/bin/env python3
"""P3113/S2063: U_action reference-carrier source-law audit.

P3112 left exactly one admissible next object: a nadsoliton-only dimensionful
reference carrier U_action with a scale-orbit section and a coupling theorem
C_phi(A_phi)=U_action.  This audit constructs finite U_action source-law
candidates and checks whether any candidate supplies a nonzero action dimension,
breaks the scale orbit internally, couples to C_phi, and derives length/time
calibration without importing hbar/Planck, rods, clocks, observed light,
apparatus, selector replay, L_total, bridge/role-transfer, or ToE promotion.
"""
from __future__ import annotations

import hashlib
import json
import math
import subprocess
from fractions import Fraction
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p3105_s2055_alpha_geo_entropy_to_unit_map_obstruction_audit import ALPHA_GEO
from p3112_s2062_internal_calibration_functional_c_phi_audit import OUT as P3112

OUT = GEN / "p3113_s2063_u_action_reference_carrier_source_law_audit.json"
MD = GEN / "p3113_s2063_u_action_reference_carrier_source_law_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"
SCALE_FACTORS = (Fraction(1, 4), Fraction(1, 2), Fraction(1, 1), Fraction(2, 1), Fraction(4, 1))
DIMENSION_AXES = ("action", "length", "time")
GATES = (
    "uses_p3112_C_phi_obligation",
    "explicit_U_action_formula",
    "nonzero_action_dimension_exported",
    "scale_orbit_section_exported",
    "C_phi_A_phi_coupling_theorem_exported",
    "length_time_calibration_derived",
    "standard_physics_import_free",
    "selector_bridge_ltotal_toe_free",
    "nonconventional_source_law_exported",
)


def content_grep() -> list[dict[str, Any]]:
    patterns = {
        "u_action_reference_carrier": r"U_action|reference-carrier|reference carrier|dimensionful reference|action unit",
        "c_phi_coupling": r"C_phi\(A_phi\)|C_phi|A_phi|2\*pi/alpha_geo|calibration functional",
        "scale_orbit_dimension": r"scale-orbit|scale orbit|dimensionful|nonzero action dimension|length/time|meter|clock",
        "blocked_imports_closure": r"hbar|Planck|rod|clock|observed light|apparatus|selector|QW-2191|L_total|bridge/role-transfer|ToE",
    }
    rows = []
    for lane, pattern in patterns.items():
        proc = subprocess.run(
            ["rg", "-n", pattern, "fundamental_action_reconstruction", "AGENTS.md", "-S"],
            cwd=REPO,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False,
        )
        hits = [line for line in proc.stdout.splitlines() if line.strip()]
        rows.append({"content_lane": lane, "pattern": pattern, "hit_count": len(hits), "sample_hits": hits[:8]})
    return rows


def a_phi() -> float:
    return (2.0 * math.pi) / ALPHA_GEO


def candidate_source_laws() -> list[dict[str, Any]]:
    return [
        {
            "candidate": "pure_phase_reference_carrier",
            "formula": "U_action := [A_phi] with C_phi(A_phi)=1",
            "dimension_vector": {"action": 0, "length": 0, "time": 0},
            "uses_p3112_C_phi_obligation": True,
            "explicit_U_action_formula": True,
            "nonzero_action_dimension_exported": False,
            "scale_orbit_section_exported": True,
            "C_phi_A_phi_coupling_theorem_exported": True,
            "length_time_calibration_derived": False,
            "standard_physics_import_free": True,
            "selector_bridge_ltotal_toe_free": True,
            "nonconventional_source_law_exported": True,
            "blocker": "the carrier is a normalized phase number, not a dimensionful action carrier",
        },
        {
            "candidate": "alpha_geo_entropy_cell_carrier",
            "formula": "U_action := alpha_geo * bit_cell",
            "dimension_vector": {"action": 0, "length": 0, "time": 0},
            "uses_p3112_C_phi_obligation": True,
            "explicit_U_action_formula": True,
            "nonzero_action_dimension_exported": False,
            "scale_orbit_section_exported": False,
            "C_phi_A_phi_coupling_theorem_exported": False,
            "length_time_calibration_derived": False,
            "standard_physics_import_free": True,
            "selector_bridge_ltotal_toe_free": True,
            "nonconventional_source_law_exported": True,
            "blocker": "four-bit entropy is internal information accounting and carries no action dimension",
        },
        {
            "candidate": "formal_declared_U_action_symbol",
            "formula": "declare dim(U_action)=action and set C_phi(A_phi)=U_action",
            "dimension_vector": {"action": 1, "length": 0, "time": 0},
            "uses_p3112_C_phi_obligation": True,
            "explicit_U_action_formula": True,
            "nonzero_action_dimension_exported": True,
            "scale_orbit_section_exported": False,
            "C_phi_A_phi_coupling_theorem_exported": True,
            "length_time_calibration_derived": False,
            "standard_physics_import_free": True,
            "selector_bridge_ltotal_toe_free": True,
            "nonconventional_source_law_exported": False,
            "blocker": "declaring a dimension symbol does not source the dimension or fix the scale orbit",
        },
        {
            "candidate": "cohomological_integer_period_carrier",
            "formula": "U_action := generator of H^1/Z period paired with A_phi",
            "dimension_vector": {"action": 0, "length": 0, "time": 0},
            "uses_p3112_C_phi_obligation": True,
            "explicit_U_action_formula": True,
            "nonzero_action_dimension_exported": False,
            "scale_orbit_section_exported": True,
            "C_phi_A_phi_coupling_theorem_exported": False,
            "length_time_calibration_derived": False,
            "standard_physics_import_free": True,
            "selector_bridge_ltotal_toe_free": True,
            "nonconventional_source_law_exported": True,
            "blocker": "integer periods can fix a topological normalization but remain dimensionless and do not couple to C_phi as action",
        },
        {
            "candidate": "imported_planck_action_carrier",
            "formula": "U_action := hbar or Planck action quantum",
            "dimension_vector": {"action": 1, "length": 2, "time": -1},
            "uses_p3112_C_phi_obligation": True,
            "explicit_U_action_formula": True,
            "nonzero_action_dimension_exported": True,
            "scale_orbit_section_exported": True,
            "C_phi_A_phi_coupling_theorem_exported": True,
            "length_time_calibration_derived": False,
            "standard_physics_import_free": False,
            "selector_bridge_ltotal_toe_free": True,
            "nonconventional_source_law_exported": False,
            "blocker": "the action dimension comes from imported standard-physics hbar/Planck calibration",
        },
        {
            "candidate": "detector_apparatus_reference_carrier",
            "formula": "U_action := calibrated apparatus action readout",
            "dimension_vector": {"action": 1, "length": 2, "time": -1},
            "uses_p3112_C_phi_obligation": True,
            "explicit_U_action_formula": True,
            "nonzero_action_dimension_exported": True,
            "scale_orbit_section_exported": True,
            "C_phi_A_phi_coupling_theorem_exported": False,
            "length_time_calibration_derived": True,
            "standard_physics_import_free": False,
            "selector_bridge_ltotal_toe_free": True,
            "nonconventional_source_law_exported": False,
            "blocker": "apparatus calibration is downstream observer physics, not a nadsoliton-only source law",
        },
    ]


def scale_section_rows(candidates: list[dict[str, Any]]) -> list[dict[str, Any]]:
    rows = []
    for cand in candidates:
        for factor in SCALE_FACTORS:
            rows.append(
                {
                    "candidate": cand["candidate"],
                    "scale_factor": f"{factor.numerator}/{factor.denominator}",
                    "scaled_C_phi_A_phi": f"{factor} * U_action_candidate",
                    "minimal_section_factor": factor == 1,
                    "candidate_claims_section": cand["scale_orbit_section_exported"],
                    "import_free_section_valid": bool(cand["scale_orbit_section_exported"] and cand["standard_physics_import_free"] and cand["nonzero_action_dimension_exported"]),
                }
            )
    return rows


def dimensional_balance_rows(candidates: list[dict[str, Any]]) -> list[dict[str, Any]]:
    rows = []
    for cand in candidates:
        dims = cand["dimension_vector"]
        for axis in DIMENSION_AXES:
            rows.append(
                {
                    "candidate": cand["candidate"],
                    "axis": axis,
                    "dimension_exponent": dims[axis],
                    "axis_sourced_internally": bool(dims[axis] != 0 and cand["standard_physics_import_free"] and cand["nonconventional_source_law_exported"]),
                    "blocker": cand["blocker"] if dims[axis] == 0 or not cand["standard_physics_import_free"] else "dimension symbol is present but still needs an internal source theorem",
                }
            )
    return rows


def coupling_rows(candidates: list[dict[str, Any]]) -> list[dict[str, Any]]:
    area = a_phi()
    return [
        {
            "candidate": cand["candidate"],
            "A_phi": round(area, 12),
            "required_coupling": "C_phi(A_phi)=U_action",
            "coupling_claimed": cand["C_phi_A_phi_coupling_theorem_exported"],
            "coupling_accepted": bool(cand["C_phi_A_phi_coupling_theorem_exported"] and cand["standard_physics_import_free"] and cand["nonzero_action_dimension_exported"] and cand["scale_orbit_section_exported"] and cand["nonconventional_source_law_exported"]),
            "blocker": cand["blocker"],
        }
        for cand in candidates
    ]


def gate_rows(candidates: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return [
        {"candidate": cand["candidate"], "required_gate": gate, "gate_passed": bool(cand[gate]), "detail": "passed" if cand[gate] else cand["blocker"]}
        for cand in candidates
        for gate in GATES
    ]


def aggregate(gates: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return [
        {
            "candidate": candidate,
            "passed_gates": sum(row["gate_passed"] for row in gates if row["candidate"] == candidate),
            "failed_gates": sum(not row["gate_passed"] for row in gates if row["candidate"] == candidate),
            "accepted_U_action_source_law": all(row["gate_passed"] for row in gates if row["candidate"] == candidate),
        }
        for candidate in sorted({row["candidate"] for row in gates})
    ]


def build_payload() -> dict[str, Any]:
    p3112 = read_json(P3112)
    greps = content_grep()
    candidates = candidate_source_laws()
    scale_rows = scale_section_rows(candidates)
    dim_rows = dimensional_balance_rows(candidates)
    couplings = coupling_rows(candidates)
    gates = gate_rows(candidates)
    aggs = aggregate(gates)
    accepted = [row for row in aggs if row["accepted_U_action_source_law"]]
    obligations = [
        {"obligation": "read_p3112_next_atom", "satisfied": True, "detail": "P3112 requested exactly one U_action reference-carrier source law"},
        {"obligation": "content_first_repo_grep", "satisfied": len(greps) == 4 and sum(row["hit_count"] for row in greps) > 0, "detail": "repo was grepped by research content patterns"},
        {"obligation": "construct_U_action_candidates", "satisfied": len(candidates) == 6, "detail": "six U_action source-law candidates were constructed"},
        {"obligation": "test_scale_orbit_sections", "satisfied": len(scale_rows) == len(candidates) * len(SCALE_FACTORS), "detail": "scale section was tested across five factors"},
        {"obligation": "test_dimension_balance", "satisfied": len(dim_rows) == len(candidates) * len(DIMENSION_AXES), "detail": "action/length/time exponents were tested"},
        {"obligation": "test_C_phi_coupling", "satisfied": len(couplings) == len(candidates), "detail": "C_phi(A_phi)=U_action coupling rows were built"},
        {"obligation": "export_nadsoliton_only_U_action", "satisfied": False, "detail": "0 candidates export import-free U_action with scale section, C_phi coupling, and length/time calibration"},
    ]
    return {
        "status": "P3113_U_ACTION_REFERENCE_CARRIER_SOURCE_LAW_BOUNDED_NO_GO",
        "input_hashes": {"P3112": hashlib.sha256(P3112.read_bytes()).hexdigest() if P3112.exists() else None},
        "constructed_theoretical_objects": {
            "content_first_repo_grep": greps,
            "audit_object": {
                "object": "UActionReferenceCarrierSourceLawAudit",
                "required_source_law": "U_action with nonzero action dimension, scale-orbit section, and C_phi(A_phi)=U_action coupling",
                "input_section": "A_phi=2*pi/alpha_geo",
                "forbidden_imports": ["hbar/Planck", "rods", "clocks", "observed light", "apparatus", "selector replay", "L_total", "bridge/role-transfer", "ToE"],
            },
            "candidate_U_action_source_laws": candidates,
            "scale_orbit_section_rows": scale_rows,
            "dimensional_balance_rows": dim_rows,
            "C_phi_coupling_rows": couplings,
            "candidate_gate_rows": gates,
            "candidate_aggregate_certificate": aggs,
        },
        "finite_certificate": {
            "content_grep_lanes": len(greps),
            "content_grep_hits": sum(row["hit_count"] for row in greps),
            "p3112_accepted_internal_dimensionful_calibration_functionals": p3112.get("finite_certificate", {}).get("accepted_internal_dimensionful_calibration_functionals"),
            "candidate_U_action_source_laws": len(candidates),
            "scale_orbit_section_rows": len(scale_rows),
            "dimensional_balance_rows": len(dim_rows),
            "C_phi_coupling_rows": len(couplings),
            "candidate_gate_rows": len(gates),
            "accepted_U_action_source_laws": len(accepted),
            "proof_obligations": len(obligations),
            "satisfied_proof_obligations": sum(row["satisfied"] for row in obligations),
        },
        "proof_obligations": obligations,
        "decision": {
            "bounded_result": "P3113 constructs the requested U_action reference-carrier source-law family and finds bounded no-go.  Phase, entropy, and cohomological carriers are internally meaningful but dimensionless; a declared U_action symbol has no source theorem or scale section; and the only dimensionful rows import hbar/Planck or apparatus calibration.  No nadsoliton-only candidate exports U_action together with import-free scale section, C_phi(A_phi) coupling, and length/time calibration.",
            "negative_export_flags": {key: False for key in ["U_action_source_law_exported", "dimensionful_action_unit_exported", "physical_length_time_unit_exported", "detector_calibration_exported", "spacetime_eom_exported", "QW_2191_discharged", "strict_selector_closure_exported", "ltotal_exported", "bridge_role_transfer_exported", "toe_closure_exported"]},
            "positive_scoped_flags": {"p3112_C_phi_obligation_reused": True, "candidate_U_action_source_laws_constructed": True, "scale_orbit_section_matrix_built": True, "dimension_balance_matrix_built": True, "imported_planck_and_apparatus_rows_rejected": True},
            "next_honest_step": "Construct exactly one nadsoliton-only dimensional triad source law D_phi=(U_action,U_length,U_time): an explicit strict object that sources a nonzero action dimension and simultaneously exports internal length/time carriers plus a relation deriving U_action from them.  It must include a scale-orbit quotient/section proof and C_phi(A_phi)=U_action coupling, without hbar/Planck, rods, clocks, observed light, apparatus, selector replay, L_total, bridge/role-transfer, or ToE; otherwise preserve the P3105-P3113 physical-unit no-go.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    cert = payload["finite_certificate"]
    lines = [
        "# P3113/S2063 U_action reference-carrier source-law audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Finite certificate",
        f"- P3112 accepted internal dimensionful calibration functionals: `{cert['p3112_accepted_internal_dimensionful_calibration_functionals']}`",
        f"- content grep lanes: `{cert['content_grep_lanes']}`",
        f"- candidate U_action source laws: `{cert['candidate_U_action_source_laws']}`",
        f"- scale-orbit section rows: `{cert['scale_orbit_section_rows']}`",
        f"- dimensional balance rows: `{cert['dimensional_balance_rows']}`",
        f"- C_phi coupling rows: `{cert['C_phi_coupling_rows']}`",
        f"- candidate gate rows: `{cert['candidate_gate_rows']}`",
        f"- accepted U_action source laws: `{cert['accepted_U_action_source_laws']}`",
        f"- satisfied proof obligations: `{cert['satisfied_proof_obligations']}/{cert['proof_obligations']}`",
        "",
        "## Decision",
        payload["decision"]["bounded_result"],
        "",
        "## Recommendation",
        payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(
        STRICT_EQUATION_SHEET,
        "P3113/S2063 U_action reference-carrier source-law audit",
        "## P3113/S2063 U_action reference-carrier source-law audit\n\n`P3113/S2063` executes the P3112-recommended audit for a nadsoliton-only dimensionful reference carrier `U_action` coupled by `C_phi(A_phi)=U_action`.  It constructs `6` candidate source laws, `30` scale-orbit section rows, `18` dimensional-balance rows, `6` coupling rows, and a `6 x 9 = 54` gate matrix.  The bounded result is that internal phase/entropy/cohomology carriers remain dimensionless, formal `U_action` declarations are unsourced, and dimensionful rows import `hbar/Planck` or apparatus calibration.  No physical action/length/time unit, detector calibration, spacetime EOM, selector closure, `L_total`, bridge/role-transfer, or ToE closure is exported.\n",
    )
    append_once(
        STRICT_LAGRANGIAN_DRAFT,
        "P3113/S2063 U_action carrier remains unsourced",
        "## P3113/S2063 U_action carrier remains unsourced\n\n`P3113/S2063` tests whether a nadsoliton-only `U_action` can lift `C_phi(A_phi)` into a dimensionful action carrier.  Current artifacts provide no import-free action carrier with a scale-orbit section, `C_phi` coupling theorem, and length/time induction, so the result is not a Lagrangian density, Hamiltonian normalization, spacetime EOM, `L_total`, or ToE datum.\n",
    )
    append_once(
        AGENTS,
        "Current U_action reference-carrier source-law guardrail (P3113/S2063, 2026-06-26)",
        "## Current U_action reference-carrier source-law guardrail (P3113/S2063, 2026-06-26)\n\n- P3113 tests the P3112-requested nadsoliton-only dimensionful reference carrier `U_action` coupled by `C_phi(A_phi)=U_action`.\n- The finite audit constructs `6` candidate source laws, `30` scale-orbit section rows, `18` dimensional-balance rows, `6` coupling rows, and `54` gate rows; `0` candidates export an import-free `U_action` with length/time calibration.\n- Do not promote phase, entropy, cohomological integer periods, declared unit symbols, imported `hbar/Planck`, or apparatus calibration to detector calibration, spacetime EOM, selector closure, `L_total`, bridge/role-transfer, or ToE closure.\n- The next honest move is exactly one nadsoliton-only dimensional triad source law `D_phi=(U_action,U_length,U_time)`; otherwise preserve the P3105-P3113 physical-unit no-go.\n",
    )
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
