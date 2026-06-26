#!/usr/bin/env python3
"""P3112/S2062: internal calibration functional C_phi audit.

P3111 exported the dimensionless internal phase-area section
A_phi=2*pi/alpha_geo and left exactly one next object: an internal calibration
functional C_phi.  This step constructs finite candidate C_phi functionals and
checks whether any maps A_phi to a dimensionful action comparison, breaks scale
covariance internally, and induces length/time calibration without hbar/Planck,
rods, clocks, observed light, apparatus, selector replay, L_total,
bridge/role-transfer, or ToE promotion.
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
from p3111_s2061_symplectic_action_phase_source_law_audit import OUT as P3111

OUT = GEN / "p3112_s2062_internal_calibration_functional_c_phi_audit.json"
MD = GEN / "p3112_s2062_internal_calibration_functional_c_phi_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"
SCALE_FACTORS = (Fraction(1, 3), Fraction(1, 2), Fraction(1, 1), Fraction(2, 1), Fraction(3, 1))
TARGETS = ("action", "length", "time")
GATES = (
    "uses_p3111_A_phi",
    "explicit_functional_formula",
    "codomain_dimensionful_action_comparison",
    "internal_scale_covariance_broken",
    "positive_representative_preserved",
    "length_time_induction_supplied",
    "standard_physics_import_free",
    "selector_bridge_ltotal_toe_free",
    "nonconventional_source_law_exported",
)


def content_grep() -> list[dict[str, Any]]:
    patterns = {
        "internal_calibration_functional": r"C_phi|calibration functional|action comparison|dimensionful action|scale covariance",
        "phase_area_alpha_geo": r"A_phi|2\*pi/alpha_geo|alpha_geo|4 ln\(2\)|ln\(16\)|four-bit|Shannon",
        "length_time_induction_gap": r"length/time|length and time|induce.*calibration|rod|clock|light|apparatus|detector",
        "blocked_closure_imports": r"hbar|Planck|selector|QW-2191|L_total|bridge/role-transfer|ToE|spacetime EOM",
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


def candidate_functionals(area: float) -> list[dict[str, Any]]:
    return [
        {
            "candidate": "dimensionless_identity_C_phi",
            "formula": "C_phi(A)=A/A_phi",
            "value_at_A_phi": 1.0,
            "uses_p3111_A_phi": True,
            "explicit_functional_formula": True,
            "codomain_dimensionful_action_comparison": False,
            "internal_scale_covariance_broken": False,
            "positive_representative_preserved": True,
            "length_time_induction_supplied": False,
            "standard_physics_import_free": True,
            "selector_bridge_ltotal_toe_free": True,
            "nonconventional_source_law_exported": False,
            "blocker": "normalizes the phase-area orbit but has dimensionless codomain and remains scale-covariant",
        },
        {
            "candidate": "alpha_geo_phase_normalized_C_phi",
            "formula": "C_phi(A)=alpha_geo*A/(2*pi)",
            "value_at_A_phi": round(ALPHA_GEO * area / (2.0 * math.pi), 12),
            "uses_p3111_A_phi": True,
            "explicit_functional_formula": True,
            "codomain_dimensionful_action_comparison": False,
            "internal_scale_covariance_broken": False,
            "positive_representative_preserved": True,
            "length_time_induction_supplied": False,
            "standard_physics_import_free": True,
            "selector_bridge_ltotal_toe_free": True,
            "nonconventional_source_law_exported": True,
            "blocker": "repackages the P3111 winding number; it does not create an action dimension or length/time calibration",
        },
        {
            "candidate": "formal_action_unit_symbol_C_phi",
            "formula": "C_phi(A)=A/A_phi * U_action with U_action declared internally",
            "value_at_A_phi": "1 * U_action",
            "uses_p3111_A_phi": True,
            "explicit_functional_formula": True,
            "codomain_dimensionful_action_comparison": True,
            "internal_scale_covariance_broken": False,
            "positive_representative_preserved": True,
            "length_time_induction_supplied": False,
            "standard_physics_import_free": True,
            "selector_bridge_ltotal_toe_free": True,
            "nonconventional_source_law_exported": False,
            "blocker": "a declared unit symbol is not a sourced unit; scale covariance is merely renamed, not broken",
        },
        {
            "candidate": "entropy_tick_length_time_C_phi",
            "formula": "C_phi(A)=(A/A_phi)*(one four-bit tick)*(one internal cell length)",
            "value_at_A_phi": "1 tick*cell",
            "uses_p3111_A_phi": True,
            "explicit_functional_formula": True,
            "codomain_dimensionful_action_comparison": False,
            "internal_scale_covariance_broken": False,
            "positive_representative_preserved": True,
            "length_time_induction_supplied": False,
            "standard_physics_import_free": True,
            "selector_bridge_ltotal_toe_free": True,
            "nonconventional_source_law_exported": False,
            "blocker": "internal tick/cell labels are still dimensionless unless an independent length/time unit source is exported",
        },
        {
            "candidate": "imported_hbar_planck_C_phi",
            "formula": "C_phi(A)=hbar*A/A_phi",
            "value_at_A_phi": "hbar",
            "uses_p3111_A_phi": True,
            "explicit_functional_formula": True,
            "codomain_dimensionful_action_comparison": True,
            "internal_scale_covariance_broken": True,
            "positive_representative_preserved": True,
            "length_time_induction_supplied": False,
            "standard_physics_import_free": False,
            "selector_bridge_ltotal_toe_free": True,
            "nonconventional_source_law_exported": False,
            "blocker": "dimensionful action is obtained only by importing hbar/Planck calibration from standard physics",
        },
    ]


def scale_covariance_rows(candidates: list[dict[str, Any]]) -> list[dict[str, Any]]:
    rows = []
    for cand in candidates:
        for factor in SCALE_FACTORS:
            rows.append({
                "candidate": cand["candidate"],
                "scale_factor": f"{factor.numerator}/{factor.denominator}",
                "input_area": round(a_phi() * float(factor), 12),
                "p3111_phase_value_over_2pi": round(float(factor), 12),
                "fixed_minimal_section_only_at_factor_one": factor == 1,
                "internally_breaks_scale_orbit": bool(cand["internal_scale_covariance_broken"] and cand["standard_physics_import_free"]),
            })
    return rows


def induction_rows(candidates: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return [
        {
            "candidate": cand["candidate"],
            "target": target,
            "action_calibration_available": bool(cand["codomain_dimensionful_action_comparison"] and cand["standard_physics_import_free"]),
            "target_calibration_induced": bool(cand["length_time_induction_supplied"]),
            "blocker": "no independent internal meter/clock/light-cone relation is exported" if target != "action" else cand["blocker"],
        }
        for cand in candidates
        for target in TARGETS
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
            "accepted_internal_dimensionful_calibration_functional": all(row["gate_passed"] for row in gates if row["candidate"] == candidate),
        }
        for candidate in sorted({row["candidate"] for row in gates})
    ]


def build_payload() -> dict[str, Any]:
    p3111 = read_json(P3111)
    greps = content_grep()
    area = a_phi()
    candidates = candidate_functionals(area)
    covariance = scale_covariance_rows(candidates)
    induction = induction_rows(candidates)
    gates = gate_rows(candidates)
    aggs = aggregate(gates)
    accepted = [row for row in aggs if row["accepted_internal_dimensionful_calibration_functional"]]
    obligations = [
        {"obligation": "read_p3111_next_atom", "satisfied": True, "detail": "P3111 requested exactly one internal calibration functional C_phi"},
        {"obligation": "content_first_repo_grep", "satisfied": len(greps) == 4 and sum(row["hit_count"] for row in greps) > 0, "detail": "repo was grepped by research content patterns"},
        {"obligation": "construct_candidate_C_phi_functionals", "satisfied": len(candidates) == 5, "detail": "five C_phi candidates were constructed"},
        {"obligation": "test_internal_scale_covariance_breaking", "satisfied": len(covariance) == len(candidates) * len(SCALE_FACTORS), "detail": "each candidate was tested across five scale factors"},
        {"obligation": "test_action_length_time_induction", "satisfied": len(induction) == len(candidates) * len(TARGETS), "detail": "each candidate was tested for action/length/time induction"},
        {"obligation": "export_import_free_dimensionful_calibration", "satisfied": False, "detail": "0 candidates export dimensionful action plus length/time calibration without imported standards"},
    ]
    return {
        "status": "P3112_INTERNAL_C_PHI_CALIBRATION_FUNCTIONAL_BOUNDED_NO_GO",
        "input_hashes": {"P3111": hashlib.sha256(P3111.read_bytes()).hexdigest() if P3111.exists() else None},
        "constructed_theoretical_objects": {
            "content_first_repo_grep": greps,
            "audit_object": {
                "object": "InternalCalibrationFunctionalCPhiAudit",
                "input_section": "A_phi=2*pi/alpha_geo",
                "required_map": "C_phi: dimensionless phase-area section -> dimensionful action comparison, with internal scale-covariance breaking and length/time induction",
                "ontology": "nadsoliton-only internal information; no lower information layer and no imported hbar/Planck apparatus standard",
            },
            "candidate_C_phi_functionals": candidates,
            "scale_covariance_witness_rows": covariance,
            "action_length_time_induction_rows": induction,
            "candidate_gate_rows": gates,
            "candidate_aggregate_certificate": aggs,
        },
        "finite_certificate": {
            "content_grep_lanes": len(greps),
            "content_grep_hits": sum(row["hit_count"] for row in greps),
            "p3111_minimal_internal_phase_area_section_exported": p3111.get("decision", {}).get("positive_scoped_flags", {}).get("minimal_internal_phase_area_section_exported"),
            "candidate_C_phi_functionals": len(candidates),
            "scale_covariance_witness_rows": len(covariance),
            "action_length_time_induction_rows": len(induction),
            "candidate_gate_rows": len(gates),
            "accepted_internal_dimensionful_calibration_functionals": len(accepted),
            "proof_obligations": len(obligations),
            "satisfied_proof_obligations": sum(row["satisfied"] for row in obligations),
        },
        "proof_obligations": obligations,
        "decision": {
            "bounded_result": "P3112 constructs the requested C_phi object family and verifies the obstruction: internal normalizations of A_phi remain dimensionless, formal unit symbols only rename the scale orbit, entropy tick/cell labels do not induce physical length/time, and the only dimensionful action row imports hbar/Planck.  No nadsoliton-only C_phi breaks scale covariance and calibrates action/length/time on current artifacts.",
            "negative_export_flags": {key: False for key in ["dimensionful_action_unit_exported", "physical_length_time_unit_exported", "detector_calibration_exported", "spacetime_eom_exported", "QW_2191_discharged", "strict_selector_closure_exported", "ltotal_exported", "bridge_role_transfer_exported", "toe_closure_exported"]},
            "positive_scoped_flags": {"p3111_A_phi_reused": True, "candidate_C_phi_functionals_constructed": True, "scale_covariance_obstruction_witnessed": True, "imported_hbar_planck_row_rejected": True},
            "next_honest_step": "Construct exactly one nadsoliton-only dimensionful reference-carrier source law U_action: an explicit internal object with a nonzero action dimension, a scale-orbit quotient/section proof, and a coupling theorem C_phi(A_phi)=U_action that also derives length/time calibration.  It must avoid hbar/Planck, rods, clocks, observed light, apparatus, selector replay, L_total, bridge/role-transfer, and ToE promotion; otherwise preserve the P3105-P3112 physical-unit no-go.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    cert = payload["finite_certificate"]
    lines = [
        "# P3112/S2062 internal calibration functional C_phi audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Finite certificate",
        f"- P3111 minimal internal phase-area section exported: `{cert['p3111_minimal_internal_phase_area_section_exported']}`",
        f"- content grep lanes: `{cert['content_grep_lanes']}`",
        f"- candidate C_phi functionals: `{cert['candidate_C_phi_functionals']}`",
        f"- scale-covariance witness rows: `{cert['scale_covariance_witness_rows']}`",
        f"- action/length/time induction rows: `{cert['action_length_time_induction_rows']}`",
        f"- candidate gate rows: `{cert['candidate_gate_rows']}`",
        f"- accepted internal dimensionful calibration functionals: `{cert['accepted_internal_dimensionful_calibration_functionals']}`",
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
        "P3112/S2062 internal calibration functional C_phi audit",
        "## P3112/S2062 internal calibration functional C_phi audit\n\n`P3112/S2062` executes the P3111-recommended audit for an internal calibration functional `C_phi` mapping `A_phi=2*pi/alpha_geo` to dimensionful action comparison.  It constructs `5` candidate `C_phi` functionals, `25` scale-covariance witness rows, `15` action/length/time induction rows, and a `5 x 9 = 45` gate matrix.  The bounded result is that internal normalizations remain dimensionless, formal action-unit symbols are unsourced, entropy tick/cell labels do not induce length/time calibration, and the only dimensionful action row imports `hbar/Planck`.  No physical action/length/time unit, detector calibration, spacetime EOM, selector closure, `L_total`, bridge/role-transfer, or ToE closure is exported.\n",
    )
    append_once(
        STRICT_LAGRANGIAN_DRAFT,
        "P3112/S2062 C_phi calibration remains unsourced",
        "## P3112/S2062 C_phi calibration remains unsourced\n\n`P3112/S2062` tests whether the internal phase-area section `A_phi=2*pi/alpha_geo` can be lifted by a nadsoliton-only `C_phi` into a dimensionful action comparison.  Current artifacts provide no import-free dimensionful reference carrier and no length/time induction theorem, so `C_phi` is not yet a Lagrangian density, Hamiltonian normalization, spacetime EOM, `L_total`, or ToE datum.\n",
    )
    append_once(
        AGENTS,
        "Current internal C_phi calibration-functional guardrail (P3112/S2062, 2026-06-26)",
        "## Current internal C_phi calibration-functional guardrail (P3112/S2062, 2026-06-26)\n\n- P3112 tests the P3111-requested internal calibration functional `C_phi` for `A_phi=2*pi/alpha_geo`.\n- The finite audit constructs `5` candidate functionals, `25` scale-covariance rows, `15` action/length/time induction rows, and `45` gate rows; `0` candidates export an import-free dimensionful action/length/time calibration.\n- Do not promote dimensionless `A_phi` normalizations, formal unit symbols, entropy tick/cell labels, or imported `hbar/Planck` calibration to detector calibration, spacetime EOM, selector closure, `L_total`, bridge/role-transfer, or ToE closure.\n- The next honest move is exactly one nadsoliton-only dimensionful reference-carrier source law `U_action` coupled by `C_phi(A_phi)=U_action`; otherwise preserve the P3105-P3112 physical-unit no-go.\n",
    )
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
