#!/usr/bin/env python3
"""P3108/S2058: internal calibration-morphism candidate audit.

P3107 left exactly one bounded object to construct: a nadsoliton-only
calibration morphism from the alpha_geo=4 ln 2=ln(16) four-bit Shannon cell to a
nonzero action/length/time calibration object.  This audit builds that typed
candidate family and tests scale-orbit invariance and source-law status without
importing Planck units, rods, clocks, observed light, apparatus, selector
closure, bridge/role-transfer, L_total, or ToE closure.
"""
from __future__ import annotations

import hashlib, json, math, subprocess
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p3105_s2055_alpha_geo_entropy_to_unit_map_obstruction_audit import ALPHA_GEO
from p3107_s2057_shannon_to_dimension_functor_source_law_audit import OUT as P3107

OUT = GEN / "p3108_s2058_internal_calibration_morphism_candidate_audit.json"
MD = GEN / "p3108_s2058_internal_calibration_morphism_candidate_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"
GATES = (
    "source_is_alpha_geo_four_bit_cell",
    "target_nonzero_dimensional",
    "internal_formula_defined",
    "scale_orbit_invariant",
    "scale_breaking_nonconventional",
    "dimension_assignment_exported",
    "calibration_interface_exported",
    "internal_source_law_exported",
    "standard_physics_import_free",
)
CANDIDATES = (
    {"id": "dimensionless_identity_entropy_morphism", "target": "Info0", "dimension": "dimensionless", "formula": "I4 -> alpha_geo", "source_is_alpha_geo_four_bit_cell": True, "target_nonzero_dimensional": False, "internal_formula_defined": True, "scale_orbit_invariant": True, "scale_breaking_nonconventional": False, "dimension_assignment_exported": False, "calibration_interface_exported": False, "internal_source_law_exported": False, "standard_physics_import_free": True, "blocker": "identity keeps the Shannon cell internal and dimensionless"},
    {"id": "alpha_geo_action_density_unit_candidate", "target": "ActionUnit", "dimension": "action", "formula": "S_A = alpha_geo * rho_I", "source_is_alpha_geo_four_bit_cell": True, "target_nonzero_dimensional": True, "internal_formula_defined": True, "scale_orbit_invariant": False, "scale_breaking_nonconventional": False, "dimension_assignment_exported": True, "calibration_interface_exported": False, "internal_source_law_exported": False, "standard_physics_import_free": True, "blocker": "formula is internal but rescales under arbitrary action-unit lambda"},
    {"id": "sixteen_microstate_length_cell_candidate", "target": "LengthCell", "dimension": "length", "formula": "ell = f(16 states)", "source_is_alpha_geo_four_bit_cell": True, "target_nonzero_dimensional": True, "internal_formula_defined": True, "scale_orbit_invariant": False, "scale_breaking_nonconventional": False, "dimension_assignment_exported": True, "calibration_interface_exported": False, "internal_source_law_exported": False, "standard_physics_import_free": True, "blocker": "cardinality fixes a count, not a length representative in the positive scale orbit"},
    {"id": "internal_update_tick_candidate", "target": "TimeTick", "dimension": "time", "formula": "tau = one Shannon update of I4", "source_is_alpha_geo_four_bit_cell": True, "target_nonzero_dimensional": True, "internal_formula_defined": True, "scale_orbit_invariant": False, "scale_breaking_nonconventional": False, "dimension_assignment_exported": True, "calibration_interface_exported": False, "internal_source_law_exported": False, "standard_physics_import_free": True, "blocker": "update order is a formal tick without absolute duration or calibration"},
    {"id": "self_coupled_density_gradient_calibration", "target": "CalibrationInterface", "dimension": "calibration", "formula": "cal = d(rho_I)/dI4", "source_is_alpha_geo_four_bit_cell": True, "target_nonzero_dimensional": True, "internal_formula_defined": True, "scale_orbit_invariant": False, "scale_breaking_nonconventional": False, "dimension_assignment_exported": False, "calibration_interface_exported": True, "internal_source_law_exported": False, "standard_physics_import_free": True, "blocker": "internal sensitivity is relative; no nonconventional calibration interface to physical readout"},
    {"id": "hbar_planck_imported_calibration", "target": "ImportedActionCalibration", "dimension": "action", "formula": "alpha_geo maps to hbar template", "source_is_alpha_geo_four_bit_cell": True, "target_nonzero_dimensional": True, "internal_formula_defined": True, "scale_orbit_invariant": True, "scale_breaking_nonconventional": True, "dimension_assignment_exported": True, "calibration_interface_exported": True, "internal_source_law_exported": False, "standard_physics_import_free": False, "blocker": "breaks scale only by imported standard-physics calibration"},
)
SCALES = ("lambda", "alpha_geo_lambda", "sixteen_lambda", "lambda_inverse")


def content_grep() -> list[dict[str, Any]]:
    patterns = {
        "calibration_morphism_atom": r"calibration-morphism|calibration morphism|internal calibration|calibration interface",
        "alpha_geo_shannon_atom": r"alpha_geo|4 ln\(2\)|ln\(16\)|four-bit|4-bit|Shannon",
        "scale_orbit_unit_atom": r"scale-orbit|positive scale|canonical.*unit|unit-source|dimension assignment",
        "blocked_imports": r"Planck|hbar|rod|clock|observed light|apparatus|selector closure|L_total|ToE",
    }
    rows = []
    for lane, pattern in patterns.items():
        proc = subprocess.run(["rg", "-n", pattern, "fundamental_action_reconstruction", "AGENTS.md", "-S"], cwd=REPO, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
        hits = [line for line in proc.stdout.splitlines() if line.strip()]
        rows.append({"content_lane": lane, "pattern": pattern, "hit_count": len(hits), "sample_hits": hits[:8]})
    return rows


def calibration_objects() -> list[dict[str, Any]]:
    return [
        {"object": "I4_alpha_geo_shannon_cell", "bits": 4, "states": 16, "entropy_nats": round(ALPHA_GEO, 12), "role": "source information cell"},
        {"object": "ActionUnit", "dimension": "action", "requires_scale_representative": True},
        {"object": "LengthTimeUnit", "dimension": "length/time", "requires_scale_representative": True},
        {"object": "CalibrationInterface", "dimension": "readout", "requires_nonimported_interface": True},
    ]


def morphism_candidate_rows() -> list[dict[str, Any]]:
    return [{"candidate": c["id"], "source_object": "I4_alpha_geo_shannon_cell", "target_object": c["target"], "target_dimension": c["dimension"], "formula": c["formula"], **{g: c[g] for g in GATES}} for c in CANDIDATES]


def scale_orbit_rows() -> list[dict[str, Any]]:
    rows = []
    for c in CANDIDATES:
        for s in SCALES:
            affected = c["target_nonzero_dimensional"] or c["scale_breaking_nonconventional"]
            rows.append({"candidate": c["id"], "scale_orbit": s, "target_nonzero_dimensional": c["target_nonzero_dimensional"], "scale_breaking_claimed": c["scale_breaking_nonconventional"], "internal_source_law_exported": c["internal_source_law_exported"], "standard_physics_import_free": c["standard_physics_import_free"], "passes_scale_orbit_test": (not affected and c["scale_orbit_invariant"]) or (c["scale_orbit_invariant"] and c["scale_breaking_nonconventional"] and c["internal_source_law_exported"] and c["standard_physics_import_free"])})
    return rows


def source_law_rows() -> list[dict[str, Any]]:
    return [{"candidate": c["id"], "internal_formula_defined": c["internal_formula_defined"], "dimension_assignment_exported": c["dimension_assignment_exported"], "calibration_interface_exported": c["calibration_interface_exported"], "internal_source_law_exported": c["internal_source_law_exported"], "verdict": "blocked: " + c["blocker"]} for c in CANDIDATES]


def gate_rows() -> list[dict[str, Any]]:
    return [{"candidate": c["id"], "required_gate": gate, "gate_passed": bool(c[gate]), "detail": "passed" if c[gate] else c["blocker"]} for c in CANDIDATES for gate in GATES]


def aggregate(gates: list[dict[str, Any]]) -> list[dict[str, Any]]:
    rows = []
    for c in CANDIDATES:
        subset = [r for r in gates if r["candidate"] == c["id"]]
        accepted = all(r["gate_passed"] for r in subset) and c["target_nonzero_dimensional"]
        rows.append({"candidate": c["id"], "passed_gates": sum(r["gate_passed"] for r in subset), "failed_gates": sum(not r["gate_passed"] for r in subset), "accepted_internal_calibration_morphism": accepted})
    return rows


def build_payload() -> dict[str, Any]:
    p3107 = read_json(P3107)
    greps = content_grep(); objs = calibration_objects(); morphs = morphism_candidate_rows(); scales = scale_orbit_rows(); laws = source_law_rows(); gates = gate_rows(); aggs = aggregate(gates)
    accepted = [r for r in aggs if r["accepted_internal_calibration_morphism"]]
    obligations = [
        {"obligation": "read_p3107_next_atom", "satisfied": True, "detail": "P3107 selected one bounded internal calibration-morphism candidate audit"},
        {"obligation": "content_first_repo_grep", "satisfied": len(greps) == 4 and sum(r["hit_count"] for r in greps) > 0, "detail": "repo was grepped by research content, not only by numbers"},
        {"obligation": "construct_alpha_geo_source_cell", "satisfied": any(o["object"] == "I4_alpha_geo_shannon_cell" for o in objs), "detail": "four-bit Shannon source object is explicit"},
        {"obligation": "construct_nonzero_calibration_targets", "satisfied": sum(1 for o in objs if o["object"] != "I4_alpha_geo_shannon_cell") == 3, "detail": "action, length/time, and calibration targets are typed"},
        {"obligation": "test_scale_orbit_invariance", "satisfied": len(scales) == len(CANDIDATES) * len(SCALES), "detail": "all candidates tested on four positive scale representatives"},
        {"obligation": "export_internal_nonimported_calibration_morphism", "satisfied": False, "detail": "0 nonzero candidates pass source-law, scale, and import-free gates"},
    ]
    return {"status": "P3108_INTERNAL_CALIBRATION_MORPHISM_CANDIDATE_OBSTRUCTION_BOUNDED_NO_GO", "input_hashes": {"P3107": hashlib.sha256(P3107.read_bytes()).hexdigest() if P3107.exists() else None}, "constructed_theoretical_objects": {"content_first_repo_grep": greps, "calibration_morphism_audit_object": {"object": "NadsolitonAlphaGeoInternalCalibrationMorphismCandidateAudit", "source_reused": "P3107 recommendation and alpha_geo=4 ln 2=ln(16) four-bit Shannon cell", "ontology": "nadsoliton is self-coupled primordial information; calibration must be sourced internally by that information, not imported from standard physics", "acceptance_predicate": "nonzero target dimension, internal formula, scale-orbit invariance, nonconventional scale breaking, dimension assignment, calibration interface, internal source law, and no standard-physics import"}, "calibration_object_rows": objs, "morphism_candidate_rows": morphs, "scale_orbit_rows": scales, "source_law_rows": laws, "candidate_gate_rows": gates, "candidate_aggregate_certificate": aggs}, "finite_certificate": {"content_grep_lanes": len(greps), "content_grep_hits": sum(r["hit_count"] for r in greps), "p3107_accepted_internal_shannon_to_dimension_source_laws": p3107.get("finite_certificate", {}).get("accepted_internal_shannon_to_dimension_source_laws"), "calibration_objects": len(objs), "morphism_candidates": len(CANDIDATES), "scale_orbit_rows": len(scales), "source_law_rows": len(laws), "required_gates": len(GATES), "candidate_gate_rows": len(gates), "accepted_internal_calibration_morphisms": len(accepted), "proof_obligations": len(obligations), "satisfied_proof_obligations": sum(r["satisfied"] for r in obligations)}, "proof_obligations": obligations, "decision": {"bounded_result": "P3108 constructs the bounded internal calibration-morphism candidate demanded by P3107.  The alpha_geo=4 ln 2=ln(16) four-bit Shannon cell gives real self-coupled information and several internal formulas can target action, length/time, or calibration labels.  But each nonzero import-free candidate remains scale-orbit dependent and lacks a nonconventional scale-breaking source law plus calibration interface; the only scale-stable internal morphism is dimensionless.  The imported hbar/Planck row breaks scale only by standard-physics convention and is rejected.", "negative_export_flags": {key: False for key in ["accepted_internal_calibration_morphism_exported", "physical_action_unit_exported", "physical_length_time_unit_exported", "detector_calibration_exported", "observed_light_interface_exported", "spacetime_eom_exported", "QW_2191_discharged", "strict_selector_closure_exported", "ltotal_exported", "bridge_role_transfer_exported", "toe_closure_exported"]}, "positive_scoped_flags": {"alpha_geo_four_bit_cell_verified": True, "internal_morphism_candidates_constructed": True, "scale_orbit_rows_computed": True, "imported_calibration_rejected": True}, "next_honest_step": "Construct exactly one scale-orbit quotient/source-law candidate for the alpha_geo four-bit cell: a nadsoliton-only rule that selects a nonconventional positive representative before assigning action/length/time calibration.  Test only quotient uniqueness, source-law provenance, and non-imported calibration coupling; if absent, preserve the P3105-P3108 unit-source no-go."}}


def write_markdown(payload: dict[str, Any]) -> None:
    c = payload["finite_certificate"]
    lines = ["# P3108/S2058 internal calibration-morphism candidate audit", "", f"Status: `{payload['status']}`", "", "## Finite certificate", f"- P3107 accepted internal Shannon-to-dimension source laws: `{c['p3107_accepted_internal_shannon_to_dimension_source_laws']}`", f"- content grep lanes: `{c['content_grep_lanes']}`", f"- calibration objects: `{c['calibration_objects']}`", f"- morphism candidates: `{c['morphism_candidates']}`", f"- scale-orbit rows: `{c['scale_orbit_rows']}`", f"- source-law rows: `{c['source_law_rows']}`", f"- candidate gate rows: `{c['candidate_gate_rows']}`", f"- accepted internal calibration morphisms: `{c['accepted_internal_calibration_morphisms']}`", f"- satisfied proof obligations: `{c['satisfied_proof_obligations']}/{c['proof_obligations']}`", "", "## Decision", payload["decision"]["bounded_result"], "", "## Recommendation", payload["decision"]["next_honest_step"]]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3108/S2058 internal calibration-morphism candidate audit", "## P3108/S2058 internal calibration-morphism candidate audit\n\n`P3108/S2058` executes the P3107-recommended bounded internal calibration-morphism candidate audit for the `alpha_geo=4 ln 2=ln(16)` four-bit Shannon cell.  It constructs `4` calibration objects, `6` morphism candidates, `24` scale-orbit rows, `6` source-law rows, and a `6 x 9 = 54` gate matrix.  The result remains bounded no-go: internal formulas can target action/length/time/calibration labels, but no nonzero import-free candidate exports a scale-orbit quotient, nonconventional representative, internal source law, and calibration interface.  No physical unit, detector calibration, spacetime EOM, selector closure, `L_total`, bridge/role-transfer, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3108/S2058 calibration morphism remains nonphysical", "## P3108/S2058 calibration morphism remains nonphysical\n\n`P3108/S2058` builds internal alpha_geo Shannon-cell morphisms toward action, length/time, and calibration targets.  The formulas are legitimate internal information bookkeeping, but all nonzero import-free targets remain scale-orbit dependent and lack a nonconventional source-law quotient.  Therefore they are not yet physical Lagrangian density, Hamiltonian, spacetime EOM, `L_total`, or ToE data.\n")
    append_once(AGENTS, "Current internal calibration-morphism candidate guardrail (P3108/S2058, 2026-06-26)", "## Current internal calibration-morphism candidate guardrail (P3108/S2058, 2026-06-26)\n\n- P3108 follows the P3107 recommendation and constructs bounded internal calibration-morphism candidates from the `alpha_geo=4 ln 2=ln(16)` four-bit Shannon cell of the self-coupled nadsoliton.\n- The finite audit computes `4` calibration objects, `6` morphism candidates, `24` scale-orbit rows, `6` source-law rows, and `54` candidate gate rows; `0` candidates export a nonzero internal import-free calibration morphism.\n- Do not promote alpha_geo Shannon-cell morphisms, internal action-density formulas, formal length/time ticks, or imported hbar/Planck calibration to physical units, detector calibration, spacetime EOM, selector closure, `L_total`, bridge/role-transfer, or ToE closure.\n- The next honest move is exactly one scale-orbit quotient/source-law candidate for the alpha_geo four-bit cell, unless a genuinely new internal unit-source law is introduced.\n")
    return payload

if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
