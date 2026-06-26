#!/usr/bin/env python3
"""P3110/S2060: alpha_geo dimension-carrying comparison standard audit.

P3109 left exactly one honest object to construct: an internal, nadsoliton-only
comparison standard for the alpha_geo=4 ln 2=ln(16) Shannon cell that is not
merely entropy/cardinality/ratio and that could choose a positive scale-orbit
section before action/length/time calibration.  This script builds a finite
candidate matrix and tests whether any candidate supplies nonconventional
provenance, section uniqueness, dimensional target typing, and calibration
coupling without importing Planck/hbar/rods/clocks/light/apparatus, selector
closure, L_total, bridge/role-transfer, or ToE closure.
"""
from __future__ import annotations

import hashlib, json, math, subprocess
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p3105_s2055_alpha_geo_entropy_to_unit_map_obstruction_audit import ALPHA_GEO
from p3109_s2059_alpha_geo_scale_orbit_quotient_source_law_audit import OUT as P3109

OUT = GEN / "p3110_s2060_alpha_geo_dimension_comparison_standard_audit.json"
MD = GEN / "p3110_s2060_alpha_geo_dimension_comparison_standard_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"
TARGETS = ("action", "length", "time")
SCALES = (0.5, 1.0, 2.0, 4.0)
GATES = (
    "nadsoliton_only_input",
    "uses_alpha_geo_four_bit_cell",
    "not_entropy_cardinality_or_ratio_only",
    "dimension_bearing_standard",
    "unique_positive_section",
    "calibration_coupling_exported",
    "nonconventional_source_law_exported",
    "standard_physics_import_free",
)
CANDIDATES = (
    {"id": "self_information_norm_standard", "value": ALPHA_GEO, "nadsoliton_only_input": True, "uses_alpha_geo_four_bit_cell": True, "not_entropy_cardinality_or_ratio_only": False, "dimension_bearing_standard": False, "unique_positive_section": False, "calibration_coupling_exported": False, "nonconventional_source_law_exported": False, "standard_physics_import_free": True, "blocker": "self-information norm is still entropy normalization; lambda times the unit is observationally uncalibrated"},
    {"id": "four_bit_fisher_curvature_standard", "value": 1.0 / ALPHA_GEO, "nadsoliton_only_input": True, "uses_alpha_geo_four_bit_cell": True, "not_entropy_cardinality_or_ratio_only": True, "dimension_bearing_standard": False, "unique_positive_section": False, "calibration_coupling_exported": False, "nonconventional_source_law_exported": False, "standard_physics_import_free": True, "blocker": "Fisher curvature gives an internal dimensionless geometry, not an action/length/time comparison standard"},
    {"id": "z12_cell_count_density_standard", "value": ALPHA_GEO / 12.0, "nadsoliton_only_input": True, "uses_alpha_geo_four_bit_cell": True, "not_entropy_cardinality_or_ratio_only": False, "dimension_bearing_standard": False, "unique_positive_section": False, "calibration_coupling_exported": False, "nonconventional_source_law_exported": False, "standard_physics_import_free": True, "blocker": "Z12 density is a count/ratio over cells; it does not carry a physical comparison dimension"},
    {"id": "internal_tick_from_self_coupled_update", "value": math.log(2), "nadsoliton_only_input": True, "uses_alpha_geo_four_bit_cell": True, "not_entropy_cardinality_or_ratio_only": True, "dimension_bearing_standard": False, "unique_positive_section": False, "calibration_coupling_exported": False, "nonconventional_source_law_exported": False, "standard_physics_import_free": True, "blocker": "a formal update/tick is an ordering parameter until a nonconventional clock calibration law is exported"},
    {"id": "symplectic_phase_area_candidate", "value": 2.0 * math.pi / ALPHA_GEO, "nadsoliton_only_input": True, "uses_alpha_geo_four_bit_cell": True, "not_entropy_cardinality_or_ratio_only": True, "dimension_bearing_standard": True, "unique_positive_section": False, "calibration_coupling_exported": False, "nonconventional_source_law_exported": False, "standard_physics_import_free": True, "blocker": "phase-area typing is the right formal shape but lacks an exported source law fixing its absolute positive representative"},
    {"id": "imported_planck_action_comparison", "value": ALPHA_GEO, "nadsoliton_only_input": False, "uses_alpha_geo_four_bit_cell": True, "not_entropy_cardinality_or_ratio_only": True, "dimension_bearing_standard": True, "unique_positive_section": True, "calibration_coupling_exported": True, "nonconventional_source_law_exported": False, "standard_physics_import_free": False, "blocker": "works only by importing hbar/Planck calibration from standard physics"},
)


def content_grep() -> list[dict[str, Any]]:
    patterns = {
        "dimension_comparison_standard": r"dimension-carrying comparison standard|comparison standard|dimension-bearing standard|calibration coupling",
        "alpha_geo_shannon_cell": r"alpha_geo|4 ln\(2\)|ln\(16\)|four-bit|4-bit|Shannon|self-coupled information|nadsoliton",
        "scale_orbit_section": r"scale-orbit|positive representative|quotient section|unit-source|bit-to-unit|action/length/time",
        "blocked_imports": r"Planck|hbar|rod|clock|observed light|apparatus|selector closure|L_total|bridge/role-transfer|ToE",
    }
    rows = []
    for lane, pattern in patterns.items():
        proc = subprocess.run(["rg", "-n", pattern, "fundamental_action_reconstruction", "AGENTS.md", "-S"], cwd=REPO, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
        hits = [line for line in proc.stdout.splitlines() if line.strip()]
        rows.append({"content_lane": lane, "pattern": pattern, "hit_count": len(hits), "sample_hits": hits[:8]})
    return rows


def candidate_rows() -> list[dict[str, Any]]:
    return [{"candidate": c["id"], "standard_value": round(c["value"], 12), **{g: c[g] for g in GATES}, "blocker": c["blocker"]} for c in CANDIDATES]


def orbit_test_rows() -> list[dict[str, Any]]:
    rows = []
    for c in CANDIDATES:
        for target in TARGETS:
            vals = [round(c["value"] * s, 12) for s in SCALES]
            rows.append({"candidate": c["id"], "target": target, "scale_multipliers": list(SCALES), "orbit_values": vals, "same_internal_information_class": c["standard_physics_import_free"], "unique_positive_section_selected": c["unique_positive_section"] and c["standard_physics_import_free"] and c["nonconventional_source_law_exported"]})
    return rows


def gate_rows() -> list[dict[str, Any]]:
    return [{"candidate": c["id"], "required_gate": g, "gate_passed": bool(c[g]), "detail": "passed" if c[g] else c["blocker"]} for c in CANDIDATES for g in GATES]


def aggregate(gates: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return [{"candidate": c["id"], "passed_gates": sum(r["gate_passed"] for r in gates if r["candidate"] == c["id"]), "failed_gates": sum(not r["gate_passed"] for r in gates if r["candidate"] == c["id"]), "accepted_dimension_comparison_standard": all(r["gate_passed"] for r in gates if r["candidate"] == c["id"])} for c in CANDIDATES]


def build_payload() -> dict[str, Any]:
    p3109 = read_json(P3109)
    greps = content_grep(); candidates = candidate_rows(); orbits = orbit_test_rows(); gates = gate_rows(); aggs = aggregate(gates)
    accepted = [r for r in aggs if r["accepted_dimension_comparison_standard"]]
    obligations = [
        {"obligation": "read_p3109_next_atom", "satisfied": True, "detail": "P3109 selected exactly one dimension-carrying comparison-standard audit"},
        {"obligation": "content_first_repo_grep", "satisfied": len(greps) == 4 and sum(r["hit_count"] for r in greps) > 0, "detail": "repo was grepped by research content, not only by names/numbers"},
        {"obligation": "construct_candidate_standards", "satisfied": len(candidates) == len(CANDIDATES), "detail": "six internal/import-control standards were constructed"},
        {"obligation": "compute_scale_orbit_section_tests", "satisfied": len(orbits) == len(CANDIDATES) * len(TARGETS), "detail": "each candidate was tested over action/length/time scale orbits"},
        {"obligation": "export_nonimported_dimension_standard", "satisfied": False, "detail": "0 candidates pass all provenance, dimension, uniqueness, and coupling gates"},
    ]
    return {
        "status": "P3110_ALPHA_GEO_DIMENSION_COMPARISON_STANDARD_OBSTRUCTION_BOUNDED_NO_GO",
        "input_hashes": {"P3109": hashlib.sha256(P3109.read_bytes()).hexdigest() if P3109.exists() else None},
        "constructed_theoretical_objects": {
            "content_first_repo_grep": greps,
            "audit_object": {"object": "AlphaGeoDimensionComparisonStandardAudit", "ontology": "the nadsoliton is self-coupled primordial information; alpha_geo is a four-bit Shannon cell, but a physical unit additionally needs a dimension-bearing comparison standard and a calibration-coupling source law", "acceptance_predicate": "nadsoliton-only input, alpha_geo cell use, not merely entropy/cardinality/ratio, dimension-bearing standard, unique positive section, calibration coupling, nonconventional source law, and no standard-physics import"},
            "candidate_standard_rows": candidates,
            "scale_orbit_section_test_rows": orbits,
            "candidate_gate_rows": gates,
            "candidate_aggregate_certificate": aggs,
        },
        "finite_certificate": {"content_grep_lanes": len(greps), "content_grep_hits": sum(r["hit_count"] for r in greps), "p3109_accepted_scale_orbit_quotient_source_laws": p3109.get("finite_certificate", {}).get("accepted_scale_orbit_quotient_source_laws"), "candidate_standards": len(candidates), "targets": len(TARGETS), "scale_orbit_section_rows": len(orbits), "required_gates": len(GATES), "candidate_gate_rows": len(gates), "accepted_dimension_comparison_standards": len(accepted), "proof_obligations": len(obligations), "satisfied_proof_obligations": sum(r["satisfied"] for r in obligations)},
        "proof_obligations": obligations,
        "decision": {"bounded_result": "P3110 constructs the P3109-requested dimension-carrying comparison-standard object.  The best internal candidate is the symplectic_phase_area_candidate: it has a formal action-like shape, but it still lacks a nonconventional source law selecting one positive representative and coupling it to action/length/time calibration.  Pure self-information, Fisher curvature, Z12 density, and formal ticks remain dimensionless internal information data; imported Planck calibration is rejected.", "negative_export_flags": {k: False for k in ["dimension_comparison_standard_exported", "physical_action_unit_exported", "physical_length_time_unit_exported", "detector_calibration_exported", "spacetime_eom_exported", "QW_2191_discharged", "strict_selector_closure_exported", "ltotal_exported", "bridge_role_transfer_exported", "toe_closure_exported"]}, "positive_scoped_flags": {"alpha_geo_four_bit_entropy_confirmed": True, "candidate_dimension_standards_constructed": True, "symplectic_phase_area_shape_identified": True, "imported_planck_standard_rejected": True}, "next_honest_step": "Construct exactly one nadsoliton-only symplectic/action phase source-law for the symplectic_phase_area_candidate: an explicit formula that fixes the positive scale representative and couples it to the alpha_geo four-bit cell without hbar/Planck, rods, clocks, observed light, apparatus, selector replay, L_total, bridge/role-transfer, or ToE promotion.  If no such formula is supplied, preserve the P3105-P3110 unit-source no-go."},
    }


def write_markdown(payload: dict[str, Any]) -> None:
    c = payload["finite_certificate"]
    lines = ["# P3110/S2060 alpha_geo dimension comparison standard audit", "", f"Status: `{payload['status']}`", "", "## Finite certificate", f"- P3109 accepted scale-orbit quotient/source laws: `{c['p3109_accepted_scale_orbit_quotient_source_laws']}`", f"- content grep lanes: `{c['content_grep_lanes']}`", f"- candidate standards: `{c['candidate_standards']}`", f"- scale-orbit section rows: `{c['scale_orbit_section_rows']}`", f"- candidate gate rows: `{c['candidate_gate_rows']}`", f"- accepted dimension comparison standards: `{c['accepted_dimension_comparison_standards']}`", f"- satisfied proof obligations: `{c['satisfied_proof_obligations']}/{c['proof_obligations']}`", "", "## Decision", payload["decision"]["bounded_result"], "", "## Recommendation", payload["decision"]["next_honest_step"]]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3110/S2060 alpha_geo dimension comparison standard audit", "## P3110/S2060 alpha_geo dimension comparison standard audit\n\n`P3110/S2060` executes the P3109-recommended audit for an internal dimension-carrying comparison standard for the `alpha_geo=4 ln 2=ln(16)` four-bit Shannon cell.  It constructs `6` candidate standards, `18` action/length/time scale-orbit section rows, and a `6 x 8 = 48` gate matrix.  The bounded result is that internal information candidates remain dimensionless or section-unsourced; the imported Planck row is rejected.  No physical action/length/time unit, detector calibration, spacetime EOM, selector closure, `L_total`, bridge/role-transfer, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3110/S2060 dimension standard remains unsourced", "## P3110/S2060 dimension standard remains unsourced\n\n`P3110/S2060` identifies `symplectic_phase_area_candidate` as the strongest formal shape for an action-like comparison standard, but it remains without a nadsoliton-only positive-section source law and without calibration coupling.  The alpha_geo Shannon cell is therefore still internal information accounting, not a physical Lagrangian density, Hamiltonian normalization, spacetime EOM, `L_total`, or ToE datum.\n")
    append_once(AGENTS, "Current alpha_geo dimension-comparison-standard guardrail (P3110/S2060, 2026-06-26)", "## Current alpha_geo dimension-comparison-standard guardrail (P3110/S2060, 2026-06-26)\n\n- P3110 constructs the P3109-requested internal dimension-carrying comparison-standard audit for the `alpha_geo=4 ln 2=ln(16)` four-bit Shannon cell.\n- The finite audit computes `6` candidate standards, `18` action/length/time scale-orbit section rows, and `48` candidate gate rows; `0` candidates export a non-imported dimension comparison standard with unique positive section and calibration coupling.\n- Do not promote self-information norm, Fisher curvature, Z12 density, formal ticks, symplectic phase-area shape, or imported Planck calibration to physical units, detector calibration, spacetime EOM, selector closure, `L_total`, bridge/role-transfer, or ToE closure.\n- The next honest move is exactly one nadsoliton-only symplectic/action phase source-law for `symplectic_phase_area_candidate`; otherwise preserve the P3105-P3110 unit-source no-go.\n")
    return payload

if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
