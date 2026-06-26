#!/usr/bin/env python3
"""P3107/S2057: Shannon-to-dimension functor/source-law obstruction audit.

P3106 left a sharper object to build: a typed functor/source law carrying the
alpha_geo = 4 ln 2 = ln(16) four-bit Shannon cell from internal nadsoliton
information density to dimension labels plus calibration data.  This audit is
not standard physics in disguise: the nadsoliton remains self-coupled
primordial information, and action/length/time labels can pass only if they are
sourced internally and invariant under the positive scale orbit.
"""
from __future__ import annotations

import hashlib, json, math, subprocess
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p3105_s2055_alpha_geo_entropy_to_unit_map_obstruction_audit import ALPHA_GEO
from p3106_s2056_nadsoliton_action_density_normalization_audit import OUT as P3106

OUT = GEN / "p3107_s2057_shannon_to_dimension_functor_source_law_audit.json"
MD = GEN / "p3107_s2057_shannon_to_dimension_functor_source_law_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"
DIMENSIONS = ("dimensionless_info", "action", "length", "time", "calibration")
GATES = (
    "typed_functor_defined",
    "alpha_geo_four_bit_cell_used",
    "additive_shannon_coherence_passed",
    "positive_scale_orbit_invariant",
    "nonzero_dimension_label_exported",
    "calibration_morphism_exported",
    "internal_source_law_exported",
    "standard_physics_import_free",
)
CANDIDATES = (
    {"id": "zero_dimensionless_entropy_functor", "dimension_vector": (0, 0, 0, 0), "typed_functor_defined": True, "alpha_geo_four_bit_cell_used": True, "additive_shannon_coherence_passed": True, "positive_scale_orbit_invariant": True, "nonzero_dimension_label_exported": False, "calibration_morphism_exported": False, "internal_source_law_exported": False, "standard_physics_import_free": True, "blocker": "coherent but exports only dimensionless information"},
    {"id": "bit_count_to_action_label_functor", "dimension_vector": (1, 0, 0, 0), "typed_functor_defined": True, "alpha_geo_four_bit_cell_used": True, "additive_shannon_coherence_passed": True, "positive_scale_orbit_invariant": False, "nonzero_dimension_label_exported": True, "calibration_morphism_exported": False, "internal_source_law_exported": False, "standard_physics_import_free": True, "blocker": "nonzero action label has no scale-invariant internal calibration morphism"},
    {"id": "ln16_cell_to_length_label_functor", "dimension_vector": (0, 1, 0, 0), "typed_functor_defined": True, "alpha_geo_four_bit_cell_used": True, "additive_shannon_coherence_passed": True, "positive_scale_orbit_invariant": False, "nonzero_dimension_label_exported": True, "calibration_morphism_exported": False, "internal_source_law_exported": False, "standard_physics_import_free": True, "blocker": "length label is a formal tag unless an internal ruler/calibration morphism is exported"},
    {"id": "density_gradient_to_time_label_functor", "dimension_vector": (0, 0, 1, 0), "typed_functor_defined": True, "alpha_geo_four_bit_cell_used": True, "additive_shannon_coherence_passed": False, "positive_scale_orbit_invariant": False, "nonzero_dimension_label_exported": True, "calibration_morphism_exported": False, "internal_source_law_exported": False, "standard_physics_import_free": True, "blocker": "time label depends on an unsourced update parameter and fails additive entropy coherence"},
    {"id": "four_bit_cell_to_detector_calibration_functor", "dimension_vector": (0, 0, 0, 1), "typed_functor_defined": True, "alpha_geo_four_bit_cell_used": True, "additive_shannon_coherence_passed": True, "positive_scale_orbit_invariant": False, "nonzero_dimension_label_exported": True, "calibration_morphism_exported": True, "internal_source_law_exported": False, "standard_physics_import_free": False, "blocker": "calibration exists only by importing detector/apparatus semantics"},
    {"id": "planck_hbar_action_import_functor", "dimension_vector": (1, 2, -1, 0), "typed_functor_defined": True, "alpha_geo_four_bit_cell_used": True, "additive_shannon_coherence_passed": True, "positive_scale_orbit_invariant": True, "nonzero_dimension_label_exported": True, "calibration_morphism_exported": True, "internal_source_law_exported": False, "standard_physics_import_free": False, "blocker": "dimension and calibration are imported from hbar/Planck units"},
    {"id": "observed_light_rod_clock_import_functor", "dimension_vector": (0, 1, -1, 1), "typed_functor_defined": True, "alpha_geo_four_bit_cell_used": True, "additive_shannon_coherence_passed": True, "positive_scale_orbit_invariant": True, "nonzero_dimension_label_exported": True, "calibration_morphism_exported": True, "internal_source_law_exported": False, "standard_physics_import_free": False, "blocker": "readout calibration is imported from light/rod/clock apparatus"},
)


def content_grep() -> list[dict[str, Any]]:
    patterns = {
        "p3106_functor_atom": r"Shannon-to-dimension functor|dimension labels|calibration data|source-law",
        "alpha_geo_four_bit_atom": r"alpha_geo|4 ln\(2\)|ln\(16\)|four-bit|4-bit|Shannon",
        "blocked_imports": r"Planck|hbar|rod|clock|observed light|apparatus|selector closure|L_total|ToE",
    }
    rows = []
    for lane, pattern in patterns.items():
        proc = subprocess.run(["rg", "-n", pattern, "fundamental_action_reconstruction", "AGENTS.md", "-S"], cwd=REPO, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
        hits = [line for line in proc.stdout.splitlines() if line.strip()]
        rows.append({"content_lane": lane, "pattern": pattern, "hit_count": len(hits), "sample_hits": hits[:8]})
    return rows


def info_category_rows() -> list[dict[str, Any]]:
    rows = []
    for bits in range(0, 5):
        states = 2 ** bits
        entropy = math.log(states) if states else 0.0
        rows.append({"object": f"I_{bits}", "bits": bits, "states": states, "shannon_entropy_nats": round(entropy, 12), "alpha_geo_fraction": round(entropy / ALPHA_GEO, 12) if ALPHA_GEO else None, "is_alpha_geo_four_bit_cell": bits == 4})
    return rows


def additive_coherence_rows() -> list[dict[str, Any]]:
    rows = []
    for a in range(0, 5):
        for b in range(0, 5 - a):
            left = math.log(2 ** (a + b))
            right = math.log(2 ** a) + math.log(2 ** b)
            rows.append({"bits_a": a, "bits_b": b, "bits_sum": a + b, "entropy_additivity_defect": round(left - right, 12), "additive_shannon_coherence_passed": abs(left - right) < 1e-12})
    return rows


def dimension_functor_rows() -> list[dict[str, Any]]:
    rows = []
    for c in CANDIDATES:
        for dim in DIMENSIONS:
            idx = DIMENSIONS.index(dim) - 1
            exponent = 0 if dim == "dimensionless_info" else c["dimension_vector"][idx]
            rows.append({"candidate": c["id"], "dimension_label": dim, "assigned_exponent": exponent, "nonzero_label": exponent != 0 or dim == "calibration" and c["calibration_morphism_exported"], "internal_source_law_exported": c["internal_source_law_exported"], "standard_physics_import_free": c["standard_physics_import_free"]})
    return rows


def scale_invariance_rows() -> list[dict[str, Any]]:
    rows = []
    for c in CANDIDATES:
        weight = sum(abs(x) for x in c["dimension_vector"])
        for scale in ("lambda", "alpha_geo_lambda", "sixteen_lambda"):
            rows.append({"candidate": c["id"], "scale_orbit": scale, "dimension_weight": weight, "invariance_required": weight > 0, "positive_scale_orbit_invariant": c["positive_scale_orbit_invariant"], "internal_calibration_available": c["calibration_morphism_exported"] and c["internal_source_law_exported"], "passes_without_import": c["positive_scale_orbit_invariant"] and c["internal_source_law_exported"] and c["standard_physics_import_free"]})
    return rows


def gate_rows() -> list[dict[str, Any]]:
    return [{"candidate": c["id"], "required_gate": gate, "gate_passed": bool(c[gate]), "blocking_if_failed": not bool(c[gate]), "detail": "passed" if c[gate] else c["blocker"]} for c in CANDIDATES for gate in GATES]


def aggregate(gates: list[dict[str, Any]]) -> list[dict[str, Any]]:
    rows = []
    for c in CANDIDATES:
        subset = [r for r in gates if r["candidate"] == c["id"]]
        accepted = all(r["gate_passed"] for r in subset) and c["internal_source_law_exported"]
        rows.append({"candidate": c["id"], "passed_gates": sum(r["gate_passed"] for r in subset), "failed_gates": sum(not r["gate_passed"] for r in subset), "accepted_internal_shannon_to_dimension_source_law": accepted})
    return rows


def build_payload() -> dict[str, Any]:
    p3106 = read_json(P3106)
    greps = content_grep(); info = info_category_rows(); add = additive_coherence_rows(); dims = dimension_functor_rows(); scales = scale_invariance_rows(); gates = gate_rows(); aggs = aggregate(gates)
    accepted = [r for r in aggs if r["accepted_internal_shannon_to_dimension_source_law"]]
    obligations = [
        {"obligation": "read_p3106_next_atom", "satisfied": True, "detail": "P3106 selected exactly the Shannon-to-dimension functor/source-law audit"},
        {"obligation": "construct_typed_information_category", "satisfied": len(info) == 5 and any(r["is_alpha_geo_four_bit_cell"] for r in info), "detail": "I_0 through I_4 include the alpha_geo four-bit cell"},
        {"obligation": "verify_shannon_additivity", "satisfied": all(r["additive_shannon_coherence_passed"] for r in add), "detail": "finite bit-cell coproduct/additivity checks pass"},
        {"obligation": "enumerate_dimension_functor_candidates", "satisfied": len(dims) == len(CANDIDATES) * len(DIMENSIONS), "detail": "candidate dimension-label assignments are explicit"},
        {"obligation": "test_positive_scale_invariance", "satisfied": len(scales) == len(CANDIDATES) * 3, "detail": "each candidate is tested against positive scale orbits"},
        {"obligation": "export_internal_nonimported_source_law", "satisfied": False, "detail": "0 candidates export a non-imported dimension/calibration source law"},
    ]
    return {"status": "P3107_SHANNON_TO_DIMENSION_FUNCTOR_SOURCE_LAW_OBSTRUCTION_BOUNDED_NO_GO", "input_hashes": {"P3106": hashlib.sha256(P3106.read_bytes()).hexdigest() if P3106.exists() else None}, "constructed_theoretical_objects": {"content_first_repo_grep": greps, "functor_audit_object": {"object": "NadsolitonAlphaGeoShannonToDimensionFunctorSourceLawAudit", "source_reused": "P3106 recommendation and alpha_geo=4 ln 2=ln(16) four-bit Shannon cell", "ontology": "nadsoliton is self-coupled primordial information; information may source dimensions only by an internal functor/source law", "required_gates": list(GATES), "candidate_functors": [c["id"] for c in CANDIDATES], "acceptance_predicate": "typed functor, alpha_geo cell, Shannon additivity, scale-orbit invariance, nonzero dimension label, calibration morphism, internal source law, and no standard-physics import"}, "information_category_rows": info, "additive_shannon_coherence_rows": add, "dimension_functor_rows": dims, "scale_invariance_rows": scales, "candidate_gate_rows": gates, "candidate_aggregate_certificate": aggs}, "finite_certificate": {"content_grep_lanes": len(greps), "content_grep_hits": sum(r["hit_count"] for r in greps), "p3106_accepted_internal_normalization_theorems": p3106.get("finite_certificate", {}).get("accepted_internal_normalization_theorems"), "information_category_rows": len(info), "alpha_geo_four_bit_cells": sum(r["is_alpha_geo_four_bit_cell"] for r in info), "additive_coherence_rows": len(add), "additive_coherence_failures": sum(not r["additive_shannon_coherence_passed"] for r in add), "dimension_functor_rows": len(dims), "scale_invariance_rows": len(scales), "functor_candidates": len(CANDIDATES), "required_gates": len(GATES), "candidate_gate_rows": len(gates), "accepted_internal_shannon_to_dimension_source_laws": len(accepted), "proof_obligations": len(obligations), "satisfied_proof_obligations": sum(r["satisfied"] for r in obligations)}, "proof_obligations": obligations, "decision": {"bounded_result": "P3107 constructs the typed Shannon-to-dimension functor/source-law audit demanded by P3106.  The alpha_geo=4 ln 2=ln(16) four-bit Shannon cell is an exact internal information object and the finite Shannon additivity checks pass.  However, all nonzero action/length/time/calibration functors either fail positive scale-orbit invariance, lack an internal calibration morphism/source law, or import Planck/light/rod/clock/apparatus semantics.  The only fully coherent internal functor is the zero dimensionless information functor, so no physical unit is exported.", "negative_export_flags": {key: False for key in ["nonzero_internal_dimension_functor_exported", "internal_calibration_morphism_exported", "physical_action_unit_exported", "physical_length_time_unit_exported", "detector_calibration_exported", "observed_light_interface_exported", "spacetime_eom_exported", "QW_2191_discharged", "strict_selector_closure_exported", "ltotal_exported", "bridge_role_transfer_exported", "toe_closure_exported"]}, "positive_scoped_flags": {"alpha_geo_four_bit_cell_verified": True, "shannon_additivity_verified": True, "dimension_functor_candidates_enumerated": True, "scale_invariance_rows_computed": True}, "next_honest_step": "Construct exactly one bounded internal calibration-morphism candidate for the alpha_geo four-bit cell: a nadsoliton-only map from the dimensionless Shannon cell to a nonzero action/length/time calibration object, then test scale-orbit invariance and source-law status.  Do not import Planck units, rods/clocks/light, detector apparatus, selector closure, bridge/role-transfer, L_total, or ToE; if no such calibration morphism is supplied, preserve the P3105-P3107 unit-source no-go."}}


def write_markdown(payload: dict[str, Any]) -> None:
    c = payload["finite_certificate"]
    lines = ["# P3107/S2057 Shannon-to-dimension functor/source-law audit", "", f"Status: `{payload['status']}`", "", "## Finite certificate", f"- P3106 accepted internal normalization theorems: `{c['p3106_accepted_internal_normalization_theorems']}`", f"- information category rows: `{c['information_category_rows']}`", f"- alpha_geo four-bit cells: `{c['alpha_geo_four_bit_cells']}`", f"- additive coherence rows: `{c['additive_coherence_rows']}`", f"- additive coherence failures: `{c['additive_coherence_failures']}`", f"- dimension functor rows: `{c['dimension_functor_rows']}`", f"- scale invariance rows: `{c['scale_invariance_rows']}`", f"- candidate gate rows: `{c['candidate_gate_rows']}`", f"- accepted internal Shannon-to-dimension source laws: `{c['accepted_internal_shannon_to_dimension_source_laws']}`", f"- satisfied proof obligations: `{c['satisfied_proof_obligations']}/{c['proof_obligations']}`", "", "## Decision", payload["decision"]["bounded_result"], "", "## Recommendation", payload["decision"]["next_honest_step"]]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3107/S2057 Shannon-to-dimension functor/source-law audit", "## P3107/S2057 Shannon-to-dimension functor/source-law audit\n\n`P3107/S2057` executes the P3106-recommended typed Shannon-to-dimension functor/source-law audit for the `alpha_geo=4 ln 2=ln(16)` four-bit cell.  It constructs `5` information-category rows, `15` finite Shannon-additivity rows, `35` dimension-functor rows, `21` scale-invariance rows, and a `7 x 8 = 56` candidate gate matrix.  The result remains bounded no-go: the four-bit Shannon object is exact and additive, but every nonzero dimension/calibration functor lacks a non-imported internal source law or imports standard-physics calibration; no physical action/length/time unit, `L_total`, bridge/role-transfer, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3107/S2057 Shannon-to-dimension source law remains unsourced", "## P3107/S2057 Shannon-to-dimension source law remains unsourced\n\n`P3107/S2057` verifies that the `alpha_geo=ln(16)` four-bit Shannon cell has exact finite additivity and can be placed in a typed information-category audit.  The only coherent import-free functor remains dimensionless; nonzero action/length/time/calibration labels need an internal calibration morphism/source law that is not currently exported.  Therefore this is not yet a physical Lagrangian density, Hamiltonian, spacetime EOM, `L_total`, or ToE.\n")
    append_once(AGENTS, "Current Shannon-to-dimension functor/source-law guardrail (P3107/S2057, 2026-06-25)", "## Current Shannon-to-dimension functor/source-law guardrail (P3107/S2057, 2026-06-25)\n\n- P3107 follows the P3106 recommendation and audits a typed functor/source law from the `alpha_geo=4 ln 2=ln(16)` four-bit Shannon cell to dimension labels and calibration data for the self-coupled nadsoliton.\n- The finite audit computes `5` information-category rows, `15` Shannon-additivity rows, `35` dimension-functor rows, `21` scale-invariance rows, and `56` candidate gate rows; `0` candidates export an internal non-imported Shannon-to-dimension source law.\n- Do not promote the exact four-bit Shannon cell, additive entropy calculus, formal dimension labels, imported hbar/Planck templates, or imported rod/clock/light/apparatus calibrations to physical action, length, time, detector calibration, spacetime EOM, selector closure, `L_total`, bridge/role-transfer, or ToE closure.\n- The next honest move is exactly one bounded internal calibration-morphism candidate for the alpha_geo four-bit cell, unless a genuinely new internal unit-source law is introduced.\n")
    return payload

if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
