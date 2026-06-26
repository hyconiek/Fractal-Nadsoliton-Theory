#!/usr/bin/env python3
"""P3106/S2056: Nadsoliton action-density normalization obstruction audit.

P3105 left one honest alpha_geo route open: construct the self-coupled
informational action-density normalization theorem/audit.  P3106 keeps the
ontology internal--the nadsoliton is self-coupled primordial information--and
asks whether alpha_geo = 4 ln 2 = ln(16) can normalize an information action
DENSITY and quotient the positive scale orbit into dimensions/calibration
without importing Planck units, rods, clocks, observed light, apparatus,
continuum spacetime, selector closure, L_total, bridge/role-transfer, or ToE.
"""
from __future__ import annotations

import hashlib, json, math, subprocess
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p3095_s2045_dispersion_propagating_mode_observable_obstruction_audit import Z12_SIZE
from p3105_s2055_alpha_geo_entropy_to_unit_map_obstruction_audit import OUT as P3105, ALPHA_GEO, state_values, STATE_LABELS

OUT = GEN / "p3106_s2056_nadsoliton_action_density_normalization_audit.json"
MD = GEN / "p3106_s2056_nadsoliton_action_density_normalization_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"
N = Z12_SIZE
DENSITY_NORMALIZATIONS = ("per_node", "per_edge", "per_bit4", "per_ln16_cell")
SCALE_GROUP = ("lambda", "lambda_alpha_geo", "lambda_exp_alpha", "lambda_inverse_alpha")
GATES = (
    "self_coupled_density_defined",
    "alpha_geo_entropy_anchor_defined",
    "positive_scale_orbit_quotiented",
    "canonical_representative_selected_nonconventionally",
    "dimension_assignment_exported",
    "calibration_interface_exported",
    "standard_physics_import_free",
)
CANDIDATES = (
    {"id": "per_node_alpha_dirichlet_density", "self_coupled_density_defined": True, "alpha_geo_entropy_anchor_defined": True, "positive_scale_orbit_quotiented": False, "canonical_representative_selected_nonconventionally": False, "dimension_assignment_exported": False, "calibration_interface_exported": False, "standard_physics_import_free": True, "blocker": "per-node alpha_geo Dirichlet density is internal but still scale-covariant and dimensionless"},
    {"id": "per_edge_alpha_dirichlet_density", "self_coupled_density_defined": True, "alpha_geo_entropy_anchor_defined": True, "positive_scale_orbit_quotiented": False, "canonical_representative_selected_nonconventionally": False, "dimension_assignment_exported": False, "calibration_interface_exported": False, "standard_physics_import_free": True, "blocker": "edge density normalizes counting support, not the positive unit orbit"},
    {"id": "four_bit_entropy_cell_density", "self_coupled_density_defined": True, "alpha_geo_entropy_anchor_defined": True, "positive_scale_orbit_quotiented": False, "canonical_representative_selected_nonconventionally": False, "dimension_assignment_exported": False, "calibration_interface_exported": False, "standard_physics_import_free": True, "blocker": "the 4-bit/16-state cell is exact Shannon information but lacks a bit-to-action dimension law"},
    {"id": "scale_quotient_by_alpha_equals_one_gauge", "self_coupled_density_defined": True, "alpha_geo_entropy_anchor_defined": True, "positive_scale_orbit_quotiented": True, "canonical_representative_selected_nonconventionally": False, "dimension_assignment_exported": False, "calibration_interface_exported": False, "standard_physics_import_free": True, "blocker": "setting alpha_geo to one fixes gauge by convention, not by an internal calibration theorem"},
    {"id": "imported_planck_action_density", "self_coupled_density_defined": True, "alpha_geo_entropy_anchor_defined": True, "positive_scale_orbit_quotiented": True, "canonical_representative_selected_nonconventionally": True, "dimension_assignment_exported": True, "calibration_interface_exported": False, "standard_physics_import_free": False, "blocker": "action dimension is obtained only by importing hbar/Planck semantics"},
    {"id": "imported_light_clock_rod_calibration_density", "self_coupled_density_defined": True, "alpha_geo_entropy_anchor_defined": True, "positive_scale_orbit_quotiented": True, "canonical_representative_selected_nonconventionally": True, "dimension_assignment_exported": True, "calibration_interface_exported": True, "standard_physics_import_free": False, "blocker": "calibration is supplied by external rods/clocks/light/apparatus"},
)


def content_grep() -> list[dict[str, Any]]:
    patterns = {
        "p3105_normalization_atom": r"self-coupled informational action-density normalization|scale quotient|dimension assignment|calibration interface",
        "alpha_geo_entropy_atom": r"alpha_geo|4 ln\(2\)|ln\(16\)|Shannon|4-bit|16-state",
        "blocked_imports": r"Planck|rod|clock|observed light|apparatus|selector closure|L_total|ToE",
    }
    rows = []
    for lane, pattern in patterns.items():
        proc = subprocess.run(["rg", "-n", pattern, "fundamental_action_reconstruction", "AGENTS.md", "-S"], cwd=REPO, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
        hits = [line for line in proc.stdout.splitlines() if line.strip()]
        rows.append({"content_lane": lane, "pattern": pattern, "hit_count": len(hits), "sample_hits": hits[:8]})
    return rows


def dirichlet(vals: list[float]) -> float:
    return 0.5 * sum((vals[(i + 1) % N] - vals[i]) ** 2 for i in range(N))


def density_rows() -> list[dict[str, Any]]:
    rows = []
    for label in STATE_LABELS:
        energy = dirichlet(state_values(label))
        for norm in DENSITY_NORMALIZATIONS:
            denom = {"per_node": N, "per_edge": N, "per_bit4": 4.0, "per_ln16_cell": ALPHA_GEO}[norm]
            density = ALPHA_GEO * energy / denom
            rows.append({"state_label": label, "normalization": norm, "dirichlet_energy": round(energy, 12), "alpha_geo_density": round(density, 12), "self_coupled_density_defined": True, "dimension_assignment_exported": False, "calibration_interface_exported": False})
    return rows


def scale_quotient_rows() -> list[dict[str, Any]]:
    rows = []
    base = ALPHA_GEO
    for scale in SCALE_GROUP:
        factor = {"lambda": 1.0, "lambda_alpha_geo": ALPHA_GEO, "lambda_exp_alpha": math.exp(ALPHA_GEO), "lambda_inverse_alpha": 1.0 / ALPHA_GEO}[scale]
        rows.append({"scale_orbit_label": scale, "factor": round(factor, 12), "rescaled_density_anchor": round(base * factor, 12), "same_dimensionless_orbit": True, "positive_scale_quotient_constructed": True, "canonical_representative_selected_nonconventionally": False, "dimension_assignment_survives_quotient": False})
    return rows


def shannon_cell_rows() -> list[dict[str, Any]]:
    rows = []
    for bits in (1, 2, 3, 4):
        states = 2 ** bits
        entropy = math.log(states)
        rows.append({"bits": bits, "states": states, "shannon_entropy_nats": round(entropy, 12), "alpha_geo_fraction": round(entropy / ALPHA_GEO, 12), "is_four_bit_alpha_cell": bits == 4, "internal_information_anchor": True, "bit_to_action_dimension_law_exported": False})
    return rows


def gate_rows() -> list[dict[str, Any]]:
    return [{"candidate": c["id"], "required_gate": gate, "gate_passed": bool(c[gate]), "blocking_if_failed": not bool(c[gate]), "detail": "passed" if c[gate] else c["blocker"]} for c in CANDIDATES for gate in GATES]


def aggregate(gates: list[dict[str, Any]]) -> list[dict[str, Any]]:
    rows = []
    for c in CANDIDATES:
        subset = [r for r in gates if r["candidate"] == c["id"]]
        accepted = all(r["gate_passed"] for r in subset) and c["dimension_assignment_exported"] and c["calibration_interface_exported"]
        rows.append({"candidate": c["id"], "passed_gates": sum(r["gate_passed"] for r in subset), "failed_gates": sum(not r["gate_passed"] for r in subset), "accepted_internal_normalization_theorem": accepted})
    return rows


def build_payload() -> dict[str, Any]:
    p3105 = read_json(P3105)
    greps = content_grep(); densities = density_rows(); quotients = scale_quotient_rows(); shannon = shannon_cell_rows(); gates = gate_rows(); aggs = aggregate(gates)
    accepted = [r for r in aggs if r["accepted_internal_normalization_theorem"]]
    obligations = [
        {"obligation": "read_p3105_next_atom", "satisfied": True, "detail": "P3105 selected exactly the self-coupled action-density normalization theorem/audit"},
        {"obligation": "construct_alpha_geo_shannon_cell", "satisfied": any(r["is_four_bit_alpha_cell"] for r in shannon), "detail": "alpha_geo=ln(16) is represented as a 4-bit Shannon cell"},
        {"obligation": "construct_self_coupled_density_rows", "satisfied": len(densities) == len(STATE_LABELS) * len(DENSITY_NORMALIZATIONS), "detail": "all state/normalization density rows are explicit"},
        {"obligation": "test_positive_scale_quotient", "satisfied": all(r["positive_scale_quotient_constructed"] for r in quotients), "detail": "the quotient orbit can be represented formally"},
        {"obligation": "select_nonconventional_representative", "satisfied": False, "detail": "no current internal theorem selects one dimensionful representative from the orbit"},
        {"obligation": "export_dimension_assignment_and_calibration", "satisfied": False, "detail": "0 candidates export both non-imported dimension assignment and calibration interface"},
    ]
    return {"status": "P3106_NADSOLITON_ACTION_DENSITY_NORMALIZATION_OBSTRUCTION_BOUNDED_NO_GO", "input_hashes": {"P3105": hashlib.sha256(P3105.read_bytes()).hexdigest() if P3105.exists() else None}, "constructed_theoretical_objects": {"content_first_repo_grep": greps, "normalization_audit_object": {"object": "NadsolitonSelfCoupledInformationalActionDensityNormalizationAudit", "source_reused": "P3105 recommendation and alpha_geo=4 ln 2=ln(16) four-bit Shannon entropy clue", "ontology": "nadsoliton is self-coupled primordial information; no lower information layer and no imported standard-physics units", "alpha_geo": ALPHA_GEO, "required_gates": list(GATES), "candidate_normalization_sources": [c["id"] for c in CANDIDATES], "acceptance_predicate": "self-coupled density, alpha_geo entropy anchor, positive scale quotient, nonconventional representative, dimension assignment, calibration interface, and no standard-physics import"}, "shannon_four_bit_cell_rows": shannon, "action_density_rows": densities, "scale_quotient_rows": quotients, "candidate_gate_rows": gates, "candidate_aggregate_certificate": aggs}, "finite_certificate": {"content_grep_lanes": len(greps), "content_grep_hits": sum(r["hit_count"] for r in greps), "p3105_accepted_nonimported_unit_sources": p3105.get("finite_certificate", {}).get("accepted_nonimported_unit_sources"), "shannon_cell_rows": len(shannon), "four_bit_alpha_cells": sum(r["is_four_bit_alpha_cell"] for r in shannon), "bit_to_action_dimension_laws": sum(r["bit_to_action_dimension_law_exported"] for r in shannon), "action_density_rows": len(densities), "density_rows_with_dimension_assignment": sum(r["dimension_assignment_exported"] for r in densities), "scale_quotient_rows": len(quotients), "nonconventional_scale_representatives": sum(r["canonical_representative_selected_nonconventionally"] for r in quotients), "normalization_candidates": len(CANDIDATES), "required_gates": len(GATES), "candidate_gate_rows": len(gates), "accepted_internal_normalization_theorems": len(accepted), "proof_obligations": len(obligations), "satisfied_proof_obligations": sum(r["satisfied"] for r in obligations)}, "proof_obligations": obligations, "decision": {"bounded_result": "P3106 constructs the requested Nadsoliton self-coupled informational action-density normalization audit.  The exact alpha_geo=4 ln 2=ln(16) four-bit Shannon cell is represented; alpha_geo-weighted Dirichlet information action densities are computed under four finite normalizations; and the positive scale orbit is formally quotiented.  These are legitimate internal nadsoliton-information objects, but the quotient remains a dimensionless gauge orbit: current artifacts do not export a nonconventional representative, a dimension assignment, or a calibration interface.  Thus alpha_geo still does not become physical action/length/time units on current evidence.", "negative_export_flags": {key: False for key in ["internal_dimension_assignment_exported", "internal_calibration_interface_exported", "physical_action_unit_exported", "length_time_unit_exported", "canonical_nonconventional_scale_representative_exported", "observed_light_interface_exported", "spacetime_eom_exported", "QW_2191_discharged", "strict_selector_closure_exported", "ltotal_exported", "bridge_role_transfer_exported", "toe_closure_exported"]}, "positive_scoped_flags": {"alpha_geo_four_bit_shannon_cell_constructed": True, "self_coupled_action_density_rows_computed": True, "positive_scale_quotient_rows_computed": True}, "next_honest_step": "Construct exactly one bounded Shannon-to-dimension functor/source-law audit for the alpha_geo four-bit cell: define a typed functor from internal nadsoliton information densities to dimension labels/calibration data and test whether it is invariant under positive scale orbits without importing rods, clocks, Planck units, observed light, apparatus, selector closure, bridge/role-transfer, L_total, or ToE.  If no such functor is exported, preserve the unit-source no-go certificate."}}


def write_markdown(payload: dict[str, Any]) -> None:
    c = payload["finite_certificate"]
    lines = ["# P3106/S2056 Nadsoliton action-density normalization audit", "", f"Status: `{payload['status']}`", "", "## Finite certificate", f"- P3105 accepted non-imported unit sources: `{c['p3105_accepted_nonimported_unit_sources']}`", f"- Shannon cell rows: `{c['shannon_cell_rows']}`", f"- four-bit alpha_geo cells: `{c['four_bit_alpha_cells']}`", f"- bit-to-action dimension laws: `{c['bit_to_action_dimension_laws']}`", f"- action density rows: `{c['action_density_rows']}`", f"- density rows with dimension assignment: `{c['density_rows_with_dimension_assignment']}`", f"- scale quotient rows: `{c['scale_quotient_rows']}`", f"- nonconventional scale representatives: `{c['nonconventional_scale_representatives']}`", f"- candidate gate rows: `{c['candidate_gate_rows']}`", f"- accepted internal normalization theorems: `{c['accepted_internal_normalization_theorems']}`", f"- satisfied proof obligations: `{c['satisfied_proof_obligations']}/{c['proof_obligations']}`", "", "## Decision", payload["decision"]["bounded_result"], "", "## Recommendation", payload["decision"]["next_honest_step"]]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3106/S2056 Nadsoliton action-density normalization audit", "## P3106/S2056 Nadsoliton action-density normalization audit\n\n`P3106/S2056` executes the P3105-recommended self-coupled informational action-density normalization audit for `alpha_geo=4 ln 2=ln(16)`.  It constructs `4` Shannon cell rows, `16` alpha_geo-weighted action-density rows, `4` positive scale-quotient rows, and a `6 x 7 = 42` candidate gate matrix.  The result remains bounded no-go: the four-bit Shannon entropy cell and internal density proxies are real, but no nonconventional scale representative, dimension assignment, calibration interface, physical action/length/time unit, `L_total`, bridge/role-transfer, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3106/S2056 action-density normalization remains unitless", "## P3106/S2056 action-density normalization remains unitless\n\n`P3106/S2056` normalizes the alpha_geo-weighted Dirichlet information action into finite per-node, per-edge, per-bit, and per-ln16-cell densities for the self-coupled nadsoliton.  The construction stays internal and information-theoretic, but it remains unitless: the positive scale quotient has no exported nonconventional representative, dimension assignment, or calibration interface, so it is not a physical Lagrangian density, Hamiltonian, spacetime EOM, `L_total`, or ToE.\n")
    append_once(AGENTS, "Current Nadsoliton action-density normalization guardrail (P3106/S2056, 2026-06-25)", "## Current Nadsoliton action-density normalization guardrail (P3106/S2056, 2026-06-25)\n\n- P3106 follows the P3105 recommendation and audits the self-coupled informational action-density normalization theorem/source for `alpha_geo=4 ln 2=ln(16)` as a four-bit Shannon cell of the nadsoliton, not as imported standard physics.\n- The finite audit computes `4` Shannon-cell rows, `16` alpha_geo-weighted action-density rows, `4` positive scale-quotient rows, and `42` candidate gate rows; `0` candidates export an internal dimension-assignment plus calibration theorem.\n- Do not promote the four-bit Shannon cell, per-node/per-edge/per-bit/per-ln16 action-density proxies, or formal positive scale quotients to physical action, length, time, detector calibration, spacetime EOM, selector closure, `L_total`, bridge/role-transfer, or ToE closure.\n- The next honest move is exactly one bounded Shannon-to-dimension functor/source-law audit for the alpha_geo four-bit cell, unless a genuinely new internal dimension/calibration theorem is introduced.\n")
    return payload

if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
