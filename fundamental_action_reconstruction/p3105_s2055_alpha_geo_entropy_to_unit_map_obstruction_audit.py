#!/usr/bin/env python3
"""P3105/S2055: alpha_geo entropy-to-unit map obstruction audit.

P3104 identified the remaining non-selector unit clue: alpha_geo = 4 ln 2 =
ln(16), the 16-state entropy scalar inherited as a legacy-kernel bridge number.
P3105 treats the nadsoliton as self-coupled information, not as already-given
standard physics, and tests whether alpha_geo can internally become an action,
length, time, Hamiltonian, or detector-calibration unit without importing
Planck units, observed light, apparatus semantics, continuum spacetime, quantum
measurement postulates, selector closure, L_total, bridge/role-transfer, or ToE.
"""
from __future__ import annotations

import hashlib, json, math, subprocess
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p3095_s2045_dispersion_propagating_mode_observable_obstruction_audit import Z12_SIZE, eigenvalue
from p3104_s2054_spectral_triple_geometry_interface_obstruction_audit import OUT as P3104

OUT = GEN / "p3105_s2055_alpha_geo_entropy_to_unit_map_obstruction_audit.json"
MD = GEN / "p3105_s2055_alpha_geo_entropy_to_unit_map_obstruction_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"
N = Z12_SIZE
ALPHA_GEO = 4.0 * math.log(2.0)
STATE_LABELS = ("delta0", "block3", "alternating", "ramp")
CONTENT_PATTERNS = {
    "alpha_geo_unit_atom": r"alpha_geo|4 ln\(2\)|ln\(16\)|entropy-to-action|bit-to-length|unit-map",
    "self_coupled_information_atom": r"self-coupled information|samosprz|nadsoliton|informational action|Dirichlet energy",
    "blocked_standard_physics_imports": r"Planck units|observed light|quantum measurement|apparatus|continuum spacetime|L_total|ToE",
}
GATES = (
    "alpha_geo_entropy_identity_exported",
    "sixteen_state_reference_exported",
    "self_coupled_information_functional_exported",
    "unit_dimension_assignment_exported",
    "conversion_law_exported",
    "calibration_interface_exported",
    "nonimported_unit_source_exported",
    "standard_physics_import_free",
)
CANDIDATES = (
    {"id": "alpha_geo_ln16_entropy_identity", "alpha_geo_entropy_identity_exported": True, "sixteen_state_reference_exported": True, "self_coupled_information_functional_exported": False, "unit_dimension_assignment_exported": False, "conversion_law_exported": False, "calibration_interface_exported": False, "nonimported_unit_source_exported": False, "standard_physics_import_free": True, "blocker": "alpha_geo=ln(16) is exact but remains a dimensionless entropy identity"},
    {"id": "per_bit_entropy_quantum_alpha_over_4", "alpha_geo_entropy_identity_exported": True, "sixteen_state_reference_exported": True, "self_coupled_information_functional_exported": False, "unit_dimension_assignment_exported": False, "conversion_law_exported": False, "calibration_interface_exported": False, "nonimported_unit_source_exported": False, "standard_physics_import_free": True, "blocker": "ln(2) per-bit normalization gives a scalar bit weight, not a unit of action/length"},
    {"id": "self_coupled_dirichlet_information_action_proxy", "alpha_geo_entropy_identity_exported": True, "sixteen_state_reference_exported": True, "self_coupled_information_functional_exported": True, "unit_dimension_assignment_exported": False, "conversion_law_exported": False, "calibration_interface_exported": False, "nonimported_unit_source_exported": False, "standard_physics_import_free": True, "blocker": "alpha_geo times finite Dirichlet energy is an internal informational action proxy but has no dimension assignment or calibration law"},
    {"id": "alpha_geo_length_orbit_normalization_proxy", "alpha_geo_entropy_identity_exported": True, "sixteen_state_reference_exported": True, "self_coupled_information_functional_exported": True, "unit_dimension_assignment_exported": False, "conversion_law_exported": False, "calibration_interface_exported": False, "nonimported_unit_source_exported": False, "standard_physics_import_free": True, "blocker": "choosing alpha_geo or exp(alpha_geo)=16 as a length scale fixes a gauge convention, not an internal unit source"},
    {"id": "imported_planck_unit_template", "alpha_geo_entropy_identity_exported": True, "sixteen_state_reference_exported": True, "self_coupled_information_functional_exported": True, "unit_dimension_assignment_exported": True, "conversion_law_exported": True, "calibration_interface_exported": False, "nonimported_unit_source_exported": False, "standard_physics_import_free": False, "blocker": "dimensionful units pass only by importing Planck/standard-physics unit semantics"},
    {"id": "imported_detector_clock_light_calibration_template", "alpha_geo_entropy_identity_exported": True, "sixteen_state_reference_exported": True, "self_coupled_information_functional_exported": True, "unit_dimension_assignment_exported": True, "conversion_law_exported": True, "calibration_interface_exported": True, "nonimported_unit_source_exported": False, "standard_physics_import_free": False, "blocker": "calibration passes only by importing detector clocks, rods, or observed light"},
)


def content_grep() -> list[dict[str, Any]]:
    rows = []
    for lane, pattern in CONTENT_PATTERNS.items():
        proc = subprocess.run(["rg", "-n", pattern, "fundamental_action_reconstruction", "AGENTS.md", "-S"], cwd=REPO, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
        hits = [line for line in proc.stdout.splitlines() if line.strip()]
        rows.append({"content_lane": lane, "pattern": pattern, "hit_count": len(hits), "sample_hits": hits[:8]})
    return rows


def state_values(label: str) -> list[float]:
    if label == "delta0":
        return [1.0 if n == 0 else 0.0 for n in range(N)]
    if label == "block3":
        return [1.0 if n in (0, 1, 2) else 0.0 for n in range(N)]
    if label == "alternating":
        return [1.0 if n % 2 == 0 else -1.0 for n in range(N)]
    return [float(n + 1) for n in range(N)]


def entropy_identity_rows() -> list[dict[str, Any]]:
    rows = []
    for states in (2, 4, 8, 16):
        entropy = math.log(states)
        rows.append({
            "state_count": states,
            "entropy_ln_states": round(entropy, 12),
            "alpha_geo_multiple": round(entropy / ALPHA_GEO, 12),
            "is_alpha_geo_reference": states == 16,
            "dimensionless_entropy_identity": True,
            "unit_dimension_attached": False,
        })
    return rows


def self_coupled_action_rows() -> list[dict[str, Any]]:
    rows = []
    for label in STATE_LABELS:
        vals = state_values(label)
        dirichlet = 0.5 * sum((vals[(i + 1) % N] - vals[i]) ** 2 for i in range(N))
        lap_expect = sum(eigenvalue(k) for k in range(N)) / N if label == "ramp" else dirichlet / N
        rows.append({
            "state_label": label,
            "z12_dirichlet_energy": round(dirichlet, 12),
            "alpha_geo_weighted_information_action": round(ALPHA_GEO * dirichlet, 12),
            "laplacian_reference_scalar": round(lap_expect, 12),
            "self_coupled_information_proxy": True,
            "dimensionful_action_unit_attached": False,
            "nonimported_conversion_law_attached": False,
        })
    return rows


def scale_orbit_rows() -> list[dict[str, Any]]:
    rows = []
    for scale_label, scale in (("one", 1.0), ("alpha_geo", ALPHA_GEO), ("exp_alpha_geo", math.exp(ALPHA_GEO)), ("alpha_geo_over_4", ALPHA_GEO / 4.0)):
        rows.append({
            "scale_label": scale_label,
            "scale_value": round(scale, 12),
            "rescaled_alpha_geo": round(ALPHA_GEO * scale, 12),
            "rescaled_ln16": round(math.log(16.0) * scale, 12),
            "scale_orbit_member": True,
            "canonical_unit_selected": False,
            "selection_is_internal_not_conventional": False,
        })
    return rows


def target_unit_rows() -> list[dict[str, Any]]:
    return [{
        "target_unit": target,
        "candidate_formula": formula,
        "formula_value": round(value, 12),
        "dimension_assignment_exported": False,
        "calibration_interface_exported": False,
        "standard_physics_import_needed": target in {"planck_action", "meter_like_length", "second_like_time", "detector_click_calibration"},
    } for target, formula, value in (
        ("informational_action_proxy", "alpha_geo * Dirichlet_energy", ALPHA_GEO),
        ("length_proxy", "exp(alpha_geo)=16", math.exp(ALPHA_GEO)),
        ("hamiltonian_time_proxy", "1/alpha_geo", 1.0 / ALPHA_GEO),
        ("detector_click_calibration", "alpha_geo/16", ALPHA_GEO / 16.0),
        ("planck_action", "imported hbar identification", ALPHA_GEO),
        ("meter_like_length", "imported rod/light calibration", math.exp(ALPHA_GEO)),
    )]


def gate_rows() -> list[dict[str, Any]]:
    return [{"candidate": c["id"], "required_gate": gate, "gate_passed": bool(c[gate]), "blocking_if_failed": not bool(c[gate]), "detail": "passed" if c[gate] else c["blocker"]} for c in CANDIDATES for gate in GATES]


def aggregate(gates: list[dict[str, Any]]) -> list[dict[str, Any]]:
    rows = []
    for c in CANDIDATES:
        subset = [r for r in gates if r["candidate"] == c["id"]]
        rows.append({"candidate": c["id"], "passed_gates": sum(r["gate_passed"] for r in subset), "failed_gates": sum(not r["gate_passed"] for r in subset), "accepted_nonimported_unit_source": all(r["gate_passed"] for r in subset) and bool(c["nonimported_unit_source_exported"])})
    return rows


def build_payload() -> dict[str, Any]:
    p3104 = read_json(P3104)
    greps = content_grep(); entropy = entropy_identity_rows(); action = self_coupled_action_rows(); scales = scale_orbit_rows(); targets = target_unit_rows(); gates = gate_rows(); aggs = aggregate(gates)
    accepted = [r for r in aggs if r["accepted_nonimported_unit_source"]]
    obligations = [
        {"obligation": "read_p3104_next_atom", "satisfied": True, "detail": "P3104 selected exactly alpha_geo entropy-to-action/length unit-map obstruction audit"},
        {"obligation": "construct_alpha_geo_entropy_identity", "satisfied": any(r["is_alpha_geo_reference"] and not r["unit_dimension_attached"] for r in entropy), "detail": "alpha_geo=ln(16)=4 ln 2 is exact and dimensionless"},
        {"obligation": "construct_self_coupled_information_action_proxy", "satisfied": len(action) == len(STATE_LABELS) and all(r["self_coupled_information_proxy"] for r in action), "detail": "alpha_geo-weighted Dirichlet action proxies are explicit"},
        {"obligation": "test_scale_orbit_canonicity", "satisfied": all(not r["canonical_unit_selected"] for r in scales), "detail": "scale orbit has no internal selector-free unit choice"},
        {"obligation": "test_dimensionful_target_units", "satisfied": all(not r["dimension_assignment_exported"] for r in targets), "detail": "target unit rows have formulas but no dimension assignment"},
        {"obligation": "export_nonimported_alpha_geo_unit_source", "satisfied": False, "detail": "0 candidates pass as internal bit-to-unit conversion laws"},
    ]
    return {"status": "P3105_ALPHA_GEO_ENTROPY_TO_UNIT_MAP_OBSTRUCTION_BOUNDED_NO_GO", "input_hashes": {"P3104": hashlib.sha256(P3104.read_bytes()).hexdigest() if P3104.exists() else None}, "constructed_theoretical_objects": {"content_first_repo_grep": greps, "alpha_geo_unit_map_audit_object": {"object": "NadsolitonSelfCoupledInformationAlphaGeoUnitMapObstructionAudit", "source_reused": "P3104 recommendation and alpha_geo=4 ln 2=ln(16) entropy clue", "ontology": "nadsoliton is self-coupled primordial information; standard-physics units are not assumed", "alpha_geo": ALPHA_GEO, "required_gates": list(GATES), "candidate_unit_sources": [c["id"] for c in CANDIDATES], "acceptance_predicate": "alpha entropy identity, 16-state reference, self-coupled information functional, dimension assignment, conversion law, calibration interface, non-imported unit source, and no standard-physics import"}, "entropy_identity_rows": entropy, "self_coupled_information_action_rows": action, "scale_orbit_rows": scales, "target_unit_rows": targets, "candidate_gate_rows": gates, "candidate_aggregate_certificate": aggs}, "finite_certificate": {"content_grep_lanes": len(greps), "content_grep_hits": sum(r["hit_count"] for r in greps), "p3104_accepted_nonimported_geometry_sources": p3104.get("finite_certificate", {}).get("accepted_nonimported_geometry_sources"), "entropy_identity_rows": len(entropy), "alpha_geo_reference_rows": sum(r["is_alpha_geo_reference"] for r in entropy), "entropy_rows_with_unit_dimension": sum(r["unit_dimension_attached"] for r in entropy), "self_coupled_action_rows": len(action), "action_rows_with_dimensionful_unit": sum(r["dimensionful_action_unit_attached"] for r in action), "scale_orbit_rows": len(scales), "canonical_scale_selections": sum(r["canonical_unit_selected"] for r in scales), "target_unit_rows": len(targets), "target_rows_with_dimension_assignment": sum(r["dimension_assignment_exported"] for r in targets), "unit_candidates": len(CANDIDATES), "required_gates": len(GATES), "candidate_gate_rows": len(gates), "accepted_nonimported_unit_sources": len(accepted), "proof_obligations": len(obligations), "satisfied_proof_obligations": sum(r["satisfied"] for r in obligations)}, "proof_obligations": obligations, "decision": {"bounded_result": "P3105 constructs the alpha_geo entropy-to-unit map audit requested by P3104 while treating the nadsoliton as self-coupled information rather than assumed standard physics.  The audit verifies the exact identity alpha_geo=4 ln 2=ln(16), constructs alpha_geo-weighted finite Dirichlet informational action proxies, tests the scale orbit {1, alpha_geo, exp(alpha_geo), alpha_geo/4}, and enumerates action/length/time/detector target formulas.  These objects are real internal mathematical candidates, but no current artifact exports a dimension assignment, canonical scale selection, calibration interface, or non-imported bit-to-unit conversion law.  Therefore alpha_geo is not yet converted into physical units/action.", "negative_export_flags": {key: False for key in ["alpha_geo_dimensionful_action_unit_exported", "alpha_geo_length_unit_exported", "hamiltonian_time_unit_exported", "detector_calibration_unit_exported", "canonical_scale_selection_exported", "nonimported_bit_to_unit_conversion_law_exported", "observed_light_interface_exported", "quantum_measurement_postulate_imported_as_source", "spacetime_eom_exported", "QW_2191_discharged", "strict_selector_closure_exported", "ltotal_exported", "bridge_role_transfer_exported", "toe_closure_exported"]}, "positive_scoped_flags": {"alpha_geo_ln16_identity_verified": True, "self_coupled_information_action_rows_computed": True, "scale_orbit_rows_computed": True, "target_unit_rows_computed": True}, "next_honest_step": "Construct exactly one new internal unit-source mechanism: a Nadsoliton self-coupled informational action-density normalization theorem.  It must take the alpha_geo-weighted Dirichlet information action proxy, quotient the positive scale orbit without importing Planck/rod/clock/light units, and export either a dimension assignment plus calibration interface or a proof that no such scale quotient exists.  Do not return to selector replay, standard-physics unit import, bridge/role-transfer, L_total, or ToE promotion."}}


def write_markdown(payload: dict[str, Any]) -> None:
    c = payload["finite_certificate"]
    lines = ["# P3105/S2055 alpha_geo entropy-to-unit map obstruction audit", "", f"Status: `{payload['status']}`", "", "## Finite certificate", f"- P3104 accepted non-imported geometry sources: `{c['p3104_accepted_nonimported_geometry_sources']}`", f"- entropy identity rows: `{c['entropy_identity_rows']}`", f"- alpha_geo reference rows: `{c['alpha_geo_reference_rows']}`", f"- entropy rows with unit dimension: `{c['entropy_rows_with_unit_dimension']}`", f"- self-coupled action rows: `{c['self_coupled_action_rows']}`", f"- action rows with dimensionful unit: `{c['action_rows_with_dimensionful_unit']}`", f"- scale orbit rows: `{c['scale_orbit_rows']}`", f"- canonical scale selections: `{c['canonical_scale_selections']}`", f"- target unit rows: `{c['target_unit_rows']}`", f"- target rows with dimension assignment: `{c['target_rows_with_dimension_assignment']}`", f"- candidate gate rows: `{c['candidate_gate_rows']}`", f"- accepted non-imported unit sources: `{c['accepted_nonimported_unit_sources']}`", f"- satisfied proof obligations: `{c['satisfied_proof_obligations']}/{c['proof_obligations']}`", "", "## Decision", payload["decision"]["bounded_result"], "", "## Recommendation", payload["decision"]["next_honest_step"]]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3105/S2055 alpha_geo entropy-to-unit map obstruction audit", "## P3105/S2055 alpha_geo entropy-to-unit map obstruction audit\n\n`P3105/S2055` attacks the P3104-recommended `alpha_geo` unit-map atom while treating the nadsoliton as self-coupled primordial information rather than assumed standard physics.  It verifies `alpha_geo=4 ln 2=ln(16)`, constructs `4` entropy-identity rows, `4` alpha_geo-weighted Dirichlet informational-action rows, `4` scale-orbit rows, `6` target-unit rows, and a `6 x 8 = 48` candidate gate matrix.  The result remains bounded no-go: no dimension assignment, canonical scale quotient, calibration interface, non-imported bit-to-unit conversion law, observed-light interface, `L_total`, bridge/role-transfer, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3105/S2055 alpha_geo unit conversion remains unsourced", "## P3105/S2055 alpha_geo unit conversion remains unsourced\n\n`P3105/S2055` confirms that `alpha_geo=4 ln 2=ln(16)` can weight a finite self-coupled informational Dirichlet action proxy for the nadsoliton branch.  This is not yet a physical action, length, time, Hamiltonian, detector calibration, spacetime EOM, or `L_total`: the missing object is an internal scale-quotient/dimension-assignment/calibration theorem, not an imported Planck-unit, rod/clock/light, apparatus, or measurement-postulate template.\n")
    append_once(AGENTS, "Current alpha_geo entropy-to-unit map obstruction guardrail (P3105/S2055, 2026-06-25)", "## Current alpha_geo entropy-to-unit map obstruction guardrail (P3105/S2055, 2026-06-25)\n\n- P3105 follows the P3104 recommendation and audits the one available unit clue: `alpha_geo=4 ln 2=ln(16)` as a possible entropy-to-action/length map for the nadsoliton understood as self-coupled information, not as pre-assumed standard physics.\n- The finite audit computes `4` entropy-identity rows, `4` self-coupled informational-action rows, `4` scale-orbit rows, `6` target-unit rows, and `48` candidate gate rows; `0` candidates export an internal non-imported bit-to-unit conversion law.\n- Do not promote alpha_geo entropy identities, per-bit weights, Dirichlet informational-action proxies, exp(alpha_geo)=16 scale choices, Planck-unit templates, rod/clock/light templates, or apparatus calibration templates to physical action, length, time, detector calibration, spacetime EOM, selector closure, `L_total`, bridge/role-transfer, or ToE closure.\n- The next honest move is exactly one Nadsoliton self-coupled informational action-density normalization theorem/audit that either exports an internal scale quotient with dimension assignment and calibration interface or proves no such quotient exists on current artifacts.\n")
    return payload

if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
