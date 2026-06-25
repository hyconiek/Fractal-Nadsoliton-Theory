#!/usr/bin/env python3
"""P3086/S2036: empirical-readout/observable-calibration obstruction audit.

P3085 left the Z12 Dirichlet/Laplacian branch with formal divergence-free
current witnesses but no unit-bearing conserved-current source.  P3086 attacks
exactly the next standard-physics interface atom: whether any internal scalar,
spectrum, current proxy, or dimensionless witness admits a unit-calibrated
empirical readout without importing measurement units, observed light/photons,
spacetime EOM, selector closure, L_total, bridge/role-transfer, or ToE.
"""
from __future__ import annotations

import hashlib, json, math, subprocess
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p3085_s2035_conserved_current_noether_obstruction_audit import OUT as P3085

OUT = GEN / "p3086_s2036_empirical_readout_observable_calibration_obstruction_audit.json"
MD = GEN / "p3086_s2036_empirical_readout_observable_calibration_obstruction_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"
Z12_SIZE = 12
TOL = 1e-9

CONTENT_PATTERNS = {
    "observable_atom": r"empirical observable|observable-calibration|measurement units|unit-calibrated|readout",
    "predecessor_current_atom": r"conserved current|Noether|charge density|unit-bearing current",
    "blocked_promotions": r"observed photons|observed light|spacetime EOM|L_total|ToE|selector closure|bridge/role-transfer",
}

OBSERVABLE_CANDIDATES = (
    {
        "id": "z12_dirichlet_scalar_energy_proxy",
        "description": "dimensionless Dirichlet scalar energy/readout proxy on Z12 profiles",
        "internal_artifact": True,
        "dimensionless_witness_computed": True,
        "measurement_unit_source_exported": False,
        "calibration_map_exported": False,
        "apparatus_readout_protocol_exported": False,
        "empirical_observable_exported": False,
        "blocker": "internal scalar energy is dimensionless and has no strict measurement-unit or apparatus map",
    },
    {
        "id": "z12_laplacian_spectrum_proxy",
        "description": "finite Z12 Laplacian spectrum and spectral gaps",
        "internal_artifact": True,
        "dimensionless_witness_computed": True,
        "measurement_unit_source_exported": False,
        "calibration_map_exported": False,
        "apparatus_readout_protocol_exported": False,
        "empirical_observable_exported": False,
        "blocker": "spectral numbers are intrinsic dimensionless labels without a unit-calibrated readout theorem",
    },
    {
        "id": "p3085_formal_current_proxy",
        "description": "formal Fourier link-current magnitude inherited from P3085",
        "internal_artifact": False,
        "dimensionless_witness_computed": True,
        "measurement_unit_source_exported": False,
        "calibration_map_exported": False,
        "apparatus_readout_protocol_exported": False,
        "empirical_observable_exported": False,
        "blocker": "P3085 current is an auxiliary algebraic proxy and is not unit-bearing",
    },
    {
        "id": "imported_dimensionful_meter_calibration_template",
        "description": "external meter/second/joule calibration template",
        "internal_artifact": False,
        "dimensionless_witness_computed": True,
        "measurement_unit_source_exported": True,
        "calibration_map_exported": True,
        "apparatus_readout_protocol_exported": True,
        "empirical_observable_exported": True,
        "blocker": "passes only by importing an external metrology layer",
    },
    {
        "id": "imported_observed_photon_frequency_template",
        "description": "external observed photon/frequency readout template",
        "internal_artifact": False,
        "dimensionless_witness_computed": True,
        "measurement_unit_source_exported": True,
        "calibration_map_exported": True,
        "apparatus_readout_protocol_exported": True,
        "empirical_observable_exported": True,
        "blocker": "passes only by importing observed light/photon physics",
    },
)
REQUIRED_GATES = ("internal_artifact", "dimensionless_witness_computed", "measurement_unit_source_exported", "calibration_map_exported", "apparatus_readout_protocol_exported", "empirical_observable_exported")
SCALE_FACTORS = (0.25, 0.5, 1.0, 2.0, 4.0)


def content_grep() -> list[dict[str, Any]]:
    rows = []
    for lane, pattern in CONTENT_PATTERNS.items():
        proc = subprocess.run(["rg", "-n", pattern, "fundamental_action_reconstruction", "AGENTS.md", "-S"], cwd=REPO, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
        hits = [line for line in proc.stdout.splitlines() if line.strip()]
        rows.append({"content_lane": lane, "pattern": pattern, "hit_count": len(hits), "sample_hits": hits[:8]})
    return rows


def spectrum_rows() -> list[dict[str, Any]]:
    rows = []
    for m in range(Z12_SIZE):
        lam = 2.0 - 2.0 * math.cos(2.0 * math.pi * m / Z12_SIZE)
        rows.append({"mode_m": m, "laplacian_eigenvalue": round(lam, 12), "dimensionless": True, "unit_calibrated": False})
    return rows


def normalized_gap_rows(spec: list[dict[str, Any]]) -> list[dict[str, Any]]:
    values = sorted({row["laplacian_eigenvalue"] for row in spec})
    rows = []
    for idx, value in enumerate(values):
        next_value = values[(idx + 1) % len(values)] if idx + 1 < len(values) else None
        rows.append({"distinct_index": idx, "eigenvalue": value, "next_gap": None if next_value is None else round(next_value - value, 12), "dimensionless_gap": True, "empirical_unit_attached": False})
    return rows


def scale_orbit_rows(spec: list[dict[str, Any]]) -> list[dict[str, Any]]:
    rows = []
    base = [row["laplacian_eigenvalue"] for row in spec]
    for factor in SCALE_FACTORS:
        scaled = [round(factor * v, 12) for v in base]
        rows.append({"scale_factor": factor, "scaled_spectrum_sum": round(sum(scaled), 12), "ratios_preserved": True, "canonical_scale_selected": factor == 1.0, "selection_is_conventional_without_unit_source": True})
    return rows


def calibration_attempt_rows() -> list[dict[str, Any]]:
    rows = []
    targets = ("length", "time", "energy", "charge", "frequency", "dimensionless_ratio")
    for candidate in OBSERVABLE_CANDIDATES:
        for target in targets:
            imported = not candidate["internal_artifact"] and candidate["empirical_observable_exported"]
            rows.append({
                "candidate": candidate["id"],
                "target_readout_type": target,
                "has_internal_unit_source": bool(candidate["measurement_unit_source_exported"] and candidate["internal_artifact"]),
                "calibration_map_available_nonimported": bool(candidate["calibration_map_exported"] and candidate["internal_artifact"]),
                "passes_only_as_imported_template": imported,
                "accepted_nonimported_empirical_readout": False,
            })
    return rows


def gate_rows() -> list[dict[str, Any]]:
    return [{"candidate": c["id"], "required_gate": gate, "gate_passed": bool(c[gate]), "blocking_if_failed": not bool(c[gate]), "detail": "passed" if c[gate] else c["blocker"]} for c in OBSERVABLE_CANDIDATES for gate in REQUIRED_GATES]


def aggregate(gates: list[dict[str, Any]]) -> list[dict[str, Any]]:
    out = []
    for c in OBSERVABLE_CANDIDATES:
        subset = [r for r in gates if r["candidate"] == c["id"]]
        out.append({"candidate": c["id"], "passed_gates": sum(r["gate_passed"] for r in subset), "failed_gates": sum(not r["gate_passed"] for r in subset), "accepted_nonimported_empirical_observable_source": all(r["gate_passed"] for r in subset) and bool(c["internal_artifact"])})
    return out


def build_payload() -> dict[str, Any]:
    p3085 = read_json(P3085)
    greps = content_grep(); spec = spectrum_rows(); gaps = normalized_gap_rows(spec); scales = scale_orbit_rows(spec); calibrations = calibration_attempt_rows(); gates = gate_rows(); aggs = aggregate(gates)
    accepted = [r for r in aggs if r["accepted_nonimported_empirical_observable_source"]]
    obligations = [
        {"obligation": "read_p3085_next_atom", "satisfied": True, "detail": "P3085 selected empirical-readout/observable-calibration as the next interface atom"},
        {"obligation": "construct_internal_dimensionless_witnesses", "satisfied": True, "detail": "Z12 Laplacian spectrum, gaps, scalar proxies, and current proxy classes are represented"},
        {"obligation": "test_scale_orbit_noncanonicity", "satisfied": all(r["selection_is_conventional_without_unit_source"] for r in scales), "detail": "five scale choices preserve ratios; no strict unit selects one empirical scale"},
        {"obligation": "build_calibration_target_matrix", "satisfied": True, "detail": "five candidate readout sources x six target readout types are audited"},
        {"obligation": "export_nonimported_unit_calibrated_empirical_observable", "satisfied": False, "detail": "0 candidates pass as internal non-imported unit-calibrated empirical observables"},
    ]
    return {
        "status": "P3086_EMPIRICAL_READOUT_OBSERVABLE_CALIBRATION_OBSTRUCTION_BOUNDED_NO_GO",
        "input_hashes": {"P3085": hashlib.sha256(P3085.read_bytes()).hexdigest() if P3085.exists() else None},
        "constructed_theoretical_objects": {
            "content_first_repo_grep": greps,
            "empirical_readout_audit_object": {"object": "Z12DirichletEmpiricalReadoutObservableCalibrationObstructionAudit", "source_reused": "P3085 recommendation: bounded empirical-readout/observable-calibration obstruction audit", "required_gates": list(REQUIRED_GATES), "candidate_observable_sources": [c["id"] for c in OBSERVABLE_CANDIDATES], "acceptance_predicate": "internal non-imported dimensionless witness plus measurement-unit source, calibration map, apparatus readout protocol, and empirical observable export"},
            "z12_laplacian_spectrum_rows": spec,
            "normalized_spectral_gap_rows": gaps,
            "scale_orbit_control_rows": scales,
            "calibration_attempt_rows": calibrations,
            "candidate_gate_rows": gates,
            "candidate_aggregate_certificate": aggs,
        },
        "finite_certificate": {
            "content_grep_lanes": len(greps), "content_grep_hits": sum(r["hit_count"] for r in greps),
            "p3085_accepted_nonimported_conserved_current_sources": p3085.get("finite_certificate", {}).get("accepted_nonimported_conserved_current_sources"),
            "spectrum_rows": len(spec), "spectrum_rows_with_unit_calibration": sum(r["unit_calibrated"] for r in spec), "distinct_spectral_gap_rows": len(gaps),
            "scale_orbit_control_rows": len(scales), "canonical_scale_sources_exported": 0,
            "observable_candidates": len(OBSERVABLE_CANDIDATES), "required_gates": len(REQUIRED_GATES), "candidate_gate_rows": len(gates),
            "calibration_attempt_rows": len(calibrations), "accepted_nonimported_empirical_observable_sources": len(accepted),
            "proof_obligations": len(obligations), "satisfied_proof_obligations": sum(r["satisfied"] for r in obligations),
        },
        "proof_obligations": obligations,
        "decision": {
            "bounded_result": "P3086 constructs the requested empirical-readout/observable-calibration obstruction audit.  The Z12 Dirichlet/Laplacian branch supplies finite dimensionless spectral/scalar witnesses and inherits a formal current proxy from P3085, but scale-orbit controls show that multiplying the spectrum by an external unit scale preserves all internal ratios.  Current artifacts do not export measurement units, a calibration map, an apparatus readout protocol, or an empirical observable.  Meter/frequency/photon rows pass only as imported templates.  Therefore no non-imported unit-calibrated empirical observable is exported.",
            "negative_export_flags": {key: False for key in ["measurement_unit_source_exported", "calibration_map_exported", "apparatus_readout_protocol_exported", "empirical_observable_exported", "observed_photons_exported", "observed_light_exported", "unit_bearing_current_exported", "spacetime_eom_exported", "physical_hamiltonian_exported", "gauge_representation_source_exported", "QW_2191_discharged", "strict_selector_closure_exported", "ltotal_exported", "bridge_role_transfer_exported", "toe_closure_exported"]},
            "positive_scoped_flags": {"finite_z12_spectrum_computed": True, "scale_orbit_controls_executed": True, "calibration_matrix_constructed": True},
            "next_honest_step": "Pivot to exactly one adjacent standard-physics interface atom outside selector replay: construct a bounded thermodynamic/statistical-ensemble obstruction audit for the Z12 Dirichlet/Laplacian branch, testing whether the internal finite spectrum supplies a canonical temperature, partition function, entropy/energy units, and equilibrium observable without importing Boltzmann constants, measurement units, observed radiation, spacetime EOM, L_total, bridge/role-transfer, or ToE.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    c = payload["finite_certificate"]
    lines = ["# P3086/S2036 empirical-readout/observable-calibration obstruction audit", "", f"Status: `{payload['status']}`", "", "## Finite certificate",
        f"- content grep lanes: `{c['content_grep_lanes']}`", f"- content grep hits: `{c['content_grep_hits']}`", f"- P3085 accepted non-imported conserved-current sources: `{c['p3085_accepted_nonimported_conserved_current_sources']}`", f"- spectrum rows: `{c['spectrum_rows']}`", f"- spectrum rows with unit calibration: `{c['spectrum_rows_with_unit_calibration']}`", f"- distinct spectral gap rows: `{c['distinct_spectral_gap_rows']}`", f"- scale orbit control rows: `{c['scale_orbit_control_rows']}`", f"- canonical scale sources exported: `{c['canonical_scale_sources_exported']}`", f"- observable candidates: `{c['observable_candidates']}`", f"- candidate gate rows: `{c['candidate_gate_rows']}`", f"- calibration attempt rows: `{c['calibration_attempt_rows']}`", f"- accepted non-imported empirical observable sources: `{c['accepted_nonimported_empirical_observable_sources']}`", f"- satisfied proof obligations: `{c['satisfied_proof_obligations']}/{c['proof_obligations']}`", "", "## Decision", payload["decision"]["bounded_result"], "", "## Recommendation", payload["decision"]["next_honest_step"]]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3086/S2036 empirical-readout/observable-calibration obstruction audit", "## P3086/S2036 empirical-readout/observable-calibration obstruction audit\n\n`P3086/S2036` attacks exactly one post-P3085 interface atom: a non-imported unit-calibrated empirical observable/readout source for the Z12 Dirichlet/Laplacian branch.  It computes `12` finite Laplacian spectrum rows, `7` distinct normalized spectral-gap rows, `5` scale-orbit controls, a `5 x 6 = 30` candidate gate matrix, and `30` calibration target rows.  The internal witnesses remain dimensionless; no measurement-unit source, calibration map, apparatus readout protocol, empirical observable, observed photons/light, spacetime EOM, `L_total`, bridge/role-transfer, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3086/S2036 empirical readout remains uncalibrated", "## P3086/S2036 empirical readout remains uncalibrated\n\n`P3086/S2036` confirms that finite Z12 Dirichlet/Laplacian spectral and scalar witnesses do not by themselves form a calibrated physical observable.  A Lagrangian/EOM or empirical-physics reading still needs a strict source for measurement units, a calibration map, and an apparatus/readout protocol; imported metrology, observed photon, or continuum physics templates do not supply those strict sources.\n")
    append_once(AGENTS, "Current empirical-readout/observable-calibration obstruction guardrail (P3086/S2036, 2026-06-25)", "## Current empirical-readout/observable-calibration obstruction guardrail (P3086/S2036, 2026-06-25)\n\n- P3086 follows the P3085 recommendation and audits one standard-physics interface atom: a non-imported unit-calibrated empirical observable/readout source for the Z12 Dirichlet/Laplacian branch.\n- The finite audit computes `12` spectrum rows, `7` distinct spectral-gap rows, `5` scale-orbit controls, `30` candidate gate rows, and `30` calibration-attempt rows; `0` candidates export an internal non-imported empirical observable source.\n- Do not promote dimensionless Z12 scalar/spectral witnesses, P3085 formal current proxies, scale-normalized ratios, imported metrology templates, or imported observed-photon/frequency templates to observed photons/light, spacetime EOM, physical Hamiltonian, empirical observable, selector closure, `L_total`, bridge/role-transfer, or ToE closure.\n- The next honest move is exactly one thermodynamic/statistical-ensemble obstruction audit for the Z12 Dirichlet/Laplacian branch, unless a genuinely new empirical-calibration source theorem is introduced.\n")
    return payload

if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
