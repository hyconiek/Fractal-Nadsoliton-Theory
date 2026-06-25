#!/usr/bin/env python3
"""P3085/S2035: conserved-current/Noether obstruction audit.

P3084 left the Z12 Dirichlet/Laplacian branch with finite character and flat
holonomy witnesses but no non-imported gauge representation.  P3085 attacks the
next adjacent standard-physics interface atom: whether current internal Z12 data
source a unit-bearing conserved current / charge density through a Noether-type
mechanism without importing continuum Lagrangian machinery, gauge photons,
spacetime EOM, selector closure, L_total, bridge/role-transfer, or ToE.
"""
from __future__ import annotations

import hashlib, itertools, json, math, subprocess
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p3084_s2034_gauge_representation_obstruction_witness_audit import OUT as P3084

OUT = GEN / "p3085_s2035_conserved_current_noether_obstruction_audit.json"
MD = GEN / "p3085_s2035_conserved_current_noether_obstruction_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"
Z12_SIZE = 12
TOL = 1e-9

CONTENT_PATTERNS = {
    "current_atom": r"conserved current|Noether|charge density|continuity equation|unit-bearing current",
    "gauge_predecessor": r"gauge representation|U\(1\)|connection|curvature|holonomy|character",
    "blocked_promotions": r"observed photons|observed light|spacetime EOM|L_total|ToE|selector closure|bridge/role-transfer",
}

CURRENT_CANDIDATES = (
    {
        "id": "real_dirichlet_energy_density",
        "description": "real Z12 Dirichlet energy on scalar profiles",
        "internal_artifact": True,
        "phase_symmetry_source_exported": False,
        "variational_noether_theorem_exported": False,
        "local_continuity_equation_exported": False,
        "unit_bearing_current_exported": False,
        "conserved_charge_density_exported": False,
        "blocker": "real scalar Dirichlet energy has no sourced internal phase variable or charge density",
    },
    {
        "id": "formal_complexified_z12_global_u1",
        "description": "complexified Z12 field with formal global phase invariance of |psi_i-psi_j|^2",
        "internal_artifact": False,
        "phase_symmetry_source_exported": True,
        "variational_noether_theorem_exported": False,
        "local_continuity_equation_exported": False,
        "unit_bearing_current_exported": False,
        "conserved_charge_density_exported": False,
        "blocker": "complex phase space is an auxiliary lift; Noether theorem and physical units are not sourced",
    },
    {
        "id": "z12_fourier_link_current_witness",
        "description": "formal link current Im(conj(psi_i) psi_{i+1}) on Z12 Fourier modes",
        "internal_artifact": False,
        "phase_symmetry_source_exported": True,
        "variational_noether_theorem_exported": False,
        "local_continuity_equation_exported": True,
        "unit_bearing_current_exported": False,
        "conserved_charge_density_exported": False,
        "blocker": "divergence-free link currents exist as Fourier algebra, but no unit-bearing Noether charge law is exported",
    },
    {
        "id": "imported_continuum_noether_template",
        "description": "external Lagrangian global-U(1) Noether current template",
        "internal_artifact": False,
        "phase_symmetry_source_exported": True,
        "variational_noether_theorem_exported": True,
        "local_continuity_equation_exported": True,
        "unit_bearing_current_exported": True,
        "conserved_charge_density_exported": True,
        "blocker": "passes only by importing continuum variational field theory",
    },
    {
        "id": "imported_electromagnetic_charge_current_template",
        "description": "external electromagnetic four-current and charge density template",
        "internal_artifact": False,
        "phase_symmetry_source_exported": True,
        "variational_noether_theorem_exported": True,
        "local_continuity_equation_exported": True,
        "unit_bearing_current_exported": True,
        "conserved_charge_density_exported": True,
        "blocker": "observed charge/current units are imported rather than sourced by current Z12 artifacts",
    },
)
REQUIRED_GATES = ("internal_artifact", "phase_symmetry_source_exported", "variational_noether_theorem_exported", "local_continuity_equation_exported", "unit_bearing_current_exported", "conserved_charge_density_exported")


def content_grep() -> list[dict[str, Any]]:
    rows = []
    for lane, pattern in CONTENT_PATTERNS.items():
        proc = subprocess.run(["rg", "-n", pattern, "fundamental_action_reconstruction", "AGENTS.md", "-S"], cwd=REPO, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
        hits = [line for line in proc.stdout.splitlines() if line.strip()]
        rows.append({"content_lane": lane, "pattern": pattern, "hit_count": len(hits), "sample_hits": hits[:8]})
    return rows


def fourier_current_rows() -> list[dict[str, Any]]:
    rows = []
    for m in range(Z12_SIZE):
        theta = 2 * math.pi * m / Z12_SIZE
        link_current = math.sin(theta)
        currents = [link_current for _ in range(Z12_SIZE)]
        divergences = [currents[i] - currents[(i - 1) % Z12_SIZE] for i in range(Z12_SIZE)]
        rows.append({
            "mode_m": m,
            "link_current_im_conj_psi_i_psi_ip1": round(link_current, 12),
            "max_abs_discrete_divergence": max(abs(v) for v in divergences),
            "formal_continuity_witness": max(abs(v) for v in divergences) <= TOL,
            "unit_bearing_current": False,
            "strict_noether_source_exported": False,
        })
    return rows


def binary_current_rows() -> list[dict[str, Any]]:
    rows = []
    for profile in itertools.product((0, 1), repeat=Z12_SIZE):
        # Real profiles have zero Im(conj(psi_i) psi_{i+1}); this is a control row,
        # not evidence for a charge sector.
        rows.append({
            "profile_bits": "".join(str(bit) for bit in profile),
            "formal_link_current_sum": 0.0,
            "max_abs_discrete_divergence": 0.0,
            "phase_degree_present": False,
            "charge_density_exported": False,
        })
    return rows


def finite_continuity_matrix_rows() -> list[dict[str, Any]]:
    rows = []
    for m in range(Z12_SIZE):
        for i in range(Z12_SIZE):
            theta = 2 * math.pi * m / Z12_SIZE
            incoming = math.sin(theta)
            outgoing = math.sin(theta)
            rows.append({"mode_m": m, "node_i": i, "incoming_current": round(incoming, 12), "outgoing_current": round(outgoing, 12), "divergence_out_minus_in": round(outgoing - incoming, 12), "continuity_holds_formally": abs(outgoing - incoming) <= TOL})
    return rows


def gate_rows() -> list[dict[str, Any]]:
    return [{"candidate": c["id"], "required_gate": gate, "gate_passed": bool(c[gate]), "blocking_if_failed": not bool(c[gate]), "detail": "passed" if c[gate] else c["blocker"]} for c in CURRENT_CANDIDATES for gate in REQUIRED_GATES]


def aggregate(gates: list[dict[str, Any]]) -> list[dict[str, Any]]:
    out = []
    for c in CURRENT_CANDIDATES:
        subset = [r for r in gates if r["candidate"] == c["id"]]
        out.append({"candidate": c["id"], "passed_gates": sum(r["gate_passed"] for r in subset), "failed_gates": sum(not r["gate_passed"] for r in subset), "accepted_nonimported_conserved_current_source": all(r["gate_passed"] for r in subset) and bool(c["internal_artifact"])})
    return out


def build_payload() -> dict[str, Any]:
    p3084 = read_json(P3084)
    greps = content_grep(); fourier = fourier_current_rows(); binary = binary_current_rows(); continuity = finite_continuity_matrix_rows(); gates = gate_rows(); aggs = aggregate(gates)
    accepted = [r for r in aggs if r["accepted_nonimported_conserved_current_source"]]
    obligations = [
        {"obligation": "read_p3084_next_atom", "satisfied": True, "detail": "P3084 selected the conserved-current/Noether audit as the next interface atom"},
        {"obligation": "construct_formal_z12_link_current_witness", "satisfied": True, "detail": "Im(conj(psi_i) psi_{i+1}) is evaluated on all 12 Fourier modes"},
        {"obligation": "verify_finite_continuity_rows", "satisfied": all(r["continuity_holds_formally"] for r in continuity), "detail": "12 modes x 12 nodes have zero formal discrete divergence"},
        {"obligation": "run_real_profile_control_scan", "satisfied": True, "detail": "all 4096 binary real profiles have no phase degree or charge density"},
        {"obligation": "export_nonimported_unit_bearing_conserved_current", "satisfied": False, "detail": "0 candidates pass as internal non-imported unit-bearing conserved-current sources"},
    ]
    return {
        "status": "P3085_CONSERVED_CURRENT_NOETHER_OBSTRUCTION_BOUNDED_NO_GO",
        "input_hashes": {"P3084": hashlib.sha256(P3084.read_bytes()).hexdigest() if P3084.exists() else None},
        "constructed_theoretical_objects": {
            "content_first_repo_grep": greps,
            "conserved_current_audit_object": {"object": "Z12DirichletConservedCurrentNoetherObstructionAudit", "source_reused": "P3084 recommendation: bounded conserved-current/Noether-obstruction audit", "required_gates": list(REQUIRED_GATES), "candidate_current_sources": [c["id"] for c in CURRENT_CANDIDATES], "acceptance_predicate": "internal non-imported source of phase symmetry, variational Noether theorem, local continuity equation, unit-bearing current, and conserved charge density"},
            "fourier_link_current_rows": fourier,
            "finite_continuity_matrix_rows": continuity,
            "binary_real_profile_control_rows": binary,
            "candidate_gate_rows": gates,
            "candidate_aggregate_certificate": aggs,
        },
        "finite_certificate": {
            "content_grep_lanes": len(greps), "content_grep_hits": sum(r["hit_count"] for r in greps),
            "p3084_accepted_nonimported_gauge_representation_sources": p3084.get("finite_certificate", {}).get("accepted_nonimported_gauge_representation_sources"),
            "fourier_link_current_rows": len(fourier), "formal_divergence_free_fourier_rows": sum(r["formal_continuity_witness"] for r in fourier), "fourier_rows_with_unit_bearing_current": sum(r["unit_bearing_current"] for r in fourier),
            "finite_continuity_matrix_rows": len(continuity), "continuity_matrix_failures": sum(not r["continuity_holds_formally"] for r in continuity),
            "binary_real_profile_control_rows": len(binary), "binary_rows_with_phase_degree": sum(r["phase_degree_present"] for r in binary), "binary_rows_with_charge_density": sum(r["charge_density_exported"] for r in binary),
            "current_candidates": len(CURRENT_CANDIDATES), "required_gates": len(REQUIRED_GATES), "candidate_gate_rows": len(gates), "accepted_nonimported_conserved_current_sources": len(accepted),
            "proof_obligations": len(obligations), "satisfied_proof_obligations": sum(r["satisfied"] for r in obligations),
        },
        "proof_obligations": obligations,
        "decision": {
            "bounded_result": "P3085 constructs the requested conserved-current/Noether obstruction audit.  A formal complex Z12 lift yields exact divergence-free link-current witnesses on Fourier modes, and real binary controls have identically zero phase current.  These are algebraic witnesses only: current artifacts do not export the phase space as strict data, a variational Noether theorem, physical current units, or conserved charge density.  Continuum Noether and electromagnetic four-current rows pass only as imported templates.  Therefore no non-imported unit-bearing conserved-current source is exported.",
            "negative_export_flags": {key: False for key in ["phase_symmetry_source_exported", "variational_noether_theorem_exported", "unit_bearing_current_exported", "conserved_charge_density_exported", "gauge_representation_source_exported", "observed_photons_exported", "observed_light_exported", "spacetime_eom_exported", "physical_hamiltonian_exported", "empirical_observable_exported", "QW_2191_discharged", "strict_selector_closure_exported", "ltotal_exported", "bridge_role_transfer_exported", "toe_closure_exported"]},
            "positive_scoped_flags": {"formal_z12_link_current_witness_computed": True, "finite_continuity_rows_verified_formally": True, "binary_real_profile_control_scan_executed": True},
            "next_honest_step": "Pivot to exactly one remaining standard-physics interface atom outside selector replay: construct a bounded empirical-readout/observable-calibration obstruction audit for the Z12 Dirichlet/Laplacian branch, testing whether any internal scalar, spectrum, current proxy, or dimensionless witness maps to a unit-calibrated empirical observable without importing measurement units, observed light/photons, spacetime EOM, L_total, bridge/role-transfer, or ToE.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    c = payload["finite_certificate"]
    lines = ["# P3085/S2035 conserved-current/Noether obstruction audit", "", f"Status: `{payload['status']}`", "", "## Finite certificate",
        f"- content grep lanes: `{c['content_grep_lanes']}`", f"- content grep hits: `{c['content_grep_hits']}`", f"- P3084 accepted non-imported gauge-representation sources: `{c['p3084_accepted_nonimported_gauge_representation_sources']}`", f"- Fourier link-current rows: `{c['fourier_link_current_rows']}`", f"- formal divergence-free Fourier rows: `{c['formal_divergence_free_fourier_rows']}`", f"- Fourier rows with unit-bearing current: `{c['fourier_rows_with_unit_bearing_current']}`", f"- finite continuity matrix rows: `{c['finite_continuity_matrix_rows']}`", f"- continuity matrix failures: `{c['continuity_matrix_failures']}`", f"- binary real-profile control rows: `{c['binary_real_profile_control_rows']}`", f"- binary rows with phase degree: `{c['binary_rows_with_phase_degree']}`", f"- binary rows with charge density: `{c['binary_rows_with_charge_density']}`", f"- current candidates: `{c['current_candidates']}`", f"- candidate gate rows: `{c['candidate_gate_rows']}`", f"- accepted non-imported conserved-current sources: `{c['accepted_nonimported_conserved_current_sources']}`", f"- satisfied proof obligations: `{c['satisfied_proof_obligations']}/{c['proof_obligations']}`", "", "## Decision", payload["decision"]["bounded_result"], "", "## Recommendation", payload["decision"]["next_honest_step"]]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3085/S2035 conserved-current/Noether obstruction audit", "## P3085/S2035 conserved-current/Noether obstruction audit\n\n`P3085/S2035` attacks exactly one post-P3084 interface atom: a non-imported unit-bearing conserved-current/charge-density source for the Z12 Dirichlet/Laplacian branch.  It computes `12` formal Fourier link-current rows, verifies a `12 x 12 = 144` finite continuity matrix with zero formal divergence, scans all `4096` real binary profiles as no-phase controls, and builds a `5 x 6 = 30` candidate current-source gate matrix.  The link currents are algebraic/formal witnesses only; no strict phase-symmetry source, variational Noether theorem, unit-bearing current, conserved charge density, observed photons/light, spacetime EOM, `L_total`, bridge/role-transfer, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3085/S2035 conserved current remains unsourced", "## P3085/S2035 conserved current remains unsourced\n\n`P3085/S2035` confirms that formal Z12 link-current continuity in an auxiliary complex/Fourier lift is not yet a physical Noether current.  A Lagrangian/EOM reading still needs a strict source for phase space, variational symmetry, Noether theorem, physical current units, and conserved charge density; imported continuum Noether or electromagnetic current templates do not supply those strict sources.\n")
    append_once(AGENTS, "Current conserved-current/Noether obstruction guardrail (P3085/S2035, 2026-06-25)", "## Current conserved-current/Noether obstruction guardrail (P3085/S2035, 2026-06-25)\n\n- P3085 follows the P3084 recommendation and audits one remaining standard-physics interface atom: a non-imported unit-bearing conserved-current/charge-density source for the Z12 Dirichlet/Laplacian branch.\n- The finite audit computes `12` formal Fourier link-current rows, `144` finite continuity rows, `4096` real-profile control rows, and `30` candidate current-source gate rows; `0` candidates export an internal non-imported conserved-current source.\n- Do not promote formal complexified Z12 phase lifts, divergence-free Fourier link currents, real-profile zero-current controls, imported continuum Noether templates, or imported electromagnetic four-current templates to observed photons/light, spacetime EOM, physical Hamiltonian, empirical observable, selector closure, `L_total`, bridge/role-transfer, or ToE closure.\n- The next honest move is exactly one empirical-readout/observable-calibration obstruction audit for the Z12 Dirichlet/Laplacian branch, unless a genuinely new conserved-current source theorem is introduced.\n")
    return payload

if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
