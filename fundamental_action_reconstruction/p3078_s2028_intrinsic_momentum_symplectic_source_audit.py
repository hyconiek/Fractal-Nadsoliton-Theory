#!/usr/bin/env python3
"""P3078/S2028: intrinsic momentum/symplectic-source audit.

P3077 isolated the exact obstruction preventing the internal Z12 Dirichlet
source from becoming an internally sourced second-order/wave object: the formal
Hamiltonian lift imports momentum, symplectic pairing, kinetic normalization,
time units, and Lorentzian embedding.  P3078 executes the next bounded test by
constructing a finite candidate inventory for intrinsic momentum or
antisymmetric two-form sources already available in current nadsoliton/Z12
artifacts.  It deliberately avoids selector replay and does not promote formal
phase-space mathematics to observed light or standard-model physics.
"""
from __future__ import annotations

import hashlib, json, subprocess
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p3077_s2027_second_order_lift_obstruction_table import OUT as P3077

OUT = GEN / "p3078_s2028_intrinsic_momentum_symplectic_source_audit.json"
MD = GEN / "p3078_s2028_intrinsic_momentum_symplectic_source_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"
Z12_SIZE = 12
CONTENT_PATTERNS = {
    "momentum_symbols": r"\bpi\b|momentum|conjugate momentum|phase-space",
    "symplectic_two_form": r"symplectic|antisymmetric 2-form|two-form|Poisson|canonical bracket",
    "dirichlet_laplacian_source": r"Dirichlet|Laplacian|lambda_j|E_D\(rho\)",
    "no_physics_promotion": r"observed light|gauge photons|spacetime EOM|unit-bearing time|L_total|ToE|selector closure",
}

CANDIDATE_SOURCES = (
    {
        "id": "z12_shift_derivative_current",
        "object_type": "antisymmetric nearest-neighbor difference operator",
        "formula": "D = S - S^{-1}",
        "internally_defined_on_z12": True,
        "antisymmetric_pairing_available": True,
        "conjugate_pi_variable_exported": False,
        "nondegenerate_symplectic_form_exported": False,
        "hamiltonian_coupling_to_dirichlet_source_exported": False,
        "time_unit_exported": False,
        "blocker": "a skew operator on configuration space is not a canonical rho/pi phase-space two-form",
    },
    {
        "id": "fourier_quadrature_phase_pairing",
        "object_type": "modewise sine/cosine quadrature pairing",
        "formula": "(cos_j, sin_j) quadrature pairs for j=1..5 plus Nyquist handling",
        "internally_defined_on_z12": True,
        "antisymmetric_pairing_available": True,
        "conjugate_pi_variable_exported": False,
        "nondegenerate_symplectic_form_exported": False,
        "hamiltonian_coupling_to_dirichlet_source_exported": False,
        "time_unit_exported": False,
        "blocker": "quadrature labels provide representation bookkeeping, not an exported momentum source law",
    },
    {
        "id": "gradient_flow_velocity_proxy",
        "object_type": "first-order smoothing velocity proxy",
        "formula": "v = -L rho",
        "internally_defined_on_z12": True,
        "antisymmetric_pairing_available": False,
        "conjugate_pi_variable_exported": False,
        "nondegenerate_symplectic_form_exported": False,
        "hamiltonian_coupling_to_dirichlet_source_exported": False,
        "time_unit_exported": False,
        "blocker": "dissipative velocity proxy is not independent conjugate momentum and has no symplectic conservation law",
    },
    {
        "id": "chiral_bispectrum_sign_torsor",
        "object_type": "orientation-sensitive pseudoscalar sign marker",
        "formula": "sign(Im(B_{1,5}))",
        "internally_defined_on_z12": True,
        "antisymmetric_pairing_available": False,
        "conjugate_pi_variable_exported": False,
        "nondegenerate_symplectic_form_exported": False,
        "hamiltonian_coupling_to_dirichlet_source_exported": False,
        "time_unit_exported": False,
        "blocker": "signed torsor data is not a rho/pi phase-space variable and lacks the P2719/P2721 coupling/localizer theorem",
    },
    {
        "id": "formal_imported_canonical_phase_space",
        "object_type": "external canonical cotangent bundle ansatz",
        "formula": "T*R^12 with omega = sum_j d rho_j wedge d pi_j",
        "internally_defined_on_z12": False,
        "antisymmetric_pairing_available": True,
        "conjugate_pi_variable_exported": True,
        "nondegenerate_symplectic_form_exported": True,
        "hamiltonian_coupling_to_dirichlet_source_exported": True,
        "time_unit_exported": False,
        "blocker": "works only as imported continuum mechanics structure, not a strict nadsoliton/Z12 export",
    },
)

REQUIRED_GATES = (
    "internally_defined_on_z12",
    "conjugate_pi_variable_exported",
    "nondegenerate_symplectic_form_exported",
    "hamiltonian_coupling_to_dirichlet_source_exported",
    "time_unit_exported",
)


def content_grep() -> list[dict[str, Any]]:
    rows = []
    for lane, pattern in CONTENT_PATTERNS.items():
        proc = subprocess.run(["rg", "-n", pattern, "fundamental_action_reconstruction", "AGENTS.md", "-S"], cwd=REPO, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
        hits = [line for line in proc.stdout.splitlines() if line.strip()]
        rows.append({"content_lane": lane, "pattern": pattern, "hit_count": len(hits), "sample_hits": hits[:8]})
    return rows


def gate_rows() -> list[dict[str, Any]]:
    rows = []
    for candidate in CANDIDATE_SOURCES:
        for gate in REQUIRED_GATES:
            passed = bool(candidate[gate])
            rows.append({
                "candidate_source": candidate["id"],
                "required_gate": gate,
                "gate_passed": passed,
                "blocking_if_failed": not passed,
                "detail": "passed" if passed else candidate["blocker"],
            })
    return rows


def z12_skew_operator_rows() -> list[dict[str, Any]]:
    rows = []
    for j in range(Z12_SIZE):
        # D=S-S^-1 has Fourier eigenvalue 2i sin(2*pi*j/12).  This is a real
        # skew configuration-space current, not a nondegenerate rho/pi two-form.
        rows.append({
            "mode_j": j,
            "skew_current_eigenvalue_label": f"2*i*sin(2*pi*{j}/12)",
            "zero_skew_eigenvalue": j in (0, 6),
            "can_pair_as_canonical_momentum": False,
            "reason": "zero mode/Nyquist degeneracy and no independent pi source law",
        })
    return rows


def aggregate(gates: list[dict[str, Any]]) -> list[dict[str, Any]]:
    out = []
    for candidate in CANDIDATE_SOURCES:
        subset = [r for r in gates if r["candidate_source"] == candidate["id"]]
        out.append({
            "candidate_source": candidate["id"],
            "passed_gates": sum(1 for r in subset if r["gate_passed"]),
            "failed_gates": sum(1 for r in subset if not r["gate_passed"]),
            "accepted_intrinsic_momentum_symplectic_source": all(r["gate_passed"] for r in subset),
        })
    return out


def build_payload() -> dict[str, Any]:
    p3077 = read_json(P3077)
    greps = content_grep()
    gates = gate_rows()
    skew_rows = z12_skew_operator_rows()
    aggs = aggregate(gates)
    accepted = [a for a in aggs if a["accepted_intrinsic_momentum_symplectic_source"]]
    proof_obligations = [
        {"obligation": "content_first_grep_for_momentum_and_symplectic_sources", "satisfied": True, "detail": "searched momentum, symplectic/two-form, Dirichlet/Laplacian, and no-promotion lanes"},
        {"obligation": "construct_finite_candidate_source_inventory", "satisfied": True, "detail": "five candidate source classes are audited, including the strongest formal imported cotangent ansatz"},
        {"obligation": "execute_candidate_gate_matrix", "satisfied": True, "detail": "5 candidates x 5 required gates = 25 exact gate rows"},
        {"obligation": "execute_z12_skew_current_modal_table", "satisfied": True, "detail": "12 rows compute the internal skew-current eigenvalue labels and degeneracy blockers"},
        {"obligation": "export_accepted_intrinsic_phase_space_source", "satisfied": False, "detail": "no candidate passes all internal pi, nondegenerate symplectic, Hamiltonian-coupling, and time-unit gates"},
    ]
    return {
        "status": "P3078_INTRINSIC_MOMENTUM_SYMPLECTIC_SOURCE_NOT_EXPORTED_BOUNDED_NO_GO",
        "input_hashes": {"P3077": hashlib.sha256(P3077.read_bytes()).hexdigest() if P3077.exists() else None},
        "constructed_theoretical_objects": {
            "content_first_repo_grep": greps,
            "intrinsic_source_audit_object": {
                "object": "Z12IntrinsicMomentumSymplecticSourceAudit",
                "source_reused": "P3077 formal Hamiltonian lift obstruction boundary",
                "candidate_sources": [c["id"] for c in CANDIDATE_SOURCES],
                "required_gates": list(REQUIRED_GATES),
                "acceptance_predicate": "current artifacts internally define an independent pi variable, a nondegenerate antisymmetric two-form, a Hamiltonian coupling to E_D(rho), and unit-bearing time",
            },
            "candidate_gate_rows": gates,
            "z12_skew_current_modal_rows": skew_rows,
            "candidate_aggregate_certificate": aggs,
        },
        "finite_certificate": {
            "content_grep_lanes": len(greps),
            "content_grep_hits": sum(r["hit_count"] for r in greps),
            "p3077_accepted_internal_second_order_wave_rows": p3077.get("finite_certificate", {}).get("accepted_internal_second_order_wave_rows"),
            "candidate_sources": len(CANDIDATE_SOURCES),
            "required_gates": len(REQUIRED_GATES),
            "candidate_gate_rows": len(gates),
            "z12_skew_current_modal_rows": len(skew_rows),
            "z12_skew_current_zero_rows": sum(1 for r in skew_rows if r["zero_skew_eigenvalue"]),
            "accepted_intrinsic_momentum_symplectic_sources": len(accepted),
            "formal_imported_candidate_passed_gates": next(a["passed_gates"] for a in aggs if a["candidate_source"] == "formal_imported_canonical_phase_space"),
            "proof_obligations": len(proof_obligations),
            "satisfied_proof_obligations": sum(1 for o in proof_obligations if o["satisfied"]),
        },
        "proof_obligations": proof_obligations,
        "decision": {
            "bounded_result": "P3078 constructs a finite intrinsic momentum/symplectic-source audit.  Z12 supplies useful configuration-space skew currents and Fourier quadrature bookkeeping, and the imported cotangent ansatz is mathematically coherent, but no current artifact exports an independent internal pi variable, nondegenerate symplectic two-form, Hamiltonian coupling to the Dirichlet source, and unit-bearing time together.  The P3077 wave promotion therefore remains blocked.",
            "negative_export_flags": {k: False for k in ["intrinsic_momentum_source_exported", "intrinsic_symplectic_form_exported", "unit_bearing_time_exported", "internally_sourced_wave_eom_exported", "observed_light_exported", "gauge_photon_sector_exported", "spacetime_eom_exported", "empirical_physics_exported", "QW_2191_discharged", "strict_selector_closure_exported", "ltotal_exported", "bridge_role_transfer_exported", "toe_closure_exported"]},
            "positive_scoped_flags": {"z12_skew_current_table_constructed": True, "candidate_source_gate_matrix_executed": True, "second_order_wave_promotion_frozen_on_current_artifacts": True},
            "next_honest_step": "Freeze the Dirichlet-to-wave promotion on current artifacts and pivot to one different non-selector typed object: a bounded light-cone/causal-order compatibility audit for the internal Z12 dispersion data, testing whether any current artifact supplies a metric signature, finite propagation cone, or unit-normalized causal order without importing spacetime EOM.  If it does not, preserve the smoothing-only interpretation and do not claim observed light, gauge photons, empirical physics, selector closure, L_total, bridge/role-transfer, or ToE.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    c = payload["finite_certificate"]
    lines = [
        "# P3078/S2028 intrinsic momentum/symplectic-source audit", "", f"Status: `{payload['status']}`", "", "## Finite certificate",
        f"- content grep lanes: `{c['content_grep_lanes']}`",
        f"- content grep hits: `{c['content_grep_hits']}`",
        f"- P3077 accepted internal second-order wave rows: `{c['p3077_accepted_internal_second_order_wave_rows']}`",
        f"- candidate sources: `{c['candidate_sources']}`",
        f"- required gates: `{c['required_gates']}`",
        f"- candidate gate rows: `{c['candidate_gate_rows']}`",
        f"- Z12 skew-current modal rows: `{c['z12_skew_current_modal_rows']}`",
        f"- Z12 skew-current zero rows: `{c['z12_skew_current_zero_rows']}`",
        f"- accepted intrinsic momentum/symplectic sources: `{c['accepted_intrinsic_momentum_symplectic_sources']}`",
        f"- formal imported candidate passed gates: `{c['formal_imported_candidate_passed_gates']}`",
        f"- satisfied proof obligations: `{c['satisfied_proof_obligations']}/{c['proof_obligations']}`",
        "", "## Decision", payload["decision"]["bounded_result"], "", "## Recommendation", payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3078/S2028 intrinsic momentum/symplectic-source audit", "## P3078/S2028 intrinsic momentum/symplectic-source audit\n\n`P3078/S2028` constructs `Z12IntrinsicMomentumSymplecticSourceAudit` after the `P3077` second-order lift obstruction.  It audits `5` candidate momentum/symplectic source classes, `5` required gates, `25` gate rows, and `12` Z12 skew-current modal rows.  Internal Z12 skew currents and Fourier quadrature pairings are real configuration-space structures, but they do not export an independent `pi`, a nondegenerate rho/pi symplectic form, a Hamiltonian coupling theorem to `E_D(rho)`, or unit-bearing time.  The imported cotangent ansatz remains formal/imported, so no internally sourced wave EOM, observed light, gauge photons, empirical physics, selector closure, `L_total`, bridge/role-transfer, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3078/S2028 no intrinsic phase-space source", "## P3078/S2028 no intrinsic phase-space source\n\n`P3078/S2028` freezes the P3077 promotion attempt on current artifacts: the repository has Z12 skew-current and quadrature bookkeeping, but no strict internal source for `pi`, nondegenerate symplectic two-form, Hamiltonian coupling to the Dirichlet energy, or unit-bearing time.  Therefore `H = 1/2*pi^2 + E_D(rho)` remains an imported formal lift rather than a sourced Lagrangian/Hamiltonian EOM.\n")
    append_once(AGENTS, "Current intrinsic momentum/symplectic-source audit guardrail (P3078/S2028, 2026-06-24)", "## Current intrinsic momentum/symplectic-source audit guardrail (P3078/S2028, 2026-06-24)\n\n- P3078 follows the P3077 recommendation and audits current nadsoliton/Z12 artifacts for an intrinsic `pi` variable or nondegenerate symplectic/two-form source coupled to the Dirichlet source.\n- The finite matrix has `5` candidate source classes, `25` gate rows, and `12` Z12 skew-current modal rows; no candidate exports all required internal phase-space, Hamiltonian-coupling, and time-unit gates.\n- Do not promote Z12 skew currents, Fourier quadrature bookkeeping, or the imported cotangent ansatz to observed light, gauge photons, spacetime EOM, empirical physics, `QW-2191` discharge, selector closure, `L_total`, bridge/role-transfer, or ToE closure.\n- Freeze Dirichlet-to-wave promotion on current artifacts.  A next honest non-selector move may test one different typed object, such as a bounded light-cone/causal-order compatibility audit for internal Z12 dispersion data, without importing spacetime EOM.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
