#!/usr/bin/env python3
"""P3038/S1988: integrated viscous-retarded-chiral selector candidate.

This is the direct follow-up to P3037: construct one integrated candidate
operator rather than forcing the selector into a familiar prechosen schema.  The
operator combines memory/viscosity, c-retardation anisotropy, an inversion-odd
chiral projection, and a unit/readout coupling slot.  It produces a nonzero
finite branch-separating score, but the sign/localizer, unit, and readout
coupling are still not sourced, so it remains an obstruction certificate.
"""
from __future__ import annotations

import hashlib, json, math
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p3023_s1973_kernel_dissipation_time_order_candidate_obstruction import N, k_strict
from p3037_s1987_selector_mechanism_hint_sheaf_acceptance_matrix import OUT as P3037

OUT = GEN / "p3038_s1988_viscous_retarded_chiral_selector_candidate_obstruction.json"
MD = GEN / "p3038_s1988_viscous_retarded_chiral_selector_candidate_obstruction.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

MEMORY_LAMBDA = 0.35
RETARD_RHO = 0.25
RETARD_OMEGA = 0.18575
UNIT_COUPLING_SLOT = 1.0
TOL = 1e-12


def kernel_vector() -> list[float]:
    return [k_strict(label) for label in range(1, N + 1)]


def memory_viscosity_trace(values: list[float], lam: float = MEMORY_LAMBDA) -> list[float]:
    trace = []
    state = 0.0
    for value in values:
        state = lam * value + (1.0 - lam) * state
        trace.append(state)
    return trace


def retardation_anisotropy(trace: list[float], rho: float = RETARD_RHO, omega: float = RETARD_OMEGA) -> list[float]:
    return [math.cos(omega * (1.0 + rho * value)) - math.cos(omega * (1.0 - rho * value)) for value in trace]


def chiral_projector() -> list[float]:
    # This is intentionally an inversion-odd probe, not a sourced orientation law.
    return [math.sin(2.0 * math.pi * label / N) for label in range(1, N + 1)]


def integrated_density(branch: int, unit_coupling: float = UNIT_COUPLING_SLOT) -> list[float]:
    values = kernel_vector()
    trace = memory_viscosity_trace(values)
    delta = retardation_anisotropy(trace)
    chi = chiral_projector()
    return [branch * unit_coupling * d * c for d, c in zip(delta, chi)]


def branch_score(branch: int) -> float:
    return sum(integrated_density(branch))


def build_matrix() -> dict[str, Any]:
    plus_score = branch_score(+1)
    minus_score = branch_score(-1)
    separation = plus_score - minus_score
    candidate_features = {
        "damping_memory_feedback_hint": True,
        "retardation_anisotropy_anchor_hint": True,
        "inversion_odd_or_chiral_sign_hint": True,
        "nonpremise_source_localizer_needed": True,
        "coupling_polarity_needed": True,
        "readout_unit_coupling_needed": True,
    }
    operator_rows = [
        {
            "row": "memory_viscosity_trace",
            "formula": "M_i = lambda*K_i + (1-lambda)*M_{i-1}",
            "finite_nonzero": any(abs(x) > TOL for x in memory_viscosity_trace(kernel_vector())),
            "accepted_as_source": False,
            "failure": "memory smoothing is a candidate dynamics, not a strict exported viscosity law",
        },
        {
            "row": "c_retardation_anisotropy_split",
            "formula": "Delta_i = cos(omega*(1+rho*M_i)) - cos(omega*(1-rho*M_i))",
            "finite_nonzero": any(abs(x) > TOL for x in retardation_anisotropy(memory_viscosity_trace(kernel_vector()))),
            "accepted_as_source": False,
            "failure": "anisotropic path split is inserted as a candidate slot, not sourced path geometry",
        },
        {
            "row": "inversion_odd_chiral_projection",
            "formula": "chi_i = sin(2*pi*i/12)",
            "finite_nonzero": any(abs(x) > TOL for x in chiral_projector()),
            "accepted_as_source": False,
            "failure": "the sine orientation is an inversion-odd probe but no nonpremise sign/localizer theorem selects its polarity",
        },
        {
            "row": "unit_readout_coupling_slot",
            "formula": "U_readout = 1 placeholder",
            "finite_nonzero": UNIT_COUPLING_SLOT != 0,
            "accepted_as_source": False,
            "failure": "unit/readout coupling is a dimensionless placeholder after P3036, not a physical unit theorem",
        },
    ]
    obligations = [
        {"obligation": "integrated_candidate_operator_constructed", "satisfied": True, "detail": "memory, retardation anisotropy, chiral projection, and unit/readout slot are multiplied into one density"},
        {"obligation": "all_p3037_hint_features_present", "satisfied": all(candidate_features.values()), "detail": "candidate covers the six P3037 hint features"},
        {"obligation": "finite_branch_score_nonzero", "satisfied": abs(plus_score) > TOL and abs(minus_score) > TOL, "detail": "both branch scores are nonzero and opposite"},
        {"obligation": "branch_separation_computable", "satisfied": abs(separation) > TOL, "detail": "plus and minus branch scores differ"},
        {"obligation": "nonpremise_chiral_sign_localizer", "satisfied": False, "detail": "chi_i is an inserted inversion-odd probe; no strict source/localizer fixes its sign"},
        {"obligation": "sourced_retardation_path_anisotropy", "satisfied": False, "detail": "rho and parallel/perpendicular path split are candidate slots, not exported geometry"},
        {"obligation": "physical_unit_readout_coupling", "satisfied": False, "detail": "unit coupling remains dimensionless after P3036"},
        {"obligation": "nonproxy_variational_or_readout_theorem", "satisfied": False, "detail": "no theorem installs the density into L_total or a classical readout row"},
    ]
    return {
        "object": "ViscousRetardedChiralReadoutSelector_CandidateObstructionMatrix",
        "operator_formula": "S_b = sum_i b * U * [cos(omega*(1+rho*M_i))-cos(omega*(1-rho*M_i))] * sin(2*pi*i/12)",
        "parameters": {"memory_lambda": MEMORY_LAMBDA, "retard_rho": RETARD_RHO, "retard_omega": RETARD_OMEGA, "unit_coupling_slot": UNIT_COUPLING_SLOT},
        "candidate_features": candidate_features,
        "operator_rows": operator_rows,
        "branch_scores": {"plus": round(plus_score, 15), "minus": round(minus_score, 15), "separation": round(separation, 15)},
        "proof_obligations": obligations,
        "finite_certificate": {
            "operator_rows": len(operator_rows),
            "finite_nonzero_operator_rows": sum(1 for row in operator_rows if row["finite_nonzero"]),
            "accepted_operator_source_rows": sum(1 for row in operator_rows if row["accepted_as_source"]),
            "p3037_features_present": sum(1 for value in candidate_features.values() if value),
            "branch_score_nonzero": abs(plus_score) > TOL and abs(minus_score) > TOL,
            "branch_separation_abs": round(abs(separation), 15),
            "proof_obligations": len(obligations),
            "satisfied_proof_obligations": sum(1 for row in obligations if row["satisfied"]),
            "strict_selector_mechanism_exported": all(row["satisfied"] for row in obligations) and all(row["accepted_as_source"] for row in operator_rows),
        },
    }


def build_payload() -> dict[str, Any]:
    read_json(P3037)
    matrix = build_matrix()
    return {
        "status": "P3038_VISCOUS_RETARDED_CHIRAL_SELECTOR_CANDIDATE_OBSTRUCTION_NO_SELECTOR_EXPORT",
        "input_hashes": {"P3037": hashlib.sha256(P3037.read_bytes()).hexdigest() if P3037.exists() else None},
        "constructed_theoretical_objects": matrix,
        "finite_certificate": matrix["finite_certificate"],
        "decision": {
            "bounded_no_go": "The integrated candidate operator is real and finite: it combines memory/viscosity, c-retardation anisotropy, an inversion-odd chiral projection, and a unit/readout coupling slot, yielding nonzero opposite branch scores.  This is progress beyond hint collection, but the sign/localizer, path anisotropy, physical unit coupling, and nonproxy variational/readout theorem remain unsourced.  Therefore no strict selector mechanism is exported.",
            "negative_export_flags": {k: False for k in ["strict_selector_mechanism_exported", "qw2191_discharged", "strict_selector_branch_source_exported", "observed_physics_exported", "unit_bearing_action_eom_hamiltonian_exported", "ltotal_exported", "bridge_closure_exported", "role_transfer_exported", "toe_closure_exported"]},
            "next_honest_step": "Do not promote the branch-separating score alone to selector closure. The next proof-grade move should attack exactly one missing source premise of this integrated operator: either a nonpremise chiral sign/localizer theorem, a sourced retardation path-anisotropy theorem, or a physical unit/readout coupling theorem.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    c = payload["finite_certificate"]
    scores = payload["constructed_theoretical_objects"]["branch_scores"]
    lines = [
        "# P3038/S1988 viscous-retarded-chiral selector candidate obstruction", "",
        f"Status: `{payload['status']}`", "", "## Finite certificate",
        f"- operator rows: `{c['operator_rows']}`",
        f"- finite nonzero operator rows: `{c['finite_nonzero_operator_rows']}`",
        f"- accepted operator source rows: `{c['accepted_operator_source_rows']}`",
        f"- P3037 features present: `{c['p3037_features_present']}`",
        f"- plus branch score: `{scores['plus']}`",
        f"- minus branch score: `{scores['minus']}`",
        f"- branch separation abs: `{c['branch_separation_abs']}`",
        f"- satisfied proof obligations: `{c['satisfied_proof_obligations']}/{c['proof_obligations']}`",
        f"- strict selector mechanism exported: `{c['strict_selector_mechanism_exported']}`", "",
        "## Decision", payload["decision"]["bounded_no_go"], "", "## Recommendation", payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3038/S1988 viscous-retarded-chiral selector candidate obstruction", "## P3038/S1988 viscous-retarded-chiral selector candidate obstruction\n\n`P3038/S1988` executes the P3037 recommendation by constructing one integrated candidate operator: a memory/viscosity trace, c-retardation anisotropy split, inversion-odd chiral projection, and unit/readout coupling slot.  The finite score is genuinely branch-separating (`+` and `-` branches get nonzero opposite scores), so this is stronger than a hint sheaf alone.  The obstruction is also explicit: the chiral sign/localizer, retardation path anisotropy, physical unit/readout coupling, and nonproxy variational/readout theorem are still inserted candidate slots rather than exported strict sources.  No `QW-2191` discharge, selector closure, observed physics, `L_total`, bridge/role transfer, or ToE closure follows.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3038/S1988 integrated selector candidate `L_total` guard", "## P3038/S1988 integrated selector candidate `L_total` guard\n\n`P3038/S1988` adds no physical `L_total` term.  Although the integrated density gives a finite branch-separating score, its chiral sign, path anisotropy, unit/readout coupling, and variational/readout insertion remain unsourced candidate slots.\n")
    append_once(AGENTS, "Current integrated viscous-retarded-chiral selector candidate guardrail (P3038/S1988, 2026-06-23)", "## Current integrated viscous-retarded-chiral selector candidate guardrail (P3038/S1988, 2026-06-23)\n\n- P3038 constructs one integrated candidate operator combining memory/viscosity, c-retardation anisotropy, an inversion-odd chiral projection, and a unit/readout coupling slot.\n- The finite computation gives nonzero opposite branch scores, but all source-critical components remain candidate slots: nonpremise chiral sign/localizer, sourced path anisotropy, physical unit/readout coupling, and nonproxy variational/readout insertion are not exported.\n- Do not promote the branch-separating score alone to `QW-2191` discharge, selector closure, observed physics, unit-bearing action/EOM, `L_total`, bridge/role-transfer, or ToE closure.\n- A next move may attack exactly one missing source premise of this integrated operator: chiral sign/localizer theorem, retardation path-anisotropy theorem, or physical unit/readout coupling theorem.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
