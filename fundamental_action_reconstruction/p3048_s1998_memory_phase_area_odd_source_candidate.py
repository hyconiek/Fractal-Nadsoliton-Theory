#!/usr/bin/env python3
"""P3048/S1998: memory phase-space area odd-source candidate.

P3047 says the next admissible memory-lag move must supply one concrete nonzero
strict inversion-odd source value coupled to the memory-lag sign, or pivot away.
P3048 supplies a concrete candidate value without promoting it to closure: the
oriented triangle-area pseudoscalar of the (K_i, M_i) phase curve,

    A_s = 1/2 sum_i det[[K_i, M_i, 1], [K_{i+s}, M_{i+s}, 1], [K_{i+2s}, M_{i+2s}, 1]].

This is not the P3044 pair commutator sum; it is a three-point oriented area
functional.  It is finite and inversion-odd under s -> -s, so it is a real new
computational hint.  It still does not export strict selector closure because
the orientation of the phase curve and the K/M receiver provenance are not a
nonpremise nadsoliton source law or a unit-bearing selector/readout theorem.
"""
from __future__ import annotations

import hashlib, json
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p3023_s1973_kernel_dissipation_time_order_candidate_obstruction import N
from p3038_s1988_viscous_retarded_chiral_selector_candidate_obstruction import memory_viscosity_trace
from p3044_s1994_memory_lag_commutator_source_candidate import kernel_vector
from p3047_s1997_memory_lag_polarity_source_law_obstruction import OUT as P3047

OUT = GEN / "p3048_s1998_memory_phase_area_odd_source_candidate.json"
MD = GEN / "p3048_s1998_memory_phase_area_odd_source_candidate.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"
TOL = 1e-12
LAGS = list(range(1, N // 2 + 1))


def triangle_area(k: list[float], m: list[float], lag: int) -> float:
    total = 0.0
    for i in range(N):
        j = (i + lag) % N
        ell = (i + 2 * lag) % N
        total += 0.5 * (
            k[i] * (m[j] - m[ell])
            + k[j] * (m[ell] - m[i])
            + k[ell] * (m[i] - m[j])
        )
    return total


def area_rows() -> list[dict[str, Any]]:
    k = kernel_vector()
    m = memory_viscosity_trace(k)
    rows = []
    for lag in LAGS:
        value = triangle_area(k, m, lag)
        inv_value = triangle_area(k, m, -lag)
        rows.append({
            "lag": lag,
            "oriented_area": round(value, 15),
            "inverted_lag_area": round(inv_value, 15),
            "finite_nonzero": abs(value) > TOL,
            "inversion_odd_verified": abs(value + inv_value) < TOL,
            "candidate_sign": 1 if value > TOL else (-1 if value < -TOL else 0),
            "accepted_strict_odd_source_value": False,
            "failure": "oriented phase-area value is nonzero and odd, but phase-curve orientation and K/M provenance are not strict source laws",
        })
    return rows


def source_acceptance_rows(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return [
        {
            "criterion": "concrete_formula_supplied",
            "satisfied": True,
            "detail": "A_s is an explicit three-point oriented area functional, not just an abstract odd slot",
        },
        {
            "criterion": "finite_nonzero_signed_value",
            "satisfied": any(row["finite_nonzero"] for row in rows),
            "detail": "at least one lag has a nonzero oriented area",
        },
        {
            "criterion": "inversion_odd_representation",
            "satisfied": all(row["inversion_odd_verified"] for row in rows),
            "detail": "s -> -s sends A_s to -A_s on every finite row",
        },
        {
            "criterion": "strict_nadsoliton_source_provenance",
            "satisfied": False,
            "detail": "the area is built from sampled K and P3038 memory trace, not exported by a strict source theorem",
        },
        {
            "criterion": "nonconventional_phase_curve_orientation",
            "satisfied": False,
            "detail": "the sign uses the chosen cyclic orientation of the (K,M) phase curve",
        },
        {
            "criterion": "unique_coupling_polarity_to_p3046",
            "satisfied": False,
            "detail": "no theorem selects one P3046 coupling polarity from A_s without importing orientation convention",
        },
        {
            "criterion": "unit_bearing_selector_readout_installation",
            "satisfied": False,
            "detail": "no QW-2191 discharge, physical readout, or action/EOM insertion is exported",
        },
    ]


def build_matrix() -> dict[str, Any]:
    read_json(P3047)
    rows = area_rows()
    acceptance = source_acceptance_rows(rows)
    obligations = [
        {"obligation": "p3047_missing_odd_value_targeted", "satisfied": True, "detail": "P3048 supplies a concrete candidate for the missing inversion-odd source value"},
        {"obligation": "three_point_area_not_pair_commutator_replay", "satisfied": True, "detail": "A_s uses three phase-curve points and oriented triangle area"},
        {"obligation": "finite_nonzero_odd_value", "satisfied": any(row["finite_nonzero"] for row in rows) and all(row["inversion_odd_verified"] for row in rows), "detail": "nonzero rows flip sign under lag inversion"},
        {"obligation": "strict_nadsoliton_source_law", "satisfied": False, "detail": "no theorem exports this area as a primordial nadsoliton source law"},
        {"obligation": "nonconventional_orientation_and_coupling", "satisfied": False, "detail": "phase-curve orientation and P3046 coupling polarity remain unselected"},
        {"obligation": "selector_readout_or_ltotal_installation", "satisfied": False, "detail": "no selector/readout, unit-bearing action, EOM, or L_total installation follows"},
    ]
    return {
        "object": "MemoryPhaseSpaceArea_OddSourceCandidateMatrix",
        "formula": "A_s=1/2*sum_i det[[K_i,M_i,1],[K_{i+s},M_{i+s},1],[K_{i+2s},M_{i+2s},1]]",
        "area_rows": rows,
        "source_acceptance_rows": acceptance,
        "proof_obligations": obligations,
        "finite_certificate": {
            "area_rows": len(rows),
            "finite_nonzero_area_rows": sum(1 for row in rows if row["finite_nonzero"]),
            "inversion_odd_rows": sum(1 for row in rows if row["inversion_odd_verified"]),
            "accepted_strict_odd_source_value_rows": sum(1 for row in rows if row["accepted_strict_odd_source_value"]),
            "source_acceptance_criteria": len(acceptance),
            "satisfied_source_acceptance_criteria": sum(1 for row in acceptance if row["satisfied"]),
            "proof_obligations": len(obligations),
            "satisfied_proof_obligations": sum(1 for row in obligations if row["satisfied"]),
            "strict_inversion_odd_source_value_exported": False,
        },
    }


def build_payload() -> dict[str, Any]:
    matrix = build_matrix()
    return {
        "status": "P3048_MEMORY_PHASE_AREA_ODD_SOURCE_CANDIDATE_BOUNDED_NO_EXPORT",
        "input_hashes": {"P3047": hashlib.sha256(P3047.read_bytes()).hexdigest() if P3047.exists() else None},
        "constructed_theoretical_objects": matrix,
        "finite_certificate": matrix["finite_certificate"],
        "decision": {
            "bounded_no_go": "P3048 supplies a concrete nonzero inversion-odd candidate value: the oriented triangle area of the (K,M) memory phase curve.  This is real computational progress beyond an abstract odd slot, but it remains a candidate receiver/source value rather than a strict source law because phase-curve orientation, K/M provenance, coupling polarity, and unit-bearing selector/readout installation are not exported.",
            "negative_export_flags": {k: False for k in ["strict_inversion_odd_source_value_exported", "strict_nadsoliton_source_law_exported", "coupling_polarity_selected", "selector_readout_coupling_exported", "qw2191_discharged", "strict_selector_closure_exported", "observed_physics_exported", "unit_bearing_action_eom_hamiltonian_exported", "ltotal_exported", "bridge_closure_exported", "role_transfer_exported", "toe_closure_exported"]},
            "next_honest_step": "Do not promote the oriented phase-area sign alone.  A next proof-grade move may attack exactly one missing premise for this candidate: strict provenance of the (K,M) phase curve as nadsoliton source data, or a nonconventional orientation/coupling-polarity theorem tying A_s to the P3046 torsor.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    c = payload["finite_certificate"]
    lines = [
        "# P3048/S1998 memory phase-space area odd-source candidate", "",
        f"Status: `{payload['status']}`", "", "## Finite certificate",
        f"- area rows: `{c['area_rows']}`",
        f"- finite nonzero area rows: `{c['finite_nonzero_area_rows']}`",
        f"- inversion-odd rows: `{c['inversion_odd_rows']}`",
        f"- accepted strict odd-source rows: `{c['accepted_strict_odd_source_value_rows']}`",
        f"- source acceptance criteria: `{c['satisfied_source_acceptance_criteria']}/{c['source_acceptance_criteria']}`",
        f"- satisfied proof obligations: `{c['satisfied_proof_obligations']}/{c['proof_obligations']}`",
        f"- strict inversion-odd source value exported: `{c['strict_inversion_odd_source_value_exported']}`", "",
        "## Decision", payload["decision"]["bounded_no_go"], "", "## Recommendation", payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3048/S1998 memory phase-space area odd-source candidate", "## P3048/S1998 memory phase-space area odd-source candidate\n\n`P3048/S1998` responds to the P3047 missing-object requirement by constructing a concrete inversion-odd candidate value: the three-point oriented triangle area of the `(K_i,M_i)` memory phase curve.  The finite audit finds nonzero area rows and exact sign flip under lag inversion, so this is a real computational hint beyond an abstract odd-source slot.  It still does not export a strict source law: the `(K,M)` phase-curve provenance, nonconventional orientation, P3046 coupling polarity, selector/readout installation, unit-bearing action/EOM, `L_total`, bridge/role transfer, and ToE closure remain unexported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3048/S1998 memory phase-area `L_total` guard", "## P3048/S1998 memory phase-area `L_total` guard\n\n`P3048/S1998` adds no physical `L_total` term.  The oriented `(K,M)` phase-area pseudoscalar is finite, nonzero, and inversion-odd, but it lacks strict phase-curve provenance, nonconventional orientation/coupling polarity, and unit-bearing variational/action/EOM installation.\n")
    append_once(AGENTS, "Current memory phase-area odd-source candidate guardrail (P3048/S1998, 2026-06-23)", "## Current memory phase-area odd-source candidate guardrail (P3048/S1998, 2026-06-23)\n\n- P3048 constructs one concrete inversion-odd candidate value after P3047: the three-point oriented triangle area of the `(K_i,M_i)` memory phase curve.\n- The finite area rows are nonzero and flip sign under lag inversion, but phase-curve provenance, nonconventional orientation, unique P3046 coupling polarity, and selector/readout installation remain unexported.\n- Do not promote oriented phase-area signs, cyclic phase-curve orientation, or area-score winners to `QW-2191` discharge, selector closure, observed physics, unit-bearing action/EOM, `L_total`, bridge/role-transfer, or ToE closure.\n- A next move may attack exactly one missing premise for this candidate: strict `(K,M)` phase-curve provenance or nonconventional orientation/coupling-polarity theorem; otherwise preserve bounded no-export.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
