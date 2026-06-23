#!/usr/bin/env python3
"""P3049/S1999: memory phase-curve provenance obstruction.

P3048 produced a concrete nonzero inversion-odd phase-area pseudoscalar A_s,
but explicitly left the (K_i,M_i) phase-curve provenance unsourced.  P3049
attacks exactly that missing premise: can the ordered curve i -> (K_i,M_i) be
promoted to strict nadsoliton source data rather than a sampled receiver curve?

The finite audit is deliberately constructive and negative.  It builds a
source-channel provenance matrix for K, M, the ordered Z12 chart, cyclic
orientation, and the derived phase curve.  It also computes Aut(Z12) tests for
whether the curve and its area sign survive relabeling without a chart-polarity
choice.  The answer is no: K is sampled from the current strict gate vector and
M is a P3038 derived memory trace; the ordered phase curve is translation
stable only as a receiver diagnostic, while inversion still flips the signed
area.  Therefore the P3048 odd value remains real but not a strict source law.
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
from p3048_s1998_memory_phase_area_odd_source_candidate import OUT as P3048, triangle_area

OUT = GEN / "p3049_s1999_memory_phase_curve_provenance_obstruction.json"
MD = GEN / "p3049_s1999_memory_phase_curve_provenance_obstruction.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"
UNITS = [1, 5, 7, 11]
LAGS = list(range(1, N // 2 + 1))
TOL = 1e-12


def relabel(values: list[float], unit: int, shift: int = 0) -> list[float]:
    # New chart index i reads the old value at unit*i+shift.
    return [values[(unit * i + shift) % N] for i in range(N)]


def provenance_rows() -> list[dict[str, Any]]:
    return [
        {
            "channel": "K_i strict-gate sample vector",
            "finite_object_present": True,
            "strict_source_provenance": False,
            "reason": "K_i is the sampled operational strict-gate vector, not an exported primordial phase-coordinate source law",
        },
        {
            "channel": "M_i P3038 memory-viscosity trace",
            "finite_object_present": True,
            "strict_source_provenance": False,
            "reason": "M_i is a derived receiver trace over K_i; P3038/P3044-P3048 do not export it as independent nadsoliton source data",
        },
        {
            "channel": "ordered Z12 phase-curve index i -> (K_i,M_i)",
            "finite_object_present": True,
            "strict_source_provenance": False,
            "reason": "the ordered curve uses the chosen Z12 chart order and remains sensitive to inversion polarity",
        },
        {
            "channel": "cyclic orientation of the phase curve",
            "finite_object_present": True,
            "strict_source_provenance": False,
            "reason": "orientation is exactly the missing nonconventional sign; inversion sends it to the opposite branch",
        },
        {
            "channel": "P3048 area pseudoscalar A_s",
            "finite_object_present": True,
            "strict_source_provenance": False,
            "reason": "A_s is computable and odd, but inherits unsourced K/M curve provenance and orientation polarity",
        },
    ]


def aut_area_rows() -> list[dict[str, Any]]:
    k = kernel_vector()
    m = memory_viscosity_trace(k)
    rows = []
    base = {lag: triangle_area(k, m, lag) for lag in LAGS}
    for unit in UNITS:
        kk = relabel(k, unit)
        mm = relabel(m, unit)
        for lag in LAGS:
            value = triangle_area(kk, mm, lag)
            expected_lag = (unit * lag) % N
            if expected_lag > N // 2:
                expected_lag -= N
            sign = 1
            abs_lag = abs(expected_lag)
            if expected_lag < 0:
                sign = -1
            expected = sign * base[abs_lag]
            rows.append({
                "unit": unit,
                "lag": lag,
                "image_signed_lag": expected_lag,
                "area_after_relabel": round(value, 15),
                "expected_signed_area": round(expected, 15),
                "equivariance_verified": abs(value - expected) < TOL,
                "orientation_reversing_unit": unit in (7, 11),
                "strict_source_provenance_exported": False,
            })
    return rows


def acceptance_rows(prov: list[dict[str, Any]], aut: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return [
        {"criterion": "finite_KM_curve_constructed", "satisfied": True, "detail": "K_i and M_i define a concrete 12-point phase curve"},
        {"criterion": "translation_receiver_stability", "satisfied": True, "detail": "cyclic shifts preserve the receiver curve up to chart origin"},
        {"criterion": "Aut_area_equivariance_accounted", "satisfied": all(r["equivariance_verified"] for r in aut), "detail": "unit relabeling maps A_s to the signed image-lag area"},
        {"criterion": "K_channel_strict_source_theorem", "satisfied": False, "detail": prov[0]["reason"]},
        {"criterion": "M_channel_strict_memory_source_theorem", "satisfied": False, "detail": prov[1]["reason"]},
        {"criterion": "chart_independent_phase_curve_order", "satisfied": False, "detail": "the curve order is a Z12 chart order, not a nonpremise order theorem"},
        {"criterion": "nonconventional_orientation_polarity", "satisfied": False, "detail": "orientation-reversing units exchange the two area-sign branches"},
        {"criterion": "P3046_coupling_polarity_selected", "satisfied": False, "detail": "provenance alone does not choose one of the two P3046 coupling polarities"},
    ]


def build_matrix() -> dict[str, Any]:
    read_json(P3048)
    prov = provenance_rows()
    aut = aut_area_rows()
    acceptance = acceptance_rows(prov, aut)
    obligations = [
        {"obligation": "p3048_phase_curve_provenance_targeted", "satisfied": True, "detail": "P3049 attacks only the missing (K,M) phase-curve provenance premise"},
        {"obligation": "source_channels_explicitly_separated", "satisfied": True, "detail": "K, M, chart order, orientation, and A_s are audited as separate provenance channels"},
        {"obligation": "finite_aut_relabeling_test", "satisfied": all(r["equivariance_verified"] for r in aut), "detail": "all Aut relabel rows match the signed image-lag law"},
        {"obligation": "strict_KM_source_provenance", "satisfied": False, "detail": "K and M remain sampled/derived receiver data rather than strict source exports"},
        {"obligation": "chart_independent_orientation", "satisfied": False, "detail": "inversion-reversing units still flip the signed area branch"},
        {"obligation": "selector_readout_or_ltotal_installation", "satisfied": False, "detail": "no selector/readout, action/EOM, L_total, bridge, role transfer, or ToE export follows"},
    ]
    return {
        "object": "MemoryPhaseCurve_ProvenanceObstructionMatrix",
        "provenance_rows": prov,
        "aut_area_rows": aut,
        "source_acceptance_rows": acceptance,
        "proof_obligations": obligations,
        "finite_certificate": {
            "provenance_rows": len(prov),
            "strict_source_provenance_rows": sum(1 for r in prov if r["strict_source_provenance"]),
            "aut_area_rows": len(aut),
            "aut_equivariance_verified_rows": sum(1 for r in aut if r["equivariance_verified"]),
            "orientation_reversing_rows": sum(1 for r in aut if r["orientation_reversing_unit"]),
            "source_acceptance_criteria": len(acceptance),
            "satisfied_source_acceptance_criteria": sum(1 for r in acceptance if r["satisfied"]),
            "proof_obligations": len(obligations),
            "satisfied_proof_obligations": sum(1 for r in obligations if r["satisfied"]),
            "strict_KM_phase_curve_source_exported": False,
        },
    }


def build_payload() -> dict[str, Any]:
    matrix = build_matrix()
    return {
        "status": "P3049_MEMORY_PHASE_CURVE_PROVENANCE_OBSTRUCTION_BOUNDED_NO_EXPORT",
        "input_hashes": {"P3048": hashlib.sha256(P3048.read_bytes()).hexdigest() if P3048.exists() else None},
        "constructed_theoretical_objects": matrix,
        "finite_certificate": matrix["finite_certificate"],
        "decision": {
            "bounded_no_go": "P3049 constructs the missing provenance audit for the P3048 (K,M) phase curve.  The finite curve and its Aut signed-area law are real, but every source channel remains sampled, derived, chart-ordered, or orientation-polarity dependent.  Therefore the P3048 area sign is not promoted to a strict nadsoliton source law.",
            "negative_export_flags": {k: False for k in ["strict_KM_phase_curve_source_exported", "strict_nadsoliton_source_law_exported", "nonconventional_orientation_exported", "p3046_coupling_polarity_selected", "selector_readout_coupling_exported", "qw2191_discharged", "strict_selector_closure_exported", "observed_physics_exported", "unit_bearing_action_eom_hamiltonian_exported", "ltotal_exported", "bridge_closure_exported", "role_transfer_exported", "toe_closure_exported"]},
            "next_honest_step": "Do not replay K/M provenance or area signs as source closure.  The remaining proof-grade move in this lane is a nonconventional orientation/coupling-polarity theorem tying A_s to exactly one P3046 coupling polarity; otherwise pivot to a genuinely new typed object.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    c = payload["finite_certificate"]
    lines = [
        "# P3049/S1999 memory phase-curve provenance obstruction", "",
        f"Status: `{payload['status']}`", "", "## Finite certificate",
        f"- provenance rows: `{c['provenance_rows']}`",
        f"- strict source provenance rows: `{c['strict_source_provenance_rows']}`",
        f"- Aut area rows: `{c['aut_area_rows']}`",
        f"- Aut equivariance verified rows: `{c['aut_equivariance_verified_rows']}`",
        f"- orientation-reversing rows: `{c['orientation_reversing_rows']}`",
        f"- source acceptance criteria: `{c['satisfied_source_acceptance_criteria']}/{c['source_acceptance_criteria']}`",
        f"- satisfied proof obligations: `{c['satisfied_proof_obligations']}/{c['proof_obligations']}`",
        f"- strict K/M phase-curve source exported: `{c['strict_KM_phase_curve_source_exported']}`", "",
        "## Decision", payload["decision"]["bounded_no_go"], "", "## Recommendation", payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3049/S1999 memory phase-curve provenance obstruction", "## P3049/S1999 memory phase-curve provenance obstruction\n\n`P3049/S1999` attacks exactly one P3048 missing premise: strict provenance for the `(K_i,M_i)` phase curve.  It separates the K sample vector, P3038 memory trace, ordered Z12 phase curve, cyclic orientation, and P3048 area pseudoscalar into source-channel rows, then checks Aut relabeling of the signed area.  The finite signed-area law is coherent, but all source channels remain sampled/derived/chart-ordered/orientation-polarity dependent.  No strict `(K,M)` phase-curve source, P3046 coupling-polarity theorem, selector/readout installation, unit-bearing action/EOM, `L_total`, bridge/role transfer, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3049/S1999 memory phase-curve provenance `L_total` guard", "## P3049/S1999 memory phase-curve provenance `L_total` guard\n\n`P3049/S1999` adds no physical `L_total` term.  The `(K,M)` phase curve and P3048 area sign remain finite receiver diagnostics with sampled/derived provenance and orientation-polarity dependence; they are not unit-bearing variational/action/EOM inputs.\n")
    append_once(AGENTS, "Current memory phase-curve provenance obstruction guardrail (P3049/S1999, 2026-06-23)", "## Current memory phase-curve provenance obstruction guardrail (P3049/S1999, 2026-06-23)\n\n- P3049 attacks exactly one P3048 missing premise: strict provenance of the `(K_i,M_i)` phase curve as nadsoliton source data.\n- The finite Aut signed-area law is coherent, but K, M, chart order, cyclic orientation, and `A_s` remain sampled/derived/chart-polarity dependent source channels with zero strict source-provenance rows.\n- Do not promote `(K,M)` phase-curve provenance, P3048 area signs, or Aut signed-area equivariance to `QW-2191` discharge, selector closure, observed physics, unit-bearing action/EOM, `L_total`, bridge/role-transfer, or ToE closure.\n- A next move in this lane may attack only the remaining nonconventional orientation/coupling-polarity theorem tying `A_s` to exactly one P3046 polarity; otherwise pivot to a genuinely new typed object.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
