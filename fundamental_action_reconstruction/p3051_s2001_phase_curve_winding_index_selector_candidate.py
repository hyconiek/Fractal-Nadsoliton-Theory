#!/usr/bin/env python3
"""P3051/S2001: phase-curve winding-index selector candidate.

After P3050, the P3048 area-sign lane is bounded unless a genuinely new strict
orientation-law/source object appears.  P3051 pivots to a different typed object
on the same finite receiver data: the global winding/turning index of the closed
(K_i,M_i) phase curve around its centroid.

This is not another area-score or lag-polarity map.  It computes the total
unwrapped angular turn of the 12-point phase curve.  The finite result is a real
new topological hint: the base curve has winding +1; Aut units 1 and 5 preserve
it, units 7 and 11 reverse it; all translations preserve the value.  The honest
obstruction is unchanged at the source level: the centroided phase curve and its
orientation are still sampled/derived receiver data, not a strict nadsoliton
orientation law or selector/readout installation.
"""
from __future__ import annotations

import hashlib, json, math
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p3023_s1973_kernel_dissipation_time_order_candidate_obstruction import N
from p3038_s1988_viscous_retarded_chiral_selector_candidate_obstruction import memory_viscosity_trace
from p3044_s1994_memory_lag_commutator_source_candidate import kernel_vector
from p3046_s1996_memory_lag_torsor_coupling_polarity_audit import INVERSION_UNITS, UNITS
from p3050_s2000_area_sign_to_p3046_polarity_coupling_audit import OUT as P3050

OUT = GEN / "p3051_s2001_phase_curve_winding_index_selector_candidate.json"
MD = GEN / "p3051_s2001_phase_curve_winding_index_selector_candidate.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"
TOL = 1e-9


def relabel(values: list[float], unit: int = 1, shift: int = 0) -> list[float]:
    return [values[(unit * i + shift) % N] for i in range(N)]


def wrap_delta(delta: float) -> float:
    return (delta + math.pi) % (2 * math.pi) - math.pi


def winding_index(xs: list[float], ys: list[float]) -> dict[str, Any]:
    cx = sum(xs) / len(xs)
    cy = sum(ys) / len(ys)
    angles = [math.atan2(y - cy, x - cx) for x, y in zip(xs, ys)]
    deltas = [wrap_delta(angles[(i + 1) % len(angles)] - angles[i]) for i in range(len(angles))]
    total_turn = sum(deltas)
    winding = total_turn / (2 * math.pi)
    rounded = round(winding)
    return {
        "centroid": [round(cx, 15), round(cy, 15)],
        "total_turn_radians": round(total_turn, 15),
        "winding": round(winding, 15),
        "integer_winding": int(rounded) if abs(winding - rounded) < TOL else None,
        "nonzero_integer_winding": abs(winding - rounded) < TOL and rounded != 0,
        "turn_deltas": [round(d, 15) for d in deltas],
    }


def base_winding_row() -> dict[str, Any]:
    k = kernel_vector()
    m = memory_viscosity_trace(k)
    row = winding_index(k, m)
    row.update({
        "row": "base_phase_curve_winding",
        "strict_orientation_source_exported": False,
        "failure": "nonzero winding is computed from sampled K and derived M receiver curve, not strict source data",
    })
    return row


def aut_translation_rows() -> list[dict[str, Any]]:
    k = kernel_vector()
    m = memory_viscosity_trace(k)
    rows = []
    base = winding_index(k, m)["integer_winding"]
    for unit in UNITS:
        expected = -base if unit in INVERSION_UNITS else base
        for shift in range(N):
            kk = relabel(k, unit, shift)
            mm = relabel(m, unit, shift)
            w = winding_index(kk, mm)
            rows.append({
                "unit": unit,
                "shift": shift,
                "integer_winding": w["integer_winding"],
                "expected_integer_winding": expected,
                "translation_stable": w["integer_winding"] == expected,
                "orientation_reversing_unit": unit in INVERSION_UNITS,
                "strict_source_exported": False,
            })
    return rows


def source_acceptance_rows(base: dict[str, Any], aut_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return [
        {"criterion": "global_winding_index_constructed", "satisfied": True, "detail": "total unwrapped turn of the closed (K,M) curve is explicitly computed"},
        {"criterion": "nonzero_integer_winding", "satisfied": base["nonzero_integer_winding"], "detail": "base winding is +1"},
        {"criterion": "translation_stability", "satisfied": all(r["translation_stable"] for r in aut_rows), "detail": "all 48 unit/shift relabel rows match the expected signed winding"},
        {"criterion": "inversion_odd_orientation_signal", "satisfied": any(r["orientation_reversing_unit"] and r["integer_winding"] == -base["integer_winding"] for r in aut_rows), "detail": "units 7 and 11 reverse the winding sign"},
        {"criterion": "strict_KM_phase_curve_source", "satisfied": False, "detail": "P3049 already blocks strict source provenance for the sampled/derived phase curve"},
        {"criterion": "nonconventional_winding_orientation_law", "satisfied": False, "detail": "choosing +1 rather than -1 remains a chart-orientation polarity choice"},
        {"criterion": "P3046_or_selector_coupling_selected", "satisfied": False, "detail": "P3050 blocks unique area-sign-to-polarity selection, and winding sign has the same two-branch obstruction"},
        {"criterion": "unit_bearing_action_eom_installation", "satisfied": False, "detail": "no variational/action/EOM row is exported"},
    ]


def build_matrix() -> dict[str, Any]:
    read_json(P3050)
    base = base_winding_row()
    aut_rows = aut_translation_rows()
    acceptance = source_acceptance_rows(base, aut_rows)
    obligations = [
        {"obligation": "new_typed_object_not_area_replay", "satisfied": True, "detail": "uses global winding/turning index rather than A_s area score or lag winner"},
        {"obligation": "finite_nonzero_topological_hint", "satisfied": base["nonzero_integer_winding"], "detail": "base closed phase curve has winding +1"},
        {"obligation": "aut_translation_witness_matrix", "satisfied": all(r["translation_stable"] for r in aut_rows), "detail": "48 unit/shift rows match preserving/reversing expectations"},
        {"obligation": "strict_winding_orientation_source", "satisfied": False, "detail": "no strict source theorem selects the winding orientation sign"},
        {"obligation": "selector_readout_or_p3046_coupling", "satisfied": False, "detail": "no unique P3046 polarity or selector/readout installation follows"},
        {"obligation": "ltotal_bridge_role_toe_installation", "satisfied": False, "detail": "no L_total, bridge, role-transfer, or ToE export follows"},
    ]
    return {
        "object": "PhaseCurveWindingIndex_SelectorCandidateAudit",
        "base_winding_row": base,
        "aut_translation_rows": aut_rows,
        "source_acceptance_rows": acceptance,
        "proof_obligations": obligations,
        "finite_certificate": {
            "base_integer_winding": base["integer_winding"],
            "base_nonzero_integer_winding": base["nonzero_integer_winding"],
            "aut_translation_rows": len(aut_rows),
            "translation_stable_rows": sum(1 for r in aut_rows if r["translation_stable"]),
            "orientation_reversing_rows": sum(1 for r in aut_rows if r["orientation_reversing_unit"]),
            "strict_source_exported_rows": sum(1 for r in aut_rows if r["strict_source_exported"]),
            "source_acceptance_criteria": len(acceptance),
            "satisfied_source_acceptance_criteria": sum(1 for r in acceptance if r["satisfied"]),
            "proof_obligations": len(obligations),
            "satisfied_proof_obligations": sum(1 for r in obligations if r["satisfied"]),
            "strict_winding_selector_source_exported": False,
        },
    }


def build_payload() -> dict[str, Any]:
    matrix = build_matrix()
    return {
        "status": "P3051_PHASE_CURVE_WINDING_INDEX_SELECTOR_CANDIDATE_BOUNDED_NO_EXPORT",
        "input_hashes": {"P3050": hashlib.sha256(P3050.read_bytes()).hexdigest() if P3050.exists() else None},
        "constructed_theoretical_objects": matrix,
        "finite_certificate": matrix["finite_certificate"],
        "decision": {
            "bounded_no_go": "P3051 constructs a genuinely different finite object: the global winding/turning index of the (K,M) phase curve.  The base curve has winding +1, all translations preserve the signed value, and Aut inversion units reverse it.  This is a real topological orientation hint, but it remains receiver-level because K/M provenance, nonconventional orientation law, P3046 coupling, selector/readout installation, and unit-bearing action/EOM are absent.",
            "negative_export_flags": {k: False for k in ["strict_winding_selector_source_exported", "strict_KM_phase_curve_source_exported", "nonconventional_orientation_theorem_exported", "p3046_coupling_polarity_selected", "selector_readout_coupling_exported", "qw2191_discharged", "strict_selector_closure_exported", "observed_physics_exported", "unit_bearing_action_eom_hamiltonian_exported", "ltotal_exported", "bridge_closure_exported", "role_transfer_exported", "toe_closure_exported"]},
            "next_honest_step": "Do not promote winding +1 as selector closure.  A next proof-grade move needs a strict source theorem for the phase-curve orientation/winding sign or an independent typed object outside sampled K/M receiver geometry; otherwise preserve the P3048-P3051 bounded no-export certificate.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    c = payload["finite_certificate"]
    lines = [
        "# P3051/S2001 phase-curve winding-index selector candidate", "",
        f"Status: `{payload['status']}`", "", "## Finite certificate",
        f"- base integer winding: `{c['base_integer_winding']}`",
        f"- base nonzero integer winding: `{c['base_nonzero_integer_winding']}`",
        f"- Aut/translation rows: `{c['aut_translation_rows']}`",
        f"- translation-stable rows: `{c['translation_stable_rows']}`",
        f"- orientation-reversing rows: `{c['orientation_reversing_rows']}`",
        f"- strict source exported rows: `{c['strict_source_exported_rows']}`",
        f"- source acceptance criteria: `{c['satisfied_source_acceptance_criteria']}/{c['source_acceptance_criteria']}`",
        f"- satisfied proof obligations: `{c['satisfied_proof_obligations']}/{c['proof_obligations']}`",
        f"- strict winding selector source exported: `{c['strict_winding_selector_source_exported']}`", "",
        "## Decision", payload["decision"]["bounded_no_go"], "", "## Recommendation", payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3051/S2001 phase-curve winding-index selector candidate", "## P3051/S2001 phase-curve winding-index selector candidate\n\n`P3051/S2001` pivots beyond the P3050 area-sign coupling replay by constructing the global winding/turning index of the closed `(K_i,M_i)` phase curve around its centroid.  The base curve has integer winding `+1`; all `48` Aut/translation rows match the expected preserving/reversing signed winding, with units `7` and `11` reversing the sign.  This is a real topological orientation hint, but it remains sampled/derived receiver geometry: no strict phase-curve source theorem, nonconventional winding-orientation law, P3046 coupling, selector/readout installation, unit-bearing action/EOM, `L_total`, bridge/role transfer, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3051/S2001 phase-curve winding-index `L_total` guard", "## P3051/S2001 phase-curve winding-index `L_total` guard\n\n`P3051/S2001` adds no physical `L_total` term.  The winding `+1` of the sampled `(K,M)` receiver curve is a finite topological diagnostic, but no strict orientation source, P3046 polarity coupling, or unit-bearing variational/action/EOM input is exported.\n")
    append_once(AGENTS, "Current phase-curve winding-index selector candidate guardrail (P3051/S2001, 2026-06-23)", "## Current phase-curve winding-index selector candidate guardrail (P3051/S2001, 2026-06-23)\n\n- P3051 introduces a different typed object after the P3050 area-sign lane: the global winding/turning index of the closed `(K_i,M_i)` phase curve around its centroid.\n- The finite curve has winding `+1`, translations preserve it, and Aut inversion units reverse it, but this remains sampled/derived receiver geometry without strict orientation-source provenance or selector/readout coupling.\n- Do not promote winding `+1`, phase-curve turning direction, or Aut signed-winding rows to `QW-2191` discharge, selector closure, observed physics, unit-bearing action/EOM, `L_total`, bridge/role-transfer, or ToE closure.\n- A next move requires a strict source theorem for the phase-curve orientation/winding sign or an independent typed object outside sampled K/M receiver geometry; otherwise preserve the P3048-P3051 bounded no-export certificate.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
