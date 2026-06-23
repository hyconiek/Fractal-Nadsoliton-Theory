#!/usr/bin/env python3
"""P3052/S2002: phase-curve winding stability-margin certificate.

P3051 found a real winding +1 topological hint but left it receiver-level.  P3052
makes the next honest step more proof-like: it asks whether the winding hint is
fragile or whether the closed (K_i,M_i) curve has a positive finite stability
margin.  This is not a selector closure attempt; it is a robustness certificate
for the P3051 object.

The constructed object is a centroid-to-edge clearance matrix plus deterministic
Fourier-mode perturbation witness table.  The minimum centroid-to-edge clearance
is positive, and all tested perturbations with amplitudes up to 0.8 of that
clearance preserve the winding.  This strengthens the finite topological hint,
but it still does not export a strict source theorem for the phase-curve
orientation, P3046 polarity, selector/readout installation, or L_total.
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
from p3051_s2001_phase_curve_winding_index_selector_candidate import OUT as P3051, relabel, winding_index

OUT = GEN / "p3052_s2002_phase_curve_winding_stability_margin_certificate.json"
MD = GEN / "p3052_s2002_phase_curve_winding_stability_margin_certificate.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"
FRACTIONS = [0.1, 0.25, 0.4, 0.8]
MODES = [1, 2, 3, 4, 5]
SIGNS = [-1, 1]


def segment_distance(px: float, py: float, ax: float, ay: float, bx: float, by: float) -> float:
    dx = bx - ax
    dy = by - ay
    length_sq = dx * dx + dy * dy
    if length_sq == 0:
        return math.hypot(px - ax, py - ay)
    t = max(0.0, min(1.0, ((px - ax) * dx + (py - ay) * dy) / length_sq))
    qx = ax + t * dx
    qy = ay + t * dy
    return math.hypot(px - qx, py - qy)


def phase_curve() -> tuple[list[float], list[float]]:
    k = kernel_vector()
    return k, memory_viscosity_trace(k)


def clearance_rows() -> list[dict[str, Any]]:
    xs, ys = phase_curve()
    cx = sum(xs) / N
    cy = sum(ys) / N
    rows = []
    for i in range(N):
        j = (i + 1) % N
        distance = segment_distance(cx, cy, xs[i], ys[i], xs[j], ys[j])
        rows.append({
            "edge": [i, j],
            "centroid": [round(cx, 15), round(cy, 15)],
            "distance_to_centroid": round(distance, 15),
            "positive_clearance": distance > 0,
        })
    return rows


def perturb(xs: list[float], ys: list[float], amplitude: float, mode: int) -> tuple[list[float], list[float]]:
    x2 = []
    y2 = []
    for i, (x, y) in enumerate(zip(xs, ys)):
        theta = 2 * math.pi * mode * i / N
        x2.append(x + amplitude * math.cos(theta))
        y2.append(y + amplitude * math.sin(theta))
    return x2, y2


def perturbation_rows(min_clearance: float) -> list[dict[str, Any]]:
    xs, ys = phase_curve()
    base = winding_index(xs, ys)["integer_winding"]
    rows = []
    for fraction in FRACTIONS:
        for mode in MODES:
            for sign in SIGNS:
                amplitude = sign * fraction * min_clearance
                px, py = perturb(xs, ys, amplitude, mode)
                w = winding_index(px, py)
                rows.append({
                    "fraction_of_clearance": fraction,
                    "mode": mode,
                    "amplitude_sign": sign,
                    "amplitude": round(amplitude, 15),
                    "integer_winding": w["integer_winding"],
                    "winding_preserved": w["integer_winding"] == base,
                    "strict_source_exported": False,
                })
    return rows


def aut_rows() -> list[dict[str, Any]]:
    xs, ys = phase_curve()
    base = winding_index(xs, ys)["integer_winding"]
    rows = []
    for unit in UNITS:
        expected = -base if unit in INVERSION_UNITS else base
        rx = relabel(xs, unit, 0)
        ry = relabel(ys, unit, 0)
        w = winding_index(rx, ry)
        rows.append({
            "unit": unit,
            "integer_winding": w["integer_winding"],
            "expected_integer_winding": expected,
            "aut_signed_winding_verified": w["integer_winding"] == expected,
            "orientation_reversing_unit": unit in INVERSION_UNITS,
        })
    return rows


def source_acceptance_rows(base: dict[str, Any], clearance: list[dict[str, Any]], perturbations: list[dict[str, Any]], aut: list[dict[str, Any]]) -> list[dict[str, Any]]:
    min_clearance = min(row["distance_to_centroid"] for row in clearance)
    return [
        {"criterion": "winding_index_reused_from_p3051", "satisfied": base["nonzero_integer_winding"], "detail": "base winding remains +1"},
        {"criterion": "positive_centroid_edge_clearance", "satisfied": min_clearance > 0, "detail": f"minimum clearance is {min_clearance}"},
        {"criterion": "deterministic_perturbation_stability", "satisfied": all(row["winding_preserved"] for row in perturbations), "detail": "all Fourier-mode perturbations preserve winding"},
        {"criterion": "aut_signed_winding_verified", "satisfied": all(row["aut_signed_winding_verified"] for row in aut), "detail": "Aut units preserve or reverse as expected"},
        {"criterion": "strict_winding_source_theorem", "satisfied": False, "detail": "robustness does not make sampled receiver geometry a strict source"},
        {"criterion": "nonconventional_orientation_law", "satisfied": False, "detail": "robust +1 still requires an orientation-source theorem"},
        {"criterion": "p3046_polarity_or_selector_coupling", "satisfied": False, "detail": "no unique P3046 polarity or selector readout is installed"},
        {"criterion": "unit_bearing_variational_installation", "satisfied": False, "detail": "no action/EOM/Hamiltonian row is exported"},
    ]


def build_matrix() -> dict[str, Any]:
    read_json(P3051)
    xs, ys = phase_curve()
    base = winding_index(xs, ys)
    clearance = clearance_rows()
    min_clearance = min(row["distance_to_centroid"] for row in clearance)
    perturbations = perturbation_rows(min_clearance)
    aut = aut_rows()
    acceptance = source_acceptance_rows(base, clearance, perturbations, aut)
    obligations = [
        {"obligation": "p3051_winding_stability_targeted", "satisfied": True, "detail": "P3052 strengthens only the robustness status of the P3051 winding object"},
        {"obligation": "positive_clearance_certificate", "satisfied": min_clearance > 0, "detail": "centroid has positive distance to every polygon edge"},
        {"obligation": "finite_perturbation_witness_table", "satisfied": all(row["winding_preserved"] for row in perturbations), "detail": "all deterministic perturbation rows preserve winding"},
        {"obligation": "aut_signed_winding_boundary", "satisfied": all(row["aut_signed_winding_verified"] for row in aut), "detail": "Aut preserving/reversing boundary is explicit"},
        {"obligation": "strict_orientation_source_theorem", "satisfied": False, "detail": "no strict theorem sources the winding orientation sign"},
        {"obligation": "selector_ltotal_bridge_toe_installation", "satisfied": False, "detail": "no selector closure, L_total, bridge, role transfer, or ToE export follows"},
    ]
    return {
        "object": "PhaseCurveWinding_StabilityMarginCertificate",
        "base_winding_row": base,
        "clearance_rows": clearance,
        "perturbation_rows": perturbations,
        "aut_signed_winding_rows": aut,
        "source_acceptance_rows": acceptance,
        "proof_obligations": obligations,
        "finite_certificate": {
            "base_integer_winding": base["integer_winding"],
            "clearance_rows": len(clearance),
            "positive_clearance_rows": sum(1 for row in clearance if row["positive_clearance"]),
            "minimum_centroid_edge_clearance": min_clearance,
            "perturbation_rows": len(perturbations),
            "winding_preserved_perturbation_rows": sum(1 for row in perturbations if row["winding_preserved"]),
            "aut_signed_winding_rows": len(aut),
            "aut_signed_winding_verified_rows": sum(1 for row in aut if row["aut_signed_winding_verified"]),
            "source_acceptance_criteria": len(acceptance),
            "satisfied_source_acceptance_criteria": sum(1 for row in acceptance if row["satisfied"]),
            "proof_obligations": len(obligations),
            "satisfied_proof_obligations": sum(1 for row in obligations if row["satisfied"]),
            "strict_winding_source_theorem_exported": False,
        },
    }


def build_payload() -> dict[str, Any]:
    matrix = build_matrix()
    return {
        "status": "P3052_PHASE_CURVE_WINDING_STABILITY_MARGIN_CERTIFICATE_BOUNDED_NO_EXPORT",
        "input_hashes": {"P3051": hashlib.sha256(P3051.read_bytes()).hexdigest() if P3051.exists() else None},
        "constructed_theoretical_objects": matrix,
        "finite_certificate": matrix["finite_certificate"],
        "decision": {
            "bounded_no_go": "P3052 strengthens P3051 from a single winding witness to a finite stability-margin certificate: the centroid has positive clearance from every edge, all deterministic Fourier perturbation rows preserve winding, and Aut signed-winding behavior is explicit.  This is stronger evidence for a robust receiver-level orientation hint, but robustness is not a strict source theorem, selector/readout coupling, or L_total installation.",
            "negative_export_flags": {k: False for k in ["strict_winding_source_theorem_exported", "nonconventional_orientation_theorem_exported", "p3046_coupling_polarity_selected", "selector_readout_coupling_exported", "qw2191_discharged", "strict_selector_closure_exported", "observed_physics_exported", "unit_bearing_action_eom_hamiltonian_exported", "ltotal_exported", "bridge_closure_exported", "role_transfer_exported", "toe_closure_exported"]},
            "next_honest_step": "Do not promote robust winding to selector closure.  The next proof-grade move must supply a strict source theorem explaining why this robust winding orientation is the nadsoliton selector sign, or pivot to an independent typed object outside sampled K/M geometry.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    c = payload["finite_certificate"]
    lines = [
        "# P3052/S2002 phase-curve winding stability-margin certificate", "",
        f"Status: `{payload['status']}`", "", "## Finite certificate",
        f"- base integer winding: `{c['base_integer_winding']}`",
        f"- clearance rows: `{c['clearance_rows']}`",
        f"- positive clearance rows: `{c['positive_clearance_rows']}`",
        f"- minimum centroid-edge clearance: `{c['minimum_centroid_edge_clearance']}`",
        f"- perturbation rows: `{c['perturbation_rows']}`",
        f"- winding-preserved perturbation rows: `{c['winding_preserved_perturbation_rows']}`",
        f"- Aut signed-winding rows: `{c['aut_signed_winding_rows']}`",
        f"- Aut signed-winding verified rows: `{c['aut_signed_winding_verified_rows']}`",
        f"- source acceptance criteria: `{c['satisfied_source_acceptance_criteria']}/{c['source_acceptance_criteria']}`",
        f"- satisfied proof obligations: `{c['satisfied_proof_obligations']}/{c['proof_obligations']}`",
        f"- strict winding source theorem exported: `{c['strict_winding_source_theorem_exported']}`", "",
        "## Decision", payload["decision"]["bounded_no_go"], "", "## Recommendation", payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3052/S2002 phase-curve winding stability-margin certificate", "## P3052/S2002 phase-curve winding stability-margin certificate\n\n`P3052/S2002` strengthens the P3051 winding hint with a finite stability-margin certificate.  The centroid has positive distance to every closed phase-curve edge, the minimum centroid-edge clearance is positive, all deterministic Fourier-mode perturbations up to the audited clearance fractions preserve winding, and Aut signed-winding rows verify the preserving/reversing boundary.  This is robust receiver-level evidence only: no strict winding-source theorem, nonconventional orientation theorem, P3046 polarity selection, selector/readout installation, unit-bearing action/EOM, `L_total`, bridge/role transfer, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3052/S2002 winding stability-margin `L_total` guard", "## P3052/S2002 winding stability-margin `L_total` guard\n\n`P3052/S2002` adds no physical `L_total` term.  Positive clearance and perturbation-stable winding strengthen the receiver diagnostic, but do not supply a strict orientation source, P3046 polarity coupling, or unit-bearing variational/action/EOM input.\n")
    append_once(AGENTS, "Current phase-curve winding stability-margin guardrail (P3052/S2002, 2026-06-23)", "## Current phase-curve winding stability-margin guardrail (P3052/S2002, 2026-06-23)\n\n- P3052 strengthens P3051 by adding a centroid-edge clearance matrix and deterministic perturbation witness table for the `(K_i,M_i)` winding `+1` diagnostic.\n- The audited winding is robust under the finite perturbation table and has explicit Aut signed-winding behavior, but this remains receiver-level geometry without strict orientation-source provenance or selector/readout coupling.\n- Do not promote robust winding, positive clearance, perturbation stability, or Aut signed-winding rows to `QW-2191` discharge, selector closure, observed physics, unit-bearing action/EOM, `L_total`, bridge/role-transfer, or ToE closure.\n- A next move requires a strict source theorem explaining why this robust winding orientation is the nadsoliton selector sign, or an independent typed object outside sampled K/M geometry.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
