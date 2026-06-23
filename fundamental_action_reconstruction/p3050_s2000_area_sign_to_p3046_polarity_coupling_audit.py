#!/usr/bin/env python3
"""P3050/S2000: area-sign to P3046 coupling-polarity audit.

P3049 left exactly one admissible follow-up in the P3048 lane: test whether the
P3048 signed area A_s supplies a nonconventional orientation/coupling-polarity
theorem tying the area sign to exactly one P3046 polarity.  P3050 attacks only
that premise.

The new object is the finite area-sign torsor for nonzero P3048 rows, crossed
with the P3046 coupling-polarity torsor {+,-}.  Inversion sends A_s to -A_s and
exchanges the two coupling-polarity rows.  The computation proves the precise
obstruction: every nonzero area row admits two equivariant sign-to-polarity
maps, not one.  The lag-6 neutral row has zero area and cannot select either
polarity.  Thus A_s has the right odd representation type, but no strict law in
current artifacts chooses the polarity convention.
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
from p3046_s1996_memory_lag_torsor_coupling_polarity_audit import INVERSION_UNITS, UNITS
from p3048_s1998_memory_phase_area_odd_source_candidate import triangle_area
from p3049_s1999_memory_phase_curve_provenance_obstruction import OUT as P3049

OUT = GEN / "p3050_s2000_area_sign_to_p3046_polarity_coupling_audit.json"
MD = GEN / "p3050_s2000_area_sign_to_p3046_polarity_coupling_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"
LAGS = list(range(1, N // 2 + 1))
AREA_SIGNS = [-1, 1]
POLARITIES = ["-", "+"]
TOL = 1e-12


def sign_of(value: float) -> int:
    if value > TOL:
        return 1
    if value < -TOL:
        return -1
    return 0


def sign_action(unit: int, sign: int) -> int:
    if sign == 0:
        return 0
    return -sign if unit in INVERSION_UNITS else sign


def polarity_action(unit: int, polarity: str) -> str:
    if unit not in INVERSION_UNITS:
        return polarity
    return "+" if polarity == "-" else "-"


def candidate_maps() -> list[dict[str, Any]]:
    return [
        {"map_name": "positive_area_selects_positive_p3046_polarity", "map": {-1: "-", 1: "+"}},
        {"map_name": "positive_area_selects_negative_p3046_polarity", "map": {-1: "+", 1: "-"}},
    ]


def is_equivariant(mapping: dict[int, str]) -> bool:
    for unit in UNITS:
        for sign in AREA_SIGNS:
            lhs = mapping[sign_action(unit, sign)]
            rhs = polarity_action(unit, mapping[sign])
            if lhs != rhs:
                return False
    return True


def area_sign_rows() -> list[dict[str, Any]]:
    k = kernel_vector()
    m = memory_viscosity_trace(k)
    rows = []
    for lag in LAGS:
        value = triangle_area(k, m, lag)
        inv_value = triangle_area(k, m, -lag)
        sign = sign_of(value)
        rows.append({
            "lag": lag,
            "area": round(value, 15),
            "inverted_area": round(inv_value, 15),
            "area_sign": sign,
            "nonzero_area_sign": sign != 0,
            "inversion_odd_verified": abs(value + inv_value) < TOL,
            "can_enter_polarity_torsor": sign != 0,
        })
    return rows


def coupling_rows(area_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    rows = []
    for area_row in area_rows:
        if not area_row["can_enter_polarity_torsor"]:
            rows.append({
                "lag": area_row["lag"],
                "map_name": "neutral_area_no_polarity_map",
                "area_sign": 0,
                "aut_equivariant": False,
                "polarity_selected_by_current_artifacts": False,
                "accepted_as_orientation_coupling_theorem": False,
                "failure": "zero area sign cannot choose a P3046 coupling polarity",
            })
            continue
        for candidate in candidate_maps():
            mapping = candidate["map"]
            rows.append({
                "lag": area_row["lag"],
                "map_name": candidate["map_name"],
                "area_sign": area_row["area_sign"],
                "selected_polarity_for_this_area_sign": mapping[area_row["area_sign"]],
                "positive_area_maps_to": mapping[1],
                "negative_area_maps_to": mapping[-1],
                "aut_equivariant": is_equivariant(mapping),
                "polarity_selected_by_current_artifacts": False,
                "accepted_as_orientation_coupling_theorem": False,
                "failure": "equivariant map exists, but the opposite equivariant map exists too; no strict theorem selects this polarity convention",
            })
    return rows


def acceptance_rows(area_rows: list[dict[str, Any]], couplings: list[dict[str, Any]]) -> list[dict[str, Any]]:
    nonzero_lags = [r for r in area_rows if r["nonzero_area_sign"]]
    nonzero_couplings = [r for r in couplings if r["area_sign"] != 0]
    return [
        {"criterion": "concrete_area_sign_torsor_present", "satisfied": bool(nonzero_lags), "detail": "P3048 supplies nonzero signed area rows"},
        {"criterion": "inversion_odd_area_sign", "satisfied": all(r["inversion_odd_verified"] for r in area_rows), "detail": "A_s changes sign under lag inversion"},
        {"criterion": "finite_equivariant_area_to_polarity_maps", "satisfied": all(r["aut_equivariant"] for r in nonzero_couplings), "detail": "both sign-to-polarity maps are Aut-equivariant on nonzero rows"},
        {"criterion": "unique_p3046_polarity_selected", "satisfied": False, "detail": "each nonzero area sign admits two opposite equivariant couplings"},
        {"criterion": "zero_lag_neutral_row_eliminated_by_theorem", "satisfied": False, "detail": "lag 6 remains a neutral A_s=0 row rather than a strict exclusion theorem"},
        {"criterion": "nonconventional_orientation_law", "satisfied": False, "detail": "choosing positive-area -> positive-polarity is a convention unless a strict orientation law is exported"},
        {"criterion": "selector_readout_installation", "satisfied": False, "detail": "no QW-2191/readout row is installed by the conditional coupling pair"},
    ]


def build_matrix() -> dict[str, Any]:
    read_json(P3049)
    areas = area_sign_rows()
    couplings = coupling_rows(areas)
    acceptance = acceptance_rows(areas, couplings)
    obligations = [
        {"obligation": "p3049_remaining_coupling_premise_targeted", "satisfied": True, "detail": "P3050 attacks only the A_s -> P3046 polarity theorem left by P3049"},
        {"obligation": "area_sign_torsor_constructed", "satisfied": any(r["nonzero_area_sign"] for r in areas), "detail": "nonzero P3048 area signs form an inversion-odd torsor"},
        {"obligation": "all_nonzero_equivariant_maps_enumerated", "satisfied": sum(1 for r in couplings if r["area_sign"] != 0 and r["aut_equivariant"]) == 2 * sum(1 for r in areas if r["nonzero_area_sign"]), "detail": "two equivariant maps are found for each nonzero area-sign row"},
        {"obligation": "unique_p3046_coupling_polarity", "satisfied": False, "detail": "the two maps are opposite polarity conventions and current artifacts select neither"},
        {"obligation": "nonconventional_orientation_theorem", "satisfied": False, "detail": "no theorem makes positive area the strict orientation branch"},
        {"obligation": "selector_readout_or_ltotal_installation", "satisfied": False, "detail": "no selector/readout, action/EOM, L_total, bridge, role transfer, or ToE export follows"},
    ]
    return {
        "object": "AreaSignToP3046Polarity_CouplingAudit",
        "area_sign_rows": areas,
        "coupling_rows": couplings,
        "source_acceptance_rows": acceptance,
        "proof_obligations": obligations,
        "finite_certificate": {
            "area_sign_rows": len(areas),
            "nonzero_area_sign_rows": sum(1 for r in areas if r["nonzero_area_sign"]),
            "neutral_area_rows": sum(1 for r in areas if not r["nonzero_area_sign"]),
            "coupling_rows": len(couplings),
            "nonzero_candidate_coupling_rows": sum(1 for r in couplings if r["area_sign"] != 0),
            "aut_equivariant_nonzero_coupling_rows": sum(1 for r in couplings if r["area_sign"] != 0 and r["aut_equivariant"]),
            "polarity_selected_rows": sum(1 for r in couplings if r["polarity_selected_by_current_artifacts"]),
            "accepted_orientation_coupling_rows": sum(1 for r in couplings if r["accepted_as_orientation_coupling_theorem"]),
            "source_acceptance_criteria": len(acceptance),
            "satisfied_source_acceptance_criteria": sum(1 for r in acceptance if r["satisfied"]),
            "proof_obligations": len(obligations),
            "satisfied_proof_obligations": sum(1 for r in obligations if r["satisfied"]),
            "p3046_coupling_polarity_selected": False,
        },
    }


def build_payload() -> dict[str, Any]:
    matrix = build_matrix()
    return {
        "status": "P3050_AREA_SIGN_TO_P3046_POLARITY_COUPLING_AUDIT_BOUNDED_NO_EXPORT",
        "input_hashes": {"P3049": hashlib.sha256(P3049.read_bytes()).hexdigest() if P3049.exists() else None},
        "constructed_theoretical_objects": matrix,
        "finite_certificate": matrix["finite_certificate"],
        "decision": {
            "bounded_no_go": "P3050 constructs the exact A_s-to-P3046 coupling-polarity audit left by P3049.  Nonzero area signs have the right inversion-odd type and every nonzero row admits Aut-equivariant maps to the P3046 polarity torsor, but there are always two opposite maps.  The neutral lag-6 row supplies no sign.  Therefore no unique nonconventional orientation/coupling-polarity theorem is exported.",
            "negative_export_flags": {k: False for k in ["p3046_coupling_polarity_selected", "nonconventional_orientation_theorem_exported", "selector_readout_coupling_exported", "qw2191_discharged", "strict_selector_closure_exported", "observed_physics_exported", "unit_bearing_action_eom_hamiltonian_exported", "ltotal_exported", "bridge_closure_exported", "role_transfer_exported", "toe_closure_exported"]},
            "next_honest_step": "Do not replay A_s sign-to-polarity maps as selector closure.  This lane is bounded unless a genuinely new strict orientation-law/source object selects one of the two maps; otherwise pivot to a different typed object or preserve the P3048-P3050 no-export certificate.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    c = payload["finite_certificate"]
    lines = [
        "# P3050/S2000 area-sign to P3046 coupling-polarity audit", "",
        f"Status: `{payload['status']}`", "", "## Finite certificate",
        f"- area-sign rows: `{c['area_sign_rows']}`",
        f"- nonzero area-sign rows: `{c['nonzero_area_sign_rows']}`",
        f"- neutral area rows: `{c['neutral_area_rows']}`",
        f"- coupling rows: `{c['coupling_rows']}`",
        f"- nonzero candidate coupling rows: `{c['nonzero_candidate_coupling_rows']}`",
        f"- Aut-equivariant nonzero coupling rows: `{c['aut_equivariant_nonzero_coupling_rows']}`",
        f"- polarity-selected rows: `{c['polarity_selected_rows']}`",
        f"- accepted orientation-coupling rows: `{c['accepted_orientation_coupling_rows']}`",
        f"- source acceptance criteria: `{c['satisfied_source_acceptance_criteria']}/{c['source_acceptance_criteria']}`",
        f"- satisfied proof obligations: `{c['satisfied_proof_obligations']}/{c['proof_obligations']}`",
        f"- P3046 coupling polarity selected: `{c['p3046_coupling_polarity_selected']}`", "",
        "## Decision", payload["decision"]["bounded_no_go"], "", "## Recommendation", payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3050/S2000 area-sign to P3046 coupling-polarity audit", "## P3050/S2000 area-sign to P3046 coupling-polarity audit\n\n`P3050/S2000` attacks the remaining P3049 premise: a nonconventional orientation/coupling-polarity theorem tying the P3048 signed area `A_s` to exactly one P3046 polarity.  The finite audit constructs the area-sign torsor and enumerates sign-to-polarity maps.  Every nonzero area row has two Aut-equivariant maps to the P3046 polarity torsor, while the neutral lag-6 row has no sign.  Thus the representation type is real, but no unique polarity-selection theorem, selector/readout installation, unit-bearing action/EOM, `L_total`, bridge/role transfer, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3050/S2000 area-sign polarity `L_total` guard", "## P3050/S2000 area-sign polarity `L_total` guard\n\n`P3050/S2000` adds no physical `L_total` term.  The `A_s` sign can be conditionally coupled to P3046 polarity in two opposite equivariant ways, and the neutral lag-6 row has no sign; no unique nonconventional orientation law or unit-bearing variational/action/EOM input is exported.\n")
    append_once(AGENTS, "Current area-sign to P3046 coupling-polarity guardrail (P3050/S2000, 2026-06-23)", "## Current area-sign to P3046 coupling-polarity guardrail (P3050/S2000, 2026-06-23)\n\n- P3050 attacks exactly one P3049 remaining premise: a nonconventional orientation/coupling-polarity theorem tying `A_s` to exactly one P3046 polarity.\n- Each nonzero P3048 area-sign row admits two Aut-equivariant sign-to-polarity maps, and the neutral lag-6 row has no sign; no strict theorem selects one polarity convention.\n- Do not promote `A_s` sign-to-polarity maps, positive-area conventions, or lag-area winners to `QW-2191` discharge, selector closure, observed physics, unit-bearing action/EOM, `L_total`, bridge/role-transfer, or ToE closure.\n- This lane is bounded unless a genuinely new strict orientation-law/source object selects one of the two maps; otherwise pivot to a different typed object or preserve the P3048-P3050 no-export certificate.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
