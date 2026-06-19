#!/usr/bin/env python3
"""P2941/S1891: P2938 carrier extremal-selector polarity obstruction.

P2940 quantified the selector burden for the exact P2938 carrier [1,2,2,2,2].
P2941 constructs the next theorem-object candidate rather than another carrier:
use the carrier's own finite values on each U(12)-orbit to define extremal
representative selectors (min-value and max-value).

The finite calculation is deliberately two-sided.  It shows that the min-value
rule gives a unique representative on every audited orbit, while the max-value
rule has a tie on the unit orbit.  However, the min/max polarity is not sourced
by a strict theorem, and the selected representatives are not U(12)-compatible
as global sections.  Therefore the construction is a useful conditional
selector skeleton, not strict provenance for P2938 and not a damping source.
"""
from __future__ import annotations

import hashlib
import json
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p2940_s1890_p2938_carrier_aut_orbit_selector_burden import MODULUS, OUT as P2940, UNITS, carrier_value, orbit_rows

OUT = GEN / "p2941_s1891_p2938_carrier_extremal_selector_polarity_obstruction.json"
MD = GEN / "p2941_s1891_p2938_carrier_extremal_selector_polarity_obstruction.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"


def selector_rows() -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for orbit_row in orbit_rows():
        orbit = orbit_row["orbit"]
        values = {int(k): v for k, v in orbit_row["carrier_values"].items()}
        for polarity in ["min", "max"]:
            extremal = min(values.values()) if polarity == "min" else max(values.values())
            selected = [node for node in orbit if values[node] == extremal]
            rows.append({
                "orbit": orbit,
                "polarity": polarity,
                "extremal_value": extremal,
                "selected_nodes": selected,
                "unique_selector": len(selected) == 1,
                "selector_node": selected[0] if len(selected) == 1 else None,
            })
    return rows


def compatibility_rows(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    comp: list[dict[str, Any]] = []
    for row in rows:
        if not row["unique_selector"]:
            continue
        s = row["selector_node"]
        orbit = row["orbit"]
        for u in UNITS:
            moved = (u * s) % MODULUS
            comp.append({
                "orbit": orbit,
                "polarity": row["polarity"],
                "unit": u,
                "selector_node": s,
                "moved_selector_node": moved,
                "remains_selected": moved == s,
                "carrier_value_at_moved_node": carrier_value(moved),
            })
    return comp


def acceptance_rows(sel: list[dict[str, Any]], comp: list[dict[str, Any]]) -> list[dict[str, Any]]:
    min_rows = [row for row in sel if row["polarity"] == "min"]
    max_rows = [row for row in sel if row["polarity"] == "max"]
    return [
        {
            "criterion": "carrier_extremal_selector_skeleton_constructed",
            "satisfied": True,
            "evidence": "min-value and max-value selectors are computed on every P2940 orbit",
        },
        {
            "criterion": "min_polarity_unique_on_all_orbits",
            "satisfied": all(row["unique_selector"] for row in min_rows),
            "evidence": f"{sum(row['unique_selector'] for row in min_rows)} of {len(min_rows)} min rows are unique",
        },
        {
            "criterion": "max_polarity_unique_on_all_orbits",
            "satisfied": all(row["unique_selector"] for row in max_rows),
            "evidence": f"{sum(row['unique_selector'] for row in max_rows)} of {len(max_rows)} max rows are unique; unit orbit ties remain",
        },
        {
            "criterion": "global_section_U12_compatible",
            "satisfied": all(row["remains_selected"] for row in comp),
            "evidence": f"{sum(not row['remains_selected'] for row in comp)} U(12)-motion defects occur on unique selector rows",
        },
        {
            "criterion": "strict_polarity_source_theorem_exported",
            "satisfied": False,
            "evidence": "current artifacts do not source the choice of min rather than max as a strict nadsoliton theorem",
        },
        {
            "criterion": "delta_eta_and_beta_eta_coupling_exported",
            "satisfied": False,
            "evidence": "the extremal selector skeleton is not coupled to the damping delta/eta or beta/eta source obligations",
        },
    ]


def build_payload(_: dict[str, Any]) -> dict[str, Any]:
    sel = selector_rows()
    comp = compatibility_rows(sel)
    criteria = acceptance_rows(sel, comp)
    min_unique = all(row["unique_selector"] for row in sel if row["polarity"] == "min")
    max_ties = [row for row in sel if row["polarity"] == "max" and not row["unique_selector"]]
    motion_defects = [row for row in comp if not row["remains_selected"]]
    return {
        "status": "P2941_P2938_CARRIER_EXTREMAL_SELECTOR_POLARITY_OBSTRUCTION_NO_STRICT_PROVENANCE",
        "input_hashes": {"P2940": hashlib.sha256(P2940.read_bytes()).hexdigest() if P2940.exists() else None},
        "constructed_theoretical_objects": {
            "candidate_object": "P2938_Carrier_Value_Extremal_Orbit_Selector_Skeleton",
            "selector_rows": sel,
            "u12_motion_compatibility_rows": comp,
            "acceptance_rows": criteria,
        },
        "selector_certificate": {
            "orbit_selector_rows": len(sel),
            "min_polarity_unique_on_all_orbits": min_unique,
            "max_polarity_tie_count": len(max_ties),
            "u12_motion_rows": len(comp),
            "u12_motion_defect_count": len(motion_defects),
            "strict_polarity_source_theorem_exported": False,
            "strict_provenance_theorem_exported": False,
            "accepted_strict_source": False,
        },
        "decision": {
            "positive_witnesses": {
                "finite_extremal_selector_skeleton_exists": True,
                "min_polarity_repairs_orbit_representative_ties_conditionally": min_unique,
            },
            "negative_export_flags": {
                "strict_selector_source_exported": False,
                "strict_aut_breaking_prime_coordinate_source_law_exported": False,
                "strict_prime_log_value_source_exported": False,
                "strict_delta_eta_source_exported": False,
                "strict_beta_eta_coupling_theorem_exported": False,
                "strict_damping_beta_eta_source_packet_exported": False,
                "nonproxy_ltotal_exported": False,
                "eom_hamiltonian_exported": False,
                "bridge_closure_exported": False,
                "role_transfer_exported": False,
                "toe_closure_exported": False,
            },
            "reason": "P2941 constructs a conditional extremal-selector skeleton from the exact P2938 carrier.  The min-value polarity is finite and unique on all audited orbits, but choosing min over max is an unsourced polarity convention, max has a unit-orbit tie, and the resulting sections have U(12)-motion defects.  Thus no strict provenance or damping coupling theorem is exported.",
            "next_honest_step": "A next admissible move must supply a strict polarity/source theorem choosing the min-value extremal selector and coupling it to delta/eta plus beta/eta, or pivot to a genuinely new typed object outside the P2938 extremal-selector lane.  Otherwise preserve the P2929-P2941 no-new-live-frontier boundary.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    cert = payload["selector_certificate"]
    lines = [
        "# P2941/S1891 P2938 carrier extremal-selector polarity obstruction",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Selector certificate",
        f"- orbit selector rows: `{cert['orbit_selector_rows']}`",
        f"- min polarity unique on all orbits: `{cert['min_polarity_unique_on_all_orbits']}`",
        f"- max polarity tie count: `{cert['max_polarity_tie_count']}`",
        f"- U(12)-motion rows: `{cert['u12_motion_rows']}`",
        f"- U(12)-motion defect count: `{cert['u12_motion_defect_count']}`",
        f"- strict polarity source theorem exported: `{cert['strict_polarity_source_theorem_exported']}`",
        f"- accepted strict source: `{cert['accepted_strict_source']}`",
        "",
        "## Boundary",
        payload["decision"]["reason"],
        "",
        "## Recommendation",
        payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload(read_json(P2940))
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2941/S1891 P2938 extremal-selector polarity obstruction", "## P2941/S1891 P2938 extremal-selector polarity obstruction\n\n`P2941/S1891` constructs the next theorem-object candidate for P2938 provenance: an extremal orbit-selector skeleton using carrier values on the P2940 `U(12)` orbits.  The min-value polarity is unique on all audited orbits, but min-versus-max remains an unsourced polarity convention, max has a unit-orbit tie, and the unique selected sections have `U(12)`-motion defects.  Therefore the skeleton does not export strict provenance, delta/eta source, beta/eta coupling, strict `L_p`, nonproxy `L_total`, EOM, Hamiltonian, bridge closure, role transfer, or ToE closure.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2941/S1891 extremal-selector `L_total` guard", "## P2941/S1891 extremal-selector `L_total` guard\n\n`P2941/S1891` gives a conditional min-extremal selector skeleton for the P2938 carrier, but no strict polarity/source theorem and no damping coupling are exported.  The skeleton cannot enter nonproxy `L_total`, EOM, Hamiltonian, bridge closure, role transfer, or ToE.\n")
    append_once(AGENTS, "Current P2938 extremal-selector polarity guardrail (P2941/S1891, 2026-06-19)", "## Current P2938 extremal-selector polarity guardrail (P2941/S1891, 2026-06-19)\n\n- P2941 constructs an extremal orbit-selector skeleton from the exact P2938 carrier values on the P2940 `U(12)` orbits.\n- The min-value polarity is finite and unique on all audited orbits, but choosing min rather than max is not strictly sourced; max has a unit-orbit tie and the selected sections have `U(12)`-motion defects.\n- Do not promote this skeleton to selector closure, strict provenance, strict `L_p`, delta/eta source, beta/eta coupling, nonproxy `L_total`, bridge closure, role transfer, or ToE.\n- A next admissible move must export a strict polarity/source theorem plus damping coupling for this exact skeleton, or pivot to a genuinely new typed object outside the P2938 extremal-selector lane; otherwise preserve the P2929-P2941 no-new-live-frontier boundary.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
