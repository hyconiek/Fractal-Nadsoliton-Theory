#!/usr/bin/env python3
"""P2942/S1892: extremal-selector value-order orientation source gate.

P2941 built a conditional extremal selector skeleton for the P2938 carrier.
The missing theorem object was not another selector variant, but a strict source
for the polarity/order orientation: why should the carrier order be read as
"min" rather than "max"?

P2942 makes that missing object finite and executable.  It scans affine
value-order orientations q(v)=a*v+b over a bounded rational/integer grid and
checks the induced argmin selector on every P2940 orbit.  Positive slopes
reproduce the P2941 min selector and are unique on every orbit; negative slopes
reverse polarity and inherit the max unit-orbit tie; zero slopes collapse all
value information.  Therefore the only remaining proof-grade gap in this
micro-lane is a strict theorem sourcing the positive orientation sign and then
coupling it to delta/eta and beta/eta.
"""
from __future__ import annotations

import hashlib
import json
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p2940_s1890_p2938_carrier_aut_orbit_selector_burden import orbit_rows
from p2941_s1891_p2938_carrier_extremal_selector_polarity_obstruction import OUT as P2941

OUT = GEN / "p2942_s1892_extremal_selector_value_order_orientation_source_gate.json"
MD = GEN / "p2942_s1892_extremal_selector_value_order_orientation_source_gate.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

SLOPES = list(range(-3, 4))
INTERCEPTS = list(range(-2, 3))


def induced_selector_for_orbit(orbit_row: dict[str, Any], slope: int, intercept: int) -> dict[str, Any]:
    orbit = orbit_row["orbit"]
    values = {int(k): v for k, v in orbit_row["carrier_values"].items()}
    scores = {node: slope * values[node] + intercept for node in orbit}
    minimum = min(scores.values())
    selected = [node for node in orbit if scores[node] == minimum]
    return {
        "orbit": orbit,
        "carrier_values": values,
        "slope": slope,
        "intercept": intercept,
        "orientation_class": "positive_min" if slope > 0 else "negative_max" if slope < 0 else "zero_collapse",
        "scores": scores,
        "selected_nodes": selected,
        "unique_selector": len(selected) == 1,
        "selector_node": selected[0] if len(selected) == 1 else None,
    }


def orientation_scan_rows() -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    orbits = orbit_rows()
    for slope in SLOPES:
        for intercept in INTERCEPTS:
            orbit_results = [induced_selector_for_orbit(row, slope, intercept) for row in orbits]
            rows.append({
                "slope": slope,
                "intercept": intercept,
                "orientation_class": "positive_min" if slope > 0 else "negative_max" if slope < 0 else "zero_collapse",
                "unique_orbit_count": sum(1 for row in orbit_results if row["unique_selector"]),
                "all_orbits_unique": all(row["unique_selector"] for row in orbit_results),
                "tie_orbits": [row["orbit"] for row in orbit_results if not row["unique_selector"]],
                "orbit_results": orbit_results,
            })
    return rows


def source_candidate_rows() -> list[dict[str, Any]]:
    return [
        {
            "candidate": "positive_affine_value_order_q(v)=a*v+b_with_a>0",
            "finite_effect": "reproduces the P2941 min selector and is unique on all audited orbits",
            "accepted_strict_source": False,
            "blocker": "the positive sign of a is not sourced by a strict nadsoliton theorem",
        },
        {
            "candidate": "negative_affine_value_order_q(v)=a*v+b_with_a<0",
            "finite_effect": "reverses to the max selector and inherits the unit-orbit tie",
            "accepted_strict_source": False,
            "blocker": "opposite polarity remains equally formal without an orientation source",
        },
        {
            "candidate": "zero_slope_value_order_q(v)=b",
            "finite_effect": "collapses all carrier value distinctions and ties every non-singleton orbit",
            "accepted_strict_source": False,
            "blocker": "no selector information remains",
        },
        {
            "candidate": "external_order_or_minimization_convention",
            "finite_effect": "can name the positive slope but imports a convention",
            "accepted_strict_source": False,
            "blocker": "not a strict nadsoliton-sourced polarity theorem",
        },
        {
            "candidate": "future_strict_positive_orientation_source_theorem",
            "finite_effect": "would pass this gate only if it fixes a>0 nonconventionally and couples to damping",
            "accepted_strict_source": False,
            "blocker": "schema only; theorem not exported in current artifacts",
        },
    ]


def acceptance_rows(scan: list[dict[str, Any]]) -> list[dict[str, Any]]:
    positive = [row for row in scan if row["orientation_class"] == "positive_min"]
    negative = [row for row in scan if row["orientation_class"] == "negative_max"]
    zero = [row for row in scan if row["orientation_class"] == "zero_collapse"]
    return [
        {
            "criterion": "affine_value_order_orientation_space_computed",
            "satisfied": True,
            "evidence": f"{len(scan)} affine orientation rows were audited over slopes {SLOPES} and intercepts {INTERCEPTS}",
        },
        {
            "criterion": "positive_orientation_conditionally_selects_unique_min_skeleton",
            "satisfied": all(row["all_orbits_unique"] for row in positive),
            "evidence": f"{sum(row['all_orbits_unique'] for row in positive)} of {len(positive)} positive rows select uniquely on all orbits",
        },
        {
            "criterion": "negative_orientation_rejected_by_unit_orbit_tie",
            "satisfied": all(not row["all_orbits_unique"] for row in negative),
            "evidence": f"{sum(not row['all_orbits_unique'] for row in negative)} of {len(negative)} negative rows have a tie",
        },
        {
            "criterion": "zero_orientation_rejected_by_collapse",
            "satisfied": all(not row["all_orbits_unique"] for row in zero),
            "evidence": f"{sum(not row['all_orbits_unique'] for row in zero)} of {len(zero)} zero-slope rows collapse distinctions",
        },
        {
            "criterion": "strict_positive_orientation_source_theorem_exported",
            "satisfied": False,
            "evidence": "the finite scan identifies a sufficient sign condition a>0, but no strict artifact sources that sign",
        },
        {
            "criterion": "delta_eta_and_beta_eta_coupling_exported",
            "satisfied": False,
            "evidence": "the value-order orientation gate is not coupled to delta/eta or beta/eta source theorems",
        },
    ]


def build_payload(_: dict[str, Any]) -> dict[str, Any]:
    scan = orientation_scan_rows()
    criteria = acceptance_rows(scan)
    candidates = source_candidate_rows()
    positive = [row for row in scan if row["orientation_class"] == "positive_min"]
    negative = [row for row in scan if row["orientation_class"] == "negative_max"]
    zero = [row for row in scan if row["orientation_class"] == "zero_collapse"]
    return {
        "status": "P2942_EXTREMAL_SELECTOR_VALUE_ORDER_ORIENTATION_SOURCE_GATE_NO_STRICT_POLARITY_SOURCE",
        "input_hashes": {"P2941": hashlib.sha256(P2941.read_bytes()).hexdigest() if P2941.exists() else None},
        "constructed_theoretical_objects": {
            "candidate_object": "Strict_Positive_Value_Order_Orientation_Source_For_P2938_Extremal_Selector",
            "affine_orientation_scan_rows": scan,
            "source_candidate_rows": candidates,
            "acceptance_rows": criteria,
        },
        "orientation_certificate": {
            "scan_row_count": len(scan),
            "positive_orientation_row_count": len(positive),
            "negative_orientation_row_count": len(negative),
            "zero_orientation_row_count": len(zero),
            "positive_rows_unique_on_all_orbits": sum(row["all_orbits_unique"] for row in positive),
            "negative_rows_unique_on_all_orbits": sum(row["all_orbits_unique"] for row in negative),
            "zero_rows_unique_on_all_orbits": sum(row["all_orbits_unique"] for row in zero),
            "strict_positive_orientation_source_theorem_exported": False,
            "strict_provenance_theorem_exported": False,
            "accepted_strict_source": False,
        },
        "decision": {
            "positive_witnesses": {
                "finite_gate_identifies_sufficient_orientation_sign": "a>0",
                "positive_orientation_conditionally_recovers_P2941_min_skeleton": True,
            },
            "negative_export_flags": {
                "strict_positive_orientation_source_theorem_exported": False,
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
            "reason": "P2942 turns the P2941 polarity gap into a finite value-order orientation gate.  Positive affine orientations conditionally recover a unique min-selector skeleton, while negative orientations inherit the max unit-orbit tie and zero orientations collapse distinctions.  The sign a>0 is therefore the exact missing strict source atom; it is not exported here and is not coupled to damping.",
            "next_honest_step": "A next admissible move must export one strict positive value-order orientation source theorem for a>0 and prove delta/eta plus beta/eta coupling for that theorem, or pivot to a genuinely new typed object outside the P2938 polarity/orientation lane.  Otherwise preserve the P2929-P2942 no-new-live-frontier boundary.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    cert = payload["orientation_certificate"]
    lines = [
        "# P2942/S1892 extremal-selector value-order orientation source gate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Orientation certificate",
        f"- scan rows: `{cert['scan_row_count']}`",
        f"- positive orientation rows: `{cert['positive_orientation_row_count']}`",
        f"- negative orientation rows: `{cert['negative_orientation_row_count']}`",
        f"- zero orientation rows: `{cert['zero_orientation_row_count']}`",
        f"- positive rows unique on all orbits: `{cert['positive_rows_unique_on_all_orbits']}`",
        f"- negative rows unique on all orbits: `{cert['negative_rows_unique_on_all_orbits']}`",
        f"- zero rows unique on all orbits: `{cert['zero_rows_unique_on_all_orbits']}`",
        f"- strict positive orientation source theorem exported: `{cert['strict_positive_orientation_source_theorem_exported']}`",
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
    payload = build_payload(read_json(P2941))
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2942/S1892 extremal-selector value-order orientation source gate", "## P2942/S1892 extremal-selector value-order orientation source gate\n\n`P2942/S1892` turns the P2941 min/max polarity gap into a finite value-order orientation gate by scanning affine orientations `q(v)=a*v+b`.  Positive slopes conditionally recover the unique P2941 min-selector skeleton on all audited orbits; negative slopes reverse to the max polarity and inherit the unit-orbit tie; zero slopes collapse the value distinctions.  Thus the exact missing strict source atom is a theorem sourcing the positive orientation sign `a>0`, followed by delta/eta and beta/eta coupling.  No strict provenance, strict `L_p`, nonproxy `L_total`, EOM, Hamiltonian, bridge closure, role transfer, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2942/S1892 value-order orientation `L_total` guard", "## P2942/S1892 value-order orientation `L_total` guard\n\n`P2942/S1892` identifies `a>0` in an affine value-order orientation `q(v)=a*v+b` as the exact missing polarity source for the conditional P2941 min-selector skeleton.  Since no strict positive-orientation source theorem or damping coupling is exported, the gate cannot enter nonproxy `L_total`, EOM, Hamiltonian, bridge closure, role transfer, or ToE.\n")
    append_once(AGENTS, "Current extremal-selector value-order orientation guardrail (P2942/S1892, 2026-06-19)", "## Current extremal-selector value-order orientation guardrail (P2942/S1892, 2026-06-19)\n\n- P2942 scans affine value-order orientations `q(v)=a*v+b` for the P2941 extremal-selector skeleton.\n- Positive slopes conditionally recover the unique min-selector skeleton; negative slopes reverse to the max selector and inherit the unit-orbit tie; zero slopes collapse distinctions.\n- The exact missing strict source atom is a nonconventional theorem sourcing the positive orientation sign `a>0`, plus delta/eta and beta/eta coupling.\n- Do not promote P2942 to selector closure, strict provenance, strict `L_p`, damping source, nonproxy `L_total`, bridge closure, role transfer, or ToE; next move must export that source theorem/coupling or pivot outside the P2938 polarity/orientation lane.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
