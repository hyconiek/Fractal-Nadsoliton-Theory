#!/usr/bin/env python3
"""P2901/S1851: explicit defect-placement source-law candidate gate.

P2900 found no current artifact exporting the missing object.  P2901 therefore
constructs the nearest explicit theoretical candidate instead of searching again:
a 9/5 carrier-coupled defect-placement law parameterized by a signed phase
source (b, sigma) on the free Z12 torsor.

The construction is intentionally acceptance-gated.  It proves that the formula
has the right *shape* if a nonimported strict signed phase source is supplied,
but it also computes that the current construction itself still has 24 imported
parameter choices (12 basepoints times 2 polarities).  Hence it is a conditional
candidate/schema, not a strict source export or ToE closure.
"""
from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2803_s1753_meringer_scd_import_decode_validation_certificate import sha
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN

P2900 = GEN / "p2900_s1850_defect_placement_source_law_inventory_no_go.json"
OUT = GEN / "p2901_s1851_explicit_defect_placement_source_law_candidate_gate.json"
MD = GEN / "p2901_s1851_explicit_defect_placement_source_law_candidate_gate.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"
N = 12
CARRIER_STEP = 5  # the named 9/5 carrier offset, represented as a Z12 step.


def z12() -> range:
    return range(N)


def source_signal(basepoint: int, polarity: int) -> list[float]:
    """A concrete signed phase profile whose zero/edge anchors a defect.

    The profile is not claimed to be strict-sourced.  It is a candidate formula:
    I_{b,sigma}(n)=sigma*sin(2*pi*(n-b)/12).  If b and sigma were exported by a
    strict law, the positive first 9/5 edge would be (b, b+5) and the reversed
    polarity would place the opposite signed edge.
    """
    return [polarity * math.sin(2 * math.pi * ((n - basepoint) % N) / N) for n in z12()]


def defect_edge(basepoint: int, polarity: int) -> tuple[int, int]:
    if polarity not in (-1, 1):
        raise ValueError("polarity must be +/-1")
    return (basepoint % N, (basepoint + polarity * CARRIER_STEP) % N)


def density_support(basepoint: int, polarity: int) -> dict[str, Any]:
    tail, head = defect_edge(basepoint, polarity)
    return {
        "formal_density_name": "rho_9_5_defect[b,sigma]",
        "support_edge": [tail, head],
        "carrier_offset_mod_12": (head - tail) % N,
        "opposite_orientation_offset_mod_12": (tail - head) % N,
        "formal_expression": "sigma * delta_{(b, b + sigma*5)} * q_9_5",
        "unit_status": "formal_symbolic_only_no_unit_bearing_Ltotal_export",
    }


def candidate_rows() -> list[dict[str, Any]]:
    rows = []
    for b in z12():
        for sigma in (-1, 1):
            signal = source_signal(b, sigma)
            rows.append({
                "basepoint": b,
                "polarity": sigma,
                "defect_edge": list(defect_edge(b, sigma)),
                "signal_zero_at_basepoint": abs(signal[b]) < 1e-12,
                "signal_at_forward_step_5": signal[(b + CARRIER_STEP) % N],
                "density_support": density_support(b, sigma),
            })
    return rows


def translate_row(row: dict[str, Any], t: int) -> tuple[int, int, tuple[int, int]]:
    b = (row["basepoint"] + t) % N
    sigma = row["polarity"]
    return b, sigma, defect_edge(b, sigma)


def orbit_representatives(rows: list[dict[str, Any]]) -> dict[str, Any]:
    unseen = {(row["basepoint"], row["polarity"]) for row in rows}
    orbits = []
    while unseen:
        b, sigma = min(unseen)
        orbit = {((b + t) % N, sigma) for t in z12()}
        unseen -= orbit
        orbits.append(sorted([{"basepoint": x, "polarity": y, "edge": list(defect_edge(x, y))} for x, y in orbit], key=lambda r: r["basepoint"]))
    return {"orbit_count": len(orbits), "orbit_sizes": [len(o) for o in orbits], "orbits": orbits}


def build_payload(p2900: dict[str, Any]) -> dict[str, Any]:
    rows = candidate_rows()
    orbits = orbit_representatives(rows)
    unique_unimported_choice_count = 0
    imported_parameter_count = len(rows)
    has_correct_formal_shape = all(row["signal_zero_at_basepoint"] and row["density_support"]["carrier_offset_mod_12"] in (5, 7) for row in rows)
    return {
        "status": "P2901_EXPLICIT_DEFECT_PLACEMENT_SOURCE_LAW_CANDIDATE_CONDITIONAL_NO_CLOSURE",
        "input_hashes": {"P2900": sha(P2900)},
        "constructed_theoretical_object": {
            "object_name": "conditional signed-phase defect-placement source-law schema",
            "formula": "I_{b,sigma}(n)=sigma*sin(2*pi*(n-b)/12); D_{b,sigma}=(b,b+sigma*5); rho_{9/5}=sigma*delta_{D_{b,sigma}}*q_{9/5}",
            "intended_role": "If a strict nonimported law supplied (b,sigma), this formula would compute a basepoint, polarity, one directed defect placement, and a formal 9/5 carrier-density support.",
            "constructed_rows": rows,
            "translation_orbits": orbits,
        },
        "acceptance_matrix": {
            "p2900_rechecked": p2900.get("status") == "P2900_DEFECT_PLACEMENT_SOURCE_LAW_INVENTORY_NO_GO_NO_CLOSURE",
            "formal_candidate_constructed": True,
            "has_correct_formal_defect_and_9_over_5_shape": has_correct_formal_shape,
            "candidate_parameter_count": imported_parameter_count,
            "translation_orbit_count": orbits["orbit_count"],
            "translation_orbit_sizes": orbits["orbit_sizes"],
            "nonimported_basepoint_supplied_by_formula": False,
            "nonimported_polarity_supplied_by_formula": False,
            "strict_source_theorem_exported": False,
            "unique_unimported_choice_count": unique_unimported_choice_count,
            "accepted_as_missing_object": False,
        },
        "decision": {
            "positive_witnesses": {
                "explicit_formula_schema_exists": True,
                "computes_defect_edge_given_basepoint_and_polarity": True,
                "computes_formal_9_over_5_density_support_given_parameters": True,
            },
            "negative_export_flags": {
                "strict_defect_placement_source_law_exported": False,
                "nonimported_basepoint_or_polarity_law_exported": False,
                "coupling_to_9_over_5_variational_density_exported": False,
                "nonimported_9_over_5_variational_chain_rule_exported": False,
                "localized_unit_bearing_action_density_exported": False,
                "strict_damping_compression_bridge_exported": False,
                "full_kernel_bridge_exported": False,
                "role_transfer_exported": False,
                "ltotal_exported": False,
                "eom_closure_exported": False,
                "hamiltonian_closure_exported": False,
                "toe_closure_exported": False,
            },
            "reason": "P2901 constructs the missing-object-shaped formula explicitly, but the formula is a 24-member family indexed by imported (basepoint, polarity). Translation acts freely on each polarity family, leaving two length-12 orbits and no internally selected member. Therefore P2901 is a conditional schema/readiness witness, not a strict source theorem or variational L_total construction.",
            "next_honest_step": "The next proof-grade move should attack exactly the remaining premise exposed by P2901: supply a nonimported strict law that selects one (b,sigma) pair, then prove the unit-bearing variational coupling of rho_9/5 into L_total. If no such selector/source law is supplied, pivot outside torsor/basepoint/scalar-score/relation/defect/support/orbit/Fourier/inventory families or preserve no-new-live-frontier.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    acc = payload["acceptance_matrix"]
    lines = [
        "# P2901/S1851 explicit defect-placement source-law candidate gate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Constructed object",
        f"`{payload['constructed_theoretical_object']['formula']}`",
        "",
        "## Finite acceptance gate",
        f"- candidate parameter count: `{acc['candidate_parameter_count']}`",
        f"- translation orbit count: `{acc['translation_orbit_count']}`",
        f"- translation orbit sizes: `{acc['translation_orbit_sizes']}`",
        f"- unique unimported choices: `{acc['unique_unimported_choice_count']}`",
        f"- accepted as missing object: `{acc['accepted_as_missing_object']}`",
        "",
        "## Boundary",
        payload["decision"]["reason"],
        "",
        "## Recommendation",
        payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload(read_json(P2900))
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2901/S1851 explicit defect-placement source-law candidate gate", "## P2901/S1851 explicit defect-placement source-law candidate gate\n\n`P2901/S1851` constructs an explicit conditional missing-object-shaped schema `I_{b,sigma}(n)=sigma*sin(2*pi*(n-b)/12)`, `D_{b,sigma}=(b,b+sigma*5)`, and formal `rho_{9/5}=sigma*delta_{D_{b,sigma}}*q_{9/5}`.  The finite gate finds `24` parameterized candidates split into two translation orbits of size `12`, with `0` internally selected unimported `(b,sigma)` choices.  Thus the construction is a readiness/schema witness only: it does not export a strict defect-placement source law, unit-bearing localized action density, variational chain rule, `L_total`, EOM, Hamiltonian, role transfer, bridge closure, or ToE closure.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2901/S1851 conditional defect-placement schema `L_total` guard", "## P2901/S1851 conditional defect-placement schema `L_total` guard\n\n`P2901/S1851` gives a formal `9/5` support-density schema only after imported `(basepoint, polarity)` data are supplied.  Because no nonimported source theorem selects that pair and no unit-bearing variational coupling theorem is proved, the schema is not a nonproxy `L_total`, EOM, or Hamiltonian construction.\n")
    append_once(AGENTS, "Current explicit defect-placement source-law candidate guardrail (P2901/S1851, 2026-06-19)", "## Current explicit defect-placement source-law candidate guardrail (P2901/S1851, 2026-06-19)\n\n- P2901 constructs a conditional signed-phase defect-placement schema `I_{b,sigma}(n)=sigma*sin(2*pi*(n-b)/12)`, `D_{b,sigma}=(b,b+sigma*5)`, and formal `rho_{9/5}` support.\n- The finite gate has `24` parameterized candidates, two translation orbits of size `12`, and `0` unimported internally selected `(basepoint, polarity)` choices.\n- Do not promote the conditional formula schema, imported `(b,sigma)` parameters, formal support density, chosen carrier edge, or translation-orbit representatives to strict phase/origin sourcehood, strict damping/compression bridge, role transfer, `L_total`, EOM, Hamiltonian, or ToE closure.\n- A next admissible proof-grade move must supply a nonimported strict law selecting one `(b,sigma)` pair and then a unit-bearing `rho_9/5 -> L_total` variational coupling theorem, pivot to a genuinely different typed object outside torsor/basepoint/scalar-score/relation/defect/support/orbit/Fourier/inventory constructions, or preserve no-new-live-frontier.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
