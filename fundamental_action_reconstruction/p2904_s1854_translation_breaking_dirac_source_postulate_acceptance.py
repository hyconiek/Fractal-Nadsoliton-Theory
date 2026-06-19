#!/usr/bin/env python3
"""P2904/S1854: translation-breaking Dirac source postulate acceptance.

P2903 proved that translation-neutral internal data cannot select the pointed
signed axiom A(b,sigma).  P2904 therefore constructs the minimal object that
would pass that specific obstruction: a signed Dirac source Xi on Z12 with one
nonzero value.  This object is genuinely translation-breaking and couples
computably to A(0,+), but it is a *new source postulate/candidate* rather than a
strict nadsoliton-derived theorem.  The audit records the positive construction
and the remaining provenance/physical-unit blockers without promoting closure.
"""
from __future__ import annotations

import json
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2803_s1753_meringer_scd_import_decode_validation_certificate import sha
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN

P2903 = GEN / "p2903_s1853_internal_pointed_axiom_fixed_point_obstruction.json"
OUT = GEN / "p2904_s1854_translation_breaking_dirac_source_postulate_acceptance.json"
MD = GEN / "p2904_s1854_translation_breaking_dirac_source_postulate_acceptance.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"
N = 12
CARRIER_STEP = 5
UNIT_SYMBOL = "U_9_5"


def dirac_source(basepoint: int = 0, sign: int = 1) -> list[int]:
    return [sign if n == basepoint % N else 0 for n in range(N)]


def translate_signal(signal: list[int], shift: int) -> list[int]:
    return [signal[(n - shift) % N] for n in range(N)]


def stabilizer(signal: list[int]) -> list[int]:
    return [shift for shift in range(N) if translate_signal(signal, shift) == signal]


def orbit(signal: list[int]) -> list[list[int]]:
    seen = []
    for shift in range(N):
        translated = translate_signal(signal, shift)
        if translated not in seen:
            seen.append(translated)
    return seen


def coupling_to_pointed_axiom(signal: list[int]) -> dict[str, Any]:
    nonzero = [(idx, value) for idx, value in enumerate(signal) if value != 0]
    if len(nonzero) != 1:
        return {"coupling_defined": False, "reason": "source must have exactly one nonzero signed support point"}
    basepoint, value = nonzero[0]
    sigma = 1 if value > 0 else -1
    edge = [basepoint, (basepoint + sigma * CARRIER_STEP) % N]
    return {
        "coupling_defined": True,
        "selected_axiom": {"basepoint": basepoint, "polarity": sigma, "name": f"A({basepoint},{'+' if sigma > 0 else '-'})"},
        "defect_edge": edge,
        "density_template": f"{sigma}*{UNIT_SYMBOL}*delta_edge({edge[0]},{edge[1]})*q_9_5({edge[0]},{edge[1]})",
        "coupling_theorem": "A singleton signed Dirac source Xi with support b and sign sigma maps uniquely to A(b,sigma) and then to D=(b,b+sigma*5).",
    }


def build_payload(p2903: dict[str, Any]) -> dict[str, Any]:
    source = dirac_source(0, 1)
    stab = stabilizer(source)
    orb = orbit(source)
    coupling = coupling_to_pointed_axiom(source)
    return {
        "status": "P2904_TRANSLATION_BREAKING_DIRAC_SOURCE_POSTULATE_ACCEPTED_AS_CANDIDATE_NO_STRICT_PROVENANCE",
        "input_hashes": {"P2903": sha(P2903)},
        "constructed_theoretical_objects": {
            "source_object": "Xi_{0,+}: signed Dirac source on Z12",
            "source_values": source,
            "nonzero_value": 1,
            "support": [0],
            "stabilizer_under_translation": stab,
            "translation_orbit_size": len(orb),
            "coupling_to_pointed_axiom": coupling,
        },
        "acceptance_matrix": {
            "p2903_rechecked": p2903.get("status") == "P2903_INTERNAL_POINTED_AXIOM_FIXED_POINT_OBSTRUCTION_NO_STRICT_SOURCE",
            "translation_breaking_source_constructed": True,
            "nonzero_computed_value": source[0],
            "unique_signed_support_count": 1,
            "translation_stabilizer_size": len(stab),
            "translation_orbit_size": len(orb),
            "coupling_to_one_pointed_axiom_constructed": coupling["coupling_defined"],
            "selected_basepoint": coupling["selected_axiom"]["basepoint"],
            "selected_polarity": coupling["selected_axiom"]["polarity"],
            "defect_edge": coupling["defect_edge"],
            "passes_p2903_fixed_point_obstruction": True,
            "strict_nadsoliton_provenance_exported": False,
            "unit_bearing_nonproxy_ltotal_coupling_exported": False,
        },
        "decision": {
            "positive_witnesses": {
                "genuinely_translation_breaking_candidate_constructed": True,
                "nonzero_signed_value_computed": True,
                "unique_A_b_sigma_coupling_constructed": True,
                "defect_edge_and_density_template_computed": True,
            },
            "negative_export_flags": {
                "strict_nadsoliton_source_provenance_exported": False,
                "strict_defect_placement_source_law_exported": False,
                "unit_bearing_strict_density_exported": False,
                "nonproxy_ltotal_exported": False,
                "eom_closure_exported": False,
                "hamiltonian_closure_exported": False,
                "strict_damping_compression_bridge_exported": False,
                "role_transfer_exported": False,
                "toe_closure_exported": False,
            },
            "reason": "P2904 constructs the minimal translation-breaking object demanded by P2903: Xi_{0,+} has one nonzero signed value, trivial translation stabilizer, orbit size 12, and a unique coupling to A(0,+), D=(0,5), and the symbolic rho_9/5 template.  This passes the fixed-point obstruction as a candidate/source postulate.  It still does not prove strict nadsoliton provenance for why Xi_{0,+} rather than a translate or sign flip is exported, and it does not lift U_9_5 to a unit-bearing nonproxy L_total coupling.  Therefore closure remains blocked.",
            "next_honest_step": "The next proof-grade move should audit provenance of the new Xi_{0,+} source candidate: either derive its support and sign from a strict nadsoliton asymmetry/chiral/defect-generation theorem, or prove that no current artifact sources Xi rather than its 23 translated/sign-flipped alternatives.  Do not spend the next step on more translation-neutral selectors or more pointed templates.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    acc = payload["acceptance_matrix"]
    lines = [
        "# P2904/S1854 translation-breaking Dirac source postulate acceptance",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Constructed source candidate",
        "- `Xi_{0,+}`: signed Dirac source on `Z12`",
        f"- source values: `{payload['constructed_theoretical_objects']['source_values']}`",
        f"- nonzero computed value: `{acc['nonzero_computed_value']}`",
        f"- translation stabilizer size: `{acc['translation_stabilizer_size']}`",
        f"- translation orbit size: `{acc['translation_orbit_size']}`",
        "",
        "## Coupling",
        f"- selected basepoint: `{acc['selected_basepoint']}`",
        f"- selected polarity: `{acc['selected_polarity']}`",
        f"- defect edge: `{acc['defect_edge']}`",
        f"- coupling constructed: `{acc['coupling_to_one_pointed_axiom_constructed']}`",
        "",
        "## Boundary",
        payload["decision"]["reason"],
        "",
        "## Recommendation",
        payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload(read_json(P2903))
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2904/S1854 translation-breaking Dirac source postulate acceptance", "## P2904/S1854 translation-breaking Dirac source postulate acceptance\n\n`P2904/S1854` constructs the minimal translation-breaking object demanded by `P2903`: a signed Dirac source candidate `Xi_{0,+}` on `Z12` with one nonzero value, trivial translation stabilizer, orbit size `12`, and a unique coupling to `A(0,+)`, defect edge `D=(0,5)`, and the symbolic `rho_9/5` template.  This passes the fixed-point obstruction only as a new source candidate/postulate.  It does not prove strict nadsoliton provenance for choosing `Xi_{0,+}` over its translated/sign-flipped alternatives and does not lift `U_9_5` to a unit-bearing nonproxy `L_total` coupling; no EOM, Hamiltonian, bridge, role transfer, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2904/S1854 Dirac source candidate `L_total` guard", "## P2904/S1854 Dirac source candidate `L_total` guard\n\n`P2904/S1854` supplies a translation-breaking signed Dirac source candidate and computes its symbolic `rho_9/5` support, but it is still a provenance-unproved source postulate and keeps `U_9_5` symbolic.  Therefore it does not export unit-bearing nonproxy `L_total`, EOM, Hamiltonian closure, role transfer, or ToE closure.\n")
    append_once(AGENTS, "Current translation-breaking Dirac source candidate guardrail (P2904/S1854, 2026-06-19)", "## Current translation-breaking Dirac source candidate guardrail (P2904/S1854, 2026-06-19)\n\n- P2904 constructs a minimal translation-breaking signed Dirac source candidate `Xi_{0,+}` with one nonzero value, trivial translation stabilizer, orbit size `12`, and unique coupling to `A(0,+)` and `D=(0,5)`.\n- Treat this as a source candidate/postulate only: current artifacts do not prove strict nadsoliton provenance for choosing `Xi_{0,+}` over translated/sign-flipped alternatives and do not lift `U_9_5` to a unit-bearing nonproxy `L_total` coupling.\n- Do not promote `Xi_{0,+}`, its chosen coordinate `0`, its positive sign, symbolic `rho_9/5`, or `U_9_5` to strict sourcehood, nonproxy `L_total`, EOM, Hamiltonian, bridge closure, role transfer, or ToE closure without a provenance theorem and unit-bearing coupling theorem.\n- A next admissible proof-grade move must audit/prove strict provenance for `Xi_{0,+}` or prove no current artifact sources it over its `23` translated/sign-flipped alternatives; otherwise pivot outside this source-candidate lane or preserve no-new-live-frontier.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
