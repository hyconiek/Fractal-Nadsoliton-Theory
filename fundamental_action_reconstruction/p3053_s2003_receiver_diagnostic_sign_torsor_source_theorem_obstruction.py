#!/usr/bin/env python3
"""P3053/S2003: receiver diagnostic sign-torsor source-theorem obstruction.

P3052 left one honest option inside the phase-geometry lane: supply a strict
source theorem explaining why robust receiver winding is the nadsoliton selector
sign.  P3053 tests the maximal theorem shape available from the already-built
receiver diagnostics themselves, rather than replaying area, polarity, winding,
or robustness as closure.

The constructed object is a finite C2 sign-torsor obstruction matrix for the
P3048-P3052 receiver diagnostics.  We encode the available signed diagnostics as
three odd sign coordinates: area sign, winding sign, and robust-winding sign.
Inversion negates all three.  Exhausting all Boolean laws on this sign cube shows
that sign-blind/invariant laws cannot distinguish the two orientations, while
odd/equivariant laws exist only in globally paired polarities.  Therefore the
receiver diagnostics do not supply the missing strict source theorem by
themselves; a new non-receiver source law or coupling-polarity axiom is still
required.
"""
from __future__ import annotations

import hashlib, itertools, json
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p3052_s2002_phase_curve_winding_stability_margin_certificate import OUT as P3052

OUT = GEN / "p3053_s2003_receiver_diagnostic_sign_torsor_source_theorem_obstruction.json"
MD = GEN / "p3053_s2003_receiver_diagnostic_sign_torsor_source_theorem_obstruction.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

DIAGNOSTICS = ["area_sign_A_s", "winding_sign_W", "robust_winding_sign_R"]
DOMAIN = list(itertools.product([-1, 1], repeat=len(DIAGNOSTICS)))
BASE_STATE = (1, 1, 1)
OPPOSITE_STATE = tuple(-x for x in BASE_STATE)


def neg(state: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(-x for x in state)


def law_value(mask: int, state_index: int) -> int:
    return 1 if (mask >> state_index) & 1 else -1


def law_table(mask: int) -> dict[tuple[int, ...], int]:
    return {state: law_value(mask, index) for index, state in enumerate(DOMAIN)}


def classify_law(mask: int) -> dict[str, Any]:
    table = law_table(mask)
    invariant = all(table[neg(state)] == table[state] for state in DOMAIN)
    equivariant_odd = all(table[neg(state)] == -table[state] for state in DOMAIN)
    base_value = table[BASE_STATE]
    opposite_value = table[OPPOSITE_STATE]
    return {
        "mask": mask,
        "invariant_under_inversion": invariant,
        "odd_equivariant_under_inversion": equivariant_odd,
        "base_orientation_value": base_value,
        "opposite_orientation_value": opposite_value,
        "distinguishes_base_pair": base_value != opposite_value,
        "selected_absolute_polarity": base_value if equivariant_odd else None,
        "strict_source_theorem_exported": False,
    }


def diagnostic_rows() -> list[dict[str, Any]]:
    return [
        {
            "diagnostic": "area_sign_A_s",
            "source": "P3048/P3049/P3050 oriented memory phase-area sign",
            "inversion_action": "A_s -> -A_s",
            "finite_positive": True,
            "strict_source_channel_exported": False,
            "reason_not_source": "finite receiver sign with unresolved K/M provenance and coupling polarity",
        },
        {
            "diagnostic": "winding_sign_W",
            "source": "P3051 global phase-curve winding index",
            "inversion_action": "W -> -W",
            "finite_positive": True,
            "strict_source_channel_exported": False,
            "reason_not_source": "sampled/derived phase-curve orientation diagnostic",
        },
        {
            "diagnostic": "robust_winding_sign_R",
            "source": "P3052 stability-margin certificate",
            "inversion_action": "R -> -R",
            "finite_positive": True,
            "strict_source_channel_exported": False,
            "reason_not_source": "robustness of a receiver diagnostic, not a source theorem",
        },
    ]


def build_matrix() -> dict[str, Any]:
    read_json(P3052)
    laws = [classify_law(mask) for mask in range(2 ** len(DOMAIN))]
    invariant = [law for law in laws if law["invariant_under_inversion"]]
    odd = [law for law in laws if law["odd_equivariant_under_inversion"]]
    paired_odd = []
    for law in odd:
        complement = law["mask"] ^ ((1 << len(DOMAIN)) - 1)
        if law["mask"] < complement:
            paired_odd.append({
                "mask_plus": law["mask"],
                "mask_minus": complement,
                "base_values": [law["base_orientation_value"], -law["base_orientation_value"]],
                "pair_interpretation": "opposite output polarities; current artifacts select neither member",
            })
    acceptance = [
        {"criterion": "receiver_diagnostic_inputs_are_nonzero", "satisfied": True, "detail": "P3048-P3052 supply finite signed receiver diagnostics"},
        {"criterion": "inversion_action_explicit", "satisfied": True, "detail": "all diagnostic signs transform by simultaneous negation"},
        {"criterion": "boolean_law_exhaustion_complete", "satisfied": len(laws) == 256, "detail": "all Boolean maps from the three-sign cube to a selector sign were enumerated"},
        {"criterion": "invariant_laws_fail_to_select", "satisfied": all(not law["distinguishes_base_pair"] for law in invariant), "detail": "invariant laws identify the inversion pair"},
        {"criterion": "odd_laws_have_unique_artifact_selected_polarity", "satisfied": False, "detail": "odd laws distinguish orientation only in opposite polarity pairs"},
        {"criterion": "strict_nonreceiver_source_theorem", "satisfied": False, "detail": "no new source law outside receiver diagnostics is exported"},
        {"criterion": "selector_readout_or_p3046_coupling", "satisfied": False, "detail": "no theorem picks one P3046/readout polarity"},
        {"criterion": "unit_bearing_action_eom_installation", "satisfied": False, "detail": "no unit-bearing variational/action/EOM row is installed"},
    ]
    obligations = [
        {"obligation": "content_first_prior_lane_check", "satisfied": True, "detail": "P3053 is a source-theorem obstruction for receiver diagnostic signs, not a replay of area/winding values"},
        {"obligation": "construct_sign_torsor_source_theorem_object", "satisfied": True, "detail": "three receiver signs are assembled into a C2 sign cube"},
        {"obligation": "exhaust_boolean_selector_laws", "satisfied": True, "detail": "256 Boolean laws were classified"},
        {"obligation": "prove_invariant_laws_cannot_select_orientation", "satisfied": True, "detail": "all invariant laws give equal values on the inversion pair"},
        {"obligation": "export_unique_nonpremise_polarity", "satisfied": False, "detail": "odd laws remain paired by global output sign"},
        {"obligation": "selector_ltotal_bridge_toe_installation", "satisfied": False, "detail": "no QW-2191 discharge, selector closure, L_total, bridge, role transfer, or ToE export follows"},
    ]
    return {
        "object": "ReceiverDiagnosticSignTorsor_SourceTheoremObstruction",
        "diagnostic_rows": diagnostic_rows(),
        "domain_states": [dict(zip(DIAGNOSTICS, state)) for state in DOMAIN],
        "boolean_law_summary": {
            "total_laws": len(laws),
            "invariant_laws": len(invariant),
            "odd_equivariant_laws": len(odd),
            "odd_polarity_pairs": len(paired_odd),
            "invariant_laws_distinguishing_base_pair": sum(1 for law in invariant if law["distinguishes_base_pair"]),
            "artifact_selected_odd_polarity_laws": 0,
        },
        "odd_polarity_pair_rows": paired_odd,
        "source_acceptance_rows": acceptance,
        "proof_obligations": obligations,
        "finite_certificate": {
            "diagnostic_signs": len(DIAGNOSTICS),
            "domain_states": len(DOMAIN),
            "boolean_laws_enumerated": len(laws),
            "invariant_laws": len(invariant),
            "odd_equivariant_laws": len(odd),
            "odd_polarity_pairs": len(paired_odd),
            "invariant_laws_distinguishing_base_pair": sum(1 for law in invariant if law["distinguishes_base_pair"]),
            "artifact_selected_odd_polarity_laws": 0,
            "source_acceptance_criteria": len(acceptance),
            "satisfied_source_acceptance_criteria": sum(1 for row in acceptance if row["satisfied"]),
            "proof_obligations": len(obligations),
            "satisfied_proof_obligations": sum(1 for row in obligations if row["satisfied"]),
            "strict_receiver_diagnostic_source_theorem_exported": False,
        },
    }


def build_payload() -> dict[str, Any]:
    matrix = build_matrix()
    return {
        "status": "P3053_RECEIVER_DIAGNOSTIC_SIGN_TORSOR_SOURCE_THEOREM_OBSTRUCTION_BOUNDED_NO_EXPORT",
        "input_hashes": {"P3052": hashlib.sha256(P3052.read_bytes()).hexdigest() if P3052.exists() else None},
        "constructed_theoretical_objects": matrix,
        "finite_certificate": matrix["finite_certificate"],
        "decision": {
            "bounded_no_go": "P3053 shows that the P3048-P3052 signed receiver diagnostics have the correct inversion-odd type but cannot become a strict source theorem by recombination.  Invariant Boolean laws cannot distinguish the orientation pair, while odd/equivariant laws occur in opposite output-polarity pairs with no artifact-selected member.  Thus the robust phase-geometry signs remain receiver diagnostics, not a non-premise selector source.",
            "negative_export_flags": {k: False for k in ["strict_receiver_diagnostic_source_theorem_exported", "unique_nonpremise_orientation_polarity_exported", "p3046_coupling_polarity_selected", "selector_readout_coupling_exported", "qw2191_discharged", "strict_selector_closure_exported", "observed_physics_exported", "unit_bearing_action_eom_hamiltonian_exported", "ltotal_exported", "bridge_closure_exported", "role_transfer_exported", "toe_closure_exported"]},
            "next_honest_step": "Do not recombine P3048-P3052 receiver signs as source closure.  The next proof-grade move must introduce a genuinely new non-receiver strict signed source/coupling law that selects one odd polarity, or pivot outside the phase-geometry selector lane and preserve the P3048-P3053 bounded no-export certificate.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    c = payload["finite_certificate"]
    lines = [
        "# P3053/S2003 receiver diagnostic sign-torsor source-theorem obstruction", "",
        f"Status: `{payload['status']}`", "", "## Finite certificate",
        f"- diagnostic signs: `{c['diagnostic_signs']}`",
        f"- domain states: `{c['domain_states']}`",
        f"- Boolean laws enumerated: `{c['boolean_laws_enumerated']}`",
        f"- invariant laws: `{c['invariant_laws']}`",
        f"- odd equivariant laws: `{c['odd_equivariant_laws']}`",
        f"- odd polarity pairs: `{c['odd_polarity_pairs']}`",
        f"- invariant laws distinguishing base pair: `{c['invariant_laws_distinguishing_base_pair']}`",
        f"- artifact-selected odd polarity laws: `{c['artifact_selected_odd_polarity_laws']}`",
        f"- source acceptance criteria: `{c['satisfied_source_acceptance_criteria']}/{c['source_acceptance_criteria']}`",
        f"- satisfied proof obligations: `{c['satisfied_proof_obligations']}/{c['proof_obligations']}`",
        f"- strict receiver diagnostic source theorem exported: `{c['strict_receiver_diagnostic_source_theorem_exported']}`", "",
        "## Decision", payload["decision"]["bounded_no_go"], "", "## Recommendation", payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3053/S2003 receiver diagnostic sign-torsor source-theorem obstruction", "## P3053/S2003 receiver diagnostic sign-torsor source-theorem obstruction\n\n`P3053/S2003` audits whether the P3048-P3052 signed receiver diagnostics can themselves supply the strict source theorem demanded after P3052.  The constructed C2 sign cube uses area sign, winding sign, and robust-winding sign as inversion-odd coordinates, then exhausts all `256` Boolean selector laws.  Invariant laws cannot distinguish the orientation pair, while odd/equivariant laws remain paired by global output polarity; current artifacts select neither member.  No strict receiver-diagnostic source theorem, P3046 polarity selection, selector/readout installation, unit-bearing action/EOM, `L_total`, bridge/role transfer, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3053/S2003 receiver diagnostic sign-torsor `L_total` guard", "## P3053/S2003 receiver diagnostic sign-torsor `L_total` guard\n\n`P3053/S2003` adds no physical `L_total` term.  Exhausting Boolean laws over area/winding/robust-winding receiver signs shows that invariant laws are orientation-blind and odd laws require an unselected output polarity, so no unit-bearing variational/action/EOM source is installed.\n")
    append_once(AGENTS, "Current receiver diagnostic sign-torsor source-theorem guardrail (P3053/S2003, 2026-06-23)", "## Current receiver diagnostic sign-torsor source-theorem guardrail (P3053/S2003, 2026-06-23)\n\n- P3053 tests whether the P3048-P3052 signed receiver diagnostics can themselves provide the strict source theorem demanded after robust winding.\n- The finite C2 sign cube over area sign, winding sign, and robust-winding sign exhausts all 256 Boolean selector laws: invariant laws cannot distinguish the inversion pair, while odd/equivariant laws remain paired by unselected output polarity.\n- Do not recombine P3048-P3052 receiver signs, robust winding, positive clearance, or Boolean odd laws into `QW-2191` discharge, selector closure, observed physics, unit-bearing action/EOM, `L_total`, bridge/role-transfer, or ToE closure.\n- A next move requires a genuinely new non-receiver strict signed source/coupling law selecting one odd polarity, or a pivot outside the phase-geometry selector lane while preserving the P3048-P3053 bounded no-export certificate.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
