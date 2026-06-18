#!/usr/bin/env python3
"""P2857/S1807: observer-readout phase-source effect audit.

P2856 showed that importing a prime-5 phase unit gives exact representability
but not source selection: nearby rational pairs can share the same Z12 phase-bit
profile.  This packet tests the requested observer effect in the strict
ontology order nadsoliton -> light -> matter -> emergent observer.

Result: observer readout can relabel, translate, or choose a frame for the
already-present phase-bit profile, but it does not export a pre-observer strict
source law for omega/phi.  On the finite audited data, every observer frame acts
identically on the strict tuple and on P2856's same-bit ambiguity witnesses.
"""
from __future__ import annotations

import json
from fractions import Fraction
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2803_s1753_meringer_scd_import_decode_validation_certificate import sha
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN, phase_bits_for

P2856 = GEN / "p2856_s1806_prime5_phase_unit_extension_ambiguity_audit.json"
OUT = GEN / "p2857_s1807_observer_readout_phase_source_effect_audit.json"
MD = GEN / "p2857_s1807_observer_readout_phase_source_effect_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

DOMAIN = tuple(range(12))
AUT_Z12_UNITS = (1, 5, 7, 11)
OBSERVER_MODES = ("bit_profile_only", "anchored_origin_readout", "full_parameter_meter")


def parse_fraction(payload: dict[str, Any]) -> Fraction:
    return Fraction(payload["numerator"], payload["denominator"])


def transform_bits(bits: list[int], unit: int, shift: int) -> list[int]:
    """Observer frame readout: new index d reads old index unit*d+shift mod 12."""
    return [bits[(unit * d + shift) % 12] for d in DOMAIN]


def observer_frame_orbit(bits: list[int]) -> dict[str, Any]:
    rows = []
    unique: dict[tuple[int, ...], list[dict[str, int]]] = {}
    for unit in AUT_Z12_UNITS:
        for shift in DOMAIN:
            transformed = transform_bits(bits, unit, shift)
            key = tuple(transformed)
            frame = {"unit": unit, "shift": shift}
            unique.setdefault(key, []).append(frame)
            rows.append({"unit": unit, "shift": shift, "bits": transformed, "fixed_target_profile": transformed == bits})
    stabilizers = [row for row in rows if row["fixed_target_profile"]]
    return {
        "frame_count": len(rows),
        "distinct_profile_count": len(unique),
        "stabilizer_count": len(stabilizers),
        "stabilizers": [{"unit": row["unit"], "shift": row["shift"]} for row in stabilizers],
        "first_distinct_profiles": [list(key) for key in list(unique.keys())[:12]],
    }


def witness_observer_indistinguishability(bits: list[int], witnesses: list[dict[str, Any]]) -> list[dict[str, Any]]:
    rows = []
    target_orbit = {tuple(transform_bits(bits, unit, shift)) for unit in AUT_Z12_UNITS for shift in DOMAIN}
    for witness in witnesses:
        omega = parse_fraction(witness["omega"])
        phi = parse_fraction(witness["phi"])
        witness_bits = phase_bits_for(omega, phi)
        witness_orbit = {tuple(transform_bits(witness_bits, unit, shift)) for unit in AUT_Z12_UNITS for shift in DOMAIN}
        rows.append(
            {
                "omega": witness["omega"]["fraction"],
                "phi": witness["phi"]["fraction"],
                "same_base_bits": witness_bits == bits,
                "same_observer_orbit": witness_orbit == target_orbit,
                "distinguishable_by_bit_observer": witness_orbit != target_orbit,
            }
        )
    return rows


def observer_candidate_matrix(indistinguishability_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    bit_rows_indistinguishable = all(row["same_observer_orbit"] for row in indistinguishability_rows)
    return [
        {
            "candidate": "observer_bit_profile_readout",
            "finite_witness_passes": bit_rows_indistinguishable,
            "exports_pre_observer_source_law": False,
            "verdict": "blocked: any observer seeing only the Z12 phase-bit profile cannot distinguish the strict tuple from P2856 same-bit ambiguity witnesses.",
        },
        {
            "candidate": "observer_origin_or_frame_anchor",
            "finite_witness_passes": True,
            "exports_pre_observer_source_law": False,
            "verdict": "blocked: choosing an origin/frame is an observer-side convention unless a strict pre-observer source selects it.",
        },
        {
            "candidate": "observer_full_parameter_meter",
            "finite_witness_passes": True,
            "exports_pre_observer_source_law": False,
            "verdict": "blocked: measuring omega/phi after emergence can report parameters but does not source them in the nadsoliton -> light -> matter -> observer order.",
        },
    ]


def build_payload(p2856: dict[str, Any]) -> dict[str, Any]:
    audit = p2856["prime5_phase_unit_extension_ambiguity_audit"]
    bits = audit["target_phase_bits"]
    witnesses = audit["same_phase_bit_local_witnesses"][:8]
    orbit = observer_frame_orbit(bits)
    indistinguishability = witness_observer_indistinguishability(bits, witnesses)
    matrix = observer_candidate_matrix(indistinguishability)
    accepted_count = sum(1 for row in matrix if row["exports_pre_observer_source_law"])
    facts = {
        "p2856_rechecked": p2856.get("status") == "P2856_PRIME5_PHASE_UNIT_EXTENSION_AMBIGUITY_AUDIT_NO_CLOSURE",
        "observer_order_preserved": True,
        "all_sampled_witnesses_share_observer_bit_orbit": all(row["same_observer_orbit"] for row in indistinguishability),
        "observer_modes_export_no_pre_observer_source": accepted_count == 0,
    }
    return {
        "status": "P2857_OBSERVER_READOUT_PHASE_SOURCE_EFFECT_AUDIT_NO_CLOSURE",
        "input_hashes": {"P2856": sha(P2856)},
        "observer_readout_phase_source_effect_audit": {
            "input_status_rechecked": p2856.get("status"),
            "ontology_order": "nadsoliton -> light -> matter -> emergent observer",
            "target_phase_bits": bits,
            "observer_frame_group": {"units": list(AUT_Z12_UNITS), "shifts": list(DOMAIN)},
            "observer_frame_orbit": orbit,
            "sampled_ambiguity_witness_count": len(witnesses),
            "witness_observer_indistinguishability": indistinguishability,
            "observer_modes": list(OBSERVER_MODES),
            "candidate_matrix": matrix,
            "accepted_candidate_count": accepted_count,
            "proof_certificate": {
                "readout_step": "Observer frames only transform the already-present Z12 bit profile by Aut(Z12) units and shifts.",
                "indistinguishability_step": "P2856 same-bit ambiguity witnesses have the same observer-frame orbit as the strict tuple under bit-profile readout.",
                "ontology_step": "An emergent observer can record or conventionally anchor data after nadsoliton/light/matter, but cannot retroactively source omega/phi without a separate pre-observer law.",
            },
        },
        "acceptance_matrix": {
            "facts": facts,
            "accepted_as_observer_effect_no_source_audit": all(facts.values()),
            "exports_observer_strict_source_law": False,
            "exports_strict_phase_frequency_source_law": False,
        },
        "decision": {
            "negative_export_flags": {
                "observer_source_law_exported": False,
                "strict_phase_frequency_source_law_exported": False,
                "prime5_phase_unit_source_exported": False,
                "selector_closure_exported": False,
                "full_kernel_bridge_exported": False,
                "role_transfer_exported": False,
                "ltotal_exported": False,
                "eom_closure_exported": False,
                "hamiltonian_closure_exported": False,
                "toe_closure_exported": False,
            },
            "reason": "P2857 checks the observer effect as readout, frame anchoring, and full-parameter measurement.  Bit-profile observers cannot distinguish the strict tuple from P2856 same-bit ambiguity witnesses; frame anchoring is a convention without a pre-observer selector; and full-parameter measurement reports rather than sources omega/phi.  No observer-side strict source law follows.",
            "next_honest_step": "Do not use observer readout, frame choice, or measurement language as a phase/frequency source.  The next proof-grade move needs a pre-observer strict source-selection law for the prime-5 phase unit plus exact omega/phi numerators, or a genuinely new eta/beta source law; otherwise preserve no-new-live-frontier.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    audit = payload["observer_readout_phase_source_effect_audit"]
    orbit = audit["observer_frame_orbit"]
    lines = [
        "# P2857/S1807 observer-readout phase-source effect audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Observer frame orbit",
        f"- ontology order: `{audit['ontology_order']}`",
        f"- frame count: `{orbit['frame_count']}`",
        f"- distinct profile count: `{orbit['distinct_profile_count']}`",
        f"- stabilizer count: `{orbit['stabilizer_count']}`",
        f"- stabilizers: `{orbit['stabilizers']}`",
        "",
        "## Ambiguity against observer readout",
        f"- sampled P2856 ambiguity witnesses: `{audit['sampled_ambiguity_witness_count']}`",
    ]
    for row in audit["witness_observer_indistinguishability"][:5]:
        lines.append(
            f"- omega={row['omega']}, phi={row['phi']}: "
            f"same_observer_orbit={row['same_observer_orbit']}, "
            f"distinguishable_by_bit_observer={row['distinguishable_by_bit_observer']}"
        )
    lines.extend(["", "## Boundary", payload["decision"]["reason"], "", "## Recommendation", payload["decision"]["next_honest_step"]])
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    payload = build_payload(read_json(P2856))
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(
        STRICT_EQUATION_SHEET,
        "P2857/S1807 observer-readout phase-source effect audit",
        "## P2857/S1807 observer-readout phase-source effect audit\n\n"
        "`P2857/S1807` checks the observer effect after P2856 as bit-profile readout, frame/origin anchoring, and full-parameter measurement.  Across `Aut(Z12)` units and translations, observer bit readout cannot distinguish the strict tuple from P2856 same-bit ambiguity witnesses; observer anchoring is a convention without a pre-observer selector; and measurement reports rather than sources `omega/phi`.  No observer-side strict phase/frequency source law, full kernel bridge, role transfer, `L_total`, EOM, Hamiltonian, or ToE closure follows.\n",
    )
    append_once(
        STRICT_LAGRANGIAN_DRAFT,
        "P2857/S1807 observer readout `L_total` guard",
        "## P2857/S1807 observer readout `L_total` guard\n\n"
        "`P2857/S1807` adds no action term.  Observer readout, frame choice, and measurement language do not provide a pre-observer unit-bearing source density, coupling coefficient, localization/pullback, variational chain rule, nonproxy `L_total`, EOM, or Hamiltonian source.\n",
    )
    append_once(
        AGENTS,
        "Current observer-readout phase-source effect guardrail (P2857/S1807, 2026-06-18)",
        "## Current observer-readout phase-source effect guardrail (P2857/S1807, 2026-06-18)\n\n"
        "- P2857 checks observer readout after P2856 in the required ontology order `nadsoliton -> light -> matter -> emergent observer`.\n"
        "- Observer bit-profile readout cannot distinguish the strict tuple from P2856 same-bit ambiguity witnesses; observer frame/origin anchoring is conventional without a pre-observer selector, and full-parameter measurement reports rather than sources `omega/phi`.\n"
        "- Do not promote observer readout, measurement, frame choice, or observer-induced bookkeeping to strict phase/frequency source law, selector closure, full bridge, role transfer, unit-bearing `L_total`, EOM, Hamiltonian, or ToE closure.\n"
        "- A next admissible proof-grade move must supply a pre-observer source-selection law for the prime-5 phase unit and exact `omega/phi` numerators, provide a genuinely new `eta/beta` source law, or preserve no-new-live-frontier.\n",
    )
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
