#!/usr/bin/env python3
"""Scratch probe: phase-resonance / Pythagorean-limma audit for strict alpha -> eta.

This packet checks a proposed physical reading of the exact radical payload

    q^5 = 256/243 = 2^8/3^5.

What survives arithmetic audit:

1. 256/243 is exactly the Pythagorean limma and exactly the fifth power of the
   radical correction q=(256/243)^(1/5).
2. A phase-labelled Z_12 / Fourier reading can conditionally motivate labelled
   channels, but only after a generator/orientation or source phase is fixed.
3. The common claim that 256/243 is the 12-fifth closure residual is false: the
   12-fifth / 7-octave closure residual is the Pythagorean comma
   3^12/2^19, not the limma.

So resonance is a useful candidate source of labelled channels, but the honest
next obligation is a strict phase/generator/source theorem.  The limma identity
alone does not discharge QW-2191 or prove eta=9/5.
"""
from __future__ import annotations

import json
import math
from fractions import Fraction
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
BINARY_RESCALE = HERE / "bridge_strict_alpha_binary_rescale_quotient_report.json"
OUT_JSON = HERE / "bridge_strict_alpha_phase_resonance_limmatic_audit_report.json"
OUT_MD = HERE / "bridge_strict_alpha_phase_resonance_limmatic_audit_report.md"

LIMMA = Fraction(256, 243)
PYTHAGOREAN_COMMA = Fraction(3**12, 2**19)
PYTHAGOREAN_APOTOME = Fraction(3**7, 2**11)
PYTHAGOREAN_WHOLE_TONE = Fraction(9, 8)
PERFECT_FOURTH = Fraction(4, 3)
PERFECT_FIFTH = Fraction(3, 2)
Z12_AUT_UNITS = [1, 5, 7, 11]
Z12_FIFTH_STEP = 7
Z12_NODE_COUNT = 12
STRICT_TARGET_ETA = 9.0 / 5.0
NAD12_SUPPORT_SIZE = 12
ALPHA_SCALE = 4.0


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def fraction_label(value: Fraction) -> str:
    return f"{value.numerator}/{value.denominator}"


def eta_from_correction(correction: float) -> float:
    return math.log(NAD12_SUPPORT_SIZE * correction) / math.log(ALPHA_SCALE)


def z12_orbit(value: int) -> list[int]:
    return sorted({(unit * value) % Z12_NODE_COUNT for unit in Z12_AUT_UNITS})


def fifth_cycle_mod_12() -> list[int]:
    current = 0
    out = []
    for _ in range(Z12_NODE_COUNT):
        out.append(current)
        current = (current + Z12_FIFTH_STEP) % Z12_NODE_COUNT
    return out


def interval_identity_rows() -> list[dict[str, Any]]:
    rows = []
    identities = [
        ("limma", LIMMA, "2^8/3^5"),
        ("pythagorean_comma_12_fifths_vs_7_octaves", PYTHAGOREAN_COMMA, "3^12/2^19"),
        ("apotome", PYTHAGOREAN_APOTOME, "3^7/2^11"),
        ("whole_tone", PYTHAGOREAN_WHOLE_TONE, "9/8"),
        ("perfect_fourth", PERFECT_FOURTH, "4/3"),
    ]
    for name, value, closed_form in identities:
        rows.append(
            {
                "name": name,
                "fraction": fraction_label(value),
                "closed_form": closed_form,
                "decimal": float(value),
                "log2_size": math.log2(float(value)),
            }
        )
    return rows


def main() -> None:
    binary_rescale_report = load_json(BINARY_RESCALE)
    q_radical = float(LIMMA) ** (1.0 / 5.0)
    eta = eta_from_correction(q_radical)
    fifth_cycle = fifth_cycle_mod_12()
    limma_from_fourth_minus_two_tones = PERFECT_FOURTH / (PYTHAGOREAN_WHOLE_TONE**2)
    limma_times_apotome = LIMMA * PYTHAGOREAN_APOTOME
    twelve_fifth_closure = PERFECT_FIFTH**12 / (2**7)

    report = {
        "status": "OPEN_STRICT_ALPHA_PHASE_RESONANCE_LIMMATIC_AUDIT_NO_RESONANCE_SELECTOR_THEOREM",
        "result_kind": "SCRATCH_STRICT_ALPHA_PHASE_RESONANCE_LIMMATIC_AUDIT_PROBE__NOT_A_THEOREM",
        "source_reports": {
            "binary_rescale_quotient": str(BINARY_RESCALE.relative_to(ROOT)),
        },
        "limma_identity_check": {
            "q_radical_power_5": fraction_label(LIMMA),
            "q_radical": q_radical,
            "q_radical_closed_form": "(256/243)^(1/5)",
            "eta_from_q_radical": eta,
            "eta_residual_vs_9_5": eta - STRICT_TARGET_ETA,
            "limma_is_2_8_over_3_5": LIMMA == Fraction(2**8, 3**5),
            "limma_from_fourth_div_two_whole_tones": fraction_label(limma_from_fourth_minus_two_tones),
            "limma_from_fourth_identity_ok": limma_from_fourth_minus_two_tones == LIMMA,
            "limma_times_apotome_equals_whole_tone": fraction_label(limma_times_apotome),
            "limma_apotome_identity_ok": limma_times_apotome == PYTHAGOREAN_WHOLE_TONE,
        },
        "twelve_node_resonance_audit": {
            "z12_fifth_step": Z12_FIFTH_STEP,
            "fifth_step_is_generator_mod_12": math.gcd(Z12_FIFTH_STEP, Z12_NODE_COUNT) == 1,
            "fifth_cycle_mod_12": fifth_cycle,
            "fifth_cycle_visits_all_nodes": sorted(fifth_cycle) == list(range(Z12_NODE_COUNT)),
            "aut_z12_units": Z12_AUT_UNITS,
            "orbit_of_fifth_step_under_aut_z12": z12_orbit(Z12_FIFTH_STEP),
            "fifth_step_equivalent_to_generator_only_after_generator_gauge_choice": True,
        },
        "closure_residual_discriminator": {
            "claim_under_audit": "256/243 is the 12-fifth closure residual on a 12-node/octave ring",
            "actual_12_fifth_closure_residual": fraction_label(twelve_fifth_closure),
            "actual_12_fifth_closure_name": "Pythagorean comma",
            "limma_fraction": fraction_label(LIMMA),
            "limma_equals_12_fifth_closure_residual": LIMMA == twelve_fifth_closure,
            "ratio_limmatic_to_comma": fraction_label(LIMMA / twelve_fifth_closure),
            "verdict": "REJECT_12_FIFTH_CLOSURE_EQUALS_LIMMA; ACCEPT_LIMMA_AS_PYTHAGOREAN_INTERVAL_IDENTITY",
        },
        "interval_rows": interval_identity_rows(),
        "phase_labelled_selector_interpretation": {
            "supported_conditionally": True,
            "content": "A fixed phase source/generator on Z_12 would make channels labelled by phase distance, which aligns with the labelled-selector probes.",
            "strict_core_blocker": "N462 boundary: typed Z_12 plus Aut(Z_12) alone does not canonically choose a generator/orientation.",
            "not_supported_without_extra_premise": "Phase labels do not by themselves prove strict-core selector closure unless a strict phase source or explicit generator/orientation premise is supplied.",
            "binary_rescale_context_result_kind": binary_rescale_report["result_kind"],
        },
        "candidate_interpretation": {
            "supported_by_this_probe": bool(
                LIMMA == Fraction(2**8, 3**5)
                and limma_from_fourth_minus_two_tones == LIMMA
                and limma_times_apotome == PYTHAGOREAN_WHOLE_TONE
                and twelve_fifth_closure == PYTHAGOREAN_COMMA
                and LIMMA != twelve_fifth_closure
                and sorted(fifth_cycle) == list(range(Z12_NODE_COUNT))
            ),
            "content": "The limma/resonance reading is arithmetically real, but the specific 12-fifth closure claim is false and phase-labelled selection remains premise-based.",
            "why_this_is_more_proof_like": "It separates exact interval identities, Z_12 phase-cycle facts, and the false 12-fifth closure identification by machine-checkable fractions.",
            "why_this_is_not_enough": "It does not derive a strict phase resonance selector, a canonical generator/orientation, or eta=9/5 from strict nadsoliton geometry.",
        },
        "hard_limits": [
            "No identity K_legacy_ont == K_strict_gate is exported.",
            "No legacy physical-role transfer is licensed.",
            "The 12-fifth closure residual is 3^12/2^19, not 256/243.",
            "No theorem derives a highest-resonance phase selector from strict nadsoliton geometry.",
            "No Aut(Z_12)-invariant generator/orientation fixing is claimed; phase labelling remains premise-based under the N462 boundary.",
            "No selector theorem chooses (2,2,2,1,1) from strict-core resonance alone.",
            "No theorem derives eta=9/5 without adopting the branch model, phase/generator premise, and selector convention as extra premises.",
            "No QW-2191 selector discharge and no ToE closure are claimed.",
        ],
        "next_honest_step": "Build a strict phase-resonance selector only if a generator/source-phase premise is explicitly supplied or derived; otherwise keep the limma identity as a suggestive but non-closing resonance shadow.",
    }

    OUT_JSON.write_text(json.dumps(report, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    OUT_MD.write_text(
        "# Scratch strict-alpha phase resonance limmatic audit probe\n\n"
        "Status: limma/resonance audit for eta=9/5; no resonance selector theorem.\n\n"
        f"- Accepted identity: `q^5={fraction_label(LIMMA)}=2^8/3^5`, the Pythagorean limma; eta residual `{eta-STRICT_TARGET_ETA:.3e}`.\n"
        f"- Rejected claim: 12-fifth closure residual is `{fraction_label(twelve_fifth_closure)}` (Pythagorean comma), not `{fraction_label(LIMMA)}`.\n"
        f"- Z12 check: fifth step `{Z12_FIFTH_STEP}` is a generator and visits all 12 nodes, but phase labels need a generator/source premise.\n"
        "- Honest read: resonance can motivate labelled channels only conditionally; N462 blocks free strict-core generator/orientation canonicity.\n"
        "- No false pass: no resonance selector theorem, no QW-2191 discharge, no ToE closure.\n",
        encoding="utf-8",
    )
    print(OUT_JSON)
    print(OUT_MD)


if __name__ == "__main__":
    main()
