#!/usr/bin/env python3
"""Scratch probe: C12 source-translation no-go after chiral orientation selection.

The chiral bispectrum audit showed a conditional gain: if a handedness convention
is allowed, C12-invariant / reflection-chiral bispectrum signs can distinguish
the two d5 orientations.  This probe asks the remaining source question.

Result: once an orientation is fixed, the 12 possible sources are one free C12
translation orbit.  Any C12-translation-invariant score is constant on that orbit
and cannot choose the source.  The computation verifies this with complete
Fourier power data and complete bispectrum signatures B(k,l)=F_k F_l conj(F_{k+l})
for all k,l in Z12: both are constant over all 12 source shifts for a fixed
orientation.  Raw Fourier phases vary over sources, but using them requires an
origin/phase reference, exactly the missing non-strict premise.
"""
from __future__ import annotations

import cmath
import json
import math
from fractions import Fraction
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
CHIRAL_BISPECTRUM = HERE / "bridge_strict_alpha_c12_chiral_bispectrum_orientation_report.json"
OUT_JSON = HERE / "bridge_strict_alpha_c12_source_translation_no_go_report.json"
OUT_MD = HERE / "bridge_strict_alpha_c12_source_translation_no_go_report.md"

Z12_NODE_COUNT = 12
SUPPORT_SIZE = 5
DISTANCE_SELECTED = 5
TARGET_BINARY_EXPONENT = 8
DENOMINATOR = 3
STRICT_TARGET_ETA = 9.0 / 5.0
NAD12_SUPPORT_SIZE = 12
ALPHA_SCALE = 4.0
FORWARD_ASSIGNMENT = (2, 2, 2, 1, 1)
ROUND_DIGITS = 9


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def eta_from_product(product: Fraction, branch_count: int) -> float:
    correction = float(product) ** (1.0 / branch_count)
    return math.log(NAD12_SUPPORT_SIZE * correction) / math.log(ALPHA_SCALE)


def ordered_d5_support(source: int, orientation: int) -> tuple[int, ...]:
    if orientation not in {-1, 1}:
        raise ValueError("orientation must be +/-1")
    step = (orientation * DISTANCE_SELECTED) % Z12_NODE_COUNT
    return tuple((source + index * step) % Z12_NODE_COUNT for index in range(SUPPORT_SIZE))


def value_configuration(source: int, orientation: int, assignment: tuple[int, ...] = FORWARD_ASSIGNMENT) -> tuple[int, ...]:
    values = [0] * Z12_NODE_COUNT
    for node, value in zip(ordered_d5_support(source, orientation), assignment):
        values[node] = value
    return tuple(values)


def translate_config(config: tuple[int, ...], shift: int) -> tuple[int, ...]:
    out = [0] * Z12_NODE_COUNT
    for node, value in enumerate(config):
        out[(node + shift) % Z12_NODE_COUNT] = value
    return tuple(out)


def dft(config: tuple[int, ...], mode: int) -> complex:
    return sum(value * cmath.exp(-2j * math.pi * mode * node / Z12_NODE_COUNT) for node, value in enumerate(config))


def bispectrum(config: tuple[int, ...], left_mode: int, right_mode: int) -> complex:
    return dft(config, left_mode) * dft(config, right_mode) * dft(config, (left_mode + right_mode) % Z12_NODE_COUNT).conjugate()


def rounded(value: float) -> float:
    result = round(value, ROUND_DIGITS)
    return 0.0 if result == -0.0 else result


def rounded_complex(value: complex) -> dict[str, float]:
    return {"real": rounded(value.real), "imag": rounded(value.imag)}


def phase_turns(value: complex) -> float | None:
    if abs(value) < 1e-12:
        return None
    phase = (math.atan2(value.imag, value.real) / (2.0 * math.pi)) % 1.0
    return rounded(0.0 if round(phase, ROUND_DIGITS) in {1.0, -0.0} else phase)


def power_signature(config: tuple[int, ...]) -> dict[str, float]:
    return {str(mode): rounded(abs(dft(config, mode)) ** 2) for mode in range(Z12_NODE_COUNT)}


def full_bispectrum_signature(config: tuple[int, ...]) -> dict[str, dict[str, float]]:
    return {
        f"{left},{right}": rounded_complex(bispectrum(config, left, right))
        for left in range(Z12_NODE_COUNT)
        for right in range(Z12_NODE_COUNT)
    }


def raw_phase_signature(config: tuple[int, ...]) -> dict[str, float | None]:
    return {str(mode): phase_turns(dft(config, mode)) for mode in range(Z12_NODE_COUNT)}


def orientation_rows(orientation: int) -> list[dict[str, Any]]:
    rows = []
    for source in range(Z12_NODE_COUNT):
        config = value_configuration(source, orientation)
        rows.append(
            {
                "source": source,
                "orientation": orientation,
                "ordered_support": list(ordered_d5_support(source, orientation)),
                "value_configuration": list(config),
                "power_signature": power_signature(config),
                "full_bispectrum_signature": full_bispectrum_signature(config),
                "raw_phase_signature": raw_phase_signature(config),
            }
        )
    return rows


def unique_signature_count(rows: list[dict[str, Any]], key: str) -> int:
    return len({json.dumps(row[key], sort_keys=True) for row in rows})


def source_orbit_audit(orientation: int) -> dict[str, Any]:
    base = value_configuration(0, orientation)
    orbit = {translate_config(base, shift) for shift in range(Z12_NODE_COUNT)}
    stabilizer = [shift for shift in range(Z12_NODE_COUNT) if translate_config(base, shift) == base]
    rows = orientation_rows(orientation)
    phase_mode_1_values = sorted({row["raw_phase_signature"]["1"] for row in rows})
    phase_mode_5_values = sorted({row["raw_phase_signature"]["5"] for row in rows})
    return {
        "orientation": orientation,
        "computed_translation_orbit_size": len(orbit),
        "computed_translation_stabilizer": stabilizer,
        "computed_translation_stabilizer_size": len(stabilizer),
        "row_count": len(rows),
        "power_signature_unique_count_over_sources": unique_signature_count(rows, "power_signature"),
        "full_bispectrum_signature_unique_count_over_sources": unique_signature_count(rows, "full_bispectrum_signature"),
        "raw_phase_signature_unique_count_over_sources": unique_signature_count(rows, "raw_phase_signature"),
        "raw_phase_mode_1_unique_values": phase_mode_1_values,
        "raw_phase_mode_5_unique_values": phase_mode_5_values,
        "representative_power_signature": rows[0]["power_signature"],
        "representative_bispectrum_sample": {
            key: rows[0]["full_bispectrum_signature"][key]
            for key in ("1,1", "1,2", "1,5", "2,3", "5,5")
        },
    }


def translation_covariance_audit() -> dict[str, Any]:
    max_power_error = 0.0
    max_bispectrum_error = 0.0
    max_dft_covariance_error = 0.0
    checked_power = 0
    checked_bispectrum = 0
    checked_dft = 0
    for orientation in (-1, 1):
        base = value_configuration(0, orientation)
        base_dft = {mode: dft(base, mode) for mode in range(Z12_NODE_COUNT)}
        base_power = {mode: abs(base_dft[mode]) ** 2 for mode in range(Z12_NODE_COUNT)}
        base_bispectrum = {(left, right): bispectrum(base, left, right) for left in range(Z12_NODE_COUNT) for right in range(Z12_NODE_COUNT)}
        for shift in range(Z12_NODE_COUNT):
            shifted = translate_config(base, shift)
            for mode in range(Z12_NODE_COUNT):
                coeff = dft(shifted, mode)
                expected_coeff = cmath.exp(-2j * math.pi * mode * shift / Z12_NODE_COUNT) * base_dft[mode]
                max_dft_covariance_error = max(max_dft_covariance_error, abs(coeff - expected_coeff))
                max_power_error = max(max_power_error, abs((abs(coeff) ** 2) - base_power[mode]))
                checked_dft += 1
                checked_power += 1
            for left in range(Z12_NODE_COUNT):
                for right in range(Z12_NODE_COUNT):
                    shifted_bispectrum = bispectrum(shifted, left, right)
                    max_bispectrum_error = max(max_bispectrum_error, abs(shifted_bispectrum - base_bispectrum[(left, right)]))
                    checked_bispectrum += 1
    return {
        "checked_dft_cases": checked_dft,
        "checked_power_cases": checked_power,
        "checked_bispectrum_cases": checked_bispectrum,
        "max_dft_covariance_error": max_dft_covariance_error,
        "max_power_invariance_error": max_power_error,
        "max_bispectrum_invariance_error": max_bispectrum_error,
        "dft_translation_law": "F_k(T_s x)=exp(-2πik s/12)*F_k(x)",
        "bispectrum_translation_law": "B(k,l;T_s x)=B(k,l;x) because shift phases cancel: exp(-2πi(k+l-(k+l))s/12)=1",
    }


def main() -> None:
    chiral = load_json(CHIRAL_BISPECTRUM)
    product = Fraction(2 ** TARGET_BINARY_EXPONENT, DENOMINATOR ** SUPPORT_SIZE)
    audits = {str(orientation): source_orbit_audit(orientation) for orientation in (-1, 1)}
    covariance = translation_covariance_audit()

    report = {
        "status": "OPEN_STRICT_ALPHA_C12_SOURCE_TRANSLATION_NO_GO_NO_STRICT_SELECTOR_DISCHARGE",
        "result_kind": "SCRATCH_STRICT_ALPHA_C12_SOURCE_TRANSLATION_NO_GO_PROBE__NOT_A_THEOREM",
        "source_reports": {
            "c12_chiral_bispectrum_orientation": str(CHIRAL_BISPECTRUM.relative_to(ROOT)),
        },
        "previous_chiral_bispectrum_replay": {
            "result_kind": chiral["result_kind"],
            "orientation_separating_pair_count": chiral["c12_chiral_orientation_scan"]["orientation_separating_pair_count"],
            "source_degeneracy_minus": chiral["source_degeneracy_after_c12_chiral_phase"]["-1"]["source_counts_per_signature"],
            "source_degeneracy_plus": chiral["source_degeneracy_after_c12_chiral_phase"]["1"]["source_counts_per_signature"],
        },
        "target_identity_replay": {
            "q_power_product": f"{product.numerator}/{product.denominator}",
            "eta_residual_vs_9_5": eta_from_product(product, SUPPORT_SIZE) - STRICT_TARGET_ETA,
            "forward_assignment": list(FORWARD_ASSIGNMENT),
        },
        "finite_no_go_statement": {
            "premise": "A source selector score S is C12-translation-invariant on a fixed-orientation d5 source orbit: S(T_s x)=S(x).",
            "transitivity_certificate": "For each fixed orientation, the 12 source placements are one free C12 translation orbit with stabilizer size 1.",
            "conclusion": "Every C12-translation-invariant score is constant on the fixed-orientation source orbit and cannot choose a unique source.",
            "what_can_select": "A unique source requires an origin/phase reference, boundary, defect, or another strict source-localizing term that is not translation-invariant.",
        },
        "source_orbit_audits_by_orientation": audits,
        "translation_covariance_and_invariance_audit": covariance,
        "selector_consequence": {
            "what_is_gained_from_previous_chiral_probe": "Conditional handedness can distinguish orientation.",
            "what_this_probe_rules_out": "Complete translation-invariant Fourier power and bispectrum data cannot distinguish the source within a fixed orientation.",
            "why_raw_phase_is_not_enough": "Raw Fourier phases vary with source, but their values are defined relative to an origin/phase reference not derived here.",
        },
        "candidate_interpretation": {
            "supported_by_this_probe": True,
            "content": "After a conditional chiral orientation discriminator, the remaining source ambiguity is a C12 translation-invariant no-go.",
            "why_this_is_more_proof_like": "The transitive-orbit argument proves constancy for every C12-invariant score, and the complete power/bispectrum scan verifies the Fourier invariant class used here.",
            "why_this_is_not_enough": "The probe identifies the missing source-localizing premise; it does not derive one from strict nadsoliton geometry.",
            "status": "candidate-supported-but-source-localizer-not-derived",
        },
        "hard_limits": [
            "No identity K_legacy_ont == K_strict_gate is exported.",
            "No legacy physical-role transfer is licensed.",
            "C12-translation-invariant scores cannot choose the source representative within a fixed-orientation source orbit.",
            "Complete Fourier power and complete bispectrum signatures are translation-invariant and remain constant over the 12 source shifts.",
            "Raw Fourier phase variation requires an origin/phase reference; it is not a strict-core source selector by itself.",
            "No theorem derives the origin, boundary, defect, external phase, chirality, or source-localizing reference from strict nadsoliton geometry.",
            "No theorem derives branch count m=5, total exponent n=8, denominator 3, or binary-rescale quotient.",
            "No Aut(Z_12)/N462/QW-2191 selector obstruction is discharged by this C12 source no-go audit.",
            "No theorem derives eta=9/5 without the branch model, denominator/quotient convention, support premise, ledger selector, orientation premise, assignment premise, chirality premise, and source-localizer premise.",
            "No QW-2191 discharge and no ToE closure are claimed.",
            "No ToE closure is claimed.",
        ],
        "next_honest_step": "Look for a strict source-localizing term that breaks C12 translations, or prove no-go for a specified class of source-localizing candidates that preserve strict-core symmetries.",
    }

    OUT_JSON.write_text(json.dumps(report, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    OUT_MD.write_text(
        "# Scratch strict-alpha C12 source-translation no-go probe\n\n"
        "Status: C12 translation-invariant source no-go after chiral orientation audit; no strict selector discharge.\n\n"
        f"- Orientation -1 source orbit size/stabilizer: `{audits['-1']['computed_translation_orbit_size']}` / `{audits['-1']['computed_translation_stabilizer_size']}`.\n"
        f"- Orientation +1 source orbit size/stabilizer: `{audits['1']['computed_translation_orbit_size']}` / `{audits['1']['computed_translation_stabilizer_size']}`.\n"
        f"- Full bispectrum signature counts over sources: `-1 -> {audits['-1']['full_bispectrum_signature_unique_count_over_sources']}`, `+1 -> {audits['1']['full_bispectrum_signature_unique_count_over_sources']}`.\n"
        f"- Translation covariance max errors: DFT `{covariance['max_dft_covariance_error']:.3e}`, power `{covariance['max_power_invariance_error']:.3e}`, bispectrum `{covariance['max_bispectrum_invariance_error']:.3e}`.\n"
        f"- Target replay: `q^5={product.numerator}/{product.denominator}`, eta residual `{eta_from_product(product, SUPPORT_SIZE)-STRICT_TARGET_ETA:.3e}`.\n"
        "- Honest read: translation-invariant Fourier/bispectrum data cannot choose source; raw phase needs an origin reference.\n"
        "- No false pass: no strict source-localizer theorem, no QW-2191 discharge, no ToE closure.\n",
        encoding="utf-8",
    )
    print(OUT_JSON)
    print(OUT_MD)


if __name__ == "__main__":
    main()
