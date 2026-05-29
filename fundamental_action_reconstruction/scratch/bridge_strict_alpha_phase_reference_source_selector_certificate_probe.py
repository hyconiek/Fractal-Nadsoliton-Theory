#!/usr/bin/env python3
"""Scratch probe: phase-reference source selector certificate after C12 no-go.

The C12 source-translation no-go showed that translation-invariant Fourier power
and bispectrum data cannot choose a source.  This probe records the exact
conditional escape hatch: if a phase origin/reference is supplied, raw Fourier
phase is not translation-invariant and can recover the source.

For a fixed orientation and nonzero Fourier mode k coprime to 12,
F_k(T_s x)=exp(-2πi k s/12)F_k(x).  Therefore, after calibrating the source-0
phase for that orientation, the phase difference recovers s modulo 12.  A
handedness convention can first choose orientation via the chiral bispectrum sign,
and then a phase-origin convention can recover source.  This is a constructive
selector certificate, but it is explicitly premise-based: the origin/phase
reference and handedness are not derived from strict nadsoliton geometry here.
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
SOURCE_NO_GO = HERE / "bridge_strict_alpha_c12_source_translation_no_go_report.json"
OUT_JSON = HERE / "bridge_strict_alpha_phase_reference_source_selector_certificate_report.json"
OUT_MD = HERE / "bridge_strict_alpha_phase_reference_source_selector_certificate_report.md"

Z12_NODE_COUNT = 12
SUPPORT_SIZE = 5
DISTANCE_SELECTED = 5
TARGET_BINARY_EXPONENT = 8
DENOMINATOR = 3
STRICT_TARGET_ETA = 9.0 / 5.0
NAD12_SUPPORT_SIZE = 12
ALPHA_SCALE = 4.0
FORWARD_ASSIGNMENT = (2, 2, 2, 1, 1)
COPRIME_SOURCE_MODES = (1, 5, 7, 11)
ORIENTATION_BISPECTRUM_PAIR = (1, 5)
ROUND_DIGITS = 12


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


def dft(config: tuple[int, ...], mode: int) -> complex:
    return sum(value * cmath.exp(-2j * math.pi * mode * node / Z12_NODE_COUNT) for node, value in enumerate(config))


def bispectrum(config: tuple[int, ...], left_mode: int, right_mode: int) -> complex:
    return dft(config, left_mode) * dft(config, right_mode) * dft(config, (left_mode + right_mode) % Z12_NODE_COUNT).conjugate()


def rounded(value: float) -> float:
    result = round(value, ROUND_DIGITS)
    return 0.0 if result == -0.0 else result


def phase_turns(value: complex) -> float:
    if abs(value) < 1e-12:
        raise ValueError("phase is undefined for a near-zero Fourier coefficient")
    phase = (math.atan2(value.imag, value.real) / (2.0 * math.pi)) % 1.0
    return rounded(0.0 if round(phase, ROUND_DIGITS) in {1.0, -0.0} else phase)


def signed_orientation_from_bispectrum(config: tuple[int, ...]) -> int:
    marker = bispectrum(config, *ORIENTATION_BISPECTRUM_PAIR)
    if abs(marker.imag) < 1e-10:
        raise ValueError("chiral bispectrum marker has zero imaginary part")
    return -1 if marker.imag > 0 else 1


def calibrated_phase_references() -> dict[str, dict[str, float]]:
    return {
        str(orientation): {
            str(mode): phase_turns(dft(value_configuration(0, orientation), mode))
            for mode in COPRIME_SOURCE_MODES
        }
        for orientation in (-1, 1)
    }


def recover_source_from_phase(phase: float, reference_phase: float, mode: int) -> int:
    inverse = pow(mode, -1, Z12_NODE_COUNT)
    raw_shift = ((reference_phase - phase) % 1.0) * Z12_NODE_COUNT * inverse
    return int(round(raw_shift)) % Z12_NODE_COUNT


def source_recovery_rows() -> list[dict[str, Any]]:
    references = calibrated_phase_references()
    rows = []
    for orientation in (-1, 1):
        for source in range(Z12_NODE_COUNT):
            config = value_configuration(source, orientation)
            predicted_orientation = signed_orientation_from_bispectrum(config)
            per_mode: dict[str, Any] = {}
            for mode in COPRIME_SOURCE_MODES:
                observed_phase = phase_turns(dft(config, mode))
                reference_phase = references[str(predicted_orientation)][str(mode)]
                recovered_source = recover_source_from_phase(observed_phase, reference_phase, mode)
                per_mode[str(mode)] = {
                    "observed_phase_turns": observed_phase,
                    "reference_phase_turns": reference_phase,
                    "recovered_source": recovered_source,
                    "matches_actual_source": recovered_source == source,
                }
            rows.append(
                {
                    "actual_source": source,
                    "actual_orientation": orientation,
                    "predicted_orientation_from_chiral_bispectrum": predicted_orientation,
                    "orientation_matches_actual": predicted_orientation == orientation,
                    "ordered_support": list(ordered_d5_support(source, orientation)),
                    "per_coprime_mode_source_recovery": per_mode,
                    "all_modes_recover_source": all(packet["matches_actual_source"] for packet in per_mode.values()),
                }
            )
    return rows


def main() -> None:
    source_no_go = load_json(SOURCE_NO_GO)
    product = Fraction(2 ** TARGET_BINARY_EXPONENT, DENOMINATOR ** SUPPORT_SIZE)
    references = calibrated_phase_references()
    rows = source_recovery_rows()
    modes_all_work = {
        str(mode): all(row["per_coprime_mode_source_recovery"][str(mode)]["matches_actual_source"] for row in rows)
        for mode in COPRIME_SOURCE_MODES
    }
    predicted_pairs = {
        (row["predicted_orientation_from_chiral_bispectrum"], row["per_coprime_mode_source_recovery"]["1"]["recovered_source"])
        for row in rows
    }

    report = {
        "status": "OPEN_STRICT_ALPHA_PHASE_REFERENCE_SOURCE_SELECTOR_CERTIFICATE_PREMISE_BASED_NO_STRICT_DISCHARGE",
        "result_kind": "SCRATCH_STRICT_ALPHA_PHASE_REFERENCE_SOURCE_SELECTOR_CERTIFICATE_PROBE__NOT_A_THEOREM",
        "source_reports": {
            "c12_source_translation_no_go": str(SOURCE_NO_GO.relative_to(ROOT)),
        },
        "previous_source_no_go_replay": {
            "result_kind": source_no_go["result_kind"],
            "source_orbit_minus_size": source_no_go["source_orbit_audits_by_orientation"]["-1"]["computed_translation_orbit_size"],
            "source_orbit_plus_size": source_no_go["source_orbit_audits_by_orientation"]["1"]["computed_translation_orbit_size"],
            "full_bispectrum_minus_signature_count": source_no_go["source_orbit_audits_by_orientation"]["-1"]["full_bispectrum_signature_unique_count_over_sources"],
            "full_bispectrum_plus_signature_count": source_no_go["source_orbit_audits_by_orientation"]["1"]["full_bispectrum_signature_unique_count_over_sources"],
        },
        "target_identity_replay": {
            "q_power_product": f"{product.numerator}/{product.denominator}",
            "eta_residual_vs_9_5": eta_from_product(product, SUPPORT_SIZE) - STRICT_TARGET_ETA,
            "forward_assignment": list(FORWARD_ASSIGNMENT),
        },
        "phase_recovery_formula": {
            "translation_law": "F_k(T_s x)=exp(-2πik s/12)*F_k(x)",
            "required_modes": list(COPRIME_SOURCE_MODES),
            "formula": "s = inv(k mod 12) * 12 * (phase_ref_k(orientation)-phase_obs_k) mod 12",
            "why_coprime": "k must be invertible modulo 12 for one Fourier phase to recover all 12 sources without aliasing.",
            "calibrated_source_0_phase_references": references,
        },
        "conditional_selector_scan": {
            "row_count": len(rows),
            "rows": rows,
            "all_orientations_recovered_by_chiral_bispectrum": all(row["orientation_matches_actual"] for row in rows),
            "all_rows_recovered_by_all_coprime_modes": all(row["all_modes_recover_source"] for row in rows),
            "source_recovery_success_by_mode": modes_all_work,
            "unique_predicted_orientation_source_pairs_using_mode_1": len(predicted_pairs),
            "expected_pair_count": 2 * Z12_NODE_COUNT,
        },
        "selector_consequence": {
            "what_is_sufficient": "Handedness/chiral bispectrum plus calibrated phase origin for a coprime Fourier mode recovers orientation and source for all 24 anchored representatives.",
            "what_is_not_derived": "The handedness convention, source-0 phase calibration, and origin/phase reference are external premises in this packet.",
            "relation_to_no_go": "This does not contradict the C12 no-go because raw phase is not translation-invariant; the selector works by supplying the missing translation-breaking reference.",
        },
        "candidate_interpretation": {
            "supported_by_this_probe": True,
            "content": "A phase-origin premise is sufficient to recover source exactly after conditional chiral orientation selection.",
            "why_this_is_more_proof_like": "The probe gives the exact modular inverse recovery formula and verifies every source/orientation row for all coprime modes 1,5,7,11.",
            "why_this_is_not_enough": "No strict theorem supplies the phase origin or handedness; the certificate is premise-based rather than a strict selector discharge.",
            "status": "candidate-supported-but-phase-origin-premise-not-derived",
        },
        "hard_limits": [
            "No identity K_legacy_ont == K_strict_gate is exported.",
            "No legacy physical-role transfer is licensed.",
            "This is a sufficiency certificate under explicit handedness and phase-origin premises, not a strict-core selector theorem.",
            "The calibrated source-0 phase references are imported selector data; they are not derived from strict nadsoliton geometry here.",
            "Raw Fourier phase is not C12-translation-invariant, so using it requires breaking translation symmetry with an origin/phase reference.",
            "No theorem derives the origin, boundary, defect, external phase, chirality, handedness, or source-localizing reference from strict nadsoliton geometry.",
            "No theorem derives branch count m=5, total exponent n=8, denominator 3, or binary-rescale quotient.",
            "No Aut(Z_12)/N462/QW-2191 selector obstruction is discharged by this phase-reference source certificate.",
            "No theorem derives eta=9/5 without the branch model, denominator/quotient convention, support premise, ledger selector, orientation premise, assignment premise, chirality premise, and phase-origin premise.",
            "No QW-2191 discharge and no ToE closure are claimed.",
            "No ToE closure is claimed.",
        ],
        "next_honest_step": "Try to derive a strict internal phase-origin/source-localizing term, or prove that no admissible strict-core term can provide the needed calibrated phase reference.",
    }

    OUT_JSON.write_text(json.dumps(report, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    OUT_MD.write_text(
        "# Scratch strict-alpha phase-reference source selector certificate probe\n\n"
        "Status: premise-based phase-origin sufficiency certificate; no strict selector discharge.\n\n"
        f"- Rows checked: `{len(rows)}`; all chiral orientations recovered: `{report['conditional_selector_scan']['all_orientations_recovered_by_chiral_bispectrum']}`.\n"
        f"- All coprime modes recover all sources: `{report['conditional_selector_scan']['all_rows_recovered_by_all_coprime_modes']}` with modes `{list(COPRIME_SOURCE_MODES)}`.\n"
        f"- Unique predicted orientation/source pairs using mode 1: `{len(predicted_pairs)}` / `{2 * Z12_NODE_COUNT}`.\n"
        f"- Target replay: `q^5={product.numerator}/{product.denominator}`, eta residual `{eta_from_product(product, SUPPORT_SIZE)-STRICT_TARGET_ETA:.3e}`.\n"
        "- Honest read: phase origin plus handedness is sufficient to select source/orientation, but those are imported premises.\n"
        "- No false pass: no strict phase-origin theorem, no QW-2191 discharge, no ToE closure.\n",
        encoding="utf-8",
    )
    print(OUT_JSON)
    print(OUT_MD)


if __name__ == "__main__":
    main()
