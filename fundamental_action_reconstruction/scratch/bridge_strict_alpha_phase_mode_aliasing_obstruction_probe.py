#!/usr/bin/env python3
"""Scratch probe: Fourier phase-mode aliasing obstruction for source recovery.

The phase-reference source selector certificate showed that coprime Fourier modes
k=1,5,7,11 recover all 12 source positions once an origin/phase reference is
supplied.  This probe makes the complementary arithmetic obstruction explicit.

For mode k, the translation law is F_k(T_s x)=exp(-2πi k s/12)F_k(x), so a
calibrated phase can determine only k*s modulo 12.  The source is unique iff k
is invertible modulo 12.  If gcd(k,12)=g>1, the phase aliases sources into g
indistinguishable candidates and only determines s modulo 12/g.  This is a
finite arithmetic obstruction, not a numerical accident.
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
PHASE_SOURCE = HERE / "bridge_strict_alpha_phase_reference_source_selector_certificate_report.json"
OUT_JSON = HERE / "bridge_strict_alpha_phase_mode_aliasing_obstruction_report.json"
OUT_MD = HERE / "bridge_strict_alpha_phase_mode_aliasing_obstruction_report.md"

Z12_NODE_COUNT = 12
SUPPORT_SIZE = 5
DISTANCE_SELECTED = 5
TARGET_BINARY_EXPONENT = 8
DENOMINATOR = 3
STRICT_TARGET_ETA = 9.0 / 5.0
NAD12_SUPPORT_SIZE = 12
ALPHA_SCALE = 4.0
FORWARD_ASSIGNMENT = (2, 2, 2, 1, 1)
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


def rounded(value: float) -> float:
    result = round(value, ROUND_DIGITS)
    return 0.0 if result in {-0.0, 1.0} else result


def phase_turns(value: complex) -> float:
    if abs(value) < 1e-12:
        raise ValueError("phase is undefined for a near-zero Fourier coefficient")
    return rounded((math.atan2(value.imag, value.real) / (2.0 * math.pi)) % 1.0)


def phase_reference(orientation: int, mode: int) -> float:
    return phase_turns(dft(value_configuration(0, orientation), mode))


def phase_residue_mod_12(source: int, orientation: int, mode: int) -> int:
    observed = phase_turns(dft(value_configuration(source, orientation), mode))
    reference = phase_reference(orientation, mode)
    return int(round(((reference - observed) % 1.0) * Z12_NODE_COUNT)) % Z12_NODE_COUNT


def candidate_sources_from_phase_residue(residue: int, mode: int) -> list[int]:
    mode_mod = mode % Z12_NODE_COUNT
    if mode_mod == 0:
        return list(range(Z12_NODE_COUNT))
    divisor = math.gcd(mode_mod, Z12_NODE_COUNT)
    if residue % divisor != 0:
        return []
    reduced_mode = mode_mod // divisor
    reduced_residue = residue // divisor
    reduced_modulus = Z12_NODE_COUNT // divisor
    inverse = pow(reduced_mode, -1, reduced_modulus)
    source_mod_class = (inverse * reduced_residue) % reduced_modulus
    return [source_mod_class + step * reduced_modulus for step in range(divisor)]


def mode_orientation_rows(mode: int, orientation: int) -> list[dict[str, Any]]:
    rows = []
    for source in range(Z12_NODE_COUNT):
        residue = phase_residue_mod_12(source, orientation, mode)
        candidates = candidate_sources_from_phase_residue(residue, mode)
        rows.append(
            {
                "source": source,
                "orientation": orientation,
                "mode": mode,
                "phase_residue_k_times_source_mod_12": residue,
                "candidate_sources_from_phase": candidates,
                "candidate_count": len(candidates),
                "actual_source_in_candidates": source in candidates,
                "unique_recovery": candidates == [source],
            }
        )
    return rows


def mode_aliasing_summary(mode: int, orientation: int) -> dict[str, Any]:
    rows = mode_orientation_rows(mode, orientation)
    divisor = math.gcd(mode, Z12_NODE_COUNT)
    expected_phase_classes = Z12_NODE_COUNT // divisor if mode != 0 else 1
    unique_residues = sorted({row["phase_residue_k_times_source_mod_12"] for row in rows})
    unique_candidate_sets = sorted({tuple(row["candidate_sources_from_phase"]) for row in rows})
    return {
        "mode": mode,
        "orientation": orientation,
        "gcd_mode_12": divisor,
        "expected_phase_class_count": expected_phase_classes,
        "expected_alias_size": divisor,
        "observed_phase_residue_count": len(unique_residues),
        "observed_phase_residues": unique_residues,
        "unique_candidate_set_count": len(unique_candidate_sets),
        "candidate_set_sizes": sorted({len(candidate_set) for candidate_set in unique_candidate_sets}),
        "all_actual_sources_in_candidate_sets": all(row["actual_source_in_candidates"] for row in rows),
        "all_sources_uniquely_recovered": all(row["unique_recovery"] for row in rows),
        "rows": rows,
    }


def all_mode_aliasing_summary() -> dict[str, Any]:
    return {
        str(orientation): {str(mode): mode_aliasing_summary(mode, orientation) for mode in range(Z12_NODE_COUNT)}
        for orientation in (-1, 1)
    }


def main() -> None:
    phase_source = load_json(PHASE_SOURCE)
    product = Fraction(2 ** TARGET_BINARY_EXPONENT, DENOMINATOR ** SUPPORT_SIZE)
    summaries = all_mode_aliasing_summary()
    source_complete_modes = [
        mode for mode in range(Z12_NODE_COUNT)
        if all(summaries[str(orientation)][str(mode)]["all_sources_uniquely_recovered"] for orientation in (-1, 1))
    ]
    aliasing_modes = [mode for mode in range(Z12_NODE_COUNT) if mode not in source_complete_modes]
    expected_counts_by_mode = {str(mode): (Z12_NODE_COUNT // math.gcd(mode, Z12_NODE_COUNT) if mode != 0 else 1) for mode in range(Z12_NODE_COUNT)}

    report = {
        "status": "OPEN_STRICT_ALPHA_PHASE_MODE_ALIASING_OBSTRUCTION_NO_STRICT_SELECTOR_DISCHARGE",
        "result_kind": "SCRATCH_STRICT_ALPHA_PHASE_MODE_ALIASING_OBSTRUCTION_PROBE__NOT_A_THEOREM",
        "source_reports": {
            "phase_reference_source_selector_certificate": str(PHASE_SOURCE.relative_to(ROOT)),
        },
        "previous_phase_source_certificate_replay": {
            "result_kind": phase_source["result_kind"],
            "all_rows_recovered_by_all_coprime_modes": phase_source["conditional_selector_scan"]["all_rows_recovered_by_all_coprime_modes"],
            "source_recovery_success_by_mode": phase_source["conditional_selector_scan"]["source_recovery_success_by_mode"],
            "unique_predicted_orientation_source_pairs_using_mode_1": phase_source["conditional_selector_scan"]["unique_predicted_orientation_source_pairs_using_mode_1"],
        },
        "target_identity_replay": {
            "q_power_product": f"{product.numerator}/{product.denominator}",
            "eta_residual_vs_9_5": eta_from_product(product, SUPPORT_SIZE) - STRICT_TARGET_ETA,
            "forward_assignment": list(FORWARD_ASSIGNMENT),
        },
        "aliasing_theorem_statement": {
            "translation_law": "F_k(T_s x)=exp(-2πik s/12)*F_k(x)",
            "phase_observable": "A calibrated phase observes k*s modulo 12.",
            "unique_source_condition": "A single mode recovers source uniquely iff gcd(k,12)=1.",
            "aliasing_condition": "If gcd(k,12)=g>1, the phase leaves g candidate sources and only determines s modulo 12/g.",
        },
        "phase_mode_aliasing_scan": {
            "summaries_by_orientation_and_mode": summaries,
            "expected_phase_class_count_by_mode": expected_counts_by_mode,
            "source_complete_modes": source_complete_modes,
            "aliasing_modes": aliasing_modes,
            "coprime_modes_mod_12": [mode for mode in range(Z12_NODE_COUNT) if math.gcd(mode, Z12_NODE_COUNT) == 1],
            "all_observed_counts_match_gcd_formula": all(
                summaries[str(orientation)][str(mode)]["observed_phase_residue_count"] == expected_counts_by_mode[str(mode)]
                for orientation in (-1, 1)
                for mode in range(Z12_NODE_COUNT)
            ),
            "all_alias_sizes_match_gcd_formula": all(
                summaries[str(orientation)][str(mode)]["candidate_set_sizes"] == [math.gcd(mode, Z12_NODE_COUNT)]
                for orientation in (-1, 1)
                for mode in range(Z12_NODE_COUNT)
            ),
        },
        "selector_consequence": {
            "what_is_refined": "The previous phase-origin certificate is minimal at the single-mode level: the mode must be coprime to 12.",
            "what_fails_for_non_coprime_modes": "Non-coprime modes are phase-aliased and leave multiple source candidates even with a calibrated origin.",
            "why_this_is_not_a_strict_selector": "The probe assumes a phase origin/reference; it only classifies which Fourier modes would work after that premise is supplied.",
        },
        "candidate_interpretation": {
            "supported_by_this_probe": True,
            "content": "Coprime Fourier modes are exactly the source-complete phase modes; all other modes have gcd-controlled aliasing.",
            "why_this_is_more_proof_like": "The probe states and exhaustively verifies the modular arithmetic criterion gcd(k,12)=1 across every source/orientation row.",
            "why_this_is_not_enough": "The phase origin and handedness premises remain external; the probe does not derive them from strict nadsoliton geometry.",
            "status": "candidate-supported-but-phase-origin-still-not-derived",
        },
        "hard_limits": [
            "No identity K_legacy_ont == K_strict_gate is exported.",
            "No legacy physical-role transfer is licensed.",
            "This is a phase-mode aliasing classification under an explicit phase-origin premise, not a strict-core selector theorem.",
            "Non-coprime modes do not rescue source selection; they leave gcd(k,12) aliased candidates.",
            "The phase origin/reference and handedness premises are still imported selector data and are not derived here.",
            "No theorem derives the origin, boundary, defect, external phase, chirality, handedness, or source-localizing reference from strict nadsoliton geometry.",
            "No theorem derives branch count m=5, total exponent n=8, denominator 3, or binary-rescale quotient.",
            "No Aut(Z_12)/N462/QW-2191 selector obstruction is discharged by this phase-mode aliasing audit.",
            "No theorem derives eta=9/5 without the branch model, denominator/quotient convention, support premise, ledger selector, orientation premise, assignment premise, chirality premise, and phase-origin premise.",
            "No QW-2191 discharge and no ToE closure are claimed.",
            "No ToE closure is claimed.",
        ],
        "next_honest_step": "Use the gcd(k,12)=1 criterion as a constraint on any proposed strict phase-origin term, or prove that no admissible strict-core source can privilege a coprime phase mode.",
    }

    OUT_JSON.write_text(json.dumps(report, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    OUT_MD.write_text(
        "# Scratch strict-alpha phase-mode aliasing obstruction probe\n\n"
        "Status: phase-mode aliasing classification under phase-origin premise; no strict selector discharge.\n\n"
        f"- Source-complete modes: `{source_complete_modes}`; aliasing modes: `{aliasing_modes}`.\n"
        f"- Expected phase class counts by mode: `{expected_counts_by_mode}`.\n"
        f"- All observed counts match gcd formula: `{report['phase_mode_aliasing_scan']['all_observed_counts_match_gcd_formula']}`.\n"
        f"- Target replay: `q^5={product.numerator}/{product.denominator}`, eta residual `{eta_from_product(product, SUPPORT_SIZE)-STRICT_TARGET_ETA:.3e}`.\n"
        "- Honest read: a single phase mode recovers source iff gcd(k,12)=1; non-coprime modes alias sources.\n"
        "- No false pass: no strict phase-origin theorem, no QW-2191 discharge, no ToE closure.\n",
        encoding="utf-8",
    )
    print(OUT_JSON)
    print(OUT_MD)


if __name__ == "__main__":
    main()
