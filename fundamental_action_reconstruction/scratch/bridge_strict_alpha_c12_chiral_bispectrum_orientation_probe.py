#!/usr/bin/env python3
"""Scratch probe: C12-invariant chiral bispectrum orientation audit.

The D12 Fourier phase-reference obstruction showed that Fourier phase is useful
only after a phase reference/handedness is supplied.  This probe asks the next
finite question: if we weaken the symmetry from D12 to rotations C12 (no fixed
origin, but a handedness is allowed), can phase information distinguish the two
d5 orientations?

Answer: yes, conditionally.  The bispectrum B(k,l)=F_k F_l conj(F_{k+l}) is
translation-invariant, so it is constant over the 12 source shifts for a fixed
orientation.  For the anchored d5 ledger, several bispectrum imaginary parts
flip sign under reflection and therefore separate the two orientation orbits.
But that is a chiral / handedness selector, not a D12-invariant selector, and it
still cannot choose the source within the 12 rotations.  Thus the computation
identifies exactly what extra premise buys: orientation only, not source, and
only after a handedness convention is imported or derived.
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
PHASE_OBSTRUCTION = HERE / "bridge_strict_alpha_d12_fourier_phase_reference_obstruction_report.json"
OUT_JSON = HERE / "bridge_strict_alpha_c12_chiral_bispectrum_orientation_report.json"
OUT_MD = HERE / "bridge_strict_alpha_c12_chiral_bispectrum_orientation_report.md"

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
SELECTED_BISPECTRUM_PAIRS = ((1, 1), (1, 2), (1, 5), (2, 3), (5, 5))


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


def rounded_complex(value: complex) -> dict[str, float]:
    return {"real": rounded(value.real), "imag": rounded(value.imag)}


def pair_key(pair: tuple[int, int]) -> str:
    return f"{pair[0]},{pair[1]}"


def all_rows() -> list[dict[str, Any]]:
    rows = []
    for source in range(Z12_NODE_COUNT):
        for orientation in (-1, 1):
            config = value_configuration(source, orientation)
            rows.append(
                {
                    "source": source,
                    "orientation": orientation,
                    "ordered_support": list(ordered_d5_support(source, orientation)),
                    "value_configuration": list(config),
                    "bispectrum": {
                        pair_key(pair): rounded_complex(bispectrum(config, *pair))
                        for pair in SELECTED_BISPECTRUM_PAIRS
                    },
                }
            )
    return rows


def orientation_bispectrum_summary(rows: list[dict[str, Any]]) -> dict[str, Any]:
    summary: dict[str, Any] = {}
    for pair in SELECTED_BISPECTRUM_PAIRS:
        key = pair_key(pair)
        by_orientation: dict[str, Any] = {}
        for orientation in (-1, 1):
            packets = {
                json.dumps(row["bispectrum"][key], sort_keys=True)
                for row in rows
                if row["orientation"] == orientation
            }
            packet = json.loads(next(iter(packets))) if len(packets) == 1 else None
            by_orientation[str(orientation)] = {
                "unique_value_count_over_12_sources": len(packets),
                "representative_value": packet,
                "imaginary_sign": 0 if packet is None or packet["imag"] == 0 else (1 if packet["imag"] > 0 else -1),
            }
        minus = by_orientation["-1"]["representative_value"]
        plus = by_orientation["1"]["representative_value"]
        summary[key] = {
            "by_orientation": by_orientation,
            "real_parts_match_under_reflection": minus is not None and plus is not None and minus["real"] == plus["real"],
            "imaginary_parts_flip_under_reflection": minus is not None and plus is not None and minus["imag"] == -plus["imag"] and minus["imag"] != 0,
            "orientation_separating": by_orientation["-1"]["imaginary_sign"] == -by_orientation["1"]["imaginary_sign"] and by_orientation["-1"]["imaginary_sign"] != 0,
        }
    return summary


def source_degeneracy_summary(rows: list[dict[str, Any]]) -> dict[str, Any]:
    orientation_classes = {}
    for orientation in (-1, 1):
        sources_by_signature: dict[str, list[int]] = {}
        for row in rows:
            if row["orientation"] != orientation:
                continue
            signature = json.dumps(row["bispectrum"], sort_keys=True)
            sources_by_signature.setdefault(signature, []).append(row["source"])
        orientation_classes[str(orientation)] = {
            "signature_count": len(sources_by_signature),
            "source_counts_per_signature": sorted(len(sources) for sources in sources_by_signature.values()),
            "representative_sources": list(sources_by_signature.values())[0],
        }
    return orientation_classes


def main() -> None:
    phase_obstruction = load_json(PHASE_OBSTRUCTION)
    product = Fraction(2 ** TARGET_BINARY_EXPONENT, DENOMINATOR ** SUPPORT_SIZE)
    rows = all_rows()
    orientation_summary = orientation_bispectrum_summary(rows)
    source_summary = source_degeneracy_summary(rows)
    separating_pairs = [key for key, packet in orientation_summary.items() if packet["orientation_separating"]]

    report = {
        "status": "OPEN_STRICT_ALPHA_C12_CHIRAL_BISPECTRUM_ORIENTATION_AUDIT_NO_STRICT_SELECTOR_DISCHARGE",
        "result_kind": "SCRATCH_STRICT_ALPHA_C12_CHIRAL_BISPECTRUM_ORIENTATION_PROBE__NOT_A_THEOREM",
        "source_reports": {
            "d12_fourier_phase_reference_obstruction": str(PHASE_OBSTRUCTION.relative_to(ROOT)),
        },
        "previous_phase_obstruction_replay": {
            "result_kind": phase_obstruction["result_kind"],
            "all_power_spectra_constant_on_orbit": phase_obstruction["fourier_orbit_scan"]["all_power_spectra_constant_on_orbit"],
            "dft_covariance_max_error": phase_obstruction["d12_covariance_law_audit"]["max_error"],
            "phase_reference_status": phase_obstruction["candidate_interpretation"]["status"],
        },
        "target_identity_replay": {
            "q_power_product": f"{product.numerator}/{product.denominator}",
            "eta_residual_vs_9_5": eta_from_product(product, SUPPORT_SIZE) - STRICT_TARGET_ETA,
            "forward_assignment": list(FORWARD_ASSIGNMENT),
        },
        "bispectrum_definition": {
            "formula": "B(k,l)=F_k * F_l * conjugate(F_{k+l mod 12})",
            "translation_invariance_reason": "Under a shift by s, F_k gains exp(-2πik s/12), so the three factors cancel exactly.",
            "reflection_behavior": "Reflection conjugates/reindexes Fourier modes; the audited chiral imaginary parts flip sign between the two d5 orientations.",
            "selected_pairs": [list(pair) for pair in SELECTED_BISPECTRUM_PAIRS],
        },
        "c12_chiral_orientation_scan": {
            "row_count": len(rows),
            "orientation_bispectrum_summary": orientation_summary,
            "orientation_separating_pair_count": len(separating_pairs),
            "orientation_separating_pairs": separating_pairs,
            "all_selected_pairs_translation_invariant_over_sources": all(
                packet["by_orientation"][str(orientation)]["unique_value_count_over_12_sources"] == 1
                for packet in orientation_summary.values()
                for orientation in (-1, 1)
            ),
            "all_selected_pairs_reflection_chiral": all(packet["imaginary_parts_flip_under_reflection"] for packet in orientation_summary.values()),
        },
        "source_degeneracy_after_c12_chiral_phase": source_summary,
        "selector_consequence": {
            "what_is_gained": "A C12-invariant but reflection-chiral bispectrum can distinguish the two d5 orientations after a handedness convention is allowed.",
            "what_remains_unselected": "The source/translation index remains 12-fold degenerate because the bispectrum is translation-invariant by construction.",
            "why_this_is_not_D12_invariant": "The separating observable is the sign of an imaginary bispectrum component, which flips under reflection.",
        },
        "candidate_interpretation": {
            "supported_by_this_probe": True,
            "content": "Chiral Fourier phase invariants can supply an orientation discriminator, but only in a C12 + handedness setting and not as a full source selector.",
            "why_this_is_more_proof_like": "The probe uses the exact bispectrum shift-cancellation identity and a complete 24-row finite scan to isolate orientation versus source selection power.",
            "why_this_is_not_enough": "No strict theorem supplies the handedness convention, and translation invariance intentionally leaves the source unresolved.",
            "status": "candidate-supported-but-handedness-and-source-not-derived",
        },
        "hard_limits": [
            "No identity K_legacy_ont == K_strict_gate is exported.",
            "No legacy physical-role transfer is licensed.",
            "The bispectrum orientation discriminator is C12-invariant but reflection-chiral; it is not D12-invariant.",
            "Using the imaginary bispectrum sign requires a handedness convention or strict chirality source that is not derived here.",
            "The source/translation representative remains 12-fold degenerate under translation-invariant bispectrum data.",
            "No theorem derives the origin, orientation, external phase, chirality, or handedness reference from strict nadsoliton geometry.",
            "No theorem derives branch count m=5, total exponent n=8, denominator 3, or binary-rescale quotient.",
            "No Aut(Z_12)/N462/QW-2191 selector obstruction is discharged by this chiral bispectrum audit.",
            "No theorem derives eta=9/5 without the branch model, denominator/quotient convention, support premise, ledger selector, source/orientation premise, assignment premise, and chiral phase-reference premise.",
            "No QW-2191 discharge and no ToE closure are claimed.",
            "No ToE closure is claimed.",
        ],
        "next_honest_step": "Either derive a strict chirality/handedness source and a separate source-localizing term, or prove a no-go for translation-invariant phase invariants selecting the source representative.",
    }

    OUT_JSON.write_text(json.dumps(report, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    OUT_MD.write_text(
        "# Scratch strict-alpha C12 chiral bispectrum orientation probe\n\n"
        "Status: C12-invariant / reflection-chiral Fourier phase audit; no strict selector discharge.\n\n"
        f"- Rows: `{len(rows)}`; separating bispectrum pairs: `{separating_pairs}`.\n"
        f"- All selected pairs translation-invariant over sources: `{report['c12_chiral_orientation_scan']['all_selected_pairs_translation_invariant_over_sources']}`.\n"
        f"- Source degeneracy per orientation: `{source_summary}`.\n"
        f"- Target replay: `q^5={product.numerator}/{product.denominator}`, eta residual `{eta_from_product(product, SUPPORT_SIZE)-STRICT_TARGET_ETA:.3e}`.\n"
        "- Honest read: chiral bispectrum phase can distinguish orientation only after handedness is allowed, and it leaves source unresolved.\n"
        "- No false pass: no strict chirality/source theorem, no QW-2191 discharge, no ToE closure.\n",
        encoding="utf-8",
    )
    print(OUT_JSON)
    print(OUT_MD)


if __name__ == "__main__":
    main()
