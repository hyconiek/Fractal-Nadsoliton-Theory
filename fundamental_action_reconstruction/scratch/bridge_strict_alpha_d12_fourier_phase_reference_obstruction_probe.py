#!/usr/bin/env python3
"""Scratch probe: Fourier phase-reference obstruction on the anchored d5 orbit.

The D12-invariant no-go says invariant scores are constant on the 24 anchored
source/orientation representatives.  This probe checks the most natural next
escape route: use Fourier phase as the selector.  The finite result is sharp but
conditional.  Fourier magnitudes / powers are D12-invariant and therefore remain
constant on the full orbit, while Fourier phases are only D12-covariant: shifts
rotate phase and reflections conjugate modes.  Therefore phase can label a
source/orientation only after an origin and handedness/phase reference have been
provided.  That is useful resonance bookkeeping, but it is not a strict-core
selector theorem by itself.
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
NO_GO = HERE / "bridge_strict_alpha_d12_invariant_selector_no_go_report.json"
OUT_JSON = HERE / "bridge_strict_alpha_d12_fourier_phase_reference_obstruction_report.json"
OUT_MD = HERE / "bridge_strict_alpha_d12_fourier_phase_reference_obstruction_report.md"

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


def dihedral_action_on_config(config: tuple[int, ...], shift: int, reflect: bool) -> tuple[int, ...]:
    out = [0] * Z12_NODE_COUNT
    for node, value in enumerate(config):
        target = ((-node if reflect else node) + shift) % Z12_NODE_COUNT
        out[target] = value
    return tuple(out)


def dft(config: tuple[int, ...], mode: int) -> complex:
    return sum(value * cmath.exp(-2j * math.pi * mode * node / Z12_NODE_COUNT) for node, value in enumerate(config))


def rounded_complex(value: complex) -> dict[str, float]:
    real = round(value.real, ROUND_DIGITS)
    imag = round(value.imag, ROUND_DIGITS)
    if real == -0.0:
        real = 0.0
    if imag == -0.0:
        imag = 0.0
    return {"real": real, "imag": imag}


def rounded_power(value: complex) -> float:
    power = round((value.real * value.real) + (value.imag * value.imag), ROUND_DIGITS)
    return 0.0 if power == -0.0 else power


def phase_turns(value: complex) -> float | None:
    if abs(value) < 1e-12:
        return None
    phase = (math.atan2(value.imag, value.real) / (2.0 * math.pi)) % 1.0
    rounded = round(phase, ROUND_DIGITS)
    return 0.0 if rounded in {1.0, -0.0} else rounded


def coefficient_packet(config: tuple[int, ...]) -> dict[str, Any]:
    packet: dict[str, Any] = {}
    for mode in range(Z12_NODE_COUNT):
        coeff = dft(config, mode)
        packet[str(mode)] = {
            "coefficient": rounded_complex(coeff),
            "power": rounded_power(coeff),
            "phase_turns": phase_turns(coeff),
        }
    return packet


def all_anchored_rows() -> list[dict[str, Any]]:
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
                    "fourier_coefficients": coefficient_packet(config),
                }
            )
    return rows


def unique_counts_by_mode(rows: list[dict[str, Any]], field: str) -> dict[str, int]:
    return {
        str(mode): len({json.dumps(row["fourier_coefficients"][str(mode)][field], sort_keys=True) for row in rows})
        for mode in range(Z12_NODE_COUNT)
    }


def representative_power_spectrum(rows: list[dict[str, Any]]) -> dict[str, float]:
    first = rows[0]["fourier_coefficients"]
    return {str(mode): first[str(mode)]["power"] for mode in range(Z12_NODE_COUNT)}


def phase_samples_by_mode(rows: list[dict[str, Any]], modes: tuple[int, ...] = (1, 5, 7, 11)) -> dict[str, list[dict[str, Any]]]:
    samples: dict[str, list[dict[str, Any]]] = {}
    for mode in modes:
        by_phase: dict[str, dict[str, Any]] = {}
        for row in rows:
            phase = row["fourier_coefficients"][str(mode)]["phase_turns"]
            key = json.dumps(phase)
            by_phase.setdefault(
                key,
                {
                    "phase_turns": phase,
                    "example_source": row["source"],
                    "example_orientation": row["orientation"],
                },
            )
        samples[str(mode)] = list(by_phase.values())[:12]
    return samples


def transformation_law_audit() -> dict[str, Any]:
    base = value_configuration(0, -1)
    base_coeffs = {mode: dft(base, mode) for mode in range(Z12_NODE_COUNT)}
    max_error = 0.0
    worst_case: dict[str, Any] | None = None
    checked = 0
    for shift in range(Z12_NODE_COUNT):
        for reflect in (False, True):
            transformed = dihedral_action_on_config(base, shift, reflect)
            for mode in range(Z12_NODE_COUNT):
                actual = dft(transformed, mode)
                phase = cmath.exp(-2j * math.pi * mode * shift / Z12_NODE_COUNT)
                source_mode = (-mode) % Z12_NODE_COUNT if reflect else mode
                expected = phase * base_coeffs[source_mode]
                error = abs(actual - expected)
                checked += 1
                if error > max_error:
                    max_error = error
                    worst_case = {"shift": shift, "reflect": reflect, "mode": mode, "error": error}
    return {
        "checked_cases": checked,
        "max_error": max_error,
        "worst_case": worst_case,
        "law": "DFT_k(shift/reflection.x)=exp(-2πik shift/12)*DFT_{-k if reflected else k}(x)",
    }


def main() -> None:
    no_go = load_json(NO_GO)
    product = Fraction(2 ** TARGET_BINARY_EXPONENT, DENOMINATOR ** SUPPORT_SIZE)
    rows = all_anchored_rows()
    power_counts = unique_counts_by_mode(rows, "power")
    phase_counts = unique_counts_by_mode(rows, "phase_turns")
    transform_audit = transformation_law_audit()

    report = {
        "status": "OPEN_STRICT_ALPHA_D12_FOURIER_PHASE_REFERENCE_OBSTRUCTION_NO_STRICT_SELECTOR_DISCHARGE",
        "result_kind": "SCRATCH_STRICT_ALPHA_D12_FOURIER_PHASE_REFERENCE_OBSTRUCTION_PROBE__NOT_A_THEOREM",
        "source_reports": {
            "d12_invariant_selector_no_go": str(NO_GO.relative_to(ROOT)),
        },
        "previous_no_go_replay": {
            "result_kind": no_go["result_kind"],
            "computed_orbit_size": no_go["d12_orbit_replay"]["computed_orbit_size"],
            "computed_stabilizer_size": no_go["d12_orbit_replay"]["computed_stabilizer_size"],
            "invariant_feature_packets": no_go["invariant_feature_audit"]["unique_invariant_feature_packet_count"],
        },
        "target_identity_replay": {
            "q_power_product": f"{product.numerator}/{product.denominator}",
            "eta_residual_vs_9_5": eta_from_product(product, SUPPORT_SIZE) - STRICT_TARGET_ETA,
            "forward_assignment": list(FORWARD_ASSIGNMENT),
        },
        "fourier_orbit_scan": {
            "row_count": len(rows),
            "representative_power_spectrum_modes_0_to_11": representative_power_spectrum(rows),
            "power_unique_counts_by_mode": power_counts,
            "all_power_spectra_constant_on_orbit": all(count == 1 for count in power_counts.values()),
            "phase_unique_counts_by_mode": phase_counts,
            "nontrivial_phase_modes": [mode for mode, count in phase_counts.items() if count > 1],
            "phase_samples_by_mode": phase_samples_by_mode(rows),
        },
        "d12_covariance_law_audit": transform_audit,
        "phase_reference_obstruction": {
            "invariant_part": "Fourier powers |F_k|^2 are D12-invariant and are constant on all 24 anchored representatives.",
            "covariant_part": "Fourier phases are D12-covariant: they rotate under shifts and conjugate under reflections, exactly as the audited DFT covariance law records.",
            "selector_consequence": "A phase selector can label source/orientation only after an origin and handedness/phase reference are supplied; the phase reference is extra selector data unless derived strictly.",
        },
        "candidate_interpretation": {
            "supported_by_this_probe": True,
            "content": "Fourier resonance phase is useful labelled-channel bookkeeping, but without a strict phase reference it cannot be a D12-invariant source/orientation selector.",
            "why_this_is_more_proof_like": "The probe combines the exact DFT covariance law with a complete 24-row orbit scan showing invariant powers and covariant phases.",
            "why_this_is_not_enough": "The computation does not export an internal origin, handedness, or phase-reference theorem from strict nadsoliton geometry.",
            "status": "candidate-supported-but-phase-reference-not-derived",
        },
        "hard_limits": [
            "No identity K_legacy_ont == K_strict_gate is exported.",
            "No legacy physical-role transfer is licensed.",
            "Fourier power/magnitude data are D12-invariant and cannot uniquely select the anchored representative.",
            "Fourier phase data are D12-covariant, not D12-invariant; using phase to select source/orientation requires an extra reference frame/handedness premise.",
            "No theorem derives the origin, orientation, external phase, or handedness reference from strict nadsoliton geometry.",
            "No theorem derives branch count m=5, total exponent n=8, denominator 3, or binary-rescale quotient.",
            "No Aut(Z_12)/N462/QW-2191 selector obstruction is discharged by this phase-reference audit.",
            "No theorem derives eta=9/5 without the branch model, denominator/quotient convention, support premise, ledger selector, source/orientation premise, assignment premise, and phase-reference premise.",
            "No QW-2191 discharge and no ToE closure are claimed.",
            "No ToE closure is claimed.",
        ],
        "next_honest_step": "Either derive a strict internal phase reference/handedness source, or prove a broader no-go for a specified class of D12-covariant phase selectors without external frame data.",
    }

    OUT_JSON.write_text(json.dumps(report, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    OUT_MD.write_text(
        "# Scratch strict-alpha D12 Fourier phase-reference obstruction probe\n\n"
        "Status: Fourier phase audit for the anchored d5 orbit; no strict selector discharge.\n\n"
        f"- Orbit rows: `{len(rows)}`; all Fourier power spectra constant: `{all(count == 1 for count in power_counts.values())}`.\n"
        f"- Nontrivial phase modes: `{[mode for mode, count in phase_counts.items() if count > 1]}`.\n"
        f"- DFT covariance max error: `{transform_audit['max_error']:.3e}` over `{transform_audit['checked_cases']}` checked cases.\n"
        f"- Target replay: `q^5={product.numerator}/{product.denominator}`, eta residual `{eta_from_product(product, SUPPORT_SIZE)-STRICT_TARGET_ETA:.3e}`.\n"
        "- Honest read: Fourier phase can label source/orientation only after an origin plus handedness/phase reference is supplied.\n"
        "- No false pass: no strict phase-reference theorem, no QW-2191 discharge, no ToE closure.\n",
        encoding="utf-8",
    )
    print(OUT_JSON)
    print(OUT_MD)


if __name__ == "__main__":
    main()
