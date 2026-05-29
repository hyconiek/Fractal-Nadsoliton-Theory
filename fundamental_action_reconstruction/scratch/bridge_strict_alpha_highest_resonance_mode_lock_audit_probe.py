#!/usr/bin/env python3
"""Scratch probe: highest-resonance mode-lock audit for d5 support.

The user-level physics intuition says that the nadsoliton tends toward highest
resonance.  This probe makes that statement finite and checkable on Z_12 for
5-node supports.

Honest finite result:
- If the resonance channel is explicitly locked to Fourier mode k=5 (or its
  conjugate k=7), maximizing coherent single-mode power selects exactly the
  fifth-step supports {a, a+5, a+10, a+15, a+20} mod 12, up to translation.
- If "highest resonance" is left unqualified over all modes k=1..6, the global
  winner is the Nyquist/parity mode k=6, not the d5 orbit.
- If one maximizes the full coprime-mode packet {1,5,7,11}, the maximizer set is
  broad and does not select d5.

Thus highest resonance is a useful selector only after a mode-lock/source
premise chooses the fifth resonance channel.  This is consistent with the
ontology guardrail: the nadsoliton can be treated as primordial information
that self-records its formation channel, but this probe does not derive that
self-recorded k=5 channel from strict geometry.
"""
from __future__ import annotations

import json
import math
from itertools import combinations
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
OUT_JSON = HERE / "bridge_strict_alpha_highest_resonance_mode_lock_audit_report.json"
OUT_MD = HERE / "bridge_strict_alpha_highest_resonance_mode_lock_audit_report.md"

N = 12
SUPPORT_SIZE = 5
TARGET_MODE = 5
CONJUGATE_TARGET_MODE = N - TARGET_MODE
COPRIME_MODES = [1, 5, 7, 11]
INDEPENDENT_MODES = [1, 2, 3, 4, 5, 6]
FLOAT_TOL = 1e-9
TARGET_Q_POWER = "256/243"
TARGET_ETA = "9/5"


def support_list() -> list[tuple[int, ...]]:
    return list(combinations(range(N), SUPPORT_SIZE))


def canonical_support(support: tuple[int, ...] | list[int]) -> list[int]:
    return sorted(int(x) % N for x in support)


def fifth_step_supports() -> list[list[int]]:
    supports = {
        tuple(canonical_support([(start + TARGET_MODE * step) % N for step in range(SUPPORT_SIZE)]))
        for start in range(N)
    }
    return [list(row) for row in sorted(supports)]


def contiguous_supports() -> list[list[int]]:
    supports = {
        tuple(canonical_support([(start + step) % N for step in range(SUPPORT_SIZE)]))
        for start in range(N)
    }
    return [list(row) for row in sorted(supports)]


def coherent_sum_components(support: tuple[int, ...], mode: int) -> tuple[float, float]:
    real = sum(math.cos(2.0 * math.pi * mode * node / N) for node in support)
    imag = sum(math.sin(2.0 * math.pi * mode * node / N) for node in support)
    return real, imag


def coherent_power(support: tuple[int, ...], mode: int) -> float:
    real, imag = coherent_sum_components(support, mode)
    return real * real + imag * imag


def rounded(value: float) -> float:
    return round(value, 12)


def mode_scan(mode: int, supports: list[tuple[int, ...]]) -> dict[str, Any]:
    rows = [(coherent_power(support, mode), support) for support in supports]
    maximum = max(power for power, _support in rows)
    winners = [canonical_support(support) for power, support in rows if abs(power - maximum) <= FLOAT_TOL]
    return {
        "mode": mode,
        "maximum_power": rounded(maximum),
        "winner_count": len(winners),
        "winners": winners,
    }


def packet_scan(name: str, modes: list[int], supports: list[tuple[int, ...]]) -> dict[str, Any]:
    rows = [(sum(coherent_power(support, mode) for mode in modes), support) for support in supports]
    maximum = max(power for power, _support in rows)
    winners = [canonical_support(support) for power, support in rows if abs(power - maximum) <= FLOAT_TOL]
    return {
        "packet": name,
        "modes": modes,
        "maximum_power_sum": rounded(maximum),
        "winner_count": len(winners),
        "winners": winners,
    }


def global_single_mode_scan(mode_scans: list[dict[str, Any]]) -> dict[str, Any]:
    maximum = max(scan["maximum_power"] for scan in mode_scans)
    winning_modes = [scan["mode"] for scan in mode_scans if abs(scan["maximum_power"] - maximum) <= FLOAT_TOL]
    winning_supports = []
    for scan in mode_scans:
        if scan["mode"] in winning_modes:
            for support in scan["winners"]:
                if support not in winning_supports:
                    winning_supports.append(support)
    return {
        "modes_checked": [scan["mode"] for scan in mode_scans],
        "global_maximum_power": maximum,
        "winning_modes": winning_modes,
        "winner_count": len(winning_supports),
        "winners": sorted(winning_supports),
        "is_d5_selector": sorted(winning_supports) == fifth_step_supports(),
    }


def support_set(rows: list[list[int]]) -> set[tuple[int, ...]]:
    return {tuple(row) for row in rows}


def main() -> None:
    supports = support_list()
    d5_supports = fifth_step_supports()
    contiguous = contiguous_supports()
    scans = [mode_scan(mode, supports) for mode in INDEPENDENT_MODES]
    scan_by_mode = {scan["mode"]: scan for scan in scans}
    target_scan = scan_by_mode[TARGET_MODE]
    conjugate_scan = mode_scan(CONJUGATE_TARGET_MODE, supports)
    global_scan = global_single_mode_scan(scans)
    coprime_packet = packet_scan("coprime_modes_1_5_7_11", COPRIME_MODES, supports)

    target_winners = support_set(target_scan["winners"])
    d5_set = support_set(d5_supports)
    contiguous_set = support_set(contiguous)
    mode1_set = support_set(scan_by_mode[1]["winners"])
    coprime_set = support_set(coprime_packet["winners"])

    report: dict[str, Any] = {
        "result_kind": "SCRATCH_STRICT_ALPHA_HIGHEST_RESONANCE_MODE_LOCK_AUDIT_PROBE__NOT_A_THEOREM",
        "status": "candidate-supported-only-with-fifth-mode-lock",
        "finite_model": {
            "ring": "Z_12",
            "support_size": SUPPORT_SIZE,
            "support_count": len(supports),
            "coherent_power": "P_k(S)=|sum_{x in S} exp(2*pi*i*k*x/12)|^2",
            "target_q_power": TARGET_Q_POWER,
            "target_eta": TARGET_ETA,
        },
        "mode_locked_certificate": {
            "target_mode": TARGET_MODE,
            "conjugate_target_mode": CONJUGATE_TARGET_MODE,
            "target_mode_scan": target_scan,
            "conjugate_target_mode_scan": conjugate_scan,
            "fifth_step_supports": d5_supports,
            "target_mode_winners_equal_fifth_step_supports": target_winners == d5_set,
            "conjugate_mode_winners_equal_fifth_step_supports": support_set(conjugate_scan["winners"]) == d5_set,
            "interpretation": "A fifth-mode resonance lock turns maximum coherent power into a finite d5 support selector, up to translation/conjugate orientation.",
        },
        "unqualified_highest_resonance_audit": {
            "single_mode_scans_k_1_to_6": scans,
            "global_single_mode_scan": global_scan,
            "mode_1_winners_are_contiguous_supports": mode1_set == contiguous_set,
            "mode_5_winners_are_d5_supports": target_winners == d5_set,
            "unqualified_highest_resonance_selects_d5": global_scan["is_d5_selector"],
            "obstruction": "Without a fifth-mode lock, the largest single-mode coherence over k=1..6 is the k=6 parity/Nyquist resonance, not the d5 orbit.",
        },
        "coprime_packet_audit": {
            "packet_scan": coprime_packet,
            "contains_all_d5_supports": d5_set.issubset(coprime_set),
            "winner_count": coprime_packet["winner_count"],
            "is_unique_d5_selector": coprime_set == d5_set,
            "obstruction": "Maximizing the whole coprime packet {1,5,7,11} leaves many supports tied and does not isolate the fifth-step orbit.",
        },
        "candidate_interpretation": {
            "honest_gain": "The phrase highest resonance can be made into a d5 support selector if it means highest coherent power in the self-recorded fifth Fourier channel.",
            "honest_limit": "The finite scan also shows that unqualified highest resonance is too broad/ambiguous and can select the Nyquist parity channel instead.",
            "ontology_guardrail": "The nadsoliton is treated as primordial information that may contain information about its own formation/existence; the fifth-mode record is not introduced as a deeper layer below the nadsoliton.",
            "selector_status": "conditional finite resonance certificate, not QW-2191 discharge",
        },
        "hard_limits": [
            "No identity K_legacy_ont == K_strict_gate is asserted or used.",
            "No legacy physical-role claim is transferred onto K_strict_gate.",
            "No theorem derives the fifth-mode resonance lock from strict nadsoliton geometry.",
            "No theorem derives source, phase, orientation, branch count m=5, n=8, denominator 3, or ledger selector from this support scan.",
            "Unqualified highest resonance is explicitly not enough: it selects the k=6 parity/Nyquist channel in this finite audit.",
            "No Aut(Z_12)/N462/QW-2191 selector obstruction is discharged.",
            "No QW-2191 discharge and no ToE closure are claimed.",
            "No ToE closure is claimed.",
        ],
        "next_honest_step": "Try to derive a self-recorded fifth-mode lock from a strict nadsoliton formation/existence record, or prove that highest-resonance selection necessarily needs an extra mode-lock premise.",
    }

    OUT_JSON.write_text(json.dumps(report, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    OUT_MD.write_text(
        "# Scratch strict-alpha highest-resonance mode-lock audit probe\n\n"
        "Status: d5 support selected only under an explicit fifth-mode resonance lock; no strict selector discharge.\n\n"
        f"- Supports scanned: `{len(supports)}`.\n"
        f"- Mode `k=5` winner count: `{target_scan['winner_count']}`; equals d5 supports: `{target_winners == d5_set}`.\n"
        f"- Conjugate mode `k=7` winner count: `{conjugate_scan['winner_count']}`; equals d5 supports: `{support_set(conjugate_scan['winners']) == d5_set}`.\n"
        f"- Unqualified single-mode highest resonance winning modes: `{global_scan['winning_modes']}`; selects d5: `{global_scan['is_d5_selector']}`.\n"
        f"- Coprime packet winner count: `{coprime_packet['winner_count']}`; unique d5 selector: `{coprime_set == d5_set}`.\n"
        f"- Target replay: `q^5={TARGET_Q_POWER}`, eta `{TARGET_ETA}`.\n"
        "- Honest read: highest resonance supports d5 only after a fifth-mode/source lock is supplied.\n"
        "- No false pass: no fifth-mode theorem, no QW-2191 discharge, no ToE closure.\n",
        encoding="utf-8",
    )
    print(OUT_JSON)
    print(OUT_MD)


if __name__ == "__main__":
    main()
