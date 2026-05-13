#!/usr/bin/env python3
"""Shared policy-gate helpers for local-only rerun checkpoints."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class PolicyBand:
    delta_min: float
    delta_max: float


def in_policy_band(delta: float, band: PolicyBand) -> bool:
    return band.delta_min <= delta <= band.delta_max


def gate_decision(delta: float, band: PolicyBand) -> dict:
    allowed = in_policy_band(delta, band)
    return {
        "delta": round(float(delta), 6),
        "gate_status": "PASS_GATE_IN_BAND" if allowed else "FAIL_POLICY_BAND_VIOLATION",
        "allow_rerun": allowed,
    }
