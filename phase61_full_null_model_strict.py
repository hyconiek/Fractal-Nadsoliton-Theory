#!/usr/bin/env python3
"""
Phase 61 strict repair:
- removes hard-coded observation,
- computes observation and nulls under one identical protocol.
"""

from __future__ import annotations

import json
import os
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from qw1660_strict_common import (
    cross_mfdfa_q0,
    empirical_p_lower,
    empirical_p_upper,
    phase_randomized_signal,
    read_strain,
    scales_for_n,
)


ROOT = Path(__file__).resolve().parent
RAW = ROOT / "raw_strain_unfiltered"
OUT = ROOT / "QW_1660_v61_FullNullModel_strict.json"


def main() -> None:
    gps = 1266965117
    n_samples = int(os.getenv("PHASE61_STRICT_N_SAMPLES", "524288"))
    n_trials = int(os.getenv("PHASE61_STRICT_N_TRIALS", "240"))
    seed = int(os.getenv("PHASE61_STRICT_SEED", "6101"))

    rng = np.random.default_rng(seed)

    x_h1 = read_strain(RAW, "H1", gps, n_samples=n_samples)
    x_l1 = read_strain(RAW, "L1", gps, n_samples=n_samples)
    n0 = min(len(x_h1), len(x_l1))
    x_h1 = x_h1[:n0]
    x_l1 = x_l1[:n0]
    scales = scales_for_n(n0)

    h_obs = cross_mfdfa_q0(x_h1, x_l1, scales)

    null_phase = []
    for _ in range(n_trials):
        s1 = phase_randomized_signal(x_h1, rng)
        s2 = phase_randomized_signal(x_l1, rng)
        null_phase.append(cross_mfdfa_q0(s1, s2, scales))
    null_phase = np.array(null_phase, dtype=float)

    m = float(np.mean(null_phase))
    s = float(np.std(null_phase))
    p_lo = empirical_p_lower(h_obs, null_phase)
    p_hi = empirical_p_upper(h_obs, null_phase)
    p_two = float(min(1.0, 2.0 * min(p_lo, p_hi)))
    z = float((h_obs - m) / (s + 1e-12))

    verdict = (
        "STRICT_NULL_REJECTED"
        if p_two <= 0.01
        else "STRICT_NULL_NOT_REJECTED"
    )

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "protocol": {
            "gps": gps,
            "n_samples": n0,
            "n_trials": n_trials,
            "seed": seed,
            "estimator": "cross_mfdfa_q0",
            "observation_computed_from_same_protocol": True,
        },
        "observation": {"h_cross_obs": float(h_obs)},
        "null_phase_randomized": {
            "mean": m,
            "std": s,
            "min": float(np.min(null_phase)),
            "max": float(np.max(null_phase)),
            "p_lower": p_lo,
            "p_upper": p_hi,
            "p_two_sided": p_two,
            "z_score": z,
        },
        "verdict": verdict,
    }
    OUT.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(out, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

