#!/usr/bin/env python3
"""
Phase 63 strict repair:
- canonical phase-only whitening,
- null-calibrated whitened cross-H estimate.
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
    exact_phase_only_whiten,
    phase_randomized_signal,
    read_strain,
    scales_for_n,
)


ROOT = Path(__file__).resolve().parent
RAW = ROOT / "raw_strain_unfiltered"
OUT = ROOT / "QW_1660_v63_Whitening_strict.json"


def main() -> None:
    gps = 1266965117
    n_samples = int(os.getenv("PHASE63_STRICT_N_SAMPLES", "524288"))
    n_trials = int(os.getenv("PHASE63_STRICT_N_TRIALS", "200"))
    seed = int(os.getenv("PHASE63_STRICT_SEED", "6301"))

    rng = np.random.default_rng(seed)

    x_h1 = read_strain(RAW, "H1", gps, n_samples=n_samples)
    x_l1 = read_strain(RAW, "L1", gps, n_samples=n_samples)
    n0 = min(len(x_h1), len(x_l1))
    x_h1 = x_h1[:n0]
    x_l1 = x_l1[:n0]
    scales = scales_for_n(n0)

    h_raw = cross_mfdfa_q0(x_h1, x_l1, scales)

    w_h1 = exact_phase_only_whiten(x_h1)
    w_l1 = exact_phase_only_whiten(x_l1)
    h_white = cross_mfdfa_q0(w_h1, w_l1, scales)

    null_white = []
    for _ in range(n_trials):
        s1 = phase_randomized_signal(w_h1, rng)
        s2 = phase_randomized_signal(w_l1, rng)
        null_white.append(cross_mfdfa_q0(s1, s2, scales))
    null_white = np.array(null_white, dtype=float)

    m = float(np.mean(null_white))
    s = float(np.std(null_white))
    p_lo = empirical_p_lower(h_white, null_white)
    p_hi = empirical_p_upper(h_white, null_white)
    p_two = float(min(1.0, 2.0 * min(p_lo, p_hi)))
    z = float((h_white - m) / (s + 1e-12))

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "protocol": {
            "gps": gps,
            "n_samples": n0,
            "n_trials": n_trials,
            "seed": seed,
            "whitening": "exact_phase_only_whiten",
            "estimator": "cross_mfdfa_q0",
        },
        "results": {
            "h_cross_raw": float(h_raw),
            "h_cross_whitened": float(h_white),
        },
        "null_whitened_phase_randomized": {
            "mean": m,
            "std": s,
            "p_lower": p_lo,
            "p_upper": p_hi,
            "p_two_sided": p_two,
            "z_score": z,
        },
        "verdict": "WHITENED_CROSS_H_CALIBRATED",
    }
    OUT.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(out, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
