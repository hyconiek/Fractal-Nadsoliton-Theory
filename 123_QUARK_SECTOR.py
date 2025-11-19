#!/usr/bin/env python3
# Author: Krzysztof Żuchowski

"""
BADANIE 123: QUARK SECTOR FROM OCTAVE TOPOLOGY

Goal: Provide a runnable framework that maps octave winding numbers and
topological coupling to quark flavors (u,d,s,c,b,t), computes a first-pass
mass prediction using the same composite-Higgs quantization:

    m_i = |w_i| * c * <H> * A_i * color_factor

This script is intentionally conservative: it produces a transparent
calculation and a JSON report. The mapping and amplification factors
are hypotheses and intended to be refined in follow-up work.

Author: AI Research Agent
Date: 15 listopada 2025
"""

import numpy as np
import json
from datetime import datetime

# -----------------------------
# Fixed parameters (from 117/118)
# -----------------------------
VEV_GeV = 246.0
M_E_OBS = 0.5109989e-3  # GeV (electron)

# Winding numbers (absolute) from Badanie 117 (indexed same way as OCTAVES)
OCTAVES_EFFECTIVE = np.array([1,3,4,6,7,9,10,12])
WINDING_NUMBERS = np.array([
    0.015410,  # d=1
    0.035010,  # d=3
    0.448359,  # d=4
    0.090475,  # d=6
    0.093498,  # d=7
    0.141299,  # d=9
    0.346460,  # d=10
    0.175617   # d=12
])

# Observed (PDG) quark masses (MSbar approximate) in GeV for comparison
PDG_QUARK_MASSES = {
    'u': 2.2e-3,   # GeV (approx)
    'd': 4.7e-3,
    's': 95e-3,
    'c': 1.27,
    'b': 4.18,
    't': 172.76
}

# Hypothetical mapping: assign octaves (indices) to quark flavors.
# This mapping is a hypothesis; we expose it so it can be changed.
# We try to place light quarks on small windings and heavy quarks on larger ones.
QUARK_OCTAVE_MAPPING = {
    'u': 0,   # d=1 (very small winding) → light
    'd': 1,   # d=3
    's': 3,   # d=6
    'c': 6,   # d=10
    'b': 7,   # d=12
    't': 2    # d=4 (we put top on a large winding slot to test)
}

# Color factor: quarks carry color triplet — include factor ~3 as naive multiplicity
COLOR_FACTOR = 3.0

# Compute coupling constant c from electron as in Badanie 118
w_e = WINDING_NUMBERS[0]
if w_e > 0:
    c_coupling = M_E_OBS / (w_e * VEV_GeV)
else:
    c_coupling = 1.0e-4

# Simple amplification model for quarks: eigenvalue-projection style
# Here we use a lightweight heuristic: amplification grows with winding and
# proximity to the strongest eigenmode. This is a placeholder to be replaced
# by a dynamical eigenmode calculation.

def compute_quark_amplification(octave_idx):
    # Base amplification proportional to sqrt(|winding|)
    w = WINDING_NUMBERS[octave_idx]
    base = 1.0 + 10.0 * np.sqrt(abs(w))

    # Additional factor: heavier quarks have larger bound-state amplification
    # simple heuristic based on octave index distance from electron (idx 0)
    dist = abs(octave_idx - 0)
    dist_factor = 1.0 + 0.5 * dist

    return base * dist_factor


def predict_quark_masses(mapping):
    results = {}
    for q, idx in mapping.items():
        w = abs(WINDING_NUMBERS[idx])
        A = compute_quark_amplification(idx)
        m = w * c_coupling * VEV_GeV * A * COLOR_FACTOR
        results[q] = {
            'octave_index': int(idx),
            'octave_d': int(OCTAVES_EFFECTIVE[idx]),
            'winding': float(w),
            'amplification_A': float(A),
            'mass_GeV_predicted': float(m),
            'mass_GeV_observed_ref': float(PDG_QUARK_MASSES.get(q, np.nan)),
            'ratio_to_ref': float(m / PDG_QUARK_MASSES.get(q, np.nan)) if q in PDG_QUARK_MASSES else None
        }
    return results


def main():
    print("BADANIE 123: Quark sector — running predictions")
    print(f"Date: {datetime.now().isoformat()}")
    print(f"Using VEV = {VEV_GeV} GeV, coupling c={c_coupling:.3e}")

    results = predict_quark_masses(QUARK_OCTAVE_MAPPING)

    print("\nPredicted quark masses (first pass):")
    for q in ['u','d','s','c','b','t']:
        r = results[q]
        print(f" {q}: octave d={r['octave_d']:2d}, w={r['winding']:.6f}, A={r['amplification_A']:.3f}, m_pred={r['mass_GeV_predicted']:.6e} GeV, ref={r['mass_GeV_observed_ref']:.6e} GeV, ratio={r['ratio_to_ref']:.2f}")

    # Save JSON report
    report = {
        'study': 123,
        'title': 'Quark sector from octave topology',
        'date': datetime.now().isoformat(),
        'parameters': {
            'VEV_GeV': VEV_GeV,
            'c_coupling': float(c_coupling),
            'color_factor': COLOR_FACTOR,
            'winding_numbers': WINDING_NUMBERS.tolist(),
            'octaves': OCTAVES_EFFECTIVE.tolist()
        },
        'mapping': QUARK_OCTAVE_MAPPING,
        'predictions': results
    }

    with open('report_123_quark_sector.json', 'w') as f:
        json.dump(report, f, indent=2)

    print('\nReport saved to report_123_quark_sector.json')

if __name__ == '__main__':
    main()
