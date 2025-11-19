#!/usr/bin/env python3
# Author: Krzysztof Żuchowski

"""
BADANIE 123.1: Quark sector — Fit amplification scale

This script takes the eigenmode-based amplification scores and finds a global
scaling factor `A_scale` such that `A_i = 1 + A_scale * proj_norm[i]`
minimizes the log-ratio error for quark mass predictions.

Outputs: `report_123_quark_fit.json` with the best-fit scale and predictions.
"""

import numpy as np
import itertools
import json
from datetime import datetime
from scipy.linalg import eigh
from scipy.optimize import minimize_scalar

# --- constants (from repository) ---
VEV_GeV = 246.0
M_E_OBS = 0.5109989e-3
OCTAVES_EFFECTIVE = np.array([1,3,4,6,7,9,10,12])
WINDING_NUMBERS = np.array([
    0.015410, 0.035010, 0.448359, 0.090475,
    0.093498, 0.141299, 0.346460, 0.175617
])
ALPHA_GEO = 2.77
BETA_TORS = 0.01
OMEGA = 2*np.pi/8
PHI = 0.5236
COLOR_FACTOR = 3.0

PDG_QUARK_MASSES = {
    'u': 2.2e-3,
    'd': 4.7e-3,
    's': 95e-3,
    'c': 1.27,
    'b': 4.18,
    't': 172.76
}
QUARKS = ['u','d','s','c','b','t']

# compute c coupling from electron winding
w_e = WINDING_NUMBERS[0]
c_coupling = M_E_OBS / (w_e * VEV_GeV) if w_e>0 else 1e-4

# candidate octaves: choose top-6 by winding magnitude
sorted_idx = np.argsort(np.abs(WINDING_NUMBERS))[::-1]
candidate_idxs = list(sorted_idx[:6])

# build coupling matrix S (8x8) as in other scripts
n = len(OCTAVES_EFFECTIVE)
S = np.zeros((n,n))
for i in range(n):
    for j in range(n):
        d = abs(i-j)+1
        K = ALPHA_GEO * np.cos(OMEGA * d + PHI) / (1.0 + BETA_TORS * d)
        S[i,j] = K
# symmetrize
S = (S + S.T)/2.0

# eigen-decomposition
eigvals, eigvecs = eigh(S)
# sort descending by abs
idxs = np.argsort(-np.abs(eigvals))
eigvals = eigvals[idxs]
eigvecs = eigvecs[:,idxs]

# compute projection-based amplification raw scores for each octave index
proj_scores = []
for i in range(n):
    p = np.abs(eigvecs[i,:])**2
    score = float(np.sum(np.abs(eigvals) * p))
    proj_scores.append(score)
proj_scores = np.array(proj_scores)
# normalize proj_scores
proj_norm = proj_scores / (np.max(proj_scores) + 1e-12)

# --- Best mapping from previous analysis (123_QUARK_ANALYSIS.py) ---
# This mapping is fixed for the fitting of A_scale
BEST_MAPPING = {'u': 3, 'd': 4, 's': 5, 'c': 7, 'b': 6, 't': 2}

def calculate_error(A_scale):
    preds = {}
    for q in QUARKS:
        idx = BEST_MAPPING[q]
        w = abs(WINDING_NUMBERS[idx])
        A = 1.0 + A_scale * proj_norm[idx]
        m = w * c_coupling * VEV_GeV * A * COLOR_FACTOR
        preds[q] = m
    
    err = 0.0
    count = 0
    for q in QUARKS:
        obs = PDG_QUARK_MASSES.get(q)
        if obs and obs > 0:
            ratio = preds[q] / obs
            err += (np.log(max(ratio, 1e-6)))**2
            count += 1
    return err / max(count, 1)

def main():
    print("BADANIE 123.1: Quark sector — Fitting amplification scale")
    print(f"Date: {datetime.now().isoformat()}")
    print(f"Using VEV = {VEV_GeV} GeV, coupling c={c_coupling:.3e}")
    print(f"Fixed mapping: {BEST_MAPPING}")

    # Perform scalar minimization to find best A_scale
    # Search range for A_scale: [0.1, 1000.0]
    res = minimize_scalar(calculate_error, bounds=(0.1, 1000.0), method='bounded')
    
    best_A_scale = res.x
    min_error = res.fun

    print(f"\nBest fit A_scale: {best_A_scale:.3f}")
    print(f"Minimum error metric: {min_error:.3f}")

    # Calculate predictions with best_A_scale
    final_preds = {}
    final_amplifications = {}
    for q in QUARKS:
        idx = BEST_MAPPING[q]
        w = abs(WINDING_NUMBERS[idx])
        A = 1.0 + best_A_scale * proj_norm[idx]
        m = w * c_coupling * VEV_GeV * A * COLOR_FACTOR
        final_preds[q] = m
        final_amplifications[q] = A

    print("\nPredicted quark masses (with best-fit A_scale):")
    for q in QUARKS:
        r_pred = final_preds[q]
        r_obs = PDG_QUARK_MASSES.get(q)
        ratio = r_pred / r_obs if r_obs else np.nan
        print(f" {q}: m_pred={r_pred:.6e} GeV, ref={r_obs:.6e} GeV, ratio={ratio:.2f}, A={final_amplifications[q]:.3f}")

    # Save JSON report
    report = {
        'study': 'Badanie 123.1: Quark sector - Fit amplification scale',
        'date': datetime.now().isoformat(),
        'parameters': {
            'VEV_GeV': VEV_GeV,
            'c_coupling': float(c_coupling),
            'color_factor': COLOR_FACTOR,
            'winding_numbers': WINDING_NUMBERS.tolist(),
            'octaves': OCTAVES_EFFECTIVE.tolist(),
            'proj_norm': {int(i): float(proj_norm[i]) for i in range(n)}
        },
        'fixed_mapping': BEST_MAPPING,
        'fit_results': {
            'best_A_scale': float(best_A_scale),
            'min_error_metric': float(min_error)
        },
        'predictions_with_fit_scale': final_preds,
        'amplifications_with_fit_scale': final_amplifications
    }

    with open('report_123_quark_fit.json', 'w') as f:
        json.dump(report, f, indent=2)

    print('\nReport saved to report_123_quark_fit.json')

if __name__ == '__main__':
    main()
