#!/usr/bin/env python3
# Author: Krzysztof Żuchowski

"""
BADANIE 123: Quark sector analysis — permutation search + eigenmode amplifications

Performs:
 - Permutation search mapping quark flavors to candidate octaves (top6 by |winding|)
 - Computes amplification factors A_i from eigenmode projections of the coupling matrix S
 - Evaluates predicted masses and reports best mappings for a few amplification scale factors

Outputs: `report_123_quark_analysis.json`
"""

import numpy as np
import itertools
import json
from datetime import datetime
from scipy.linalg import eigh

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
# pick first 6 indices
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
    # projection onto top modes weighted by eigenvalue magnitude
    p = np.abs(eigvecs[i,:])**2  # fraction of that octave in each eigenmode
    score = float(np.sum(np.abs(eigvals) * p))
    proj_scores.append(score)
proj_scores = np.array(proj_scores)
# normalize proj_scores
proj_norm = proj_scores / (np.max(proj_scores) + 1e-12)

# amplification model: A_i = 1 + scale * proj_norm[i]
scale_candidates = [5.0, 10.0, 20.0, 50.0]

# Permutation search among candidate_idxs for mapping to 6 quarks
perms = list(itertools.permutations(candidate_idxs, 6))
print(f"Candidate octave indices (top6 by |winding|): {candidate_idxs}")
print(f"Total permutations: {len(perms)}")

results = []
for scale in scale_candidates:
    best = None
    best_metric = 1e18
    for perm in perms:
        # construct mapping quark->octave_idx
        mapping = {q: perm[i] for i,q in enumerate(QUARKS)}
        # compute predicted masses
        preds = {}
        for q in QUARKS:
            idx = mapping[q]
            w = abs(WINDING_NUMBERS[idx])
            A = 1.0 + scale * proj_norm[idx]
            m = w * c_coupling * VEV_GeV * A * COLOR_FACTOR
            preds[q] = m
        # compute error metric: sum((log(pred/obs))^2) for quarks with obs>0
        err = 0.0
        count = 0
        for q in QUARKS:
            obs = PDG_QUARK_MASSES.get(q)
            if obs and obs>0:
                ratio = preds[q]/obs
                # penalize huge mismatches but avoid zero issues
                err += (np.log(max(ratio,1e-6)))**2
                count += 1
        metric = err / max(count,1)
        if metric < best_metric:
            best_metric = metric
            best = {
                'scale': scale,
                'metric': float(metric),
                'mapping': {q: int(mapping[q]) for q in QUARKS},
                'predictions': {q: float(preds[q]) for q in QUARKS}
            }
    results.append(best)

# select overall best by metric
overall_best = min(results, key=lambda r: r['metric'])

# also produce a mapping where t is assigned to the maximum winding index forcefully
max_w_idx = int(sorted_idx[0])

# prepare final report
report = {
    'date': datetime.now().isoformat(),
    'c_coupling': float(c_coupling),
    'candidate_idxs': [int(i) for i in candidate_idxs],
    'proj_scores': {int(i): float(proj_scores[i]) for i in range(n)},
    'proj_norm': {int(i): float(proj_norm[i]) for i in range(n)},
    'scale_candidates': scale_candidates,
    'results_per_scale': results,
    'overall_best': overall_best,
    'notes': (
        'Amplifications computed from eigenmode projections of S. ' 
        'Search restricted to top-6 octaves by |winding| (720 permutations). '
        'This is a first-pass analysis; further physics (QCD renormalization, '
        'color dynamics, or alternative amplification mechanisms) may be required for top/b mass scale.'
    )
}

with open('report_123_quark_analysis.json','w') as f:
    json.dump(report,f,indent=2)

print('\nAnalysis complete. Report saved to report_123_quark_analysis.json')
print('Best overall mapping (scale, metric):', overall_best['scale'], overall_best['metric'])
print('Mapping (quark -> octave_idx):', overall_best['mapping'])
print('Predictions (GeV):')
for q in QUARKS:
    print(f"  {q}: {overall_best['predictions'][q]:.6e} | ref {PDG_QUARK_MASSES[q]:.6e}")
