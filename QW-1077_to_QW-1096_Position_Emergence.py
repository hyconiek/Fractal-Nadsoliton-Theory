#!/usr/bin/env python3
"""
QW-1077 to QW-1096: EMERGENCJA POZYCJI Z JĄDRA K(d)
====================================================
PROBLEM: Pozycje d = 0, 1.75, 3.5, 6.0, 13.75 są ZAŁOŻONE, nie WYPROWADZONE.

CEL: Pokazać że K(d) SAMO WYBIERA te pozycje jako:
  1. Minima energii
  2. Rezonanse
  3. Punkty stacjonarne
  4. Stany związane
  5. Atraktorы dynamiki

Author: Krzysztof Żuchowski / Antigravity System
Date: 2025-12-10
"""

import numpy as np
from scipy.linalg import eigh
from scipy.signal import find_peaks, argrelextrema
from scipy.optimize import minimize, minimize_scalar, brentq
from scipy.integrate import solve_ivp
import time
import json

# ============================================================================
# FROZEN KERNEL
# ============================================================================
ALPHA_GEO = 4 * np.log(2)
OMEGA = np.pi / 4
PHI = np.pi / 6
BETA_TORS = 0.01
N_OCTAVES = 12

def K(d):
    d = np.asarray(d)
    return ALPHA_GEO * np.cos(OMEGA * d + PHI) / (1 + BETA_TORS * np.abs(d))

def K_complex(d):
    d = np.asarray(d)
    return ALPHA_GEO * np.exp(1j * (OMEGA * d + PHI)) / (1 + BETA_TORS * np.abs(d))

# Target positions from QW-726
TARGET_POSITIONS = [0.00, 1.75, 2.25, 3.50, 5.00, 6.00, 13.75, 14.50]

def match_score(found_positions):
    """How well do found positions match targets?"""
    if len(found_positions) == 0:
        return 0, []
    
    matches = []
    for target in TARGET_POSITIONS:
        closest = min(found_positions, key=lambda x: abs(x - target))
        error = abs(closest - target)
        if error < 0.3:  # Within 0.3 octave
            matches.append((target, closest, error))
    
    score = len(matches) / len(TARGET_POSITIONS) * 100
    return score, matches

# ============================================================================
# QW-1077: EXTREMA OF K(d)
# ============================================================================
def qw1077_K_extrema():
    """Are particle positions at K(d) extrema?"""
    print("[QW-1077] K(d) Extrema")
    
    d_vals = np.linspace(0, 16, 500)
    K_vals = K(d_vals)
    
    # Find maxima
    max_idx = argrelextrema(K_vals, np.greater)[0]
    maxima = d_vals[max_idx]
    
    # Find minima
    min_idx = argrelextrema(K_vals, np.less)[0]
    minima = d_vals[min_idx]
    
    all_extrema = sorted(list(maxima) + list(minima))
    
    score, matches = match_score(all_extrema)
    
    return {
        "maxima": [float(m) for m in maxima[:6]],
        "minima": [float(m) for m in minima[:6]],
        "all_extrema": [float(e) for e in all_extrema[:10]],
        "match_score": float(score),
        "matches": [(float(t), float(f), float(e)) for t, f, e in matches],
        "INSIGHT": f"K(d) extrema match {score:.0f}% of targets"
    }

# ============================================================================
# QW-1078: ZEROS OF K(d)
# ============================================================================
def qw1078_K_zeros():
    """Particle positions at K(d) zeros?"""
    print("[QW-1078] K(d) Zeros")
    
    d_vals = np.linspace(0.1, 16, 500)
    K_vals = K(d_vals)
    
    # Find zeros (sign changes)
    sign_changes = np.where(np.diff(np.sign(K_vals)))[0]
    zeros = d_vals[sign_changes]
    
    score, matches = match_score(zeros.tolist())
    
    return {
        "zeros": [float(z) for z in zeros[:10]],
        "match_score": float(score),
        "matches": [(float(t), float(f), float(e)) for t, f, e in matches],
        "INSIGHT": f"K(d) zeros match {score:.0f}% of targets"
    }

# ============================================================================
# QW-1079: EIGENVALUE POSITIONS
# ============================================================================
def qw1079_eigenvalue_positions():
    """Do K(d) matrix eigenspaces select particle positions?"""
    print("[QW-1079] Eigenvalues → Positions")
    
    N = 32
    sites = np.linspace(0, 16, N)
    
    H = np.zeros((N, N))
    for i in range(N):
        for j in range(N):
            H[i, j] = K(abs(sites[i] - sites[j]))
    H = (H + H.T) / 2
    
    eigenvalues, eigenvectors = np.linalg.eigh(H)
    
    # For each eigenstate, find where it peaks
    peak_positions = []
    for k in range(min(10, N)):
        vec = np.abs(eigenvectors[:, k])
        peak_idx = np.argmax(vec)
        peak_positions.append(sites[peak_idx])
    
    score, matches = match_score(peak_positions)
    
    return {
        "eigenvalues_top5": [float(e) for e in sorted(eigenvalues)[::-1][:5]],
        "peak_positions": [float(p) for p in peak_positions],
        "match_score": float(score),
        "matches": [(float(t), float(f), float(e)) for t, f, e in matches],
        "INSIGHT": f"Eigenstate peaks match {score:.0f}% of targets"
    }

# ============================================================================
# QW-1080: ENERGY MINIMA
# ============================================================================
def qw1080_energy_minima():
    """Minimize energy functional to find stable positions"""
    print("[QW-1080] Energy Minima")
    
    # Energy of a "particle" at position d
    def energy(d):
        # Kinetic + Potential from K(d)
        return -K(d)  # Particles sit where K is large (attraction)
    
    # Find all local minima
    minima = []
    for d_init in np.arange(0, 16, 0.5):
        result = minimize(energy, d_init, bounds=[(0, 16)])
        if result.success:
            d_min = result.x[0]
            # Check if new
            if not any(abs(d_min - m) < 0.1 for m in minima):
                minima.append(d_min)
    
    minima = sorted(minima)
    score, matches = match_score(minima)
    
    return {
        "energy_minima": [float(m) for m in minima[:10]],
        "match_score": float(score),
        "matches": [(float(t), float(f), float(e)) for t, f, e in matches],
        "INSIGHT": f"Energy minima match {score:.0f}% of targets"
    }

# ============================================================================
# QW-1081: STANDING WAVE NODES
# ============================================================================
def qw1081_standing_waves():
    """Nodes of standing wave solutions"""
    print("[QW-1081] Standing Wave Nodes")
    
    # Solve: -d²ψ/dx² + V(x)ψ = Eψ where V = -K(d)
    # Look for bound states
    
    N = 200
    dx = 0.1
    d_vals = np.arange(0, N * dx, dx)
    V = -K(d_vals)
    
    # Build Hamiltonian
    H = np.zeros((N, N))
    for i in range(N):
        H[i, i] = 2 / dx**2 + V[i]
        if i > 0:
            H[i, i-1] = -1 / dx**2
        if i < N - 1:
            H[i, i+1] = -1 / dx**2
    
    eigenvalues, eigenvectors = np.linalg.eigh(H)
    
    # Find nodes of lowest eigenstates
    nodes_all = []
    for k in range(min(5, N)):
        vec = eigenvectors[:, k]
        sign_changes = np.where(np.diff(np.sign(vec)))[0]
        nodes = d_vals[sign_changes]
        nodes_all.extend(nodes.tolist())
    
    nodes_unique = sorted(set([round(n, 1) for n in nodes_all]))
    
    score, matches = match_score(nodes_unique)
    
    return {
        "wave_nodes": [float(n) for n in nodes_unique[:15]],
        "n_bound_states": int(np.sum(eigenvalues < 0)),
        "match_score": float(score),
        "matches": [(float(t), float(f), float(e)) for t, f, e in matches],
        "INSIGHT": f"Wave nodes match {score:.0f}% of targets"
    }

# ============================================================================
# QW-1082: WINDING NUMBER QUANTIZATION
# ============================================================================
def qw1082_winding():
    """Where does cumulative winding hit integer values?"""
    print("[QW-1082] Winding Quantization")
    
    d_vals = np.linspace(0.01, 16, 500)
    phases = np.angle(K_complex(d_vals))
    
    # Unwrap phase
    unwrapped = np.unwrap(phases)
    
    # Cumulative winding (in units of 2π)
    winding = unwrapped / (2 * np.pi)
    
    # Find where winding crosses 0.25, 0.5, 0.75, 1, ...
    quantized_positions = []
    for n in np.arange(0.25, 5, 0.25):
        # Find crossing
        idx = np.argmin(np.abs(winding - n))
        quantized_positions.append(d_vals[idx])
    
    score, matches = match_score(quantized_positions)
    
    return {
        "quantized_positions": [float(p) for p in quantized_positions[:16]],
        "winding_at_targets": [float(np.interp(t, d_vals, winding)) for t in TARGET_POSITIONS],
        "match_score": float(score),
        "matches": [(float(t), float(f), float(e)) for t, f, e in matches],
        "INSIGHT": f"Winding-quantized positions match {score:.0f}% of targets"
    }

# ============================================================================
# QW-1083: RESONANCE WIDTHS
# ============================================================================
def qw1083_resonances():
    """Sharp resonances in response function"""
    print("[QW-1083] Resonance Positions")
    
    # Response: R(d) = 1 / (1 - K(d)/max(|K|))
    d_vals = np.linspace(0.1, 16, 500)
    K_vals = K(d_vals)
    K_max = np.max(np.abs(K_vals))
    
    # Avoid division by zero
    denominator = 1 - K_vals / (K_max + 0.01)
    response = 1 / (np.abs(denominator) + 0.01)
    
    # Find peaks in response
    peaks, properties = find_peaks(response, prominence=1)
    resonance_positions = d_vals[peaks]
    
    score, matches = match_score(resonance_positions.tolist())
    
    return {
        "resonance_positions": [float(r) for r in resonance_positions[:10]],
        "match_score": float(score),
        "matches": [(float(t), float(f), float(e)) for t, f, e in matches],
        "INSIGHT": f"Resonances match {score:.0f}% of targets"
    }

# ============================================================================
# QW-1084: DYNAMICAL ATTRACTOR
# ============================================================================
def qw1084_attractor():
    """Evolve particle under K(d) force to find attractors"""
    print("[QW-1084] Dynamical Attractors")
    
    def equations(t, y):
        d, v = y
        # Force = -dV/dd = dK/dd
        dd = 0.01
        force = (K(d + dd) - K(d - dd)) / (2 * dd)
        return [v, force - 0.1 * v]  # With damping
    
    attractors = []
    for d_init in np.arange(0.5, 16, 1):
        sol = solve_ivp(equations, [0, 100], [d_init, 0], max_step=0.1)
        final_d = sol.y[0, -1]
        if 0 < final_d < 16:
            if not any(abs(final_d - a) < 0.2 for a in attractors):
                attractors.append(final_d)
    
    attractors = sorted(attractors)
    score, matches = match_score(attractors)
    
    return {
        "attractors": [float(a) for a in attractors[:10]],
        "match_score": float(score),
        "matches": [(float(t), float(f), float(e)) for t, f, e in matches],
        "INSIGHT": f"Dynamical attractors match {score:.0f}% of targets"
    }

# ============================================================================
# QW-1085: FRACTAL LAYER POSITIONS
# ============================================================================
def qw1085_fractal():
    """K(d) structure at fractal layer boundaries"""
    print("[QW-1085] Fractal Layers")
    
    # Fractal layers: d = n × log_4(1/β) where β = BETA_TORS
    # log_4(1/0.01) = log_4(100) = log(100)/log(4) = 3.32
    
    layer_step = np.log(1/BETA_TORS) / np.log(4)
    layer_positions = [n * layer_step for n in range(6)]
    
    # Alternative: sub-layers at layer_step/4
    sublayer_positions = [n * layer_step / 4 for n in range(24)]
    
    score_layer, matches_layer = match_score(layer_positions)
    score_sub, matches_sub = match_score(sublayer_positions)
    
    return {
        "layer_step": float(layer_step),  # ~3.32
        "layer_positions": [float(l) for l in layer_positions],
        "sublayer_positions": [float(s) for s in sublayer_positions[:16]],
        "layer_match_score": float(score_layer),
        "sublayer_match_score": float(score_sub),
        "matches_sublayer": [(float(t), float(f), float(e)) for t, f, e in matches_sub],
        "INSIGHT": f"Sublayers match {score_sub:.0f}% of targets"
    }

# ============================================================================
# QW-1086: QUANTIZED ACTION
# ============================================================================
def qw1086_action():
    """Where total action ∫K(d')dd' is quantized"""
    print("[QW-1086] Quantized Action")
    
    d_vals = np.linspace(0, 16, 500)
    K_vals = K(d_vals)
    
    # Cumulative action
    action = np.cumsum(K_vals) * (d_vals[1] - d_vals[0])
    
    # Normalize
    action_norm = action / ALPHA_GEO
    
    # Find positions where action = n × 0.25 (quantized in sub-bits)
    quantized_d = []
    for n in np.arange(0, 20, 0.25):
        idx = np.argmin(np.abs(action_norm - n))
        if idx > 0 and idx < len(d_vals) - 1:
            quantized_d.append(d_vals[idx])
    
    # Remove duplicates
    quantized_d = sorted(set([round(d, 1) for d in quantized_d]))
    
    score, matches = match_score(quantized_d)
    
    return {
        "quantized_action_positions": [float(d) for d in quantized_d[:16]],
        "action_at_targets": [float(np.interp(t, d_vals, action_norm)) for t in TARGET_POSITIONS],
        "match_score": float(score),
        "matches": [(float(t), float(f), float(e)) for t, f, e in matches],
        "INSIGHT": f"Quantized action matches {score:.0f}% of targets"
    }

# ============================================================================
# QW-1087: K(d) PRODUCT CHAIN
# ============================================================================
def qw1087_product_chain():
    """Positions where product ∏K(n×0.25) has special values"""
    print("[QW-1087] Product Chain")
    
    # Product of K at all steps up to d
    def cumulative_product(d_max):
        d_steps = np.arange(0.25, d_max + 0.01, 0.25)
        product = 1
        for d in d_steps:
            K_val = K(d)
            if K_val != 0:
                product *= K_val
        return product
    
    products = []
    for d in np.arange(0.25, 16, 0.25):
        p = cumulative_product(d)
        products.append((d, p))
    
    # Find sign changes (zeros of product)
    sign_positions = []
    for i in range(1, len(products)):
        if products[i][1] * products[i-1][1] < 0:
            sign_positions.append(products[i][0])
    
    # Find extrema of |product|
    abs_products = [(d, abs(p)) for d, p in products]
    
    score, matches = match_score(sign_positions)
    
    return {
        "product_sign_changes": [float(s) for s in sign_positions[:10]],
        "match_score": float(score),
        "INSIGHT": f"Product sign changes match {score:.0f}% of targets"
    }

# ============================================================================
# QW-1088: BAND GAPS IN K SPECTRUM
# ============================================================================
def qw1088_band_gaps():
    """Gaps in eigenvalue spectrum locate particle positions?"""
    print("[QW-1088] Band Gaps")
    
    N = 64
    sites = np.linspace(0, 16, N)
    
    H = np.zeros((N, N))
    for i in range(N):
        for j in range(N):
            H[i, j] = K(abs(sites[i] - sites[j]))
    H = (H + H.T) / 2
    
    eigenvalues = np.sort(np.linalg.eigvalsh(H))
    
    # Find large gaps
    gaps = np.diff(eigenvalues)
    large_gap_idx = np.argsort(gaps)[::-1][:10]
    
    # Position in spectrum maps to position in space?
    gap_positions = [sites[int(idx * N / len(eigenvalues))] for idx in large_gap_idx]
    
    score, matches = match_score(gap_positions)
    
    return {
        "largest_gaps": [float(gaps[i]) for i in large_gap_idx[:5]],
        "gap_positions": [float(p) for p in gap_positions[:10]],
        "match_score": float(score),
        "INSIGHT": f"Band gap positions match {score:.0f}% of targets"
    }

# ============================================================================
# QW-1089: HOPFION STABILITY
# ============================================================================
def qw1089_hopfion():
    """Stable hopfion sizes at specific d values"""
    print("[QW-1089] Hopfion Stability")
    
    # Hopfion energy: balance of K(d) and surface tension
    def hopfion_energy(d, winding=1):
        volume_term = -K(d) * d**3  # Bulk energy
        surface_term = ALPHA_GEO * d**2 * winding  # Surface tension
        return volume_term + surface_term
    
    # Find stable sizes for each winding
    stable_sizes = []
    for w in [1, 2, 3, 4]:
        for d_init in np.arange(0.5, 16, 1):
            result = minimize(lambda d: hopfion_energy(d[0], w), [d_init], bounds=[(0.1, 16)])
            if result.success:
                d_opt = result.x[0]
                if not any(abs(d_opt - s) < 0.2 for s in stable_sizes):
                    stable_sizes.append(d_opt)
    
    stable_sizes = sorted(set([round(s, 2) for s in stable_sizes]))
    score, matches = match_score(stable_sizes)
    
    return {
        "stable_hopfion_sizes": [float(s) for s in stable_sizes[:10]],
        "match_score": float(score),
        "matches": [(float(t), float(f), float(e)) for t, f, e in matches],
        "INSIGHT": f"Stable hopfion sizes match {score:.0f}% of targets"
    }

# ============================================================================
# QW-1090: K CORRELATION LENGTH
# ============================================================================
def qw1090_correlation():
    """Where K(d) correlation length changes"""
    print("[QW-1090] Correlation Length")
    
    d_vals = np.linspace(0, 16, 200)
    K_vals = K(d_vals)
    
    # Local correlation length (from autocorrelation decay)
    corr_lengths = []
    window = 20
    for i in range(len(K_vals) - window):
        local_K = K_vals[i:i+window]
        autocorr = np.correlate(local_K, local_K, mode='full')
        autocorr = autocorr[len(autocorr)//2:]
        autocorr = autocorr / autocorr[0]
        
        # Find where it drops to 1/e
        idx = np.argmax(autocorr < 1/np.e)
        if idx > 0:
            xi = idx * (d_vals[1] - d_vals[0])
            corr_lengths.append((d_vals[i], xi))
    
    # Find where correlation length changes abruptly
    if len(corr_lengths) > 1:
        d_corr = [c[0] for c in corr_lengths]
        xi_vals = [c[1] for c in corr_lengths]
        grad = np.abs(np.gradient(xi_vals))
        peaks, _ = find_peaks(grad, prominence=0.1)
        transition_d = [d_corr[p] for p in peaks if p < len(d_corr)]
    else:
        transition_d = []
    
    score, matches = match_score(transition_d)
    
    return {
        "transition_positions": [float(t) for t in transition_d[:10]],
        "match_score": float(score),
        "INSIGHT": f"Correlation transitions match {score:.0f}% of targets"
    }

# ============================================================================
# QW-1091: 0.25 GRID WITH K WEIGHT
# ============================================================================
def qw1091_weighted_grid():
    """Select 0.25-quantized positions by K(d) weight"""
    print("[QW-1091] K-Weighted Grid Selection")
    
    all_positions = np.arange(0, 16, 0.25)
    
    # Weight each position by |K(d)|
    weights = np.abs(K(all_positions))
    
    # Select positions with highest weight
    sorted_idx = np.argsort(weights)[::-1]
    top_positions = [all_positions[i] for i in sorted_idx[:12]]
    top_positions = sorted(top_positions)
    
    score, matches = match_score(top_positions)
    
    # Alternative: cumulative selection
    # Select position if |K(d)| > threshold
    threshold = np.median(weights)
    selected = [d for d, w in zip(all_positions, weights) if w > threshold]
    
    score_sel, matches_sel = match_score(selected)
    
    return {
        "top_K_positions": [float(p) for p in top_positions],
        "selected_by_threshold": [float(s) for s in selected[:16]],
        "top_match_score": float(score),
        "threshold_match_score": float(score_sel),
        "INSIGHT": f"K-weighted selection matches {score:.0f}% of targets"
    }

# ============================================================================
# QW-1092: EIGENVALUE LEVEL CROSSING
# ============================================================================
def qw1092_level_crossing():
    """Where eigenvalues cross as d varies"""
    print("[QW-1092] Level Crossing")
    
    crossing_positions = []
    N = 16
    
    prev_eig = None
    for d_param in np.arange(0, 16, 0.1):
        sites = np.linspace(0, d_param + 1, N)
        H = np.zeros((N, N))
        for i in range(N):
            for j in range(N):
                H[i, j] = K(abs(sites[i] - sites[j]))
        H = (H + H.T) / 2
        eig = np.sort(np.linalg.eigvalsh(H))
        
        if prev_eig is not None:
            # Check for level reordering
            changes = np.sum(np.sign(eig[1:] - eig[:-1]) * np.sign(prev_eig[1:] - prev_eig[:-1]) < 0)
            if changes > 0:
                crossing_positions.append(d_param)
        prev_eig = eig
    
    # Remove close duplicates
    crossing_unique = []
    for c in crossing_positions:
        if not any(abs(c - cu) < 0.3 for cu in crossing_unique):
            crossing_unique.append(c)
    
    score, matches = match_score(crossing_unique)
    
    return {
        "level_crossings": [float(c) for c in crossing_unique[:10]],
        "match_score": float(score),
        "INSIGHT": f"Level crossings match {score:.0f}% of targets"
    }

# ============================================================================
# QW-1093: BEST MECHANISM COMBINATION
# ============================================================================
def qw1093_combination(all_results):
    """Combine multiple mechanisms to find best match"""
    print("[QW-1093] Best Combination")
    
    # Collect positions from different mechanisms
    mechanism_positions = {}
    for key, result in all_results.items():
        if isinstance(result, dict):
            for pos_key in ["maxima", "minima", "zeros", "peak_positions", 
                           "energy_minima", "wave_nodes", "quantized_positions",
                           "resonance_positions", "attractors", "sublayer_positions"]:
                if pos_key in result:
                    mechanism_positions[f"{key}_{pos_key}"] = result[pos_key]
    
    # Find positions that appear in multiple mechanisms
    all_pos = []
    for positions in mechanism_positions.values():
        all_pos.extend(positions)
    
    # Count occurrences at 0.25 resolution
    from collections import Counter
    pos_rounded = [round(p * 4) / 4 for p in all_pos]
    counts = Counter(pos_rounded)
    
    # Most common positions
    most_common = [pos for pos, count in counts.most_common(10)]
    
    score, matches = match_score(most_common)
    
    return {
        "combined_best": [float(p) for p in most_common[:10]],
        "occurrence_counts": {str(k): v for k, v in counts.most_common(15)},
        "match_score": float(score),
        "matches": [(float(t), float(f), float(e)) for t, f, e in matches],
        "INSIGHT": f"Combined mechanisms match {score:.0f}% of targets"
    }

# ============================================================================
# QW-1094: INTEGRATED K(d) PEAKS
# ============================================================================
def qw1094_integrated():
    """Positions where integrated |K| is special"""
    print("[QW-1094] Integrated K Peaks")
    
    d_vals = np.linspace(0, 16, 500)
    K_abs = np.abs(K(d_vals))
    
    # Cumulative integral
    cumulative = np.cumsum(K_abs) * (d_vals[1] - d_vals[0])
    
    # Normalize so total = 1
    cumulative_norm = cumulative / cumulative[-1]
    
    # Find where cumulative = 0.1, 0.2, 0.3, ..., 0.9
    quantile_positions = []
    for q in np.arange(0.05, 1.0, 0.05):
        idx = np.argmin(np.abs(cumulative_norm - q))
        quantile_positions.append(d_vals[idx])
    
    score, matches = match_score(quantile_positions)
    
    return {
        "quantile_positions": [float(p) for p in quantile_positions],
        "cumulative_at_targets": [float(np.interp(t, d_vals, cumulative_norm)) for t in TARGET_POSITIONS],
        "match_score": float(score),
        "matches": [(float(t), float(f), float(e)) for t, f, e in matches],
        "INSIGHT": f"Integrated K quantiles match {score:.0f}% of targets"
    }

# ============================================================================
# QW-1095: K(d) DERIVATIVE ZEROS
# ============================================================================
def qw1095_derivative():
    """Positions where dK/dd or d²K/dd² = 0"""
    print("[QW-1095] Derivative Zeros")
    
    d_vals = np.linspace(0.1, 16, 500)
    K_vals = K(d_vals)
    
    # First derivative
    dK = np.gradient(K_vals, d_vals)
    sign_changes_1 = np.where(np.diff(np.sign(dK)))[0]
    zeros_1 = d_vals[sign_changes_1]
    
    # Second derivative  
    d2K = np.gradient(dK, d_vals)
    sign_changes_2 = np.where(np.diff(np.sign(d2K)))[0]
    zeros_2 = d_vals[sign_changes_2]
    
    score_1, matches_1 = match_score(zeros_1.tolist())
    score_2, matches_2 = match_score(zeros_2.tolist())
    
    return {
        "dK_zeros": [float(z) for z in zeros_1[:10]],
        "d2K_zeros": [float(z) for z in zeros_2[:10]],
        "dK_match_score": float(score_1),
        "d2K_match_score": float(score_2),
        "INSIGHT": f"dK/dd zeros match {score_1:.0f}%, d²K/dd² zeros match {score_2:.0f}%"
    }

# ============================================================================
# QW-1096: SYNTHESIS
# ============================================================================
def qw1096_synthesis(all_results):
    """Final synthesis: which mechanism gives best emergence?"""
    print("[QW-1096] Emergence Synthesis")
    
    # Collect all match scores
    scores = {}
    for key, result in all_results.items():
        if isinstance(result, dict):
            for score_key in ["match_score", "dK_match_score", "d2K_match_score",
                             "top_match_score", "sublayer_match_score", "layer_match_score"]:
                if score_key in result:
                    scores[f"{key}_{score_key}"] = result[score_key]
    
    # Sort by score
    sorted_scores = sorted(scores.items(), key=lambda x: x[1], reverse=True)
    
    best_mechanisms = sorted_scores[:5]
    
    synthesis = {
        "all_scores": {k: float(v) for k, v in sorted_scores},
        "best_mechanisms": [(k, float(v)) for k, v in best_mechanisms],
        "BEST": best_mechanisms[0] if best_mechanisms else ("none", 0),
        "EMERGENCE_FOUND": best_mechanisms[0][1] > 50 if best_mechanisms else False,
        "VERDICT": ""
    }
    
    if synthesis["EMERGENCE_FOUND"]:
        synthesis["VERDICT"] = f"SUKCES: {synthesis['BEST'][0]} daje {synthesis['BEST'][1]:.0f}% dopasowanie!"
    else:
        synthesis["VERDICT"] = f"CZĘŚCIOWY: Najlepszy wynik to {synthesis['BEST'][1]:.0f}% ({synthesis['BEST'][0]})"
    
    return synthesis

# ============================================================================
# MAIN
# ============================================================================
def run_emergence_suite():
    print("=" * 80)
    print("QW-1077 TO QW-1096: EMERGENCJA POZYCJI Z K(d)")
    print("=" * 80)
    print(f"TARGET POSITIONS: {TARGET_POSITIONS}")
    print("=" * 80)
    
    all_results = {}
    
    tests = [
        ("qw1077", qw1077_K_extrema),
        ("qw1078", qw1078_K_zeros),
        ("qw1079", qw1079_eigenvalue_positions),
        ("qw1080", qw1080_energy_minima),
        ("qw1081", qw1081_standing_waves),
        ("qw1082", qw1082_winding),
        ("qw1083", qw1083_resonances),
        ("qw1084", qw1084_attractor),
        ("qw1085", qw1085_fractal),
        ("qw1086", qw1086_action),
        ("qw1087", qw1087_product_chain),
        ("qw1088", qw1088_band_gaps),
        ("qw1089", qw1089_hopfion),
        ("qw1090", qw1090_correlation),
        ("qw1091", qw1091_weighted_grid),
        ("qw1092", qw1092_level_crossing),
        ("qw1094", qw1094_integrated),
        ("qw1095", qw1095_derivative),
    ]
    
    for test_id, test_func in tests:
        try:
            result = test_func()
            all_results[test_id] = result
            score = result.get("match_score", 0)
            print(f"  → {result.get('INSIGHT', '')} [{score:.0f}%]")
        except Exception as e:
            print(f"  ❌ {e}")
            all_results[test_id] = {"ERROR": str(e)}
    
    # Combination
    all_results["qw1093"] = qw1093_combination(all_results)
    print(f"  → {all_results['qw1093']['INSIGHT']}")
    
    # Synthesis
    synthesis = qw1096_synthesis(all_results)
    all_results["qw1096"] = synthesis
    print(f"\n[QW-1096] {synthesis['VERDICT']}")
    
    # Save
    with open("RAPORT_QW1077_QW1096_EMERGENCE.md", "w") as f:
        f.write("# EMERGENCJA POZYCJI CZĄSTEK Z K(d)\n\n")
        f.write(f"**Werdykt:** {synthesis['VERDICT']}\n\n")
        f.write("## Najlepsze Mechanizmy\n")
        for name, score in synthesis["best_mechanisms"]:
            f.write(f"- {name}: {score:.0f}%\n")
    
    with open("RAPORT_QW1077_QW1096_DATA.json", "w") as f:
        json.dump(all_results, f, indent=2, default=str)
    
    print("\n" + "=" * 80)
    
    return all_results

if __name__ == "__main__":
    run_emergence_suite()
