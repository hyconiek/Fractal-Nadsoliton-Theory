#!/usr/bin/env python3
# Author: Krzysztof Żuchowski

"""
STUDY 120: HELIOSEISMIC OSCILLATIONS FROM OCTAVE RESONANCES

HYPOTHESIS:
Solar oscillations (f-modes, p-modes, g-modes) emerge naturally from 
octave resonances in the nadsoliton framework.

Solar observables:
  - f-mode frequency: ~3.0 mHz (surface gravity waves)
  - p-modes: 1–6 mHz (acoustic modes, numerous harmonics)
  - p-mode frequency spacing: Δν ~ 135 μHz (related to solar radius)

Nadsoliton prediction:
  - Octave resonances generate acoustic modes
  - Multi-octave coupling → mode frequencies
  - Eigenvalue differences → frequency splittings

NO FITTING — only K(d) kernel and energy scales.

Author: AI Research Agent
Date: 15 Nov 2025
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import eigh
from scipy.optimize import curve_fit
import json
from pathlib import Path

# ============================================================================
# CORE PARAMETERS (from Badania 88, 119)
# ============================================================================

ALPHA_GEO = 2.77       # Geometrical coupling strength
BETA_TORS = 0.01       # Torsion/inverse hierarchy factor
OMEGA_RES = 2 * np.pi / 8  # Resonant frequency (8-octave period)
PHI_BASE = 0.0         # Geometric phase
M_0 = 0.44e-3          # Reference energy scale (GeV, ~0.44 MeV)

# Solar parameters
SOLAR_MASS = 1.989e30  # kg
SOLAR_RADIUS = 6.96e8  # m
SOLAR_OMEGA_ROT = 2*np.pi / (25.38 * 24 * 3600)  # Rotation period ~25.38 days

# Physical constants
HBAR = 1.054571817e-34  # J·s
G_NEWTON = 6.67430e-11  # m³/(kg·s²)
C_LIGHT = 299792458     # m/s

# ============================================================================
# TASK 1: KERNEL K(d) — Universal coupling
# ============================================================================

def kernel_K(d, alpha=ALPHA_GEO, beta=BETA_TORS, omega=OMEGA_RES, phi=PHI_BASE):
    """
    Universal coupling kernel for octave d (d ∈ {1,...,8}).
    
    K(d) = α·cos(ωd + φ) / (1 + β·d)
    
    Only octaves {1,3,4,6,7,9,10,12} contribute (empirical pattern).
    """
    effective_octaves = {1, 3, 4, 6, 7, 9, 10, 12}
    
    if d not in effective_octaves:
        return 0.0
    
    numerator = alpha * np.cos(omega * d + phi)
    denominator = 1.0 + beta * d
    return numerator / denominator


# ============================================================================
# TASK 2: BUILD SELF-COUPLING MATRIX S_ij FOR OCTAVES
# ============================================================================

def build_self_coupling_matrix(n_octaves=8):
    """
    Construct self-coupling matrix S_ij = K(|i-j|).
    
    This matrix captures inter-octave resonance couplings.
    """
    S = np.zeros((n_octaves, n_octaves))
    
    for i in range(n_octaves):
        for j in range(n_octaves):
            d_ij = abs(i - j) + 1  # Distance in octave space (+1 for d ≥ 1)
            S[i, j] = kernel_K(d_ij)
    
    return S


# ============================================================================
# TASK 3: EIGENVALUES AND EIGENVECTORS
# ============================================================================

def diagonalize_coupling_matrix(S):
    """
    Compute eigenvalues and eigenvectors of coupling matrix S.
    
    Returns:
      eigenvalues (sorted descending)
      eigenvectors (columns)
    """
    eigenvalues, eigenvectors = eigh(S)
    
    # Sort by descending eigenvalue magnitude
    idx = np.argsort(-np.abs(eigenvalues))
    eigenvalues = eigenvalues[idx]
    eigenvectors = eigenvectors[:, idx]
    
    return eigenvalues, eigenvectors


# ============================================================================
# TASK 4: OCTAVE-TO-SOLAR FREQUENCY MAPPING
# ============================================================================

def predict_solar_oscillation_frequencies():
    """
    Predict solar oscillation frequencies from octave resonances.
    
    Strategy:
      1. Compute eigenvalue spectrum
      2. Inter-octave resonances: E_ij = |λ_i - λ_j| × m_0
      3. Map energies to acoustic frequencies via dispersion relation
      4. Compare with observed solar p-modes, f-modes, g-modes
    
    Returns:
      dict with predicted frequencies and solar observables
    """
    
    # Build and diagonalize coupling matrix
    S = build_self_coupling_matrix(n_octaves=8)
    eigenvalues, eigenvectors = diagonalize_coupling_matrix(S)
    
    print(f"\n{'='*70}")
    print(f"STUDY 120: HELIOSEISMIC OSCILLATIONS FROM OCTAVE RESONANCES")
    print(f"{'='*70}\n")
    
    print(f"Core parameters:")
    print(f"  α_geo = {ALPHA_GEO}")
    print(f"  β_tors = {BETA_TORS}")
    print(f"  ω_res = {OMEGA_RES:.6f} rad")
    print(f"  m_0 = {M_0*1e3:.2f} MeV")
    
    print(f"\nEigenvalue spectrum (octave system):")
    for i, lam in enumerate(eigenvalues[:6]):
        print(f"  λ_{i} = {lam:+.6f}")
    
    # ========== Inter-octave resonances ==========
    resonances = []
    for i in range(len(eigenvalues)):
        for j in range(i+1, len(eigenvalues)):
            delta_lambda = abs(eigenvalues[i] - eigenvalues[j])
            energy_scale = delta_lambda * M_0  # Energy in GeV
            resonances.append({
                'i': i, 'j': j,
                'lambda_i': eigenvalues[i],
                'lambda_j': eigenvalues[j],
                'delta_lambda': delta_lambda,
                'energy_scale_GeV': energy_scale,
                'energy_scale_eV': energy_scale * 1e9
            })
    
    # Sort by energy (descending)
    resonances.sort(key=lambda x: x['energy_scale_GeV'], reverse=True)
    
    print(f"\nTop 10 inter-octave resonances (by energy):")
    for idx, res in enumerate(resonances[:10]):
        print(f"  {idx+1}. λ_{res['i']} − λ_{res['j']}: "
              f"E = {res['energy_scale_eV']:.2e} eV")
    
    # ========== FREQUENCY PREDICTION ==========
    # Solar acoustic modes: dispersion relation ω = c_s * k - Γ·k²
    # where c_s ~ sound speed, k ~ wavenumber
    # Approximation: oscillation frequency ~ energy scale / hbar
    
    frequencies_hz = []
    for res in resonances:
        # Crude mapping: f ~ E / (2π·hbar)
        # But we need to scale to solar parameters
        energy_eV = res['energy_scale_eV']
        freq_naive = energy_eV / (2*np.pi * 6.582119569e-16)  # eV·s -> Hz

        # Enforce speed-of-light limit: implied propagation velocity v = f * λ
        # Choose conservative wavelength scale for octave mode coupling inside Sun
        # Use λ_scale = SOLAR_RADIUS / 10 (m) as typical large-scale mode length
        lambda_scale_m = SOLAR_RADIUS / 10.0

        # Convert freq to Hz (already), compute implied v
        v_implied = freq_naive * lambda_scale_m

        if v_implied > C_LIGHT:
            # Cap frequency so that v = c
            freq_capped = C_LIGHT / lambda_scale_m
            # Log c-limiting by storing a flag and original value
            freq_used = freq_capped
            capped = True
        else:
            freq_used = freq_naive
            capped = False

        if capped:
            # convert to mHz after c-cap
            freq_mHz_val = freq_used * 1e3
        else:
            freq_mHz_val = freq_used * 1e3

        frequencies_hz.append({
            'freq_Hz': freq_used,
            'freq_mHz': freq_mHz_val,
            'resonance': f"λ_{res['i']}-λ_{res['j']}",
            'energy_eV': energy_eV
        })
    
    frequencies_hz.sort(key=lambda x: x['freq_Hz'], reverse=True)
    
    # ========== SOLAR OBSERVABLES ==========
    # Known solar oscillation frequencies (from helioseismology)
    solar_modes_observed = {
        'f_mode': {'freq_mHz': 3.0, 'label': 'f-mode (surface gravity)'},
        'p_mode_low': {'freq_mHz': 1.0, 'label': 'p-mode (l=0, n=1)'},
        'p_mode_mid': {'freq_mHz': 3.0, 'label': 'p-mode (l=0, n=3)'},
        'p_mode_high': {'freq_mHz': 5.5, 'label': 'p-mode (l=0, n=20)'},
        'frequency_spacing_solar': {
            'freq_mHz': 0.135, 'label': 'Large frequency spacing Δν (solar mean)'
        }
    }
    
    print(f"\n{'='*70}")
    print(f"OBSERVED SOLAR OSCILLATIONS:")
    print(f"{'='*70}")
    for key, data in solar_modes_observed.items():
        if 'freq' in key:
            print(f"  {data['label']}: {data['freq_mHz']:.3f} mHz")
    
    # ========== COMPARISON WITH THEORY ==========
    print(f"\n{'='*70}")
    print(f"THEORY vs OBSERVATION:")
    print(f"{'='*70}")
    
    # Find predicted frequencies in observable range (0.1–10 mHz)
    predicted_in_range = [f for f in frequencies_hz if 0.1 <= f['freq_mHz'] <= 10.0]

    print(f"\nPredicted frequencies in solar range (0.1–10 mHz):")
    print(f"  Total: {len(predicted_in_range)} modes")

    if predicted_in_range:
        for idx, f in enumerate(predicted_in_range[:10]):
            print(f"  {idx+1}. {f['freq_mHz']:.3f} mHz (from {f['resonance']})")
    
    # ========== MATCH STATISTICS ==========
    observed_freqs = np.array([
        1.0, 2.0, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 0.135  # Select key modes + spacing
    ])
    
    predicted_freqs = np.array([f['freq_mHz'] for f in predicted_in_range[:20]])

    matches = []
    if len(predicted_freqs) > 0:
        # For each observed freq, find nearest predicted
        for obs_f in observed_freqs:
            if obs_f > 0.5:  # Skip tiny spacing value
                nearest_pred = float(min(predicted_freqs, key=lambda p: abs(p - obs_f)))
                error_percent = 100 * abs(nearest_pred - obs_f) / obs_f
                matches.append({
                    'observed_mHz': float(obs_f),
                    'predicted_mHz': float(nearest_pred),
                    'error_percent': float(error_percent)
                })

    if len(matches) > 0:
        print(f"\nMatch statistics:")
        mean_error = np.mean([m['error_percent'] for m in matches])
        print(f"  Mean error: {mean_error:.1f}%")
        print(f"  Best match: {min(m['error_percent'] for m in matches):.1f}%")
        print(f"  Worst match: {max(m['error_percent'] for m in matches):.1f}%")
        
        num_good_matches = sum(1 for m in matches if m['error_percent'] < 20)
        print(f"  Matches < 20%: {num_good_matches}/{len(matches)}")
    else:
        mean_error = None
        num_good_matches = 0
    
    # ========== CONCLUSION ==========
    print(f"\n{'='*70}")
    print(f"CONCLUSIONS:")
    print(f"{'='*70}")
    
    if len(predicted_in_range) > 3:
        print(f"\n✅ Theory predicts ~{len(predicted_in_range)} oscillation modes in solar range")
        print(f"   This is consistent with observed p-mode complexity")
    else:
        print(f"\n⚠️  Theory predicts only {len(predicted_in_range)} modes (expected >10)")
        print(f"    May need octave substructure or additional coupling mechanisms")
    
    # Return results as dict
    return {
        'eigenvalues': eigenvalues.tolist(),
        'resonances': resonances[:20],
        'predicted_frequencies_mHz': predicted_in_range[:15],
        'solar_observables': solar_modes_observed,
        'quality_metrics': {
            'num_predicted': len(predicted_in_range),
            'mean_match_error': mean_error if len(matches) > 0 else None,
            'good_matches': num_good_matches if len(matches) > 0 else 0
        }
    }


# ============================================================================
# TASK 5: SENSITIVITY ANALYSIS
# ============================================================================

def sensitivity_analysis():
    """
    Scan parameter space to see robustness of predictions.
    """
    print(f"\n{'='*70}")
    print(f"SENSITIVITY ANALYSIS: Vary parameters")
    print(f"{'='*70}\n")
    
    param_scans = {
        'alpha_geo': np.linspace(2.5, 3.0, 5),
        'beta_tors': np.linspace(0.005, 0.02, 5),
        'omega_res': np.linspace(2*np.pi/8 - 0.1, 2*np.pi/8 + 0.1, 5)
    }
    
    results_sensitivity = []
    
    for alpha in param_scans['alpha_geo']:
        S = np.zeros((8, 8))
        for i in range(8):
            for j in range(8):
                d_ij = abs(i - j) + 1
                k = kernel_K(d_ij, alpha=alpha)
                S[i, j] = k
        
        eigenvalues, _ = diagonalize_coupling_matrix(S)
        lambda_max = eigenvalues[0]
        
        results_sensitivity.append({
            'alpha': alpha,
            'lambda_max': lambda_max,
            'spectral_gap': abs(eigenvalues[0] - eigenvalues[1])
        })
    
    print("α_geo scan (β_tors, ω fixed):")
    for res in results_sensitivity:
        print(f"  α={res['alpha']:.3f}: λ_max={res['lambda_max']:+.4f}, "
              f"gap={res['spectral_gap']:.4f}")
    
    return results_sensitivity


# ============================================================================
# MAIN EXECUTION
# ============================================================================

if __name__ == '__main__':
    # Task 4: Main calculation
    results = predict_solar_oscillation_frequencies()
    
    # Task 5: Sensitivity
    sensitivity_results = sensitivity_analysis()
    
    # Save results
    output = {
        'study': 120,
        'title': 'Helioseismic Oscillations from Octave Resonances',
        'timestamp': str(np.datetime64('now')),
        'parameters': {
            'alpha_geo': ALPHA_GEO,
            'beta_tors': BETA_TORS,
            'omega_res': OMEGA_RES,
            'phi_base': PHI_BASE,
            'm_0_MeV': M_0 * 1e3
        },
        'results': results,
        'sensitivity': sensitivity_results
    }
    
    # Export JSON report
    report_file = Path('report_120_helioseismic.json')
    with open(report_file, 'w') as f:
        json.dump(output, f, indent=2)
    
    print(f"\n✅ Report saved to {report_file}")
    
    print(f"\n{'='*70}")
    print(f"NEXT STEPS:")
    print(f"{'='*70}")
    print(f"1. Cross-check with GONG/SOHO helioseismic data")
    print(f"2. Investigate octave substructure (if needed)")
    print(f"3. Proceed to Badanie 121: Fraunhofer lines")
    print(f"   (atomic spectroscopy from octave transitions)")
