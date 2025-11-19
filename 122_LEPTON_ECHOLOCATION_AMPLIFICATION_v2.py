#!/usr/bin/env python3
# Author: Krzysztof Żuchowski

"""
STUDY 122 REVISED: LEPTON MASS HIERARCHY — ENHANCED ECHOLOCATION

The naive echolocation attempt gave 99% error because the echo amplification
was too weak. This version uses STRONGER coupling mechanisms:

1. TOPOLOGICAL WINDING PHASE: phase ~ W(d_i, d_j) instead of just ω·Δd
2. RESONANT TUNNELING: coupling ~ exp(-barriers between octaves)
3. CHIRAL FEEDBACK: multiplicative amplification across octave chains

Expected outcome:
  If mechanism works → m_μ/m_e ≈ 207, m_τ/m_μ ≈ 16.8
"""

import numpy as np
import json
from pathlib import Path
from scipy.special import gamma, erf

# ============================================================================
# CORE PARAMETERS
# ============================================================================

ALPHA_GEO = 2.77
BETA_TORS = 0.01
OMEGA_RES = 2 * np.pi / 8
PHI_BASE = 0.0
M_0 = 0.44e-3  # GeV

# Observed ratios
M_E_OBS = 0.5109989e-3
M_MU_OBS = 105.6583745e-3
M_TAU_OBS = 1776.86e-3
RATIO_MU_E_OBS = M_MU_OBS / M_E_OBS
RATIO_TAU_MU_OBS = M_TAU_OBS / M_MU_OBS

def kernel_K(d):
    effective_octaves = {1, 3, 4, 6, 7, 9, 10, 12}
    if d not in effective_octaves:
        return 0.0
    numerator = ALPHA_GEO * np.cos(OMEGA_RES * d + PHI_BASE)
    denominator = 1.0 + BETA_TORS * d
    return numerator / denominator


# ============================================================================
# REVISED AMPLIFICATION: STRONGER MECHANISMS
# ============================================================================

def compute_topological_winding_number(d1, d2):
    """
    Compute topological winding number W(d1, d2).
    
    Interpretation: how many "times" the field winds around between octaves
    d1 and d2. Related to Berry phase and adiabatic parameter.
    
    Formula (approximate): W ~ Σ_d ω(d) = ω * |d2-d1| * (structure factor)
    
    We use a proxy: W ~ (d2-d1) * α_geo (coupling strength)
    """
    distance = abs(d2 - d1)
    # Winding number ~ distance × coupling strength
    W = distance * ALPHA_GEO
    return W


def resonant_tunneling_coupling(d_lepton, d_echo):
    """
    Quantum tunneling through octave potential barriers.
    
    Tunneling probability: T ~ exp(-2κ·L)
    where κ ~ inverse penetration depth
    
    For octaves: effective barrier ~ 1/|K(d)|
    
    Result: coupling via tunneling ~ T ~ exp(-β·|K(d)|^{-1})
    This gives EXPONENTIAL enhancement for small K(d) values!
    """
    
    k_min = 0.1  # Minimum coupling strength
    barriers = []
    
    for d in range(min(d_lepton, d_echo), max(d_lepton, d_echo) + 1):
        k = kernel_K(d)
        if abs(k) < k_min:
            barrier = 1.0 / k_min  # Large barrier
        else:
            barrier = 1.0 / abs(k)
        barriers.append(barrier)
    
    total_barrier = sum(barriers)
    # Tunneling decay: T ~ exp(-total_barrier / normalization)
    T = np.exp(-total_barrier / 10.0)  # 10 = normalization
    
    return max(T, 0.01)  # Floor at 1%


def compute_amplified_mass(d_lepton):
    """
    Enhanced echolocation with THREE amplification mechanisms:
    1. Topological winding number
    2. Resonant tunneling
    3. Chiral feedback (multiplicative cascade)
    """
    
    # Direct coupling (baseline)
    k_direct = kernel_K(d_lepton)
    m_baseline = abs(k_direct) * M_0
    
    # MECHANISM 1: Topological phases
    echo_from_e = compute_topological_winding_number(d_lepton, 1)  # d=1 is electron
    phase_1 = echo_from_e / 10.0  # Normalized phase
    topological_factor = 1.0 + np.sin(phase_1) * 2.0  # sin can give -1 to +1
    
    # MECHANISM 2: Resonant tunneling
    tunnel_prob = resonant_tunneling_coupling(d_lepton, 1)
    
    # MECHANISM 3: Chiral feedback cascade
    # Path: lepton (d) ← d-3 ← d-6 ← ... ← d=1 (electron)
    chiral_cascade = 1.0
    d_current = d_lepton
    while d_current > 1:
        d_next = max(1, d_current - 3)
        k_step = kernel_K(d_current)
        k_target = kernel_K(d_next)
        
        # Chiral coupling ~ (K_target / K_source)
        if abs(k_step) > 0.01:
            ratio = k_target / k_step
            chiral_cascade *= (1.0 + abs(ratio))
        
        d_current = d_next
    
    # TOTAL AMPLIFICATION
    total_amplification = topological_factor * (1.0 + tunnel_prob * 10.0) * chiral_cascade
    
    m_final = m_baseline * total_amplification
    
    return {
        'mass_GeV': m_final,
        'mass_MeV': m_final * 1e3,
        'amplification': total_amplification,
        'baseline_MeV': m_baseline * 1e3,
        'components': {
            'topological_factor': topological_factor,
            'tunneling_factor': 1.0 + tunnel_prob * 10.0,
            'chiral_cascade': chiral_cascade
        }
    }


# ============================================================================
# MAIN CALCULATION
# ============================================================================

print(f"\n{'='*70}")
print(f"STUDY 122 REVISED: ENHANCED ECHOLOCATION AMPLIFICATION")
print(f"{'='*70}\n")

print(f"Observed lepton mass ratios:")
print(f"  m_μ/m_e = {RATIO_MU_E_OBS:.4f}")
print(f"  m_τ/m_μ = {RATIO_TAU_MU_OBS:.4f}\n")

results = {}
for lepton, d in [('e', 1), ('μ', 4), ('τ', 7)]:
    print(f"Lepton: {lepton} (octave d={d})")
    print(f"{'-'*60}")
    
    result = compute_amplified_mass(d)
    results[lepton] = result
    
    print(f"  Direct coupling: K({d}) = {kernel_K(d):.4f}")
    print(f"  Baseline mass: {result['baseline_MeV']:.4f} MeV")
    print(f"  Total amplification factor: {result['amplification']:.4f}x")
    print(f"  AMPLIFIED MASS: {result['mass_MeV']:.4f} MeV")
    print(f"  Components:")
    print(f"    - Topological: {result['components']['topological_factor']:.4f}x")
    print(f"    - Tunneling: {result['components']['tunneling_factor']:.4f}x")
    print(f"    - Chiral cascade: {result['components']['chiral_cascade']:.4f}x")
    print()

# Calculate mass ratios
m_e = results['e']['mass_GeV']
m_mu = results['μ']['mass_GeV']
m_tau = results['τ']['mass_GeV']

ratio_mu_e_theory = m_mu / m_e
ratio_tau_mu_theory = m_tau / m_mu

error_mu_e = 100 * abs(ratio_mu_e_theory - RATIO_MU_E_OBS) / RATIO_MU_E_OBS
error_tau_mu = 100 * abs(ratio_tau_mu_theory - RATIO_TAU_MU_OBS) / RATIO_TAU_MU_OBS

print(f"\n{'='*70}")
print(f"AMPLIFIED PREDICTIONS vs OBSERVATION:")
print(f"{'='*70}\n")

print(f"{'Ratio':<15} {'Theory':<15} {'Observed':<15} {'Error %':<12}")
print(f"{'-'*60}")
print(f"{'m_μ/m_e':<15} {ratio_mu_e_theory:<15.4f} {RATIO_MU_E_OBS:<15.4f} {error_mu_e:<12.2f}%")
print(f"{'m_τ/m_μ':<15} {ratio_tau_mu_theory:<15.4f} {RATIO_TAU_MU_OBS:<15.4f} {error_tau_mu:<12.2f}%")

mean_error = (error_mu_e + error_tau_mu) / 2

print(f"\n{'='*70}")
print(f"RESULT:")
print(f"{'='*70}\n")

if error_mu_e < 10 and error_tau_mu < 10:
    print(f"✅ STRONG SUCCESS")
    print(f"   Both ratios within 10% error!")
elif error_mu_e < 20 and error_tau_mu < 20:
    print(f"🟡 MODERATE SUCCESS")
    print(f"   Getting closer — amplification mechanisms working")
elif mean_error < 50:
    print(f"⚠️  PARTIAL SUCCESS")
    print(f"   Amplification direction is correct but magnitude off")
else:
    print(f"❌ NEEDS MORE WORK")
    print(f"   Mean error: {mean_error:.1f}%")

# Save JSON
output = {
    'study': 122,
    'title': 'Lepton Mass Hierarchy (Revised with Enhanced Echolocation)',
    'parameters': {
        'alpha_geo': ALPHA_GEO,
        'beta_tors': BETA_TORS,
        'omega_res': OMEGA_RES,
        'm_0_MeV': M_0 * 1e3
    },
    'predictions': {
        'electron_MeV': results['e']['mass_MeV'],
        'muon_MeV': results['μ']['mass_MeV'],
        'tau_MeV': results['τ']['mass_MeV'],
        'ratio_mu_e': ratio_mu_e_theory,
        'ratio_tau_mu': ratio_tau_mu_theory,
        'error_mu_e_percent': error_mu_e,
        'error_tau_mu_percent': error_tau_mu,
        'mean_error_percent': mean_error
    },
    'observed': {
        'electron_MeV': M_E_OBS * 1e3,
        'muon_MeV': M_MU_OBS * 1e3,
        'tau_MeV': M_TAU_OBS * 1e3,
        'ratio_mu_e': RATIO_MU_E_OBS,
        'ratio_tau_mu': RATIO_TAU_MU_OBS
    }
}

with open(Path('report_122_enhanced_echolocation.json'), 'w') as f:
    json.dump(output, f, indent=2)

print(f"\n✅ Report saved to report_122_enhanced_echolocation.json")
