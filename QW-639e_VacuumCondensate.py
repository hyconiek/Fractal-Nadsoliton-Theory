#!/usr/bin/env python3
"""
QW-639e: ELECTRON MASS - HYBRID OCTAVE-LAYER MODEL (CORRECTED)
Purpose: Proper octave-layer separation with correct β scaling
Key: β^N means DOWN from Planck (damping), not up
Date: 2025-12-06
"""

import numpy as np
from scipy.linalg import eigh
from scipy.stats import entropy

print("="*80)
print("QW-639e: HYBRID OCTAVE-LAYER MODEL (CORRECTED)")
print("="*80)
print("Insight: Warstwy określają SPADEK od Plancka, nie wzrost")
print("="*80)

# Constants
ALPHA_GEO = 4 * np.log(2)
OMEGA = np.pi / 4
PHI = np.pi / 6
BETA_TORS = 0.01
M_PLANCK_GeV = 1.2209e19

# Experimental
M_ELECTRON_MeV = 0.511
M_MUON_MeV = 105.66
M_TAU_MeV = 1776.86

kappa = ALPHA_GEO / (OMEGA * PHI)
print(f"\nκ = {kappa:.4f}")
print(f"β = {BETA_TORS}")

# ============================================================================
# KEY REALIZATION
# ============================================================================

print("\n" + "="*80)
print("CORRECTED UNDERSTANDING:")
print("="*80)
print("• Wszystkie LEPTONY są w WARST WIE 10 (cząstki)")
print("• Różnią się OKTAWĄ (1, 4, 7)")
print("• Warstwa 10: β^10 = 10^-20 (spadek od Plancka)")
print("• Oktawa N: Amplifikacja poprzez NIELINIOWY rezonans")

# Vacuum Hamiltonian
N_octaves = 12
def K(d):
    return ALPHA_GEO * np.cos(OMEGA * d + PHI) / (1 + BETA_TORS * d)

H_vac = np.zeros((N_octaves, N_octaves))
for i in range(N_octaves):
    for j in range(N_octaves):
        H_vac[i, j] = -K(abs(i - j))

evals, evecs = eigh(H_vac)

# Base parameters
C_stability = 12.027675
N_layer = 10  # All leptons at same layer
frac_damp = BETA_TORS ** N_layer

# ============================================================================
# HYPOTHESIS: OKTAWA określa KONDENSACJĘ PRÓŻNI
# ============================================================================

print("\n" + "="*80)
print("NOWA HIPOTEZA: Vacuum Condensate Density")
print("="*80)
print("⟨H⟩(N) = Gęstość kondensatu próżni w oktawie N")
print("Im wyższa oktawa → większa gęstość → większa masa")

# Higgs VEV scaling (to be determined)
def vacuum_VEV(N_octave):
    """
    Vacuum expectation value grows with octave
    Hypothesis: ⟨H⟩ ∝ exp(α × N)
    """
    # Solve for α that gives muon/electron ratio
    # m_μ/m_e = ⟨H⟩(4)/⟨H⟩(1) × (κ^4/κ^1) × (A_res ratio)
    # 207 ≈ exp(α × 3) × smaller terms
    
    # Rough estimate: exp(α×3) ~ 100 → α ~ ln(100)/3 ~ 1.5
    alpha_condensate = 1.5
    return np.exp(alpha_condensate * N_octave)

print(f"\nVacuum condensate scaling:")
print(f"  ⟨H⟩(1) = {vacuum_VEV(1):.2f} (electron)")
print(f"  ⟨H⟩(4) = {vacuum_VEV(4):.2f} (muon)")
print(f"  ⟨H⟩(7) = {vacuum_VEV(7):.2f} (tau)")
print(f"  Ratio μ/e: {vacuum_VEV(4)/vacuum_VEV(1):.1f}")
print(f"  Ratio τ/e: {vacuum_VEV(7)/vacuum_VEV(1):.1f}")

# ============================================================================
# PREDICTIONS
# ============================================================================

def compute_mass_hybrid(N_octave):
    """Unified formula with vacuum condensate"""
    # Octave amplification
    oct_amp = kappa ** (N_octave / 12)
    
    # Resonance
    psi = np.zeros(N_octaves)
    psi[N_octave] = 1.0
    A_res = abs(np.dot(psi, evecs[:, 0]))
    
    # Processing intensity (base)
    psi_dist = np.exp(-0.5 * ((np.arange(N_octaves) - N_octave) / 0.8)**2)
    psi_dist /= np.sum(psi_dist)
    S = entropy(psi_dist, base=2)
    I_proc = (S * 0.1 / C_stability)
    
    # VACUUM CONDENSATE (key mechanism!)
    VEV = vacuum_VEV(N_octave)
    
    # Total mass
    m_GeV = M_PLANCK_GeV * 1 * oct_amp * A_res * frac_damp * I_proc * VEV
    return m_GeV * 1000, {'oct_amp': oct_amp, 'A_res': A_res,<br/>                           'I_proc': I_proc, 'VEV': VEV}

print("\n" + "="*80)
print("PREDICTIONS WITH VACUUM CONDENSATE")
print("="*80)

leptons = [
    ('Electron', 1, M_ELECTRON_MeV),
    ('Muon', 4, M_MUON_MeV),
    ('Tau', 7, M_TAU_MeV)
]

for name, N_oct, m_exp in leptons:
    m_pred, components = compute_mass_hybrid(N_oct)
    error = abs(m_pred - m_exp) / m_exp * 100
    
    print(f"\n{name} (Octave {N_oct}):")
    print(f"  κ^({N_oct}/12) = {components['oct_amp']:.4f}")
    print(f"  A_res = {components['A_res']:.4f}")
    print(f"  ⟨H⟩({N_oct}) = {components['VEV']:.2f}")
    print(f"  Pred: {m_pred:.2f} MeV")
    print(f"  Exp:  {m_exp:.2f} MeV")
    print(f"  Error: {error:.1f}%")

print("\n" + "="*80)
print("FINALNA FORMUŁA:")
print("="*80)
print("\n┌────────────────────────────────────────────────────────┐")
print("│ m = m_Planck × |W| × κ^(N/12) × A_res × β^10 ×       │")
print("│     I_proc × ⟨H⟩(N)                                   │")
print("└────────────────────────────────────────────────────────┘")
print("\nwhere:")
print("  ⟨H⟩(N) = exp(α × N) - Vacuum condensate density")
print("  α ≈ 1.5 - Condensation rate")

print("\n🔑 KLUCZOWA IDEA:")
print("  Wszystkie leptony: Warstwa 10 (cząstki)")
print("  Różne oktawy → różna gęstość próżni → różne masy!")
print("="*80)
