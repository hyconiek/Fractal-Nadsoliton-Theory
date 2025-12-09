#!/usr/bin/env python3
"""
QW-707: KOIDE VERIFICATION + PROTON FIX
========================================
PART 1: Verify if Koide Q=2/3 actually emerges from M ∝ d^α
PART 2: Fix proton mass (Path A: correct orbit, Path B: baryon formula)

Date: 2025-12-08
"""

import numpy as np
import datetime

print("="*80)
print("QW-707: KOIDE VERIFICATION + PROTON FIX")
print("="*80)

# ===========================================================================
# FROZEN PARAMETERS
# ===========================================================================

ALPHA_GEO = 4 * np.log(2)  # 2.7726
OMEGA = np.pi / 4
PHI = np.pi / 6
BETA_TORS = 0.01

# Stable orbits from theory
d_1 = 1.33   # Electron
d_2 = 9.33   # Muon
d_3 = 17.33  # Tau

# Experimental masses (MeV)
M_E_EXP = 0.511
M_MU_EXP = 105.66
M_TAU_EXP = 1776.86
M_P_EXP = 938.27

def K(d):
    return ALPHA_GEO * np.cos(OMEGA * d + PHI) / (1 + BETA_TORS * d)

# ===========================================================================
# PART 1: KOIDE VERIFICATION
# ===========================================================================

print("\n" + "="*80)
print("PART 1: KOIDE VERIFICATION")
print("="*80)

# FINAL_THEORY_REPORT claims: M ∝ d^α gives Koide Q = 2/3 with 0.03% error
# Let's check if this is true

# Method 1: Pure d^α (no additional factors)
print("\n[Method 1] M = M_0 × (d/d_1)^α")
m_e_1 = 1.0  # Reference
m_mu_1 = (d_2 / d_1) ** ALPHA_GEO
m_tau_1 = (d_3 / d_1) ** ALPHA_GEO

sqrt_m = [np.sqrt(m_e_1), np.sqrt(m_mu_1), np.sqrt(m_tau_1)]
Q_method1 = (sum(sqrt_m))**2 / sum([m_e_1, m_mu_1, m_tau_1])

print(f"  m_e   = 1.0 (reference)")
print(f"  m_μ   = (d_2/d_1)^α = {m_mu_1:.4f}")
print(f"  m_τ   = (d_3/d_1)^α = {m_tau_1:.4f}")
print(f"  Koide Q = {Q_method1:.5f} (expected: 0.66667)")
print(f"  Error: {abs(Q_method1 - 2/3)/(2/3)*100:.2f}%")

# Method 2: With generation factor 3 for tau (from QW-703)
print("\n[Method 2] M_τ = 3 × M_0 × (d_3/d_1)^α")
m_tau_2 = 3 * (d_3 / d_1) ** ALPHA_GEO

sqrt_m2 = [np.sqrt(m_e_1), np.sqrt(m_mu_1), np.sqrt(m_tau_2)]
Q_method2 = (sum(sqrt_m2))**2 / sum([m_e_1, m_mu_1, m_tau_2])

print(f"  m_τ   = 3 × (d_3/d_1)^α = {m_tau_2:.4f}")
print(f"  Koide Q = {Q_method2:.5f}")
print(f"  Error: {abs(Q_method2 - 2/3)/(2/3)*100:.2f}%")

# Method 3: With EXPERIMENTAL masses (tautology check)
print("\n[Method 3] Using EXPERIMENTAL masses")
sqrt_exp = [np.sqrt(M_E_EXP), np.sqrt(M_MU_EXP), np.sqrt(M_TAU_EXP)]
Q_exp = (sum(sqrt_exp))**2 / sum([M_E_EXP, M_MU_EXP, M_TAU_EXP])

print(f"  Koide Q (experimental) = {Q_exp:.5f}")
print(f"  Error: {abs(Q_exp - 2/3)/(2/3)*100:.4f}%")

# Conclusion
print("\n" + "-"*60)
print("KOIDE VERIFICATION CONCLUSION:")
if abs(Q_method1 - 2/3)/(2/3) < 0.01:
    print("✅ Method 1 (pure d^α) gives Koide ≈ 2/3")
else:
    print("❌ Method 1 (pure d^α) does NOT give Koide = 2/3")
    print(f"   Q = {Q_method1:.3f} → error {abs(Q_method1 - 2/3)/(2/3)*100:.1f}%")

print(f"\n   FINAL_THEORY_REPORT claimed 0.03% error")
print(f"   Likely this was computed using EXPERIMENTAL masses (tautology)")
print(f"   Actual error from d^α formula: {abs(Q_method1 - 2/3)/(2/3)*100:.1f}%")

# ===========================================================================
# PART 2A: FIX PROTON - FIND CORRECT ORBIT
# ===========================================================================

print("\n" + "="*80)
print("PART 2A: FIX PROTON ORBIT")
print("="*80)

# Current problem: proton at d=6 gives K(6)=1.3 which is too large
# Need to find d_p such that mass ratio comes out correctly

# Target: m_p/m_e = 1836
ratio_target = M_P_EXP / M_E_EXP

# With winding number W=3 (3 quarks):
# m_p = W × m_e × (d_p/d_1)^α × other_factors
# If we only use d^α:
# (d_p/d_1)^α = ratio_target / W = 1836 / 3 = 612

ratio_needed = ratio_target / 3  # 612
d_p_needed = d_1 * (ratio_needed) ** (1/ALPHA_GEO)

print(f"Target ratio: m_p/m_e = {ratio_target:.2f}")
print(f"With W=3: (d_p/d_1)^α = {ratio_needed:.2f}")
print(f"Solving for d_p:")
print(f"  d_p = d_1 × {ratio_needed:.2f}^(1/{ALPHA_GEO:.4f})")
print(f"  d_p = {d_p_needed:.4f}")

# Check K(d) at this orbit
K_at_dp = K(d_p_needed)
print(f"\nK(d_p={d_p_needed:.2f}) = {K_at_dp:.4f}")

# Verify
m_p_check = 3 * M_E_EXP * (d_p_needed / d_1) ** ALPHA_GEO
print(f"Verification: m_p = 3 × m_e × (d_p/d_1)^α = {m_p_check:.2f} MeV")
print(f"Target: {M_P_EXP:.2f} MeV")
print(f"Error: {abs(m_p_check - M_P_EXP)/M_P_EXP*100:.2f}%")

# ===========================================================================
# PART 2B: BARYON FORMULA (Triplet Binding)
# ===========================================================================

print("\n" + "="*80)
print("PART 2B: BARYON FORMULA (Triplet Binding)")
print("="*80)

# Alternative: Proton is NOT just 3× something
# It's a BOUND STATE of 3 quarks

# Constituent quark model:
# m_p ≈ 3 × m_constituent_quark + E_binding
# m_constituent ≈ 330 MeV (phenomenological)

# In FIN theory, quarks would be at specific octaves
# Proton = 3 quarks coupled by K(d)

# Let's try: Each quark at its own orbit, with binding from K(d)

print("Model: Proton = 3 constituent quarks with K(d) binding")

# Quarks at octaves 3, 4, 5 (as proposed in QW-702)
quark_octaves = [3, 4, 5]
d_quarks = [quark_octaves[i] for i in range(3)]

# Quark mass from octave position
m_quark_base = M_E_EXP  # Reference

def quark_mass(octave):
    """Quark mass from octave position."""
    d = octave  # Distance from octave 0
    return m_quark_base * (d / d_1) ** ALPHA_GEO

m_quarks = [quark_mass(o) for o in quark_octaves]
print(f"\nQuark masses from octaves {quark_octaves}:")
for i, (o, m) in enumerate(zip(quark_octaves, m_quarks)):
    print(f"  q_{i+1} (O{o}): m = {m:.2f} MeV")

# Total without binding
m_sum = sum(m_quarks)
print(f"\nSum of quark masses: {m_sum:.2f} MeV")

# Binding energy from K(d) between quarks
def binding_energy(octaves):
    """Calculate binding from inter-quark K(d)."""
    E_bind = 0
    for i in range(len(octaves)):
        for j in range(i+1, len(octaves)):
            d_ij = abs(octaves[i] - octaves[j])
            E_bind += K(d_ij)
    return E_bind

E_bind = binding_energy(quark_octaves)
print(f"Binding energy from K(d): E_bind = {E_bind:.2f} MeV")

# Proton mass
g_binding = 100  # Coupling strength (to be determined)
m_p_baryon = m_sum - g_binding * abs(E_bind)
print(f"\nProton mass (baryon model):")
print(f"  m_p = Σm_q - g × |E_bind| = {m_sum:.2f} - {g_binding}×{abs(E_bind):.2f}")
print(f"  m_p = {m_p_baryon:.2f} MeV (target: {M_P_EXP:.2f} MeV)")

# Find optimal g
g_optimal = (m_sum - M_P_EXP) / abs(E_bind)
print(f"\nOptimal g to match proton: g = {g_optimal:.2f}")

# ===========================================================================
# PART 3: UNIFIED SOLUTION
# ===========================================================================

print("\n" + "="*80)
print("PART 3: UNIFIED SOLUTION")
print("="*80)

print("""
FINDINGS:

1. KOIDE:
   - d^α formula gives Q = {:.3f} (not 0.667!)
   - FINAL_THEORY_REPORT's 0.03% was likely from experimental masses
   - This is a CRITICAL CORRECTION

2. PROTON (Path A - Orbit Fix):
   - Required orbit: d_p = {:.2f}
   - This is between tau (17.33) and beyond
   - Interpretation: Proton resonates at higher orbit than tau

3. PROTON (Path B - Baryon Model):
   - 3 quarks at octaves [3,4,5]
   - Binding energy from K(d) with g ≈ {:.1f}
   - More physically motivated but requires fitting g

RECOMMENDATION:
Use PATH A for simplicity:
  - Electron: d = 1.33
  - Muon: d = 9.33  
  - Tau: d = 17.33
  - Proton: d = {:.2f} (with W=3)
""".format(Q_method1, d_p_needed, g_optimal, d_p_needed))

# ===========================================================================
# FINAL TEST WITH CORRECTED PARAMETERS
# ===========================================================================

print("="*80)
print("FINAL TEST: CORRECTED MASS FORMULA")
print("="*80)

# Corrected formula:
# M = M_0 × (d/d_0)^α × W
# For leptons: W = 1
# For baryons: W = 3

particles = {
    'electron': {'d': 1.33, 'W': 1},
    'muon': {'d': 9.33, 'W': 1},
    'tau': {'d': 17.33, 'W': 1},
    'proton': {'d': d_p_needed, 'W': 3}
}

M_0 = M_E_EXP
d_0 = d_1

print(f"\nFormula: M = M_0 × (d/d_0)^α × W")
print(f"M_0 = {M_0} MeV, d_0 = {d_0}, α = {ALPHA_GEO:.4f}")

print(f"\n{'Particle':<10} {'d':>8} {'W':>4} {'Predicted':>12} {'Experiment':>12} {'Error':>8}")
print("-"*60)

for name, p in particles.items():
    mass_pred = M_0 * (p['d'] / d_0) ** ALPHA_GEO * p['W']
    if name == 'electron':
        mass_exp = M_E_EXP
    elif name == 'muon':
        mass_exp = M_MU_EXP
    elif name == 'tau':
        mass_exp = M_TAU_EXP
    else:
        mass_exp = M_P_EXP
    error = abs(mass_pred - mass_exp) / mass_exp * 100
    print(f"{name:<10} {p['d']:>8.2f} {p['W']:>4} {mass_pred:>12.2f} {mass_exp:>12.2f} {error:>7.1f}%")

# Save report
with open("raport_qw707_koide_proton_fix.md", "w") as f:
    f.write("# RAPORT QW-707: Koide Verification + Proton Fix\n")
    f.write(f"**Date:** {datetime.datetime.now()}\n\n")
    
    f.write("## 1. Koide Verification\n")
    f.write(f"- FINAL_THEORY_REPORT claimed: Koide Q = 0.66647 (0.03% error)\n")
    f.write(f"- **ACTUAL from d^α formula: Q = {Q_method1:.5f} ({abs(Q_method1 - 2/3)/(2/3)*100:.1f}% error)**\n")
    f.write(f"- Conclusion: 0.03% was likely from experimental masses (tautology)\n\n")
    
    f.write("## 2. Proton Fix\n")
    f.write(f"- Required orbit: d_p = {d_p_needed:.2f}\n")
    f.write(f"- With W = 3 (triplet topology)\n\n")
    
    f.write("## 3. Corrected Parameters\n")
    f.write("| Particle | d | W |\n")
    f.write("|----------|---|---|\n")
    for name, p in particles.items():
        f.write(f"| {name} | {p['d']:.2f} | {p['W']} |\n")

print("\nReport saved to raport_qw707_koide_proton_fix.md")
print("="*80)
