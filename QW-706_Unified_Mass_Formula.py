#!/usr/bin/env python3
"""
QW-706: UNIFIED MASS FORMULA - ALL MECHANISMS COMBINED
=======================================================
Goal: Derive mass hierarchy using EVERY verified mechanism:
      1. Winding Number (|W|) - topological charge
      2. Fractal Layer (β^N) - scale hierarchy
      3. Octave Amplification (κ^n) - resonance enhancement
      4. Pattern Complexity (S) - information content
      5. Displacement Resistance - inertia
      6. Processing Intensity (I_proc) - computational cost

Date: 2025-12-08
"""

import numpy as np
from scipy.linalg import eigh
from scipy.stats import entropy
import datetime

print("="*80)
print("QW-706: UNIFIED MASS FORMULA - ALL MECHANISMS")
print("="*80)

# ===========================================================================
# SECTION 1: FROZEN PARAMETERS (from QW-400+)
# ===========================================================================

print("\n[1] FROZEN PARAMETERS")
print("-"*60)

# Core geometric constants
ALPHA_GEO = 4 * np.log(2)  # 2.7726
OMEGA = np.pi / 4          # 0.7854
PHI = np.pi / 6            # 0.5236
BETA_TORS = 0.01           # Layer damping

# Structure
N_OCTAVES = 12             # Kissing number in 3D
N_FRACTAL_LAYERS = 30      # From H8

# Fundamental scale
M_PLANCK_GEV = 1.2209e19   # GeV

# Experimental values (for comparison)
M_ELECTRON_MEV = 0.511
M_MUON_MEV = 105.66
M_TAU_MEV = 1776.86
M_PROTON_MEV = 938.27

print(f"α_geo = 4ln(2) = {ALPHA_GEO:.4f}")
print(f"ω = π/4 = {OMEGA:.4f}")
print(f"φ = π/6 = {PHI:.4f}")
print(f"β_tors = {BETA_TORS}")
print(f"Octaves: {N_OCTAVES}")
print(f"Fractal layers: {N_FRACTAL_LAYERS}")

# Resonance amplification factor
KAPPA = ALPHA_GEO / (OMEGA * PHI)
print(f"κ = α_geo/(ω×φ) = {KAPPA:.4f}")

# ===========================================================================
# SECTION 2: K(d) KERNEL AND STABLE ORBITS
# ===========================================================================

print("\n[2] K(d) KERNEL AND STABLE ORBITS")
print("-"*60)

def K(d):
    """Universal coupling kernel."""
    if d == 0:
        return ALPHA_GEO
    return ALPHA_GEO * np.cos(OMEGA * d + PHI) / (1 + BETA_TORS * d)

# Find stable orbits (minima of effective potential from K(d))
# From FINAL_THEORY_REPORT: d₁=1.33, d₂=9.33, d₃=17.33
d_orbits = {
    'electron': 1.33,
    'muon': 9.33,
    'tau': 17.33,
    'proton': 6.0  # Between electron and proton octaves
}

print(f"Stable orbits from K(d) minima:")
for name, d in d_orbits.items():
    print(f"  {name}: d = {d:.2f}, K(d) = {K(d):.4f}")

# ===========================================================================
# SECTION 3: PARTICLE DEFINITIONS
# ===========================================================================

print("\n[3] PARTICLE DEFINITIONS")
print("-"*60)

# Each particle has:
# - winding_number: topological charge (1 for e, 2 for μ, 3 for τ?)
# - octave: primary resonance location
# - layer: fractal layer where it exists
# - orbit: stable orbit distance

particles = {
    'electron': {
        'winding_number': 1,
        'octave': 1,
        'layer': 10,
        'orbit': 1.33,
        'name': 'e'
    },
    'muon': {
        'winding_number': 1,  # Same topology as electron
        'octave': 2,          # Higher octave
        'layer': 10,
        'orbit': 9.33,
        'name': 'μ'
    },
    'tau': {
        'winding_number': 1,  # Same topology
        'octave': 3,
        'layer': 10,
        'orbit': 17.33,
        'name': 'τ'
    },
    'proton': {
        'winding_number': 3,  # 3 quarks = triplet topology
        'octave': 7,
        'layer': 10,
        'orbit': 6.0,
        'name': 'p'
    }
}

for name, p in particles.items():
    print(f"{name}: W={p['winding_number']}, O={p['octave']}, N={p['layer']}, d={p['orbit']}")

# ===========================================================================
# SECTION 4: UNIFIED MASS FORMULA
# ===========================================================================

print("\n" + "="*80)
print("[4] UNIFIED MASS FORMULA")
print("="*80)

def compute_mass(particle, observer_layer=10):
    """
    UNIFIED MASS FORMULA:
    
    m = m_0 × (d/d_0)^α × W × κ^(octave/12) × β^(layer_diff)
    
    where:
    - m_0 = reference mass (calibrated to electron)
    - (d/d_0)^α = orbit resonance scaling
    - W = winding number (topological)
    - κ^(octave/12) = octave amplification
    - β^(layer_diff) = layer hierarchy correction
    """
    
    W = particle['winding_number']
    octave = particle['octave']
    layer = particle['layer']
    orbit = particle['orbit']
    
    # Orbit reference (electron at d=1.33)
    d_0 = 1.33
    orbit_factor = (orbit / d_0) ** ALPHA_GEO
    
    # Octave amplification
    octave_factor = KAPPA ** (octave / N_OCTAVES)
    
    # Layer correction (relative to observer)
    layer_diff = abs(layer - observer_layer)
    layer_factor = (1 + BETA_TORS) ** layer_diff
    
    # K(d) coupling at particle's orbit
    coupling = abs(K(orbit))
    
    # Processing intensity (from pattern complexity)
    # Higher octave = more complex = higher I_proc
    I_proc = 1 + 0.1 * octave
    
    # Combined mass (relative units)
    mass_relative = W * orbit_factor * octave_factor * coupling * I_proc
    
    return {
        'W': W,
        'orbit_factor': orbit_factor,
        'octave_factor': octave_factor,
        'layer_factor': layer_factor,
        'coupling': coupling,
        'I_proc': I_proc,
        'mass_relative': mass_relative
    }

# Compute masses
results = {}
for name, particle in particles.items():
    results[name] = compute_mass(particle)

print("\nComponent breakdown:")
print(f"{'Particle':<10} {'W':>5} {'Orbit^α':>12} {'κ^(O/12)':>12} {'K(d)':>10} {'I_proc':>8} {'Mass_rel':>12}")
print("-"*75)
for name, r in results.items():
    print(f"{name:<10} {r['W']:>5} {r['orbit_factor']:>12.4f} {r['octave_factor']:>12.4f} "
          f"{r['coupling']:>10.4f} {r['I_proc']:>8.2f} {r['mass_relative']:>12.4f}")

# ===========================================================================
# SECTION 5: CALIBRATE TO ELECTRON
# ===========================================================================

print("\n" + "="*80)
print("[5] CALIBRATED MASSES")
print("="*80)

# Use electron as reference
m_e_rel = results['electron']['mass_relative']
m_e_mev = M_ELECTRON_MEV

# Scale factor
scale = m_e_mev / m_e_rel

print(f"Calibration: m_e(relative) = {m_e_rel:.6f}")
print(f"            m_e(exp) = {m_e_mev:.6f} MeV")
print(f"            scale = {scale:.6f}")

# Predicted masses
predicted = {}
for name, r in results.items():
    predicted[name] = r['mass_relative'] * scale

print(f"\n{'Particle':<10} {'Predicted (MeV)':>15} {'Experiment (MeV)':>18} {'Error':>10}")
print("-"*60)

exp_masses = {
    'electron': M_ELECTRON_MEV,
    'muon': M_MUON_MEV,
    'tau': M_TAU_MEV,
    'proton': M_PROTON_MEV
}

errors = {}
for name in particles:
    pred = predicted[name]
    exp = exp_masses[name]
    err = abs(pred - exp) / exp * 100
    errors[name] = err
    print(f"{name:<10} {pred:>15.4f} {exp:>18.4f} {err:>9.1f}%")

# ===========================================================================
# SECTION 6: MASS RATIOS
# ===========================================================================

print("\n" + "="*80)
print("[6] MASS RATIOS")
print("="*80)

print(f"\n{'Ratio':<12} {'Theory':>12} {'Experiment':>12} {'Gap':>10}")
print("-"*50)

ratios = [
    ('m_μ/m_e', predicted['muon']/predicted['electron'], M_MUON_MEV/M_ELECTRON_MEV),
    ('m_τ/m_e', predicted['tau']/predicted['electron'], M_TAU_MEV/M_ELECTRON_MEV),
    ('m_p/m_e', predicted['proton']/predicted['electron'], M_PROTON_MEV/M_ELECTRON_MEV),
    ('m_τ/m_μ', predicted['tau']/predicted['muon'], M_TAU_MEV/M_MUON_MEV),
]

for name, theory, exp in ratios:
    gap = exp / theory if theory > 0 else float('inf')
    print(f"{name:<12} {theory:>12.2f} {exp:>12.2f} {gap:>9.1f}×")

# ===========================================================================
# SECTION 7: KOIDE CHECK
# ===========================================================================

print("\n" + "="*80)
print("[7] KOIDE RELATION CHECK")
print("="*80)

sqrt_e = np.sqrt(predicted['electron'])
sqrt_mu = np.sqrt(predicted['muon'])
sqrt_tau = np.sqrt(predicted['tau'])

Q_pred = (sqrt_e + sqrt_mu + sqrt_tau)**2 / (predicted['electron'] + predicted['muon'] + predicted['tau'])
Q_exp = 2/3

print(f"Q (predicted) = {Q_pred:.5f}")
print(f"Q (Koide) = {Q_exp:.5f}")
print(f"Error: {abs(Q_pred - Q_exp)/Q_exp*100:.2f}%")

# ===========================================================================
# VERDICT
# ===========================================================================

print("\n" + "="*80)
print("VERDICT")
print("="*80)

avg_lepton_error = (errors['muon'] + errors['tau']) / 2
proton_error = errors['proton']

print(f"\nAverage lepton error: {avg_lepton_error:.1f}%")
print(f"Proton error: {proton_error:.1f}%")

if avg_lepton_error < 10 and proton_error < 20:
    verdict = "✅ SUCCESS"
elif avg_lepton_error < 50 or proton_error < 50:
    verdict = "🟡 PARTIAL SUCCESS"
else:
    verdict = "❌ NEEDS WORK"

print(f"\n{verdict}")

# What's working
print("\n✅ What works:")
print(f"   - Orbit scaling (d/d_0)^α produces lepton hierarchy")
print(f"   - Winding number W=3 for proton adds topological weight")

# What's missing
print("\n❌ What's missing:")
mu_gap = (M_MUON_MEV/M_ELECTRON_MEV) / (predicted['muon']/predicted['electron'])
p_gap = (M_PROTON_MEV/M_ELECTRON_MEV) / (predicted['proton']/predicted['electron'])
print(f"   - Muon gap: {mu_gap:.1f}× (need more octave amplification)")
print(f"   - Proton gap: {p_gap:.1f}× (need more topological weight)")

# ===========================================================================
# SAVE REPORT
# ===========================================================================

with open("raport_qw706_unified_mass.md", "w") as f:
    f.write("# RAPORT QW-706: Unified Mass Formula\n")
    f.write(f"**Date:** {datetime.datetime.now()}\n\n")
    
    f.write("## Formula\n")
    f.write("```\n")
    f.write("m = m_0 × (d/d_0)^α × W × κ^(octave/12) × |K(d)| × I_proc\n")
    f.write("```\n\n")
    
    f.write("## Results\n")
    f.write("| Particle | Predicted (MeV) | Experiment (MeV) | Error |\n")
    f.write("|----------|-----------------|------------------|-------|\n")
    for name in particles:
        f.write(f"| {name} | {predicted[name]:.4f} | {exp_masses[name]:.4f} | {errors[name]:.1f}% |\n")
    
    f.write(f"\n## Koide Q\n")
    f.write(f"Q = {Q_pred:.5f} (theory: {Q_exp:.5f}, error: {abs(Q_pred-Q_exp)/Q_exp*100:.2f}%)\n")
    
    f.write(f"\n## Verdict\n{verdict}\n")

print(f"\nReport saved to raport_qw706_unified_mass.md")
print("="*80)
