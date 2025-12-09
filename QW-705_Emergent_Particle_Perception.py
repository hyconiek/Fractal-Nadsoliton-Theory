#!/usr/bin/env python3
"""
QW-705: EMERGENT PARTICLE PERCEPTION
=====================================
Hypothesis: Particles are NOT fundamental objects, but stable patterns
            that the Nadsoliton "sees" when observing itself.
            
Key Insight: We (observers) ARE part of the Nadsoliton.
             "Particles" are how local correlations APPEAR from inside.
             Mass = resistance to displacement FROM INTERNAL perspective.

Goal: Define "emergent observer" mathematically and derive how mass
      emerges as a property of internal perception.

Date: 2025-12-08
"""

import numpy as np
from scipy.linalg import eigh
from scipy.stats import entropy
import datetime

print("="*80)
print("QW-705: EMERGENT PARTICLE PERCEPTION")
print("="*80)
print("Hypothesis: Particles are emergent perceptions of self-interacting Nadsoliton")
print("="*80)

# ===========================================================================
# SECTION 1: NADSOLITON AS SELF-INTERACTING PROCESS
# ===========================================================================

print("\n[1] NADSOLITON: ONE SELF-INTERACTING PROCESS")
print("-"*60)

# Frozen parameters
N_OCTAVES = 12
ALPHA_GEO = 4 * np.log(2)
BETA_TORS = 0.01
OMEGA = np.pi / 4
PHI = np.pi / 6

def K(d):
    """Coupling kernel - how octaves interact with each other."""
    if d == 0:
        return ALPHA_GEO
    return ALPHA_GEO * np.cos(OMEGA * d + PHI) / (1 + BETA_TORS * d)

# Build coupling matrix
H_coupling = np.zeros((N_OCTAVES, N_OCTAVES))
for i in range(N_OCTAVES):
    for j in range(N_OCTAVES):
        H_coupling[i, j] = K(abs(i - j))

print(f"Nadsoliton: {N_OCTAVES} octaves coupled by K(d)")
print(f"Total couplings: {N_OCTAVES**2}")

# Eigenstructure - these are the "natural modes" of the system
eigenvalues, eigenvectors = eigh(H_coupling)
print(f"\nNatural modes (eigenvalues):")
for i, ev in enumerate(eigenvalues[:4]):
    print(f"  Mode {i}: λ = {ev:.4f}")

# ===========================================================================
# SECTION 2: DEFINE "EMERGENT OBSERVER"
# ===========================================================================

print("\n" + "="*80)
print("[2] EMERGENT OBSERVER DEFINITION")
print("-"*60)

# Key insight: An "observer" is a localized subsystem within the Nadsoliton
# that can record/process correlations with other parts.

# Observer = subset of octaves that form a coherent pattern
OBSERVER_OCTAVES = [5, 6, 7]  # Central octaves (our scale)
print(f"Observer: localized in octaves {OBSERVER_OCTAVES}")

# Observer's "state" is the projection of full system onto observer's subspace
def observer_projection(full_state, observer_octaves):
    """Project full Nadsoliton state onto observer's subspace."""
    proj = np.zeros_like(full_state)
    for o in observer_octaves:
        proj[o] = full_state[o]
    return proj / (np.linalg.norm(proj) + 1e-10)

# Observer sees REDUCED density matrix
def observer_reduced_density(full_density, observer_octaves):
    """Trace out non-observer degrees of freedom."""
    n_obs = len(observer_octaves)
    rho_reduced = np.zeros((n_obs, n_obs), dtype=complex)
    for i, oi in enumerate(observer_octaves):
        for j, oj in enumerate(observer_octaves):
            rho_reduced[i, j] = full_density[oi, oj]
    # Normalize
    trace = np.trace(rho_reduced)
    if np.abs(trace) > 1e-10:
        rho_reduced /= trace
    return rho_reduced

print("Observer perceives REDUCED state (traced over non-local octaves)")

# ===========================================================================
# SECTION 3: "PARTICLE" AS STABLE CORRELATION PATTERN
# ===========================================================================

print("\n" + "="*80)
print("[3] PARTICLE = STABLE CORRELATION PATTERN")
print("-"*60)

# A "particle" is NOT an object, but a stable pattern of correlations
# that the observer repeatedly perceives.

# Create a localized excitation (candidate "particle")
def create_particle_pattern(center_octave, width=1.5):
    """Create Gaussian-localized pattern in octave space."""
    pattern = np.exp(-0.5 * ((np.arange(N_OCTAVES) - center_octave) / width)**2)
    return pattern / np.linalg.norm(pattern)

# Electron pattern: localized at octave 1
electron_pattern = create_particle_pattern(1, width=1.0)
print(f"'Electron' pattern: peak at octave 1")

# Muon pattern: localized at octave 2 with different structure
muon_pattern = create_particle_pattern(2, width=1.2)
print(f"'Muon' pattern: peak at octave 2")

# Proton pattern: spread across 3 octaves (triplet)
proton_pattern = create_particle_pattern(7, width=1.5)
print(f"'Proton' pattern: peak at octave 7")

# ===========================================================================
# SECTION 4: WHAT OBSERVER "SEES"
# ===========================================================================

print("\n" + "="*80)
print("[4] WHAT OBSERVER 'SEES'")
print("-"*60)

# Observer at octaves 5,6,7 "sees" other patterns through correlations
# The "visibility" depends on overlap with K(d) propagator

def observer_sees_pattern(observer_octaves, particle_pattern):
    """Calculate how strongly observer perceives a pattern."""
    # Correlation through K(d)
    visibility = 0.0
    for o_obs in observer_octaves:
        for o_part in range(N_OCTAVES):
            visibility += particle_pattern[o_part] * K(abs(o_obs - o_part))
    return visibility

vis_electron = observer_sees_pattern(OBSERVER_OCTAVES, electron_pattern)
vis_muon = observer_sees_pattern(OBSERVER_OCTAVES, muon_pattern)
vis_proton = observer_sees_pattern(OBSERVER_OCTAVES, proton_pattern)

print(f"Observer visibility:")
print(f"  'Electron': {vis_electron:.4f}")
print(f"  'Muon':     {vis_muon:.4f}")
print(f"  'Proton':   {vis_proton:.4f}")

# ===========================================================================
# SECTION 5: MASS = RESISTANCE TO DISPLACEMENT (FROM OBSERVER'S VIEW)
# ===========================================================================

print("\n" + "="*80)
print("[5] MASS = RESISTANCE TO DISPLACEMENT")
print("-"*60)

# Key insight: Mass is NOT stored energy, but how hard it is to MOVE a pattern
# from the observer's perspective.

# "Moving" means shifting the pattern in octave space
# Resistance = energy cost of displacement

def displacement_energy(pattern, shift):
    """Calculate energy cost to shift pattern by 'shift' octaves."""
    shifted = np.roll(pattern, shift)
    # Energy = change in K-weighted self-energy
    E_original = np.dot(pattern, H_coupling @ pattern)
    E_shifted = np.dot(shifted, H_coupling @ shifted)
    return abs(E_shifted - E_original)

def effective_mass(pattern):
    """
    Mass = d²E/dx² at equilibrium (curvature of energy landscape).
    Approximated as energy cost for unit shift.
    """
    E_plus = displacement_energy(pattern, 1)
    E_minus = displacement_energy(pattern, -1)
    # Second derivative approximation
    mass = (E_plus + E_minus) / 2
    return mass

m_electron = effective_mass(electron_pattern)
m_muon = effective_mass(muon_pattern)
m_proton = effective_mass(proton_pattern)

print(f"Effective mass (resistance to displacement):")
print(f"  'Electron': m_e = {m_electron:.6f}")
print(f"  'Muon':     m_μ = {m_muon:.6f}")
print(f"  'Proton':   m_p = {m_proton:.6f}")

print(f"\nMass ratios:")
print(f"  m_μ/m_e = {m_muon/m_electron:.2f}")
print(f"  m_p/m_e = {m_proton/m_electron:.2f}")

print(f"\nExperimental ratios:")
print(f"  m_μ/m_e = 207")
print(f"  m_p/m_e = 1836")

# ===========================================================================
# SECTION 6: COMPLEXITY = MASS
# ===========================================================================

print("\n" + "="*80)
print("[6] COMPLEXITY OF PATTERN → MASS")
print("-"*60)

# Alternative: Mass scales with COMPLEXITY of pattern (information content)

def pattern_complexity(pattern):
    """Measure information complexity of pattern."""
    # Shannon entropy
    prob = np.abs(pattern)**2
    prob = prob / np.sum(prob)
    S = entropy(prob, base=2)
    return S

def pattern_energy(pattern):
    """Total K-weighted energy of pattern."""
    return np.dot(pattern, H_coupling @ pattern)

def mass_from_complexity(pattern):
    """
    Mass = (Pattern Complexity) × (Coupling Strength)
    
    Interpretation: More complex patterns require more "computational cycles"
    to maintain against diffusion. This is perceived as "mass" by observer.
    """
    comp = pattern_complexity(pattern)
    energy = abs(pattern_energy(pattern))
    return comp * energy

mass_e = mass_from_complexity(electron_pattern)
mass_mu = mass_from_complexity(muon_pattern)
mass_p = mass_from_complexity(proton_pattern)

print(f"Pattern complexity:")
print(f"  S(electron) = {pattern_complexity(electron_pattern):.4f} bits")
print(f"  S(muon)     = {pattern_complexity(muon_pattern):.4f} bits")
print(f"  S(proton)   = {pattern_complexity(proton_pattern):.4f} bits")

print(f"\nMass from complexity:")
print(f"  m_e = {mass_e:.4f}")
print(f"  m_μ = {mass_mu:.4f}")
print(f"  m_p = {mass_p:.4f}")

print(f"\nRatios:")
print(f"  m_μ/m_e = {mass_mu/mass_e:.2f}")
print(f"  m_p/m_e = {mass_p/mass_e:.2f}")

# ===========================================================================
# SECTION 7: HIERARCHY FROM OBSERVER'S LAYER
# ===========================================================================

print("\n" + "="*80)
print("[7] HIERARCHY FROM OBSERVER'S LAYER")
print("-"*60)

# Key insight: Mass hierarchy depends on WHERE the observer is located
# Different layers see different ratios!

def mass_as_seen_from_layer(pattern, observer_layer):
    """
    Mass depends on observer's layer relative to particle's layer.
    
    Mechanism: Patterns further from observer appear "heavier" because
    correlations are weaker (harder to couple to / move).
    """
    pattern_layer = np.argmax(pattern)  # Peak location
    layer_distance = abs(observer_layer - pattern_layer)
    
    # Base mass from complexity
    base_mass = mass_from_complexity(pattern)
    
    # Layer scaling: further patterns appear heavier
    layer_factor = (1 + BETA_TORS * layer_distance) ** ALPHA_GEO
    
    return base_mass * layer_factor

# Observer at layer 6 (central)
obs_layer = 6

mass_e_seen = mass_as_seen_from_layer(electron_pattern, obs_layer)
mass_mu_seen = mass_as_seen_from_layer(muon_pattern, obs_layer)
mass_p_seen = mass_as_seen_from_layer(proton_pattern, obs_layer)

print(f"Observer at octave {obs_layer}")
print(f"\nMass AS SEEN BY OBSERVER:")
print(f"  m_e = {mass_e_seen:.4f}")
print(f"  m_μ = {mass_mu_seen:.4f}")
print(f"  m_p = {mass_p_seen:.4f}")

print(f"\nRatios AS SEEN:")
print(f"  m_μ/m_e = {mass_mu_seen/mass_e_seen:.2f}")
print(f"  m_p/m_e = {mass_p_seen/mass_e_seen:.2f}")

# ===========================================================================
# VERDICT
# ===========================================================================

print("\n" + "="*80)
print("VERDICT: EMERGENT PARTICLE PERCEPTION")
print("="*80)

print("""
KEY INSIGHTS:

1. PARTICLES ARE NOT FUNDAMENTAL
   They are stable correlation patterns in the Nadsoliton.
   
2. OBSERVER IS PART OF THE SYSTEM
   We observe from INSIDE, not outside.
   
3. MASS = RESISTANCE TO DISPLACEMENT
   From internal perspective, moving a pattern requires energy.
   This appears as "mass" (inertia).
   
4. COMPLEXITY CONTRIBUTES TO MASS
   More complex patterns are harder to maintain/move.
   s
5. HIERARCHY DEPENDS ON OBSERVER POSITION
   Patterns far from observer appear heavier.
""")

# Check if this explains hierarchy
ratio_mu_e = mass_mu_seen / mass_e_seen
ratio_p_e = mass_p_seen / mass_e_seen

print(f"Results vs Experiment:")
print(f"| Ratio | Theory | Experiment | Gap |")
print(f"|-------|--------|------------|-----|")
print(f"| m_μ/m_e | {ratio_mu_e:.1f} | 207 | {207/ratio_mu_e:.1f}× |")
print(f"| m_p/m_e | {ratio_p_e:.1f} | 1836 | {1836/ratio_p_e:.1f}× |")

if ratio_mu_e > 10 and ratio_p_e > 100:
    print("\n✅ PARTIAL SUCCESS: Mechanism produces hierarchical ratios")
else:
    print("\n🟡 Mechanism exists but does not produce large enough ratios")
    print("   Need: Better localization model or additional winding numbers")

print("\n" + "="*80)
print("NEXT STEPS:")
print("- Integrate winding number (topological complexity)")
print("- Use fractal layers (30 layers) for deeper hierarchy")
print("- Calculate processing intensity from pattern dynamics")
print("="*80)

# Save report
with open("raport_qw705_emergent_particle_perception.md", "w") as f:
    f.write("# RAPORT QW-705: Emergent Particle Perception\n")
    f.write(f"**Date:** {datetime.datetime.now()}\n\n")
    f.write("## Key Insight\n")
    f.write("Particles are NOT fundamental objects but stable patterns\n")
    f.write("that the Nadsoliton 'sees' when observing itself.\n\n")
    f.write("## Mass Mechanism\n")
    f.write("Mass = Resistance to displacement from internal perspective\n")
    f.write("     = Pattern complexity × Coupling strength × Layer factor\n\n")
    f.write("## Results\n")
    f.write(f"- m_μ/m_e = {ratio_mu_e:.1f} (exp: 207)\n")
    f.write(f"- m_p/m_e = {ratio_p_e:.1f} (exp: 1836)\n\n")
    f.write("## Verdict\n")
    f.write("Mechanism produces hierarchy but not at experimental scale.\n")
    f.write("Requires integration with winding numbers and fractal layers.\n")

print("Report saved to raport_qw705_emergent_particle_perception.md")
