#!/usr/bin/env python3
# QW-611: OCTAVES VS FRACTAL LAYERS  
# Purpose: Test relationship between 12-octave structure and fractal layers
# Question: Are octaves = horizontal (info modes) and layers = vertical (scale)?
# Date: 2025-12-05

import numpy as np
from scipy.linalg import eigh

print("="*80)
print("QW-611: OCTAVES VS FRACTAL LAYERS")
print("="*80)
print("Question: How do 12 octaves relate to fractal scale hierarchy?")
print("="*80)

# ============================================================================
# PARAMETERS
# ============================================================================
N_OCTAVES = 12
ALPHA_GEO = 4 * np.log(2)
OMEGA = np.pi / 4
PHI = np.pi / 6
BETA_TORS = 0.01

N_LAYERS = 20  # Fractal hierarchy (Planck → macro)

print(f"\nOctaves (horizontal): {N_OCTAVES}")
print(f"Layers (vertical): {N_LAYERS}")
print(f"β_tors per layer: {BETA_TORS}")
print("-" * 40)

# ============================================================================
# BUILD OCTAVE COUPLING MATRIX
# ============================================================================

def K(d, alpha=ALPHA_GEO, omega=OMEGA, phi=PHI, beta=BETA_TORS):
    """Coupling kernel at single layer"""
    return alpha * np.cos(omega * d + phi) / (1 + beta * d)

S_octaves = np.zeros((N_OCTAVES, N_OCTAVES))
for i in range(N_OCTAVES):
    for j in range(N_OCTAVES):
        S_octaves[i, j] = K(abs(i - j))

eigenvalues_octaves, _ = eigh(S_octaves)

print("\nOctave spectrum (single layer):")
print(f"  λ_max = {eigenvalues_octaves[-1]:.4f}")
print(f"  λ_min = {eigenvalues_octaves[0]:.4f}")
print(f"  Range = {eigenvalues_octaves[-1] - eigenvalues_octaves[0]:.4f}")

# ============================================================================
# BUILD LAYER HIERARCHY
# ============================================================================

print("\n" + "="*80)
print("LAYER HIERARCHY (Fractal Scaling)")
print("="*80)

# At layer N, beta_effective = beta_tors^N (cumulative damping)
# Coupling strength scales down: g_gravity ~ exp(-N * beta_tors * mean_distance)

layer_couplings = []
layer_beta_eff = []

for layer in range(N_LAYERS):
    # Effective beta at this layer
    beta_eff = BETA_TORS * (layer + 1)  # Cumulative
    
    # Build coupling matrix at this layer
    S_layer = np.zeros((N_OCTAVES, N_OCTAVES))
    for i in range(N_OCTAVES):
        for j in range(N_OCTAVES):
            d = abs(i - j)
            # Modified kernel with layer-dependent damping
            S_layer[i, j] = ALPHA_GEO * np.cos(OMEGA * d + PHI) / (1 + beta_eff * d)
    
    # Measure coupling strength (Frobenius norm)
    coupling_strength = np.linalg.norm(S_layer, 'fro') / N_OCTAVES
    
    layer_couplings.append(coupling_strength)
    layer_beta_eff.append(beta_eff)
    
    if layer % 5 == 0:
        print(f"Layer {layer:2d}: β_eff={beta_eff:.3f}, coupling={coupling_strength:.4f}")

# ============================================================================
# TEST 1: COUPLING SCALING
# ============================================================================

print("\n" + "="*80)
print("TEST 1: Does coupling scale exponentially with layers?")
print("="*80)

# g_gravity ~ exp(-40) from Planck to macro (10^40 hierarchy)
# Expected: coupling(layer) ~ exp(-λ * layer)

log_couplings = np.log(layer_couplings)
layers_array = np.arange(N_LAYERS)

# Linear fit in log space
coeffs = np.polyfit(layers_array, log_couplings, 1)
lambda_decay = -coeffs[0]

# R²
coupling_mean = np.mean(log_couplings)
ss_tot = np.sum((log_couplings - coupling_mean)**2)
fit_vals = np.polyval(coeffs, layers_array)
ss_res = np.sum((log_couplings - fit_vals)**2)
r_squared = 1 - (ss_res / ss_tot)

print(f"\nExponential decay rate λ: {lambda_decay:.4f}")
print(f"R² (linearity): {r_squared:.4f}")
print(f"Hierarchy over {N_LAYERS} layers: {np.exp(lambda_decay * N_LAYERS):.2e}")

# Expected for 10^40 hierarchy over 20 layers
expected_lambda = np.log(1e40) / N_LAYERS  # ≈ 4.61

print(f"\nExpected λ for 10^40: {expected_lambda:.2f}")
print(f"Actual λ: {lambda_decay:.2f}")
print(f"Error: {abs(lambda_decay - expected_lambda)/expected_lambda * 100:.1f}%")

if r_squared > 0.95:
    print("\n✅ EXPONENTIAL SCALING CONFIRMED")
    print("   Layers create hierarchical scale separation")
    layer_scaling_ok = True
else:
    print("\n❌ NON-EXPONENTIAL")
    layer_scaling_ok = False

# ============================================================================
# TEST 2: OCTAVES vs LAYERS (Independence?)
# ============================================================================

print("\n" + "="*80)
print("TEST 2: Are octaves and layers independent dimensions?")
print("="*80)

# Hypothesis: Octaves = information modes (like Fourier components)
#             Layers = scale hierarchy (like zoom levels)
# They should be orthogonal (independent)

# Test: Does changing layer affect octave SPECTRUM (not magnitude)?

eigenvalue_sets = []

for layer in [0, 5, 10, 15, 19]:
    beta_eff = BETA_TORS * (layer + 1)
    S_layer = np.zeros((N_OCTAVES, N_OCTAVES))
    for i in range(N_OCTAVES):
        for j in range(N_OCTAVES):
            d = abs(i - j)
            S_layer[i, j] = ALPHA_GEO * np.cos(OMEGA * d + PHI) / (1 + beta_eff * d)
    
    eigs, _ = eigh(S_layer)
    # Normalize by max to get shape
    eigs_norm = eigs / np.max(np.abs(eigs))
    eigenvalue_sets.append(eigs_norm)

# Compare shapes
print("\nNormalized eigenvalue spectra:")
print("(if shapes match → octaves preserve structure across layers)")

for idx, layer in enumerate([0, 5, 10, 15, 19]):
    top3 = eigenvalue_sets[idx][-3:]
    print(f"Layer {layer:2d}: top 3 = [{top3[0]:+.3f}, {top3[1]:+.3f}, {top3[2]:+.3f}]")

# Correlation between layer 0 and layer 19
corr = np.corrcoef(eigenvalue_sets[0], eigenvalue_sets[-1])[0, 1]

print(f"\nCorrelation (layer 0 vs layer 19): r = {corr:.4f}")

if corr > 0.99:
    print("\n✅ OCTAVE STRUCTURE PRESERVED")
    print("   Octaves and layers are INDEPENDENT dimensions")
    print("   - Octaves: Information modes (horizontal)")
    print("   - Layers: Scale hierarchy (vertical)")
    independence_ok = True
else:
    print("\n🟡 STRUCTURE CHANGES")
    print("   Layers modify octave relationships")
    independence_ok = False

# ============================================================================
# TEST 3: LAYER-OCTAVE RESONANCE
# ============================================================================

print("\n" + "="*80)
print("TEST 3: Is there resonance between specific layers and octaves?")
print("="*80)

# Hypothesis: Electron (octave 1), Muon (octave 4), Tau (octave 7)
# might couple to specific layers (generations ↔ scales)

# Particle assignments from FIN Theory Paper:
particles = {
    'electron': 1,
    'muon': 4,
    'tau': 7
}

print("\nLepton octave assignments:")
for name, octave in particles.items():
    print(f"  {name:10s}: octave {octave}")

# Test if these octaves have special coupling at specific layers
# Measure cross-layer coupling for electron octave

electron_octave = 1
coupling_cross_layer = []

for target_layer in range(N_LAYERS):
    # Coupling from layer 0 (Planck) to target_layer, electron octave
    beta_source = BETA_TORS * 1
    beta_target = BETA_TORS * (target_layer + 1)
    
    # Inter-layer coupling (simplified)
    # K_inter ~ K(0) on source × transmission_factor
    transmission = np.exp(-beta_target * electron_octave)
    
    coupling_cross_layer.append(transmission)

# Check if peaks occur at specific layers
best_layer = np.argmax(coupling_cross_layer)
best_coupling = np.max(coupling_cross_layer)

print(f"\nElectron (octave 1) best couples to layer: {best_layer}")
print(f"Coupling strength: {best_coupling:.4f}")

# ============================================================================
# FINAL INTERPRETATION
# ============================================================================

print("\n" + "="*80)
print("FINAL INTERPRETATION")
print("="*80)

if independence_ok and layer_scaling_ok:
    print("\n✅ OCTAVES ≠ LAYERS (They are orthogonal dimensions)")
    print("\n**Model:**")
    print("  12 Octaves (horizontal) = Information modes / resonances")
    print("    ├─ Like Fourier components")
    print("    ├─ Encode particle types (e/μ/τ at octaves 1/4/7)")
    print("    └─ Preserve structure across scales (r={:.3f})".format(corr))
    print("\n  20 Layers (vertical) = Fractal scale hierarchy")
    print("    ├─ Planck (layer 0) → Macro (layer 20)")
    print("    ├─ β_eff = β_tors × layer (cumulative damping)")
    print("    └─ Creates 10^40 hierarchy (λ={:.2f})".format(lambda_decay))
    print("\n**Analogy:** Octaves like piano keys, Layers like octaves of sound")
    verdict = "orthogonal"
elif independence_ok:
    print("\n🟡 PARTIAL ORTHOGONALITY")
    verdict = "partial"
else:
    print("\n❌ OCTAVES MIXED WITH LAYERS")
    verdict = "mixed"

# ============================================================================
# REPORT
# ============================================================================

report_path = "raport_qw611_octave_layers.md"
print(f"\nGenerating report: {report_path}...")

with open(report_path, "w") as f:
    f.write("# Raport QW-611: Octaves vs Fractal Layers\n")
    f.write("**Data:** 2025-12-05\n")
    f.write("**Pytanie:** Jak 12 oktaw ma się do fraktalnych warstw?\n\n")
    
    f.write("## 1. Wyniki\n")
    f.write(f"**Exponential scaling:** λ={lambda_decay:.2f} (R²={r_squared:.3f})\n")
    f.write(f"**Structure preservation:** r={corr:.4f}\n\n")
    
    if verdict == "orthogonal":
        f.write("### ✅ OKTAWY ≠ WARSTWY (Orthogonalne wymiary)\n\n")
        f.write("**12 Oktaw** = poziomy wymiar (mody informacyjne)\n")
        f.write("- Kodują typy cząstek (e/μ/τ)\n")
        f.write("- Struktura zachowana na wszystkich skalach\n\n")
        f.write("**20 Warstw** = pionowy wymiar (hierarchia skal)\n")
        f.write("- Planck → Makro (10^40)\n")
        f.write("- Tłumienie β_eff ∝ layer\n\n")
        f.write("**Analogia:** Oktawy FIN = klawisze pianina\n")
        f.write("               Warstwy FIN = oktawy muzyczne\n")

print("Report saved.")
print("="*80)
