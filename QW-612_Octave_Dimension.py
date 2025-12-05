#!/usr/bin/env python3
# QW-612: OCTAVE CORRELATION DIMENSION (Proper Version)
# Purpose: Test if d_eff emerges from 12-octave correlation structure  
# Replaces QW-609 which incorrectly used spatial hopfion grid
# Proper approach: Measure correlation dimension from octave network
# Date: 2025-12-05

import numpy as np
from scipy.linalg import eigh

print("="*80)
print("QW-612: OCTAVE CORRELATION DIMENSION (Proper Octave-Based)")
print("="*80)
print("Question: What dimension emerges from 12-octave correlation structure?")
print("Hypothesis H1: Space emerges from information correlation")
print("="*80)

# ============================================================================
# PARAMETERS
# ============================================================================
N_OCTAVES = 12
ALPHA_GEO = 4 * np.log(2)
OMEGA = np.pi / 4
PHI = np.pi / 6
BETA_TORS = 0.01

print(f"\nNetwork: {N_OCTAVES} octaves")
print(f"Testing correlation dimension from eigenmode structure")
print("-" * 40)

# ============================================================================
# BUILD COUPLING MATRIX
# ============================================================================

def K(d):
    """Coupling kernel"""
    return ALPHA_GEO * np.cos(OMEGA * d + PHI) / (1 + BETA_TORS * d)

S = np.zeros((N_OCTAVES, N_OCTAVES))
for i in range(N_OCTAVES):
    for j in range(N_OCTAVES):
        S[i, j] = K(abs(i - j))

eigenvalues, eigenvectors = eigh(S)

print("\nCoupling matrix eigenvalues:")
print(f"  λ_max = {eigenvalues[-1]:.4f}")
print(f"  λ_min = {eigenvalues[0]:.4f}")
print(f"  Range = {eigenvalues[-1] - eigenvalues[0]:.4f}")

# ============================================================================
# METHOD 1: CORRELATION DIMENSION FROM OCTAVE SPACING
# ============================================================================

print("\n" + "="*80)
print("METHOD 1: Correlation Function C(d) Scaling")
print("="*80)

# Measure correlation C(d) = <ψ(i) ψ(i+d)> for ground state
ground_state = eigenvectors[:, -1]  # Highest eigenvalue

C_d = []
distances = list(range(1, N_OCTAVES))

for d in distances:
    # Correlation at distance d
    corr = 0
    count = 0
    for i in range(N_OCTAVES - d):
        corr += ground_state[i] * ground_state[i + d]
        count += 1
    C_d.append(corr / count if count > 0 else 0)

C_d = np.array(C_d)

print("\nCorrelation function C(d):")
for d, c in zip(distances[:5], C_d[:5]):
    print(f"  d={d}: C={c:+.4f}")

# Fit power law: C(d) ~ d^(-α)
# Correlation dimension: d_c relates to α
valid = np.abs(C_d) > 1e-6
if np.sum(valid) > 3:
    log_d = np.log(np.array(distances)[valid])
    log_C = np.log(np.abs(C_d[valid]))
    
    coeffs = np.polyfit(log_d, log_C, 1)
    alpha_corr = -coeffs[0]
    
    # R²
    log_C_fit = np.polyval(coeffs, log_d)
    ss_tot = np.sum((log_C - np.mean(log_C))**2)
    ss_res = np.sum((log_C - log_C_fit)**2)
    r_squared = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0
    
    print(f"\nPower law: C(d) ~ d^({coeffs[0]:.3f})")
    print(f"R² = {r_squared:.3f}")
    
    # Interpretation depends on embedding
    # For 1D chain embedded in d_eff dimensions:
    # C(d) ~ d^(-(d_eff-1)) approximately
    d_eff_method1 = 1 + alpha_corr
    
    print(f"\nEstimated d_eff (method 1): {d_eff_method1:.2f}")
else:
    d_eff_method1 = np.nan
    r_squared = 0

# ============================================================================
# METHOD 2: PARTICIPATION RATIO
# ============================================================================

print("\n" + "="*80)
print("METHOD 2: Participation Ratio (Spreading)")
print("="*80)

# Participation ratio PR = 1 / Σ|ψ_i|^4
# Measures how many octaves participate in eigenmode
# For d_eff dimensional space: PR ~ N^(d_eff/(d_eff+1))

PR_values = []

for idx in range(N_OCTAVES):
    psi = eigenvectors[:, idx]
    PR = 1.0 / np.sum(psi**4)
    PR_values.append(PR)

PR_mean = np.mean(PR_values)

print(f"\nMean Participation Ratio: {PR_mean:.2f} / {N_OCTAVES}")
print(f"Fraction: {PR_mean / N_OCTAVES:.2f}")

# For fractal dimension:
# PR ~ N^β where β = d/(d+1)
# Solve: β = PR/N → d = β/(1-β)

beta_pr = PR_mean / N_OCTAVES
if beta_pr < 1:
    d_eff_method2 = beta_pr / (1 - beta_pr)
    print(f"\nEstimated d_eff (method 2): {d_eff_method2:.2f}")
else:
    d_eff_method2 = np.nan

# ============================================================================
# METHOD 3: SPECTRAL DIMENSION FROM DENSITY OF STATES
# ============================================================================

print("\n" + "="*80)
print("METHOD 3: Spectral Dimension (Eigenvalue Density)")
print("="*80)

# For d-dimensional system: density of states ρ(E) ~ E^(d/2 - 1)
# Count eigenvalues in energy bins

# Normalize eigenvalues to [0, 1]
eigs_norm = (eigenvalues - eigenvalues[0]) / (eigenvalues[-1] - eigenvalues[0])

# Bin counting (simplified)
n_bins = 5
bins = np.linspace(0, 1, n_bins + 1)
counts, _ = np.histogram(eigs_norm, bins=bins)

# Bin centers
bin_centers = (bins[:-1] + bins[1:]) / 2

# Filter non-zero counts
valid_bins = counts > 0
if np.sum(valid_bins) > 2:
    log_E = np.log(bin_centers[valid_bins] + 1e-6)
    log_rho = np.log(counts[valid_bins])
    
    # Fit: log(ρ) = (d/2 - 1) * log(E) + const
    coeffs_dos = np.polyfit(log_E, log_rho, 1)
    spectral_exponent = coeffs_dos[0]
    
    # d = 2*(exponent + 1)
    d_eff_method3 = 2 * (spectral_exponent + 1)
    
    print(f"Density of states exponent: {spectral_exponent:.3f}")
    print(f"Estimated d_eff (method 3): {d_eff_method3:.2f}")
else:
    d_eff_method3 = np.nan

# ============================================================================
# FINAL VERDICT
# ============================================================================

print("\n" + "="*80)
print("FINAL RESULTS: Emergent Dimension from 12-Octave Network")
print("="*80)

methods = {
    'Correlation C(d)': d_eff_method1,
    'Participation Ratio': d_eff_method2,
    'Spectral Density': d_eff_method3
}

print("\n| Method | d_eff |")
print("|--------|-------|")
for name, d in methods.items():
    if not np.isnan(d):
        print(f"| {name:20s} | {d:5.2f} |")
    else:
        print(f"| {name:20s} | N/A   |")

# Average (excluding NaN)
valid_d = [d for d in methods.values() if not np.isnan(d)]
if len(valid_d) > 0:
    d_mean = np.mean(valid_d)
    d_std = np.std(valid_d)
    
    print(f"\nConsensus d_eff: {d_mean:.2f} ± {d_std:.2f}")
    
    # Compare to known results
    print("\n" + "="*80)
    print("COMPARISON")
    print("="*80)
    print(f"This test (octave network): d ≈ {d_mean:.1f}")
    print(f"QW-171 (entanglement entropy): d ≈ 2.6")
    print(f"QW-208 (compactification): d ≈ 2.6")
    print(f"QW-384 (holographic lift): ψ₃D = ψ₁D ⊗ ψ₁D ⊗ ψ₁D")
    
    if 2.4 < d_mean < 2.8:
        print("\n✅ CONSISTENT WITH d ≈ 2.6 (Fractal)")
        print("   Octave network generates fractal spatial dimension!")
        verdict = "confirmed"
    elif 2.8 < d_mean < 3.2:
        print("\n✅ d ≈ 3.0 from Octaves")
        print("   12-octave structure emergently 3D!")
        verdict = "3d"
    else:
        print(f"\n🟡 d ≈ {d_mean:.1f} (Unexpected)")
        verdict = "other"
else:
    verdict = "failed"

# ============================================================================
# REPORT
# ============================================================================

report_path = "raport_qw612_octave_dimension.md"
print(f"\nGenerating report: {report_path}...")

with open(report_path, "w") as f:
    f.write("# Raport QW-612: Octave Correlation Dimension\n")
    f.write("**Data:** 2025-12-05\n")
    f.write("**Hipoteza:** H1 - Przestrzeń emerguje z korelacji oktaw\n\n")
    
    f.write("## 1. Metodologia\n")
    f.write("Proper octave-based approach (replaces QW-609 spatial error):\n")
    f.write("- Method 1: Correlation function C(d) scaling\n")
    f.write("- Method 2: Participation ratio of eigenmodes\n")
    f.write("- Method 3: Spectral density of states\n\n")
    
    if valid_d:
        f.write("## 2. Wyniki\n")
        for name, d in methods.items():
            if not np.isnan(d):
                f.write(f"**{name}:** d_eff = {d:.2f}\n")
        f.write(f"\n**Consensus:** d = {d_mean:.2f} ± {d_std:.2f}\n\n")
        
        f.write("## 3. Interpretacja\n")
        if verdict == "confirmed":
            f.write("### ✅ ZGODNOŚĆ z d≈2.6 (Fractal)\n")
            f.write("12-oktawowa sieć generuje fraktalny wymiar przestrzenny!\n")
            f.write("Potwierdza QW-171, QW-208.\n")
        elif verdict == "3d":
            f.write("### ✅ d≈3.0 z Oktaw\n")
            f.write("Struktura 12-oktaw jest emergentnie 3D!\n")
        else:
            f.write(f"### 🟡 d≈{d_mean:.1f}\n")

print("Report saved.")
print("="*80)
