#!/usr/bin/env python3
# QW-608: SPECTRAL ORIGIN OF CONSTANTS
# Purpose: Test if fundamental constants emerge from eigenvalue spectrum
# Hypothesis H7: α_geo, ω, φ, β are geometric/spectral, not phenomenological
# Date: 2025-12-05

import numpy as np
from scipy.linalg import eigh

print("="*80)
print("QW-608: SPECTRAL ORIGIN OF CONSTANTS")
print("="*80)
print("Hypothesis H7: Constants emerge from network spectrum")
print("="*80)

# ============================================================================
# FROZEN PARAMETERS (from theory)
# ============================================================================
ALPHA_GEO = 4 * np.log(2)  # ≈ 2.7726
OMEGA = np.pi / 4           # ≈ 0.7854
PHI = np.pi / 6             # ≈ 0.5236
BETA_TORS = 0.01

N_OCTAVES = 12

print(f"\nTarget constants (from theory):")
print(f"  α_geo  = {ALPHA_GEO:.6f}")
print(f"  ω      = {OMEGA:.6f} ({OMEGA/np.pi:.4f}π)")
print(f"  φ      = {PHI:.6f} ({PHI/np.pi:.4f}π)")
print(f"  β_tors = {BETA_TORS:.6f}")
print("-" * 40)

# ============================================================================
# BUILD COUPLING MATRIX
# ============================================================================

def K(d, alpha=ALPHA_GEO, omega=OMEGA, phi=PHI, beta=BETA_TORS):
    """Universal coupling kernel"""
    return alpha * np.cos(omega * d + phi) / (1 + beta * d)

def build_matrix(N, alpha=ALPHA_GEO, omega=OMEGA, phi=PHI, beta=BETA_TORS):
    """Build self-coupling matrix S_ij = K(|i-j|)"""
    S = np.zeros((N, N))
    for i in range(N):
        for j in range(N):
            S[i, j] = K(abs(i - j), alpha, omega, phi, beta)
    return S

S = build_matrix(N_OCTAVES)

print(f"\nCoupling matrix S ({N_OCTAVES}×{N_OCTAVES}):")
print(f"Trace(S) = {np.trace(S):.6f}")
print(f"Frobenius norm = {np.linalg.norm(S, 'fro'):.6f}")

# ============================================================================
# COMPUTE EIGENVALUES
# ============================================================================

eigenvalues, eigenvectors = eigh(S)

print(f"\nEigenvalue spectrum:")
for i, ev in enumerate(eigenvalues):
    print(f"  λ_{i:2d} = {ev:+.6f}")

print(f"\nSpectral properties:")
print(f"  λ_max = {eigenvalues[-1]:.6f}")
print(f"  λ_min = {eigenvalues[0]:.6f}")
print(f"  Range = {eigenvalues[-1] - eigenvalues[0]:.6f}")
print(f"  Mean  = {np.mean(eigenvalues):.6f}")
print(f"  Std   = {np.std(eigenvalues):.6f}")

# ============================================================================
# TEST 1: α_geo from Eigenvalue Ratios
# ============================================================================

print("\n" + "="*80)
print("TEST 1: α_geo from Spectral Ratios")
print("="*80)

# Test various ratio formulas
candidates_alpha = {
    'λ_max/λ_min (abs)': abs(eigenvalues[-1] / eigenvalues[0]),
    'λ_max - λ_min': eigenvalues[-1] - eigenvalues[0],
    '(λ_max - λ_min) / N': (eigenvalues[-1] - eigenvalues[0]) / N_OCTAVES,
    'Mean(λ)': np.mean(eigenvalues),
    'Std(λ)': np.std(eigenvalues),
    'Trace(S) / N': np.trace(S) / N_OCTAVES,
    'sqrt(Trace(S²) / N)': np.sqrt(np.trace(S @ S) / N_OCTAVES),
}

print("\nCandidates for α_geo:")
print("| Formula | Value | Error vs 2.7726 |")
print("|---------|-------|-----------------|")

best_alpha = None
best_alpha_error = float('inf')

for name, value in candidates_alpha.items():
    error = abs(value - ALPHA_GEO) / ALPHA_GEO * 100
    print(f"| {name:25s} | {value:7.4f} | {error:6.2f}% |")
    if error < best_alpha_error:
        best_alpha_error = error
        best_alpha = (name, value)

print(f"\n✅ Best match: {best_alpha[0]}")
print(f"   Value: {best_alpha[1]:.6f}")
print(f"   Error: {best_alpha_error:.2f}%")

# ============================================================================
# TEST 2: ω from Eigenvalue Spacing
# ============================================================================

print("\n" + "="*80)
print("TEST 2: ω from Eigenvalue Spacing")
print("="*80)

# Compute gaps between consecutive eigenvalues
gaps = np.diff(eigenvalues)
gap_mean = np.mean(gaps)
gap_std = np.std(gaps)
gap_max = np.max(gaps)

# Test if ω relates to gap structure
candidates_omega = {
    'Mean gap': gap_mean,
    'Max gap': gap_max,
    'Std gap': gap_std,
    'π × mean gap': np.pi * gap_mean,
    'arctan(gap_max)': np.arctan(gap_max),
    'log(gap_max)': np.log(gap_max) if gap_max > 0 else np.nan,
}

print("\nCandidates for ω (target: 0.7854 = π/4):")
print("| Formula | Value | Error |")
print("|---------|-------|-------|")

best_omega = None
best_omega_error = float('inf')

for name, value in candidates_omega.items():
    if np.isnan(value):
        continue
    error = abs(value - OMEGA) / OMEGA * 100
    print(f"| {name:20s} | {value:7.4f} | {error:6.2f}% |")
    if error < best_omega_error:
        best_omega_error = error
        best_omega = (name, value)

print(f"\n{'✅' if best_omega_error < 10 else '🟡'} Best match: {best_omega[0]}")
print(f"   Value: {best_omega[1]:.6f}")
print(f"   Error: {best_omega_error:.2f}%")

# ============================================================================
# TEST 3: φ from Phase of Complex Eigenvalues
# ============================================================================

print("\n" + "="*80)
print("TEST 3: φ from Eigenvalue Phases")
print("="*80)

# For real symmetric matrix, eigenvalues are real
# But we can test if φ relates to ratios
candidates_phi = {
    'atan(λ_1 / λ_0)': np.arctan(eigenvalues[1] / eigenvalues[0]) if eigenvalues[0] != 0 else np.nan,
    'atan(λ_mid / λ_max)': np.arctan(eigenvalues[N_OCTAVES//2] / eigenvalues[-1]),
    'Mean(λ_neg) / Mean(λ_pos)': np.mean(eigenvalues[eigenvalues < 0]) / np.mean(eigenvalues[eigenvalues > 0]) if np.any(eigenvalues > 0) else np.nan,
    'π / (λ_max - λ_min)': np.pi / (eigenvalues[-1] - eigenvalues[0]),
}

print("\nCandidates for φ (target: 0.5236 = π/6):")
print("| Formula | Value | Error |")
print("|---------|-------|-------|")

best_phi = None
best_phi_error = float('inf')

for name, value in candidates_phi.items():
    if np.isnan(value):
        continue
    error = abs(value - PHI) / PHI * 100
    print(f"| {name:30s} | {value:7.4f} | {error:6.2f}% |")
    if error < best_phi_error:
        best_phi_error = error
        best_phi = (name, value)

print(f"\n{'✅' if best_phi_error < 10 else '🟡'} Best match: {best_phi[0]}")
print(f"   Value: {best_phi[1]:.6f}")
print(f"   Error: {best_phi_error:.2f}%")

# ============================================================================
# TEST 4: β from Small Eigenvalue Corrections
# ============================================================================

print("\n" + "="*80)
print("TEST 4: β_tors from Spectral Corrections")
print("="*80)

candidates_beta = {
    '1/λ_max': 1.0 / eigenvalues[-1],
    '1/(N × λ_max)': 1.0 / (N_OCTAVES * eigenvalues[-1]),
    'gap_min / gap_max': np.min(gaps) / np.max(gaps) if gap_max > 0 else np.nan,
    'Smallest |λ|': np.min(np.abs(eigenvalues)),
}

print("\nCandidates for β_tors (target: 0.01):")
print("| Formula | Value | Error |")
print("|---------|-------|-------|")

best_beta = None
best_beta_error = float('inf')

for name, value in candidates_beta.items():
    if np.isnan(value):
        continue
    error = abs(value - BETA_TORS) / BETA_TORS * 100
    print(f"| {name:20s} | {value:9.6f} | {error:6.2f}% |")
    if error < best_beta_error:
        best_beta_error = error
        best_beta = (name, value)

print(f"\n{'✅' if best_beta_error < 10 else '🟡'} Best match: {best_beta[0]}")
print(f"   Value: {best_beta[1]:.6f}")
print(f"   Error: {best_beta_error:.2f}%")

# ============================================================================
# FINAL VERDICT
# ============================================================================

print("\n" + "="*80)
print("FINAL VERDICT: H7 (Constants = Geometry)")
print("="*80)

matches = 0
if best_alpha_error < 10:
    print(f"✅ α_geo: {best_alpha[0]} (error {best_alpha_error:.1f}%)")
    matches += 1
else:
    print(f"❌ α_geo: No good match (best {best_alpha_error:.1f}%)")

if best_omega_error < 10:
    print(f"✅ ω:     {best_omega[0]} (error {best_omega_error:.1f}%)")
    matches += 1
else:
    print(f"❌ ω:     No good match (best {best_omega_error:.1f}%)")

if best_phi_error < 10:
    print(f"✅ φ:     {best_phi[0]} (error {best_phi_error:.1f}%)")
    matches += 1
else:
    print(f"❌ φ:     No good match (best {best_phi_error:.1f}%)")

if best_beta_error < 10:
    print(f"✅ β:     {best_beta[0]} (error {best_beta_error:.1f}%)")
    matches += 1
else:
    print(f"❌ β:     No good match (best {best_beta_error:.1f}%)")

print(f"\nMatches: {matches}/4")

if matches >= 3:
    print("\n✅ H7 POTWIERDZONE!")
    print("   Większość stałych wynika ze spektrum sieci")
elif matches >= 2:
    print("\n🟡 H7 CZĘŚCIOWO POTWIERDZONE")
    print("   Część stałych ma pochodzenie spektralne")
else:
    print("\n❌ H7 NIEPOTWIERDZONE")
    print("   Stałe wydają się fenomenologiczne")

# ============================================================================
# REPORT
# ============================================================================

report_path = "raport_qw608_spectral.md"
print(f"\nGenerating report: {report_path}...")

with open(report_path, "w") as f:
    f.write("# Raport QW-608: Spectral Origin of Constants\n")
    f.write("**Data:** 2025-12-05\n")
    f.write("**Hipoteza:** H7 - Stałe Fizyczne to Liczby Geometryczne/Spektralne\n\n")
    
    f.write("## 1. Metodologia\n")
    f.write(f"Budowa macierzy sprzężeń S ({N_OCTAVES}×{N_OCTAVES}) z jądra K(d).\n")
    f.write("Obliczenie eigenvalues i testowanie wszelkich możliwych formuł spektralnych.\n\n")
    
    f.write("## 2. Wyniki\n")
    f.write(f"**α_geo:** {best_alpha[0]} = {best_alpha[1]:.4f} (błąd {best_alpha_error:.1f}%)\n")
    f.write(f"**ω:**     {best_omega[0]} = {best_omega[1]:.4f} (błąd {best_omega_error:.1f}%)\n")
    f.write(f"**φ:**     {best_phi[0]} = {best_phi[1]:.4f} (błąd {best_phi_error:.1f}%)\n")
    f.write(f"**β:**     {best_beta[0]} = {best_beta[1]:.6f} (błąd {best_beta_error:.1f}%)\n\n")
    
    f.write(f"## 3. Werdykt\n")
    f.write(f"**Matches:** {matches}/4\n\n")
    
    if matches >= 3:
        f.write("### ✅ H7 POTWIERDZ ONE\n")
        f.write("Fundamentalne stałe emergują ze spektrum eigenvalues sieci!\n")
    elif matches >= 2:
        f.write("### 🟡 H7 CZĘŚCIOWO POTWIERDZONE\n")
        f.write("Część stałych ma pochodzenie spektralne, ale nie wszystkie.\n")
    else:
        f.write("### ❌ H7 NIEPOTWIERDZONE\n")
        f.write("Stałe pozostają fenomenologiczne.\n")

print("Report saved.")
print("="*80)
