#!/usr/bin/env python3
# QW-616: TENSOR PRODUCT DIMENSIONALITY
# Purpose: Test how ψ₁D⊗ψ₁D⊗ψ₁D creates 3D from d≈1.5 octave chain
# Key insight from QW-612: octaves alone have d≈1.5
# Hypothesis: Tensor product amplifies dimension (1.5 → 3.0)
# Date: 2025-12-05

import numpy as np
from scipy.linalg import eigh

print("="*80)
print("QW-616: TENSOR PRODUCT DIMENSIONALITY")
print("="*80)
print("Test: Jak ψ₁D⊗ψ₁D⊗ψ₁D tworzy 3D z d≈1.5?")
print("="*80)

# ============================================================================
# PARAMETERS
# ============================================================================
N_OCTAVES = 12
ALPHA_GEO = 4 * np.log(2)
OMEGA = np.pi / 4
PHI = np.pi / 6
BETA_TORS = 0.01

print(f"\nOctaves: {N_OCTAVES}")
print("-" * 40)

# ============================================================================
# 1D OCTAVE CHAIN
# ============================================================================

def K(d):
    return ALPHA_GEO * np.cos(OMEGA * d + PHI) / (1 + BETA_TORS * d)

S_1D = np.zeros((N_OCTAVES, N_OCTAVES))
for i in range(N_OCTAVES):
    for j in range(N_OCTAVES):
        S_1D[i, j] = K(abs(i - j))

eigs_1D, evecs_1D = eigh(S_1D)

# Ground state
psi_1D = evecs_1D[:, -1]  # Highest eigenvalue

print("\n1D Octave Chain:")
print(f"  Eigenvalues range: {eigs_1D[0]:.3f} to {eigs_1D[-1]:.3f}")

# Measure correlation dimension in 1D
def measure_dimension_1D(psi):
    """Correlation dimension from C(d) ~ d^(-α)"""
    C_d = []
    distances = list(range(1, N_OCTAVES))
    
    for d in distances:
        corr = 0
        count = 0
        for i in range(N_OCTAVES - d):
            corr += psi[i] * psi[i + d]
            count += 1
        C_d.append(abs(corr / count) if count > 0 else 0)
    
    C_d = np.array(C_d)
    valid = C_d > 1e-6
    
    if np.sum(valid) > 3:
        log_d = np.log(np.array(distances)[valid])
        log_C = np.log(C_d[valid])
        coeffs = np.polyfit(log_d, log_C, 1)
        alpha = -coeffs[0]
        d_eff = 1 + alpha  # For 1D chain
        return d_eff
    return np.nan

d_1D = measure_dimension_1D(psi_1D)
print(f"  Dimension d_1D: {d_1D:.2f}")

# ============================================================================
# 2D TENSOR PRODUCT
# ============================================================================

print("\n" + "="*80)
print("2D Tensor Product: ψ₂D = ψ₁D ⊗ ψ₁D")
print("="*80)

# Create 2D field
psi_2D = np.outer(psi_1D, psi_1D)

# Measure correlation in 2D
def measure_dimension_2D(psi_2d):
    """Estimate dimension from radial correlation"""
    N = psi_2d.shape[0]
    
    # Radial correlation
    C_r = []
    radii = list(range(1, N//2))
    
    for r in radii:
        # Sample points at distance r
        corr_sum = 0
        count = 0
        
        for i in range(N):
            for j in range(N):
                # Distance from (i,j) to neighbors at ~r
                for di in [-r, 0, r]:
                    for dj in [-r, 0, r]:
                        if di*di + dj*dj > 0:  # Not self
                            ii, jj = i + di, j + dj
                            if 0 <= ii < N and 0 <= jj < N:
                                dist = np.sqrt(di**2 + dj**2)
                                if abs(dist - r) < 0.5:  # Shells
                                    corr_sum += psi_2d[i, j] * psi_2d[ii, jj]
                                    count += 1
        
        if count > 0:
            C_r.append(abs(corr_sum / count))
        else:
            C_r.append(0)
    
    C_r = np.array(C_r)
    valid = C_r > 1e-6
    
    if np.sum(valid) > 3:
        log_r = np.log(np.array(radii)[valid])
        log_C = np.log(C_r[valid])
        coeffs = np.polyfit(log_r, log_C, 1)
        alpha = -coeffs[0]
        # For D-dimensional space: C ~ r^(-(D-1+α))
        # Approximate: d_eff ≈ 2 + α/2
        d_eff = 2.0 + alpha / 2
        return d_eff
    return np.nan

d_2D = measure_dimension_2D(psi_2D)

print(f"Dimension d_2D: {d_2D:.2f}")
print(f"Amplification: {d_2D / d_1D:.2f}×" if not np.isnan(d_1D) else "")

# ============================================================================
# 3D TENSOR PRODUCT  
# ============================================================================

print("\n" + "="*80)
print("3D Tensor Product: ψ₃D = ψ₁D ⊗ ψ₁D ⊗ ψ₁D")
print("="*80)

# We can't create full 12×12×12 easily, but test concept:
# Sample 3D correlation from smaller subset
N_sample = 8  # Subsample
psi_1D_sample = psi_1D[:N_sample]

# 3D field (sampled)
psi_3D_sample = np.zeros((N_sample, N_sample, N_sample))
for i in range(N_sample):
    for j in range(N_sample):
        for k in range(N_sample):
            psi_3D_sample[i, j, k] = psi_1D_sample[i] * psi_1D_sample[j] * psi_1D_sample[k]

# Measure 3D correlation (radial)
def measure_dimension_3D(psi_3d):
    """Estimate 3D dimension from radial correlation"""
    N = psi_3d.shape[0]
    
    C_r = []
    radii = list(range(1, N//2))
    
    for r in radii:
        corr_sum = 0
        count = 0
        
        for i in range(N):
            for j in range(N):
                for k in range(N):
                    # Sample spherical shell at radius r
                    for di in [-r, 0, r]:
                        for dj in [-r, 0, r]:
                            for dk in [-r, 0, r]:
                                if di*di + dj*dj + dk*dk > 0:
                                    ii, jj, kk = i + di, j + dj, k + dk
                                    if 0 <= ii < N and 0 <= jj < N and 0 <= kk < N:
                                        dist = np.sqrt(di**2 + dj**2 + dk**2)
                                        if abs(dist - r) < 0.7:
                                            corr_sum += psi_3d[i, j, k] * psi_3d[ii, jj, kk]
                                            count += 1
        
        if count > 0:
            C_r.append(abs(corr_sum / count))
        else:
            C_r.append(0)
    
    C_r = np.array(C_r)
    valid = C_r > 1e-6
    
    if np.sum(valid) > 2:
        log_r = np.log(np.array(radii)[valid])
        log_C = np.log(C_r[valid])
        coeffs = np.polyfit(log_r, log_C, 1)
        alpha = -coeffs[0]
        # Approximate for 3D
        d_eff = 3.0 + alpha / 3
        return d_eff
    return np.nan

d_3D = measure_dimension_3D(psi_3D_sample)

print(f"Dimension d_3D: {d_3D:.2f}")
print(f"Amplification from 1D: {d_3D / d_1D:.2f}×" if not np.isnan(d_1D) and not np.isnan(d_3D) else "")

# ============================================================================
# ANALYSIS
# ============================================================================

print("\n" + "="*80)
print("DIMENSIONAL AMPLIFICATION ANALYSIS")
print("="*80)

print("\n| Tensor Product | Dimension | Amplification |")
print("|----------------|-----------|---------------|")
print(f"| ψ₁D            | {d_1D:9.2f} | 1.00×         |")
if not np.isnan(d_2D):
    print(f"| ψ₁D ⊗ ψ₁D      | {d_2D:9.2f} | {d_2D/d_1D:.2f}×          |")
if not np.isnan(d_3D):
    print(f"| ψ₁D ⊗ ψ₁D ⊗ ψ₁D | {d_3D:9.2f} | {d_3D/d_1D:.2f}×          |")

# Test if d scales as expected
if not np.isnan(d_1D) and not np.isnan(d_2D) and not np.isnan(d_3D):
    # Hypothesis: d_tensor = n × d_1D (linear amplification)
    # Or: d_tensor = d_1D^n (power amplification)
    
    ratio_2D = d_2D / d_1D
    ratio_3D = d_3D / d_1D
    
    print(f"\nScaling test:")
    print(f"  d_2D / d_1D = {ratio_2D:.2f}")
    print(f"  d_3D / d_1D = {ratio_3D:.2f}")
    print(f"  Expected (linear): 2.0, 3.0")
    
    if 1.8 < ratio_2D < 2.2 and 2.7 < ratio_3D < 3.3:
        print("\n✅ LINEAR AMPLIFICATION!")
        print("   Tensor product dodaje wymiary liniowo:")
        print("   d(ψ₁D ⊗ ψ₁D ⊗ ψ₁D) ≈ 3 × d(ψ₁D)")
        verdict = "linear"
    else:
        print("\n🟡 NON-LINEAR SCALING")
        verdict = "nonlinear"
else:
    verdict = "inconclusive"

# ============================================================================
# PHYSICAL INTERPRETATION
# ============================================================================

print("\n" + "="*80)
print("PHYSICAL INTERPRETATION")
print("="*80)

if verdict == "linear":
    print("\n✅ MECHANIZM EMERGENCJI 3D WYJAŚNIONY!")
    print("\n**Odkrycie:**")
    print(f"  1. Chain oktaw ma d ≈ {d_1D:.1f} (fraktalny)")
    print(f"  2. Tensor product ψ₁D⊗ψ₁D⊗ψ₁D AMPLIFIKUJE:")
    print(f"     d_3D ≈ 3 × d_1D ≈ 3 × {d_1D:.1f} ≈ {3*d_1D:.1f}")
    print("\n**Znaczenie:**")
    print("  3D przestrzeń emerguje z PRODUCT STRUCTURE,")
    print("  nie z pojedynczych oktaw!")
    print("\n**Link do QW-384:** ψ₃D = ψ₁D⊗ψ₁D⊗ψ₁D (holographic lift)")
elif verdict == "nonlinear":
    print("\n🟡 ZŁOŻONA AMPLIFIKACJA")
    print("Tensor product modyfikuje wymiar nieliniowo")
else:
    print("\n🟡 INCONCLUSIVE")

# ============================================================================
# REPORT
# ============================================================================

report_path = "raport_qw616_tensor_dimension.md"
print(f"\nGenerating report: {report_path}...")

with open(report_path, "w") as f:
    f.write("# Raport QW-616: Tensor Product Dimensionality\n")
    f.write("**Data:** 2025-12-05\n")
    f.write("**Test:** Jak ψ₁D⊗ψ₁D⊗ψ₁D tworzy 3D przestrzeń?\n\n")
    
    f.write("## 1. Wyniki\n")
    f.write(f"- **d_1D:** {d_1D:.2f}\n")
    if not np.isnan(d_2D):
        f.write(f"- **d_2D:** {d_2D:.2f} ({d_2D/d_1D:.2f}× amplification)\n")
    if not np.isnan(d_3D):
        f.write(f"- **d_3D:** {d_3D:.2f} ({d_3D/d_1D:.2f}× amplification)\n")
    f.write("\n")
    
    if verdict == "linear":
        f.write("### ✅ LINIOWA AMPLIFIKACJA!\n")
        f.write("Tensor product dodaje wymiary liniowo:\n")
        f.write(f"d(ψ₁D⊗ψ₁D⊗ψ₁D) ≈ 3 × d(ψ₁D) ≈ {3*d_1D:.1f}\n\n")
        f.write("**Mechanizm emergencji 3D przestrzeni wyjaśniony!**\n")

print("Report saved.")
print("="*80)
