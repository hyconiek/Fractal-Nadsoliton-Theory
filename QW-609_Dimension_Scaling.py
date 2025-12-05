#!/usr/bin/env python3
# QW-609: CORRELATION DIMENSION SCALING
# Purpose: Test if d_eff transitions from 2.6 (Planck) to 3.0 (macro)
# Hypothesis H1: 3D space emerges at large scales
# Date: 2025-12-05

import numpy as np
from scipy.ndimage import convolve

print("="*80)
print("QW-609: CORRELATION DIMENSION SCALING")
print("="*80)
print("Hypothesis H1: Does d=2.6 (fractal) → d=3.0 (smooth) at large scales?")
print("="*80)

# ============================================================================
# PARAMETERS
# ============================================================================
GRID_SIZE = 64
ALPHA_GEO = 4 * np.log(2)
OMEGA = np.pi / 4
PHI = np.pi / 6
BETA_TORS = 0.01

# Test different field densities
DENSITY_VALUES = [0.01, 0.1, 0.5, 1.0, 2.0]

print(f"\nGrid: {GRID_SIZE}^3")
print(f"Testing densities: {DENSITY_VALUES}")
print("-" * 40)

# ============================================================================
# CREATE HOPFION FIELD
# ============================================================================

def hopfion_field(grid_size, center, R=3.0, winding=1):
    x = np.linspace(-grid_size/2, grid_size/2, grid_size)
    X, Y, Z = np.meshgrid(x, x, x, indexing='ij')
    X, Y, Z = X - center[0], Y - center[1], Z - center[2]
    rho = np.sqrt(X**2 + Y**2)
    rho[rho == 0] = 1e-10
    eta = np.arctan2(Z, rho - R)
    xi = np.arctan2(Y, X)
    phase = winding * (xi + eta)
    dist_to_ring = np.sqrt((rho - R)**2 + Z**2)
    amplitude = np.tanh(dist_to_ring / 1.5)
    return amplitude * np.exp(1j * phase)

# ============================================================================
# MEASURE CORRELATION DIMENSION
# ============================================================================

def measure_correlation_dimension(psi, r_values):
    """Measure C(r) = <ψ(x) ψ(x+r)>"""
    field = np.abs(psi)
    N = field.shape[0]
    correlations = []
    
    for r in r_values:
        if r == 0:
            correlations.append(np.mean(field**2))
            continue
        
        # Sample correlation at distance r
        # Use periodic boundary
        shift_x = int(r)
        if shift_x >= N:
            correlations.append(0)
            continue
        
        shifted = np.roll(field, shift_x, axis=0)
        C_r = np.mean(field * shifted)
        correlations.append(C_r)
    
    return np.array(correlations)

results = []

print("\nTesting correlation dimension at different densities...")

for density in DENSITY_VALUES:
    print(f"\n{'='*60}")
    print(f"Density ρ = {density}")
    print('='*60)
    
    # Create field with specified density
    psi = hopfion_field(GRID_SIZE, center=(0, 0, 0), R=3.0, winding=+1)
    psi = psi * density  # Scale amplitude
    
    # Normalize
    psi = psi / (np.max(np.abs(psi)) + 1e-10)
    
    # Measure correlation
    r_values = np.arange(1, GRID_SIZE//2)
    C_r = measure_correlation_dimension(psi, r_values)
    
    # Fit power law: C(r) ~ r^-alpha
    # d_eff = D - alpha where D is embedding dimension
    valid = (C_r > 1e-10) & (r_values > 2)
    if np.sum(valid) > 5:
        log_r = np.log(r_values[valid])
        log_C = np.log(C_r[valid])
        
        # Linear fit
        coeffs = np.polyfit(log_r, log_C, 1)
        alpha = -coeffs[0]  # Exponent (C ~ r^-alpha)
        
        # Correlation dimension: d_c = D - alpha
        # For 3D embedding: d_eff = 3 - alpha
        # But in fractal: d_eff from correlation scaling
        d_eff = 3.0 + alpha  # If alpha < 0, d < 3
        
        # R²
        log_C_fit = np.polyval(coeffs, log_r)
        ss_tot = np.sum((log_C - np.mean(log_C))**2)
        ss_res = np.sum((log_C - log_C_fit)**2)
        r_squared = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0
        
        results.append({
            'density': density,
            'alpha': alpha,
            'd_eff': d_eff,
            'r_squared': r_squared
        })
        
        print(f"  Exponent α: {alpha:.3f}")
        print(f"  d_eff:      {d_eff:.3f}")
        print(f"  R²:         {r_squared:.3f}")
    else:
        print("  Not enough valid points for fit")

print("\n" + "="*80)
print("RESULTS SUMMARY")
print("="*80)

if len(results) > 0:
    print("\n| Density | α (exponent) | d_eff | R² |")
    print("|---------|--------------|-------|-----|")
    for r in results:
        print(f"| {r['density']:7.2f} | {r['alpha']:12.3f} | {r['d_eff']:5.2f} | {r['r_squared']:.3f} |")
    
    # Check if d_eff increases with density
    densities = np.array([r['density'] for r in results])
    d_values = np.array([r['d_eff'] for r in results])
    
    if len(densities) > 2:
        # Correlation between density and d_eff
        from scipy.stats import pearsonr
        corr, p_value = pearsonr(densities, d_values)
        
        print(f"\nCorrelation (density vs d_eff): r = {corr:.3f}")
        
        # Check if d approaches 3.0
        d_max = np.max(d_values)
        d_min = np.min(d_values)
        
        print(f"d_eff range: {d_min:.2f} → {d_max:.2f}")
        
        if corr > 0.7 and d_max > 2.9:
            print("\n✅ SCALE TRANSITION CONFIRMED!")
            print(f"   d emerguje → 3.0 przy wysokiej gęstości (ρ={densities[np.argmax(d_values)]:.1f})")
            print("   FIN space staje się 3D na makroskali!")
            verdict = "confirmed"
        elif corr > 0.5:
            print("\n🟡 PARTIAL TRANSITION")
            print(f"   d rośnie z gęstością (r={corr:.2f}), ale nie osiąga 3.0")
            verdict = "partial"
        else:
            print("\n❌ NO TRANSITION")
            print("   d pozostaje fractal (d≈2.6) na wszystkich skalach")
            verdict = "none"
    else:
        verdict = "insufficient"
else:
    verdict = "failed"

# ============================================================================
# REPORT
# ============================================================================

report_path = "raport_qw609_dimension.md"
print(f"\nGenerating report: {report_path}...")

with open(report_path, "w") as f:
    f.write("# Raport QW-609: Correlation Dimension Scaling\n")
    f.write("**Data:** 2025-12-05\n")
    f.write("**Hipoteza:** H1 - Czy d=2.6 przechodzi w d=3.0 na dużych skalach?\n\n")
    
    f.write("## 1. Metodologia\n")
    f.write("Test korelacji C(r) ~ r^-α przy różnych gęstościach pola.\n")
    f.write("Wymiar: d_eff = 3 + α (α<0 dla fraktalu)\n\n")
    
    if len(results) > 0:
        f.write("## 2. Wyniki\n")
        f.write("| Density | α | d_eff |\n|---------|---|-------|\n")
        for r in results:
            f.write(f"| {r['density']:.2f} | {r['alpha']:.3f} | {r['d_eff']:.2f} |\n")
        f.write("\n")
        
        f.write("## 3. Analiza\n")
        if verdict == "confirmed":
            f.write("### ✅ SCALE TRANSITION POTWIERDZONA\n")
            f.write("Wymiar emergentny rośnie z gęstością i osiąga d≈3.0!\n")
            f.write("FIN space ma **fractal microstructure** ale **smooth macrostructure**.\n")
        elif verdict == "partial":
            f.write("### 🟡 CZĘŚCIOWA TRANSITION\n")
            f.write("d rośnie z gęstością, ale nie osiąga pełnego d=3.0.\n")
        else:
            f.write("### ❌ BRAK TRANSITION\n")
            f.write("d pozostaje fractal na wszystkich skalach.\n")

print("Report saved.")
print("="*80)
