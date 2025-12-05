#!/usr/bin/env python3
# QW-618: SUPER-BALLISTIC NOISE CHECK (RK4)
# Purpose: Definitive test if super-ballistic (b≈2.4) survives noise
# Improvement: Uses RK4 integrator to eliminate Euler damping artifacts
# Date: 2025-12-05

import numpy as np

print("="*80)
print("QW-618: SUPER-BALLISTIC NOISE CHECK (RK4)")
print("="*80)
print("Test: Czy b≈2.4 przetrwa szum przy dokładniejszym całkowaniu?")
print("="*80)

# ============================================================================
# PARAMETERS
# ============================================================================
GRID_SIZE = 64  # Increased for better resolution
DX = 0.5
DT = 0.005      # Smaller timestep
STEPS = 400     # Longer evolution

ALPHA_GEO = 4 * np.log(2)
BETA_TORS = 0.01
GAMMA = 0.1     # Attracteur

NOISE_LEVELS = [0.0, 0.01, 0.05, 0.1, 0.2]

print(f"Grid: {GRID_SIZE}")
print(f"Integration: RK4 (dt={DT})")
print(f"Noise levels: {NOISE_LEVELS}")
print("-" * 40)

# ============================================================================
# COUPLING
# ============================================================================

def K(d):
    return ALPHA_GEO / (1 + BETA_TORS * d)

S = np.zeros((GRID_SIZE, GRID_SIZE))
for i in range(GRID_SIZE):
    for j in range(GRID_SIZE):
        d = abs(i - j)
        S[i, j] = K(d)

# ============================================================================
# DYNAMICS (RK4)
# ============================================================================

def derivatives(psi, t, noise_field):
    # Hamiltonian part: i dpsi/dt = - [Coupling + Nonlinear] psi
    # dpsi/dt = i * (S @ psi - gamma |psi|^2 psi) + noise
    
    coupling = np.dot(S, psi)
    nonlinear = -GAMMA * np.abs(psi)**2 * psi
    
    dpsi = 1j * (coupling + nonlinear) + noise_field
    return dpsi

results = []

for noise_level in NOISE_LEVELS:
    print(f"\nProcessing noise σ = {noise_level}...")
    
    # Initial Gaussian
    psi = np.zeros(GRID_SIZE, dtype=complex)
    center = GRID_SIZE // 2
    sigma_0 = 3.0
    x = (np.arange(GRID_SIZE) - center) * DX
    psi = np.exp(-x**2 / (2 * sigma_0**2)) + 0j
    psi = psi / np.sqrt(np.sum(np.abs(psi)**2) * DX)
    
    widths = []
    times = []
    
    # Pre-generate noise field (static thermal background or dynamic?)
    # For robust test, let's use dynamic noise refreshed every step
    
    for step in range(STEPS):
        t = step * DT
        
        # Measure
        rho = np.abs(psi)**2
        x_mean = np.sum(x * rho) * DX
        sigma_t = np.sqrt(np.sum((x - x_mean)**2 * rho) * DX)
        
        widths.append(sigma_t)
        times.append(t)
        
        # RK4 Step
        # Generate noise for this step
        if noise_level > 0:
            noise = (np.random.randn(GRID_SIZE) + 1j*np.random.randn(GRID_SIZE)) * noise_level
        else:
            noise = 0
            
        k1 = derivatives(psi, t, noise)
        k2 = derivatives(psi + 0.5*DT*k1, t + 0.5*DT, noise)
        k3 = derivatives(psi + 0.5*DT*k2, t + 0.5*DT, noise)
        k4 = derivatives(psi + DT*k3, t + DT, noise)
        
        psi = psi + (DT/6.0) * (k1 + 2*k2 + 2*k3 + k4)
        
        # Renormalize (if open system, norm might change, but for dispersion we check shape)
        # Re-normalizing keeps it stable against numerical drift
        psi = psi / np.sqrt(np.sum(np.abs(psi)**2) * DX + 1e-12)
        
    # Fit power law σ(t) ~ t^b
    widths = np.array(widths)
    times = np.array(times)
    
    # Filter initial transient
    mask = (times > 0.5) & (widths > sigma_0)
    
    if np.sum(mask) > 10:
        log_t = np.log(times[mask])
        log_w = np.log(widths[mask])
        
        coeffs = np.polyfit(log_t, log_w, 1)
        b = coeffs[0]
        
        # R2
        fit = np.polyval(coeffs, log_t)
        ss_res = np.sum((log_w - fit)**2)
        ss_tot = np.sum((log_w - np.mean(log_w))**2)
        r2 = 1 - (ss_res / ss_tot)
        
        print(f"  Exponent b: {b:.3f} (R²={r2:.3f})")
        
        results.append({'noise': noise_level, 'b': b, 'r2': r2})
    else:
        print("  Insufficient valid points for fit")

# ============================================================================
# ANALYSIS
# ============================================================================

print("\n" + "="*80)
print("FINAL RESULTS")
print("="*80)
print("\n| Noise | Exponent b | R² | Status |")
print("|-------|------------|-----|--------|")

for res in results:
    # Check if super-ballistic (b > 1)
    # Ideally b ≈ 2.0 - 2.5
    status = "✅ Robust" if res['b'] > 1.5 else ("🟡 Weak" if res['b'] > 1.1 else "❌ Normal/Sub")
    print(f"| {res['noise']:5.2f} | {res['b']:10.3f} | {res['r2']:.3f} | {status} |")

# Determine verdict
robust_count = sum(1 for r in results if r['b'] > 1.5)
if robust_count >= len(results) - 1: # Allow fail at highest noise
    verdict = "robust"
    print("\n✅ SUPER-BALLISTIC CONFIRMED ROBUST (RK4)")
    print("   Dynamika przetrwała szum dzięki lepszemu integratorowi!")
else:
    verdict = "fragile"
    print("\n❌ STILL FRAGILE OR DAMPED")

# Save report
with open("raport_qw618_superballistic_rk4.md", "w") as f:
    f.write("# Raport QW-618: Super-Ballistic Noise Check (RK4)\n")
    f.write("**Data:** 2025-12-05\n\n")
    f.write("## Wyniki (Exponent b)\n")
    f.write("| Noise | b | R² |\n|---|---|---|\n")
    for r in results:
        f.write(f"| {r['noise']} | {r['b']:.3f} | {r['r2']:.3f} |\n")
    
    if verdict == "robust":
        f.write("\n### ✅ Wnioski: Robust!\n")
        f.write("Użycie RK4 potwierdza, że b > 1.5 (super-ballistic) jest fizyczne i odporne na szum.\n")
    else:
        f.write("\n### ❌ Wnioski: Problematyczne\n")

print("Report saved.")
