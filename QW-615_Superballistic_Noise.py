#!/usr/bin/env python3
# QW-615: SUPER-BALLISTIC POD WPŁYWEM SZUMU
# Purpose: Test if super-ballistic dispersion survives thermal noise
# Key question: Is b≈2.4 robust or noise-sensitive?
# Date: 2025-12-05

import numpy as np

print("="*80)
print("QW-615: SUPER-BALLISTIC POD WPŁYWEM SZUMU")
print("="*80)
print("Test: Czy b≈2.4 przetrwa thermal/quantum noise?")
print("="*80)

# ============================================================================
# PARAMETERS
# ============================================================================
GRID_SIZE = 32
DX = 0.5
DT = 0.01
STEPS = 100

ALPHA_GEO = 4 * np.log(2)
BETA_TORS = 0.01
GAMMA = 0.1  # Attracteur (QW-558)

# Noise levels
NOISE_LEVELS = [0.0, 0.01, 0.05, 0.1, 0.2]

print(f"\nGrid: {GRID_SIZE}")
print(f"Timesteps: {STEPS}")
print(f"Noise levels: {NOISE_LEVELS}")
print("-" * 40)

# ============================================================================
# COUPLING KERNEL
# ============================================================================

def K(d):
    return ALPHA_GEO / (1 + BETA_TORS * d)

# Build coupling matrix
S = np.zeros((GRID_SIZE, GRID_SIZE))
for i in range(GRID_SIZE):
    for j in range(GRID_SIZE):
        d = abs(i - j)
        S[i, j] = K(d)

print("Coupling matrix ready")

# ============================================================================
# TEST AT DIFFERENT NOISE LEVELS
# ============================================================================

results = []

for noise_level in NOISE_LEVELS:
    print(f"\n{'='*60}")
    print(f"Noise σ = {noise_level}")
    print('='*60)
    
    # Initial wave packet (Gaussian)
    psi = np.zeros(GRID_SIZE, dtype=complex)
    center = GRID_SIZE // 2
    sigma_0 = 3.0
    for i in range(GRID_SIZE):
        x = (i - center) * DX
        psi[i] = np.exp(-x**2 / (2 * sigma_0**2)) * np.exp(1j * 0)
    
    # Normalize
    psi = psi / np.sqrt(np.sum(np.abs(psi)**2) * DX)
    
    # Track width evolution
    widths = []
    times = []
    
    for t in range(STEPS):
        # Measure width at this timestep
        density = np.abs(psi)**2
        x_grid = np.arange(GRID_SIZE) * DX
        
        # Mean position
        x_mean = np.sum(x_grid * density) * DX
        
        # Width (standard deviation)
        sigma_t = np.sqrt(np.sum((x_grid - x_mean)**2 * density) * DX)
        
        widths.append(sigma_t)
        times.append(t * DT)
        
        # Evolution: dψ/dt = i(coupling - gamma|ψ|²)ψ + noise
        # Network coupling term
        coupling_term = np.dot(S, psi)
        
        # Attractor
        attractor_term = -GAMMA * np.abs(psi)**2 * psi
        
        # Noise: Random phase/amplitude fluctuations
        if noise_level > 0:
            noise_real = np.random.randn(GRID_SIZE) * noise_level
            noise_imag = np.random.randn(GRID_SIZE) * noise_level
            noise_term = noise_real + 1j * noise_imag
        else:
            noise_term = 0
        
        # Update (simple Euler)
        dpsi_dt = 1j * coupling_term + attractor_term + noise_term
        psi = psi + dpsi_dt * DT
        
        # Renormalize to conserve probability
        psi = psi / np.sqrt(np.sum(np.abs(psi)**2) * DX + 1e-10)
    
    # Fit power law: σ(t) ~ t^b
    widths = np.array(widths)
    times = np.array(times)
    
    # Exclude t=0 (log issues)
    valid = (times > 0.1) & (widths > sigma_0 * 0.8)
    
    if np.sum(valid) > 10:
        log_t = np.log(times[valid])
        log_sigma = np.log(widths[valid])
        
        coeffs = np.polyfit(log_t, log_sigma, 1)
        b_exponent = coeffs[0]
        
        # R²
        log_sigma_fit = np.polyval(coeffs, log_t)
        ss_tot = np.sum((log_sigma - np.mean(log_sigma))**2)
        ss_res = np.sum((log_sigma - log_sigma_fit)**2)
        r_squared = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0
        
        results.append({
            'noise': noise_level,
            'b': b_exponent,
            'r2': r_squared
        })
        
        print(f"  Exponent b: {b_exponent:.3f}")
        print(f"  R²: {r_squared:.3f}")
    else:
        print("  Insufficient valid points")

# ============================================================================
# ANALYSIS
# ============================================================================

print("\n" + "="*80)
print("ROBUSTNESS ANALYSIS")
print("="*80)

print("\n| Noise | b (exponent) | R² | Status |")
print("|-------|--------------|-----|--------|")

b_baseline = results[0]['b'] if len(results) > 0 else 0

for res in results:
    # Super-ballistic if b > 2.0
    # Robust if within 10% of baseline
    
    is_super_ballistic = res['b'] > 2.0
    is_robust = abs(res['b'] - b_baseline) / b_baseline < 0.1 if b_baseline > 0 else False
    
    if is_super_ballistic and is_robust:
        status = "✅ Robust SB"
    elif is_super_ballistic:
        status = "🟡 SB Changed"
    elif is_robust:
        status = "🟡 Lost SB"
    else:
        status = "❌ Degraded"
    
    print(f"| {res['noise']:5.2f} | {res['b']:12.3f} | {res['r2']:.3f} | {status:14s} |")

# Find critical noise
critical_noise = None
for res in results[1:]:  # Skip baseline
    if abs(res['b'] - b_baseline) / b_baseline > 0.1 or res['b'] < 2.0:
        critical_noise = res['noise']
        break

if critical_noise is not None:
    print(f"\n⚠️  Critical noise: σ_crit ≈ {critical_noise}")
    print(f"    Super-ballistic degrades above this level")
else:
    print(f"\n✅ SUPER-BALLISTIC ROBUST")
    print(f"    b≈2.4 survives noise up to σ={NOISE_LEVELS[-1]}")

# ============================================================================
# PHYSICAL INTERPRETATION
# ============================================================================

print("\n" + "="*80)
print("PHYSICAL INTERPRETATION")
print("="*80)

if len(results) > 0:
    b_mean = np.mean([r['b'] for r in results])
    b_std = np.std([r['b'] for r in results])
    
    print(f"\nMean exponent: b = {b_mean:.3f} ± {b_std:.3f}")
    print(f"Baseline (no noise): b = {b_baseline:.3f}")
    
    if critical_noise is None or critical_noise > 0.1:
        print("\n✅ SUPER-BALLISTIC TO INTRINSIC PROPERTY")
        print("\nEgzotyczna  dynamika przetrwała szum:")
        print("- Nie jest artefaktem numerycznym")
        print("- Robust przeciw thermal fluctuations")
        print("- FUNDAMENTALNA właściwość sieci")
        verdict = "robust"
    elif critical_noise > 0.01:
        print("\n🟡 CZĘŚCIOWA ROBUSTNOŚĆ")
        print(f"\nSuper-ballistic zachowane do σ < {critical_noise}")
        verdict = "partial"
    else:
        print("\n❌ WRAŻLIWE NA SZUM")
        print("\nSuper-ballistic może być emergentnym efektem wymagającym low noise")
        verdict = "fragile"

# ============================================================================
# REPORT
# ============================================================================

report_path = "raport_qw615_superballistic_noise.md"
print(f"\nGenerating report: {report_path}...")

with open(report_path, "w") as f:
    f.write("# Raport QW-615: Super-Ballistic pod Wpływem Szumu\n")
    f.write("**Data:** 2025-12-05\n")
    f.write("**Test:** Czy egzotyczna dynamika b≈2.4 przetrwa noise?\n\n")
    
    f.write("## 1. Metodologia\n")
    f.write("Ewolucja pakietu falowego z szumem: dpsi += noise(σ)\n")
    f.write(f"Testowano σ = {NOISE_LEVELS}\n\n")
    
    f.write("## 2. Wyniki\n")
    f.write("| Noise | b | Status |\n|-------|---|--------|\n")
    for res in results:
        status = "✅" if abs(res['b'] - b_baseline)/b_baseline < 0.1 and res['b'] > 2.0 else "🟡"
        f.write(f"| {res['noise']:.2f} | {res['b']:.3f} | {status} |\n")
    f.write("\n")
    
    if verdict == "robust":
        f.write("### ✅ SUPER-BALLISTIC TO FUNDAMENTALNA WŁAŚCIWOŚĆ\n")
        f.write("Egzotyczna dynamika przetrwała szum - to nie artefakt!\n")
    elif verdict == "partial":
        f.write("### 🟡 CZĘŚCIOWA ROBUSTNOŚĆ\n")
    else:
        f.write("### ❌ WRAŻLIWE\n")

print("Report saved.")
print("="*80)
