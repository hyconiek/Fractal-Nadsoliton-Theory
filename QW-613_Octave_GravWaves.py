#!/usr/bin/env python3
# QW-613: OCTAVE GRAVITY WAVE PROPAGATION (Proper Version)
# Purpose: Test if mass perturbations propagate as waves in octave network
# Replaces QW-607 which used spatial Hebbian network incorrectly
# Proper approach: Mass on specific octaves, measure connectivity wave
# Date: 2025-12-05

import numpy as np

print("="*80)
print("QW-613: OCTAVE GRAVITY WAVE PROPAGATION (Proper Octave-Based)")
print("="*80)
print("Test: Do mass perturbations create waves in octave network?")
print("Hypothesis H9: Gravity = Hebbian strengthening")
print("="*80)

# ============================================================================
# NETWORK PARAMETERS
# ============================================================================
N_OCTAVES = 12
ALPHA_GEO = 4 * np.log(2)
OMEGA = np.pi / 4
PHI = np.pi / 6
BETA_TORS = 0.01

ALPHA_HEBBIAN = 0.05 
DT = 0.01
STEPS = 300

print(f"\nNetwork: {N_OCTAVES} octaves")
print(f"Hebbian α: {ALPHA_HEBBIAN}")
print(f"Timesteps: {STEPS}")
print("-" * 40)

# ============================================================================
# INITIAL CONNECTIVITY
# ============================================================================

def K(d):
    return ALPHA_GEO * np.cos(OMEGA * d + PHI) / (1 + BETA_TORS * d)

K_initial = np.zeros((N_OCTAVES, N_OCTAVES))
for i in range(N_OCTAVES):
    for j in range(N_OCTAVES):
        K_initial[i, j] = K(abs(i - j))

print("Initial coupling matrix ready")

# ============================================================================
# MASS PERTURBATION
# ============================================================================

# Start with uniform background mass
mass = np.ones(N_OCTAVES) * 0.1

# Add sudden mass perturbation at octave 6 (middle)
perturb_octave = 6
mass_increment = 2.0  # 20× normal

print(f"\nAdding mass perturbation:")
print(f"  Octave: {perturb_octave}")
print(f"  Δm: {mass_increment}")
print(f"  t_perturb: 0")

t_perturb = 0  # Instant at t=0

# ============================================================================
# HEBBIAN EVOLUTION
# ============================================================================

print("\nSimulating Hebbian dynamics...")

K_matrix = K_initial.copy()

# Track connectivity at different octave distances
track_octaves = [2, 4, 6, 8, 10]  # Distances from perturbation
K_history = {oct: [] for oct in track_octaves}
time_history = []

# Track peak times
peak_times = {}

for t in range(STEPS):
    # At t=0, add mass perturbation
    if t == t_perturb:
        mass[perturb_octave] += mass_increment
    
    # Hebbian update: ΔK ∝ m_i × m_j × K
    mass_outer = np.outer(mass, mass)
    dK = ALPHA_HEBBIAN * mass_outer * K_matrix * DT
    K_matrix += dK
    
    # Track relevant connections
    if t % 5 == 0:
        for oct in track_octaves:
            if oct < N_OCTAVES:
                # Connection strength from perturb_octave to oct
                K_history[oct].append(K_matrix[perturb_octave, oct])
        
        time_history.append(t * DT)
    
    if t % 75 == 0:
        K_sample = K_matrix[perturb_octave, 8] if 8 < N_OCTAVES else 0
        print(f"t={t:3d}: K[6→8] = {K_sample:.4f}")

print("-" * 40)

# ============================================================================
# WAVE ANALYSIS
# ============================================================================

print("\nAnalyzing wave propagation...")

results = {}

for oct in track_octaves:
    if len(K_history[oct]) == 0:
        continue
    
    K_values = np.array(K_history[oct])
    
    # Initial value (before perturbation)
    K_initial_val = K_values[0]
    
    # Find peak (maximum change)
    peak_idx = np.argmax(np.abs(K_values - K_initial_val))
    peak_time = time_history[peak_idx] if peak_idx < len(time_history) else 0
    peak_value = K_values[peak_idx]
    
    # Change from initial
    delta_K = peak_value - K_initial_val
    
    # Distance from perturbation
    distance = abs(oct - perturb_octave)
    
    results[oct] = {
        'distance': distance,
        'peak_time': peak_time,
        'delta_K': delta_K,
        'relative_change': (delta_K / K_initial_val * 100) if K_initial_val != 0 else 0
    }

print("\n| Octave | Distance | Peak Time | ΔK | Change % |")
print("|--------|----------|-----------|-----|----------|")
for oct in sorted(results.keys()):
    r = results[oct]
    print(f"| {oct:6d} | {r['distance']:8d} | {r['peak_time']:9.2f} | {r['delta_K']:+.3f} | {r['relative_change']:7.1f}% |")

# ============================================================================
# TEST: WAVE VELOCITY
# ============================================================================

print("\n" + "="*80)
print("TEST: Is there wave propagation? (t_peak ∝ distance)")
print("="*80)

# Extract distances and peak times
distances = np.array([results[oct]['distance'] for oct in results.keys() if results[oct]['distance'] > 0])
peak_times = np.array([results[oct]['peak_time'] for oct in results.keys() if results[oct]['distance'] > 0])

if len(distances) > 2:
    # Linear fit: t = d/v + t0
    coeffs = np.polyfit(distances, peak_times, 1)
    v_wave = 1.0 / coeffs[0] if coeffs[0] > 0 else np.inf
    t0 = coeffs[1]
    
    # R²
    t_mean = np.mean(peak_times)
    ss_tot = np.sum((peak_times - t_mean)**2)
    fit_vals = np.polyval(coeffs, distances)
    ss_res = np.sum((peak_times - fit_vals)**2)
    r_squared = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0
    
    print(f"\nLinear fit: t = d/{v_wave:.2f} + {t0:.2f}")
    print(f"Wave velocity: v = {v_wave:.3f} octaves/time")
    print(f"R² (linearity): {r_squared:.3f}")
    
    if r_squared > 0.8 and v_wave < np.inf:
        print("\n✅ WAVE PROPAGATION DETECTED!")
        print(f"   Mode travels at v={v_wave:.2f} octaves per unit time")
        print("   Connectivity changes propagate through network!")
        wave_detected = True
    elif r_squared > 0.5:
        print("\n🟡 WEAK WAVE PATTERN")
        print(f"   Some propagation (R²={r_squared:.2f}) but noisy")
        wave_detected = False
    else:
        print("\n❌ NO WAVE PROPAGATION")
        print(f"   Changes instantaneous or random (R²={r_squared:.2f})")
        wave_detected = False
else:
    wave_detected = False
    v_wave = 0
    r_squared = 0
    print("\nInsufficient data points")

# ============================================================================
# REPORT
# ============================================================================

report_path = "raport_qw613_octave_waves.md"
print(f"\nGenerating report: {report_path}...")

with open(report_path, "w") as f:
    f.write("# Raport QW-613: Octave Gravity Wave Propagation\n")
    f.write("**Data:** 2025-12-05\n")
    f.write("**Hipoteza:** H9 - Grawitacja = Uczenie Hebba (wave test)\n\n")
    
    f.write("## 1. Setup\n")
    f.write(f"Sieć: {N_OCTAVES} oktaw\n")
    f.write(f"Perturbacja masy: {mass_increment}× na oktawie {perturb_octave}\n")
    f.write(f"Hebbian α: {ALPHA_HEBBIAN}\n\n")
    
    f.write("## 2. Propagacja\n")
    f.write("| Oktawa | Odległość | Czas Szczytu | ΔK |\n")
    f.write("|--------|-----------|--------------|----|\n")
    for oct in sorted(results.keys()):
        r = results[oct]
        f.write(f"| {oct} | {r['distance']} | {r['peak_time']:.2f} | {r['delta_K']:+.3f} |\n")
    f.write("\n")
    
    f.write("## 3. Analiza Fali\n")
    f.write(f"**Prędkość:** v = {v_wave:.3f}\n")
    f.write(f"**Liniowość:** R² = {r_squared:.3f}\n\n")
    
    if wave_detected:
        f.write("### ✅ FALA GRAWITACYJNA W OKTAWACH!\n")
        f.write("Zmiana masy propaguje jako fala przez network oktaw.\n")
        f.write(f"Prędkość: v = {v_wave:.2f} oktaw/czas\n")
    else:
        f.write("### 🟡 SŁABY LUB BRAK SYGNAŁU\n")
        f.write("Propagacja nie jest wyraźnie falowa. Możliwe że:\n")
        f.write("- Coupling zbyt silny (instant)\n")
        f.write("- Wymaga większej sieci (N>12)\n")

print("Report saved.")
print("="*80)
