#!/usr/bin/env python3
# QW-607: GRAVITATIONAL WAVE PROPAGATION
# Purpose: Test if mass perturbations propagate as waves through Hebbian network
# Hypothesis H9: Gravity = Hebbian Learning → mass creates connectivity waves
# Date: 2025-12-05

import numpy as np

print("="*80)
print("QW-607: GRAVITATIONAL WAVE PROPAGATION")
print("="*80)
print("Test: Do mass perturbations create propagating waves in network?")
print("Hypothesis H9: Gravity = Hebbian strengthening")
print("="*80)

# ============================================================================
# NETWORK PARAMETERS
# ============================================================================
N_NODES = 100
DT = 0.01
STEPS = 500

# Hebbian learning rate
ALPHA_HEBBIAN = 0.05

print(f"Network size: {N_NODES} nodes")
print(f"Hebbian rate: {ALPHA_HEBBIAN}")
print("-" * 40)

# ============================================================================
# INITIALIZE NETWORK
# ============================================================================

# Positions in 1D (for simplicity)
positions = np.linspace(0, 10, N_NODES)

# Initial connectivity (distance-based)
K = np.zeros((N_NODES, N_NODES))
for i in range(N_NODES):
    for j in range(N_NODES):
        if i != j:
            d = abs(positions[i] - positions[j])
            K[i, j] = np.exp(-d / 2.0)  # Baseline connectivity

# Mass at each node (initially uniform)
mass = np.ones(N_NODES)

print("Initialized network with baseline connectivity K(d)~exp(-d/2)")

# ============================================================================
# MASS PERTURBATION
# ============================================================================

# Add sudden mass at center
center_idx = N_NODES // 2
mass_perturbation = 5.0  # 5× normal

print(f"\nAdding mass perturbation at node {center_idx}")
print(f"Mass increase: {mass_perturbation}×")

mass[center_idx] += mass_perturbation

# Note time of perturbation
t_perturb = 0

# ============================================================================
# HEBBIAN EVOLUTION
# ============================================================================

print("\nSimulating Hebbian dynamics...")

# Track connectivity at different distances from perturbation
distances_to_track = [1, 5, 10, 20, 30]
K_history = {d: [] for d in distances_to_track}
time_history = []

for t in range(STEPS):
    # Hebbian rule: ΔK_ij ∝ m_i × m_j × <activity correlation>
    # Simplified: ΔK_ij ∝ m_i × m_j × K_ij (reinforcement)
    
    # Compute mass outer product
    mass_matrix = np.outer(mass, mass)
    
    # Hebbian update
    dK = ALPHA_HEBBIAN * mass_matrix * K * DT
    K += dK
    
    # Track K values at specific distances from perturbation
    if t % 10 == 0:
        for dist in distances_to_track:
            target_idx = center_idx + dist
            if 0 <= target_idx < N_NODES:
                K_history[dist].append(K[center_idx, target_idx])
            else:
                K_history[dist].append(0)
        
        time_history.append(t * DT)
        
        if t % 100 == 0:
            print(f"t={t:4d}: K[center, +10] = {K[center_idx, min(center_idx+10, N_NODES-1)]:.4f}")

print("-" * 40)

# ============================================================================
# ANALYSIS: WAVE DETECTION
# ============================================================================

print("\nAnalyzing wave propagation...")

# Check if K increases propagate outward (wave-like)
# Expected: K at distance d peaks at time t ~ d/v_wave

results = {}

for dist in distances_to_track:
    if len(K_history[dist]) > 0:
        K_values = np.array(K_history[dist])
        
        # Find peak time
        peak_idx = np.argmax(K_values)
        peak_time = time_history[peak_idx] if peak_idx < len(time_history) else 0
        peak_value = K_values[peak_idx]
        
        # Initial value
        initial_value = K_values[0]
        
        # Amplification
        amplification = (peak_value - initial_value) / initial_value if initial_value > 0 else 0
        
        results[dist] = {
            'peak_time': peak_time,
            'peak_value': peak_value,
            'amplification': amplification * 100
        }

print("\n| Distance | Peak Time | Amplification |")
print("|----------|-----------|---------------|")
for dist in sorted(results.keys()):
    r = results[dist]
    print(f"| {dist:8d} | {r['peak_time']:9.2f} | {r['amplification']:12.1f}% |")

# Check if peak_time ∝ distance (wave behavior)
distances = np.array(list(results.keys()))
peak_times = np.array([results[d]['peak_time'] for d in distances])

# Filter out zero times
valid = peak_times > 0.1
if np.sum(valid) > 2:
    distances_valid = distances[valid]
    times_valid = peak_times[valid]
    
    # Linear fit: t = d/v
    coeffs = np.polyfit(distances_valid, times_valid, 1)
    v_wave = 1.0 / coeffs[0] if coeffs[0] > 0 else 0
    
    # R²
    t_mean = np.mean(times_valid)
    ss_tot = np.sum((times_valid - t_mean)**2)
    ss_res = np.sum((times_valid - (coeffs[0] * distances_valid + coeffs[1]))**2)
    r_squared = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0
    
    print(f"\nWave velocity: v = {v_wave:.3f} units/time")
    print(f"Linearity (R²): {r_squared:.4f}")
    
    if r_squared > 0.8:
        print("\n✅ GRAVITATIONAL WAVE DETECTED!")
        print(f"   Mass perturbation propagates as wave with v={v_wave:.2f}")
        print("   Peak time ∝ distance (R²={:.3f})".format(r_squared))
        wave_detected = True
    else:
        print("\n🟡 WEAK WAVE PATTERN")
        print(f"   Propagation unclear (R²={r_squared:.3f})")
        wave_detected = False
else:
    wave_detected = False
    v_wave = 0
    r_squared = 0

# ============================================================================
# REPORT
# ============================================================================

report_path = "raport_qw607_grav_waves.md"
print(f"\nGenerating report: {report_path}...")

with open(report_path, "w") as f:
    f.write("# Raport QW-607: Gravitational Wave Propagation\n")
    f.write("**Data:** 2025-12-05\n")
    f.write("**Hipoteza:** H9 - Grawitacja = Uczenie Hebba\n\n")
    
    f.write("## 1. Setup\n")
    f.write(f"- Sieć: {N_NODES} węzłów (1D)\n")
    f.write(f"- Perturbacja masy: {mass_perturbation}× w centrum\n")
    f.write(f"- Hebbian learning: ΔK ∝ m_i × m_j × K_ij\n\n")
    
    f.write("## 2. Propagacja\n")
    f.write("| Odległość | Czas Szczytu | Wzmocnienie |\n")
    f.write("|-----------|--------------|-------------|\n")
    for dist in sorted(results.keys()):
        r = results[dist]
        f.write(f"| {dist} | {r['peak_time']:.2f} | {r['amplification']:.1f}% |\n")
    f.write("\n")
    
    f.write("## 3. Analiza Fali\n")
    f.write(f"**Prędkość fali:** v = {v_wave:.3f}\n")
    f.write(f"**Liniowość (R²):** {r_squared:.4f}\n\n")
    
    if wave_detected:
        f.write("### ✅ FALA GRAWITACYJNA POTWIERDZONA!\n\n")
        f.write("Perturbacja masy propaguje jako fala przez sieć Hebbian.\n")
        f.write(f"Czas szczytgęści ∝ odległość (R²={r_squared:.3f})\n\n")
        f.write("**Implikacja dla H9:**\n")
        f.write("Grawitacja nie jest statyczną siłą,ale **dynamicznym procesem**\n")
        f.write("wzmacniania połączeń. Zmiany masy tworzą fale w sieci!\n")
    else:
        f.write("### 🟡 SŁABY SYGNAŁ FALI\n")
        f.write("Propagacja nie jest wyraźnie liniowa. Możliwe że:\n")
        f.write("- Sieć za mała\n")
        f.write("- Hebbian rate za mały/duży\n")
        f.write("- Wymaga 3D topologii\n")

print("Report saved.")
print("="*80)
