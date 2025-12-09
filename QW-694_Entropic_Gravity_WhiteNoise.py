#!/usr/bin/env python3
"""
QW-694: ENTROPIC GRAVITY FROM WHITE NOISE VACUUM
=================================================
Purpose: Test if gravitational attraction emerges from PURE WHITE NOISE
         via the Verlinde mechanism: F = T * dS/dx

Critical Question: QW-693 showed the vacuum is White Noise (ζ ≈ 0).
                   Can forces still emerge entropically from such chaos?

Method:
  1. Generate a "Thermal Vacuum" (pure white noise field).
  2. Insert two "topological defects" (mass sources) at distance R.
  3. Measure the ENTROPY GRADIENT dS/dx around each defect.
  4. Calculate entropic force F = T * dS/dx.
  5. Check if F ∝ 1/R² (Newton) or F ∝ 1/R (entropic surface law).

Expected:
  - If vacuum is structureless noise, dS/dx should be ZERO everywhere.
  - BUT: The *presence* of a defect (mass) may CREATE local order,
         which generates an entropy gradient (shadow entropy).
"""

import numpy as np
import datetime

print("="*80)
print("QW-694: ENTROPIC GRAVITY FROM WHITE NOISE VACUUM")
print("="*80)

# --- Parameters ---
N = 256          # Grid size (3D)
T_VACUUM = 1.0   # Effective temperature (arbitrary units)
N_SAMPLES = 50   # Number of noise realizations to average over
DEFECT_RADIUS = 3  # Radius of topological defect (mass)

np.random.seed(694)

# --- Step 1: Define the Thermal Vacuum (White Noise) ---
print("\n[1] Generating White Noise Vacuum...")

def generate_white_noise(size):
    """Pure white noise - no correlations."""
    return np.random.randn(size, size, size)

# --- Step 2: Insert Topological Defects ---
print("[2] Inserting Topological Defects (Masses)...")

def insert_defect(field, center, radius, amplitude=1.0):
    """
    Insert a topological defect (local order) into the noise field.
    The defect is a Gaussian "bump" of reduced entropy.
    """
    x = np.arange(field.shape[0])
    y = np.arange(field.shape[1])
    z = np.arange(field.shape[2])
    X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
    
    dist_sq = (X - center[0])**2 + (Y - center[1])**2 + (Z - center[2])**2
    
    # Defect creates LOCAL ORDER (low entropy region)
    # by replacing noise with a smooth function
    defect_profile = amplitude * np.exp(-dist_sq / (2 * radius**2))
    
    # The defect "replaces" noise with order
    field_with_defect = field.copy()
    field_with_defect += defect_profile  # Shift mean locally
    
    return field_with_defect

# --- Step 3: Measure Entropy ---
print("[3] Measuring Local Entropy...")

def local_entropy(field, center, measure_radius=5):
    """
    Measure the local variance (proxy for entropy) around a point.
    For white noise, variance = 1. For ordered regions, variance < 1.
    """
    x0, y0, z0 = center
    
    # Extract local cube
    x_min = max(0, x0 - measure_radius)
    x_max = min(field.shape[0], x0 + measure_radius)
    y_min = max(0, y0 - measure_radius)
    y_max = min(field.shape[1], y0 + measure_radius)
    z_min = max(0, z0 - measure_radius)
    z_max = min(field.shape[2], z0 + measure_radius)
    
    local_region = field[x_min:x_max, y_min:y_max, z_min:z_max]
    
    # Entropy ≈ Variance for Gaussian fields
    # (More rigorous: Shannon entropy of histogram, but variance is faster)
    return np.var(local_region)

def entropy_gradient(field, center, direction, delta=2, measure_radius=5):
    """
    Measure dS/dx along a given direction.
    """
    dir_vec = np.array(direction) / np.linalg.norm(direction)
    
    # Sample entropy at center - delta and center + delta
    pos_plus = (center + delta * dir_vec).astype(int)
    pos_minus = (center - delta * dir_vec).astype(int)
    
    # Clamp to grid
    pos_plus = np.clip(pos_plus, measure_radius, N - measure_radius - 1)
    pos_minus = np.clip(pos_minus, measure_radius, N - measure_radius - 1)
    
    S_plus = local_entropy(field, tuple(pos_plus), measure_radius)
    S_minus = local_entropy(field, tuple(pos_minus), measure_radius)
    
    dS_dx = (S_plus - S_minus) / (2 * delta)
    return dS_dx

# --- Step 4: Entropic Force Calculation ---
print("[4] Calculating Entropic Force vs Distance...")

# Place defect 1 at center
center1 = np.array([N//2, N//2, N//2])

# Scan distances
distances = np.arange(20, N//3, 10)
forces = []
forces_std = []

for R in distances:
    center2 = center1 + np.array([R, 0, 0])  # Second defect along x-axis
    
    # Average over multiple noise realizations
    force_samples = []
    
    for _ in range(N_SAMPLES):
        # Generate fresh white noise
        vacuum = generate_white_noise(N)
        
        # Insert both defects
        field = insert_defect(vacuum, center1, DEFECT_RADIUS)
        field = insert_defect(field, center2, DEFECT_RADIUS)
        
        # Measure entropy gradient at the midpoint, pointing toward defect 2
        midpoint = ((center1 + center2) / 2).astype(int)
        direction = center2 - center1
        
        dS_dx = entropy_gradient(field, midpoint, direction)
        
        # Entropic force: F = T * dS/dx
        F = T_VACUUM * dS_dx
        force_samples.append(F)
    
    mean_F = np.mean(force_samples)
    std_F = np.std(force_samples)
    forces.append(mean_F)
    forces_std.append(std_F)
    
    print(f"  R = {R:3d}  |  F = {mean_F:+.6f} ± {std_F:.6f}")

forces = np.array(forces)
forces_std = np.array(forces_std)

# --- Step 5: Fit Power Law ---
print("\n[5] Fitting Power Law: F ∝ R^n")

# Only use points where force is significantly non-zero
valid = np.abs(forces) > 1e-6

if np.sum(valid) > 3:
    log_R = np.log(distances[valid])
    log_F = np.log(np.abs(forces[valid]) + 1e-10)
    
    coeffs = np.polyfit(log_R, log_F, 1)
    n_power = coeffs[0]
    
    print(f"Power law exponent: n = {n_power:.4f}")
    print("Expected for Newton: n = -2")
    print("Expected for Verlinde (surface): n = -1")
else:
    n_power = 0
    print("⚠️ Force too weak to fit power law!")

# --- Step 6: Verdict ---
print("\n" + "="*80)
print("VERDICT")
print("="*80)

mean_force_magnitude = np.mean(np.abs(forces))
noise_level = np.mean(forces_std)

# Is there a signal above noise?
signal_to_noise = mean_force_magnitude / (noise_level + 1e-10)

print(f"Mean |F|: {mean_force_magnitude:.6f}")
print(f"Noise level: {noise_level:.6f}")
print(f"Signal-to-Noise: {signal_to_noise:.2f}")

if signal_to_noise < 1.5:
    verdict = "NOISE_DOMINATED"
    verdict_text = "❌ NO ENTROPIC FORCE DETECTED"
    explanation = "The 'force' is indistinguishable from random noise. White Noise vacuum does NOT generate gravity."
elif n_power < -1.5 and n_power > -2.5:
    verdict = "NEWTON_LIKE"
    verdict_text = "✅ NEWTONIAN GRAVITY EMERGES"
    explanation = f"Entropic force scales as F ~ 1/R^{abs(n_power):.1f}, consistent with Newton's Law!"
elif n_power < -0.5 and n_power > -1.5:
    verdict = "VERLINDE_LIKE"
    verdict_text = "✅ VERLINDE GRAVITY EMERGES"
    explanation = f"Entropic force scales as F ~ 1/R^{abs(n_power):.1f}, consistent with surface-law (Verlinde)."
else:
    verdict = "UNKNOWN"
    verdict_text = f"🟡 INCONCLUSIVE (n = {n_power:.2f})"
    explanation = "Power law does not match known gravity models."

print(f"\n{verdict_text}")
print(f"Explanation: {explanation}")

# --- Report ---
report_file = "raport_qw694_entropic_gravity_whitenoise.md"
print(f"\nSaving report to {report_file}...")

with open(report_file, "w") as f:
    f.write("# RAPORT QW-694: ENTROPIC GRAVITY FROM WHITE NOISE\n")
    f.write(f"**Date:** {datetime.datetime.now()}\n\n")
    
    f.write("## 1. Question\n")
    f.write("Can gravitational attraction emerge from a **pure White Noise vacuum** (QW-693)?\n\n")
    
    f.write("## 2. Method\n")
    f.write("- Generated 3D white noise field (N=256).\n")
    f.write("- Inserted two 'topological defects' (mass sources) at varying distances R.\n")
    f.write("- Measured entropy gradient dS/dx between defects.\n")
    f.write(f"- Calculated entropic force F = T * dS/dx (averaged over {N_SAMPLES} realizations).\n\n")
    
    f.write("## 3. Results\n")
    f.write("| R (distance) | F (entropic force) | σ (noise) |\n")
    f.write("|--------------|--------------------|-----------|\n")
    for i, R in enumerate(distances):
        f.write(f"| {R} | {forces[i]:+.6f} | {forces_std[i]:.6f} |\n")
    f.write("\n")
    
    f.write(f"**Power Law Fit:** F ∝ R^{n_power:.4f}\n\n")
    f.write(f"**Signal-to-Noise:** {signal_to_noise:.2f}\n\n")
    
    f.write("## 4. Verdict\n")
    f.write(f"### {verdict_text}\n")
    f.write(f"{explanation}\n\n")
    
    if verdict == "NOISE_DOMINATED":
        f.write("### ⚠️ IMPLICATIONS FOR FIN THEORY\n")
        f.write("If the vacuum is truly White Noise (QW-693), then:\n")
        f.write("1. **Gravity cannot be purely entropic** from the vacuum itself.\n")
        f.write("2. **Defects (masses) must create local order** that then generates entropy gradients.\n")
        f.write("3. Forces may emerge from **topological interactions**, not vacuum entropy.\n")
    elif "EMERGES" in verdict_text:
        f.write("### ✅ IMPLICATIONS FOR FIN THEORY\n")
        f.write("The entropic mechanism SURVIVES even in White Noise vacuum:\n")
        f.write("- Topological defects create 'entropy shadows'.\n")
        f.write("- The gradient between shadows generates attractive force.\n")
        f.write("- This is consistent with H6 (Forces = Entropic Gradients).\n")

print("Done.")
