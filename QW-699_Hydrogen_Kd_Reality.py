#!/usr/bin/env python3
"""
QW-699: HYDROGEN SPECTRUM FROM K(d) KERNEL - REALITY TEST
==========================================================
Purpose: RIGOROUS TEST - Can K(d) naturally produce 1/n² Rydberg spectrum?

Background:
  - PRIOR FAILURES: QW-221 (250%), QW-505 (no spectrum), QW-V81 (316%)
  - Octaves give LINEAR f ~ n, but Rydberg requires f ~ 1/n²
  - This is the ACID TEST of whether K(d) connects to reality.

Key Question:
  If theory is TRUE, then K(d) should produce 1/n² levels somewhere.
  If theory is FALSE/NUMEROLOGY, K(d) will produce linear/random levels.

Method:
  1. Find eigenvalues of K(d) operator on finite space.
  2. Check if eigenvalues follow 1/n² pattern (Rydberg).
  3. NO FITTING - pure prediction from K(d) parameters.
"""

import numpy as np
from scipy import linalg
import datetime

print("="*80)
print("QW-699: HYDROGEN SPECTRUM FROM K(d) - REALITY TEST")
print("="*80)

# --- K(d) Parameters (FIXED - no fitting) ---
ALPHA = 4 * np.log(2)  # 2.77
OMEGA = np.pi / 4
PHI = np.pi / 6
BETA = 0.1

def K(d):
    """Fixed coupling kernel."""
    if d < 0.1:
        d = 0.1  # Regularize
    return ALPHA * np.cos(OMEGA * d + PHI) / (1 + BETA * d)

# --- Step 1: Build Hamiltonian from K(d) ---
print("\n[1] Building Hamiltonian from K(d)...")

N = 20  # Number of sites (radial points)
H = np.zeros((N, N))

# Hamiltonian: H_ij = K(|i-j|) for coupling, diagonal ~ potential
for i in range(N):
    for j in range(N):
        if i != j:
            d = abs(i - j)
            H[i, j] = -K(d)  # Negative for binding
        else:
            # Diagonal: use K(0) regularized
            H[i, i] = K(0.1)

print(f"  Hamiltonian size: {N}x{N}")
print(f"  H[0,0] = {H[0,0]:.4f}")
print(f"  H[0,1] = {H[0,1]:.4f}")

# --- Step 2: Diagonalize ---
print("\n[2] Diagonalizing...")

eigenvalues, eigenvectors = linalg.eigh(H)
eigenvalues = np.sort(eigenvalues)

# Filter bound states (negative eigenvalues)
bound = eigenvalues[eigenvalues < 0]
print(f"  Total eigenvalues: {len(eigenvalues)}")
print(f"  Bound states (E < 0): {len(bound)}")

if len(bound) == 0:
    print("  ⚠️ No bound states! Theory doesn't predict binding.")
    bound = eigenvalues[:5]  # Take lowest anyway

# --- Step 3: Compare to Hydrogen 1/n² ---
print("\n[3] Comparing to Hydrogen Rydberg Series (1/n²)...")

def rydberg(n):
    """Rydberg formula: E_n = -E_0/n²"""
    return -1 / n**2

# Take first N_compare levels
N_compare = min(len(bound), 5)
E_model = bound[:N_compare]

# Normalize model to ground state
E_model_norm = E_model / abs(E_model[0])

# Rydberg normalized
E_rydberg = [rydberg(n+1) for n in range(N_compare)]
E_rydberg_norm = np.array(E_rydberg) / abs(E_rydberg[0])

print("| n | E_model (norm) | E_Rydberg (norm) | Error |")
print("|---|----------------|------------------|-------|")

errors = []
for n in range(N_compare):
    err = abs(E_model_norm[n] - E_rydberg_norm[n]) / abs(E_rydberg_norm[n]) * 100
    errors.append(err)
    print(f"| {n+1} | {E_model_norm[n]:.4f} | {E_rydberg_norm[n]:.4f} | {err:.1f}% |")

mean_error = np.mean(errors)
print(f"\nMean Error: {mean_error:.1f}%")

# --- Step 4: Check if levels are closer to LINEAR n ---
print("\n[4] Comparison: Is spectrum closer to LINEAR or 1/n²?")

E_linear = [(n+1) / 1 for n in range(N_compare)]  # Linear
E_linear_norm = np.array(E_linear) / E_linear[0]
# Match sign
E_linear_norm = -E_linear_norm / E_linear_norm[-1] * E_model_norm[-1]

linear_errors = []
for n in range(N_compare):
    err = abs(E_model_norm[n] - E_linear_norm[n]) / abs(E_linear_norm[n] + 1e-10) * 100
    linear_errors.append(err)

mean_linear_error = np.mean(linear_errors)

print(f"Mean Error vs 1/n² (Rydberg): {mean_error:.1f}%")
print(f"Mean Error vs Linear n: {mean_linear_error:.1f}%")

if mean_error < mean_linear_error and mean_error < 50:
    verdict = "✅ CLOSER TO RYDBERG"
elif mean_error < 100:
    verdict = "🟡 PARTIALLY RYDBERG-LIKE"
else:
    verdict = "❌ NOT RYDBERG (LINEAR OR RANDOM)"

# --- Step 5: Emergent Observer Check ---
print("\n[5] Emergent Observer Perspective...")
print("""
In Emergent Observer framework:
- The 1/n² pattern emerges from INTERNAL perspective
- K(d) gives EXTERNAL coupling structure
- Transformation between perspectives might give 1/n²

If K(d) → H → 1/n²: Theory connects to reality.
If K(d) → H → linear: Theory is disconnected from atomic physics.
""")

# --- Verdict ---
print("\n" + "="*80)
print("VERDICT")
print("="*80)

print(f"\n{verdict}")
print(f"Mean Rydberg Error: {mean_error:.1f}%")
print(f"Prior tests (QW-221/V81): 250-316% error")
print(f"Improvement: {'YES' if mean_error < 200 else 'NO'}")

if mean_error < 50:
    conclusion = "K(d) PRODUCES hydrogen-like spectrum! Theory connects to reality."
elif mean_error < 150:
    conclusion = "K(d) produces PARTIAL hydrogen structure. Needs refinement."
else:
    conclusion = "K(d) FAILS to produce hydrogen spectrum. Theory disconnected from atomic physics."

print(f"\n{conclusion}")

# --- Report ---
report_file = "raport_qw699_hydrogen_Kd_reality.md"
print(f"\nSaving report to {report_file}...")

with open(report_file, "w") as f:
    f.write("# RAPORT QW-699: HYDROGEN FROM K(d) - REALITY TEST\n")
    f.write(f"**Date:** {datetime.datetime.now()}\n\n")
    
    f.write("## 1. Question\n")
    f.write("Can K(d) naturally produce hydrogen-like 1/n² spectrum?\n\n")
    
    f.write("## 2. Prior Results\n")
    f.write("| Test | Error |\n")
    f.write("|------|-------|\n")
    f.write("| QW-221 | 250% |\n")
    f.write("| QW-505 | No spectrum |\n")
    f.write("| QW-V81 | 316% |\n\n")
    
    f.write("## 3. QW-699 Results\n")
    f.write("| n | E_model | E_Rydberg | Error |\n")
    f.write("|---|---------|-----------|-------|\n")
    for n in range(N_compare):
        f.write(f"| {n+1} | {E_model_norm[n]:.4f} | {E_rydberg_norm[n]:.4f} | {errors[n]:.1f}% |\n")
    f.write(f"\n**Mean Error:** {mean_error:.1f}%\n\n")
    
    f.write("## 4. Verdict\n")
    f.write(f"### {verdict}\n")
    f.write(f"{conclusion}\n\n")
    
    if mean_error > 100:
        f.write("## 5. Implications\n")
        f.write("The K(d) kernel, in its current form, does NOT produce atomic spectra.\n")
        f.write("This is a FUNDAMENTAL DISCONNECT between FIN Theory and atomic physics.\n")
        f.write("Options:\n")
        f.write("1. K(d) operates at different scale than atomic\n")
        f.write("2. Additional physics needed (Coulomb potential?)\n")
        f.write("3. Theory is incomplete for atomic regime\n")

print("Done.")
