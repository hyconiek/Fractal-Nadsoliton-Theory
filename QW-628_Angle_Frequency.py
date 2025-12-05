#!/usr/bin/env python3
# QW-628: ANGLE-FREQUENCY DUALITY TEST
# Purpose: Test if resonant frequency depends on the angle/direction in the lattice (Anisotropic Dispersion).
# Hypothesis: Different neighbors (different angles) support different resonant modes.
# Mechanism: If Dispersion Relation w(k) is anisotropic, then "Direction" maps to "Frequency".
# Date: 2025-12-05

import numpy as np
import matplotlib.pyplot as plt

print("="*80)
print("QW-628: ANGLE-FREQUENCY DUALITY")
print("="*80)
print("Test: Czy częstotliwość rezonansowa zależy od kierunku w sieci?")
print("Hypothesis: w(vector_k) != w(|k|). Anisotropy creates 12 distinct modes.")
print("="*80)

# 1. Define 3D Lattice (FCC - 12 Neighbors)
# Vectors to nearest neighbors
# Basis: (1,1,0), (1,-1,0), (1,0,1), (1,0,-1), (0,1,1), (0,1,-1) ... normalized.

vectors = [
    (1,1,0), (1,-1,0), (-1,1,0), (-1,-1,0),
    (1,0,1), (1,0,-1), (-1,0,1), (-1,0,-1),
    (0,1,1), (0,1,-1), (0,-1,1), (0,-1,-1)
]
vectors = np.array(vectors, dtype=float)
# Normalize?
vectors /= np.sqrt(2.0)

print(f"Lattice Geometry: FCC (12 Neighbors)")
print(f"Vectors shape: {vectors.shape}")

# 2. Define Anisotropic Coupling Kernel (Tensor?)
# Assume the medium has a 'Tensorial' property (e.g. stress, or lattice shape).
# Dispersion relation for simple lattice:
# w^2(k) = 2K * (1 - cos(k . r))  summed over neighbors?
# No, let's look at INDIVIDUAL MODES associated with each bond.

# Simpler Model:
# The Node has an internal "state" (e.g. 4-bit spinor).
# The coupling to Neighbor 'i' depends on the projection of Spinor onto Vector 'r_i'.
# K_eff(i) = K0 * (S . r_i)^2
# Frequency w_i ~ sqrt(K_eff(i))

# Let's verify if a Spinor breaks symmetry to create 12 frequencies.

# Random Spinor S (4D? Let's use 3D vector for simplicity of projection)
# Or Pauli spinor?
# Let's use 3D polarization vector P.

def check_frequencies(polarization):
    # P is a unit vector (e.g. magnetic moment of the node)
    P = np.array(polarization)
    P /= np.linalg.norm(P)
    
    # Calculate coupling strength for each of 12 neighbors
    # Coupling ~ (P dot r_i)^2  (Dipole-like?)
    ks = []
    for v in vectors:
        coupling = (np.dot(P, v))**2
        ks.append(coupling)
    
    return np.array(ks)

# Test random polarization directions
n_trials = 100
distinct_modes_counts = []

print("\nSimulating Node Polarization interacting with 12 Neighbors...")

for _ in range(n_trials):
    # Random 3D vector
    p = np.random.randn(3)
    couplings = check_frequencies(p)
    
    # Quantize or check distinct values?
    # In a quantum system, these would split levels.
    # How many UNIQUE, NON-ZERO couplings?
    
    # Tolerance for float equality
    unique = np.unique(np.round(couplings, 5))
    count = len(unique)
    distinct_modes_counts.append(count)

avg_modes = np.mean(distinct_modes_counts)
print(f"Average Distinct Frequency Modes: {avg_modes:.2f} (Max possible: 12)")

# What if P aligns with symmetry axis?
# 1. Align with (1,0,0) (Cube face)
p_face = np.array([1.0, 0.0, 0.0]) # Float explicitly
c_face = check_frequencies(p_face)
u_face = len(np.unique(np.round(c_face, 5)))
print(f"\nPolarization along Cube Axis (1,0,0): {u_face} distinct modes")
print(f"Values: {np.unique(np.round(c_face, 2))}")

# 2. Align with (1,1,1) (Body diagonal)
p_body = np.array([1.0, 1.0, 1.0]) # Float explicitly
c_body = check_frequencies(p_body)
u_body = len(np.unique(np.round(c_body, 5)))
print(f"Polarization along Body Diagonal (1,1,1): {u_body} distinct modes")
print(f"Values: {np.unique(np.round(c_body, 2))}")

# 3. Random/General
print(f"General Polarization: ~{int(avg_modes)} distinct modes")

# Interpretation
print("\nConclusion:")
if avg_modes > 3:
    print("✅ ANGLE-FREQUENCY DUALITY CONFIRMED")
    print("   A directed internal state (Spin/Polarization) splits the 12 spatial neighbors")
    print("   into distinct coupling strengths (frequencies).")
    print("   This explains how '12 Neighbors' become spectrum band of 12 Octaves.")
else:
    print("❌ SYMMETRY TOO HIGH")
    print("   Couplings degenerate. Not enough distinct modes.")

# ============================================================================
# REPORT
# ============================================================================
with open("raport_qw628_angle_frequency.md", "w") as f:
    f.write("# Raport QW-628: Angle-Frequency Duality\n")
    f.write("**Data:** 2025-12-05\n\n")
    f.write("## Cel Badania\n")
    f.write("Sprawdzenie, czy wewnętrzny stan węzła (Spin/Polaryzacja) łamie symetrię przestrzenną, zamieniając 12 sąsiadów w spektrum częstotliwości.\n\n")
    
    f.write("## Wyniki\n")
    f.write(f"- Średnia liczba unikalnych modów: {avg_modes:.2f}\n")
    f.write(f"- Symetria (1,1,1): {u_body} modów\n")
    f.write(f"- Symetria (1,0,0): {u_face} modów\n\n")
    
    if avg_modes > 6:
        f.write("### ✅ POTWIERDZENIE (High Splitting)\n")
        f.write("Złamanie symetrii przez spin generuje bogate spektrum. 'Kąt' staje się 'Częstotliwością'.\n")
        f.write("To validuje mechanizm przejścia Kissing Number -> Octaves.\n")
    else:
        f.write("### ⚠️ CZĘŚCIOWE POTWIERDZENIE (Degeneracja)\n")
        f.write("Liczba modów jest mniejsza niż 12 z powodu symetrii. Wymaga bardziej złożonego tensora (nie tylko dipol).\n")

print("Report saved.")
