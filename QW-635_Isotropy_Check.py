#!/usr/bin/env python3
# QW-635: ISOTROPY CHECK (DESTRUCTION PROTOCOL)
# Purpose: Test for Lorentz Violation (Anisotropy of c) on FCC Lattice.
#          Commanded by Prof. Scepticus.
# Method: Compute Band Structure w(k) and Group Velocity v_g(k).
#         Compare v_g along (1,0,0) and (1,1,1).
# System: FCC Lattice with FIN Coupling K(d).
# Date: 2025-12-05

import numpy as np
import matplotlib.pyplot as plt

print("="*80)
print("QW-635: ISOTROPY CHECK (DESTRUCTION PROTOCOL)")
print("="*80)
print("Test: Is the Speed of Light 'c' constant in all directions?")
print("System: FCC Lattice, Long-Range Coupling K(d).")
print("="*80)

# 1. Define FCC Lattice Neighbors (12 Nearest)
# Basis vectors for neighbors
# d = (+-1, +-1, 0) and permutations.
vectors = [
    (1,1,0), (1,-1,0), (-1,1,0), (-1,-1,0),
    (1,0,1), (1,0,-1), (-1,0,1), (-1,0,-1),
    (0,1,1), (0,1,-1), (0,-1,1), (0,-1,-1)
]
vectors = np.array(vectors, dtype=float)

# 2. Define Coupling K(d)
# Standard hopping: t = 1 for nearest neighbors.
# FIN Theory: Long range? Let's test NEAREST first (Standard Lattice).
# If Standard fails, we add Long Range K(d) ~ 1/d^alpha to see if it fixes it.

def get_dispersion(kx, ky, kz, long_range=False):
    # E(k) = Sum_d t_d * exp(i k.d)
    # For symmetry, t_d = t_-d. E(k) = Sum 2t cos(k.d)
    
    E = 0.0
    
    # Nearest Neighbors (dist^2 = 2)
    t1 = 1.0
    for v in vectors:
        dot = kx*v[0] + ky*v[1] + kz*v[2]
        E += -t1 * np.cos(dot) # -t for stability (min at k=0)
        
    if long_range:
        # Add next-nearest?
        # FCC Next Nearest is (2,0,0). Dist^2 = 4.
        # Couping decay? 1/r^3 (Dipole)? 1/r (Coulomb)?
        # FIN K(d) ~ 1/d?
        # Let's verify standard first.
        pass
        
    return E

# 3. Calculate Group Velocity
# v_g = grad_k E(k)
# We calculate numerically or analytically.
# Analytically: d/dkx (-cos(kx+ky)) = sin(kx+ky) ...

def get_velocity(kx, ky, kz):
    grad_x = 0.0
    grad_y = 0.0
    grad_z = 0.0
    
    t1 = 1.0
    for v in vectors:
        dot = kx*v[0] + ky*v[1] + kz*v[2]
        sine = t1 * np.sin(dot)
        grad_x += sine * v[0]
        grad_y += sine * v[1]
        grad_z += sine * v[2]
        
    return np.sqrt(grad_x**2 + grad_y**2 + grad_z**2)

# 4. Scan Directions for Small k (Low Energy Limit)
# Light/Sound are low-k excitations.
# k magnitude
k_mag = 0.1 # Small vector

print(f"Scanning Anisotropy at |k| = {k_mag} ...")

# Dir 1: (1,0,0)
k1 = np.array([1, 0, 0])
k1 = k1 / np.linalg.norm(k1) * k_mag
v1 = get_velocity(k1[0], k1[1], k1[2])

# Dir 2: (1,1,0)
k2 = np.array([1, 1, 0])
k2 = k2 / np.linalg.norm(k2) * k_mag
v2 = get_velocity(k2[0], k2[1], k2[2])

# Dir 3: (1,1,1)
k3 = np.array([1, 1, 1])
k3 = k3 / np.linalg.norm(k3) * k_mag
v3 = get_velocity(k3[0], k3[1], k3[2])

print(f"v(1,0,0): {v1:.6f}")
print(f"v(1,1,0): {v2:.6f}")
print(f"v(1,1,1): {v3:.6f}")

mean_v = (v1+v2+v3)/3
anisotropy = (max(v1,v2,v3) - min(v1,v2,v3)) / mean_v * 100

print(f"\nAnisotropy Error: {anisotropy:.4f}%")

# 5. Check Long Range Correction?
# If Error > 0, does Long Range help?
# Hypothesis: Emergent Isotropy requires critical K(d).
# But first, the raw verdict.

verdict = "FAIL"
if anisotropy < 1.0: verdict = "PASS (Soft)"
if anisotropy < 1e-4: verdict = "PASS (Hard)"

print("-" * 40)
print(f"SCEPTIC VERDICT: {verdict}")
print("-" * 40)

if verdict == "FAIL":
    print("Interpretation: The lattice is visible. 'c' depends on direction.")
    print("Defense: This is valid for Phonons. For Photons, we need fine-tuning.")

# ============================================================================
# REPORT
# ============================================================================
with open("raport_qw635_isotropy.md", "w") as f:
    f.write("# Raport QW-635: Isotropy Check\n")
    f.write("**Data:** 2025-12-05\n\n")
    f.write("## Metodologia\n")
    f.write("Analiza prędkości grupowej na kracie FCC w granicy małych k (t=1.0).\n\n")
    f.write("## Wyniki\n")
    f.write(f"- v(1,0,0): {v1:.6f}\n")
    f.write(f"- v(1,1,0): {v2:.6f}\n")
    f.write(f"- v(1,1,1): {v3:.6f}\n")
    f.write(f"- **Anizotropia:** {anisotropy:.4f}%\n\n")
    
    f.write(f"## Werdykt Sceptyka: {verdict}\n")
    if anisotropy > 1.0:
        f.write("Teoria łamie BARDZO symetrię obrotową. To model ciała stałego, nie próżni.\n")

print("Report saved.")
