#!/usr/bin/env python3
# QW-636: PARITY & CHIRALITY CHECK (SCEPTIC'S CHALLENGE)
# Purpose: Test if FIN Lattice Theory suffers from Fermion Doubling (Nielsen-Ninomiya).
#          Test if Energy Spectrum is Parity Invariant E(k) = E(-k).
#          Verify if Long-Range Coupling (1/r) evades Doubling.
# Date: 2025-12-05

import numpy as np
import matplotlib.pyplot as plt

print("="*80)
print("QW-636: PARITY & CHIRALITY CHECK")
print("="*80)
print("Test 1: Fermion Doubling (Count zero crossings in dispersion)")
print("Test 2: Parity Symmetry (Is E(k) = E(-k)?)")
print("Hamiltonian: Long-Range Hopping t(r) ~ 1/r^alpha")
print("="*80)

# 1. Define 1D Lattice for clear visualization of Doubling
# (3D is just product of 1D in simple cases).
# L nodes.
L = 1000
r_vec = np.arange(1, L//2) 

# Coupling Parameters
# Short Range: Only r=1
# Long Range: r^-alpha
alpha = 1.0 # Coulomb-like / FIN Geometric
# alpha = 2.0 (Dipole)

def compute_dispersion(k_vals, alpha_val, long_range=True):
    E_vals = []
    # E(k) = Sum_r  2 * t(r) * cos(k*r)  (Assuming symmetric hopping t(r)=t(-r))
    # If hopping is antisymmetric t(-r) = -t(r), we get i*sin(k*r).
    # FIN Theory: Coupling K(d) is symmetric on distance distance d.
    # So E(k) is cosine series.
    
    for k in k_vals:
        E = 0.0
        # Nearest Neighbor
        t1 = 1.0
        E += -2 * t1 * np.cos(k * 1)
        
        if long_range:
            # Add neighbors r=2 to L/2
            for r in r_vec[1:]:
                # t(r) ~ 1/r^alpha. Sign?
                # Usually hopping is negative (-t).
                # Alternating sign? or constant?
                # Let's assume constant sign decay.
                tr = 1.0 / (r ** alpha_val)
                E += -2 * tr * np.cos(k * r)
        
        E_vals.append(E)
    return np.array(E_vals)

k_range = np.linspace(-np.pi, np.pi, 200)

print("Computing Bands...")

# Case A: Nearest Neighbor (Standard)
E_nn = compute_dispersion(k_range, alpha, long_range=False)

# Case B: Long Range (alpha=1, FIN)
E_fin = compute_dispersion(k_range, 1.0, long_range=True)

# Case C: Long Range (alpha=2)
E_dip = compute_dispersion(k_range, 2.0, long_range=True)

# Analyze Doubling (Zero Crossings)
# Shift energy so global min is at E_min (Vacuum).
# Standard Fermion doubling usually refers to bands crossing the Fermi level multiple times.
# Let's count minima.
# Standard Cosine: Min at 0. Max at pi. Monotonic in 0..pi.
# If Dispersion is Monotonic in 0..pi, then there is NO Doubling (Single species).
# If there are extra valleys, we have doublers.

def count_minima(E_array):
    # Only check 0 to pi ( indices 100 to 200 )
    half = E_array[100:]
    deriv = np.diff(half)
    # Count sign changes of derivative
    sign_changes = 0
    for i in range(len(deriv)-1):
        if deriv[i] * deriv[i+1] < 0:
            sign_changes += 1
    return sign_changes # 0 implies monotonic (1 minimum at start)

doublers_nn = count_minima(E_nn)
doublers_fin = count_minima(E_fin)

print(f"Standard (NN) Extra Minima: {doublers_nn}")
print(f"FIN (1/r) Extra Minima: {doublers_fin}")

# Check Parity Asymmetry
# Is E(k) == E(-k)?
# Since our coupling K(|r|) depends only on modulus, t(r)=t(-r).
# So Hamiltonian is symmetric. E(k) is even function.
diff_nn = np.linalg.norm(E_nn - np.flip(E_nn))
diff_fin = np.linalg.norm(E_fin - np.flip(E_fin))

print(f"\nParity Asymmetry Check (E_L - E_R):")
print(f"Standard: {diff_nn:.6f}")
print(f"FIN: {diff_fin:.6f}")

if diff_fin < 1e-6:
    print("\n❌ RESULT: THEORY IS PARITY SYMMETRIC.")
    print("   E(k) = E(-k). Left and Right states have same energy.")
    print("   We CANNOT explain Weak Interaction (Chirality) with scalar coupling.")
    print("   Prof. Scepticus WINS this point.")
else:
    print("\n✅ RESULT: PARITY BROKEN.")

# Solution Search:
# How to break parity? t(r) must NOT be t(-r).
# We need Directional Coupling.
# t(r) = K(r) * (1 + Spin . r) ?
# Let's test "Spin-Orbit Hopper".
# t(r) = (1/r) * sign(r) ? No, hopping term usually Hermitian conjugate c_i+ c_j + c_j+ c_i.
# Requires complex phase? Peierls phase?
# If we add magnetic flux (Hopfion?), we break T-symmetry.
# What breaks P-symmetry?
# Chirality term Gamma_5. 

print("\nTrying 'Hopfion Phase' Coupling (Complex)...")
# t(r) = (1/r) * exp(i * theta * sign(r))
def compute_chiral_dispersion(k_vals):
    E_vals = []
    theta = np.pi / 2 # Max chirality
    for k in k_vals:
        E = 0.0
        # Sum over r > 0
        for r in r_vec:
            # H_ij = t. H_ji = t*. 
            # contribution = t * e^(ikr) + t* * e^(-ikr)
            tr = (1.0 / r) * np.exp(1j * theta)
            val = tr * np.exp(1j * k * r) + np.conj(tr) * np.exp(-1j * k * r)
            E += val.real
        E_vals.append(E)
    return np.array(E_vals)

E_chiral = compute_chiral_dispersion(k_range)
diff_chiral = np.linalg.norm(E_chiral - np.flip(E_chiral))
print(f"Chiral Phase Model Asymmetry: {diff_chiral:.6f}")

if diff_chiral > 1.0:
    print(">> A Complex Phase (Hopfion Flux?) CAN break Parity!")
    print("   Hypothesis: FIN Network nodes have internal phase (Hopfion).")

# ============================================================================
# REPORT
# ============================================================================
with open("raport_qw636_parity.md", "w") as f:
    f.write("# Raport QW-636: Parity Check\n")
    f.write("**Data:** 2025-12-05\n\n")
    f.write("## Wyniki\n")
    f.write(f"- Standard Parity Asymmetry: {diff_nn:.6f}\n")
    f.write(f"- FIN Scalar Asymmetry: {diff_fin:.6f}\n\n")
    f.write(f"- Chiral Phase Asymmetry: {diff_chiral:.6f}\n")
    
    if diff_fin < 1e-4:
        f.write("### ❌ PARITY CONSERVED (Scalar Model)\n")
        f.write("Podstawowy model FIN (skalarny) nie łamie parzystości. Nie wyjaśnia neutrin.\n")
    
    if diff_chiral > 1.0:
        f.write("### ✅ PARITY BROKEN (Hopfion Model)\n")
        f.write("Wprowadzenie zespolonej fazy (Flux/Winding) do wiązań łamie symetrię P.\n")
        f.write("Wymaga to uznania węzłów za 'Hopfiony' z wewnętrzną fazą.\n")

print("Report saved.")
