#!/usr/bin/env python3
# QW-634: LATTICE HYDROGEN SPECTRUM
# Purpose: Derive the Hydrogen Spectrum (Rydberg series 1/n^2) from the FCC Lattice.
#          This tests if the discrete geometry reproduces continuous atomic physics.
# System: Single electron on FCC Lattice (N=~2000) with central Coulomb potential.
# Date: 2025-12-05

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import matplotlib.pyplot as plt

print("="*80)
print("QW-634: HYDROGEN SPECTRUM ON FCC LATTICE")
print("="*80)
print("Test: Do eigenvalues follow E_n ~ -1/n^2?")
print("System: 3D FCC Lattice + Coulomb Potential V(r) ~ -1/r.")
print("="*80)

# 1. Construct FCC Lattice (Larger L for better spectrum)
L = 10 # 10x10x10. N approx 4 * 1000 = 4000.
# Coordinates need to be centered to have a "center" for the proton.
# Let's use range [-L, L].

print(f"Constructing FCC Lattice (Range -{L}..{L})...")

nodes = []
map_pos_to_idx = {}
idx_counter = 0

# Scan volume
# FCC criterion: x,y,z integers, sum even (D3 lattice).
# Scale: Distance between nearest neighbors is sqrt(2).
# We want "r" in physical units.
# Lattice constant a_lat. Let's assume natural units where 1/m_eff * 1/a^2 ...
# We just look for the RATIO of energy levels.

for x in range(-L, L+1):
    for y in range(-L, L+1):
        for z in range(-L, L+1):
            if (x + y + z) % 2 == 0:
                nodes.append((x,y,z))
                map_pos_to_idx[(x,y,z)] = idx_counter
                idx_counter += 1

N = len(nodes)
print(f"System Size N: {N}")

# 2. Build Hamiltonian H = T + V
# Kinetic T: Laplacian on FCC.
# Potential V: -1/r

print("Building Hamiltonian...")
row = []
col = []
data = []

shifts = [
    (1,1,0), (1,-1,0), (-1,1,0), (-1,-1,0),
    (1,0,1), (1,0,-1), (-1,0,1), (-1,0,-1),
    (0,1,1), (0,1,-1), (0,-1,1), (0,-1,-1)
]

# Coulomb Potential Parameters
# Regularize at r=0 to avoid infinity. V(0) = -V_max.
# Center is (0,0,0).
# Regularization: V(r) = -1 / sqrt(r^2 + a^2)
epsilon = 0.5 
Z_proton = 1.0
coupling_scale = 10.0 # Adjust to get binding

diag_kinetic = 12.0 # Degree

for i in range(N):
    x, y, z = nodes[i]
    r_sq = x*x + y*y + z*z
    r = np.sqrt(r_sq)
    
    # Potential Term
    # Note: On lattice, x,y,z are integers. Distances are sqrt(even).
    if r < 1e-6:
        V = -Z_proton * coupling_scale / epsilon # Origin
    else:
        # Distance unit? Neighbor=sqrt(2).
        # Let's scale potential so ground state is deep.
        V = -Z_proton * coupling_scale / r
    
    # Kinetic Diagonal
    idx = i
    
    # H_ii = Kinetic_ii + V_i
    # Kinetic: Laplacian = sum(neighbors) - diag.
    # Usually H = -t * sum(neighbors) + V.
    # T = -Laplacian = Diag - Sum(Off).
    
    t_hop = 1.0
    H_ii = t_hop * diag_kinetic + V
    
    row.append(i)
    col.append(i)
    data.append(H_ii)
    
    # Off-diagonal hopping
    for dx, dy, dz in shifts:
        nx, ny, nz = x+dx, y+dy, z+dz
        if (nx, ny, nz) in map_pos_to_idx:
            j = map_pos_to_idx[(nx, ny, nz)]
            
            row.append(i)
            col.append(j)
            data.append(-t_hop) # Standard hopping

H = sp.csr_matrix((data, (row, col)), shape=(N, N))

print("Diagonalizing (Lowest States)...")
# We need lowest eigenvalues (most negative). "SA" = Smallest Algebraic.
k_levels = 20
vals, vecs = spla.eigsh(H, k=k_levels, which='SA')

print("\nEnergy Levels (E_n):")
print(vals)

# Analysis: Identify Degenerate Multiplets
# Hydrogen: 
# n=1 (1s): 1 state
# n=2 (2s, 2p): 1+3 = 4 states (Degenerate in pure Coulomb)
# n=3 (3s, 3p, 3d): 1+3+5 = 9 states
# Energies: E1, E2=E1/4, E3=E1/9 ... (relative to continuum 0)

# Note: Lattice energy floor is not 0.
# Continuum starts at E_min_band?
# No, "Vacuum" is when electron is free.
# Free electron dispersion: E(k) ~ k^2. Bottom of band is 0 if shifted.
# My Hamiltonian: Diagonal = +12. Off = -1. Sum = 0 for k=0 uniform.
# So Bottom of Band is 0.
# Bound states should be < 0.

bound_states = vals[vals < 0]
print(f"Found {len(bound_states)} bound states < 0.")

if len(bound_states) > 0:
    E_ground = bound_states[0]
    print(f"Ground State E1: {E_ground:.4f}")
    
    # Check Ratios
    # Theoretical: E_n = E1 / n^2
    # E2 = E1 / 4
    # E3 = E1 / 9
    
    print("\nChecking 1/n^2 Scaling:")
    print("n | Theory (E1/n^2) | Actual | Error")
    print("-" * 40)
    
    # Group actual values by degeneracy/closeness to find "shells"
    # Find unique levels with tolerance
    unique_levels = []
    tol = 0.5 # Gap tolerance
    current_cluster = [vals[0]]
    
    for v in vals[1:]:
        if abs(v - current_cluster[0]) < tol:
            current_cluster.append(v)
        else:
            unique_levels.append(np.mean(current_cluster))
            current_cluster = [v]
    unique_levels.append(np.mean(current_cluster))
    
    for n in range(1, 4):
        theory = E_ground / (n**2)
        if n-1 < len(unique_levels):
            actual = unique_levels[n-1]
            err = abs((actual - theory)/theory) * 100
            print(f"{n} | {theory:.4f}        | {actual:.4f} | {err:.1f}%")
        else:
            print(f"{n} | {theory:.4f}        | ???    | N/A")

# ============================================================================
# REPORT
# ============================================================================
with open("raport_qw634_hydrogen_lattice.md", "w") as f:
    f.write("# Raport QW-634: Lattice Hydrogen Spectrum\n")
    f.write("**Data:** 2025-12-05\n\n")
    f.write("## Cel\n")
    f.write("Sprawdzenie, czy dyskretna krata FCC odtwarza widmo wodorowe $1/n^2$.\n\n")
    f.write("## Wyniki\n")
    f.write(f"- Ground State E1: {E_ground if len(bound_states)>0 else 'None'}\n")
    f.write("### Poziomy Energetyczne:\n")
    for v in vals:
        f.write(f"- {v:.4f}\n")

print("Report saved.")
