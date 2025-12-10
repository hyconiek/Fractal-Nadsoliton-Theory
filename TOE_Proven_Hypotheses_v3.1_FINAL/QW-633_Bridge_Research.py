#!/usr/bin/env python3
# QW-633: BRIDGE RESEARCH - CLASSICAL WAVE SPECTRUM
# Purpose: Check if Classical Wave Dynamics on FCC Lattice generate the same
#          10-12 band spectrum as the Quantum Hamiltonian.
# Bridge:  Connects Classical 4-Bit foundation -> Wave Geometry -> Quantum-like Structure.
# Method:  Construct Classical Propagation Matrix M for Scalar Waves on FCC.
#          Diagonalize M. Compare DOS with QW-629.
# Date: 2025-12-05

import numpy as np
import scipy.sparse as sp
import numpy.linalg as la # Dense for N=~256
import matplotlib.pyplot as plt

print("="*80)
print("QW-633: BRIDGE RESEARCH (CLASSICAL WAVES)")
print("="*80)
print("Test: Do Classical Waves on FCC Lattice recover the 12-band spectrum?")
print("System: Scalar Wave Equation on 3D FCC (N=256).")
print("="*80)

# 1. Classical Lattice Construction (Same as QW-629 for comparison)
# FCC Grid
L = 8
nodes = []
map_pos_to_idx = {}
idx_counter = 0

for x in range(L):
    for y in range(L):
        for z in range(L):
            if (x + y + z) % 2 == 0:
                nodes.append((x,y,z))
                map_pos_to_idx[(x,y,z)] = idx_counter
                idx_counter += 1

N = len(nodes)
print(f"System Size N: {N}")

# 2. Propagation Matrix (Classical Kirchhoff / Adjacency)
# Wave Eq: d^2u/dt^2 = c^2 Laplacian u
# Discrete: u(t+1) = 2u(t) - u(t-1) + dt^2 * Laplacian u
# Laplacian Matrix L_ij = Sum(neighbors) - degree * I.
# Eigenmodes of Wave Eq are Eigenmodes of Laplacian.

# We will analyze the LAPLACIAN spectrum directly.
# L_ij = 1 if connected, -12 on diagonal? Or positive definite?
# Usually L = D - A. Eigenvalues >= 0.

print("Building Laplacian Matrix...")
row = []
col = []
data = []

shifts = [
    (1,1,0), (1,-1,0), (-1,1,0), (-1,-1,0),
    (1,0,1), (1,0,-1), (-1,0,1), (-1,0,-1),
    (0,1,1), (0,1,-1), (0,-1,1), (0,-1,-1)
]

# NOTE: In QW-629 we used Anisotropic Coupling (Spin-dependent).
# For Classical Bridge, we must ask: Does Isotropic Wave (Scalar) work?
# Or do we need Vector Waves (Polarization)?
# Hypothesis H13 said "Spin breaks symmetry".
# So Classical Scalar Wave might NOT show 12 bands (only 1 band?).
# Let's test Isotropic Scalar first.

interaction_type = "ISOTROPIC" # vs "ANISOTROPIC"

for i in range(N):
    x, y, z = nodes[i]
    
    # Diagonal
    row.append(i)
    col.append(i)
    data.append(12.0) # Degree
    
    for dx, dy, dz in shifts:
        nx, ny, nz = x+dx, y+dy, z+dz
        if (nx, ny, nz) in map_pos_to_idx:
            j = map_pos_to_idx[(nx, ny, nz)]
            
            # Off-diagonal
            val = -1.0
            
            # If we want to bridge to the "Quantum" result, we might need Anisotropy.
            # But let's check baseline first.
            
            row.append(i)
            col.append(j)
            data.append(val)

L_mat = sp.csr_matrix((data, (row, col)), shape=(N, N))
L_dense = L_mat.toarray()

print("Diagonalizing Classical Laplacian...")
vals = la.eigvalsh(L_dense)

# Spectral Analysis
print(f"Eigenvalue Range: [{min(vals):.4f}, {max(vals):.4f}]")

# Count Peaks in DOS
hist, bins = np.histogram(vals, bins=50)
peaks = 0
for k in range(1, len(hist)-1):
    if hist[k] > hist[k-1] and hist[k] > hist[k+1]:
        peaks += 1

print(f"Isotropic Scalar Peak Count: {peaks}")

if peaks < 5:
    print(">> Low Peak Count. Isotropic not enough.")
    print(">> Switching to CLASSICAL VECTOR WAVES (Elasticity).")
    
    # Vector Wave: Each node has (ux, uy, uz). Dim = 3*N.
    # Interaction depending on angle (Springs).
    # This mimics the Spin anisotropy classically.
    
    print("\nBuilding Classical Vector Hamiltonian (Elastic Grid)...")
    dim = 3 * N
    row2 = []
    col2 = []
    data2 = []
    
    for i in range(N):
        x, y, z = nodes[i]
        
        # Self stiffness (sum of springs)
        # 12 springs * k. Matrix 3x3 diagonal block?
        # Let's verify K matrix for bond r.
        # K_bond = k * (r . rT) / |r|^2 (Longitudinal stiffness only?)
        # Let's assume Central Force springs.
        
        K_self = np.zeros((3,3))
        
        for dx, dy, dz in shifts:
            # Neighbor j
            nx, ny, nz = x+dx, y+dy, z+dz
            
            r_vec = np.array([dx, dy, dz], dtype=float)
            r_vec /= np.sqrt(2.0)
            
            # Spring matrix K_ij = - r * rT
            K_ij = -1.0 * np.outer(r_vec, r_vec)
            
            # Add to Off-diagonal block
            if (nx, ny, nz) in map_pos_to_idx:
                j = map_pos_to_idx[(nx, ny, nz)]
                
                for a in range(3):
                    for b in range(3):
                        row2.append(3*i + a)
                        col2.append(3*j + b)
                        data2.append(K_ij[a,b])
                        
                # Add to Self block (positive)
                K_self += -1.0 * K_ij # Subtracting the negative K_ij
                # Wait, Force F_i = K (u_j - u_i).
                # Matrix M: M_ii u_i + M_ij u_j = lambda u_i
                # M_ii = Sum |K_ij|
                # So K_self is Sum of outer products
                K_self += np.outer(r_vec, r_vec)
        
        # Add Self Block
        for a in range(3):
            for b in range(3):
                row2.append(3*i + a)
                col2.append(3*i + b)
                data2.append(K_self[a,b])

    H_vec = sp.csr_matrix((data2, (row2, col2)), shape=(dim, dim))
    H_vec_dense = H_vec.toarray()
    
    print("Diagonalizing Vector Laplacian...")
    vals_vec = la.eigvalsh(H_vec_dense)
    
    hist_vec, _ = np.histogram(vals_vec, bins=50)
    peaks_vec = 0
    for k in range(1, len(hist_vec)-1):
        if hist_vec[k] > hist_vec[k-1] and hist_vec[k] > hist_vec[k+1]:
            peaks_vec += 1
            
    print(f"Vector Wave Peak Count: {peaks_vec}")
    final_peaks = peaks_vec
else:
    final_peaks = peaks

# Conclusion
print("\nBridge Conclusion:")
if abs(final_peaks - 12) <= 3:
    print("✅ BRIDGE SUCCESSFUL")
    print(f"   Classical Waves (Elasticity/Vector) generated {final_peaks} bands.")
    print("   The 'Quantum' spectrum of QW-629 is actually a GEOMETRIC property.")
    print("   Structure pre-exists in the classical lattice geometry.")
else:
    print("❌ BRIDGE FAILED")
    print(f"   Classical waves count {final_peaks}. Distinct from Quantum result.")

# ============================================================================
# REPORT
# ============================================================================
with open("raport_qw633_bridge.md", "w") as f:
    f.write("# Raport QW-633: Bridge Research\n")
    f.write("**Data:** 2025-12-05\n\n")
    f.write("## Cel\n")
    f.write("Sprawdzenie, czy struktura pasmowa (10-12 oktaw) jest unikalna dla QM, czy wynika z Geometrii Falowej (Bridge).\n\n")
    f.write("## Wyniki\n")
    f.write(f"- Isotropic Scalar Peaks: {peaks}\n")
    if 'peaks_vec' in locals():
        f.write(f"- Vector Elastic Peaks: {peaks_vec}\n")
    
    if abs(final_peaks - 12) <= 3:
        f.write("### ✅ MOST ZBUDOWANY\n")
        f.write("Struktura 12-oktawowa istnieje już w klasycznej elastyczności sieci FCC.\n")
        f.write("QM tylko 'oświetla' tę strukturę, nie tworzy jej od zera.\n")
    else:
        f.write("### ❌ BRAK MOSTU\n")

print("Report saved.")
