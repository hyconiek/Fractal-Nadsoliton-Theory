#!/usr/bin/env python3
"""
QW-735: EMERGENT GRAVITY PROPAGATOR TEST
========================================
Cel: Sprawdzić, czy w losowej sieci 3D z połączeniami LOKALNYMI (Nearest Neighbors)
     funkcja Greena (odwrotność Laplasjanu) dąży do 1/r.

Jeśli TAK: Grawitacja (1/r) jest emergentną wlasnością topologii 3D.
Jeśli NIE: Grawitacja wymaga nielokalnych połączeń "z ręki".
"""

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from scipy.spatial import cKDTree
import matplotlib.pyplot as plt

def main():
    print("QW-735: EMERGENT 1/r PROPAGATOR TEST")
    print("====================================")
    
    # 1. Constants
    N_NODES = 2000
    K_NEIGHBORS = 6   # Connectivity (Kissing Number-like)
    BOX_SIZE = 10.0
    
    print(f"Network: {N_NODES} nodes in {BOX_SIZE}^3 box")
    print(f"Topology: Only {K_NEIGHBORS} nearest neighbors (LOCAL interactions only)")
    
    # 2. Generate Random Geometric Graph
    np.random.seed(735)
    positions = np.random.rand(N_NODES, 3) * BOX_SIZE
    
    print("Building adjacency matrix...")
    tree = cKDTree(positions)
    dists, indices = tree.query(positions, k=K_NEIGHBORS+1) # +1 because self is included
    
    # Construct Sparse Adjacency Matrix
    row_ind = []
    col_ind = []
    data = []
    
    for i in range(N_NODES):
        for j_idx in range(1, K_NEIGHBORS+1): # Skip self
            neighbor = indices[i, j_idx]
            # Undirected graph (symmetric)
            row_ind.append(i)
            col_ind.append(neighbor)
            data.append(1.0)
            
            row_ind.append(neighbor)
            col_ind.append(i)
            data.append(1.0)
            
    A = sp.coo_matrix((data, (row_ind, col_ind)), shape=(N_NODES, N_NODES))
    A = A.tocsr()
    
    # 3. Compute Laplacian: L = D - A
    degrees = np.array(A.sum(axis=1)).flatten()
    D = sp.diags(degrees)
    L = D - A
    
    print("Computing Pseudo-Inverse of Laplacian (Green's Function)...")
    # G = L+ (pinv)
    # Since L has a null space (constant vector), proof requires removing it or using shift
    # We use a small shift: (L + epsilon*I)^-1
    epsilon = 1e-6
    L_shifted = L + epsilon * sp.eye(N_NODES)
    
    # Invert (this is expensive, O(N^3), but okay for N=2000)
    # Using dense inversion for full analysis
    L_dense = L_shifted.toarray()
    G = np.linalg.inv(L_dense)
    
    # 4. Analyze G(r)
    print("Analyzing Propagator G(r)...")
    
    sample_pairs = 5000
    r_samples = []
    g_samples = []
    
    pairs = np.random.randint(0, N_NODES, size=(sample_pairs, 2))
    
    for i, j in pairs:
        if i == j: continue
        r = np.linalg.norm(positions[i] - positions[j])
        g_val = abs(G[i, j]) # Absolute coupling
        
        r_samples.append(r)
        g_samples.append(g_val)
        
    r_samples = np.array(r_samples)
    g_samples = np.array(g_samples)
    
    # 5. Fit Power Law: G ~ 1/r^n
    # Log-Log Fit
    mask = (r_samples > 0.5) & (r_samples < BOX_SIZE/2) # Avoid short range (lattice artifacts) and boundaries
    
    if np.sum(mask) > 100:
        log_r = np.log(r_samples[mask])
        log_g = np.log(g_samples[mask])
        
        # Linear regression
        coeffs = np.polyfit(log_r, log_g, 1)
        exponent = -coeffs[0]
        prefactor = np.exp(coeffs[1])
        
        print("\nRESULT QW-735:")
        print(f"Fit G(r) ~ 1 / r^{exponent:.4f}")
        print(f"Theoretical Expectation for 3D: n = 1.0 (Coulomb/Gravity potential)")
        
        error = abs(exponent - 1.0) * 100
        print(f"Error: {error:.2f}%")
        
        if error < 20.0:
            print("\n✅ SUCCESS: 1/r Potential EMERGES from Local Topology!")
            print("   Gravity is a statistical necessity of 3D connectivity.")
        else:
            print(f"\n❌ FAILURE: Exponent {exponent:.2f} does not match 1.0.")
            print("   (Perhaps simulation box is too small or boundaries dominate)")
            
    else:
        print("Not enough data points in valid range.")

if __name__ == "__main__":
    main()
