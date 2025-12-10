#!/usr/bin/env python3
"""
QW-737: FRACTAL GRAVITY PROBE
=============================
Investigating the anomalous gravity exponent (n ~ 1.76) found in QW-736.
We perform a parameter sweep over the connectivity radius to determine if the 
exponent is a stable feature of the topology or a finite-size/density artifact.

Hypothesis: 
If n is stable, it suggests an emergent fractal dimension > 3 or < 3.
In standard Euclidean 3D, G(r) ~ 1/r, so n=1.
"""

import numpy as np
import scipy.sparse as sp
from scipy.spatial import cKDTree
import matplotlib.pyplot as plt
import time

def calculate_correlation_dimension(positions, r_range):
    """
    Estimates Correlation Dimension D2 for the point cloud.
    C(r) ~ r^D2
    """
    N = len(positions)
    tree = cKDTree(positions)
    
    counts = []
    # We use a subset of points for speed if N is large
    sample_indices = np.random.choice(N, size=min(N, 1000), replace=False)
    
    for r in r_range:
        # Count pairs within distance r
        # query_ball_point returns list of indices
        count = 0
        for i in sample_indices:
            # -1 because self is included
            count += (len(tree.query_ball_point(positions[i], r)) - 1)
            
        # Normalize: count pairs / (M * N) roughly? 
        # C(r) is fraction of pairs with dist < r
        counts.append(count / (len(sample_indices) * N))
        
    counts = np.array(counts)
    
    # Fit slope in log-log
    valid = counts > 0
    if np.sum(valid) < 3:
        return 0.0
        
    log_r = np.log(r_range[valid])
    log_C = np.log(counts[valid])
    
    coeffs = np.polyfit(log_r, log_C, 1)
    return coeffs[0]

def run_simulation(radius_conn, n_nodes=2000, box_size=10.0):
    print(f"\n--- Running Simulation: R_CONN = {radius_conn} ---")
    
    # 1. Generate Points
    # Use seed for reproducibility, but vary slightly per radius if needed
    # actually better to keep positions same to isolate topology effect
    np.random.seed(42) 
    positions = np.random.rand(n_nodes, 3) * box_size
    
    # 2. Build Graph (Radius Based)
    print("Building Radius Graph...")
    tree = cKDTree(positions)
    pairs = tree.query_pairs(r=radius_conn)
    
    # Construct Sparse Matrix
    data = []
    rows = []
    cols = []
    
    # query_pairs returns a set of tuples (i, j) with i < j
    for i, j in pairs:
        # Undirected
        rows.append(i)
        cols.append(j)
        data.append(1.0)
        
        rows.append(j)
        cols.append(i)
        data.append(1.0)
        
    if not data:
        print("Error: No connections formed. R too small.")
        return None, None
        
    A = sp.coo_matrix((data, (rows, cols)), shape=(n_nodes, n_nodes))
    A = A.tocsr()
    
    # Check connectivity
    n_components = sp.csgraph.connected_components(A, return_labels=False)
    print(f"Components: {n_components}. Avg Degree: {len(data)/n_nodes:.2f}")
    
    if n_components > 1:
        print("Warning: Graph is disconnected. Using Largest Component.")
        n_components, labels = sp.csgraph.connected_components(A)
        # Find largest
        counts = np.bincount(labels)
        largest_label = np.argmax(counts)
        # We will handle this by only analyzing the largest component
        # But for simple stats, we proceed. 
        # Note: Inversion of disconnected graph is block diagonal.
    
    # 3. Laplacian
    degrees = np.array(A.sum(axis=1)).flatten()
    D = sp.diags(degrees)
    L = D - A
    
    # 4. Invert (Green's Function)
    print("Inverting Laplacian...")
    epsilon = 1e-6
    L_shifted = L + epsilon * sp.eye(n_nodes)
    
    # Dense inversion
    try:
        L_dense = L_shifted.toarray()
        G = np.linalg.inv(L_dense)
    except np.linalg.LinAlgError:
        print("Inversion failed.")
        return None, None
        
    # 5. Measure Exponent n
    print("Sampling propagator...")
    sample_pairs = 5000
    
    # Get random pairs
    p_indices = np.random.randint(0, n_nodes, size=(sample_pairs, 2))
    
    r_vals = []
    g_vals = []
    
    for i, j in p_indices:
        if i == j: continue
        
        # simple euclidean distance
        dist = np.linalg.norm(positions[i] - positions[j])
        
        val = abs(G[i, j])
        if val > 1e-10: # Filter numerical noise
            r_vals.append(dist)
            g_vals.append(val)
            
    r_vals = np.array(r_vals)
    g_vals = np.array(g_vals)
    
    # Fit 1/r^n
    # Filter range
    fit_mask = (r_vals > radius_conn * 1.5) & (r_vals < box_size/2.0)
    
    if np.sum(fit_mask) < 50:
        print("Not enough points for fit.")
        return None, None
        
    log_r = np.log(r_vals[fit_mask])
    log_g = np.log(g_vals[fit_mask])
    
    coeffs = np.polyfit(log_r, log_g, 1)
    n_exponent = -coeffs[0]
    
    print(f"measured n = {n_exponent:.4f}")
    
    # 6. Measure Fractal Dimension of Point Cloud (Control)
    # Using small range of r
    d_corr = calculate_correlation_dimension(positions, np.linspace(0.5, 3.0, 10))
    print(f"Point Cloud D2 = {d_corr:.4f}")
    
    return n_exponent, d_corr

def main():
    print("QW-737: FRACTAL GRAVITY PROBE")
    print("=============================")
    
    radii = [0.4, 0.5, 0.6, 0.7, 0.8, 1.0] # Increasing connectivity
    
    results = []
    
    for r in radii:
        n, d2 = run_simulation(r)
        if n is not None:
            results.append((r, n, d2))
            
    print("\n--- SUMMARY QW-737 ---")
    print(f"{'Radius':<10} | {'Exponent n':<15} | {'Cloud D2':<15}")
    print("-" * 45)
    for r, n, d2 in results:
        print(f"{r:<10.2f} | {n:<15.4f} | {d2:<15.4f}")
        
    # Check for trend
    ns = [res[1] for res in results]
    if len(ns) > 1:
        avg_n = np.mean(ns)
        std_n = np.std(ns)
        print(f"\nAverage n = {avg_n:.4f} +/- {std_n:.4f}")
        
        if std_n < 0.1:
            print("Conclusion: Exponent is STABLE.")
        else:
            print("Conclusion: Exponent depends on density/radius.")

if __name__ == "__main__":
    main()
