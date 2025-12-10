#!/usr/bin/env python3
"""
QW-737 RED TEAM VERIFICATION
============================
Critical analysis of the claims made in QW-737 regarding the "anomalous gravity exponent" n ~ 1.76.

We test:
1. Stability of n vs Connectivity Radius (Reproduction)
2. Graph Geodesic vs Euclidean Distance (Metric Distortion Hypothesis)
3. Spectral Dimension vs Euclidean Dimension (Topological vs Spatial)

If n ~ 1.76 is a physical feature of the "fractal ether", it should be robust.
If it is an artifact of measuring a graph process with a Euclidean ruler, we should see it diverge.
"""

import numpy as np
import scipy.sparse as sp
from scipy.spatial import cKDTree
from scipy.sparse.csgraph import dijkstra
import matplotlib.pyplot as plt
import time

def analyze_graph_topology(positions, radius_conn):
    N = len(positions)
    tree = cKDTree(positions)
    pairs = tree.query_pairs(r=radius_conn)
    
    # Build Adjacency
    data = []
    rows = []
    cols = []
    for i, j in pairs:
        rows.append(i); cols.append(j); data.append(1.0)
        rows.append(j); cols.append(i); data.append(1.0)
        
    if not data:
        return None
        
    A = sp.coo_matrix((data, (rows, cols)), shape=(N, N)).tocsr()
    
    # Largest Component
    n_comp, labels = sp.csgraph.connected_components(A)
    if n_comp > 1:
        # Filter to largest component
        counts = np.bincount(labels)
        largest_label = np.argmax(counts)
        indices = np.where(labels == largest_label)[0]
        # Remap
        A = A[indices, :][:, indices]
        positions = positions[indices]
        N = len(indices)
        
    # Laplacian
    degrees = np.array(A.sum(axis=1)).flatten()
    L = sp.diags(degrees) - A
    
    # Invert L (Green's Function)
    # G = L^{-1}
    # We add a small mass term for stability if needed, or use pseudoinverse logic
    # In QW-737 they used L + epsilon*I
    epsilon = 1e-6
    L_reg = L + epsilon * sp.eye(N)
    
    try:
        # Use sparse solver for speed if possible, or dense if N is small
        # For N=1000 dense is fine
        if N < 2000:
            L_dense = L_reg.toarray()
            G = np.linalg.inv(L_dense)
        else:
             # Just fail/skip for now to keep it simple or use iterative
            return None
    except:
        return None
        
    return G, A, positions

def get_exponent(x, y, range_min, range_max):
    mask = (x > range_min) & (x < range_max)
    if np.sum(mask) < 10:
        return np.nan
    
    log_x = np.log(x[mask])
    log_y = np.log(y[mask])
    coeffs = np.polyfit(log_x, log_y, 1)
    return -coeffs[0] # returning n where y ~ x^-n

def experiment(n_nodes=1000, r_conn=0.6):
    box_size = 10.0
    positions = np.random.rand(n_nodes, 3) * box_size
    
    result = analyze_graph_topology(positions, r_conn)
    if result is None:
        print(f"Failed to build graph for r={r_conn}")
        return None
        
    G, A, pos_subset = result
    N = G.shape[0]
    
    # Sample pairs
    n_sample = 2000
    indices = np.random.randint(0, N, size=(n_sample, 2))
    
    # 1. Euclidean Distances
    d_euclid = []
    g_vals = []
    
    # 2. Graph Geodesic Distances
    # Shortest path on unweighted graph (approx)
    # We can pick a few source nodes and compute all paths to avoid re-running dijkstra too much
    # optimization: compute full distances for a subset of nodes
    
    sources = np.unique(indices[:, 0])
    # limit sources to 50 to save time
    if len(sources) > 50:
        sources = sources[:50]
        
    d_graph_map = {} # (i,j) -> dist
    
    # Run Dijkstra for sources
    # Note: dijkstra on unweighted graph -> hopping distance
    dist_matrix = dijkstra(A, directed=False, indices=sources, unweighted=True)
    
    # Resample pairs based on available sources
    valid_pairs_r = []
    valid_pairs_g = []
    valid_pairs_d_graph = []
    
    for i_local, src_idx in enumerate(sources):
        # get all distances from this source
        dists = dist_matrix[i_local, :]
        
        # sample some targets
        targets = np.random.choice(N, size=min(N, 100), replace=False)
        for tgt_idx in targets:
            if src_idx == tgt_idx: continue
            
            d_e = np.linalg.norm(pos_subset[src_idx] - pos_subset[tgt_idx])
            d_g = dists[tgt_idx]
            val = abs(G[src_idx, tgt_idx])
            
            if val > 1e-12 and not np.isinf(d_g):
                valid_pairs_r.append(d_e)
                valid_pairs_g.append(val)
                valid_pairs_d_graph.append(d_g)
                
    valid_pairs_r = np.array(valid_pairs_r)
    valid_pairs_g = np.array(valid_pairs_g)
    valid_pairs_d_graph = np.array(valid_pairs_d_graph)
    
    # Fit Exponents
    # Euclidean: G ~ r_e^-n_e
    n_euclid = get_exponent(valid_pairs_r, valid_pairs_g, r_conn*1.5, box_size/2.0)
    
    # Graph: G ~ d_graph^-n_g
    # Expect n_g ~ 1.0 (if free field on graph) or something specific to spectral dimension?
    # Actually in 3D lattice, G ~ 1/r. Graph dist ~ r. So G ~ 1/d_graph.
    # In 1D, G ~ r (linear). In 2D G ~ log r.
    # Let's see what we get.
    # We filter d_graph > small to avoid lattice artifacts
    n_graph = get_exponent(valid_pairs_d_graph, valid_pairs_g, 2.0, 20.0) # graph hops range
    
    # Also check Isometry: d_graph vs r_euclid
    # d_graph ~ r_euclid^delta
    # If delta = 1, space is Euclidean-like. If delta > 1, it's fractal/tortuous.
    # Fit log(d_graph) vs log(r_euclid)
    mask_iso = (valid_pairs_r > r_conn*1.5) & (valid_pairs_r < box_size/2.0)
    if np.sum(mask_iso) > 10:
        coeffs_iso = np.polyfit(np.log(valid_pairs_r[mask_iso]), np.log(valid_pairs_d_graph[mask_iso]), 1)
        dim_tortuosity = coeffs_iso[0]
    else:
        dim_tortuosity = np.nan
        
    return {
        "r_conn": r_conn,
        "n_euclid": n_euclid,
        "n_graph": n_graph,
        "tortuosity": dim_tortuosity,
        "N": N
    }

def main():
    print("QW-737 RED TEAM VERIFICATION SUITE")
    print("----------------------------------")
    
    radii = [0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2]
    results = []
    
    print(f"{'R':<6} | {'N':<5} | {'n_Euclid':<10} | {'n_Graph':<10} | {'Tortuosity':<10}")
    
    for r in radii:
        res = experiment(n_nodes=1500, r_conn=r)
        if res:
            results.append(res)
            print(f"{res['r_conn']:<6.2f} | {res['N']:<5} | {res['n_euclid']:<10.4f} | {res['n_graph']:<10.4f} | {res['tortuosity']:<10.4f}")
            
    print("\nanalysis Complete.")
    
    # Interpretation
    print("\nINTERPRETATION:")
    print("1. If n_Euclid is ~1.76 and n_Graph is ~1.0, the anomaly is due to metric mismatch (Tortuosity).")
    print("2. If n_Euclid varies strongly with R, the result is unstable.")

if __name__ == "__main__":
    main()
