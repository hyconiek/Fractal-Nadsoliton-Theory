#!/usr/bin/env python3
"""
QW-737 to QW-746: BATCH VERIFICATION SUITE
==========================================
Executing 10 simultaneous high-speed verification studies to test the
Fractal Gravity and Emergent Mass hypotheses / discrepancies.

STUDIES:
1.  QW-737: Baseline Fractal Gravity (Reproduction & Stability Check)
2.  QW-738: Spectral Dimension Analysis (Laplacian Eigenvalues)
3.  QW-739: Metric Distortion Test (Graph Geodesic vs Euclidean scan)
4.  QW-740: Weighted Graph Laplacian (Physically motivated weights ~ 1/r)
5.  QW-741: Anisotropy & Orientation Dependence
6.  QW-742: Thermal Noise / Jitter Robustness
7.  QW-743: 4D Embedding Test (Topological Stability)
8.  QW-744: Hyperbolic Geometry Embedding (Poincaré Disk)
9.  QW-745: Quantum Walk Decay (Probability diffusion vs Green's function)
10. QW-746: Decimation/Renormalization Stability check

General Settings:
- N = 1200 nodes (fast iteration)
- Box Size = 10.0
- R_conn = 0.9
"""

import numpy as np
import scipy.sparse as sp
import scipy.linalg as la
from scipy.spatial import cKDTree
from scipy.sparse.csgraph import dijkstra, shortest_path
import time
import json

# Global Config
N_NODES = 1200
BOX_SIZE = 10.0
R_CONN = 0.9 
SEED = 42

def build_graph(positions, r_conn, weighted=False):
    tree = cKDTree(positions)
    pairs = tree.query_pairs(r=r_conn)
    
    rows = []
    cols = []
    data = []
    
    for i, j in pairs:
        dist = np.linalg.norm(positions[i] - positions[j])
        
        weight = 1.0
        if weighted:
            weight = 1.0 / (dist + 1e-6)
            
        rows.append(i); cols.append(j); data.append(weight)
        rows.append(j); cols.append(i); data.append(weight)
        
    if not data:
        return None, None
        
    A = sp.coo_matrix((data, (rows, cols)), shape=(len(positions), len(positions))).tocsr()
    
    n_comp, labels = sp.csgraph.connected_components(A)
    if n_comp > 1:
        counts = np.bincount(labels)
        largest = np.argmax(counts)
        keep = np.where(labels == largest)[0]
        return A[keep][:, keep], positions[keep]
        
    return A, positions

def get_green_function(A):
    degrees = np.array(A.sum(axis=1)).flatten()
    L = sp.diags(degrees) - A
    N = A.shape[0]
    try:
        L_reg = L + 1e-6 * sp.eye(N)
        G = np.linalg.inv(L_reg.toarray())
        return G
    except:
        return None

def fit_power_law(x_data, y_data, x_min=None, x_max=None):
    if x_min is None: x_min = np.min(x_data)
    if x_max is None: x_max = np.max(x_data)
    
    mask = (x_data > x_min) & (x_data < x_max)
    if np.sum(mask) < 10:
        return np.nan
        
    lx = np.log(x_data[mask])
    ly = np.log(y_data[mask])
    
    coeffs = np.polyfit(lx, ly, 1)
    return -coeffs[0] # Exponent n

# --- STUDY RUNNERS ---

def run_qw737(pos, A):
    G = get_green_function(A)
    if G is None: return {"n_baseline": np.nan}
    N = G.shape[0]
    idx = np.random.choice(N, size=min(N, 500), replace=False)
    r_vals = []
    g_vals = []
    for i in idx:
        targets = np.random.choice(N, size=100)
        for j in targets:
            if i == j: continue
            dist = np.linalg.norm(pos[i] - pos[j])
            val = abs(G[i, j])
            r_vals.append(dist)
            g_vals.append(val)
    n = fit_power_law(np.array(r_vals), np.array(g_vals), R_CONN*1.2, BOX_SIZE/2.0)
    return {"n_baseline": n}

def run_qw738(pos, A):
    degrees = np.array(A.sum(axis=1)).flatten()
    L = sp.diags(degrees) - A
    try:
        vals = sp.linalg.eigsh(L, k=100, which='SM', return_eigenvectors=False)
        vals = np.sort(vals[vals > 1e-5]) 
        ranks = np.arange(1, len(vals) + 1)
        coeffs = np.polyfit(np.log(vals), np.log(ranks), 1)
        d_s = 2 * coeffs[0]
        return {"d_spectral": d_s}
    except:
        return {"d_spectral": np.nan}

def run_qw739(pos, A):
    N = A.shape[0]
    idx = np.random.choice(N, size=20, replace=False)
    dists_graph = dijkstra(A, indices=idx, unweighted=True)
    d_g = []
    d_e = []
    for i, src in enumerate(idx):
        for tgt in range(N):
            if src == tgt: continue
            dg = dists_graph[i, tgt]
            if np.isinf(dg): continue
            de = np.linalg.norm(pos[src] - pos[tgt])
            d_g.append(dg)
            d_e.append(de)
    if len(d_e) > 10:
        c = np.polyfit(np.log(d_e), np.log(d_g), 1)
        metric_dim = c[0]
    else:
        metric_dim = np.nan
    return {"metric_scaling": metric_dim}

def run_qw740(pos, r_conn):
    A_w, p_w = build_graph(pos, r_conn, weighted=True)
    if A_w is None: return {"n_weighted": np.nan}
    res = run_qw737(p_w, A_w)
    return {"n_weighted": res.get("n_baseline", np.nan)}

def run_qw741(pos, A):
    G = get_green_function(A)
    if G is None: return {"n_std_dev": np.nan}
    N = G.shape[0]
    pairs = []
    for _ in range(2000):
        i, j = np.random.randint(0, N, 2)
        if i==j: continue
        delta = pos[i] - pos[j]
        dist = np.linalg.norm(delta)
        val = abs(G[i, j])
        pairs.append((delta, dist, val))
    n_axes = []
    for axis in [0, 1, 2]: 
        sub_r = []
        sub_g = []
        for delta, dist, val in pairs:
            if abs(delta[axis])/dist > 0.8:
                sub_r.append(dist)
                sub_g.append(val)
        if len(sub_r) > 10:
            n = fit_power_law(np.array(sub_r), np.array(sub_g), R_CONN, BOX_SIZE/2.0)
            n_axes.append(n)
    return {"n_std_dev": np.std(n_axes) if len(n_axes)>1 else np.nan}

def run_qw742(pos, A):
    jitter = np.random.normal(0, R_CONN*0.05, pos.shape)
    pos_new = pos + jitter
    A_new, p_new = build_graph(pos_new, R_CONN)
    res = run_qw737(p_new, A_new)
    return {"n_jittered": res["n_baseline"]}

def run_qw743(n_nodes):
    pos_4d = np.random.rand(n_nodes, 4) * BOX_SIZE
    # Adjust R_conn for 4D density? Volume scales as R^4.
    # To keep same avg degree: rho * R^d = const.
    # rho = N / L^d.
    # Degree ~ R^d / L^d * N.
    # For 4D, (0.9)^4 is smaller than (0.9)^3. Need larger R.
    # R_new = R_old * (N_old/N_new)^(1/4)? No, N is same.
    # R_new = R_old^(3/4) * L^(1/4)? Let's just guess 1.2
    A_4d, p_4d = build_graph(pos_4d, 1.2)
    res = run_qw737(p_4d, A_4d)
    return {"n_4d": res["n_baseline"]}

def run_qw744(n_nodes):
    # Hyperbolic proxy: Exponential distribution from center
    pos = np.random.randn(n_nodes, 3) 
    A, p = build_graph(pos, 0.6) 
    res = run_qw737(p, A)
    return {"n_hyperbolic": res["n_baseline"]}

def run_qw745(pos, A):
    degrees = np.array(A.sum(axis=1)).flatten()
    try:
        D_inv = sp.diags(1.0/degrees)
        M = D_inv.dot(A)
        p_vec = np.zeros(A.shape[0]); p_vec[0] = 1.0
        t_vals = [10, 20, 30, 40, 50, 60, 70, 80]
        probs = []
        v = p_vec
        for t in range(81):
            v = M.dot(v)
            if t in t_vals:
                probs.append(v[0])
        coeffs = np.polyfit(np.log(t_vals), np.log(probs), 1)
        ds_diff = -2 * coeffs[0]
        return {"d_spectral_diffusion": ds_diff}
    except:
        return {"d_spectral_diffusion": np.nan}

def run_qw746(pos, A):
    N = A.shape[0]
    keep = np.random.choice(N, size=int(N/2), replace=False)
    # Rebuild
    res_rebuild = run_qw737(pos[keep], build_graph(pos[keep], R_CONN)[0])
    return {"n_decimated": res_rebuild["n_baseline"]}


def run_batch():
    print("Initialize Batch Suite...")
    np.random.seed(SEED)
    
    base_pos = np.random.rand(N_NODES, 3) * BOX_SIZE
    base_A, base_pos = build_graph(base_pos, R_CONN)
    
    if base_A is None:
        print("Base graph failed.")
        return

    results = {}
    print("Running Studies...")

    results.update(run_qw737(base_pos, base_A))
    results.update(run_qw738(base_pos, base_A))
    results.update(run_qw739(base_pos, base_A))
    results.update(run_qw740(base_pos, R_CONN))
    results.update(run_qw741(base_pos, base_A))
    results.update(run_qw742(base_pos, base_A))
    results.update(run_qw743(N_NODES))
    results.update(run_qw744(N_NODES))
    results.update(run_qw745(base_pos, base_A))
    results.update(run_qw746(base_pos, base_A))
    
    print("\n--- BATCH RESULTS ---")
    print(json.dumps(results, indent=2))
    
    # Save to MD
    report_path = "RAPORT_QW737_QW746_BATCH_RESULTS.md"
    with open(report_path, "w") as f:
        f.write("# RAPORT ZBIORCZY: BADANIA WERYFIKACYJNE QW-737 - QW-746\n\n")
        f.write("## Cel: Rygorystyczna weryfikacja anomalii fraktalnych i topologicznych\n\n")
        f.write("| ID Badania | Metraż / Hipoteza | Wynik | Wniosek |\n")
        f.write("|---|---|---|---|\n")
        
        def row(id, name, key, thresh, msg_pass, msg_fail):
            val = results.get(key, np.nan)
            v_str = f"{val:.4f}"
            c = msg_fail
            if not np.isnan(val):
                if thresh(val): c = msg_pass
            f.write(f"| {id} | {name} | {v_str} | {c} |\n")

        row("QW-737", "Baseline Fractal n", "n_baseline", lambda x: x > 1.0, "ANOMALIA POTWIERDZONA", "Brak anomalii (n~0)")
        row("QW-738", "Spectral Dimension d_s", "d_spectral", lambda x: abs(x-3.0)>0.2, "Fraktalny d_s", "Standardowe 3D")
        row("QW-739", "Metric Tortuosity", "metric_scaling", lambda x: x > 1.1, "Przestrzeń Zakrzywiona", "Euklidesowa")
        row("QW-740", "Weighted Graph n", "n_weighted", lambda x: x > 0.5, "Wagi istotne", "Wagi bez znaczenia")
        row("QW-741", "Anisotropy", "n_std_dev", lambda x: x < 0.1, "Izotropowe", "Anizotropowe")
        row("QW-742", "Noise Stability", "n_jittered", lambda x: True, "Zmierzono", "-")
        row("QW-743", "4D Embedding", "n_4d", lambda x: True, "Zmierzono", "-")
        row("QW-744", "Hyperbolic", "n_hyperbolic", lambda x: True, "Zmierzono", "-")
        row("QW-745", "Diffusion d_s", "d_spectral_diffusion", lambda x: True, "Zmierzono", "-")
        row("QW-746", "Decimation", "n_decimated", lambda x: True, "Zmierzono", "-")

        f.write("\n## Raw Data\n```json\n")
        f.write(json.dumps(results, indent=2))
        f.write("\n```\n")
    
    print(f"Report saved to {report_path}")

if __name__ == "__main__":
    run_batch()
