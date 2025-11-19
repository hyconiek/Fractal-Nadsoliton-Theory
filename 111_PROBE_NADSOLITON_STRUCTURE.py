#!/usr/bin/env python3
# Author: Krzysztof Żuchowski

"""
111_PROBE_NADSOLITON_STRUCTURE.py
Szybkie, bez-fitowe testy diagnostyczne wewnętrznej struktury nadsolitona.
6 zadań (krótko):
- Task A: Shannon entropy top-k eigenvectors (blokowość)
- Task B: Klasteryzacja komponentów top-4 wektorów własnych
- Task C: Komutatorowy test (norma [S1,S2]) dla przesuniętych macierzy (generator proxy)
- Task D: Beta-proxy mapping: d ln lambda/d ln s oraz liczenie zmian znaku
- Task E: Czułość na perturbację fazy (phi -> phi+δ) i fidelity top eigenvector
- Task F: Masa→sprzężenie czułość (delta lambda przy ±1% alpha/beta)

Uruchom: python3 111_PROBE_NADSOLITON_STRUCTURE.py --dry-run
"""

from __future__ import annotations
import numpy as np
from scipy.linalg import eigh
import json
from datetime import datetime, timezone
import os
import argparse
from scipy.spatial.distance import cdist

OUT_JSON = "report_111_probe_nadsoliton_structure.json"
OUT_MD = "report_111_probe_nadsoliton_structure.md"

# Parametry domyślne
ALPHA = 2.7715
BETA = 0.01
OMEGA = 2*np.pi/8.0
PHI = 0.0
N_DEFAULT = 24
SCALES = np.linspace(0.5, 2.0, 25)

# Helpers
def kernel_matrix(N, alpha, beta, omega=OMEGA, phi=PHI):
    coords = np.arange(N)
    d = np.abs(coords.reshape(-1,1) - coords.reshape(1,-1))
    K = alpha * np.cos(omega * d + phi) / (1.0 + beta * d)
    return (K + K.T)/2.0

def top_eig(S, k=6):
    vals, vecs = eigh(S)
    order = np.argsort(vals)[::-1]
    return vals[order][:k], vecs[:,order][:,:k]

# Tasks

def task_A_shannon_entropy(N, alpha, beta):
    S = kernel_matrix(N, alpha, beta)
    vals, vecs = top_eig(S, k=4)
    entropies = []
    for i in range(vecs.shape[1]):
        p = np.abs(vecs[:,i])**2
        p = p / (p.sum() + 1e-16)
        H = -np.sum(p * np.log(np.maximum(p,1e-16)))
        entropies.append(float(H))
    return {"entropies": entropies, "avg_entropy": float(np.mean(entropies))}


def task_B_clustering(N, alpha, beta, k_clusters=3):
    S = kernel_matrix(N, alpha, beta)
    _, vecs = top_eig(S, k=4)
    # use absolute component magnitudes as features
    X = np.abs(vecs)
    # simple KMeans implementation (to avoid sklearn dependency)
    def simple_kmeans(X, k, maxiter=100):
        # X: (n_samples, n_features)
        n_samples = X.shape[0]
        # initialize centers randomly from samples
        rng = np.random.default_rng(0)
        centers = X[rng.choice(n_samples, size=k, replace=False)]
        labels = np.zeros(n_samples, dtype=int)
        for _ in range(maxiter):
            dists = cdist(X, centers, metric='euclidean')
            new_labels = np.argmin(dists, axis=1)
            if np.array_equal(new_labels, labels):
                break
            labels = new_labels
            for j in range(k):
                pts = X[labels == j]
                if len(pts) > 0:
                    centers[j] = pts.mean(axis=0)
        return labels, centers

    # cluster rows (components) of X
    X_rows = X
    labels, centers = simple_kmeans(X_rows, k_clusters)
    labels = labels.tolist()
    centers = centers.tolist()
    return {"labels": labels, "centers_shape": [len(centers), len(centers[0])]}


def task_C_commutator_test(N, alpha, beta, shift=1):
    S = kernel_matrix(N, alpha, beta)
    S_shift = kernel_matrix(N, alpha, beta, phi=PHI + 0.123)  # small phase shift
    comm = S.dot(S_shift) - S_shift.dot(S)
    norm = float(np.linalg.norm(comm))
    return {"comm_norm": norm}


def task_D_beta_proxy(alpha, beta, N, s_values=SCALES):
    g = []
    for s in s_values:
        S = kernel_matrix(N, alpha*s, beta)
        vals, _ = eigh(S)
        g.append(np.sort(vals)[-1])
    ln_s = np.log(s_values)
    g = np.array(g)
    dg_dlns = np.gradient(g, ln_s)
    sign_changes = int(np.sum(np.diff(np.sign(dg_dlns))!=0))
    return {"s": s_values.tolist(), "g": g.tolist(), "dg_dlns": dg_dlns.tolist(), "sign_changes": sign_changes}


def task_E_perturbation_fidelity(N, alpha, beta, delta_phi=1e-2):
    S0 = kernel_matrix(N, alpha, beta, phi=PHI)
    S1 = kernel_matrix(N, alpha, beta, phi=PHI + delta_phi)
    v0 = top_eig(S0, k=1)[1][:,0]
    v1 = top_eig(S1, k=1)[1][:,0]
    fid = float(np.abs(np.dot(np.conjugate(v0), v1)))
    return {"delta_phi": delta_phi, "fidelity": fid}


def task_F_mass_coupling_sensitivity(N, alpha, beta):
    deltas = [0.99, 1.0, 1.01]
    out = {}
    for da in deltas:
        for db in deltas:
            a = alpha * da
            b = beta * db
            S = kernel_matrix(N, a, b)
            vals, _ = eigh(S)
            out[f"a{da}_b{db}"] = float(np.sort(vals)[-1])
    return out


def run_all(dry_run=False, N=N_DEFAULT, alpha=ALPHA, beta=BETA):
    meta = {"created": datetime.now(timezone.utc).isoformat(), "script": os.path.basename(__file__)}
    results = {}
    results["A_entropy"] = task_A_shannon_entropy(N, alpha, beta)
    results["B_clustering"] = task_B_clustering(N, alpha, beta)
    results["C_commutator"] = task_C_commutator_test(N, alpha, beta)
    results["D_beta_proxy"] = task_D_beta_proxy(alpha, beta, N)
    results["E_perturb_fid"] = task_E_perturbation_fidelity(N, alpha, beta)
    results["F_mass_sens"] = task_F_mass_coupling_sensitivity(N, alpha, beta)
    out = {"meta": meta, "results": results}
    with open(OUT_JSON, "w") as f:
        json.dump(out, f, indent=2)
    # md
    with open(OUT_MD, "w") as f:
        f.write(f"# Raport 111 — Probe Nadsolitona\nData: {meta['created']}\n")
        for k,v in results.items():
            f.write(f"\n## {k}\n```")
            f.write(json.dumps(v, indent=2))
            f.write("\n```")
    return out

if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--dry-run', action='store_true')
    p.add_argument('--N', type=int, default=N_DEFAULT)
    args = p.parse_args()
    report = run_all(dry_run=args.dry_run, N=args.N)
    print('Zapisano raport:', OUT_JSON)
