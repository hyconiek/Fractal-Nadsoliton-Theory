#!/usr/bin/env python3
# Author: Krzysztof Żuchowski

"""
112_ANALYZE_NADSOLITON_DEEP.py
Głębokie, bez-fitowe analizy nadsolitona:
- Task 1: Algebraic probe (projections P_i, komutatory, ocena zamknięcia jako Lie-algebra proxy)
- Task 2: Participation Ratio (PR) scaling z N dla top-1..4 modów
- Task 3: Topological defect probe (lokalizacja modów po wprowadzeniu defektu)
- Task 4: Generator reconstruction (SVD przestrzeni S(phi) -> wymiar bazy generatorów)
- Task 5: RG landscape sweep (beta-proxy sign changes dla g(s) po scales)

Uruchomienia: python3 112_ANALYZE_NADSOLITON_DEEP.py --Ns 12 16 20 24 28 32
"""
from __future__ import annotations
import numpy as np
from scipy.linalg import eigh
import json
from datetime import datetime, timezone
import os
import argparse

OUT_JSON = 'report_112_analyze_nadsoliton_deep.json'
OUT_MD = 'report_112_analyze_nadsoliton_deep.md'

ALPHA = 2.7715
BETA = 0.01
OMEGA = 2*np.pi/8.0
PHI = 0.0
SCALES = np.linspace(0.5, 2.0, 61)

def kernel_matrix(N, alpha=ALPHA, beta=BETA, omega=OMEGA, phi=PHI):
    coords = np.arange(N)
    d = np.abs(coords.reshape(-1,1)-coords.reshape(1,-1))
    K = alpha * np.cos(omega * d + phi) / (1.0 + beta * d)
    return (K + K.T)/2.0

# Task 1: Algebraic probe
def algebraic_probe(S, top_k=4):
    vals, vecs = eigh(S)
    order = np.argsort(vals)[::-1]
    vals_k = vals[order][:top_k]
    vecs_k = vecs[:,order][:,:top_k]
    # Projections P_i = v_i v_i^T
    P = []
    for i in range(top_k):
        v = vecs_k[:,i:i+1]
        P.append((v @ v.T).astype(float))
    # Compute commutators C_ij = [P_i, P_j] and check closure onto span(P)
    topk_flat = np.stack([p.reshape(-1) for p in P], axis=1)
    # Gram matrix for projection basis
    G = np.dot(topk_flat.T, topk_flat)
    closure = {}
    f_coeffs = np.zeros((top_k, top_k, top_k))
    residual_norms = np.zeros((top_k, top_k))
    for i in range(top_k):
        for j in range(top_k):
            C = P[i] @ P[j] - P[j] @ P[i]
            c_vec = C.reshape(-1)
            # solve least squares for coefficients a in span(P): topk_flat @ a ≈ c_vec
            # use normal equations with regularization
            A = topk_flat
            rhs = np.dot(A.T, c_vec)
            try:
                coeffs = np.linalg.solve(np.dot(A.T, A) + 1e-12*np.eye(top_k), rhs)
            except np.linalg.LinAlgError:
                coeffs = np.linalg.lstsq(A, c_vec, rcond=None)[0]
            approx = A @ coeffs
            res = c_vec - approx
            residual_norms[i,j] = float(np.linalg.norm(res))
            # store coefficients mapped to basis P_k via trace normalization
            for k in range(top_k):
                f_coeffs[i,j,k] = float(coeffs[k])
    closure['residual_norms'] = residual_norms.tolist()
    closure['f_coeffs'] = f_coeffs.tolist()
    # metrics: relative residual (norm residual / norm C)
    rel_res = np.zeros((top_k, top_k))
    for i in range(top_k):
        for j in range(top_k):
            C = P[i] @ P[j] - P[j] @ P[i]
            normC = np.linalg.norm(C)
            rel_res[i,j] = float(residual_norms[i,j] / (normC + 1e-16))
    closure['rel_residual'] = rel_res.tolist()
    # summary: fraction of pairs with rel_residual < 1e-2, <1e-1
    thresh1 = np.mean(rel_res < 1e-2)
    thresh2 = np.mean(rel_res < 1e-1)
    closure['fraction_below_1e-2'] = float(thresh1)
    closure['fraction_below_1e-1'] = float(thresh2)
    return {'vals_k': vals_k.tolist(), 'closure': closure}

# Task 2: Participation Ratio scaling
def participation_ratio(vec):
    p = np.abs(vec)**2
    p = p / (p.sum() + 1e-16)
    PR = 1.0 / np.sum(p**2)
    return float(PR)

def PR_scaling(N_list, top_k=4):
    out = {}
    for N in N_list:
        S = kernel_matrix(N)
        vals, vecs = eigh(S)
        order = np.argsort(vals)[::-1]
        vecs_k = vecs[:,order][:,:top_k]
        prs = [participation_ratio(vecs_k[:,i]) for i in range(top_k)]
        out[N] = prs
    # fit power-law PR ~ c * N^alpha for each mode index
    prs_array = np.array([out[N] for N in N_list])
    logN = np.log(np.array(N_list))
    alphas = []
    fits = {}
    for mode in range(top_k):
        y = np.log(prs_array[:,mode] + 1e-12)
        coef = np.polyfit(logN, y, 1)
        alpha = float(coef[0])
        c = float(np.exp(coef[1]))
        alphas.append(alpha)
        fits[f'mode_{mode+1}'] = {'alpha': alpha, 'c': c}
    return {'prs': out, 'fits': fits}

# Task 3: Topological defect probe
def defect_probe(N, defect_width=1, defect_pos=None):
    if defect_pos is None:
        defect_pos = N//2
    S = kernel_matrix(N)
    S_def = S.copy()
    # zero out rows/cols in window defect_pos - defect_width .. + defect_width
    idx = np.arange(max(0, defect_pos-defect_width), min(N, defect_pos+defect_width+1))
    S_def[idx,:] = 0.0
    S_def[:,idx] = 0.0
    vals, vecs = eigh(S)
    vals_def, vecs_def = eigh(S_def)
    # compare top eigenvectors localization
    order = np.argsort(vals)[::-1]
    order_def = np.argsort(vals_def)[::-1]
    topv = vecs[:,order[0]]
    topv_def = vecs_def[:,order_def[0]]
    PR_orig = participation_ratio(topv)
    PR_def = participation_ratio(topv_def)
    # delta lambda
    delta_lambda = float(vals[order[0]] - vals_def[order_def[0]])
    return {'PR_orig': PR_orig, 'PR_def': PR_def, 'delta_lambda': delta_lambda}

# Task 4: Generator reconstruction via SVD of S(phi)
def generator_reconstruction(N, n_phi=9):
    phis = np.linspace(0.0, 2*np.pi, n_phi, endpoint=False)
    mats = []
    for phi in phis:
        S = kernel_matrix(N, phi=phi)
        mats.append(S.reshape(-1))
    A = np.stack(mats, axis=1)  # (N*N, n_phi)
    U, svals, Vt = np.linalg.svd(A, full_matrices=False)
    # effective dimension: number of singular values > tol
    tol = 1e-6 * svals[0]
    rank = int(np.sum(svals > tol))
    return {'n_phi': n_phi, 'singular_values': svals.tolist(), 'effective_rank': rank}

# Task 5: RG landscape sweep
def rg_landscape(N_list, s_values=np.linspace(0.2, 3.0, 141)):
    out = {}
    for N in N_list:
        gmax = []
        for s in s_values:
            S = kernel_matrix(N, alpha=ALPHA*s, beta=BETA)
            vals, _ = eigh(S)
            gmax.append(np.sort(vals)[-1])
        ln_s = np.log(s_values)
        g = np.array(gmax)
        dg = np.gradient(g, ln_s)
        sign_changes = int(np.sum(np.diff(np.sign(dg))!=0))
        out[N] = {'s_values': s_values.tolist(), 'g': g.tolist(), 'dg_dlns': dg.tolist(), 'sign_changes': sign_changes}
    return out


def run_all(N_list, out_json=OUT_JSON, out_md=OUT_MD):
    meta = {'created': datetime.now(timezone.utc).isoformat(), 'script': os.path.basename(__file__)}
    results = {}
    # Task1 per representative N (use N=24)
    S24 = kernel_matrix(24)
    results['task1_algebraic_probe_N24'] = algebraic_probe(S24, top_k=4)
    # Task2: PR scaling
    results['task2_PR_scaling'] = PR_scaling(N_list, top_k=4)
    # Task3: defect probe for each N
    results['task3_defect_probe'] = {N: defect_probe(N, defect_width=1) for N in N_list}
    # Task4: generator reconstruction for N=24
    results['task4_generator_reconstruction_N24'] = generator_reconstruction(24, n_phi=13)
    # Task5: RG landscape sweep (coarser s grid for speed)
    results['task5_rg_landscape'] = rg_landscape(N_list, s_values=np.linspace(0.5,2.5,101))
    out = {'meta': meta, 'results': results}
    with open(out_json, 'w') as f:
        json.dump(out, f, indent=2)
    # write MD with plain-language conclusions
    with open(out_md, 'w') as f:
        f.write(f"# Raport 112 — Deep Nadsoliton Probe\nData: {meta['created']}\n\n")
        # Task1 summary
        t1 = results['task1_algebraic_probe_N24']
        frac1 = t1['closure']['fraction_below_1e-2']
        frac2 = t1['closure']['fraction_below_1e-1']
        f.write("## Task 1 — Algebraic probe (N=24)\n")
        f.write("Wnioski po chłopsku:\n")
        if frac1 > 0.5:
            f.write(" - Duża zgodność z zamknięciem na bazie projekcji top-4 (wiele komutatorów dobrze wyrażalnych w bazie) → silny sygnał algebraiczny (może przypominać Lie algebra).\n")
        elif frac2 > 0.5:
            f.write(" - Częściowa zgodność (większość komutatorów wyrażalna w bazie do rzędu 10%) → sugestia struktury algebraicznej, wymaga dalszych testów.\n")
        else:
            f.write(" - Słaba zgodność z prostym zamknięciem na bazie projekcji top-4 → prawdopodobnie potrzebna większa baza (więcej trybów) lub inna baza generatorów.\n")
        f.write(f"Metryka: fraction_below_1e-2 = {frac1:.3f}, fraction_below_1e-1 = {frac2:.3f}\n\n")
        # Task2 summary
        f.write("## Task 2 — Participation Ratio scaling\n")
        pr = results['task2_PR_scaling']
        f.write("Wnioski po chłopsku:\n")
        for mode,info in pr['fits'].items():
            alpha = info['alpha']
            f.write(f" - {mode}: PR ~ N^{alpha:.3f}  (jeśli alpha≈1 → rozciągły tryb; alpha≈0 → lokalny)\n")
        f.write("Szczegóły: patrz sekcja danych w JSON\n\n")
        # Task3 summary
        f.write("## Task 3 — Topological defect probe\n")
        for N,val in results['task3_defect_probe'].items():
            f.write(f" - N={N}: PR_orig={val['PR_orig']:.3f}, PR_def={val['PR_def']:.3f}, delta_lambda={val['delta_lambda']:.6f}\n")
        f.write("Wnioski: jeśli PR_def << PR_orig i delta_lambda znaczące → lokalizowane jądro topologiczne.\n\n")
        # Task4 summary
        f.write("## Task 4 — Generator reconstruction (N=24)\n")
        gr = results['task4_generator_reconstruction_N24']
        f.write(f" - Effective rank (liczba niezależnych generatrów w przestrzeni S(phi)) = {gr['effective_rank']}\n")
        f.write("Wnioski: mała liczba generatorów sugeruje niskowymiarową algebrę operatorów; duża liczba → bardziej złożona struktura.\n\n")
        # Task5 summary
        f.write("## Task 5 — RG landscape sweep\n")
        for N,val in results['task5_rg_landscape'].items():
            f.write(f" - N={N}: sign_changes in dg/dlns = {val['sign_changes']}\n")
        f.write("Wnioski: regiony z change sign wskazują na przejścia screening/antiscreening; jeśli istnieją dla pewnych N→faza.\n\n")
        f.write("---\nPliki wygenerowane: JSON oraz MD.\n")
    return out

if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--Ns', nargs='+', type=int, default=[12,16,20,24,28,32])
    args = p.parse_args()
    report = run_all(args.Ns)
    print('Zapisano raport:', OUT_JSON)
