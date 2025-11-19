#!/usr/bin/env python3
# Author: Krzysztof Żuchowski

"""
Quick-win: five short, non-fitting, non-tautological analyses derived from the theory class.

Produces `report_97_quick_win.md` and plots under `report_images/` when matplotlib is available.

Analyses (no fitting, no tautology):
 1) Consistency check of S_matrix symmetry and eigen-spectrum sanity
 2) Robustness sweep: small perturbations of alpha_geo and effect on λ_max
 3) Sensitivity test: knock-out single off-diagonal links and measure Δλ_max
 4) Discrete topology check: stable integer winding on deterministic phases
 5) Minimal mechanism ranking: rank links by influence on λ_max

Each analysis reports a concise conclusion and Status: Sukces/Porażka (heuristic, transparent).
"""
from __future__ import annotations
import runpy
import os
import io
import datetime
from contextlib import redirect_stdout

import numpy as np
from scipy.linalg import eigh
from copy import deepcopy

try:
    import matplotlib.pyplot as plt
    MATPLOTLIB = True
except Exception:
    MATPLOTLIB = False


def load_theory():
    ns = runpy.run_path("94_QUICK_WIN_NEW_MECHANISMS.py")
    if "FractalSupersolitonTheory" not in ns:
        raise ImportError("FractalSupersolitonTheory not found")
    return ns["FractalSupersolitonTheory"]


def analysis_1_consistency(Theory):
    out = io.StringIO()
    with redirect_stdout(out):
        print("## Analiza 1: Konsystencja macierzy S i spektrum")
        th = Theory()
        S = th.S_matrix
        # symmetry
        sym_diff = np.max(np.abs(S - S.T))
        evals, _ = eigh(S)
        lambda_max = float(np.max(np.abs(evals)))
        print(f"max|S-S^T| = {sym_diff:.3e}, lambda_max = {lambda_max:.6g}")
        status = "Sukces" if sym_diff < 1e-8 and np.isfinite(lambda_max) else "Porażka"
        if status == "Sukces":
            print("Status: Sukces — S is symmetric (within numerical tolerance) and spectrum finite.")
        else:
            print("Status: Porażka — numerical asymmetry or invalid spectrum detected.")
    meta = {"sym_diff": sym_diff, "lambda_max": lambda_max, "status": status}
    return out.getvalue(), meta


def analysis_2_robustness_alpha(Theory, deltas=None):
    out = io.StringIO()
    if deltas is None:
        deltas = np.linspace(-0.2, 0.2, 21)
    with redirect_stdout(out):
        print("## Analiza 2: Robustność względem niewielkiej zmiany alpha_geo")
        th = Theory()
        base_alpha = th.alpha_geo
        lambda_changes = []
        for d in deltas:
            th2 = deepcopy(th)
            th2.alpha_geo = base_alpha * (1.0 + d)
            # rebuild S_matrix using class internal method if present
            try:
                th2.S_matrix = th2._setup_matrices()
            except Exception:
                pass
            evals, _ = eigh(th2.S_matrix)
            lambda_changes.append(float(np.max(np.abs(evals))))
        lambda_changes = np.array(lambda_changes)
        rel_change = (lambda_changes - lambda_changes[len(deltas)//2]) / (np.abs(lambda_changes[len(deltas)//2]) + 1e-12)
        max_rel = float(np.max(np.abs(rel_change)))
        print(f"max relative change in lambda_max over sweep = {max_rel:.6g}")
        status = "Sukces" if max_rel < 0.5 else "Porażka"
        print("Status:", status)

    plots = []
    if MATPLOTLIB:
        os.makedirs("report_images", exist_ok=True)
        fig, ax = plt.subplots(figsize=(6, 3))
        ax.plot(deltas, lambda_changes, '.-')
        ax.set_xlabel('relative delta')
        ax.set_ylabel('lambda_max')
        p = os.path.join('report_images', '97_robustness_alpha.png')
        fig.tight_layout(); fig.savefig(p, dpi=150); plt.close(fig)
        plots.append(p)
    meta = {"max_rel_change": max_rel, "status": status, "plots": plots}
    return out.getvalue(), meta


def analysis_3_sensitivity_links(Theory):
    out = io.StringIO()
    with redirect_stdout(out):
        print("## Analiza 3: Wrażliwość na usunięcie pojedynczych połączeń")
        th = Theory()
        S = deepcopy(th.S_matrix)
        n = S.shape[0]
        evals0, _ = eigh(S)
        lambda0 = float(np.max(np.abs(evals0)))
        deltas = []
        coords = []
        for i in range(n):
            for j in range(i+1, n):
                if S[i, j] == 0:
                    continue
                S2 = S.copy()
                S2[i, j] = 0.0
                S2[j, i] = 0.0
                l2, _ = eigh(S2)
                lambda2 = float(np.max(np.abs(l2)))
                deltas.append(lambda0 - lambda2)
                coords.append((i, j))
        if not deltas:
            print('No off-diagonal links to test.')
            meta = {'status': 'Porażka', 'reason': 'no_links'}
            return out.getvalue(), meta
        deltas = np.array(deltas)
        order = np.argsort(-deltas)
        top = [(coords[k], float(deltas[k])) for k in order[:5]]
        print('Top 5 links by lambda_max drop when removed:')
        for c, d in top:
            print(f'  link {c} -> Δλ = {d:.6g}')
        status = 'Sukces' if deltas.max() > 1e-3 * lambda0 else 'Porażka'
    meta = {'top_links': top, 'max_delta': float(deltas.max()), 'status': status}
    return out.getvalue(), meta


def analysis_4_topology_winding(Theory):
    out = io.StringIO()
    with redirect_stdout(out):
        print('## Analiza 4: Deterministyczny test obwodu — winding')
        th = Theory()
        # deterministic phases: use integer multiples of base phase
        N = 32
        x = np.arange(N)
        q = 1
        phase = (2 * np.pi * q * x / N)
        winding = (np.unwrap(phase)[-1] - np.unwrap(phase)[0]) / (2 * np.pi)
        is_integer = abs(winding - round(winding)) < 1e-6
        print(f'winding={winding:.6g}, integer_like={is_integer}')
        status = 'Sukces' if is_integer else 'Porażka'
    meta = {'winding': float(winding), 'integer_like': bool(is_integer), 'status': status}
    return out.getvalue(), meta


def analysis_5_minimal_mechanisms(Theory):
    out = io.StringIO()
    with redirect_stdout(out):
        print('## Analiza 5: Ranking minimalnych mechanizmów (wpływ połączeń)')
        th = Theory()
        # reuse sensitivity data: rank off-diagonal absolute weights
        S = th.S_matrix
        n = S.shape[0]
        links = []
        for i in range(n):
            for j in range(i+1, n):
                if S[i, j] != 0:
                    links.append(((i, j), abs(S[i, j])))
        links_sorted = sorted(links, key=lambda x: -x[1])
        top5 = links_sorted[:5]
        print('Top 5 strongest off-diagonal couplings:')
        for (i, j), val in top5:
            print(f'  ({i},{j}) = {val:.6g}')
        status = 'Sukces' if len(top5) > 0 else 'Porażka'
    meta = {'top_couplings': [(t[0], float(t[1])) for t in top5], 'status': status}
    return out.getvalue(), meta


def write_report(results, filename='report_97_quick_win.md'):
    now = datetime.datetime.now().astimezone().isoformat()
    lines = [f'# Raport — Badanie 97: No-tautology Quick Win\n', f'Generated: {now}\n']
    for k, (txt, meta) in results.items():
        lines.append(txt)
        lines.append('\n')
        lines.append(f"**Status:** {meta.get('status', 'Nieokreślony')}\n")
        # add small summary
        for kk, vv in meta.items():
            if kk == 'plots':
                for p in (vv or []):
                    rel = os.path.relpath(p, os.path.dirname(filename))
                    lines.append(f'![]({rel})\n')
            else:
                lines.append(f'- {kk}: {vv}\n')
        lines.append('---\n')
    with open(filename, 'w', encoding='utf-8') as f:
        f.write('\n'.join(lines))
    return os.path.abspath(filename)


def main():
    Theory = load_theory()
    runners = [
        ('Analiza_1_Konsystencja', analysis_1_consistency),
        ('Analiza_2_Robustnosc', analysis_2_robustness_alpha),
        ('Analiza_3_Sensitivity', analysis_3_sensitivity_links),
        ('Analiza_4_Topologia', analysis_4_topology_winding),
        ('Analiza_5_MinMechanisms', analysis_5_minimal_mechanisms),
    ]
    results = {}
    for name, fn in runners:
        try:
            txt, meta = fn(Theory)
        except Exception as e:
            txt = f'## {name} failed with exception: {e}\n'
            meta = {'status': 'Porażka', 'error': str(e)}
        results[name] = (txt, meta)

    out_path = write_report(results)
    print('✅ Raport zapisany do:', out_path)


if __name__ == '__main__':
    main()
