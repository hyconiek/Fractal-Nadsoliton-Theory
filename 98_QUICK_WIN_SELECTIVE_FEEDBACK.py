#!/usr/bin/env python3
# Author: Krzysztof Żuchowski

"""
Selective-feedback quick-win.

Takes the top-N strongest off-diagonal couplings from the theory S_matrix and applies
multiplicative feedback to those links. Sweeps feedback strength and records mass amplification
and spectral response (lambda_max). Produces `report_98_quick_win.md` and plots in report_images/.

No fitting, no tautology — purely mechanistic probes.
"""
from __future__ import annotations
import runpy
import os
import io
import datetime
from contextlib import redirect_stdout
import numpy as np
from copy import deepcopy
from scipy.linalg import eigh
from scipy.integrate import solve_ivp

try:
    import matplotlib.pyplot as plt
    MATPLOTLIB = True
except Exception:
    MATPLOTLIB = False


def load_theory():
    ns = runpy.run_path('94_QUICK_WIN_NEW_MECHANISMS.py')
    return ns.get('FractalSupersolitonTheory')


def top_links(S, N=5):
    n = S.shape[0]
    links = []
    for i in range(n):
        for j in range(i+1, n):
            links.append(((i, j), abs(S[i, j])))
    links_sorted = sorted(links, key=lambda x: -x[1])
    return [t[0] for t in links_sorted[:N]]


def apply_feedback_to_S(S, links, feedback):
    S2 = S.copy()
    for (i, j) in links:
        S2[i, j] *= (1.0 + feedback)
        S2[j, i] *= (1.0 + feedback)
    return S2


def compute_mass_amplification(S):
    # crude proxy: row-sum magnitude
    return np.abs(np.sum(S, axis=1))


def run_sweep(Theory, links, feedbacks=None):
    if feedbacks is None:
        feedbacks = np.array([0.0, 0.2, 0.5, 1.0, 2.0])
    th = Theory()
    S_base = th.S_matrix
    n = S_base.shape[0]
    results = []
    for fb in feedbacks:
        S_fb = apply_feedback_to_S(S_base, links, fb)
        mass_amp = compute_mass_amplification(S_fb)
        evals, _ = eigh(S_fb)
        lambda_max = float(np.max(np.abs(evals)))
        results.append({'feedback': float(fb), 'mass_amp': mass_amp, 'lambda_max': lambda_max})
    return feedbacks, results


def quick_simulation(S, sim_t=100.0):
    # quick network sim to test dynamics
    n = S.shape[0]
    gamma = 0.5
    def rhs(t, X):
        return -gamma * X + S.dot(np.tanh(X))
    X0 = 0.01 * np.random.randn(n)
    sol = solve_ivp(rhs, [0, sim_t], X0, max_step=0.5)
    Xf = sol.y[:, -1]
    mean_amp = float(np.mean(np.abs(Xf)))
    locked_frac = float(np.sum(np.abs(np.abs(Xf) - mean_amp) < 0.1 * (mean_amp + 1e-12)) / n)
    return {'mean_amp': mean_amp, 'locked_frac': locked_frac}


def write_report(text_blocks, meta, filename='report_98_quick_win.md'):
    now = datetime.datetime.now().astimezone().isoformat()
    lines = [f'# Raport — Badanie 98: Selective Feedback Sweep\n', f'Generated: {now}\n']
    for tb in text_blocks:
        lines.append(tb)
        lines.append('\n')
    lines.append('## Meta summary\n')
    for k, v in meta.items():
        if k == 'plots':
            for p in v:
                rel = os.path.relpath(p, os.path.dirname(filename))
                lines.append(f'![]({rel})\n')
        else:
            lines.append(f'- {k}: {v}\n')
    with open(filename, 'w', encoding='utf-8') as f:
        f.write('\n'.join(lines))
    return os.path.abspath(filename)


def main():
    Theory = load_theory()
    th = Theory()
    S = th.S_matrix
    links = top_links(S, N=5)
    text_blocks = []
    # header
    header = f'## Selektywne sprzężenie na linkach: {links} (top 5)'
    text_blocks.append(header)
    feedbacks, results = run_sweep(Theory, links)
    text_blocks.append('### Sweep results (feedback -> lambda_max)')
    for r in results:
        text_blocks.append(f"feedback={r['feedback']}: lambda_max={r['lambda_max']:.6g}")
    # mass amplification table for final feedback
    final_mass = results[-1]['mass_amp']
    text_blocks.append('### Final mass amplification (last feedback)')
    text_blocks.append(', '.join([f'{v:.4g}' for v in final_mass]))

    # quick sim on S with strongest feedback applied
    S_fb = apply_feedback_to_S(S, links, float(feedbacks[-1]))
    sim_meta = quick_simulation(S_fb)
    text_blocks.append('### Quick dynamics after strongest feedback')
    text_blocks.append(f"mean_amp={sim_meta['mean_amp']:.6g}, locked_frac={sim_meta['locked_frac']:.3g}")

    plots = []
    if MATPLOTLIB:
        os.makedirs('report_images', exist_ok=True)
        # plot lambda_max vs feedback
        fig, ax = plt.subplots(figsize=(5, 2.5))
        ax.plot(feedbacks, [r['lambda_max'] for r in results], '.-')
        ax.set_xlabel('feedback')
        ax.set_ylabel('lambda_max')
        p1 = os.path.join('report_images', '98_lambda_vs_feedback.png')
        fig.tight_layout(); fig.savefig(p1, dpi=150); plt.close(fig)
        plots.append(p1)
        # heatmap of mass amplification (octaves x feedback)
        mass_mat = np.vstack([r['mass_amp'] for r in results])
        fig, ax = plt.subplots(figsize=(6, 3))
        im = ax.imshow(mass_mat.T, aspect='auto', cmap='viridis')
        ax.set_xlabel('feedback index')
        ax.set_ylabel('octave index')
        fig.colorbar(im, ax=ax, fraction=0.046)
        p2 = os.path.join('report_images', '98_massamp_heatmap.png')
        fig.tight_layout(); fig.savefig(p2, dpi=150); plt.close(fig)
        plots.append(p2)

    meta = {'links': links, 'feedbacks': list(map(float, feedbacks)), 'plots': plots, 'sim_meta': sim_meta}
    out_path = write_report(text_blocks, meta)
    print('✅ Raport zapisany do:', out_path)


if __name__ == '__main__':
    main()
