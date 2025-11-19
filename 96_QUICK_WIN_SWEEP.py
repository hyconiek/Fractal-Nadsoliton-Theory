#!/usr/bin/env python3
# Author: Krzysztof Żuchowski

"""
Parameter sweep quick-win script based on Badanie 94/95.

Performs two sweeps:
 - Sweep A: vary an extra phase-feedback parameter to test mass-hierarchy amplification
 - Sweep B: vary damping gamma in network dynamics to search for locking regimes

Produces `report_96_quick_win.md` and plots in `report_images/`.
"""
from __future__ import annotations
import runpy
import os
import io
import datetime
from contextlib import redirect_stdout

import numpy as np
from scipy.integrate import solve_ivp

try:
    import matplotlib.pyplot as plt
    MATPLOTLIB = True
except Exception:
    MATPLOTLIB = False


def load_theory():
    ns = runpy.run_path("94_QUICK_WIN_NEW_MECHANISMS.py")
    return ns.get("FractalSupersolitonTheory")


def sweep_mass_phase(Theory, phase_range=None):
    if phase_range is None:
        phase_range = np.linspace(0.0, 2 * np.pi, 21)
    out = io.StringIO()
    with redirect_stdout(out):
        print("### Sweep A: faza sprzężenia -> wpływ na mass amplification")
        th = Theory()
        n = len(th.effective_octaves)
        ratios = []
        for phi_extra in phase_range:
            octave_phases = th.omega * th.effective_octaves + th.phi + phi_extra
            S_complex = np.zeros((n, n), dtype=complex)
            for i in range(n):
                for j in range(n):
                    if i != j:
                        phase_diff = np.exp(1j * (octave_phases[i] - octave_phases[j]))
                        S_complex[i, j] = th.S_matrix[i, j] * phase_diff
            mass_amp = np.abs(np.sum(S_complex, axis=1))
            e = mass_amp[0]
            mu = mass_amp[3] if len(mass_amp) > 3 else mass_amp[1]
            tau = mass_amp[5] if len(mass_amp) > 5 else mass_amp[-1]
            ratios.append((mu / (e + 1e-12), tau / (mu + 1e-12)))
        ratios = np.array(ratios)
        best_idx = np.argmax(ratios[:, 0])
        print(f"Best ratio mu/e = {ratios[best_idx,0]:.3f} at phi_extra={phase_range[best_idx]:.3f}")
    plots = []
    if MATPLOTLIB:
        os.makedirs("report_images", exist_ok=True)
        fig, ax = plt.subplots(figsize=(6, 3))
        ax.plot(phase_range, ratios[:, 0], label="mu/e")
        ax.plot(phase_range, ratios[:, 1], label="tau/mu")
        ax.set_xlabel("phi_extra")
        ax.set_ylabel("ratio")
        ax.legend()
        p = os.path.join("report_images", "96_sweep_mass_phase.png")
        fig.tight_layout(); fig.savefig(p, dpi=150); plt.close(fig)
        plots.append(p)
        print(f"Plot saved: {p}")

    meta = {"best_mu_e": float(ratios[best_idx, 0]), "best_phi": float(phase_range[best_idx]), "plots": plots}
    return out.getvalue(), meta


def sweep_gamma_dynamics(Theory, gammas=None):
    if gammas is None:
        gammas = np.linspace(0.1, 2.0, 20)
    out = io.StringIO()
    with redirect_stdout(out):
        print("### Sweep B: damping gamma -> network locking")
        th = Theory()
        S = th.S_matrix
        N = S.shape[0]
        locked = []
        for gamma in gammas:
            def rhs(t, X):
                return -gamma * X + S.dot(np.tanh(X))
            X0 = 0.01 * np.random.randn(N)
            sol = solve_ivp(rhs, [0, 200], X0, max_step=0.5)
            Xf = sol.y[:, -1]
            mean_amp = float(np.mean(np.abs(Xf)))
            locked_frac = float(np.sum(np.abs(np.abs(Xf) - mean_amp) < 0.1 * (mean_amp + 1e-12)) / N)
            locked.append(locked_frac)
        locked = np.array(locked)
        best_idx = int(np.argmax(locked))
        print(f"Best locked fraction {locked[best_idx]:.3f} at gamma={gammas[best_idx]:.3f}")
    plots = []
    if MATPLOTLIB:
        os.makedirs("report_images", exist_ok=True)
        fig, ax = plt.subplots(figsize=(6, 3))
        ax.plot(gammas, locked, ".-")
        ax.set_xlabel("gamma")
        ax.set_ylabel("locked_fraction")
        p = os.path.join("report_images", "96_sweep_gamma_locked.png")
        fig.tight_layout(); fig.savefig(p, dpi=150); plt.close(fig)
        plots.append(p)
        print(f"Plot saved: {p}")

    meta = {"best_locked": float(locked[best_idx]), "best_gamma": float(gammas[best_idx]), "plots": plots}
    return out.getvalue(), meta


def write_report(results, filename="report_96_quick_win.md"):
    now = datetime.datetime.now().astimezone().isoformat()
    lines = [f"# Raport — Badanie 96: Sweep Quick Win\n", f"Generated: {now}\n"]
    for k, (txt, meta) in results.items():
        lines.append(txt)
        lines.append("\n")
        lines.append(f"**Summary:**\n")
        for kk, vv in meta.items():
            if kk == 'plots':
                for p in (vv or []):
                    rel = os.path.relpath(p, os.path.dirname(filename))
                    lines.append(f"![]({rel})\n")
            else:
                lines.append(f"- {kk}: {vv}\n")
        lines.append("---\n")
    with open(filename, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))
    return os.path.abspath(filename)


def main():
    Theory = load_theory()
    results = {}
    a_txt, a_meta = sweep_mass_phase(Theory)
    results['SweepA_mass_phase'] = (a_txt, a_meta)
    b_txt, b_meta = sweep_gamma_dynamics(Theory)
    results['SweepB_gamma'] = (b_txt, b_meta)
    out_path = write_report(results)
    print("✅ Raport zapisany do:", out_path)


if __name__ == '__main__':
    main()
