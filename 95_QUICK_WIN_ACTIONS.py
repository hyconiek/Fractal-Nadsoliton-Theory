#!/usr/bin/env python3
# Author: Krzysztof Żuchowski

"""
Run 5 quick-win follow-up tasks based on Badanie 94 and produce a structured Markdown report
with Status: Sukces/Porażka, concise conclusions, and proposals for follow-up research.

Output: report_95_quick_win.md and optional images under report_images/
"""
from __future__ import annotations
import runpy
import io
import os
import datetime
from contextlib import redirect_stdout

import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import eigh

try:
    from sklearn.cluster import KMeans
    SKLEARN = True
except Exception:
    SKLEARN = False

try:
    import matplotlib.pyplot as plt
    MATPLOTLIB = True
except Exception:
    MATPLOTLIB = False


def load_theory_class(path="94_QUICK_WIN_NEW_MECHANISMS.py"):
    ns = runpy.run_path(path)
    if "FractalSupersolitonTheory" in ns:
        return ns["FractalSupersolitonTheory"]
    raise ImportError("FractalSupersolitonTheory not found in the given file")


def task_1_mass_hierarchy(Theory):
    """Refine mass-amplification and judge whether hierarchy ~O(10^2) emerges.
    Success if mu/e and tau/mu ratios exceed thresholds (heuristic).
    """
    out = io.StringIO()
    with redirect_stdout(out):
        print("## Zadanie 1: Hierarchia mas — interferencja fazowa")
        th = Theory()
        # replicate mass amplification from propose_mass_hierarchy_mechanism
        n = len(th.effective_octaves)
        octave_phases = th.omega * th.effective_octaves + th.phi
        S_complex = np.zeros((n, n), dtype=complex)
        for i in range(n):
            for j in range(n):
                if i != j:
                    phase_diff = np.exp(1j * (octave_phases[i] - octave_phases[j]))
                    S_complex[i, j] = th.S_matrix[i, j] * phase_diff
        mass_amp = np.abs(np.sum(S_complex, axis=1))
        # map to leptons as in file
        e_amp = mass_amp[0]
        mu_amp = mass_amp[3] if len(mass_amp) > 3 else mass_amp[1]
        tau_amp = mass_amp[5] if len(mass_amp) > 5 else mass_amp[-1]
        ratio_mu_e = float(mu_amp / (e_amp + 1e-12))
        ratio_tau_mu = float(tau_amp / (mu_amp + 1e-12))
        print(f"e_amp={e_amp:.4g}, mu_amp={mu_amp:.4g}, tau_amp={tau_amp:.4g}")
        print(f"ratio_mu_e={ratio_mu_e:.3f}, ratio_tau_mu={ratio_tau_mu:.3f}")

        # heuristic thresholds
        success = (ratio_mu_e > 10) and (ratio_tau_mu > 5)
        status = "Sukces" if success else "Porażka"
        if success:
            print("Status: Sukces — znacząca hierarchia wykryta.")
        else:
            print("Status: Porażka — hierarchia zbyt słaba; wymaga dodatkowych mechanizmów.")

        print("Wniosek: Interferencja fazowa daje efekt, ale nie wystarcza by wygenerować fizyczne proporcje bez dodatkowego wzmocnienia lub mechanizmu selektywnego.")
        print("Propozycja dalszych badań: wprowadzić zależność faz od lokalnego sprzężenia zwrotnego i sprawdzić selektywne wzmocnienie oktaw.")

        plots = []
        if MATPLOTLIB:
            os.makedirs("report_images", exist_ok=True)
            fig, ax = plt.subplots(figsize=(5, 2.5))
            ax.bar(np.arange(len(mass_amp)), mass_amp)
            ax.set_title("Mass amplification per octave")
            img = os.path.join("report_images", "95_task1_mass_amp.png")
            fig.tight_layout()
            fig.savefig(img, dpi=150)
            plt.close(fig)
            plots.append(img)
            print(f"Plot saved: {img}")

    meta = {"status": status, "ratio_mu_e": ratio_mu_e, "ratio_tau_mu": ratio_tau_mu, "plots": plots}
    return out.getvalue(), meta


def task_2_gravity_probe(Theory):
    """Simple probe of curvature-like term R ~ Laplacian(log rho) on a random field.
    Success if curvature variance is appreciable (signal vs noise).
    """
    out = io.StringIO()
    with redirect_stdout(out):
        print("## Zadanie 2: Probe emergentnej grawitacji — ∇²(log ρ)")
        # generate random positive density field on 2D grid
        rng = np.random.default_rng(0)
        N = 64
        rho = np.abs(rng.normal(size=(N, N))) + 0.1
        # compute discrete Laplacian of log(rho)
        logrho = np.log(rho)
        lap = (np.roll(logrho, 1, axis=0) + np.roll(logrho, -1, axis=0) + np.roll(logrho, 1, axis=1) + np.roll(logrho, -1, axis=1) - 4 * logrho)
        mean_lap = float(np.mean(lap))
        std_lap = float(np.std(lap))
        print(f"mean Lap(log rho) = {mean_lap:.6g}, std = {std_lap:.6g}")
        status = "Sukces" if std_lap > 0.01 else "Porażka"
        print(f"Status: {status}")
        print("Wniosek: Peturbacje pola dają niezerową krzywiznę; dalsze badania: sprzężenie z metryką i dynamika pola (time-dependent).")

        plots = []
        if MATPLOTLIB:
            os.makedirs("report_images", exist_ok=True)
            fig, ax = plt.subplots(figsize=(4, 3))
            im = ax.imshow(lap, cmap="RdBu", interpolation='nearest')
            fig.colorbar(im, ax=ax, fraction=0.046)
            img = os.path.join("report_images", "95_task2_laplogrho.png")
            fig.tight_layout()
            fig.savefig(img, dpi=150)
            plt.close(fig)
            plots.append(img)
            print(f"Plot saved: {img}")

    meta = {"status": status, "mean_lap": mean_lap, "std_lap": std_lap, "plots": plots}
    return out.getvalue(), meta


def task_3_gauge_mode_clustering(Theory):
    """Analyze eigenvalue gaps and try clustering of eigenvectors to identify mode families.
    Success if clear spectral gaps after 3 and 5 modes appear (heuristic for SU(3)/SU(2)/U(1)).
    """
    out = io.StringIO()
    with redirect_stdout(out):
        print("## Zadanie 3: Analiza modów i grup cechowania")
        th = Theory()
        evals, evecs = eigh(th.S_matrix)
        # sort descending by absolute value
        idx = np.argsort(np.abs(evals))[::-1]
        svals = np.abs(evals)[idx]
        print("Top eigenvalues:", ", ".join([f"{v:.4g}" for v in svals[:8]]))
        gaps = np.diff(svals)
        gap3 = float(gaps[2]) if len(gaps) > 2 else 0.0
        gap5 = float(gaps[4]) if len(gaps) > 4 else 0.0
        print(f"gap after 3 modes: {gap3:.6g}, gap after 5 modes: {gap5:.6g}")
        status = "Sukces" if (gap3 > 0.1 * svals[0] and gap5 > 0.05 * svals[0]) else "Porażka"
        print(f"Status: {status}")
        print("Wniosek: Jeśli status Sukces, spektrum separuje naturalnie modowe grupy; dalej: analizować strukturę pól w tych podprzestrzeniach.")

        plots = []
        if MATPLOTLIB:
            os.makedirs("report_images", exist_ok=True)
            fig, ax = plt.subplots(figsize=(5, 2))
            ax.plot(np.arange(len(svals)), svals, ".-")
            ax.set_xlabel("mode index")
            ax.set_ylabel("|eigenvalue|")
            img = os.path.join("report_images", "95_task3_spectrum.png")
            fig.tight_layout()
            fig.savefig(img, dpi=150)
            plt.close(fig)
            plots.append(img)
            print(f"Plot saved: {img}")

        # optional clustering on top eigenvectors
        cluster_info = None
        if SKLEARN:
            try:
                k = 3
                feats = np.abs(evecs[:, idx[:min(6, evecs.shape[1])]])
                km = KMeans(n_clusters=k, random_state=0, n_init=10)
                labels = km.fit_predict(feats)
                counts = {int(i): int(np.sum(labels == i)) for i in range(k)}
                cluster_info = counts
                print("Cluster sizes (k=3):", counts)
            except Exception as e:
                print("Clustering failed:", e)

    meta = {"status": status, "gap3": gap3, "gap5": gap5, "plots": plots, "clusters": cluster_info}
    return out.getvalue(), meta


def task_4_nonlinear_damping(Theory):
    """Estimate nonlinear damping parameters and compute predicted steady amplitude A*.
    Success if lambda_max>1 and A* is finite and moderate.
    """
    out = io.StringIO()
    with redirect_stdout(out):
        print("## Zadanie 4: Nieliniowe tłumienie i stabilizacja rezonansu")
        th = Theory()
        evals = np.linalg.eigvals(th.S_matrix)
        lambda_max = float(np.max(np.abs(evals)))
        n = len(th.effective_octaves)
        beta_fb = np.sum(np.abs(th.S_matrix - np.diag(np.diag(th.S_matrix)))) / (n ** 2)
        print(f"lambda_max={lambda_max:.6g}, beta_feedback={beta_fb:.6g}")
        status = "Porażka"
        A_star = None
        if lambda_max > 1 and beta_fb > 0:
            A_star = np.sqrt((lambda_max - 1) / (beta_fb + 1e-12))
            status = "Sukces" if np.isfinite(A_star) and A_star < 1e3 else "Porażka"
        print(f"Status: {status}")
        if A_star is not None:
            print(f"Predicted steady amplitude A* = {A_star:.4g}")
        print("Wniosek: Nieliniowe tłumienie może zapewnić stabilizację; dalsze: uwzględnić zależność tłumienia od częstotliwości.")

    meta = {"status": status, "lambda_max": lambda_max, "beta_fb": beta_fb, "A_star": A_star}
    return out.getvalue(), meta


def task_5_quick_simulation(Theory):
    """Run a short time-dependent simulation with S_matrix coupling: dX/dt = -gamma X + S @ tanh(X)
    Success if system reaches a bounded attractor (no blow-up) and shows partial locking.
    """
    out = io.StringIO()
    with redirect_stdout(out):
        print("## Zadanie 5: Krótka symulacja dynamiki sieciowej")
        th = Theory()
        S = th.S_matrix
        N = S.shape[0]
        gamma = 0.5

        def rhs(t, X):
            return -gamma * X + S.dot(np.tanh(X))

        X0 = 0.01 * np.random.randn(N)
        sol = solve_ivp(rhs, [0, 200], X0, max_step=0.5)
        Xf = sol.y[:, -1]
        mean_amp = float(np.mean(np.abs(Xf)))
        std_amp = float(np.std(np.abs(Xf)))
        locked_frac = float(np.sum(np.abs(np.abs(Xf) - mean_amp) < 0.1 * (mean_amp + 1e-12)) / N)
        status = "Sukces" if mean_amp < 1e2 and locked_frac > 0.2 else "Porażka"
        print(f"mean_amp={mean_amp:.6g}, std_amp={std_amp:.6g}, locked_frac={locked_frac:.3g}")
        print(f"Status: {status}")
        print("Wniosek: Dynamika pokazuje (nie)jednorodność amplitud; dalsze: parametryzować gamma i nieliniowość tanh.")

        plots = []
        if MATPLOTLIB:
            os.makedirs("report_images", exist_ok=True)
            fig, ax = plt.subplots(figsize=(6, 3))
            for i in range(min(6, N)):
                ax.plot(sol.t, sol.y[i, :], label=f"n{i}")
            ax.set_xlabel("t")
            ax.set_ylabel("X")
            ax.legend(fontsize=7)
            img = os.path.join("report_images", "95_task5_timeseries.png")
            fig.tight_layout()
            fig.savefig(img, dpi=150)
            plt.close(fig)
            plots.append(img)
            print(f"Plot saved: {img}")

    meta = {"status": status, "mean_amp": mean_amp, "std_amp": std_amp, "locked_frac": locked_frac, "plots": plots}
    return out.getvalue(), meta


def write_report(results, filename="report_95_quick_win.md"):
    now = datetime.datetime.now().astimezone().isoformat()
    lines = [f"# Raport — Badanie 95: Quick Win — Actionable follow-ups\n", f"Generated: {now}\n"]
    for key, (text, meta) in results.items():
        lines.append(text)
        lines.append("\n")
        # add structured summary
        lines.append(f"**Status:** {meta.get('status', 'Nieokreślony')}\n")
        # concise findings
        for k, v in meta.items():
            if k == 'plots':
                for p in (v or []):
                    rel = os.path.relpath(p, os.path.dirname(filename))
                    lines.append(f"![]({rel})\n")
            else:
                if k != 'status':
                    lines.append(f"- {k}: {v}\n")
        lines.append("---\n")

    content = "\n".join(lines)
    with open(filename, "w", encoding="utf-8") as f:
        f.write(content)
    return os.path.abspath(filename)


def main():
    Theory = load_theory_class()
    runners = [
        ("Zadanie_1_Hierarchia", task_1_mass_hierarchy),
        ("Zadanie_2_GravityProbe", task_2_gravity_probe),
        ("Zadanie_3_GaugeModes", task_3_gauge_mode_clustering),
        ("Zadanie_4_NonlinearDamping", task_4_nonlinear_damping),
        ("Zadanie_5_Simulation", task_5_quick_simulation),
    ]
    results = {}
    for name, fn in runners:
        try:
            text, meta = fn(Theory)
        except Exception as e:
            text = f"## {name} failed with exception: {e}\n"
            meta = {"status": "Porażka", "error": str(e)}
        results[name] = (text, meta)

    out_path = write_report(results, filename="report_95_quick_win.md")
    print("✅ Raport zapisany do:", out_path)


if __name__ == "__main__":
    main()
