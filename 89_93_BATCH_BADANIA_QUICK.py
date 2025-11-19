#!/usr/bin/env python3
# Author: Krzysztof Żuchowski

"""
Batch runner for studies 89-93 (quick prototypes).

Creates `report_89_93_batch.md` containing five short investigations:
 - 89: simple feedback-coupling ODE stability and time simulation
 - 90: topological winding number on a 1D ring
 - 91: algebraic/spectral symmetry search (near-degenerate eigenspaces)
 - 92: topological mass-hierarchy probe via Laplacian defect localization
 - 93: time-dependent coupled-oscillator simulation with kernel coupling

No fitting. Outputs a Markdown report and prints its path.
"""
from __future__ import annotations

import io
import os
import sys
import math
import datetime
from contextlib import redirect_stdout

import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import eigh

# plotting (optional)
try:
    import matplotlib.pyplot as plt
    MATPLOTLIB_AVAILABLE = True
except Exception:
    MATPLOTLIB_AVAILABLE = False

try:
    from sklearn.cluster import KMeans
    SKLEARN_AVAILABLE = True
except Exception:
    SKLEARN_AVAILABLE = False


def study_89_feedback_coupling():
    """Prototype feedback model: dA/dt = (-gamma + alpha)*A - beta*A^3 + feedback*tanh(A).
    We'll scan feedback strength and simulate short dynamics.
    """
    out = io.StringIO()
    with redirect_stdout(out):
        print("### Badanie 89 — Feedback coupling (ODE)\\n")
        gamma = 1.0
        alpha = 0.5
        beta = 1.0
        feedback_values = [0.0, 0.5, 1.0, 1.5]

        def rhs(t, A, fb):
            return (-gamma + alpha) * A - beta * (A ** 3) + fb * math.tanh(A)

        results = {}
        sols = {}
        t_eval = np.linspace(0, 200, 1201)
        for fb in feedback_values:
            sol = solve_ivp(lambda t, y: rhs(t, y, fb), [0, 200], [0.01], t_eval=t_eval, max_step=0.5)
            sols[fb] = sol
            A_final = float(sol.y[0, -1])
            # crude growth-rate estimate from last half
            if sol.t.size > 10:
                t1 = int(sol.t.size // 2)
                y1 = sol.y[0, t1:]
                slopes = np.diff(y1) / np.diff(sol.t[t1:])
                est_rate = float(np.mean(slopes))
            else:
                est_rate = 0.0
            results[fb] = {"A_final": A_final, "est_rate": est_rate}
            print(f"feedback={fb:.3g}: A_final={A_final:.6g}, est_rate={est_rate:.6g}")

        # conclusion
        succ_fb = [fb for fb, v in results.items() if (v["A_final"] > 1e-2 or v["est_rate"] > 1e-4)]
        if succ_fb:
            conclusion = f"Self-excitation observed for feedback values: {succ_fb}."
        else:
            conclusion = "No clear self-excitation in tested feedback range."
        print("\\nConclusion:", conclusion)

        plots = []
        if MATPLOTLIB_AVAILABLE:
            os.makedirs("report_images", exist_ok=True)
            fig, ax = plt.subplots(figsize=(6, 3))
            for fb in feedback_values:
                sol = sols[fb]
                ax.plot(sol.t, sol.y[0], label=f"fb={fb}")
            ax.set_xlabel("t")
            ax.set_ylabel("A(t)")
            ax.legend(fontsize=8)
            imgp = os.path.join("report_images", "89_feedback_timeseries.png")
            fig.tight_layout()
            fig.savefig(imgp, dpi=150)
            plt.close(fig)
            plots.append(imgp)
            print(f"Plot saved: {imgp}")
        else:
            print("matplotlib not available; skipping plots for study 89.")

    status = "Sukces" if succ_fb else "Porażka"
    meta = {"params": {"gamma": gamma, "alpha": alpha, "beta": beta}, "scan": results, "conclusion": conclusion, "plots": plots, "status": status}
    return out.getvalue(), meta


def study_90_topological_winding(N=64, q_values=(0, 1, 2, 3)):
    """Compute discrete winding number (phase total / 2pi) for rings with different charge q.
    """
    out = io.StringIO()
    with redirect_stdout(out):
        print("### Badanie 90 — Topological winding on 1D ring\n")
        results = {}
        x = np.arange(N)
        plots = []
        for q in q_values:
            phase = (2 * np.pi * q * x / N) + 0.2 * np.random.randn(N) * 0.01
            up = np.unwrap(phase)
            winding = (up[-1] - up[0]) / (2 * np.pi)
            results[q] = {"winding_est": float(winding), "total_phase": float(up[-1] - up[0])}
            print(f"q={q}: winding_est={winding:.6g}")

            if MATPLOTLIB_AVAILABLE:
                os.makedirs("report_images", exist_ok=True)
                fig, ax = plt.subplots(figsize=(6, 2))
                ax.plot(x, np.unwrap(phase), marker=".", ls="-")
                ax.set_title(f"phase (q={q})")
                ax.set_xlabel("site")
                imgp = os.path.join("report_images", f"90_phase_q{q}.png")
                fig.tight_layout()
                fig.savefig(imgp, dpi=150)
                plt.close(fig)
                plots.append(imgp)

        # concise conclusion: check integer closeness
        integer_like = {q: abs(results[q]["winding_est"] - round(results[q]["winding_est"])) < 0.1 for q in q_values}
        concl = f"Winding integer-like: {integer_like}"
        print("\nConclusion:", concl)

    status = "Sukces" if all(integer_like.values()) else "Porażka"
    meta = {"N": N, "results": results, "conclusion": concl, "plots": plots, "status": status}
    return out.getvalue(), meta


def study_91_algebraic_symmetry_search(N=30, rng_seed=42):
    """Search for near-degenerate eigenvalues and low-dimensional invariant subspaces.
    Returns counts of near-degenerate pairs and kmeans clusters on eigenvectors.
    """
    out = io.StringIO()
    with redirect_stdout(out):
        print("### Badanie 91 — Algebraic / spectral symmetry search\n")
        rng = np.random.default_rng(rng_seed)
        # build a random symmetric (Hermitian) coupling matrix with some structured block noise
        A = rng.normal(size=(N, N))
        A = 0.5 * (A + A.T)
        # add a small block to induce near-degeneracy
        for i in range(0, N, max(1, N // 5)):
            j = min(N, i + max(2, N // 10))
            A[i:j, i:j] += 0.2 * np.ones((j - i, j - i))

        evals, evecs = eigh(A)
        # find near-degenerate eigenvalues
        diffs = np.diff(evals)
        tol = 1e-3 * max(1.0, np.abs(evals).max())
        near_pairs = np.where(np.abs(diffs) < tol)[0]
        print(f"N={N}, found {len(near_pairs)} near-degenerate adjacent eigenpairs (tol={tol:.3e}).")

        cluster_summary = None
        plots = []
        if SKLEARN_AVAILABLE:
            try:
                # cluster eigenvectors in small dimension (use first few components)
                k = min(4, N)
                km = KMeans(n_clusters=k, random_state=0, n_init=5)
                # use absolute values of top-5 eigenvectors as features
                features = np.abs(evecs[:, -min(6, N):])
                labels = km.fit_predict(features)
                counts = {int(i): int(np.sum(labels == i)) for i in range(k)}
                cluster_summary = counts
                print(f"KMeans on eigenvector features -> cluster sizes: {counts}")
                if MATPLOTLIB_AVAILABLE:
                    os.makedirs("report_images", exist_ok=True)
                    fig, ax = plt.subplots(figsize=(5, 2.5))
                    ax.plot(np.arange(len(evals)), evals, ".-")
                    ax.set_xlabel("index")
                    ax.set_ylabel("eigenvalue")
                    imgp = os.path.join("report_images", "91_eigenvalues.png")
                    fig.tight_layout()
                    fig.savefig(imgp, dpi=150)
                    plt.close(fig)
                    plots.append(imgp)
            except Exception as e:
                print("KMeans clustering failed:", e)
        else:
            print("scikit-learn not available; skipping clustering.")

    concl = f"near_pairs_count={len(near_pairs)}"
    status = "Sukces" if len(near_pairs) > 0 else "Porażka"
    meta = {"N": N, "near_pairs_count": int(len(near_pairs)), "clusters": cluster_summary, "conclusion": concl, "plots": plots, "status": status}
    return out.getvalue(), meta


def study_92_laplacian_defect_localization(N=80, defect_strength=5.0):
    """Create ring Laplacian with a defect (enhanced coupling at one link) and compute eigenmodes and IPR.
    Assess whether strong localization (high IPR) creates a separation in eigenvalues (a gap).
    """
    out = io.StringIO()
    with redirect_stdout(out):
        print("### Badanie 92 — Laplacian defect and localization\n")
        # ring Laplacian
        L = np.zeros((N, N))
        for i in range(N):
            L[i, i] = 2.0
            L[i, (i - 1) % N] = -1.0
            L[i, (i + 1) % N] = -1.0
        # enhance a link between node 0 and 1 (defect) by changing off-diagonals
        L[0, 1] = L[1, 0] = -defect_strength
        L[0, 0] += (defect_strength - 1.0)
        L[1, 1] += (defect_strength - 1.0)

        evals, evecs = eigh(L)
        # inverse participation ratio for each eigenvector
        ipr = np.sum(np.abs(evecs) ** 4, axis=0)
        # find most localized mode
        idx = int(np.argmax(ipr))
        gap = float(evals[1] - evals[0]) if N > 1 else 0.0
        print(f"min eval={evals[0]:.6g}, second={evals[1]:.6g} (gap={gap:.6g}), max IPR={ipr.max():.6g} at idx={idx}")

        plots = []
        if MATPLOTLIB_AVAILABLE:
            os.makedirs("report_images", exist_ok=True)
            fig, ax = plt.subplots(figsize=(5, 2.5))
            ax.plot(np.arange(len(ipr)), ipr, ".-")
            ax.set_xlabel("mode index")
            ax.set_ylabel("IPR")
            imgp = os.path.join("report_images", "92_ipr.png")
            fig.tight_layout()
            fig.savefig(imgp, dpi=150)
            plt.close(fig)
            plots.append(imgp)

    concl = ("localized mode found" if float(ipr.max()) > 1.0 / N * 5 else "no strong localization")
    status = "Sukces" if float(ipr.max()) > 1.0 / N * 5 else "Porażka"
    meta = {"N": N, "defect_strength": defect_strength, "min_eval": float(evals[0]), "gap": gap, "max_ipr": float(ipr.max()), "conclusion": concl, "plots": plots, "status": status}
    return out.getvalue(), meta


def kernel_K(d, alpha_geo=0.8, beta_tors=0.1, omega=2.0, phi=0.3):
    return alpha_geo * np.cos(omega * d + phi) / (1.0 + beta_tors * d)


def study_93_time_dependent_network(N=24, sim_t=100.0):
    """Simulate linear-damped + nonlinear coupling dynamics on N nodes with kernel-based coupling.
    dX/dt = -gamma X + C @ tanh(X)
    """
    out = io.StringIO()
    with redirect_stdout(out):
        print("### Badanie 93 — Time-dependent coupled network simulation\n")
        gamma = 0.5
        # build coupling matrix from kernel
        idx = np.arange(N)
        C = np.zeros((N, N))
        for i in range(N):
            for j in range(N):
                d = min(abs(i - j), N - abs(i - j))
                C[i, j] = kernel_K(d)
        # normalize coupling
        C = C / (np.max(np.abs(C)) + 1e-12)

        def rhs(t, X):
            return -gamma * X + C.dot(np.tanh(X))

        X0 = 0.01 * np.random.randn(N)
        sol = solve_ivp(rhs, [0, sim_t], X0, max_step=0.5)
        X_final = sol.y[:, -1]
        amp = np.abs(X_final)
        mean_amp = float(np.mean(amp))
        std_amp = float(np.std(amp))
        # locking measure: fraction of nodes within 10% of mean
        locked = float(np.sum(np.abs(amp - mean_amp) < 0.1 * (mean_amp + 1e-12)) / N)
        print(f"sim_t={sim_t}, mean_amp={mean_amp:.6g}, std_amp={std_amp:.6g}, locked_fraction={locked:.3g}")

        plots = []
        if MATPLOTLIB_AVAILABLE:
            os.makedirs("report_images", exist_ok=True)
            fig, ax = plt.subplots(figsize=(6, 3))
            # plot a subset of node time-series
            for i in range(min(6, N)):
                ax.plot(sol.t, sol.y[i, :], label=f"n{i}")
            ax.set_xlabel("t")
            ax.set_ylabel("X")
            ax.legend(fontsize=7)
            imgp = os.path.join("report_images", "93_timeseries.png")
            fig.tight_layout()
            fig.savefig(imgp, dpi=150)
            plt.close(fig)
            plots.append(imgp)

    concl = ("partial locking" if locked > 0.4 else "no strong locking observed")
    meta = {"N": N, "sim_t": sim_t, "mean_amp": mean_amp, "std_amp": std_amp, "locked_fraction": locked, "conclusion": concl, "plots": plots}
    return out.getvalue(), meta


def write_markdown_report(results_map, filename="report_89_93_batch.md"):
    now = datetime.datetime.now().astimezone().isoformat()
    lines = []
    lines.append(f"# Raport zbiorczy: Badania 89–93\n\nGenerated: {now}\n")
    for key, (text, meta) in results_map.items():
        lines.append(text)
        # add conclusion
        if meta and isinstance(meta, dict) and meta.get("conclusion"):
            lines.append(f"**Wnioski:** {meta.get('conclusion')}\n")

        # embed plots if present
        plots = []
        if meta and isinstance(meta, dict):
            plots = meta.get("plots", []) or []
        for p in plots:
            rel = os.path.relpath(p, os.path.dirname(filename))
            lines.append(f"![]({rel})\n")

        lines.append("\n---\n")
        # include a short machine-friendly summary table
        lines.append("```")
        for k, v in (meta or {}).items():
            # avoid dumping large numpy arrays
            try:
                lines.append(f"{k}: {v}")
            except Exception:
                lines.append(f"{k}: (unrepresentable)")
        lines.append("```")
        lines.append("\n")

    content = "\n".join(lines)
    with open(filename, "w", encoding="utf-8") as f:
        f.write(content)
    return os.path.abspath(filename)


def run_all_and_report():
    runners = [
        ("89_feedback", study_89_feedback_coupling),
        ("90_topology", study_90_topological_winding),
        ("91_algebraic", study_91_algebraic_symmetry_search),
        ("92_laplacian", study_92_laplacian_defect_localization),
        ("93_dynamics", study_93_time_dependent_network),
    ]
    results = {}
    for name, fn in runners:
        try:
            text, meta = fn()
        except Exception as e:
            text = f"### {name} failed with exception: {e}\n"
            meta = {"error": str(e)}
        results[name] = (text, meta)

    out_path = write_markdown_report(results, filename="report_89_93_batch.md")
    print("✅ Raport zapisany do:", out_path)
    return results, out_path


if __name__ == "__main__":
    run_all_and_report()
