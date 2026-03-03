#!/usr/bin/env python3
"""
QW-1736: Bayesian model selection for node patterns.

Models:
- M_A: nodes at d = 2 + 3n
- M_B: nodes at d = 2 + 6n
- M_C: nodes at d = 4/3 + 4n   (analytic from omega=pi/4, phi=pi/6)
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List

import numpy as np


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1736_kernel_node_bayesian_model_selection.json"
OUT_MD = ROOT / "RAPORT_QW1736_KERNEL_NODE_BAYESIAN_MODEL_SELECTION.md"


def load_rows_1735() -> List[Dict[str, object]]:
    p = ROOT / "report_qw1735_omega_phi_from_lattice_dynamics.json"
    d = json.loads(p.read_text(encoding="utf-8"))
    return d.get("rows", [])


def extract_observed_positions(rows: List[Dict[str, object]]) -> np.ndarray:
    obs: List[float] = []
    for r in rows:
        prof = r.get("profile", {})
        y = np.array([float(prof[str(i)]) for i in range(1, 13)], dtype=float)
        ay = np.abs(y)
        thr = float(np.quantile(ay, 0.25))
        for i in range(1, 13):
            if ay[i - 1] <= thr:
                obs.append(float(i))
        # Also include explicit local minima if present.
        mins = [int(x) for x in r.get("minima_abs_positions", [])]
        for m in mins:
            obs.append(float(m))
    return np.array(obs, dtype=float)


def predicted_sequence(model: str, n_terms: int = 12) -> np.ndarray:
    if model == "M_A_2_plus_3n":
        return np.array([2.0 + 3.0 * n for n in range(n_terms)], dtype=float)
    if model == "M_B_2_plus_6n":
        return np.array([2.0 + 6.0 * n for n in range(n_terms)], dtype=float)
    if model == "M_C_4over3_plus_4n":
        return np.array([4.0 / 3.0 + 4.0 * n for n in range(n_terms)], dtype=float)
    raise ValueError(model)


def log_likelihood(obs: np.ndarray, seq: np.ndarray, sigma: float, epsilon: float) -> float:
    seq_shift = seq + epsilon
    const = -math.log(sigma * math.sqrt(2.0 * math.pi))
    ll = 0.0
    for x in obs:
        r = float(np.min(np.abs(seq_shift - x)))
        ll += const - 0.5 * (r / sigma) ** 2
    return ll


def logsumexp(a: np.ndarray) -> float:
    m = float(np.max(a))
    return m + float(np.log(np.sum(np.exp(a - m))))


def model_evidence(obs: np.ndarray, model: str, sigma: float) -> Dict[str, float]:
    seq = predicted_sequence(model, n_terms=18)
    eps_grid = np.linspace(-1.0, 1.0, 201)
    lls = np.array([log_likelihood(obs, seq, sigma, float(e)) for e in eps_grid], dtype=float)
    # Uniform prior over epsilon => average likelihood.
    log_z = logsumexp(lls) - math.log(len(eps_grid))
    idx = int(np.argmax(lls))
    return {
        "log_evidence": float(log_z),
        "best_epsilon": float(eps_grid[idx]),
        "best_log_like": float(lls[idx]),
    }


def main() -> None:
    rows = load_rows_1735()
    obs = extract_observed_positions(rows)
    if len(obs) == 0:
        raise RuntimeError("No observed quasi-zero/minima positions extracted from QW-1735.")

    sigma = 0.70  # conservative uncertainty: integer d + quasi-zero extraction noise

    models = ["M_A_2_plus_3n", "M_B_2_plus_6n", "M_C_4over3_plus_4n"]
    evid: Dict[str, Dict[str, float]] = {}
    for m in models:
        evid[m] = model_evidence(obs, m, sigma=sigma)

    logz = np.array([evid[m]["log_evidence"] for m in models], dtype=float)
    logz_norm = logz - logsumexp(logz)
    post = np.exp(logz_norm)

    post_map = {m: float(post[i]) for i, m in enumerate(models)}
    best_idx = int(np.argmax(post))
    second_idx = int(np.argsort(post)[-2])
    best_model = models[best_idx]
    second_model = models[second_idx]
    bf_best_vs_second = float(np.exp(logz[best_idx] - logz[second_idx]))

    if post[best_idx] >= 0.90 and bf_best_vs_second >= 10.0:
        verdict = "NODE_MODEL_STRONGLY_SELECTED"
    elif post[best_idx] >= 0.70 and bf_best_vs_second >= 3.0:
        verdict = "NODE_MODEL_MODERATELY_SELECTED"
    else:
        verdict = "NODE_MODEL_NOT_DECISIVE"

    compatibility_with_characteristic_priors = best_model == "M_C_4over3_plus_4n"

    output = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "protocol": {
            "source_report": "report_qw1735_omega_phi_from_lattice_dynamics.json",
            "n_observations": int(len(obs)),
            "sigma_assumed": sigma,
            "epsilon_range": [-1.0, 1.0],
            "epsilon_grid_size": 201,
        },
        "observed_positions_summary": {
            "mean": float(np.mean(obs)),
            "std": float(np.std(obs)),
            "q10": float(np.quantile(obs, 0.1)),
            "q90": float(np.quantile(obs, 0.9)),
        },
        "model_evidences": evid,
        "posterior_probabilities": post_map,
        "best_model": best_model,
        "second_model": second_model,
        "bayes_factor_best_vs_second": bf_best_vs_second,
        "compatibility_with_characteristic_priors": compatibility_with_characteristic_priors,
        "verdict": verdict,
    }
    OUT_JSON.write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1736: BAYESIAN NODE MODEL SELECTION",
        "",
        f"- Data UTC: {output['generated_utc']}",
        f"- Observations: {len(obs)}",
        f"- Sigma: {sigma:.3f}",
        f"- Best model: {best_model}",
        f"- Posterior(best): {post_map[best_model]:.4f}",
        f"- Bayes factor best/second: {bf_best_vs_second:.3f}",
        f"- Compatible with priors (pi/4, pi/6): {compatibility_with_characteristic_priors}",
        f"- Verdict: **{verdict}**",
        "",
        "## Posterior",
    ]
    for m in models:
        lines.append(f"- {m}: {post_map[m]:.6f}")
    lines.extend(["", "## Artifacts", f"- JSON: `{OUT_JSON.name}`"])
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1736] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1736] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
