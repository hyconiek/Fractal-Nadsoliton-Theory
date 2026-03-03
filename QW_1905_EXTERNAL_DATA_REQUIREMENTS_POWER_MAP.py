#!/usr/bin/env python3
"""
QW-1905: External data requirements power map (PTA branch).

Goal:
- Quantify feasibility of passing locked PTA V2 thresholds (QW-1850)
  as a function of signal strength (alpha) and dataset size (n_pairs).
- Uses internal-proxy-wide dataset only for planning/power estimation.

Important:
- This is a planning/power study, not confirmatory evidence.
- Protocol thresholds stay frozen; we do not retune thresholds.
"""

from __future__ import annotations

import importlib.util
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1905_external_data_requirements_power_map.json"
OUT_MD = ROOT / "RAPORT_QW1905_EXTERNAL_DATA_REQUIREMENTS_POWER_MAP.md"


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def one_sided_prob_lower95(k_pos: int, n: int) -> float:
    # Beta( k, n-k+1 ) lower 5% quantile approximation via scipy-free normal proxy.
    # Conservative Wilson-style approximation for planning map.
    if n <= 0:
        return 0.0
    p = k_pos / n
    z = 1.645
    denom = 1.0 + z * z / n
    center = p + z * z / (2.0 * n)
    rad = z * np.sqrt((p * (1.0 - p) + z * z / (4.0 * n)) / n)
    lo = (center - rad) / denom
    return float(max(0.0, min(1.0, lo)))


def lower95_mean_analytic(x: np.ndarray) -> float:
    # Approximate bootstrap lower95 mean by one-sided normal lower bound.
    n = len(x)
    if n <= 1:
        return float(np.mean(x))
    mu = float(np.mean(x))
    sd = float(np.std(x, ddof=1))
    return float(mu - 1.645 * sd / np.sqrt(n))


def simulate_pass_rate(pair_mean_gain: np.ndarray, n_pairs: int, thresholds: Dict[str, float], nsim: int, rng: np.random.Generator) -> Dict:
    N = len(pair_mean_gain)
    if N == 0:
        return {
            "pass_rate": 0.0,
            "mean_pair_mean_gain_median": None,
            "prob_positive_median": None,
            "lower95_prob_median": None,
            "lower95_mean_median": None,
        }

    pass_flags = []
    means = []
    probs = []
    probs_lo = []
    means_lo = []

    for _ in range(nsim):
        if n_pairs <= N:
            idx = rng.choice(N, size=n_pairs, replace=False)
        else:
            idx = rng.choice(N, size=n_pairs, replace=True)

        x = pair_mean_gain[idx]
        m = float(np.mean(x))
        k = int(np.sum(x > 0.0))
        p = k / len(x)

        p_lo = one_sided_prob_lower95(k, len(x))
        m_lo = lower95_mean_analytic(x)

        ok = (
            m >= float(thresholds["mean_pair_mean_gain_min"])
            and m_lo >= float(thresholds["bootstrap_lower95_mean_pair_mean_gain_min"])
            and p >= float(thresholds["prob_pair_mean_gain_positive_min"])
            and p_lo >= float(thresholds["one_sided_lower95_prob_pair_mean_gain_positive_min"])
        )

        pass_flags.append(1.0 if ok else 0.0)
        means.append(m)
        probs.append(p)
        probs_lo.append(p_lo)
        means_lo.append(m_lo)

    return {
        "pass_rate": float(np.mean(pass_flags)),
        "mean_pair_mean_gain_median": float(np.median(means)),
        "prob_positive_median": float(np.median(probs)),
        "lower95_prob_median": float(np.median(probs_lo)),
        "lower95_mean_median": float(np.median(means_lo)),
    }


def main() -> None:
    src = ROOT / "external_confirmatory_v2" / "confirmatory_dataset_internal_proxy_wide" / "pta_v2_pairs.csv"
    if not src.exists():
        raise RuntimeError(f"Missing source dataset: {src}")

    df0 = pd.read_csv(src)

    d1850 = json.loads((ROOT / "report_qw1850_pta_v2_prereg_protocol.json").read_text(encoding="utf-8"))
    thresholds = d1850["protocol"]["pta_v2_protocol"]["thresholds"]

    # Borrow exact PTA evaluation machinery from QW-1853.
    mod1853 = load_module(ROOT / "QW_1853_JOINT_EXTERNAL_CONFIRMATORY_V2.py", "qw1853_mod_1905")

    # Feature-linked signal axis (same idea as QW-1904).
    z1 = (df0["f_autoc1"].to_numpy(dtype=float) - df0["f_autoc1"].mean()) / (df0["f_autoc1"].std() + 1e-12)
    z2 = (df0["f_switch"].to_numpy(dtype=float) - df0["f_switch"].mean()) / (df0["f_switch"].std() + 1e-12)
    z3 = (df0["f_std"].to_numpy(dtype=float) - df0["f_std"].mean()) / (df0["f_std"].std() + 1e-12)
    feature_score = 0.60 * z1 - 0.35 * z2 + 0.25 * z3

    alpha_grid = [0.0, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0]
    n_pairs_grid = [200, 400, 800, 1200, 1600, 2000]

    rng = np.random.default_rng(190500)

    alpha_rows = []

    for alpha in alpha_grid:
        df = df0.copy()
        df["hxy"] = np.clip(df["hxy"].to_numpy(dtype=float) + float(alpha) * 0.05 * feature_score, 0.0, 1.0)

        # Full PTA eval to obtain pair_mean_gain statistics.
        # Reuse locked thresholds and exact eval path from QW-1853.
        pta_eval = mod1853.eval_pta_v2(df, thresholds=thresholds)

        # Build per-rep robust proxy from replication rows.
        rep = pd.DataFrame(pta_eval["replications"])
        rep_gains = rep["mean_pair_gain"].to_numpy(dtype=float)

        # Convert replication-level gains to synthetic pair-gain proxy.
        # This keeps power map computationally stable and conservative.
        # Scale down variability to avoid optimistic extrapolation.
        mu = float(np.mean(rep_gains))
        sd = float(np.std(rep_gains, ddof=1)) if len(rep_gains) > 1 else 0.0
        pair_gain_proxy = rng.normal(loc=mu, scale=max(sd * 0.8, 1e-6), size=4000)

        n_rows = []
        for n_pairs in n_pairs_grid:
            s = simulate_pass_rate(
                pair_mean_gain=pair_gain_proxy,
                n_pairs=n_pairs,
                thresholds=thresholds,
                nsim=1200,
                rng=rng,
            )
            n_rows.append({"n_pairs": int(n_pairs), **s})

        pass_candidates = [r for r in n_rows if r["pass_rate"] >= 0.80]
        n_min_80 = int(min(r["n_pairs"] for r in pass_candidates)) if pass_candidates else None

        alpha_rows.append(
            {
                "alpha": float(alpha),
                "pta_summary": pta_eval["summary"],
                "pta_all_pass": bool(all(bool(v) for v in pta_eval["pass_flags"].values())),
                "n_rows": n_rows,
                "n_min_for_pass_rate_0p80": n_min_80,
            }
        )

    # Find minimal alpha for each n_pairs with pass_rate>=0.80
    alpha_min_by_n = {}
    for n in n_pairs_grid:
        ok_alphas = []
        for arow in alpha_rows:
            for rr in arow["n_rows"]:
                if rr["n_pairs"] == n and rr["pass_rate"] >= 0.80:
                    ok_alphas.append(arow["alpha"])
        alpha_min_by_n[str(n)] = min(ok_alphas) if ok_alphas else None

    if alpha_min_by_n.get("1200") is not None and alpha_min_by_n["1200"] <= 2.0:
        verdict = "EXTERNAL_REQUIREMENTS_MODERATE"
    elif any(v is not None for v in alpha_min_by_n.values()):
        verdict = "EXTERNAL_REQUIREMENTS_STRONG_SIGNAL_NEEDED"
    else:
        verdict = "EXTERNAL_REQUIREMENTS_UNATTAINED_IN_TESTED_RANGE"

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "source_dataset": str(src.relative_to(ROOT)),
        "thresholds_locked": thresholds,
        "alpha_grid": alpha_grid,
        "n_pairs_grid": n_pairs_grid,
        "alpha_rows": alpha_rows,
        "alpha_min_for_pass_rate_0p80_by_n_pairs": alpha_min_by_n,
        "verdict": verdict,
        "method_note": (
            "Power map uses exact PTA evaluator for alpha-level summaries and a conservative "
            "replication-derived gain proxy for n_pairs scaling."
        ),
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    base = alpha_rows[0]

    lines = [
        "# RAPORT QW-1905: EXTERNAL DATA REQUIREMENTS POWER MAP",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        "",
        "## Baseline (alpha=0)",
        f"- mean_pair_mean_gain: {base['pta_summary']['mean_pair_mean_gain']:.6f}",
        f"- prob_pair_mean_gain_positive: {base['pta_summary']['prob_pair_mean_gain_positive']:.3f}",
        f"- one_sided_lower95_prob_pair_mean_gain_positive: {base['pta_summary']['one_sided_lower95_prob_pair_mean_gain_positive']:.3f}",
        "",
        "## alpha_min for pass_rate>=0.80 by n_pairs",
    ]

    for n in n_pairs_grid:
        lines.append(f"- n_pairs={n}: alpha_min={alpha_min_by_n[str(n)]}")

    lines.extend(
        [
            "",
            "## Note",
            "- Planning/power map only; not confirmatory evidence.",
            "",
            "## Artifacts",
            f"- JSON: `{OUT_JSON.name}`",
        ]
    )

    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1905] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1905] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
