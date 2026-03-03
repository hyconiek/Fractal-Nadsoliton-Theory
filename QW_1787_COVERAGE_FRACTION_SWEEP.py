#!/usr/bin/env python3
"""
QW-1787: Coverage-fraction sweep for stratified replications.

Purpose:
- compare replication stability across subsample fractions,
- select operational fraction for next empirical campaign stage.
"""

from __future__ import annotations

import json
from itertools import combinations
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
from scipy.special import logsumexp
from scipy.stats import linregress

from astropy.coordinates import SkyCoord
import astropy.units as u


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1787_coverage_fraction_sweep.json"
OUT_MD = ROOT / "RAPORT_QW1787_COVERAGE_FRACTION_SWEEP.md"


def load_positions(par_dir: Path) -> Dict[str, SkyCoord]:
    pos: Dict[str, SkyCoord] = {}
    for par in par_dir.rglob("*.par"):
        name = par.stem.split("_")[0]
        ra = dec = elong = elat = None
        with par.open("r", encoding="utf-8", errors="ignore") as f:
            for line in f:
                if line.startswith("RAJ"):
                    ra = line.split()[1]
                elif line.startswith("DECJ"):
                    dec = line.split()[1]
                elif line.startswith("ELONG"):
                    elong = line.split()[1]
                elif line.startswith("ELAT"):
                    elat = line.split()[1]
        try:
            if ra and dec:
                pos[name] = SkyCoord(ra, dec, unit=(u.hourangle, u.deg), frame="icrs")
            elif elong and elat:
                pos[name] = SkyCoord(
                    lon=float(elong) * u.deg,
                    lat=float(elat) * u.deg,
                    frame="barycentrictrueecliptic",
                ).transform_to("icrs")
        except Exception:
            continue
    return pos


def angular_sep(psr1: str, psr2: str, pos: Dict[str, SkyCoord]) -> float | None:
    p1 = psr1.split("_")[0]
    p2 = psr2.split("_")[0]
    if p1 in pos and p2 in pos:
        return float(pos[p1].separation(pos[p2]).deg)
    return None


def load_residuals(res_dir: Path, max_psr: int = 34) -> Dict[str, pd.DataFrame]:
    residuals: Dict[str, pd.DataFrame] = {}
    files = sorted(res_dir.rglob("*.res"))[:max_psr]
    for f in files:
        psr = f.stem
        df = pd.read_csv(
            f,
            sep=r"\s+",
            comment="#",
            names=["MJD", "res_us", "err_us"],
            usecols=[0, 1, 2],
            engine="python",
        )
        df["res_s"] = df["res_us"] * 1e-6
        residuals[psr] = df.sort_values("MJD")
    return residuals


def match_epochs(df1: pd.DataFrame, df2: pd.DataFrame, tol_days: float = 30.0) -> Tuple[np.ndarray | None, np.ndarray | None]:
    t1, x1 = df1["MJD"].values, df1["res_s"].values
    t2, x2 = df2["MJD"].values, df2["res_s"].values
    matched_x: List[float] = []
    matched_y: List[float] = []
    j = 0
    for i in range(len(t1)):
        while j < len(t2) - 1 and t2[j] < t1[i] - tol_days:
            j += 1
        if abs(t2[j] - t1[i]) <= tol_days:
            matched_x.append(float(x1[i]))
            matched_y.append(float(x2[j]))
    if len(matched_x) < 60:
        return None, None
    return np.array(matched_x, dtype=float), np.array(matched_y, dtype=float)


def cross_dfa(x: np.ndarray, y: np.ndarray, min_scale: int = 15) -> float:
    n = min(len(x), len(y))
    x = x[:n] - np.mean(x[:n])
    y = y[:n] - np.mean(y[:n])
    X = np.cumsum(x)
    Y = np.cumsum(y)
    if len(X) < 4 * min_scale:
        return float("nan")
    scales = np.unique(np.logspace(np.log10(min_scale), np.log10(n // 4), 14).astype(int))
    F = []
    for s in scales:
        k = len(X) // s
        if k < 2:
            continue
        covs = []
        for i in range(k):
            xs = X[i * s : (i + 1) * s]
            ys = Y[i * s : (i + 1) * s]
            t = np.arange(s)
            px = np.polyfit(t, xs, 1)
            py = np.polyfit(t, ys, 1)
            cov = np.mean((xs - np.polyval(px, t)) * (ys - np.polyval(py, t)))
            covs.append(cov)
        m = np.mean(covs)
        if m == 0:
            continue
        F.append(np.sign(m) * np.sqrt(abs(m)))
    if len(F) < 4:
        return float("nan")
    F = np.array(F, dtype=float)
    scl = scales[: len(F)]
    valid = F != 0
    if np.sum(valid) < 4:
        return float("nan")
    slope, *_ = linregress(np.log(scl[valid]), np.log(np.abs(F[valid])))
    return float(slope * np.sign(np.mean(F)))


def split_half_stability(x: np.ndarray, y: np.ndarray) -> float:
    n = min(len(x), len(y))
    if n < 140:
        return float("inf")
    h = n // 2
    h1 = cross_dfa(x[:h], y[:h], min_scale=12)
    h2 = cross_dfa(x[h:], y[h:], min_scale=12)
    if not np.isfinite(h1) or not np.isfinite(h2):
        return float("inf")
    return float(abs(h1 - h2))


def hellings_downs(theta_deg: np.ndarray) -> np.ndarray:
    th = np.radians(theta_deg)
    x = (1.0 - np.cos(th)) / 2.0
    return 1.5 * x * np.log(x + 1e-12) - 0.25 * x + 0.5


def loglike(theta: np.ndarray, H: np.ndarray, sigma: float, model: np.ndarray) -> float:
    return float(-0.5 * np.sum(((H - model) / sigma) ** 2))


def evidence_flat(theta: np.ndarray, H: np.ndarray, sigma: float, n_mc: int, seed: int) -> float:
    rng = np.random.default_rng(seed)
    C = rng.uniform(-1.0, 2.0, n_mc)
    ll = np.array([loglike(theta, H, sigma, np.full_like(theta, c)) for c in C], dtype=float)
    return float(logsumexp(ll) - np.log(len(ll)))


def evidence_legacy(theta: np.ndarray, H: np.ndarray, sigma: float, n_mc: int, seed: int) -> float:
    rng = np.random.default_rng(seed)
    A = rng.uniform(-2.0, 2.0, n_mc)
    C = rng.uniform(-1.0, 2.0, n_mc)
    hd = hellings_downs(theta)
    ll = np.array([loglike(theta, H, sigma, a * hd + c) for a, c in zip(A, C)], dtype=float)
    return float(logsumexp(ll) - np.log(len(ll)))


def evidence_reparam(theta: np.ndarray, H: np.ndarray, sigma: float, q_center: float, n_mc: int, seed: int) -> float:
    rng = np.random.default_rng(seed)
    A = rng.uniform(-2.0, 2.0, n_mc)
    C = rng.uniform(-1.0, 2.0, n_mc)
    q = np.clip(rng.normal(loc=q_center, scale=0.16, size=n_mc), 0.8, 2.4)
    hd0 = np.clip(hellings_downs(theta), 1e-9, None)
    ll = np.array([loglike(theta, H, sigma, a * (hd0 ** qq) + c) for a, c, qq in zip(A, C, q)], dtype=float)
    return float(logsumexp(ll) - np.log(len(ll)))


def stratified_subsample_indices(theta: np.ndarray, frac: float, rng: np.random.Generator, n_bins: int = 8) -> np.ndarray:
    bins = np.linspace(0.0, 180.0, n_bins + 1)
    idx_all = np.arange(len(theta))
    out = []
    for i in range(n_bins):
        m = (theta >= bins[i]) & (theta < bins[i + 1] if i < n_bins - 1 else theta <= bins[i + 1])
        idx = idx_all[m]
        if len(idx) == 0:
            continue
        k = max(1, int(round(frac * len(idx))))
        k = min(k, len(idx))
        out.append(rng.choice(idx, size=k, replace=False))
    if len(out) == 0:
        return np.array([], dtype=int)
    return np.sort(np.concatenate(out))


def unit_from_scale(x: float, limit: float) -> float:
    return max(0.0, min(1.0, 1.0 - x / limit))


def main() -> None:
    d1773 = json.loads((ROOT / "report_qw1773_omega_suppressed_legacy_projection.json").read_text(encoding="utf-8"))
    d1781 = json.loads((ROOT / "report_qw1781_cohort_operational_gate.json").read_text(encoding="utf-8"))
    q_center = float(d1773["projection"]["p"])
    op = d1781["operational_cohort"]
    n_match_min = int(op["n_match_min"])
    stability_max = float(op["stability_max"])

    residuals = load_residuals(ROOT / "nano15/residuals/NANOGrav15yr_PulsarTiming_v2.1.0/residuals", max_psr=34)
    positions = load_positions(ROOT / "nano15/parfiles")

    rows = []
    psr_list = list(residuals.keys())
    for p1, p2 in combinations(psr_list, 2):
        sep = angular_sep(p1, p2, positions)
        if sep is None:
            continue
        x, y = match_epochs(residuals[p1], residuals[p2], tol_days=30.0)
        if x is None:
            continue
        hxy = cross_dfa(x, y, min_scale=15)
        if not np.isfinite(hxy):
            continue
        stab = split_half_stability(x, y)
        if len(x) >= n_match_min and stab <= stability_max:
            rows.append({"theta_deg": float(sep), "hxy": float(hxy)})

    if len(rows) < 85:
        raise RuntimeError("Operational cohort too small.")

    theta_all = np.array([r["theta_deg"] for r in rows], dtype=float)
    H_all = np.array([r["hxy"] for r in rows], dtype=float)

    fractions = [0.80, 0.85, 0.90, 0.95]
    n_rep = 18
    frac_results = []

    for j, frac in enumerate(fractions):
        rng = np.random.default_rng(17870 + j)
        rep_rows = []
        for i in range(n_rep):
            idx = stratified_subsample_indices(theta_all, frac=frac, rng=rng, n_bins=8)
            if len(idx) < 70:
                continue
            th = theta_all[idx]
            hh = H_all[idx]
            sg = max(float(np.std(hh)), 1e-6)
            b0 = evidence_flat(th, hh, sg, n_mc=6000, seed=3000 + 20 * j + 10 * i + 1)
            b1 = evidence_legacy(th, hh, sg, n_mc=8000, seed=3000 + 20 * j + 10 * i + 2)
            b2 = evidence_reparam(th, hh, sg, q_center=q_center, n_mc=10000, seed=3000 + 20 * j + 10 * i + 3)
            l1 = float(b1 - b0)
            l2 = float(b2 - b0)
            rep_rows.append(
                {
                    "rep": i,
                    "n_pairs": int(len(idx)),
                    "logB_legacy": l1,
                    "logB_reparam": l2,
                    "delta_logB": float(l2 - l1),
                }
            )

        arr_rep = np.array([r["logB_reparam"] for r in rep_rows], dtype=float)
        arr_dlt = np.array([r["delta_logB"] for r in rep_rows], dtype=float)
        prob_rep = float(np.mean(arr_rep > 0.0))
        prob_dlt = float(np.mean(arr_dlt > 0.0))
        std_rep = float(np.std(arr_rep))
        std_dlt = float(np.std(arr_dlt))
        median_rep = float(np.median(arr_rep))
        median_dlt = float(np.median(arr_dlt))

        score = (
            0.35 * prob_rep
            + 0.30 * prob_dlt
            + 0.20 * unit_from_scale(std_rep, 0.22)
            + 0.15 * unit_from_scale(std_dlt, 0.25)
        )

        frac_results.append(
            {
                "fraction": frac,
                "replications": len(rep_rows),
                "prob_logB_reparam_positive": prob_rep,
                "prob_delta_logB_positive": prob_dlt,
                "std_logB_reparam": std_rep,
                "std_delta_logB": std_dlt,
                "median_logB_reparam": median_rep,
                "median_delta_logB": median_dlt,
                "selection_score": float(score),
                "pass_basic": bool(prob_rep >= 0.80 and prob_dlt >= 0.95),
            }
        )

    valid = [r for r in frac_results if r["pass_basic"]]
    if len(valid) > 0:
        recommended = max(valid, key=lambda r: r["selection_score"])
    else:
        recommended = max(frac_results, key=lambda r: r["selection_score"])

    recommended_fraction = float(recommended["fraction"])
    recommendation_strength = "STRONG" if recommended["pass_basic"] else "CONDITIONAL"

    if recommendation_strength == "STRONG" and recommended_fraction >= 0.90:
        verdict = "COVERAGE_FRACTION_SELECTION_SUPPORTED"
    elif recommendation_strength == "STRONG":
        verdict = "COVERAGE_FRACTION_SELECTION_PARTIAL"
    else:
        verdict = "COVERAGE_FRACTION_SELECTION_WEAK"

    output = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "operational_criteria": {"n_match_min": n_match_min, "stability_max": stability_max},
        "n_pairs_operational": len(rows),
        "fractions_tested": fractions,
        "results_by_fraction": frac_results,
        "recommended_fraction": recommended_fraction,
        "recommendation_strength": recommendation_strength,
        "verdict": verdict,
    }
    OUT_JSON.write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1787: COVERAGE FRACTION SWEEP",
        "",
        f"- Data UTC: {output['generated_utc']}",
        f"- Operational pairs: {len(rows)}",
        f"- Recommended fraction: {recommended_fraction:.2f} ({recommendation_strength})",
        f"- Verdict: **{verdict}**",
        "",
        "## Fraction Results",
    ]
    for r in frac_results:
        lines.append(
            "- frac={fraction:.2f} | score={selection_score:.3f} | P(reparam>0)={prob_logB_reparam_positive:.3f} "
            "| P(delta>0)={prob_delta_logB_positive:.3f} | std_reparam={std_logB_reparam:.3f} "
            "| std_delta={std_delta_logB:.3f} | pass_basic={pass_basic}".format(**r)
        )
    lines.extend(["", "## Artifacts", f"- JSON: `{OUT_JSON.name}`"])
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1787] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1787] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
