#!/usr/bin/env python3
"""
QW-1778: PTA Bayes reanalysis with reparameterized angular model.

Models:
- M0 flat: C
- M1 legacy FIN: A * HD(theta) + C
- M2 reparam FIN: A * HD(theta)^q + C, q centered by QW-1773 p
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from itertools import combinations
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
from scipy.special import logsumexp
from scipy.stats import linregress

from astropy.coordinates import SkyCoord
import astropy.units as u


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1778_pta_bayes_reparam_reanalysis.json"
OUT_MD = ROOT / "RAPORT_QW1778_PTA_BAYES_REPARAM_REANALYSIS.md"


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


def load_residuals(res_dir: Path, max_psr: int = 30) -> Dict[str, pd.DataFrame]:
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


def hellings_downs(theta_deg: np.ndarray) -> np.ndarray:
    th = np.radians(theta_deg)
    x = (1.0 - np.cos(th)) / 2.0
    return 1.5 * x * np.log(x + 1e-12) - 0.25 * x + 0.5


def loglike(theta: np.ndarray, H: np.ndarray, sigma: float, model: np.ndarray) -> float:
    return float(-0.5 * np.sum(((H - model) / sigma) ** 2))


def evidence_flat(theta: np.ndarray, H: np.ndarray, sigma: float, n_mc: int = 22000, seed: int = 1778) -> float:
    rng = np.random.default_rng(seed)
    C = rng.uniform(-1.0, 2.0, n_mc)
    ll = np.array([loglike(theta, H, sigma, np.full_like(theta, c)) for c in C], dtype=float)
    return float(logsumexp(ll) - np.log(len(ll)))


def evidence_legacy(theta: np.ndarray, H: np.ndarray, sigma: float, n_mc: int = 26000, seed: int = 1779) -> float:
    rng = np.random.default_rng(seed)
    A = rng.uniform(-2.0, 2.0, n_mc)
    C = rng.uniform(-1.0, 2.0, n_mc)
    hd = hellings_downs(theta)
    ll = []
    for a, c in zip(A, C):
        model = a * hd + c
        ll.append(loglike(theta, H, sigma, model))
    ll = np.array(ll, dtype=float)
    return float(logsumexp(ll) - np.log(len(ll)))


def evidence_reparam(theta: np.ndarray, H: np.ndarray, sigma: float, q_center: float, n_mc: int = 32000, seed: int = 1780) -> float:
    rng = np.random.default_rng(seed)
    A = rng.uniform(-2.0, 2.0, n_mc)
    C = rng.uniform(-1.0, 2.0, n_mc)
    # Prior centered near q_center from 1773 with moderate spread.
    q = np.clip(rng.normal(loc=q_center, scale=0.18, size=n_mc), 0.8, 2.4)
    hd0 = np.clip(hellings_downs(theta), 1e-9, None)
    ll = []
    for a, c, qq in zip(A, C, q):
        model = a * (hd0 ** qq) + c
        ll.append(loglike(theta, H, sigma, model))
    ll = np.array(ll, dtype=float)
    return float(logsumexp(ll) - np.log(len(ll)))


def main() -> None:
    reparam = json.loads((ROOT / "report_qw1773_omega_suppressed_legacy_projection.json").read_text(encoding="utf-8"))
    q_center = float(reparam["projection"]["p"])

    res_dir = ROOT / "nano15/residuals/NANOGrav15yr_PulsarTiming_v2.1.0/residuals"
    par_dir = ROOT / "nano15/parfiles"

    residuals = load_residuals(res_dir, max_psr=30)
    positions = load_positions(par_dir)

    psr_list = list(residuals.keys())
    pairs = list(combinations(psr_list, 2))

    theta_vals = []
    hvals = []
    for p1, p2 in pairs:
        sep = angular_sep(p1, p2, positions)
        if sep is None:
            continue
        x, y = match_epochs(residuals[p1], residuals[p2], tol_days=30.0)
        if x is None:
            continue
        hxy = cross_dfa(x, y, min_scale=15)
        if np.isfinite(hxy):
            theta_vals.append(float(sep))
            hvals.append(float(hxy))

    theta = np.array(theta_vals, dtype=float)
    H = np.array(hvals, dtype=float)
    if len(H) < 80:
        raise RuntimeError("Too few valid pair points for Bayes reanalysis.")

    sigma = float(np.std(H))
    sigma = max(sigma, 1e-6)

    logZ_flat = evidence_flat(theta, H, sigma, n_mc=22000, seed=1778)
    logZ_legacy = evidence_legacy(theta, H, sigma, n_mc=26000, seed=1779)
    logZ_reparam = evidence_reparam(theta, H, sigma, q_center=q_center, n_mc=32000, seed=1780)

    logB_legacy = float(logZ_legacy - logZ_flat)
    logB_reparam = float(logZ_reparam - logZ_flat)
    delta_reparam_vs_legacy = float(logB_reparam - logB_legacy)

    pass_reparam_gain = delta_reparam_vs_legacy >= 0.35
    pass_reparam_nonneg = logB_reparam >= -0.05
    pass_pair_count = len(H) >= 120

    if pass_reparam_gain and pass_reparam_nonneg and pass_pair_count:
        verdict = "PTA_REPARAM_BAYES_IMPROVED_SUPPORTED"
    elif pass_reparam_gain and pass_pair_count:
        verdict = "PTA_REPARAM_BAYES_IMPROVED_PARTIAL"
    else:
        verdict = "PTA_REPARAM_BAYES_IMPROVED_WEAK"

    output = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "data_summary": {
            "n_pulsars_loaded": len(psr_list),
            "n_pairs_total": len(pairs),
            "n_pairs_valid": int(len(H)),
            "sigma_crossdfa": sigma,
        },
        "model_params": {"q_center_from_1773": q_center},
        "evidence": {
            "logZ_flat": logZ_flat,
            "logZ_legacy": logZ_legacy,
            "logZ_reparam": logZ_reparam,
            "logB_legacy_vs_flat": logB_legacy,
            "logB_reparam_vs_flat": logB_reparam,
            "delta_logB_reparam_minus_legacy": delta_reparam_vs_legacy,
        },
        "pass_flags": {
            "reparam_gain_over_legacy": bool(pass_reparam_gain),
            "reparam_nonnegative_evidence": bool(pass_reparam_nonneg),
            "enough_valid_pairs": bool(pass_pair_count),
        },
        "verdict": verdict,
    }
    OUT_JSON.write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1778: PTA BAYES REPARAM REANALYSIS",
        "",
        f"- Data UTC: {output['generated_utc']}",
        f"- Valid pairs: {len(H)} / {len(pairs)}",
        f"- logB legacy vs flat: {logB_legacy:.4f}",
        f"- logB reparam vs flat: {logB_reparam:.4f}",
        f"- delta logB (reparam-legacy): {delta_reparam_vs_legacy:.4f}",
        f"- Verdict: **{verdict}**",
        "",
        "## Pass Flags",
        f"- reparam_gain_over_legacy: {pass_reparam_gain}",
        f"- reparam_nonnegative_evidence: {pass_reparam_nonneg}",
        f"- enough_valid_pairs: {pass_pair_count}",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1778] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1778] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
