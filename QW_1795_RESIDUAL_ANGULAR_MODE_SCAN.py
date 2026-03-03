#!/usr/bin/env python3
"""
QW-1795: Residual angular-mode scan after locked-model fit.

Exploratory, post-lock-in diagnostic:
- scan Legendre modes P1..P6 on residuals from QW-1794 baseline fit,
- estimate which angular mode best explains remaining structure.
"""

from __future__ import annotations

import importlib.util
import json
from itertools import combinations
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List

import numpy as np
from numpy.polynomial import legendre as npleg


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1795_residual_angular_mode_scan.json"
OUT_MD = ROOT / "RAPORT_QW1795_RESIDUAL_ANGULAR_MODE_SCAN.md"


def load_helper():
    path = ROOT / "QW_1787_COVERAGE_FRACTION_SWEEP.py"
    spec = importlib.util.spec_from_file_location("qw1787_helper", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def legendre_mode(theta_deg: np.ndarray, ell: int) -> np.ndarray:
    x = np.cos(np.radians(theta_deg))
    coeff = np.zeros(ell + 1, dtype=float)
    coeff[ell] = 1.0
    return npleg.legval(x, coeff)


def best_linear_b(resid: np.ndarray, basis: np.ndarray) -> float:
    den = float(np.dot(basis, basis))
    if den <= 1e-12:
        return 0.0
    return float(np.dot(resid, basis) / den)


def kfold_indices(n: int, k: int, seed: int) -> List[np.ndarray]:
    rng = np.random.default_rng(seed)
    idx = np.arange(n)
    rng.shuffle(idx)
    return [idx[i::k] for i in range(k)]


def main() -> None:
    helper = load_helper()
    d1793 = json.loads((ROOT / "report_qw1793_model_lockin_gate.json").read_text(encoding="utf-8"))
    d1794 = json.loads((ROOT / "report_qw1794_locked_model_residual_audit.json").read_text(encoding="utf-8"))
    d1773 = json.loads((ROOT / "report_qw1773_omega_suppressed_legacy_projection.json").read_text(encoding="utf-8"))

    q_center = float(d1773["projection"]["p"])
    q_width = float(d1793["operational_protocol"]["q_width"])
    cohort = d1793["operational_protocol"]["cohort"]
    n_match_min = int(cohort["n_match_min"])
    stability_max = float(cohort["stability_max"])

    residuals = helper.load_residuals(ROOT / "nano15/residuals/NANOGrav15yr_PulsarTiming_v2.1.0/residuals", max_psr=34)
    positions = helper.load_positions(ROOT / "nano15/parfiles")

    rows: List[Dict[str, float]] = []
    psr_list = list(residuals.keys())
    for p1, p2 in combinations(psr_list, 2):
        sep = helper.angular_sep(p1, p2, positions)
        if sep is None:
            continue
        x, y = helper.match_epochs(residuals[p1], residuals[p2], tol_days=30.0)
        if x is None:
            continue
        hxy = helper.cross_dfa(x, y, min_scale=15)
        if not np.isfinite(hxy):
            continue
        stab = helper.split_half_stability(x, y)
        if len(x) >= n_match_min and stab <= stability_max:
            rows.append({"theta_deg": float(sep), "hxy": float(hxy)})

    if len(rows) < 85:
        raise RuntimeError("Operational cohort too small.")

    theta = np.array([r["theta_deg"] for r in rows], dtype=float)
    H = np.array([r["hxy"] for r in rows], dtype=float)
    n = len(theta)

    fit = d1794["fit"]
    A = float(fit["A"])
    q = float(fit["q"])
    C = float(fit["C"])
    hd0 = np.clip(helper.hellings_downs(theta), 1e-9, None)
    pred = A * (hd0 ** q) + C
    resid = H - pred

    sse_base = float(np.sum(resid * resid))
    bic_base = float(n * np.log(max(sse_base / n, 1e-12)) + 3.0 * np.log(n))

    mode_results = []
    folds = kfold_indices(n, k=6, seed=17950)
    for ell in range(1, 7):
        basis = legendre_mode(theta, ell)
        b_hat = best_linear_b(resid, basis)
        resid_new = resid - b_hat * basis
        sse_new = float(np.sum(resid_new * resid_new))
        bic_new = float(n * np.log(max(sse_new / n, 1e-12)) + 4.0 * np.log(n))
        delta_bic = float(bic_base - bic_new)
        r2_gain = float(max(0.0, 1.0 - sse_new / max(sse_base, 1e-12)))

        # Folded generalization of the extra mode.
        fold_gains = []
        for test_idx in folds:
            train_mask = np.ones(n, dtype=bool)
            train_mask[test_idx] = False
            b_train = best_linear_b(resid[train_mask], basis[train_mask])
            sse_test_base = float(np.sum(resid[test_idx] ** 2))
            sse_test_new = float(np.sum((resid[test_idx] - b_train * basis[test_idx]) ** 2))
            gain = 1.0 - sse_test_new / max(sse_test_base, 1e-12)
            fold_gains.append(float(gain))
        cv_gain_mean = float(np.mean(fold_gains))
        cv_gain_std = float(np.std(fold_gains))

        mode_results.append(
            {
                "ell": ell,
                "b_hat": b_hat,
                "delta_bic": delta_bic,
                "r2_gain_in_sample": r2_gain,
                "cv_gain_mean": cv_gain_mean,
                "cv_gain_std": cv_gain_std,
            }
        )

    best = max(mode_results, key=lambda r: (r["delta_bic"], r["cv_gain_mean"]))
    pass_bic = best["delta_bic"] >= 6.0
    pass_cv = best["cv_gain_mean"] > 0.02 and best["cv_gain_std"] <= 0.05

    if pass_bic and pass_cv:
        verdict = "RESIDUAL_MODE_CANDIDATE_IDENTIFIED"
    elif pass_bic or pass_cv:
        verdict = "RESIDUAL_MODE_SIGNAL_PARTIAL"
    else:
        verdict = "RESIDUAL_MODE_SIGNAL_WEAK"

    output = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "cohort_size": n,
        "baseline_fit_used": {"A": A, "q": q, "C": C, "q_center": q_center, "q_width": q_width},
        "baseline_sse": sse_base,
        "baseline_bic": bic_base,
        "mode_results": mode_results,
        "best_mode": best,
        "pass_flags": {
            "bic_support": bool(pass_bic),
            "cv_support": bool(pass_cv),
        },
        "verdict": verdict,
    }
    OUT_JSON.write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1795: RESIDUAL ANGULAR MODE SCAN",
        "",
        f"- Data UTC: {output['generated_utc']}",
        f"- Cohort size: {n}",
        f"- Baseline BIC: {bic_base:.4f}",
        f"- Best mode: ell={best['ell']}, delta_BIC={best['delta_bic']:.4f}, CV_gain={best['cv_gain_mean']:.4f}±{best['cv_gain_std']:.4f}",
        f"- Verdict: **{verdict}**",
        "",
        "## Pass Flags",
        f"- bic_support: {pass_bic}",
        f"- cv_support: {pass_cv}",
        "",
        "## Mode Table",
    ]
    for r in mode_results:
        lines.append(
            f"- ell={r['ell']}: delta_BIC={r['delta_bic']:.4f}, r2_gain={r['r2_gain_in_sample']:.4f}, "
            f"cv_gain={r['cv_gain_mean']:.4f}±{r['cv_gain_std']:.4f}, b_hat={r['b_hat']:.4f}"
        )
    lines.extend(["", "## Artifacts", f"- JSON: `{OUT_JSON.name}`"])
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1795] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1795] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
