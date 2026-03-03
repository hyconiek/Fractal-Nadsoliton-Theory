#!/usr/bin/env python3
"""
QW-1788: q-prior width sweep under fixed high-coverage protocol.

Fixed protocol:
- operational cohort from QW-1781,
- stratified replications with fraction = 0.95 (QW-1787 recommendation).
"""

from __future__ import annotations

import importlib.util
import json
from itertools import combinations
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List

import numpy as np
from scipy.special import logsumexp


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1788_q_prior_width_sweep.json"
OUT_MD = ROOT / "RAPORT_QW1788_Q_PRIOR_WIDTH_SWEEP.md"


def load_helper():
    path = ROOT / "QW_1787_COVERAGE_FRACTION_SWEEP.py"
    spec = importlib.util.spec_from_file_location("qw1787_helper", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def evidence_reparam_width(helper, theta: np.ndarray, H: np.ndarray, sigma: float, q_center: float, q_width: float, n_mc: int, seed: int) -> float:
    rng = np.random.default_rng(seed)
    A = rng.uniform(-2.0, 2.0, n_mc)
    C = rng.uniform(-1.0, 2.0, n_mc)
    q = np.clip(rng.normal(loc=q_center, scale=q_width, size=n_mc), 0.8, 2.4)
    hd0 = np.clip(helper.hellings_downs(theta), 1e-9, None)
    ll = np.array(
        [helper.loglike(theta, H, sigma, a * (hd0 ** qq) + c) for a, c, qq in zip(A, C, q)],
        dtype=float,
    )
    return float(logsumexp(ll) - np.log(len(ll)))


def unit_from_scale(x: float, limit: float) -> float:
    return max(0.0, min(1.0, 1.0 - x / limit))


def main() -> None:
    helper = load_helper()

    d1773 = json.loads((ROOT / "report_qw1773_omega_suppressed_legacy_projection.json").read_text(encoding="utf-8"))
    d1781 = json.loads((ROOT / "report_qw1781_cohort_operational_gate.json").read_text(encoding="utf-8"))
    q_center = float(d1773["projection"]["p"])
    op = d1781["operational_cohort"]
    n_match_min = int(op["n_match_min"])
    stability_max = float(op["stability_max"])

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

    theta_all = np.array([r["theta_deg"] for r in rows], dtype=float)
    H_all = np.array([r["hxy"] for r in rows], dtype=float)
    frac = 0.95
    n_rep = 18
    q_widths = [0.10, 0.14, 0.16, 0.20, 0.24]

    results = []
    for j, q_width in enumerate(q_widths):
        sigma_all = max(float(np.std(H_all)), 1e-6)
        lz0 = helper.evidence_flat(theta_all, H_all, sigma_all, n_mc=9000, seed=4000 + 100 * j + 1)
        lz1 = helper.evidence_legacy(theta_all, H_all, sigma_all, n_mc=12000, seed=4000 + 100 * j + 2)
        lz2 = evidence_reparam_width(helper, theta_all, H_all, sigma_all, q_center=q_center, q_width=q_width, n_mc=15000, seed=4000 + 100 * j + 3)
        full_legacy = float(lz1 - lz0)
        full_reparam = float(lz2 - lz0)
        full_delta = float(full_reparam - full_legacy)

        rng = np.random.default_rng(17880 + j)
        rep_rows = []
        for i in range(n_rep):
            idx = helper.stratified_subsample_indices(theta_all, frac=frac, rng=rng, n_bins=8)
            if len(idx) < 80:
                continue
            th = theta_all[idx]
            hh = H_all[idx]
            sg = max(float(np.std(hh)), 1e-6)
            b0 = helper.evidence_flat(th, hh, sg, n_mc=5500, seed=5000 + 100 * j + 10 * i + 1)
            b1 = helper.evidence_legacy(th, hh, sg, n_mc=7000, seed=5000 + 100 * j + 10 * i + 2)
            b2 = evidence_reparam_width(
                helper, th, hh, sg, q_center=q_center, q_width=q_width, n_mc=9000, seed=5000 + 100 * j + 10 * i + 3
            )
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
        med_rep = float(np.median(arr_rep))
        med_dlt = float(np.median(arr_dlt))

        score = (
            0.20 * (1.0 if full_reparam > 0.0 else 0.0)
            + 0.20 * (1.0 if full_delta > 0.0 else 0.0)
            + 0.30 * prob_rep
            + 0.15 * prob_dlt
            + 0.10 * unit_from_scale(std_rep, 0.16)
            + 0.05 * unit_from_scale(std_dlt, 0.16)
        )

        results.append(
            {
                "q_width": q_width,
                "full_logB_legacy": full_legacy,
                "full_logB_reparam": full_reparam,
                "full_delta_logB": full_delta,
                "replications": len(rep_rows),
                "prob_logB_reparam_positive": prob_rep,
                "prob_delta_logB_positive": prob_dlt,
                "std_logB_reparam": std_rep,
                "std_delta_logB": std_dlt,
                "median_logB_reparam": med_rep,
                "median_delta_logB": med_dlt,
                "selection_score": float(score),
                "pass_basic": bool(full_reparam > 0.0 and full_delta > 0.0 and prob_rep >= 0.90 and prob_dlt >= 0.95),
            }
        )

    valid = [r for r in results if r["pass_basic"]]
    if len(valid) > 0:
        recommended = max(valid, key=lambda r: r["selection_score"])
        recommendation_strength = "STRONG"
    else:
        recommended = max(results, key=lambda r: r["selection_score"])
        recommendation_strength = "CONDITIONAL"

    recommended_width = float(recommended["q_width"])
    if recommendation_strength == "STRONG":
        verdict = "Q_PRIOR_WIDTH_SELECTION_SUPPORTED"
    else:
        verdict = "Q_PRIOR_WIDTH_SELECTION_PARTIAL"

    output = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "operational_criteria": {"n_match_min": n_match_min, "stability_max": stability_max},
        "protocol": {"fraction": frac, "replications": n_rep},
        "q_center": q_center,
        "q_widths_tested": q_widths,
        "results_by_width": results,
        "recommended_q_width": recommended_width,
        "recommendation_strength": recommendation_strength,
        "verdict": verdict,
    }
    OUT_JSON.write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1788: Q PRIOR WIDTH SWEEP",
        "",
        f"- Data UTC: {output['generated_utc']}",
        f"- Operational pairs: {len(rows)}",
        f"- Fixed fraction: {frac:.2f}",
        f"- Recommended q_width: {recommended_width:.3f} ({recommendation_strength})",
        f"- Verdict: **{verdict}**",
        "",
        "## Width Results",
    ]
    for r in results:
        lines.append(
            "- q_width={q_width:.3f} | score={selection_score:.3f} | full_reparam={full_logB_reparam:.3f} "
            "| full_delta={full_delta_logB:.3f} | P(reparam>0)={prob_logB_reparam_positive:.3f} "
            "| P(delta>0)={prob_delta_logB_positive:.3f} | std_reparam={std_logB_reparam:.3f} "
            "| std_delta={std_delta_logB:.3f} | pass_basic={pass_basic}".format(**r)
        )
    lines.extend(["", "## Artifacts", f"- JSON: `{OUT_JSON.name}`"])
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1788] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1788] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
