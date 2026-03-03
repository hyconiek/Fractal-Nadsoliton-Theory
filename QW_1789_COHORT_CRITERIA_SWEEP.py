#!/usr/bin/env python3
"""
QW-1789: Cohort-criteria sweep under fixed robust protocol.

Fixed protocol from QW-1787/QW-1788:
- stratified high-coverage fraction = 0.95
- q prior width = 0.20
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
OUT_JSON = ROOT / "report_qw1789_cohort_criteria_sweep.json"
OUT_MD = ROOT / "RAPORT_QW1789_COHORT_CRITERIA_SWEEP.md"


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
    q_center = float(d1773["projection"]["p"])
    q_width = 0.20
    frac = 0.95
    n_rep = 16

    residuals = helper.load_residuals(ROOT / "nano15/residuals/NANOGrav15yr_PulsarTiming_v2.1.0/residuals", max_psr=34)
    positions = helper.load_positions(ROOT / "nano15/parfiles")

    pair_rows: List[Dict[str, float]] = []
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
        pair_rows.append(
            {
                "theta_deg": float(sep),
                "hxy": float(hxy),
                "n_match": int(len(x)),
                "stability": float(stab),
            }
        )

    criteria_grid = [
        {"name": "K0_base", "n_match_min": 120, "stability_max": 0.65},
        {"name": "K1_relax_n", "n_match_min": 115, "stability_max": 0.65},
        {"name": "K2_relax_stab", "n_match_min": 120, "stability_max": 0.70},
        {"name": "K3_balanced", "n_match_min": 110, "stability_max": 0.70},
        {"name": "K4_wider", "n_match_min": 110, "stability_max": 0.75},
        {"name": "K5_widest", "n_match_min": 100, "stability_max": 0.75},
    ]

    results = []
    for j, cfg in enumerate(criteria_grid):
        subset = [
            r
            for r in pair_rows
            if r["n_match"] >= cfg["n_match_min"] and r["stability"] <= cfg["stability_max"]
        ]
        if len(subset) < 85:
            results.append(
                {
                    **cfg,
                    "n_pairs": len(subset),
                    "verdict": "TOO_SMALL",
                    "selection_score": 0.0,
                }
            )
            continue

        theta_all = np.array([r["theta_deg"] for r in subset], dtype=float)
        H_all = np.array([r["hxy"] for r in subset], dtype=float)
        sigma_all = max(float(np.std(H_all)), 1e-6)

        lz0 = helper.evidence_flat(theta_all, H_all, sigma_all, n_mc=9000, seed=6000 + 100 * j + 1)
        lz1 = helper.evidence_legacy(theta_all, H_all, sigma_all, n_mc=12000, seed=6000 + 100 * j + 2)
        lz2 = evidence_reparam_width(helper, theta_all, H_all, sigma_all, q_center=q_center, q_width=q_width, n_mc=15000, seed=6000 + 100 * j + 3)
        full_legacy = float(lz1 - lz0)
        full_reparam = float(lz2 - lz0)
        full_delta = float(full_reparam - full_legacy)

        rng = np.random.default_rng(17890 + j)
        rep_rows = []
        for i in range(n_rep):
            idx = helper.stratified_subsample_indices(theta_all, frac=frac, rng=rng, n_bins=8)
            if len(idx) < 80:
                continue
            th = theta_all[idx]
            hh = H_all[idx]
            sg = max(float(np.std(hh)), 1e-6)
            b0 = helper.evidence_flat(th, hh, sg, n_mc=5000, seed=7000 + 100 * j + 10 * i + 1)
            b1 = helper.evidence_legacy(th, hh, sg, n_mc=6500, seed=7000 + 100 * j + 10 * i + 2)
            b2 = evidence_reparam_width(
                helper, th, hh, sg, q_center=q_center, q_width=q_width, n_mc=8500, seed=7000 + 100 * j + 10 * i + 3
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

        pass_basic = full_reparam > 0.0 and full_delta > 0.0 and prob_rep >= 0.90 and prob_dlt >= 0.95
        score = (
            0.25 * unit_from_scale(max(0.0, 0.30 - full_reparam), 0.30)
            + 0.20 * unit_from_scale(max(0.0, 0.75 - full_delta), 0.75)
            + 0.25 * prob_rep
            + 0.15 * prob_dlt
            + 0.10 * unit_from_scale(std_rep, 0.16)
            + 0.05 * unit_from_scale(std_dlt, 0.16)
        )

        verdict = "SUPPORTED" if pass_basic else "PARTIAL"
        results.append(
            {
                **cfg,
                "n_pairs": len(subset),
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
                "pass_basic": bool(pass_basic),
                "selection_score": float(score),
                "verdict": verdict,
            }
        )

    valid = [r for r in results if r.get("pass_basic")]
    if len(valid) > 0:
        recommended = max(valid, key=lambda r: r["selection_score"])
        rec_strength = "STRONG"
    else:
        recommended = max(results, key=lambda r: r.get("selection_score", 0.0))
        rec_strength = "CONDITIONAL"

    if rec_strength == "STRONG":
        verdict = "COHORT_CRITERIA_SELECTION_SUPPORTED"
    else:
        verdict = "COHORT_CRITERIA_SELECTION_PARTIAL"

    output = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "protocol": {"fraction": frac, "q_width": q_width, "replications": n_rep},
        "n_pairs_precomputed": len(pair_rows),
        "criteria_grid": criteria_grid,
        "results_by_criteria": results,
        "recommended_criteria": recommended,
        "recommendation_strength": rec_strength,
        "verdict": verdict,
    }
    OUT_JSON.write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1789: COHORT CRITERIA SWEEP",
        "",
        f"- Data UTC: {output['generated_utc']}",
        f"- Precomputed pairs: {len(pair_rows)}",
        f"- Fixed protocol: frac={frac:.2f}, q_width={q_width:.2f}",
        f"- Recommendation strength: {rec_strength}",
        f"- Verdict: **{verdict}**",
        "",
        "## Criteria Results",
    ]
    for r in results:
        if r.get("verdict") == "TOO_SMALL":
            lines.append(
                "- {name}: n_pairs={n_pairs} -> TOO_SMALL".format(**r)
            )
            continue
        lines.append(
            "- {name}: n_pairs={n_pairs} | full_reparam={full_logB_reparam:.3f} | full_delta={full_delta_logB:.3f} "
            "| P(reparam>0)={prob_logB_reparam_positive:.3f} | P(delta>0)={prob_delta_logB_positive:.3f} "
            "| std_rep={std_logB_reparam:.3f} | std_delta={std_delta_logB:.3f} "
            "| score={selection_score:.3f} | pass_basic={pass_basic}".format(**r)
        )
    lines.extend(
        [
            "",
            "## Recommended",
            "- name={name}, n_match_min={n_match_min}, stability_max={stability_max}, n_pairs={n_pairs}, score={selection_score:.3f}".format(
                **recommended
            ),
            "",
            "## Artifacts",
            f"- JSON: `{OUT_JSON.name}`",
        ]
    )
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1789] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1789] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
