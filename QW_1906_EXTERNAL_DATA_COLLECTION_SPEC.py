#!/usr/bin/env python3
"""
QW-1906: External data collection specification from power-map results.

Translates QW-1905 alpha-based attainability into concrete target conditions
for truly external PTA data collection.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1906_external_data_collection_spec.json"
OUT_MD = ROOT / "RAPORT_QW1906_EXTERNAL_DATA_COLLECTION_SPEC.md"
SPEC_MD = ROOT / "external_confirmatory_v2" / "COLLECTION_SPEC_QW1906.md"


def main() -> None:
    d1905 = json.loads((ROOT / "report_qw1905_external_data_requirements_power_map.json").read_text(encoding="utf-8"))
    d1850 = json.loads((ROOT / "report_qw1850_pta_v2_prereg_protocol.json").read_text(encoding="utf-8"))

    pta_path = ROOT / "external_confirmatory_v2" / "confirmatory_dataset_internal_proxy_wide" / "pta_v2_pairs.csv"
    df = pd.read_csv(pta_path)

    thr = d1850["protocol"]["pta_v2_protocol"]["thresholds"]
    alpha_min_by_n = d1905.get("alpha_min_for_pass_rate_0p80_by_n_pairs", {})
    alpha_req = alpha_min_by_n.get("1200")
    if alpha_req is None:
        alpha_req = alpha_min_by_n.get("2000")
    if alpha_req is None:
        alpha_req = 6.0

    z1 = (df["f_autoc1"].to_numpy(dtype=float) - df["f_autoc1"].mean()) / (df["f_autoc1"].std() + 1e-12)
    z2 = (df["f_switch"].to_numpy(dtype=float) - df["f_switch"].mean()) / (df["f_switch"].std() + 1e-12)
    z3 = (df["f_std"].to_numpy(dtype=float) - df["f_std"].mean()) / (df["f_std"].std() + 1e-12)
    fs = 0.60 * z1 - 0.35 * z2 + 0.25 * z3

    y0 = df["hxy"].to_numpy(dtype=float)
    y_req = np.clip(y0 + float(alpha_req) * 0.05 * fs, 0.0, 1.0)
    dy = y_req - y0

    # Linear coupling slope between hxy and feature score.
    X = np.column_stack([np.ones_like(fs), fs])
    b0, b1 = np.linalg.lstsq(X, y0, rcond=None)[0]
    c0, c1 = np.linalg.lstsq(X, y_req, rcond=None)[0]

    corr0 = float(np.corrcoef(y0, fs)[0, 1]) if np.std(y0) > 1e-12 and np.std(fs) > 1e-12 else 0.0
    corr1 = float(np.corrcoef(y_req, fs)[0, 1]) if np.std(y_req) > 1e-12 and np.std(fs) > 1e-12 else 0.0

    baseline = d1905.get("alpha_rows", [{}])[0].get("pta_summary", {})

    spec = {
        "n_pairs_recommended_min": 1200,
        "n_pairs_preferred": 2000,
        "pta_thresholds_locked": thr,
        "minimum_signal_requirements": {
            "mean_pair_mean_gain_target_min": max(float(thr["mean_pair_mean_gain_min"]), 0.04),
            "prob_pair_mean_gain_positive_target_min": max(float(thr["prob_pair_mean_gain_positive_min"]), 0.667),
            "one_sided_lower95_prob_positive_target_min": max(float(thr["one_sided_lower95_prob_pair_mean_gain_positive_min"]), 0.60),
        },
        "feature_signal_proxy_requirements": {
            "alpha_required_for_80pct_power": float(alpha_req),
            "hxy_feature_slope_baseline": float(b1),
            "hxy_feature_slope_target": float(c1),
            "hxy_feature_corr_baseline": corr0,
            "hxy_feature_corr_target": corr1,
            "delta_hxy_mean_needed": float(np.mean(dy)),
            "delta_hxy_abs_mean_needed": float(np.mean(np.abs(dy))),
            "delta_hxy_q50_needed": float(np.quantile(dy, 0.50)),
            "delta_hxy_q90_needed": float(np.quantile(dy, 0.90)),
            "delta_hxy_positive_fraction": float(np.mean(dy > 0.0)),
        },
        "baseline_pta_state": {
            "mean_pair_mean_gain": float(baseline.get("mean_pair_mean_gain", 0.0)),
            "prob_pair_mean_gain_positive": float(baseline.get("prob_pair_mean_gain_positive", 0.0)),
            "one_sided_lower95_prob_pair_mean_gain_positive": float(
                baseline.get("one_sided_lower95_prob_pair_mean_gain_positive", 0.0)
            ),
        },
        "collection_guidance": [
            "prioritize pulsar pairs with high long-memory coherence and low sign-switching noise",
            "maximize stable cadence and long overlap windows across pair set",
            "avoid mixed-quality cohorts; enforce hard quality floor before confirmatory run",
            "keep protocol thresholds unchanged; model updates require new preregistration",
        ],
    }

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "source_reports": [
            "report_qw1905_external_data_requirements_power_map.json",
            "report_qw1850_pta_v2_prereg_protocol.json",
        ],
        "spec": spec,
        "verdict": "EXTERNAL_DATA_COLLECTION_SPEC_READY",
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1906: EXTERNAL DATA COLLECTION SPEC",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{out['verdict']}**",
        "",
        "## Minimal External Targets (PTA)",
        f"- n_pairs recommended min: {spec['n_pairs_recommended_min']}",
        f"- n_pairs preferred: {spec['n_pairs_preferred']}",
        f"- mean_pair_mean_gain >= {spec['minimum_signal_requirements']['mean_pair_mean_gain_target_min']}",
        f"- prob_pair_mean_gain_positive >= {spec['minimum_signal_requirements']['prob_pair_mean_gain_positive_target_min']}",
        f"- lower95_prob_positive >= {spec['minimum_signal_requirements']['one_sided_lower95_prob_positive_target_min']}",
        "",
        "## Feature-Signal Proxy Targets",
        f"- alpha_required_for_80pct_power: {spec['feature_signal_proxy_requirements']['alpha_required_for_80pct_power']}",
        f"- hxy-feature slope baseline -> target: {spec['feature_signal_proxy_requirements']['hxy_feature_slope_baseline']:.4f} -> {spec['feature_signal_proxy_requirements']['hxy_feature_slope_target']:.4f}",
        f"- hxy-feature corr baseline -> target: {spec['feature_signal_proxy_requirements']['hxy_feature_corr_baseline']:.4f} -> {spec['feature_signal_proxy_requirements']['hxy_feature_corr_target']:.4f}",
        f"- delta_hxy |mean| needed: {spec['feature_signal_proxy_requirements']['delta_hxy_abs_mean_needed']:.4f}",
        f"- delta_hxy positive fraction: {spec['feature_signal_proxy_requirements']['delta_hxy_positive_fraction']:.3f}",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    SPEC_MD.parent.mkdir(parents=True, exist_ok=True)
    spec_lines = [
        "# COLLECTION SPEC QW-1906",
        "",
        "This is an operational checklist for acquiring a truly external frozen PTA/GW package.",
        "",
        "## Hard Requirements",
        f"- PTA pairs >= {spec['n_pairs_recommended_min']} (preferred {spec['n_pairs_preferred']})",
        f"- Keep thresholds fixed from QW-1850 (no retuning)",
        "- Independent provider and explicit externality statement",
        "",
        "## PTA Signal Targets",
        f"- mean_pair_mean_gain >= {spec['minimum_signal_requirements']['mean_pair_mean_gain_target_min']}",
        f"- prob_pair_mean_gain_positive >= {spec['minimum_signal_requirements']['prob_pair_mean_gain_positive_target_min']}",
        f"- one_sided_lower95_prob_positive >= {spec['minimum_signal_requirements']['one_sided_lower95_prob_positive_target_min']}",
        "",
        "## Proxy Translation (planning)",
        f"- alpha_required_for_80pct_power ~= {spec['feature_signal_proxy_requirements']['alpha_required_for_80pct_power']}",
        f"- target hxy-feature slope ~= {spec['feature_signal_proxy_requirements']['hxy_feature_slope_target']:.4f}",
        f"- target hxy-feature corr ~= {spec['feature_signal_proxy_requirements']['hxy_feature_corr_target']:.4f}",
        "",
        "## Next Execution",
        "- Freeze manifest + file hashes",
        "- Run QW-1852 -> QW-1853 -> QW-1902",
        "",
    ]
    SPEC_MD.write_text("\n".join(spec_lines), encoding="utf-8")

    print(f"[QW-1906] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1906] Saved MD:   {OUT_MD.name}")
    print(f"[QW-1906] Saved spec: {SPEC_MD.name}")


if __name__ == "__main__":
    main()
