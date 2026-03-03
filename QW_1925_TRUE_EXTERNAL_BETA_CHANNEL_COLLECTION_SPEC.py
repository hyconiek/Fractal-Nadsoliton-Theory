#!/usr/bin/env python3
"""
QW-1925: True external beta-channel collection specification.

Builds a strict, executable collection spec from latest beta-channel branch:
- QW-1920 (gap diagnosis),
- QW-1921 (orthogonal observable selected),
- QW-1922 (blind intervention pass),
- QW-1924 (lambda tuning indicates data-limited bottleneck).
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1925_true_external_beta_channel_collection_spec.json"
OUT_MD = ROOT / "RAPORT_QW1925_TRUE_EXTERNAL_BETA_CHANNEL_COLLECTION_SPEC.md"
SPEC_MD = ROOT / "external_confirmatory_v2" / "COLLECTION_SPEC_QW1925_TRUE_EXTERNAL_BETA_CHANNEL.md"


def z_quantile(p: float) -> float:
    # Acklam approximation for inverse normal CDF.
    # Accurate enough for planning-level power calculations.
    a = [
        -3.969683028665376e+01,
        2.209460984245205e+02,
        -2.759285104469687e+02,
        1.383577518672690e+02,
        -3.066479806614716e+01,
        2.506628277459239e+00,
    ]
    b = [
        -5.447609879822406e+01,
        1.615858368580409e+02,
        -1.556989798598866e+02,
        6.680131188771972e+01,
        -1.328068155288572e+01,
    ]
    c = [
        -7.784894002430293e-03,
        -3.223964580411365e-01,
        -2.400758277161838e+00,
        -2.549732539343734e+00,
        4.374664141464968e+00,
        2.938163982698783e+00,
    ]
    d = [
        7.784695709041462e-03,
        3.224671290700398e-01,
        2.445134137142996e+00,
        3.754408661907416e+00,
    ]
    plow = 0.02425
    phigh = 1.0 - plow

    if p <= 0.0:
        return -1e9
    if p >= 1.0:
        return 1e9

    if p < plow:
        q = np.sqrt(-2.0 * np.log(p))
        return (
            (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5])
            / ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1.0)
        )
    if p > phigh:
        q = np.sqrt(-2.0 * np.log(1.0 - p))
        return -(
            (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5])
            / ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1.0)
        )

    q = p - 0.5
    r = q * q
    return (
        (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * q
    ) / (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1.0)


def required_n(effect: float, sigma: float, alpha_one_sided: float = 0.025, power: float = 0.90) -> int:
    z_a = z_quantile(1.0 - alpha_one_sided)
    z_b = z_quantile(power)
    n = ((z_a + z_b) * sigma / max(effect, 1e-9)) ** 2
    return int(np.ceil(n))


def main() -> None:
    d1920 = json.loads((ROOT / "report_qw1920_high_power_identifiability_interior_stability.json").read_text(encoding="utf-8"))
    d1921 = json.loads((ROOT / "report_qw1921_orthogonal_beta_observable_design.json").read_text(encoding="utf-8"))
    d1922 = json.loads((ROOT / "report_qw1922_beta_observable_blind_external_intervention.json").read_text(encoding="utf-8"))
    d1924 = json.loads((ROOT / "report_qw1924_lambda_tuning_and_transfer_retest.json").read_text(encoding="utf-8"))

    selected_obs = str(d1921["discovery"]["selected_candidate"])

    p_hold = d1922["datasets"]["primary"]["holdout"]
    s_hold = d1922["datasets"]["stress"]["holdout"]

    c_primary = float(p_hold["contrast"])
    c_stress = float(s_hold["contrast"])

    # Conservative planning target and uncertainty from bootstrap intervals.
    c_target = float(min(c_primary, c_stress))
    c_target_conservative = float(0.60 * c_target)

    # Approximate sigma from 90% CI width (q95-q05 ~= 3.29 sigma if near-normal).
    p_q05 = float(p_hold["contrast_bootstrap"]["q05"])
    p_q95 = float(p_hold["contrast_bootstrap"]["q95"])
    s_q05 = float(s_hold["contrast_bootstrap"]["q05"])
    s_q95 = float(s_hold["contrast_bootstrap"]["q95"])

    sigma_p = float((p_q95 - p_q05) / 3.29)
    sigma_s = float((s_q95 - s_q05) / 3.29)
    sigma_conservative = float(max(sigma_p, sigma_s, 0.12))

    n_req_80_eff = required_n(effect=c_target_conservative, sigma=sigma_conservative, alpha_one_sided=0.025, power=0.80)
    n_req_90_eff = required_n(effect=c_target_conservative, sigma=sigma_conservative, alpha_one_sided=0.025, power=0.90)

    # Inflate for non-idealities (cohort imbalance, metadata loss, quality exclusions),
    # then enforce strict floors for dependent pair-level data.
    inflation = 1.35
    n_req_80_adj = int(np.ceil(n_req_80_eff * inflation))
    n_req_90_adj = int(np.ceil(n_req_90_eff * inflation))
    n_req_80_adj = max(400, n_req_80_adj)
    n_req_90_adj = max(500, n_req_90_adj)

    spec = {
        "selected_beta_observable": selected_obs,
        "hard_requirements": {
            "provider_must_be_external_independent": True,
            "forbidden_provider_tokens": ["INTERNAL_PROXY", "INTERNAL"],
            "required_externality_statement_tokens": ["independent", "external"],
            "required_files": [
                "manifest_beta_channel.json",
                "beta_channel_pairs.csv",
                "intervention_events.csv",
                "protocol_freeze.json",
            ],
            "required_roles_in_manifest": [
                "beta_pairs",
                "intervention_events",
                "protocol_freeze",
            ],
        },
        "dataset_schema": {
            "beta_channel_pairs_required_columns": [
                "pair_id",
                "theta_deg",
                "hxy",
                "f_std",
                "f_autoc1",
                "f_switch",
                "f_slope",
                "intervention_id",
                "regime",
            ],
            "intervention_events_required_columns": [
                "intervention_id",
                "intervention_type",
                "source_reference",
                "start_utc",
                "end_utc",
                "is_preregistered",
                "is_exogenous",
            ],
            "regime_allowed_values": ["pre", "post"],
        },
        "minimum_signal_targets": {
            "holdout_effect_beta_min": 0.60,
            "holdout_effect_omega_max": 0.25,
            "holdout_contrast_min": 0.35,
            "holdout_contrast_bootstrap_q05_min": 0.20,
        },
        "power_targets": {
            "contrast_target_conservative": c_target_conservative,
            "sigma_conservative": sigma_conservative,
            "n_holdout_effective_iid_power_0p80": n_req_80_eff,
            "n_holdout_effective_iid_power_0p90": n_req_90_eff,
            "n_holdout_min_for_power_0p80": n_req_80_adj,
            "n_holdout_min_for_power_0p90": n_req_90_adj,
            "n_total_pairs_recommended": int(max(2 * n_req_90_adj, 1200)),
        },
        "protocol_rules": [
            "Freeze triad estimator and beta observable definition before ingest.",
            "No threshold changes after manifest freeze.",
            "Discovery/holdout split fixed by sha256(pair_id) parity.",
            "At least two exogenous intervention cohorts required.",
            "Primary analysis one-sided alpha=0.025 with bootstrap lower bound criterion.",
        ],
        "provenance": {
            "from_qw1920_verdict": d1920.get("verdict"),
            "from_qw1921_verdict": d1921.get("verdict"),
            "from_qw1922_verdict": d1922.get("verdict"),
            "from_qw1924_verdict": d1924.get("verdict"),
        },
    }

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "source_reports": [
            "report_qw1920_high_power_identifiability_interior_stability.json",
            "report_qw1921_orthogonal_beta_observable_design.json",
            "report_qw1922_beta_observable_blind_external_intervention.json",
            "report_qw1924_lambda_tuning_and_transfer_retest.json",
        ],
        "spec": spec,
        "verdict": "TRUE_EXTERNAL_BETA_CHANNEL_COLLECTION_SPEC_READY",
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1925: TRUE EXTERNAL BETA-CHANNEL COLLECTION SPEC",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{out['verdict']}**",
        f"- selected_beta_observable: `{selected_obs}`",
        "",
        "## Minimal Signal Targets",
        f"- holdout_effect_beta >= {spec['minimum_signal_targets']['holdout_effect_beta_min']}",
        f"- holdout_effect_omega <= {spec['minimum_signal_targets']['holdout_effect_omega_max']}",
        f"- holdout_contrast >= {spec['minimum_signal_targets']['holdout_contrast_min']}",
        f"- holdout_contrast_boot_q05 >= {spec['minimum_signal_targets']['holdout_contrast_bootstrap_q05_min']}",
        "",
        "## Power Targets",
        f"- conservative contrast target: {c_target_conservative:.4f}",
        f"- conservative sigma: {sigma_conservative:.4f}",
        f"- n_holdout min (80% power): {n_req_80_adj}",
        f"- n_holdout min (90% power): {n_req_90_adj}",
        f"- n_total_pairs recommended: {spec['power_targets']['n_total_pairs_recommended']}",
        "",
        "## Required Files",
        "- manifest_beta_channel.json",
        "- beta_channel_pairs.csv",
        "- intervention_events.csv",
        "- protocol_freeze.json",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    SPEC_MD.parent.mkdir(parents=True, exist_ok=True)
    spec_lines = [
        "# COLLECTION SPEC QW-1925: TRUE EXTERNAL BETA CHANNEL",
        "",
        "This is a strict operational checklist for collecting true external intervention data",
        "to resolve Stage-B beta identifiability in FIN-ToE.",
        "",
        "## Hard Requirements",
        "- External independent provider only (no INTERNAL_PROXY / INTERNAL tokens).",
        "- Manifest must include roles: beta_pairs, intervention_events, protocol_freeze.",
        "- At least two exogenous intervention cohorts.",
        "- Freeze protocol before metric execution.",
        "",
        "## File Package",
        "- manifest_beta_channel.json",
        "- beta_channel_pairs.csv",
        "- intervention_events.csv",
        "- protocol_freeze.json",
        "",
        "## Minimal Signal Targets",
        f"- holdout_effect_beta >= {spec['minimum_signal_targets']['holdout_effect_beta_min']}",
        f"- holdout_effect_omega <= {spec['minimum_signal_targets']['holdout_effect_omega_max']}",
        f"- holdout_contrast >= {spec['minimum_signal_targets']['holdout_contrast_min']}",
        f"- holdout_contrast_boot_q05 >= {spec['minimum_signal_targets']['holdout_contrast_bootstrap_q05_min']}",
        "",
        "## Power Targets",
        f"- n_holdout min (80%): {n_req_80_adj}",
        f"- n_holdout min (90%): {n_req_90_adj}",
        f"- n_total_pairs recommended: {spec['power_targets']['n_total_pairs_recommended']}",
        "",
        "## Execution After Collection",
        "- Run QW-1927 readiness gate (schema/externality/protocol lock).",
        "- Run blind intervention evaluation with frozen thresholds.",
        "- Re-run Stage-B gate with updated evidence.",
        "",
    ]
    SPEC_MD.write_text("\n".join(spec_lines), encoding="utf-8")

    print(f"[QW-1925] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1925] Saved MD:   {OUT_MD.name}")
    print(f"[QW-1925] Saved spec: {SPEC_MD.name}")


if __name__ == "__main__":
    main()
