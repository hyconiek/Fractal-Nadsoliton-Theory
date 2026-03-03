#!/usr/bin/env python3
"""
QW-1839: Joint PTA+GW confirmatory prereg protocol freeze.

Creates immutable protocol spec for next confirmatory stage after QW-1838.
"""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1839_joint_confirmatory_prereg_protocol.json"
OUT_MD = ROOT / "RAPORT_QW1839_JOINT_CONFIRMATORY_PREREG_PROTOCOL.md"


def load(name: str) -> Dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def main() -> None:
    d1823 = load("report_qw1823_quantile_score_oos.json")
    d1836 = load("report_qw1836_gw_control_calibrated_objective.json")
    d1838 = load("report_qw1838_global_reparam_readiness_gate.json")

    if d1838.get("hard_gate") != "PASS":
        raise RuntimeError("QW-1838 is not PASS; prereg freeze not allowed.")

    pta_ref = d1823.get("summary", {})
    gw_ref = d1836.get("summary", {}).get("calibrated", {})

    protocol = {
        "stage": "JOINT_CONFIRMATORY_V1",
        "frozen_utc": datetime.now(timezone.utc).isoformat(),
        "eligibility": {
            "qw1838_hard_gate_required": "PASS",
            "qw1838_readiness_required": "GLOBAL_CONDITIONAL_READY_UNDER_REPARAM_CRITERIA",
        },
        "pta_protocol": {
            "source_script": "QW_1823_QUANTILE_SCORE_OOS.py",
            "metric_primary": "mean_quantile_gain_m2_minus_m2e",
            "metric_secondary": "prob_quantile_gain_positive",
            "metric_stability": "std_quantile_gain_m2_minus_m2e",
            "thresholds": {
                "mean_quantile_gain_m2_minus_m2e_min": 0.040,
                "prob_quantile_gain_positive_min": 0.90,
                "std_quantile_gain_m2_minus_m2e_max": 0.035,
            },
            "reference_values_qw1823": {
                "mean_quantile_gain_m2_minus_m2e": pta_ref.get("mean_quantile_gain_m2_minus_m2e"),
                "prob_quantile_gain_positive": pta_ref.get("prob_quantile_gain_positive"),
                "std_quantile_gain_m2_minus_m2e": pta_ref.get("std_quantile_gain_m2_minus_m2e"),
            },
            "split_rule": "blocked_cv_by_window_idx_mod_5",
            "seed": 21923,
        },
        "gw_protocol": {
            "source_script": "QW_1836_GW_CONTROL_CALIBRATED_OBJECTIVE.py",
            "metric_primary": "calibrated_mean_auc",
            "metric_secondary": "calibrated_mean_adv",
            "metric_consistency": "calibrated_mean_control_gap",
            "thresholds": {
                "calibrated_mean_auc_min": 0.90,
                "calibrated_mean_adv_min": 0.60,
                "calibrated_mean_control_gap_max": 0.0005,
                "calibrated_prob_adv_positive_min": 0.90,
            },
            "reference_values_qw1836": {
                "calibrated_mean_auc": gw_ref.get("mean_auc"),
                "calibrated_mean_adv": gw_ref.get("mean_adv"),
                "calibrated_mean_control_gap": gw_ref.get("mean_control_gap"),
                "calibrated_prob_adv_positive": gw_ref.get("prob_adv_positive"),
            },
            "calibration_rule": "control_pair_offsets_from_train_only",
            "split_rule": "blocked_cv_by_window_idx_mod_5",
            "seed": 21924,
        },
        "joint_gate": {
            "logic": "AND",
            "require_all_pta_thresholds": True,
            "require_all_gw_thresholds": True,
            "report_both_domain_scores": True,
        },
        "anti_leakage": {
            "forbidden": [
                "fitting_offsets_on_test_folds",
                "manual_threshold_tuning_after_results",
                "changing_split_rule_posthoc",
            ],
            "allowed": [
                "pre-registered rerun with new data only",
                "adding external holdout as separate phase",
            ],
        },
    }

    canonical = json.dumps(protocol, ensure_ascii=False, sort_keys=True, separators=(",", ":"))
    protocol_hash = hashlib.sha256(canonical.encode("utf-8")).hexdigest()

    output = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "protocol": protocol,
        "protocol_sha256": protocol_hash,
        "verdict": "JOINT_CONFIRMATORY_PREREG_FROZEN",
    }
    OUT_JSON.write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")

    t_pta = protocol["pta_protocol"]["thresholds"]
    t_gw = protocol["gw_protocol"]["thresholds"]

    lines = [
        "# RAPORT QW-1839: JOINT CONFIRMATORY PREREG PROTOCOL",
        "",
        f"- Data UTC: {output['generated_utc']}",
        f"- Verdict: **{output['verdict']}**",
        f"- Protocol SHA256: `{protocol_hash}`",
        "",
        "## PTA Thresholds",
        f"- mean_quantile_gain >= {t_pta['mean_quantile_gain_m2_minus_m2e_min']}",
        f"- prob_quantile_gain_positive >= {t_pta['prob_quantile_gain_positive_min']}",
        f"- std_quantile_gain <= {t_pta['std_quantile_gain_m2_minus_m2e_max']}",
        "",
        "## GW Thresholds",
        f"- calibrated_mean_auc >= {t_gw['calibrated_mean_auc_min']}",
        f"- calibrated_mean_adv >= {t_gw['calibrated_mean_adv_min']}",
        f"- calibrated_mean_control_gap <= {t_gw['calibrated_mean_control_gap_max']}",
        f"- calibrated_prob_adv_positive >= {t_gw['calibrated_prob_adv_positive_min']}",
        "",
        "## Anti-Leakage",
        "- Offsets/calibrations train-only.",
        "- Split rule frozen (window_idx mod 5).",
        "- No post-hoc threshold tuning.",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1839] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1839] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
