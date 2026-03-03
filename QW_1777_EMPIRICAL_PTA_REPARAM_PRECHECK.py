#!/usr/bin/env python3
"""
QW-1777: Empirical PTA/GW reparameterization precheck.

Combines existing empirical outputs (phase70/71/75/76 and QW-1725/1726)
with QW-1776 readiness to determine whether the next empirical campaign
can proceed under the new parameter dictionary.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1777_empirical_pta_reparam_precheck.json"
OUT_MD = ROOT / "RAPORT_QW1777_EMPIRICAL_PTA_REPARAM_PRECHECK.md"


def load_json(path: Path) -> Dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    p70 = load_json(ROOT / "nano15/phase70_results.json")
    p71 = load_json(ROOT / "nano15/phase71_bayes_result.json")
    p75 = load_json(ROOT / "nano15/phase75_results.json")
    p76 = load_json(ROOT / "phase76_true_fin_results.json")
    q1725 = load_json(ROOT / "report_qw1725_gw_strict_cross_hurst_reanalysis.json")
    q1726 = load_json(ROOT / "report_qw1726_gw_fin_projection_retest.json")
    q1776 = load_json(ROOT / "report_qw1776_reparam_empirical_readiness_gate.json")

    # Core empirical indicators from existing pipeline
    delta_rss_phase = float(p70.get("delta_rss_phase", 0.0))
    mean_phase = float(p70.get("mean_phase", 0.0))
    logB71 = float(p71.get("logB", 0.0))
    logB_real = float(p75.get("logB_real", 0.0))
    logB_inj = float(p75.get("logB_inj", 0.0))
    sigma_q0 = float(p75.get("sigma_q0", 0.0))
    H_q0 = float(p76.get("H_q0_whitened", 0.0))

    verdict_1725 = str(q1725.get("verdict", ""))
    verdict_1726 = str(q1726.get("verdict", ""))
    readiness_1776 = str(q1776.get("readiness", ""))
    hard_gate_1776 = str(q1776.get("hard_gate", "FAIL"))

    # Preregistered precheck criteria for launching new empirical campaign
    pass_reparam_ready = (hard_gate_1776 == "PASS") and (readiness_1776 == "REPARAM_EMPIRICAL_READY")
    pass_pipeline_sensitivity = (logB_inj >= 8.0) and (0.20 <= sigma_q0 <= 0.60)
    pass_real_structure_signal = (H_q0 >= 0.15) and (delta_rss_phase >= 0.01) and (mean_phase >= 0.65)
    pass_legacy_not_robust = (verdict_1725 == "GW_CROSS_HURST_ANOMALY_NOT_ROBUST") and (verdict_1726 == "FIN_023_TO_031_PROJECTION_NOT_SUPPORTED")
    pass_bayes_headroom = (logB71 <= 0.0) and (logB_real < 0.0)

    if pass_reparam_ready and pass_pipeline_sensitivity and pass_real_structure_signal and pass_legacy_not_robust and pass_bayes_headroom:
        verdict = "EMPIRICAL_REPARAM_CAMPAIGN_READY"
    elif pass_reparam_ready and pass_pipeline_sensitivity and (pass_real_structure_signal or pass_legacy_not_robust):
        verdict = "EMPIRICAL_REPARAM_CAMPAIGN_PARTIAL_READY"
    else:
        verdict = "EMPIRICAL_REPARAM_CAMPAIGN_NOT_READY"

    output = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "inputs": {
            "phase70": p70,
            "phase71": p71,
            "phase75": p75,
            "phase76": p76,
            "q1725_verdict": verdict_1725,
            "q1726_verdict": verdict_1726,
            "q1776_readiness": readiness_1776,
            "q1776_hard_gate": hard_gate_1776,
        },
        "derived_metrics": {
            "delta_rss_phase": delta_rss_phase,
            "mean_phase": mean_phase,
            "logB71": logB71,
            "logB_real": logB_real,
            "logB_inj": logB_inj,
            "sigma_q0": sigma_q0,
            "H_q0_whitened": H_q0,
        },
        "pass_flags": {
            "reparam_ready": bool(pass_reparam_ready),
            "pipeline_sensitivity": bool(pass_pipeline_sensitivity),
            "real_structure_signal": bool(pass_real_structure_signal),
            "legacy_path_not_robust": bool(pass_legacy_not_robust),
            "bayes_headroom_for_reparam": bool(pass_bayes_headroom),
        },
        "verdict": verdict,
        "recommended_next_campaign": {
            "data_domain": "PTA + GW detectors",
            "primary_model": "reparameterized dictionary from QW-1773",
            "validation_mode": "grouped CV + OOD holdout + injection controls",
            "success_target": "logB_reparam > +3 on real data with preserved injection detectability",
        },
    }
    OUT_JSON.write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1777: EMPIRICAL PTA REPARAM PRECHECK",
        "",
        f"- Data UTC: {output['generated_utc']}",
        f"- delta_rss_phase: {delta_rss_phase:.6f}",
        f"- H_q0_whitened: {H_q0:.6f}",
        f"- logB71 / logB_real / logB_inj: {logB71:.3f} / {logB_real:.3f} / {logB_inj:.3f}",
        f"- QW-1776 readiness: {readiness_1776} ({hard_gate_1776})",
        f"- Verdict: **{verdict}**",
        "",
        "## Pass Flags",
        f"- reparam_ready: {pass_reparam_ready}",
        f"- pipeline_sensitivity: {pass_pipeline_sensitivity}",
        f"- real_structure_signal: {pass_real_structure_signal}",
        f"- legacy_path_not_robust: {pass_legacy_not_robust}",
        f"- bayes_headroom_for_reparam: {pass_bayes_headroom}",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1777] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1777] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
