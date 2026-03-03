#!/usr/bin/env python3
"""
QW-1827: GW redesign gate with prioritized failure decomposition.

Aggregates strict GW diagnostics and outputs prioritized redesign targets.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1827_gw_redesign_gate.json"
OUT_MD = ROOT / "RAPORT_QW1827_GW_REDESIGN_GATE.md"


def load(name: str) -> Dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def main() -> None:
    d1725 = load("report_qw1725_gw_strict_cross_hurst_reanalysis.json")
    d1726 = load("report_qw1726_gw_fin_projection_retest.json")
    d1826 = load("report_qw1826_gw_shared_control_identifiability.json")

    thr = d1725.get("thresholds", {})
    null_phase = d1725.get("null_phase_randomized", {})
    null_shift = d1725.get("null_circular_shift", {})
    stab = d1725.get("stability", {})
    lag = d1725.get("lag_profile", {})
    base = d1725.get("baseline", {})

    h_h1_l1 = float(base.get("h_obs_h1_l1", 0.0))
    h_h1_v1 = float(base.get("h_h1_v1", 0.0))

    checks: List[Dict[str, object]] = []

    # C1: null-model rejection quality
    pmax = float(thr.get("null_p_max", 0.01))
    c1_phase = float(null_phase.get("p_lower", 1.0)) <= pmax
    c1_shift = float(null_shift.get("p_lower", 1.0)) <= pmax
    c1_score = 0.5 * float(c1_phase) + 0.5 * float(c1_shift)
    checks.append(
        {
            "domain": "Null rejection (phase + shift)",
            "score": c1_score,
            "status": "PASS" if c1_score >= 0.75 else "FAIL",
            "metric": {
                "phase_p_lower": float(null_phase.get("p_lower", 1.0)),
                "shift_p_lower": float(null_shift.get("p_lower", 1.0)),
                "p_threshold": pmax,
            },
            "priority": "critical" if c1_score < 0.25 else "major",
        }
    )

    # C2: lag coherence
    lag_min = float(thr.get("lag_pos_corr_at_10ms_min", 0.02))
    lag_val = float(lag.get("corr_at_plus_10ms", 0.0))
    c2_pass = lag_val >= lag_min
    c2_score = min(max(lag_val / max(lag_min, 1e-12), 0.0), 1.0)
    checks.append(
        {
            "domain": "Lag coherence at +10ms",
            "score": c2_score,
            "status": "PASS" if c2_pass else "FAIL",
            "metric": {"corr_at_10ms": lag_val, "threshold": lag_min},
            "priority": "critical" if c2_score < 0.25 else "major",
        }
    )

    # C3: scale stability
    spread_max = float(thr.get("length_spread_max", 0.08))
    wstd_max = float(thr.get("window_std_max", 0.05))
    length_spread = float(stab.get("length_spread", 9.9))
    window_std = float(stab.get("window_std", 9.9))
    c3a = length_spread <= spread_max
    c3b = window_std <= wstd_max
    s3a = max(0.0, 1.0 - length_spread / max(spread_max, 1e-12))
    s3b = max(0.0, 1.0 - window_std / max(wstd_max, 1e-12))
    c3_score = 0.5 * s3a + 0.5 * s3b
    checks.append(
        {
            "domain": "Scale/window stability",
            "score": c3_score,
            "status": "PASS" if (c3a and c3b) else "FAIL",
            "metric": {
                "length_spread": length_spread,
                "length_spread_max": spread_max,
                "window_std": window_std,
                "window_std_max": wstd_max,
            },
            "priority": "major" if c3_score < 0.5 else "minor",
        }
    )

    # C4: detector consistency (simple absolute discrepancy)
    diff_det = abs(h_h1_v1 - h_h1_l1)
    diff_thr = 0.15
    c4_pass = diff_det <= diff_thr
    c4_score = max(0.0, 1.0 - diff_det / diff_thr)
    checks.append(
        {
            "domain": "Detector consistency (H1-L1 vs H1-V1)",
            "score": c4_score,
            "status": "PASS" if c4_pass else "FAIL",
            "metric": {
                "h_h1_l1": h_h1_l1,
                "h_h1_v1": h_h1_v1,
                "abs_diff": diff_det,
                "threshold": diff_thr,
            },
            "priority": "critical" if c4_score < 0.25 else "major",
        }
    )

    # C5: shared-control identifiability (from QW-1826)
    s1826 = d1826.get("summary", {})
    p1826 = d1826.get("pass_flags", {})
    c5_pass = bool(p1826.get("effect_separation") and p1826.get("quantile_nonoverlap") and p1826.get("target_prevalence_advantage"))
    c5_score = (
        0.35 * min(float(s1826.get("mean_abs_effect_size_d", 0.0)) / 0.35, 1.0)
        + 0.35 * min(float(s1826.get("mean_nonoverlap_frac", 0.0)) / 0.35, 1.0)
        + 0.30 * float(s1826.get("prob_shared_beats_control_near_target", 0.0))
    )
    checks.append(
        {
            "domain": "Shared-control identifiability (1826)",
            "score": c5_score,
            "status": "PASS" if c5_pass else "FAIL",
            "metric": {
                "mean_abs_d": float(s1826.get("mean_abs_effect_size_d", 0.0)),
                "mean_nonoverlap": float(s1826.get("mean_nonoverlap_frac", 0.0)),
                "prob_shared_beats_control": float(s1826.get("prob_shared_beats_control_near_target", 0.0)),
            },
            "priority": "critical" if c5_score < 0.4 else "major",
        }
    )

    global_score = float(sum(float(c["score"]) for c in checks) / len(checks))
    hard_gate = all(c["status"] == "PASS" for c in checks)

    if hard_gate:
        readiness = "GW_BRANCH_READY"
        recommendation = "MERGE_GW_WITH_PTA_QUANTILE_PROTOCOL"
    elif global_score >= 0.60:
        readiness = "GW_BRANCH_PARTIAL_REDESIGN"
        recommendation = "TARGET_CRITICAL_FAILURES_FIRST"
    else:
        readiness = "GW_BRANCH_REDESIGN_REQUIRED"
        recommendation = "DO_NOT_RUN_JOINT_GW_PTA_CAMPAIGN_YET"

    priority_order = sorted(checks, key=lambda c: float(c["score"]))

    action_items = [
        "Rebuild GW statistic around event-windowed, detector-coherent features (chirp-conditioned) instead of global cross-Hurst only.",
        "Introduce explicit shared-vs-unshared contrast objective and require positive near-target prevalence advantage across SNR sweep.",
        "Enforce multi-detector consistency gate (H1-L1, H1-V1, L1-V1) before any projection claim.",
        "Require lag-coherence threshold around physical delays as hard precondition for structure claims.",
    ]

    output = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "checks": checks,
        "global_score": global_score,
        "hard_gate": "PASS" if hard_gate else "FAIL",
        "readiness": readiness,
        "recommendation": recommendation,
        "priority_order": [c["domain"] for c in priority_order],
        "action_items": action_items,
    }
    OUT_JSON.write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1827: GW REDESIGN GATE",
        "",
        f"- Data UTC: {output['generated_utc']}",
        f"- Global score: {global_score:.3f}",
        f"- Hard gate: **{output['hard_gate']}**",
        f"- Readiness: **{readiness}**",
        f"- Recommendation: **{recommendation}**",
        "",
        "## Checks",
    ]
    for c in checks:
        lines.append(f"- {c['domain']}: {c['status']} | score={c['score']:.3f} | priority={c['priority']}")

    lines.extend(["", "## Priority Order"]) 
    for dom in output["priority_order"]:
        lines.append(f"- {dom}")

    lines.extend(["", "## Action Items"])
    for it in action_items:
        lines.append(f"- {it}")

    lines.extend(["", "## Artifacts", f"- JSON: `{OUT_JSON.name}`"])
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1827] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1827] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
