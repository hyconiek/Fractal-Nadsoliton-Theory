#!/usr/bin/env python3
"""
QW-1828: GW redesign target-threshold map.

Converts failures from QW-1827/1826 into quantitative improvement targets.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1828_gw_target_threshold_map.json"
OUT_MD = ROOT / "RAPORT_QW1828_GW_TARGET_THRESHOLD_MAP.md"


def load(name: str) -> Dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def ratio_needed(current: float, target: float, direction: str) -> float:
    if direction == "up":
        if current <= 0.0:
            return float("inf") if target > 0 else 1.0
        return target / current
    if direction == "down":
        if target <= 0.0:
            return float("inf")
        return current / target
    raise ValueError("direction must be up/down")


def is_met(current: float, target: float, direction: str) -> bool:
    if direction == "up":
        return current >= target
    if direction == "down":
        return current <= target
    raise ValueError("direction must be up/down")


def fmt_ratio(x: float) -> str:
    if x == float("inf"):
        return "inf"
    return f"{x:.3f}"


def main() -> None:
    d1725 = load("report_qw1725_gw_strict_cross_hurst_reanalysis.json")
    d1826 = load("report_qw1826_gw_shared_control_identifiability.json")
    d1827 = load("report_qw1827_gw_redesign_gate.json")

    thr = d1725.get("thresholds", {})
    null_phase = d1725.get("null_phase_randomized", {})
    null_shift = d1725.get("null_circular_shift", {})
    lag = d1725.get("lag_profile", {})
    stab = d1725.get("stability", {})
    base = d1725.get("baseline", {})
    s1826 = d1826.get("summary", {})

    h_h1_l1 = float(base.get("h_obs_h1_l1", 0.0))
    h_h1_v1 = float(base.get("h_h1_v1", 0.0))
    abs_diff = abs(h_h1_v1 - h_h1_l1)

    targets: List[Dict[str, object]] = []

    # 1) Null rejection targets
    p_thr = float(thr.get("null_p_max", 0.01))
    p_phase = float(null_phase.get("p_lower", 1.0))
    p_shift = float(null_shift.get("p_lower", 1.0))
    met = is_met(p_phase, p_thr, "down")
    targets.append(
        {
            "metric": "phase_null_p_lower",
            "current": p_phase,
            "target": p_thr,
            "direction": "down",
            "improvement_factor_needed": ratio_needed(p_phase, p_thr, "down"),
            "gap": p_phase - p_thr,
            "priority": "critical",
            "status": "met" if met else "missing",
        }
    )
    met = is_met(p_shift, p_thr, "down")
    targets.append(
        {
            "metric": "shift_null_p_lower",
            "current": p_shift,
            "target": p_thr,
            "direction": "down",
            "improvement_factor_needed": ratio_needed(p_shift, p_thr, "down"),
            "gap": p_shift - p_thr,
            "priority": "critical",
            "status": "met" if met else "missing",
        }
    )

    # 2) Lag coherence target
    lag_cur = float(lag.get("corr_at_plus_10ms", 0.0))
    lag_tar = float(thr.get("lag_pos_corr_at_10ms_min", 0.02))
    met = is_met(lag_cur, lag_tar, "up")
    targets.append(
        {
            "metric": "corr_at_plus_10ms",
            "current": lag_cur,
            "target": lag_tar,
            "direction": "up",
            "improvement_factor_needed": ratio_needed(lag_cur, lag_tar, "up"),
            "gap": lag_tar - lag_cur,
            "priority": "critical",
            "status": "met" if met else "missing",
        }
    )

    # 3) Scale stability targets
    ls_cur = float(stab.get("length_spread", 9.9))
    ls_tar = float(thr.get("length_spread_max", 0.08))
    ws_cur = float(stab.get("window_std", 9.9))
    ws_tar = float(thr.get("window_std_max", 0.05))
    met = is_met(ls_cur, ls_tar, "down")
    targets.append(
        {
            "metric": "length_spread",
            "current": ls_cur,
            "target": ls_tar,
            "direction": "down",
            "improvement_factor_needed": ratio_needed(ls_cur, ls_tar, "down"),
            "gap": ls_cur - ls_tar,
            "priority": "major",
            "status": "met" if met else "missing",
        }
    )
    met = is_met(ws_cur, ws_tar, "down")
    targets.append(
        {
            "metric": "window_std",
            "current": ws_cur,
            "target": ws_tar,
            "direction": "down",
            "improvement_factor_needed": ratio_needed(ws_cur, ws_tar, "down"),
            "gap": ws_cur - ws_tar,
            "priority": "major",
            "status": "met" if met else "missing",
        }
    )

    # 4) Detector consistency target
    diff_tar = 0.15
    met = is_met(abs_diff, diff_tar, "down")
    targets.append(
        {
            "metric": "abs_diff_H1L1_vs_H1V1",
            "current": abs_diff,
            "target": diff_tar,
            "direction": "down",
            "improvement_factor_needed": ratio_needed(abs_diff, diff_tar, "down"),
            "gap": abs_diff - diff_tar,
            "priority": "critical",
            "status": "met" if met else "missing",
        }
    )

    # 5) Shared-control identifiability targets
    d_cur = float(s1826.get("mean_abs_effect_size_d", 0.0))
    d_tar = 0.35
    n_cur = float(s1826.get("mean_nonoverlap_frac", 0.0))
    n_tar = 0.35
    p_cur = float(s1826.get("prob_shared_beats_control_near_target", 0.0))
    p_tar = 0.70

    met = is_met(d_cur, d_tar, "up")
    targets.append(
        {
            "metric": "mean_abs_effect_size_d",
            "current": d_cur,
            "target": d_tar,
            "direction": "up",
            "improvement_factor_needed": ratio_needed(d_cur, d_tar, "up"),
            "gap": d_tar - d_cur,
            "priority": "major",
            "status": "met" if met else "missing",
        }
    )
    met = is_met(n_cur, n_tar, "up")
    targets.append(
        {
            "metric": "mean_nonoverlap_frac",
            "current": n_cur,
            "target": n_tar,
            "direction": "up",
            "improvement_factor_needed": ratio_needed(n_cur, n_tar, "up"),
            "gap": n_tar - n_cur,
            "priority": "major",
            "status": "met" if met else "missing",
        }
    )
    met = is_met(p_cur, p_tar, "up")
    targets.append(
        {
            "metric": "prob_shared_beats_control_near_target",
            "current": p_cur,
            "target": p_tar,
            "direction": "up",
            "improvement_factor_needed": ratio_needed(p_cur, p_tar, "up"),
            "gap": p_tar - p_cur,
            "priority": "critical",
            "status": "met" if met else "missing",
        }
    )

    # Sort by priority then ratio magnitude
    priority_rank = {"critical": 0, "major": 1, "minor": 2}
    missing = [t for t in targets if t["status"] == "missing"]
    met = [t for t in targets if t["status"] == "met"]
    targets_sorted = sorted(
        missing + met,
        key=lambda t: (priority_rank.get(t["priority"], 9), -1e9 if t["improvement_factor_needed"] == float("inf") else -float(t["improvement_factor_needed"])),
    )

    readiness = d1827.get("readiness", "")

    output = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "gw_readiness_input": readiness,
        "targets": targets_sorted,
        "n_missing_targets": int(len(missing)),
        "n_met_targets": int(len(met)),
        "n_critical_targets": int(sum(1 for t in targets_sorted if t["priority"] == "critical")),
        "n_major_targets": int(sum(1 for t in targets_sorted if t["priority"] == "major")),
    }
    OUT_JSON.write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1828: GW TARGET THRESHOLD MAP",
        "",
        f"- Data UTC: {output['generated_utc']}",
        f"- GW readiness input: {readiness}",
        f"- Missing targets: {output['n_missing_targets']}",
        f"- Already met targets: {output['n_met_targets']}",
        f"- Critical targets: {output['n_critical_targets']}",
        f"- Major targets: {output['n_major_targets']}",
        "",
        "## Targets",
    ]
    for t in targets_sorted:
        lines.append(
            f"- {t['metric']}: current={t['current']:.6g}, target={t['target']:.6g}, "
            f"direction={t['direction']}, factor_needed={fmt_ratio(float(t['improvement_factor_needed']))}, "
            f"priority={t['priority']}, status={t['status']}"
        )

    lines.extend(["", "## Artifacts", f"- JSON: `{OUT_JSON.name}`"])
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1828] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1828] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
