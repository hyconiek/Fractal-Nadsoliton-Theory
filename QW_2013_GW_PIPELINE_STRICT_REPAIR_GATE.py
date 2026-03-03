#!/usr/bin/env python3
"""
QW-2013: Strict GW pipeline repair gate (Task 3).

Goal:
- execute task #3 from closure list: make GW pipeline methodologically hard,
- explicitly map legacy QW-1724 issues to post-repair controls/evidence.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List

ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw2013_gw_pipeline_strict_repair_gate.json"
OUT_MD = ROOT / "RAPORT_QW2013_GW_PIPELINE_STRICT_REPAIR_GATE.md"


def load(name: str) -> dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def main() -> None:
    q1724 = load("report_qw1724_gw_method_audit.json")
    q1725 = load("report_qw1725_gw_strict_cross_hurst_reanalysis.json")
    q2000 = load("report_qw2000_bounded_coupling_deep_audit.json")
    q2001 = load("report_qw2001_bounded_gw_triad_lockable_gate.json")

    det_gw = q2001["deterministic"]["gw"]
    det_flags = q2001["deterministic"]["flags"]
    boot = q2001["bootstrap"]
    deep = q2000["aggregate"]

    # Legacy issue -> repaired control mapping
    checks: List[Dict[str, object]] = []

    # 1) projection contradiction / null-tail control
    c1 = bool(deep["null_random_p90_p90_max"] <= 0.35 and deep["adv_constrained_max"] <= 0.10)
    checks.append(
        {
            "id": 1,
            "legacy_issue": "projection contradiction + unstable null tail",
            "repair_control": "bounded-coupling deep null audit",
            "evidence": {
                "null_random_p90_p90_max": deep["null_random_p90_p90_max"],
                "adv_constrained_max": deep["adv_constrained_max"],
            },
            "pass": c1,
        }
    )

    # 2) scale invariance / local stability
    lpr = q2001["local_deterministic_pass_rates"]
    c2 = bool(lpr["r01"] >= 0.99 and lpr["r02"] >= 0.99 and lpr["r05"] >= 0.99)
    checks.append(
        {
            "id": 2,
            "legacy_issue": "lack of scale/local invariance",
            "repair_control": "local deterministic perturbation pass-rates",
            "evidence": lpr,
            "pass": c2,
        }
    )

    # 3) inter-observatory mismatch
    c3 = bool(det_gw["control_median_gap"] <= 0.0025)
    checks.append(
        {
            "id": 3,
            "legacy_issue": "inter-observatory mismatch",
            "repair_control": "explicit control-gap threshold in lockable GW gate",
            "evidence": {"control_median_gap": det_gw["control_median_gap"], "threshold": 0.0025},
            "pass": c3,
        }
    )

    # 4) discriminative power robustness
    c4 = bool(boot["triad_pass_rate_5000_min"] >= 0.95 and det_gw["auc_h1l1_vs_ctrl"] >= 0.75)
    checks.append(
        {
            "id": 4,
            "legacy_issue": "weak identifiability / unstable significance",
            "repair_control": "bootstrap pass-rate + AUC threshold",
            "evidence": {
                "triad_pass_rate_5000_min": boot["triad_pass_rate_5000_min"],
                "auc": det_gw["auc_h1l1_vs_ctrl"],
            },
            "pass": c4,
        }
    )

    # 5) hard-coded observation risk in legacy branch
    v61 = load("QW_1660_v61_FullNullModel.json")
    legacy_hardcoded = bool("Observation" in v61)
    source_reports = set(q2001["source_reports"])
    v61_in_lockable = any("v61" in s.lower() for s in source_reports)
    c5 = bool(legacy_hardcoded and not v61_in_lockable)
    checks.append(
        {
            "id": 5,
            "legacy_issue": "possible hard-coded observation drift",
            "repair_control": "exclude legacy v61 branch from lockable source pipeline",
            "evidence": {
                "legacy_v61_has_observation_field": legacy_hardcoded,
                "v61_in_lockable_source_reports": v61_in_lockable,
            },
            "pass": c5,
        }
    )

    # 6) phase/lag interpretability debt remains? (from 1725)
    lag = q1725["lag_profile"]
    best_abs_corr = abs(float(lag["best_abs_corr"]))
    c6 = bool(best_abs_corr < 0.01)
    checks.append(
        {
            "id": 6,
            "legacy_issue": "phase/lag inconsistency",
            "repair_control": "strict lag sanity + no-lag-dependent decision rule",
            "evidence": {
                "best_abs_lag_ms": lag["best_abs_lag_ms"],
                "best_abs_corr": lag["best_abs_corr"],
                "note": "small lag-corr indicates no spurious timing-driven detection",
            },
            "pass": c6,
        }
    )

    # 7) all deterministic GW rule flags must pass
    gw_flag_names = [k for k in det_flags if k.startswith("gw_")]
    c7 = bool(all(det_flags[k] for k in gw_flag_names))
    checks.append(
        {
            "id": 7,
            "legacy_issue": "mixed GW acceptance criteria",
            "repair_control": "single explicit deterministic GW flag set",
            "evidence": {k: det_flags[k] for k in gw_flag_names},
            "pass": c7,
        }
    )

    n_pass = sum(1 for c in checks if c["pass"])
    n_total = len(checks)

    if n_pass == n_total:
        verdict = "GW_PIPELINE_STRICT_REPAIR_PASS"
        required_next = "RUN_TRUE_EXTERNAL_BLIND_CONFIRMATORY_WITH_FROZEN_LOCKABLE_PACKAGE"
    elif n_pass >= n_total - 1:
        verdict = "GW_PIPELINE_STRICT_REPAIR_PARTIAL_PASS_WITH_SMALL_DEBT"
        required_next = "CLOSE_REMAINING_DEBT_THEN_EXTERNAL_RUN"
    else:
        verdict = "GW_PIPELINE_STRICT_REPAIR_FAIL"
        required_next = "REWORK_GW_PIPELINE_BEFORE_EXTERNAL_CONFIRMATORY"

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "inputs": {
            "legacy_audit_q1724": q1724["verdict"],
            "legacy_reanalysis_q1725": q1725["verdict"],
            "deep_audit_q2000": q2000["verdict"],
            "lockable_gate_q2001": q2001["verdict"],
        },
        "checks": checks,
        "summary": {
            "n_pass": int(n_pass),
            "n_total": int(n_total),
            "pass_rate": float(n_pass / n_total),
        },
        "verdict": verdict,
        "required_next_step": required_next,
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-2013: GW PIPELINE STRICT REPAIR GATE",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- checks pass: {n_pass}/{n_total}",
        "",
        "## Legacy -> Repair Mapping",
    ]
    for c in checks:
        lines.append(f"- [{c['id']}] pass={c['pass']} | legacy: {c['legacy_issue']} | control: {c['repair_control']}")

    lines += [
        "",
        "## Required Next Step",
        f"- {required_next}",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]

    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-2013] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-2013] Saved MD:   {OUT_MD.name}")
    print(f"[QW-2013] verdict={verdict}, checks={n_pass}/{n_total}")


if __name__ == "__main__":
    main()
