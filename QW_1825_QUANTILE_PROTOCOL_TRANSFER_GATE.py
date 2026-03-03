#!/usr/bin/env python3
"""
QW-1825: Quantile-protocol transfer gate (PTA -> GW).

Combines current evidence to evaluate whether the new quantile-gated protocol
can be considered transferable across empirical domains.

Inputs:
- QW-1824: sequence branch conditional readiness under quantile gating,
- QW-1823: robust quantile OOS performance (PTA residuals),
- QW-1786: campaign robustness confirmation (PTA operational branch),
- QW-1725/QW-1726: GW robustness and projection re-test status.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1825_quantile_protocol_transfer_gate.json"
OUT_MD = ROOT / "RAPORT_QW1825_QUANTILE_PROTOCOL_TRANSFER_GATE.md"


def load(name: str) -> Dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def main() -> None:
    d1823 = load("report_qw1823_quantile_score_oos.json")
    d1824 = load("report_qw1824_quantile_gated_readiness.json")
    d1786 = load("report_qw1786_campaign_robustness_gate.json")
    d1725 = load("report_qw1725_gw_strict_cross_hurst_reanalysis.json")
    d1726 = load("report_qw1726_gw_fin_projection_retest.json")

    checks: List[Dict[str, object]] = []

    # 1) Quantile branch validity on PTA sequence benchmark
    p1823 = d1823.get("pass_flags", {})
    pass_q_pta = bool(p1823.get("quantile_gain") and p1823.get("mae_gain") and p1823.get("dispersion_control"))
    score_q_pta = 1.0 if pass_q_pta else 0.0
    checks.append(
        {
            "domain": "PTA quantile predictive branch (1823)",
            "score": score_q_pta,
            "status": "PASS" if pass_q_pta else "FAIL",
            "note": d1823.get("verdict", ""),
        }
    )

    # 2) Protocol-level decision status
    ready_1824 = d1824.get("hard_gate") == "PASS"
    score_1824 = 1.0 if ready_1824 else 0.0
    checks.append(
        {
            "domain": "Quantile-gated readiness decision (1824)",
            "score": score_1824,
            "status": "PASS" if ready_1824 else "FAIL",
            "note": d1824.get("readiness", ""),
        }
    )

    # 3) Legacy PTA operational robustness (for continuity)
    pta_campaign_pass = d1786.get("hard_gate") == "PASS"
    score_pta_campaign = float(d1786.get("global_score", 0.0)) if pta_campaign_pass else 0.0
    checks.append(
        {
            "domain": "PTA operational campaign robustness (1786)",
            "score": score_pta_campaign,
            "status": "PASS" if pta_campaign_pass else "FAIL",
            "note": d1786.get("readiness", ""),
        }
    )

    # 4) GW transfer feasibility from strict audits
    gw_not_robust = d1725.get("verdict") == "GW_CROSS_HURST_ANOMALY_NOT_ROBUST"
    gw_proj_not_supported = d1726.get("verdict") == "FIN_023_TO_031_PROJECTION_NOT_SUPPORTED"
    gw_fail_count = int(gw_not_robust) + int(gw_proj_not_supported)
    score_gw_transfer = 1.0 - 0.5 * gw_fail_count  # 1.0, 0.5, 0.0
    checks.append(
        {
            "domain": "GW transfer feasibility (1725+1726)",
            "score": score_gw_transfer,
            "status": "PASS" if score_gw_transfer >= 0.75 else "FAIL",
            "note": f"{d1725.get('verdict','')} / {d1726.get('verdict','')}",
        }
    )

    global_score = float(sum(float(c["score"]) for c in checks) / len(checks))

    pta_ready = pass_q_pta and ready_1824 and pta_campaign_pass
    gw_ready = score_gw_transfer >= 0.75

    if pta_ready and gw_ready:
        readiness = "QUANTILE_PROTOCOL_CROSS_DOMAIN_READY"
        hard_gate = "PASS"
        recommendation = "START_FULL_PTA_GW_JOINT_EMPIRICAL_RUN"
    elif pta_ready and not gw_ready:
        readiness = "QUANTILE_PROTOCOL_PTA_READY_GW_BLOCKED"
        hard_gate = "PARTIAL"
        recommendation = "RUN_GW_SPECIFIC_REDESIGN_BEFORE_JOINT_CAMPAIGN"
    else:
        readiness = "QUANTILE_PROTOCOL_NOT_READY"
        hard_gate = "FAIL"
        recommendation = "DO_NOT_FREEZE_CROSS_DOMAIN_PROTOCOL"

    blockers = []
    if not gw_ready:
        blockers.append("GW branch fails strict robustness and projection criteria (QW-1725/QW-1726).")
    if not pass_q_pta:
        blockers.append("PTA quantile branch not fully supported (QW-1823).")
    if not ready_1824:
        blockers.append("Quantile-gated readiness decision not passed (QW-1824).")

    output = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "checks": checks,
        "global_score": global_score,
        "hard_gate": hard_gate,
        "readiness": readiness,
        "recommendation": recommendation,
        "pta_ready": pta_ready,
        "gw_ready": gw_ready,
        "blockers": blockers,
    }
    OUT_JSON.write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1825: QUANTILE PROTOCOL TRANSFER GATE",
        "",
        f"- Data UTC: {output['generated_utc']}",
        f"- Global score: {global_score:.3f}",
        f"- Hard gate: **{hard_gate}**",
        f"- Readiness: **{readiness}**",
        f"- Recommendation: **{recommendation}**",
        f"- PTA ready: **{pta_ready}**",
        f"- GW ready: **{gw_ready}**",
        "",
        "## Checks",
    ]
    for c in checks:
        lines.append(f"- {c['domain']}: {c['status']} | score={c['score']:.3f} | note={c['note']}")

    lines.append("")
    lines.append("## Blockers")
    if blockers:
        for b in blockers:
            lines.append(f"- {b}")
    else:
        lines.append("- None")

    lines.extend(["", "## Artifacts", f"- JSON: `{OUT_JSON.name}`"])
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1825] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1825] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
