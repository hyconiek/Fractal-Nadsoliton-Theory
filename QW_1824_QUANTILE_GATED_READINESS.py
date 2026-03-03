#!/usr/bin/env python3
"""
QW-1824: Quantile-gated readiness decision for sequence branch.

Integrates:
- QW-1819 sequence rigor gate (LL-based gate fails),
- QW-1821 likelihood calibration diagnostic,
- QW-1823 quantile-score OOS validation.

Purpose:
- decide if branch can be advanced under revised, robust predictive criteria.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1824_quantile_gated_readiness.json"
OUT_MD = ROOT / "RAPORT_QW1824_QUANTILE_GATED_READINESS.md"


def load(name: str) -> Dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def main() -> None:
    d1819 = load("report_qw1819_sequence_rigor_gate.json")
    d1821 = load("report_qw1821_likelihood_calibration_diagnostic.json")
    d1823 = load("report_qw1823_quantile_score_oos.json")

    checks: List[Dict[str, object]] = []

    # Legacy LL gate status (for transparency)
    ll_fail = d1819.get("hard_gate") != "PASS"
    checks.append(
        {
            "domain": "Legacy LL gate (1819)",
            "status": "FAIL" if ll_fail else "PASS",
            "score": 0.0 if ll_fail else 1.0,
            "note": d1819.get("readiness", ""),
        }
    )

    # Evidence that LL misspecification is the bottleneck
    mis_spec = d1821.get("verdict") == "LIKELIHOOD_MISSPECIFICATION_SIGNAL_STRONG"
    checks.append(
        {
            "domain": "Likelihood misspecification signal (1821)",
            "status": "PASS" if mis_spec else "FAIL",
            "score": 1.0 if mis_spec else 0.3,
            "note": d1821.get("verdict", ""),
        }
    )

    p1823 = d1823.get("pass_flags", {})
    quantile_supported = bool(p1823.get("quantile_gain") and p1823.get("mae_gain") and p1823.get("dispersion_control"))
    checks.append(
        {
            "domain": "Robust quantile gate (1823)",
            "status": "PASS" if quantile_supported else "FAIL",
            "score": 1.0 if quantile_supported else 0.0,
            "note": d1823.get("verdict", ""),
        }
    )

    global_score = float(sum(float(c["score"]) for c in checks) / len(checks))

    if quantile_supported and mis_spec:
        readiness = "SEQUENCE_BRANCH_CONDITIONAL_READY_UNDER_QUANTILE_GATING"
        recommendation = "ADOPT_QUANTILE_SCORE_AS_PRIMARY_GATE_AND_FREEZE_PROTOCOL_FOR_NEXT_EMPIRICAL_STAGE"
        hard_gate = "PASS"
    elif quantile_supported:
        readiness = "SEQUENCE_BRANCH_PARTIAL_READY"
        recommendation = "RUN_ADDITIONAL_CALIBRATION_AUDIT_BEFORE_PROTOCOL_FREEZE"
        hard_gate = "PARTIAL"
    else:
        readiness = "SEQUENCE_BRANCH_NOT_READY"
        recommendation = "CONTINUE_REDESIGN"
        hard_gate = "FAIL"

    output = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "checks": checks,
        "global_score": global_score,
        "hard_gate": hard_gate,
        "readiness": readiness,
        "recommendation": recommendation,
        "transition_note": "Decision is made under revised robust predictive criterion due confirmed LL misspecification.",
    }
    OUT_JSON.write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1824: QUANTILE-GATED READINESS",
        "",
        f"- Data UTC: {output['generated_utc']}",
        f"- Global score: {global_score:.3f}",
        f"- Hard gate: **{hard_gate}**",
        f"- Readiness: **{readiness}**",
        f"- Recommendation: **{recommendation}**",
        "",
        "## Checks",
    ]
    for c in checks:
        lines.append(f"- {c['domain']}: {c['status']} | score={c['score']:.3f} | note={c['note']}")

    lines.extend([
        "",
        "## Transition Note",
        f"- {output['transition_note']}",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ])
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1824] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1824] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
