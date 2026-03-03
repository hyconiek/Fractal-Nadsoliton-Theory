#!/usr/bin/env python3
"""
QW-1860: Trust catalog for FULL_LOG evidence usage.

Converts QW-1859 risk results into actionable evidence classes.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1860_full_log_trust_catalog.json"
OUT_MD = ROOT / "RAPORT_QW1860_FULL_LOG_TRUST_CATALOG.md"


def main() -> None:
    d = json.loads((ROOT / "report_qw1859_full_logs_meta_audit.json").read_text(encoding="utf-8"))
    rows: List[Dict] = d.get("rows_sorted_by_risk", [])

    positive = []
    negative = []
    neutral = []
    excluded = []

    for r in rows:
        file = r.get("file")
        risk = float(r.get("risk_score", 0.0))
        verdict = str(r.get("verdict", ""))
        cnt = r.get("counts", {})

        has_kernel = int(cnt.get("kernel_hits", 0)) > 0
        has_freeze = int(cnt.get("freeze_hits", 0)) > 0
        has_contrad = int(r.get("contradiction_count", 0)) > 0
        frozen_negative = int(cnt.get("frozen_negative_hits", 0)) > 0

        item = {
            "file": file,
            "risk_score": risk,
            "verdict_1859": verdict,
            "kernel_hits": int(cnt.get("kernel_hits", 0)),
            "freeze_hits": int(cnt.get("freeze_hits", 0)),
            "frozen_negative_hits": int(cnt.get("frozen_negative_hits", 0)),
        }

        if has_contrad:
            item["class"] = "EXCLUDE_CONTRADICTORY"
            excluded.append(item)
        elif frozen_negative:
            item["class"] = "NEGATIVE_EVIDENCE_FROZEN_BRANCH"
            negative.append(item)
        elif has_kernel and has_freeze and risk <= 0.15:
            item["class"] = "POSITIVE_CONTEXT_EVIDENCE"
            positive.append(item)
        else:
            item["class"] = "NEUTRAL_CONTEXT"
            neutral.append(item)

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "summary": {
            "n_positive_context": len(positive),
            "n_negative_evidence": len(negative),
            "n_excluded_contradictory": len(excluded),
            "n_neutral_context": len(neutral),
        },
        "positive_context": positive,
        "negative_evidence": negative,
        "excluded_contradictory": excluded,
        "neutral_context": neutral,
        "verdict": "FULL_LOG_TRUST_CATALOG_COMPLETE",
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1860: FULL LOG TRUST CATALOG",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{out['verdict']}**",
        "",
        "## Summary",
        f"- positive context: {out['summary']['n_positive_context']}",
        f"- negative evidence: {out['summary']['n_negative_evidence']}",
        f"- excluded contradictory: {out['summary']['n_excluded_contradictory']}",
        f"- neutral context: {out['summary']['n_neutral_context']}",
        "",
        "## Negative Evidence Files",
    ]

    if not negative:
        lines.append("- none")
    else:
        for x in negative:
            lines.append(
                f"- {x['file']}: risk={x['risk_score']:.2f}, frozen_negative_hits={x['frozen_negative_hits']}"
            )

    lines += [
        "",
        "## Excluded Contradictory Files",
    ]

    if not excluded:
        lines.append("- none")
    else:
        for x in excluded:
            lines.append(
                f"- {x['file']}: risk={x['risk_score']:.2f}, verdict_1859={x['verdict_1859']}"
            )

    lines += [
        "",
        "## Positive Context Files",
    ]

    if not positive:
        lines.append("- none")
    else:
        for x in positive:
            lines.append(
                f"- {x['file']}: risk={x['risk_score']:.2f}, kernel_hits={x['kernel_hits']}, freeze_hits={x['freeze_hits']}"
            )

    lines += [
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]

    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1860] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1860] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
