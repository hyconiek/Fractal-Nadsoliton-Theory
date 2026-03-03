#!/usr/bin/env python3
"""
QW-1859: Meta-audit of all FULL_LOG*.md files.

Evaluates kernel/freeze narrative consistency and risk per full-log file.
"""

from __future__ import annotations

import json
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1859_full_logs_meta_audit.json"
OUT_MD = ROOT / "RAPORT_QW1859_FULL_LOGS_META_AUDIT.md"

QW_RE = re.compile(r"QW[-_ ]?(\d{1,4})\b", re.IGNORECASE)


def marker_counts(lines: List[str]) -> Dict:
    freeze_hits = []
    kernel_hits = []
    phi_pi6_hits = []
    phi_zero_hits = []
    beta_001_hits = []
    beta_005_hits = []
    nodes_a_hits = []
    nodes_b_hits = []
    omega_pi4_hits = []
    frozen_negative_hits = []
    frozen_positive_hits = []

    check_ok = 0
    check_fail = 0

    for i, ln in enumerate(lines, start=1):
        low = ln.lower()

        if "✅" in ln:
            check_ok += 1
        if "❌" in ln:
            check_fail += 1

        if any(k in low for k in ["frozen parameter", "kernel frozen", "zamro", "freeze", "frozen"]):
            freeze_hits.append(i)
            if any(k in low for k in ["poraż", "poraz", "kryzys", "fail", "falsyf", "nie generuje", "całkowita porażka", "calkowita porazka"]):
                frozen_negative_hits.append(i)
            if any(k in low for k in ["sukces", "success", "potwierd", "verified", "przelom", "przełom"]):
                frozen_positive_hits.append(i)

        if "kernel" in low or "kernela" in low:
            kernel_hits.append(i)

        if re.search(r"phi\s*(=|≈|~)?\s*((np\.)?pi\s*/\s*6|0\.5236|π\s*/\s*6)", low):
            phi_pi6_hits.append(i)

        if re.search(r"phi\s*(=|≈|~)?\s*0(\D|$)", low):
            phi_zero_hits.append(i)

        if re.search(r"beta[_\s-]*tors[^\n\r]{0,40}\b0\.0?1\b", low):
            beta_001_hits.append(i)

        if re.search(r"beta[_\s-]*tors[^\n\r]{0,40}\b0\.0?5\b", low):
            beta_005_hits.append(i)

        if re.search(r"omega\s*(=|≈|~)?\s*(np\.pi\s*/\s*4|pi\s*/\s*4|0\.7854)", low):
            omega_pi4_hits.append(i)

        if "2,5,8,11" in low or "2 5 8 11" in low:
            nodes_a_hits.append(i)

        if "2,8,14" in low or "2 8 14" in low:
            nodes_b_hits.append(i)

    contradictions = {
        "phi_dual_definition": len(phi_pi6_hits) > 0 and len(phi_zero_hits) > 0,
        "beta_dual_definition": len(beta_001_hits) > 0 and len(beta_005_hits) > 0,
        "node_set_conflict": len(nodes_a_hits) > 0 and len(nodes_b_hits) > 0,
    }
    contradiction_count = int(sum(1 for v in contradictions.values() if v))

    # risk score for quick ranking
    risk = 0.0
    risk += 0.35 * contradiction_count
    if len(freeze_hits) > 0 and len(frozen_negative_hits) > 0:
        risk += 0.25
    if len(freeze_hits) > 0 and len(frozen_positive_hits) > 0 and len(frozen_negative_hits) > 0:
        risk += 0.15
    if len(kernel_hits) == 0 and len(freeze_hits) > 0:
        risk += 0.10
    risk = min(1.0, risk)

    if contradiction_count >= 2:
        verdict = "FULL_LOG_INTERNAL_CONTRADICTIONS_HIGH"
    elif contradiction_count == 1:
        verdict = "FULL_LOG_INTERNAL_CONTRADICTIONS_PRESENT"
    elif len(freeze_hits) > 0 and len(frozen_negative_hits) > 0:
        verdict = "FULL_LOG_FROZEN_BRANCH_EMPIRICALLY_NEGATIVE"
    elif len(freeze_hits) > 0:
        verdict = "FULL_LOG_FROZEN_BRANCH_NO_INTERNAL_CONTRADICTION"
    else:
        verdict = "FULL_LOG_NO_STRONG_FROZEN_SIGNAL"

    return {
        "counts": {
            "freeze_hits": len(freeze_hits),
            "kernel_hits": len(kernel_hits),
            "phi_pi6_hits": len(phi_pi6_hits),
            "phi_zero_hits": len(phi_zero_hits),
            "beta_001_hits": len(beta_001_hits),
            "beta_005_hits": len(beta_005_hits),
            "nodes_a_hits": len(nodes_a_hits),
            "nodes_b_hits": len(nodes_b_hits),
            "omega_pi4_hits": len(omega_pi4_hits),
            "frozen_negative_hits": len(frozen_negative_hits),
            "frozen_positive_hits": len(frozen_positive_hits),
            "check_ok": check_ok,
            "check_fail": check_fail,
        },
        "contradictions": contradictions,
        "contradiction_count": contradiction_count,
        "risk_score": risk,
        "verdict": verdict,
        "samples": {
            "freeze_lines_head": freeze_hits[:12],
            "negative_frozen_lines_head": frozen_negative_hits[:12],
        },
    }


def main() -> None:
    full_logs = sorted([p for p in ROOT.glob("FULL_LOG*.md") if p.is_file()])

    rows = []
    all_qw_ids = set()

    for p in full_logs:
        text = p.read_text(encoding="utf-8", errors="ignore")
        lines = text.splitlines()

        qids = set(int(m.group(1)) for m in QW_RE.finditer(text))
        all_qw_ids.update(qids)

        mc = marker_counts(lines)

        row = {
            "file": p.name,
            "size": {
                "lines": len(lines),
                "chars": len(text),
            },
            "qw_coverage": {
                "unique_qw_ids": len(qids),
                "qw_min": min(qids) if qids else None,
                "qw_max": max(qids) if qids else None,
            },
            **mc,
        }
        rows.append(row)

    # risk descending
    rows_sorted = sorted(rows, key=lambda x: (-x["risk_score"], -x["contradiction_count"], -x["size"]["lines"]))

    summary = {
        "n_full_logs": len(rows),
        "n_qw_ids_union": len(all_qw_ids),
        "high_risk_files": sum(1 for r in rows if r["risk_score"] >= 0.75),
        "medium_risk_files": sum(1 for r in rows if 0.40 <= r["risk_score"] < 0.75),
        "low_risk_files": sum(1 for r in rows if r["risk_score"] < 0.40),
        "files_with_any_contradiction": sum(1 for r in rows if r["contradiction_count"] > 0),
    }

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "summary": summary,
        "rows_sorted_by_risk": rows_sorted,
        "verdict": "FULL_LOGS_META_AUDIT_COMPLETE",
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1859: FULL LOGS META AUDIT",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{out['verdict']}**",
        "",
        "## Summary",
        f"- full log files: {summary['n_full_logs']}",
        f"- union of QW ids: {summary['n_qw_ids_union']}",
        f"- high/medium/low risk: {summary['high_risk_files']} / {summary['medium_risk_files']} / {summary['low_risk_files']}",
        f"- files with contradictions: {summary['files_with_any_contradiction']}",
        "",
        "## Risk Table",
    ]

    for r in rows_sorted:
        lines.append(
            f"- {r['file']}: risk={r['risk_score']:.2f}, contradictions={r['contradiction_count']}, verdict={r['verdict']}"
        )

    lines += [
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]

    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1859] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1859] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
