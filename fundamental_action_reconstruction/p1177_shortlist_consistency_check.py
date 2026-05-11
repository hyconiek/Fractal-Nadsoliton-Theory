#!/usr/bin/env python3
"""P1177 shortlist consistency check.

Verifies shortlist integrity against P1152 summary and P1153 ranking.
"""
from __future__ import annotations
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"

def main() -> None:
    args = sys.argv[1:]
    p1152_path = Path(args[0]).resolve() if args else (GEN / "p1152_strict_candidate_registry_runner_summary.json").resolve()
    p1152 = json.loads(p1152_path.read_text(encoding="utf-8"))
    ranking_path = Path(p1152.get("ranking_summary") or (GEN / "p1153_strict_candidate_quality_ranking_summary.json")).resolve()
    p1153 = json.loads(ranking_path.read_text(encoding="utf-8"))

    results = {r.get("candidate"): r for r in p1152.get("results", [])}
    ranking = {r.get("candidate"): r for r in p1153.get("ranking", [])}
    shortlist = p1152.get("safe_margin_shortlist", [])

    issues = []
    for i, row in enumerate(shortlist):
        c = row.get("candidate")
        if c not in results:
            issues.append(f"shortlist[{i}] candidate missing in P1152 results: {c}")
            continue
        if not results[c].get("pass"):
            issues.append(f"shortlist[{i}] candidate not passed in P1152 results: {c}")
        margin = row.get("safe_region_margin")
        if not isinstance(margin, (int, float)):
            issues.append(f"shortlist[{i}] non-numeric safe_region_margin: {c}")
        rrow = ranking.get(c)
        if rrow is None:
            issues.append(f"shortlist[{i}] candidate missing in P1153 ranking: {c}")
        elif rrow.get("safe_region_margin") != margin:
            issues.append(f"shortlist[{i}] margin mismatch vs P1153: {c}")

    out = {
        "packet": "P1177",
        "as_of": "2026-05-10",
        "p1152_summary": str(p1152_path),
        "p1153_summary": str(ranking_path),
        "checked_shortlist": len(shortlist),
        "overall_pass": len(issues) == 0,
        "issues": issues,
        "note": "Consistency check only; no strict-core closure or QW-2191 discharge claim.",
    }
    out_path = GEN / "p1177_shortlist_consistency_check_summary.json"
    out_path.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[P1177] overall_pass={out['overall_pass']} checked={len(shortlist)} wrote {out_path}")
    if issues:
        raise SystemExit(1)

if __name__ == "__main__":
    main()
