#!/usr/bin/env python3
"""
QW-1855: Deep-read priority queue for kernel-freeze narrative in QW-700..1600.

Builds per-QW evidence profile and prioritizes IDs for manual scientific review.
"""

from __future__ import annotations

import json
import re
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1855_kernel_freeze_deepread_queue_700_1600.json"
OUT_MD = ROOT / "RAPORT_QW1855_KERNEL_FREEZE_DEEPREAD_QUEUE_700_1600.md"

TEXT_EXTS = {".md", ".tex", ".json", ".py"}
QW_RE = re.compile(r"QW[-_ ]?(\d{1,4})\b", re.IGNORECASE)


def is_text_file(p: Path) -> bool:
    return p.is_file() and p.suffix.lower() in TEXT_EXTS


def read_text_safe(p: Path) -> str:
    try:
        return p.read_text(encoding="utf-8", errors="ignore")
    except Exception:
        return ""


def detect_flags(t: str) -> Dict[str, bool]:
    return {
        "kernel": "kernel" in t or "kernela" in t,
        "frozen": (
            "frozen parameter set" in t
            or "kernel frozen" in t
            or ("zamro" in t and "param" in t)
        ),
        "phi_pi6": bool(re.search(r"phi\s*(=|≈|~)?\s*((np\.)?pi\s*/\s*6|0\.5236|π\s*/\s*6)", t)),
        "phi_zero": bool(re.search(r"phi\s*(=|≈|~)?\s*0(\D|$)", t)),
        "beta_001": bool(re.search(r"beta[_\s-]*tors[^\n\r]{0,40}\b0\.0?1\b", t)),
        "beta_005": bool(re.search(r"beta[_\s-]*tors[^\n\r]{0,40}\b0\.0?5\b", t)),
        "nodes_a": ("2,5,8,11" in t) or ("2 5 8 11" in t),
        "nodes_b": ("2,8,14" in t) or ("2 8 14" in t),
        "omega_pi4": bool(re.search(r"omega\s*(=|≈|~)?\s*((np\.)?pi\s*/\s*4|0\.7854|π\s*/\s*4)", t)),
    }


def main() -> None:
    files = [p for p in ROOT.rglob("*") if is_text_file(p)]

    by_qw = defaultdict(lambda: {
        "refs_any": 0,
        "refs_kernel": 0,
        "refs_frozen": 0,
        "phi_pi6": 0,
        "phi_zero": 0,
        "beta_001": 0,
        "beta_005": 0,
        "nodes_a": 0,
        "nodes_b": 0,
        "omega_pi4": 0,
        "files": [],
    })

    for p in files:
        text = read_text_safe(p)
        if not text:
            continue
        t = text.lower()

        qids = set()
        for m in QW_RE.finditer(text):
            try:
                q = int(m.group(1))
            except Exception:
                continue
            if 700 <= q <= 1600:
                qids.add(q)

        if not qids:
            continue

        flags = detect_flags(t)

        for q in qids:
            row = by_qw[q]
            row["refs_any"] += 1
            if flags["kernel"]:
                row["refs_kernel"] += 1
            if flags["frozen"]:
                row["refs_frozen"] += 1
            for k in ["phi_pi6", "phi_zero", "beta_001", "beta_005", "nodes_a", "nodes_b", "omega_pi4"]:
                if flags[k]:
                    row[k] += 1
            if len(row["files"]) < 20:
                row["files"].append(str(p.relative_to(ROOT)))

    range_ids = list(range(700, 1601))
    missing_ids = [q for q in range_ids if q not in by_qw]

    queue = []
    canonical_candidates = []

    for q in range_ids:
        if q not in by_qw:
            queue.append({
                "qw": q,
                "priority": "P0_MISSING",
                "risk_score": 1.0,
                "reason": "No trace found in scanned files for this QW id.",
            })
            continue

        r = by_qw[q]
        contradictions = {
            "phi_dual": r["phi_pi6"] > 0 and r["phi_zero"] > 0,
            "beta_dual": r["beta_001"] > 0 and r["beta_005"] > 0,
            "nodes_conflict": r["nodes_a"] > 0 and r["nodes_b"] > 0,
        }
        ccount = int(sum(1 for v in contradictions.values() if v))

        score = 0.0
        if r["refs_kernel"] == 0:
            score += 0.35
        if r["refs_frozen"] > 0 and ccount > 0:
            score += 0.45
        score += min(0.20, 0.10 * ccount)

        if score >= 0.75:
            prio = "P0_CONTRADICTORY_FREEZE"
        elif score >= 0.45:
            prio = "P1_HIGH_RISK"
        elif score >= 0.20:
            prio = "P2_MEDIUM_RISK"
        else:
            prio = "P3_LOW_RISK"

        queue.append({
            "qw": q,
            "priority": prio,
            "risk_score": round(score, 3),
            "refs_any": r["refs_any"],
            "refs_kernel": r["refs_kernel"],
            "refs_frozen": r["refs_frozen"],
            "contradictions": contradictions,
            "files_head": r["files"][:8],
        })

        if r["refs_frozen"] > 0 and ccount == 0 and r["refs_kernel"] > 0:
            canonical_candidates.append({
                "qw": q,
                "refs_frozen": r["refs_frozen"],
                "refs_kernel": r["refs_kernel"],
                "refs_any": r["refs_any"],
            })

    # sort by highest risk first, then QW ascending
    prio_order = {
        "P0_MISSING": 0,
        "P0_CONTRADICTORY_FREEZE": 1,
        "P1_HIGH_RISK": 2,
        "P2_MEDIUM_RISK": 3,
        "P3_LOW_RISK": 4,
    }
    queue.sort(key=lambda x: (prio_order.get(x["priority"], 9), -float(x.get("risk_score", 0.0)), int(x["qw"])))
    canonical_candidates.sort(key=lambda x: (-x["refs_frozen"], -x["refs_kernel"], int(x["qw"])))

    counts = defaultdict(int)
    for q in queue:
        counts[q["priority"]] += 1

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "scope": {"qw_min": 700, "qw_max": 1600, "range_size": 901},
        "summary": {
            "n_missing_ids": len(missing_ids),
            "n_with_trace": 901 - len(missing_ids),
            "priority_counts": dict(counts),
            "n_canonical_candidates": len(canonical_candidates),
        },
        "missing_ids_head": missing_ids[:120],
        "deepread_queue_head": queue[:200],
        "canonical_candidate_head": canonical_candidates[:80],
        "verdict": "DEEPREAD_QUEUE_PREPARED",
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1855: KERNEL FREEZE DEEPREAD QUEUE (QW-700..1600)",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{out['verdict']}**",
        "",
        "## Summary",
        f"- missing ids: {out['summary']['n_missing_ids']}",
        f"- ids with trace: {out['summary']['n_with_trace']}",
        f"- priority counts: {out['summary']['priority_counts']}",
        f"- canonical candidates count: {out['summary']['n_canonical_candidates']}",
        "",
        "## Top Queue (first 25)",
    ]

    for row in queue[:25]:
        lines.append(
            f"- QW-{row['qw']}: {row['priority']} | risk={row.get('risk_score')}"
        )

    lines += [
        "",
        "## Top Canonical Candidates (first 20)",
    ]

    if not canonical_candidates:
        lines.append("- none")
    else:
        for row in canonical_candidates[:20]:
            lines.append(
                f"- QW-{row['qw']}: frozen_refs={row['refs_frozen']}, kernel_refs={row['refs_kernel']}, any_refs={row['refs_any']}"
            )

    lines += [
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]

    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1855] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1855] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
