#!/usr/bin/env python3
"""
QW-1856: Conflict ledger and canonical-candidate map for kernel freeze (QW-700..1600).

Adds chronology metadata (file mtimes) to identify where freeze narrative is stable
vs contradictory.
"""

from __future__ import annotations

import json
import re
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1856_kernel_freeze_conflict_ledger_700_1600.json"
OUT_MD = ROOT / "RAPORT_QW1856_KERNEL_FREEZE_CONFLICT_LEDGER_700_1600.md"

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
        "mtimes": [],
        "files": [],
    })

    for p in ROOT.rglob("*"):
        if not is_text_file(p):
            continue
        rel = str(p.relative_to(ROOT))
        if rel.startswith(".git/"):
            continue

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
        mtime = datetime.fromtimestamp(p.stat().st_mtime, tz=timezone.utc).isoformat()

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
            row["mtimes"].append(mtime)
            if len(row["files"]) < 30:
                row["files"].append(rel)

    conflict_rows = []
    canonical_rows = []

    for q in range(700, 1601):
        if q not in by_qw:
            continue
        r = by_qw[q]

        contradictions = {
            "phi_dual": r["phi_pi6"] > 0 and r["phi_zero"] > 0,
            "beta_dual": r["beta_001"] > 0 and r["beta_005"] > 0,
            "nodes_conflict": r["nodes_a"] > 0 and r["nodes_b"] > 0,
        }
        ccount = int(sum(1 for v in contradictions.values() if v))

        item = {
            "qw": q,
            "refs_any": r["refs_any"],
            "refs_kernel": r["refs_kernel"],
            "refs_frozen": r["refs_frozen"],
            "claims": {
                "phi_pi6": r["phi_pi6"],
                "phi_zero": r["phi_zero"],
                "beta_001": r["beta_001"],
                "beta_005": r["beta_005"],
                "nodes_a": r["nodes_a"],
                "nodes_b": r["nodes_b"],
                "omega_pi4": r["omega_pi4"],
            },
            "contradictions": contradictions,
            "contradiction_count": ccount,
            "mtime_min": min(r["mtimes"]) if r["mtimes"] else None,
            "mtime_max": max(r["mtimes"]) if r["mtimes"] else None,
            "files_head": r["files"][:10],
        }

        if r["refs_frozen"] > 0 and ccount > 0:
            conflict_rows.append(item)

        if r["refs_frozen"] > 0 and ccount == 0 and r["refs_kernel"] > 0:
            canonical_rows.append(item)

    conflict_rows.sort(key=lambda x: (-x["contradiction_count"], -x["refs_frozen"], -x["refs_kernel"], x["qw"]))
    canonical_rows.sort(key=lambda x: (-x["refs_frozen"], -x["refs_kernel"], x["qw"]))

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "scope": {"qw_min": 700, "qw_max": 1600, "range_size": 901},
        "summary": {
            "n_qw_with_any_trace": len(by_qw),
            "n_conflict_rows": len(conflict_rows),
            "n_canonical_rows": len(canonical_rows),
        },
        "conflict_rows": conflict_rows[:250],
        "canonical_rows": canonical_rows[:300],
        "verdict": "CONFLICT_LEDGER_PREPARED",
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1856: KERNEL FREEZE CONFLICT LEDGER (QW-700..1600)",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{out['verdict']}**",
        "",
        "## Summary",
        f"- QW ids with trace: {out['summary']['n_qw_with_any_trace']}",
        f"- conflict rows (frozen + contradiction): {out['summary']['n_conflict_rows']}",
        f"- canonical rows (frozen + no contradiction): {out['summary']['n_canonical_rows']}",
        "",
        "## Top Conflict Rows (first 20)",
    ]

    if not conflict_rows:
        lines.append("- none")
    else:
        for r in conflict_rows[:20]:
            lines.append(
                f"- QW-{r['qw']}: c={r['contradiction_count']}, frozen={r['refs_frozen']}, kernel={r['refs_kernel']}"
            )

    lines += [
        "",
        "## Top Canonical Rows (first 20)",
    ]

    if not canonical_rows:
        lines.append("- none")
    else:
        for r in canonical_rows[:20]:
            lines.append(
                f"- QW-{r['qw']}: frozen={r['refs_frozen']}, kernel={r['refs_kernel']}, refs_any={r['refs_any']}"
            )

    lines += [
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]

    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1856] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1856] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
