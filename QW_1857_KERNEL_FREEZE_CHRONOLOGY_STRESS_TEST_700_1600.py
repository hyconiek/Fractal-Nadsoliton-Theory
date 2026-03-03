#!/usr/bin/env python3
"""
QW-1857: Chronology stress test for frozen-kernel contradictions (QW-700..1600).

For conflict QWs detected in QW-1856, checks whether contradictory claims appear
at or after first frozen-kernel declaration timestamp.
"""

from __future__ import annotations

import json
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1857_kernel_freeze_chronology_stress_test_700_1600.json"
OUT_MD = ROOT / "RAPORT_QW1857_KERNEL_FREEZE_CHRONOLOGY_STRESS_TEST_700_1600.md"

TEXT_EXTS = {".md", ".tex", ".json", ".py"}
QW_RE = re.compile(r"QW[-_ ]?(\d{1,4})\b", re.IGNORECASE)


def is_text_file(p: Path) -> bool:
    return p.is_file() and p.suffix.lower() in TEXT_EXTS


def read_text_safe(p: Path) -> str:
    try:
        return p.read_text(encoding="utf-8", errors="ignore")
    except Exception:
        return ""


def has_qw_id(text: str, q: int) -> bool:
    pat = re.compile(rf"QW[-_ ]?{q}\b", re.IGNORECASE)
    return bool(pat.search(text))


def flags(text_lower: str) -> Dict[str, bool]:
    return {
        "frozen": (
            "frozen parameter set" in text_lower
            or "kernel frozen" in text_lower
            or ("zamro" in text_lower and "param" in text_lower)
        ),
        "phi_pi6": bool(re.search(r"phi\s*(=|≈|~)?\s*((np\.)?pi\s*/\s*6|0\.5236|π\s*/\s*6)", text_lower)),
        "phi_zero": bool(re.search(r"phi\s*(=|≈|~)?\s*0(\D|$)", text_lower)),
        "beta_001": bool(re.search(r"beta[_\s-]*tors[^\n\r]{0,40}\b0\.0?1\b", text_lower)),
        "beta_005": bool(re.search(r"beta[_\s-]*tors[^\n\r]{0,40}\b0\.0?5\b", text_lower)),
        "nodes_a": ("2,5,8,11" in text_lower) or ("2 5 8 11" in text_lower),
        "nodes_b": ("2,8,14" in text_lower) or ("2 8 14" in text_lower),
    }


def ts(p: Path) -> str:
    return datetime.fromtimestamp(p.stat().st_mtime, tz=timezone.utc).isoformat()


def first_time(events: List[str]) -> Optional[str]:
    if not events:
        return None
    return sorted(events)[0]


def main() -> None:
    d1856 = json.loads((ROOT / "report_qw1856_kernel_freeze_conflict_ledger_700_1600.json").read_text(encoding="utf-8"))
    conflict_qws = [int(x["qw"]) for x in d1856.get("conflict_rows", [])]

    files = [p for p in ROOT.rglob("*") if is_text_file(p)]

    rows = []
    for q in conflict_qws:
        ev = {
            "frozen": [],
            "phi_pi6": [],
            "phi_zero": [],
            "beta_001": [],
            "beta_005": [],
            "nodes_a": [],
            "nodes_b": [],
        }
        files_head = []

        for p in files:
            text = read_text_safe(p)
            if not text:
                continue
            if not has_qw_id(text, q):
                continue

            tl = text.lower()
            fl = flags(tl)
            tstamp = ts(p)

            if len(files_head) < 20:
                files_head.append(str(p.relative_to(ROOT)))

            for k, v in fl.items():
                if v:
                    ev[k].append(tstamp)

        t_frozen = first_time(ev["frozen"])
        t_phi_pi6 = first_time(ev["phi_pi6"])
        t_phi_zero = first_time(ev["phi_zero"])
        t_beta_001 = first_time(ev["beta_001"])
        t_beta_005 = first_time(ev["beta_005"])
        t_nodes_a = first_time(ev["nodes_a"])
        t_nodes_b = first_time(ev["nodes_b"])

        def after_or_equal(a: Optional[str], b: Optional[str]) -> bool:
            if a is None or b is None:
                return False
            return a >= b

        contradiction_after_freeze = {
            "phi_dual_after_freeze": False,
            "beta_dual_after_freeze": False,
            "nodes_conflict_after_freeze": False,
        }

        if t_frozen is not None:
            if t_phi_pi6 is not None and t_phi_zero is not None and (after_or_equal(t_phi_pi6, t_frozen) or after_or_equal(t_phi_zero, t_frozen)):
                contradiction_after_freeze["phi_dual_after_freeze"] = True
            if t_beta_001 is not None and t_beta_005 is not None and (after_or_equal(t_beta_001, t_frozen) or after_or_equal(t_beta_005, t_frozen)):
                contradiction_after_freeze["beta_dual_after_freeze"] = True
            if t_nodes_a is not None and t_nodes_b is not None and (after_or_equal(t_nodes_a, t_frozen) or after_or_equal(t_nodes_b, t_frozen)):
                contradiction_after_freeze["nodes_conflict_after_freeze"] = True

        quarantine = (t_frozen is not None) and any(contradiction_after_freeze.values())

        rows.append(
            {
                "qw": q,
                "first_times": {
                    "frozen": t_frozen,
                    "phi_pi6": t_phi_pi6,
                    "phi_zero": t_phi_zero,
                    "beta_001": t_beta_001,
                    "beta_005": t_beta_005,
                    "nodes_a": t_nodes_a,
                    "nodes_b": t_nodes_b,
                },
                "contradiction_after_freeze": contradiction_after_freeze,
                "quarantine_from_canonical_freeze": quarantine,
                "files_head": files_head,
            }
        )

    n_quarantine = sum(1 for r in rows if r["quarantine_from_canonical_freeze"])

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "input_conflict_qws": conflict_qws,
        "rows": rows,
        "summary": {
            "n_conflict_qws": len(conflict_qws),
            "n_quarantine": n_quarantine,
        },
        "verdict": "CHRONOLOGY_STRESS_TEST_COMPLETE",
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1857: KERNEL FREEZE CHRONOLOGY STRESS TEST (QW-700..1600)",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{out['verdict']}**",
        f"- conflict QWs: {len(conflict_qws)}",
        f"- quarantine count: {n_quarantine}",
        "",
        "## Rows",
    ]

    for r in rows:
        lines.append(
            f"- QW-{r['qw']}: quarantine={r['quarantine_from_canonical_freeze']} | after_freeze={r['contradiction_after_freeze']}"
        )

    lines += [
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]

    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1857] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1857] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
