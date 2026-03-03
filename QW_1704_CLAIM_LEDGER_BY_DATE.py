#!/usr/bin/env python3
"""
QW-1704: Chronological claim ledger for FIN documentation.

Builds a date-ordered ledger of QW statuses from markdown evidence and
classifies each QW into D/C/M/F/P-style practical buckets.
"""

from __future__ import annotations

import json
import re
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1704_claim_ledger_by_date.json"
OUT_MD = ROOT / "RAPORT_QW1704_CLAIM_LEDGER_BY_DATE.md"

QW_RE = re.compile(r"\bQW[\s\-_]?V?(\d{1,4})\b", re.IGNORECASE)
RAPORT_QW_RE = re.compile(r"\braport[_\-]?qw[\s\-_]?(\d{1,4})\b", re.IGNORECASE)

POSITIVE_WORDS = (
    "success",
    "sukces",
    "confirmed",
    "verified",
    "pass",
    "robust",
    "breakthrough",
    "działa",
    "exact",
)

NEGATIVE_WORDS = (
    "fail",
    "failed",
    "falsified",
    "porażk",
    "nie działa",
    "inconsistent",
    "artifact",
    "critical",
    "error",
)


def extract_qw_ids(text: str):
    out = set()
    for rx in (QW_RE, RAPORT_QW_RE):
        out.update(int(m.group(1)) for m in rx.finditer(text))
    return out


def count_hits(text_lower: str, words):
    return sum(text_lower.count(w) for w in words)


def classify(pos: int, neg: int, exact: int, fit: int, taut: int) -> str:
    if neg >= max(5, int(1.2 * (pos + 1))) and (taut + fit) >= 2:
        return "F"  # practical falsification / critical
    if exact > 0 and fit == 0 and taut == 0 and pos > neg:
        return "D"  # strongest derived-like
    if fit > 0 or taut > 0:
        return "C"  # calibration-risk class
    if pos > 0 and neg > 0:
        return "M"  # mixed
    if pos > 0:
        return "M"
    return "P"  # pending/insufficient evidence


def main():
    md_files = []
    for p in ROOT.rglob("*.md"):
        if any(x in p.parts for x in (".git", ".venv", "__pycache__")):
            continue
        md_files.append(p)

    qw_docs = defaultdict(list)
    monthly = defaultdict(lambda: Counter())

    for p in md_files:
        text = p.read_text(encoding="utf-8", errors="ignore")
        lower = text.lower()
        ids = sorted(extract_qw_ids(p.name + "\n" + text))
        if not ids:
            continue
        mtime = datetime.fromtimestamp(p.stat().st_mtime, tz=timezone.utc)
        month = mtime.strftime("%Y-%m")
        pos = count_hits(lower, POSITIVE_WORDS)
        neg = count_hits(lower, NEGATIVE_WORDS)
        exact = lower.count("exact") + lower.count("0.00%")
        fit = lower.count("fitting") + lower.count("fit")
        taut = lower.count("tautolog")

        rec = {
            "path": str(p.relative_to(ROOT)),
            "mtime_utc": mtime.isoformat(),
            "positive_hits": pos,
            "negative_hits": neg,
            "exact_claims": exact,
            "fitting_mentions": fit,
            "tautology_mentions": taut,
        }
        for q in ids:
            qw_docs[q].append(rec)

        monthly[month]["docs"] += 1
        monthly[month]["positive_hits"] += pos
        monthly[month]["negative_hits"] += neg
        monthly[month]["exact_claims"] += exact
        monthly[month]["fitting_mentions"] += fit
        monthly[month]["tautology_mentions"] += taut

    ledger = []
    class_counts = Counter()
    for q, docs in sorted(qw_docs.items()):
        docs_sorted = sorted(docs, key=lambda x: x["mtime_utc"])
        pos = sum(d["positive_hits"] for d in docs_sorted)
        neg = sum(d["negative_hits"] for d in docs_sorted)
        exact = sum(d["exact_claims"] for d in docs_sorted)
        fit = sum(d["fitting_mentions"] for d in docs_sorted)
        taut = sum(d["tautology_mentions"] for d in docs_sorted)
        cls = classify(pos, neg, exact, fit, taut)
        class_counts[cls] += 1
        ledger.append(
            {
                "qw_id": q,
                "class": cls,
                "n_docs": len(docs_sorted),
                "first_seen_utc": docs_sorted[0]["mtime_utc"],
                "last_seen_utc": docs_sorted[-1]["mtime_utc"],
                "positive_hits": pos,
                "negative_hits": neg,
                "exact_claims": exact,
                "fitting_mentions": fit,
                "tautology_mentions": taut,
                "sample_docs": [d["path"] for d in docs_sorted[:3]],
            }
        )

    output = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "root": str(ROOT),
        "ledger_size": len(ledger),
        "class_counts": dict(class_counts),
        "timeline_by_month": {k: dict(v) for k, v in sorted(monthly.items())},
        "ledger_top_400": ledger[:400],
    }
    OUT_JSON.write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1704: CLAIM LEDGER BY DATE",
        "",
        f"- Data UTC: {output['generated_utc']}",
        f"- Liczba QW z udokumentowanym śladem: {len(ledger)}",
        "",
        "## Rozkład klas",
    ]
    for cls in ("D", "C", "M", "F", "P"):
        lines.append(f"- {cls}: {class_counts.get(cls, 0)}")

    lines.extend(
        [
            "",
            "## Trend miesięczny (MD)",
        ]
    )
    for month, vals in sorted(monthly.items()):
        lines.append(
            f"- {month}: docs={vals['docs']}, pos={vals['positive_hits']}, neg={vals['negative_hits']}, exact={vals['exact_claims']}, fit={vals['fitting_mentions']}, taut={vals['tautology_mentions']}"
        )

    lines.extend(
        [
            "",
            "## Próbka 20 wpisów",
        ]
    )
    for row in ledger[:20]:
        lines.append(
            f"- QW-{row['qw_id']} [{row['class']}], docs={row['n_docs']}, first={row['first_seen_utc']}, last={row['last_seen_utc']}, neg={row['negative_hits']}, fit={row['fitting_mentions']}, taut={row['tautology_mentions']}"
        )

    lines.extend(
        [
            "",
            "## Artefakty",
            f"- JSON szczegółowy: `{OUT_JSON.name}`",
        ]
    )
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"[QW-1704] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1704] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
