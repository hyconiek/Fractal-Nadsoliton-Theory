#!/usr/bin/env python3
"""
QW-1858: Dedicated audit of FULL_LOG_COMPRESSED_EXTREME_QW420_1200.md.

Focus:
- internal kernel/freeze narrative consistency,
- parameter declarations and contradictions,
- QW-range coverage and structural signals.
"""

from __future__ import annotations

import json
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List


ROOT = Path(__file__).resolve().parent
IN_FILE = ROOT / "FULL_LOG_COMPRESSED_EXTREME_QW420_1200.md"
OUT_JSON = ROOT / "report_qw1858_full_log_extreme_qw420_1200_audit.json"
OUT_MD = ROOT / "RAPORT_QW1858_FULL_LOG_EXTREME_QW420_1200_AUDIT.md"

QW_RE = re.compile(r"QW[-_ ]?(\d{1,4})\b", re.IGNORECASE)


def main() -> None:
    if not IN_FILE.exists():
        raise FileNotFoundError(f"Missing file: {IN_FILE}")

    text = IN_FILE.read_text(encoding="utf-8", errors="ignore")
    lines = text.splitlines()

    qw_ids = set()
    qw_ids_420_1200 = set()

    freeze_hits = []
    kernel_hits = []
    phi_pi6_hits = []
    phi_zero_hits = []
    beta_001_hits = []
    beta_005_hits = []
    omega_pi4_hits = []
    nodes_a_hits = []
    nodes_b_hits = []
    frozen_negative_hits = []
    frozen_positive_hits = []
    check_ok = 0
    check_fail = 0

    section_headers = []

    for i, ln in enumerate(lines, start=1):
        low = ln.lower()

        for m in QW_RE.finditer(ln):
            q = int(m.group(1))
            qw_ids.add(q)
            if 420 <= q <= 1200:
                qw_ids_420_1200.add(q)

        if ln.startswith("## "):
            section_headers.append({"line": i, "header": ln.strip()})

        if "✅" in ln:
            check_ok += 1
        if "❌" in ln:
            check_fail += 1

        if any(k in low for k in ["frozen parameter", "kernel frozen", "zamro", "freeze", "frozen"]):
            freeze_hits.append({"line": i, "text": ln[:220]})
            if any(k in low for k in ["poraż", "poraz", "kryzys", "fail", "falsyf", "nie generuje", "całkowita porażka", "calkowita porazka"]):
                frozen_negative_hits.append({"line": i, "text": ln[:220]})
            if any(k in low for k in ["sukces", "success", "potwierd", "verified", "przelom", "przełom"]):
                frozen_positive_hits.append({"line": i, "text": ln[:220]})

        if "kernel" in low or "kernela" in low:
            kernel_hits.append({"line": i, "text": ln[:220]})

        if re.search(r"phi\s*(=|≈|~)?\s*((np\.)?pi\s*/\s*6|0\.5236|π\s*/\s*6)", low):
            phi_pi6_hits.append({"line": i, "text": ln[:220]})

        if re.search(r"phi\s*(=|≈|~)?\s*0(\D|$)", low):
            phi_zero_hits.append({"line": i, "text": ln[:220]})

        if re.search(r"beta[_\s-]*tors[^\n\r]{0,40}\b0\.0?1\b", low):
            beta_001_hits.append({"line": i, "text": ln[:220]})

        if re.search(r"beta[_\s-]*tors[^\n\r]{0,40}\b0\.0?5\b", low):
            beta_005_hits.append({"line": i, "text": ln[:220]})

        if re.search(r"omega\s*(=|≈|~)?\s*(np\.pi\s*/\s*4|pi\s*/\s*4|0\.7854)", low):
            omega_pi4_hits.append({"line": i, "text": ln[:220]})

        if "2,5,8,11" in low or "2 5 8 11" in low:
            nodes_a_hits.append({"line": i, "text": ln[:220]})

        if "2,8,14" in low or "2 8 14" in low:
            nodes_b_hits.append({"line": i, "text": ln[:220]})

    contradictions = {
        "phi_dual_definition": len(phi_pi6_hits) > 0 and len(phi_zero_hits) > 0,
        "beta_dual_definition": len(beta_001_hits) > 0 and len(beta_005_hits) > 0,
        "node_set_conflict": len(nodes_a_hits) > 0 and len(nodes_b_hits) > 0,
    }

    contradiction_count = int(sum(1 for v in contradictions.values() if v))

    n_lines = len(lines)
    n_chars = len(text)

    # Heuristic verdict for this file only
    if contradiction_count >= 1 and len(frozen_negative_hits) > 0:
        verdict = "FULL_LOG_FROZEN_BRANCH_CONTRADICTORY_AND_EMPIRICALLY_NEGATIVE"
    elif contradiction_count == 0 and len(freeze_hits) > 0 and len(frozen_negative_hits) > 0:
        verdict = "FULL_LOG_FROZEN_BRANCH_EMPIRICALLY_NEGATIVE"
    elif contradiction_count == 0 and len(freeze_hits) > 0:
        verdict = "FULL_LOG_KERNEL_NARRATIVE_INTERNAL_OK"
    elif contradiction_count >= 2:
        verdict = "FULL_LOG_KERNEL_NARRATIVE_INTERNAL_CONTRADICTORY"
    else:
        verdict = "FULL_LOG_KERNEL_NARRATIVE_MIXED"

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "input_file": str(IN_FILE),
        "size": {
            "lines": n_lines,
            "chars": n_chars,
        },
        "qw_coverage": {
            "unique_qw_ids_any": len(qw_ids),
            "unique_qw_ids_420_1200": len(qw_ids_420_1200),
            "qw_420_1200_min": min(qw_ids_420_1200) if qw_ids_420_1200 else None,
            "qw_420_1200_max": max(qw_ids_420_1200) if qw_ids_420_1200 else None,
        },
        "markers": {
            "freeze_hits": len(freeze_hits),
            "kernel_hits": len(kernel_hits),
            "phi_pi6_hits": len(phi_pi6_hits),
            "phi_zero_hits": len(phi_zero_hits),
            "beta_001_hits": len(beta_001_hits),
            "beta_005_hits": len(beta_005_hits),
            "omega_pi4_hits": len(omega_pi4_hits),
            "nodes_a_hits": len(nodes_a_hits),
            "nodes_b_hits": len(nodes_b_hits),
            "frozen_negative_hits": len(frozen_negative_hits),
            "frozen_positive_hits": len(frozen_positive_hits),
            "check_ok": check_ok,
            "check_fail": check_fail,
        },
        "contradictions": contradictions,
        "contradiction_count": contradiction_count,
        "samples": {
            "freeze_head": freeze_hits[:20],
            "phi_pi6_head": phi_pi6_hits[:20],
            "phi_zero_head": phi_zero_hits[:20],
            "beta_001_head": beta_001_hits[:20],
            "beta_005_head": beta_005_hits[:20],
            "frozen_negative_head": frozen_negative_hits[:20],
            "frozen_positive_head": frozen_positive_hits[:20],
            "section_headers_head": section_headers[:50],
        },
        "verdict": verdict,
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines_out = [
        "# RAPORT QW-1858: FULL_LOG_COMPRESSED_EXTREME_QW420_1200 AUDIT",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Input: `{IN_FILE.name}`",
        f"- Verdict: **{verdict}**",
        "",
        "## Size",
        f"- lines: {n_lines}",
        f"- chars: {n_chars}",
        "",
        "## QW Coverage",
        f"- unique QW ids (any): {out['qw_coverage']['unique_qw_ids_any']}",
        f"- unique QW ids in 420..1200: {out['qw_coverage']['unique_qw_ids_420_1200']}",
        f"- min/max in 420..1200: {out['qw_coverage']['qw_420_1200_min']} / {out['qw_coverage']['qw_420_1200_max']}",
        "",
        "## Marker Counts",
        f"- freeze_hits: {out['markers']['freeze_hits']}",
        f"- kernel_hits: {out['markers']['kernel_hits']}",
        f"- phi_pi6_hits: {out['markers']['phi_pi6_hits']}",
        f"- phi_zero_hits: {out['markers']['phi_zero_hits']}",
        f"- beta_001_hits: {out['markers']['beta_001_hits']}",
        f"- beta_005_hits: {out['markers']['beta_005_hits']}",
        f"- omega_pi4_hits: {out['markers']['omega_pi4_hits']}",
        f"- nodes_a_hits: {out['markers']['nodes_a_hits']}",
        f"- nodes_b_hits: {out['markers']['nodes_b_hits']}",
        f"- frozen_negative_hits: {out['markers']['frozen_negative_hits']}",
        f"- frozen_positive_hits: {out['markers']['frozen_positive_hits']}",
        f"- checks ok/fail: {out['markers']['check_ok']} / {out['markers']['check_fail']}",
        "",
        "## Contradictions",
        f"- phi dual definition: {contradictions['phi_dual_definition']}",
        f"- beta dual definition: {contradictions['beta_dual_definition']}",
        f"- node-set conflict: {contradictions['node_set_conflict']}",
        f"- contradiction_count: {contradiction_count}",
        "",
        "## Notable Lines",
    ]

    for x in freeze_hits[:10]:
        lines_out.append(f"- freeze line {x['line']}: {x['text']}")

    lines_out += [
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]

    OUT_MD.write_text("\n".join(lines_out) + "\n", encoding="utf-8")

    print(f"[QW-1858] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1858] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
