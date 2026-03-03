#!/usr/bin/env python3
"""
QW-1854: Dedicated audit of QW-700..1600 for frozen-kernel narrative.

Purpose:
- quantify actual traceability coverage for QW IDs in [700,1600],
- detect kernel/freeze claim presence in that range,
- detect major internal contradictions (phi=pi/6 vs phi=0, beta=0.01 vs 0.05, etc.).
"""

from __future__ import annotations

import json
import re
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Set


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1854_kernel_freeze_range_audit_700_1600.json"
OUT_MD = ROOT / "RAPORT_QW1854_KERNEL_FREEZE_RANGE_AUDIT_700_1600.md"

TEXT_EXTS = {".md", ".tex", ".json", ".py"}

QW_RE = re.compile(r"QW[-_ ]?(\d{1,4})\b", re.IGNORECASE)

# Kernel/freeze narrative markers
KEYWORDS = [
    "kernel",
    "kernela",
    "nadsoliton",
    "nadsolitona",
    "beta_tors",
    "omega",
    "phi",
    "alpha_geo",
    "frozen",
    "freeze",
    "zamro",
    "zamroż",
    "zamrozone",
    "zamrożone",
]


def is_text_file(p: Path) -> bool:
    return p.is_file() and p.suffix.lower() in TEXT_EXTS


def read_text_safe(p: Path) -> str:
    try:
        return p.read_text(encoding="utf-8", errors="ignore")
    except Exception:
        return ""


def extract_qw_ids(text: str) -> Set[int]:
    out: Set[int] = set()
    for m in QW_RE.finditer(text):
        try:
            q = int(m.group(1))
        except Exception:
            continue
        if 700 <= q <= 1600:
            out.add(q)
    return out


def has_keyword(text_lower: str) -> bool:
    return any(k in text_lower for k in KEYWORDS)


def detect_claim_flags(text_lower: str) -> Dict[str, bool]:
    return {
        "frozen_parameter_set": (
            ("frozen parameter set" in text_lower)
            or ("zamro" in text_lower and "param" in text_lower)
            or ("kernel frozen" in text_lower)
        ),
        "phi_pi6": bool(re.search(r"phi\s*(=|≈|~)?\s*((np\.)?pi\s*/\s*6|0\.5236|π\s*/\s*6)", text_lower)),
        "phi_zero": bool(re.search(r"phi\s*(=|≈|~)?\s*0(\D|$)", text_lower)),
        "omega_pi4": bool(re.search(r"omega\s*(=|≈|~)?\s*((np\.)?pi\s*/\s*4|0\.7854|π\s*/\s*4)", text_lower)),
        "beta_001": bool(re.search(r"beta[_\s-]*tors[^\n\r]{0,40}\b0\.0?1\b", text_lower)),
        "beta_005": bool(re.search(r"beta[_\s-]*tors[^\n\r]{0,40}\b0\.0?5\b", text_lower)),
        "nodes_2_5_8_11": ("2,5,8,11" in text_lower) or ("2 5 8 11" in text_lower),
        "nodes_2_8_14": ("2,8,14" in text_lower) or ("2 8 14" in text_lower),
    }


def main() -> None:
    files = [p for p in ROOT.rglob("*") if is_text_file(p)]

    n_scanned = 0
    files_with_qw_range = []
    files_kernel_in_range = []

    qw_ids_all: Set[int] = set()
    qw_ids_kernel: Set[int] = set()

    claim_counter = Counter()
    claim_files: Dict[str, List[str]] = {}

    for p in files:
        rel = str(p.relative_to(ROOT))

        # keep scope to project files; skip huge binary-like export folders if needed
        if rel.startswith(".git/"):
            continue

        text = read_text_safe(p)
        if not text:
            continue

        n_scanned += 1
        text_lower = text.lower()

        qids = extract_qw_ids(text)
        if not qids:
            continue

        files_with_qw_range.append(rel)
        qw_ids_all.update(qids)

        if has_keyword(text_lower):
            flags = detect_claim_flags(text_lower)
            files_kernel_in_range.append(rel)
            qw_ids_kernel.update(qids)

            for k, v in flags.items():
                if v:
                    claim_counter[k] += 1
                    claim_files.setdefault(k, []).append(rel)

    full_range = set(range(700, 1601))
    coverage_any = len(qw_ids_all) / len(full_range)
    coverage_kernel = len(qw_ids_kernel) / len(full_range)

    contradictions = {
        "phi_dual_definition": claim_counter["phi_pi6"] > 0 and claim_counter["phi_zero"] > 0,
        "beta_dual_definition": claim_counter["beta_001"] > 0 and claim_counter["beta_005"] > 0,
        "node_set_conflict": claim_counter["nodes_2_5_8_11"] > 0 and claim_counter["nodes_2_8_14"] > 0,
    }

    contradiction_count = int(sum(1 for v in contradictions.values() if v))

    # verdict logic
    if coverage_any >= 0.95 and coverage_kernel >= 0.60 and contradiction_count == 0:
        verdict = "RANGE_700_1600_WELL_TRACED_AND_CONSISTENT"
    elif coverage_any >= 0.60:
        verdict = "RANGE_700_1600_PARTIALLY_TRACED_WITH_RISKS"
    else:
        verdict = "RANGE_700_1600_NOT_FULLY_TRACED"

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "scope": {
            "qw_min": 700,
            "qw_max": 1600,
            "range_size": len(full_range),
        },
        "scan_stats": {
            "files_scanned": n_scanned,
            "files_with_qw_range_refs": len(files_with_qw_range),
            "files_kernel_related_in_range": len(files_kernel_in_range),
        },
        "coverage": {
            "unique_qw_ids_with_any_ref": len(qw_ids_all),
            "unique_qw_ids_with_kernel_related_ref": len(qw_ids_kernel),
            "coverage_any_fraction": coverage_any,
            "coverage_kernel_fraction": coverage_kernel,
        },
        "claim_counts": dict(claim_counter),
        "contradictions": contradictions,
        "contradiction_count": contradiction_count,
        "samples": {
            "kernel_related_files_head": files_kernel_in_range[:30],
            "missing_qw_ids_head": sorted(list(full_range - qw_ids_all))[:60],
            "claim_files_head": {k: v[:10] for k, v in claim_files.items()},
        },
        "verdict": verdict,
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1854: KERNEL FREEZE RANGE AUDIT (QW-700..1600)",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        "",
        "## Scan Stats",
        f"- files_scanned: {out['scan_stats']['files_scanned']}",
        f"- files_with_qw_range_refs: {out['scan_stats']['files_with_qw_range_refs']}",
        f"- files_kernel_related_in_range: {out['scan_stats']['files_kernel_related_in_range']}",
        "",
        "## Coverage",
        (
            f"- any-ref coverage: {out['coverage']['unique_qw_ids_with_any_ref']}/"
            f"{out['scope']['range_size']} ({out['coverage']['coverage_any_fraction']:.3f})"
        ),
        (
            f"- kernel-ref coverage: {out['coverage']['unique_qw_ids_with_kernel_related_ref']}/"
            f"{out['scope']['range_size']} ({out['coverage']['coverage_kernel_fraction']:.3f})"
        ),
        "",
        "## Contradictions",
        f"- phi dual definition (pi/6 vs 0): {contradictions['phi_dual_definition']}",
        f"- beta dual definition (0.01 vs 0.05): {contradictions['beta_dual_definition']}",
        f"- node-set conflict (2,5,8,11 vs 2,8,14): {contradictions['node_set_conflict']}",
        f"- contradiction_count: {contradiction_count}",
        "",
        "## Notes",
        "- This is an evidence-range audit, not a proof of each individual QW computation.",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]

    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1854] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1854] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
