#!/usr/bin/env python3
"""
QW-1733: Closure gate for Nadsoliton -> kernel derivation.

Purpose:
1) Aggregate QW-1729, 1730, 1731, 1732 into one closure decision.
2) Separate phenomenological usefulness from first-principles closure.
3) Produce a strict gate result for "kernel derived from characteristics".
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1733_nadsoliton_kernel_closure_gate.json"
OUT_MD = ROOT / "RAPORT_QW1733_NADSOLITON_KERNEL_CLOSURE_GATE.md"


def load_json(path: Path) -> Dict:
    return json.loads(path.read_text(encoding="utf-8"))


def score_from_1731(verdict: str, flags: List[str]) -> float:
    if verdict == "NODE_NARRATIVE_COMPATIBLE":
        return 1.0
    if verdict == "NODE_NARRATIVE_PARTIALLY_INCONSISTENT":
        return 0.4
    # strongly inconsistent
    penalty = min(len(flags), 4) * 0.1
    return max(0.0, 0.2 - penalty)


def score_from_1732(verdict: str, coverage_hi: float) -> float:
    # High-quality fit region width is treated as robustness proxy.
    base = min(max(coverage_hi * 5.0, 0.0), 1.0)
    if verdict == "HYPERBOLIC_REDUCTION_ROBUST":
        return max(base, 0.8)
    if verdict == "HYPERBOLIC_REDUCTION_PLAUSIBLE_BUT_TUNED":
        return min(base, 0.5)
    return min(base, 0.2)


def main() -> None:
    p1729 = ROOT / "report_qw1729_nadsoliton_kernel_characteristics_map.json"
    p1730 = ROOT / "report_qw1730_nadsoliton_kernel_chrono_audit.json"
    p1731 = ROOT / "report_qw1731_nadsoliton_kernel_node_compatibility.json"
    p1732 = ROOT / "report_qw1732_fractal_path_to_hyperbolic_test.json"

    d1729 = load_json(p1729)
    d1730 = load_json(p1730)
    d1731 = load_json(p1731)
    d1732 = load_json(p1732)

    s_1729 = float(d1729.get("kernel_closure_score", 0.0))
    risk_1730 = int(d1730.get("risk_points", 0))
    s_1730 = max(0.0, 1.0 - risk_1730 / 10.0)
    s_1731 = score_from_1731(d1731.get("verdict", ""), d1731.get("incompatibility_flags", []))
    s_1732 = score_from_1732(
        d1732.get("verdict", ""), float(d1732.get("fit_coverage", {}).get("fraction_r2_ge_0_995", 0.0))
    )

    components = [
        {
            "domain": "Mapping Ledger (QW-1729)",
            "score": s_1729,
            "status": "PASS" if s_1729 >= 0.7 else "FAIL",
            "note": d1729.get("verdict", ""),
        },
        {
            "domain": "Chronology Consistency (QW-1730)",
            "score": s_1730,
            "status": "PASS" if s_1730 >= 0.7 else "FAIL",
            "note": d1730.get("verdict", ""),
        },
        {
            "domain": "Node Compatibility (QW-1731)",
            "score": s_1731,
            "status": "PASS" if s_1731 >= 0.7 else "FAIL",
            "note": d1731.get("verdict", ""),
        },
        {
            "domain": "Denominator Robustness (QW-1732)",
            "score": s_1732,
            "status": "PASS" if s_1732 >= 0.7 else "FAIL",
            "note": d1732.get("verdict", ""),
        },
    ]

    global_score = float(sum(c["score"] for c in components) / len(components))
    hard_pass = all(c["status"] == "PASS" for c in components)

    if hard_pass and global_score >= 0.8:
        readiness = "KERNEL_DERIVATION_CLOSED"
    elif global_score >= 0.6:
        readiness = "KERNEL_DERIVATION_PARTIAL"
    else:
        readiness = "KERNEL_DERIVATION_OPEN_NOT_CLOSED"

    output = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "components": components,
        "global_score": global_score,
        "hard_gate": "PASS" if hard_pass else "FAIL",
        "readiness": readiness,
    }
    OUT_JSON.write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1733: NADSOLITON-KERNEL CLOSURE GATE",
        "",
        f"- Data UTC: {output['generated_utc']}",
        f"- Global score: {global_score:.3f}",
        f"- Hard gate: **{output['hard_gate']}**",
        f"- Readiness: **{readiness}**",
        "",
        "## Components",
    ]
    for c in components:
        lines.append(
            f"- {c['domain']}: {c['status']} | score={c['score']:.3f} | note={c['note']}"
        )

    lines.extend(["", "## Artifacts", f"- JSON: `{OUT_JSON.name}`"])
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1733] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1733] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
