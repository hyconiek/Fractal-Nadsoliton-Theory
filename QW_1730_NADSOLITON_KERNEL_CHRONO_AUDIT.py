#!/usr/bin/env python3
"""
QW-1730: Chronological audit of Nadsoliton -> kernel claims.

Purpose:
1) Read md/tex corpus and track when core kernel-origin claims appear.
2) Detect direct internal contradictions (phi, beta_tors, node pattern).
3) Quantify chronology-based scientific risk for kernel derivation narrative.
"""

from __future__ import annotations

import json
import math
import re
from collections import defaultdict
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Pattern


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1730_nadsoliton_kernel_chrono_audit.json"
OUT_MD = ROOT / "RAPORT_QW1730_NADSOLITON_KERNEL_CHRONO_AUDIT.md"


@dataclass(frozen=True)
class ClaimPattern:
    claim_id: str
    description: str
    regex: Pattern[str]


CLAIMS: List[ClaimPattern] = [
    ClaimPattern(
        claim_id="ALPHA_4LN2",
        description="alpha_geo tied to 4*ln(2) / 2.7725",
        regex=re.compile(
            r"(4\s*\\?ln\s*2|4\s*\*?\s*ln\s*2|4ln2|alpha[_\s]*geo\s*=\s*2\.772|2\.772588)",
            flags=re.IGNORECASE,
        ),
    ),
    ClaimPattern(
        claim_id="OMEGA_PI4",
        description="omega tied to pi/4 or 2*pi/8",
        regex=re.compile(
            r"(omega\s*=\s*(?:np\.)?pi\s*/\s*4|omega\s*=\s*2\s*\*\s*(?:np\.)?pi\s*/\s*8|0\.785398)",
            flags=re.IGNORECASE,
        ),
    ),
    ClaimPattern(
        claim_id="PHI_PI6",
        description="phi tied to pi/6 / 0.5236",
        regex=re.compile(
            r"(phi\s*=\s*(?:np\.)?pi\s*/\s*6|phi\s*=\s*0\.5236|pi\s*/\s*6)",
            flags=re.IGNORECASE,
        ),
    ),
    ClaimPattern(
        claim_id="PHI_ZERO",
        description="phi fixed to 0",
        regex=re.compile(
            r"(\bphi\s*=\s*0(?:\.0+)?\b|\bvarphi\s*=\s*0(?:\.0+)?\b)",
            flags=re.IGNORECASE,
        ),
    ),
    ClaimPattern(
        claim_id="BETA_001",
        description="beta_tors fixed near 0.01",
        regex=re.compile(
            r"(beta[_\s]*tors\s*=\s*0\.01|\bbeta\s*=\s*0\.01\b|0\.01\s*\(tuned from topology\))",
            flags=re.IGNORECASE,
        ),
    ),
    ClaimPattern(
        claim_id="BETA_005",
        description="beta_tors fixed near 0.05",
        regex=re.compile(
            r"(beta[_\s]*tors\s*=\s*0\.05|\bbeta\s*=\s*0\.05\b|1\s*/\s*\(1\s*\+\s*0\.05d\))",
            flags=re.IGNORECASE,
        ),
    ),
    ClaimPattern(
        claim_id="NODES_25811",
        description="node set d=2,5,8,11",
        regex=re.compile(
            r"(d\s*=\s*2\s*,\s*5\s*,\s*8\s*,\s*11|nodes?\s+at\s+d\s*=\s*2\s*,\s*5\s*,\s*8\s*,\s*11)",
            flags=re.IGNORECASE,
        ),
    ),
    ClaimPattern(
        claim_id="NODES_2814",
        description="node set d=2,8,14",
        regex=re.compile(
            r"(2\s*,\s*8\s*,\s*14|nodes?\s+at\s+d\s*=\s*2\s*,\s*8\s*,\s*14)",
            flags=re.IGNORECASE,
        ),
    ),
    ClaimPattern(
        claim_id="HYPERBOLIC_DENOM",
        description="hyperbolic damping 1/(1+beta*d)",
        regex=re.compile(
            r"(1\s*/\s*\(\s*1\s*\+\s*beta|\(1\s*\+\s*beta[_\s]*tors\s*\*\s*d\)\s*)",
            flags=re.IGNORECASE,
        ),
    ),
    ClaimPattern(
        claim_id="PATH_SUM_DERIVATION",
        description="path/tunneling derivation of denominator",
        regex=re.compile(
            r"(path\s+integral|topological\s+tunneling|path\s+summation|many\s+paths)",
            flags=re.IGNORECASE,
        ),
    ),
]


def utc_iso_from_mtime(path: Path) -> str:
    return datetime.fromtimestamp(path.stat().st_mtime, tz=timezone.utc).isoformat()


def iter_relevant_files(root: Path) -> List[Path]:
    files: List[Path] = []
    for ext in ("*.md", "*.tex"):
        files.extend(root.rglob(ext))
    files = [p for p in files if p.is_file()]
    return sorted(files)


def detect_claims(text: str) -> List[str]:
    found = []
    for claim in CLAIMS:
        if claim.regex.search(text):
            found.append(claim.claim_id)
    return found


def summarize_claim_hits(records: List[Dict[str, object]]) -> Dict[str, Dict[str, object]]:
    by_claim: Dict[str, List[Dict[str, object]]] = defaultdict(list)
    for rec in records:
        for claim_id in rec["claims"]:
            by_claim[claim_id].append(rec)

    summary: Dict[str, Dict[str, object]] = {}
    for claim in CLAIMS:
        hits = sorted(by_claim.get(claim.claim_id, []), key=lambda r: r["mtime_iso"])
        if not hits:
            summary[claim.claim_id] = {
                "description": claim.description,
                "count_files": 0,
                "first_seen": None,
                "last_seen": None,
            }
            continue
        summary[claim.claim_id] = {
            "description": claim.description,
            "count_files": len(hits),
            "first_seen": {"mtime_iso": hits[0]["mtime_iso"], "path": hits[0]["path"]},
            "last_seen": {"mtime_iso": hits[-1]["mtime_iso"], "path": hits[-1]["path"]},
        }
    return summary


def build_contradictions(claim_summary: Dict[str, Dict[str, object]]) -> List[Dict[str, object]]:
    contradictions: List[Dict[str, object]] = []

    phi_pi6 = claim_summary["PHI_PI6"]["count_files"] > 0
    phi_zero = claim_summary["PHI_ZERO"]["count_files"] > 0
    if phi_pi6 and phi_zero:
        contradictions.append(
            {
                "id": "PHI_DUAL_DEFINITION",
                "severity": "high",
                "detail": "Both phi=pi/6 and phi=0 are present in core corpus.",
            }
        )

    beta_001 = claim_summary["BETA_001"]["count_files"] > 0
    beta_005 = claim_summary["BETA_005"]["count_files"] > 0
    if beta_001 and beta_005:
        contradictions.append(
            {
                "id": "BETA_DUAL_DEFINITION",
                "severity": "high",
                "detail": "Both beta_tors=0.01 and beta_tors=0.05 appear as fixed baselines.",
            }
        )

    nodes_25811 = claim_summary["NODES_25811"]["count_files"] > 0
    nodes_2814 = claim_summary["NODES_2814"]["count_files"] > 0
    if nodes_25811 and nodes_2814:
        contradictions.append(
            {
                "id": "NODE_PATTERN_CONFLICT",
                "severity": "high",
                "detail": "At least two incompatible node sets are simultaneously present.",
            }
        )

    # Internal arithmetic compatibility check for canonical parameters:
    # nodes solve cos(omega*d + phi)=0 => d = (pi/2 - phi + n*pi)/omega.
    omega = math.pi / 4.0
    phi = math.pi / 6.0
    first_node = (math.pi / 2.0 - phi) / omega
    spacing = math.pi / omega
    # Claimed pattern in diagrams: first node near d=2 and spacing 3.
    if abs(first_node - 2.0) > 0.4 or abs(spacing - 3.0) > 0.4:
        contradictions.append(
            {
                "id": "NODE_FORMULA_ARITHMETIC_MISMATCH",
                "severity": "critical",
                "detail": (
                    "For omega=pi/4, phi=pi/6, expected node start/spacings are "
                    f"d0={first_node:.3f}, delta_d={spacing:.3f}, not (2, 3-step)."
                ),
            }
        )

    hyperbolic = claim_summary["HYPERBOLIC_DENOM"]["count_files"] > 0
    path_deriv = claim_summary["PATH_SUM_DERIVATION"]["count_files"] > 0
    if hyperbolic and not path_deriv:
        contradictions.append(
            {
                "id": "HYPERBOLIC_NO_DERIVATION_TRACE",
                "severity": "medium",
                "detail": "Hyperbolic denominator appears without matching derivation traces in same corpus.",
            }
        )

    return contradictions


def compute_risk_points(contradictions: List[Dict[str, object]]) -> int:
    score_by_severity = {"critical": 3, "high": 2, "medium": 1, "low": 1}
    return int(sum(score_by_severity.get(c["severity"], 1) for c in contradictions))


def main() -> None:
    files = iter_relevant_files(ROOT)
    records: List[Dict[str, object]] = []

    for path in files:
        try:
            text = path.read_text(encoding="utf-8", errors="ignore")
        except OSError:
            continue
        claims = detect_claims(text)
        if not claims:
            continue
        records.append(
            {
                "path": str(path.relative_to(ROOT)),
                "mtime_iso": utc_iso_from_mtime(path),
                "claims": claims,
            }
        )

    claim_summary = summarize_claim_hits(records)
    contradictions = build_contradictions(claim_summary)
    risk_points = compute_risk_points(contradictions)

    if risk_points >= 8:
        verdict = "KERNEL_ORIGIN_CHRONOLOGY_HIGH_RISK_INCONSISTENT"
    elif risk_points >= 4:
        verdict = "KERNEL_ORIGIN_CHRONOLOGY_PARTIAL_INCONSISTENCY"
    else:
        verdict = "KERNEL_ORIGIN_CHRONOLOGY_ACCEPTABLE"

    records_sorted = sorted(records, key=lambda r: r["mtime_iso"])
    chronology_head = records_sorted[:20]
    chronology_tail = records_sorted[-20:]

    output = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "scope": {
            "files_scanned_total": len(files),
            "files_with_relevant_claims": len(records),
            "extensions": [".md", ".tex"],
        },
        "claim_summary": claim_summary,
        "contradictions": contradictions,
        "risk_points": risk_points,
        "verdict": verdict,
        "chronology_examples": {
            "earliest_20": chronology_head,
            "latest_20": chronology_tail,
        },
    }
    OUT_JSON.write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")

    lines: List[str] = [
        "# RAPORT QW-1730: NADSOLITON-KERNEL CHRONO AUDIT",
        "",
        f"- Data UTC: {output['generated_utc']}",
        f"- Files scanned: {output['scope']['files_scanned_total']}",
        f"- Files with kernel-origin claims: {output['scope']['files_with_relevant_claims']}",
        f"- Risk points: {risk_points}",
        f"- Verdict: **{verdict}**",
        "",
        "## Contradictions",
    ]
    if contradictions:
        for c in contradictions:
            lines.append(f"- [{c['severity']}] {c['id']}: {c['detail']}")
    else:
        lines.append("- None detected at this scan level.")

    lines.extend(["", "## Claim Presence"])
    for claim in CLAIMS:
        s = claim_summary[claim.claim_id]
        lines.append(f"- {claim.claim_id}: count={s['count_files']} | {claim.description}")

    lines.extend(
        [
            "",
            "## Chronology Sample (Earliest Relevant Files)",
        ]
    )
    for rec in chronology_head[:10]:
        lines.append(f"- {rec['mtime_iso']} | {rec['path']} | claims={','.join(rec['claims'])}")

    lines.extend(
        [
            "",
            "## Chronology Sample (Latest Relevant Files)",
        ]
    )
    for rec in chronology_tail[-10:]:
        lines.append(f"- {rec['mtime_iso']} | {rec['path']} | claims={','.join(rec['claims'])}")

    lines.extend(["", "## Artifacts", f"- JSON: `{OUT_JSON.name}`"])
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1730] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1730] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
