#!/usr/bin/env python3
"""
QW-2006: Deep methodological audit for pre-QW1700 campaign.

Purpose:
1) Reconstruct methodological evolution for QW<1700 from local corpus.
2) Quantify where pipeline is simulation/inference/controls vs heuristic/placeholder risk.
3) Cross-check against key prior audits and TeX-era claims.
"""

from __future__ import annotations

import json
import re
from collections import defaultdict
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional, Tuple

ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw2006_pre1700_methodology_deep_audit.json"
OUT_MD = ROOT / "RAPORT_QW2006_PRE1700_METHODOLOGY_DEEP_AUDIT_EN_PL.md"

ALLOWED_SUFFIXES = {".py", ".md", ".tex"}

BINS: List[Tuple[int, int, str]] = [
    (0, 549, "QW0000_0549_FOUNDATIONAL_HYPOTHESIS"),
    (550, 826, "QW0550_0826_RIGOR_REBOOT_MIXED"),
    (827, 1200, "QW0827_1200_DERIVATION_CAMPAIGN"),
    (1201, 1499, "QW1201_1499_TOPOLOGY_BRIDGE"),
    (1500, 1609, "QW1500_1609_GW_MODELING_PROXY"),
    (1610, 1699, "QW1610_1699_REPAIR_LIMITS_RAW_GW"),
]

QW_RE = re.compile(r"QW[-_ ]?(\d{1,4})\b", re.IGNORECASE)

MARKERS = {
    "numeric_simulation": re.compile(
        r"solve_ivp|scipy\.linalg\.eigh|dijkstra|laplace\(|np\.gradient|fft|spectral|pde|finite[-_ ]difference|simulation",
        re.IGNORECASE,
    ),
    "inference_optimization": re.compile(
        r"likelihood|posterior|bayes|mcmc|minimize\(|curve_fit|polyfit\(|importance sampling|bootstrap",
        re.IGNORECASE,
    ),
    "rigor_controls": re.compile(
        r"null|surrogate|timeslide|blind|oos|out[- ]of[- ]sample|cross[- ]validation|holdout|leakage|prereg|permutation",
        re.IGNORECASE,
    ),
    "external_data": re.compile(
        r"ligo|virgo|gwtc|gw150914|h1\b|l1\b|v1\b|nanograv|pta|raw strain|open data",
        re.IGNORECASE,
    ),
    "frozen_kernel_claim": re.compile(r"frozen|freeze|zamroz", re.IGNORECASE),
    "heuristic_or_placeholder": re.compile(
        r"heuristic|toy model|phenomenology|placeholder|hardcoded|stub|numerology|not derived|not from lagrangian|ansatz",
        re.IGNORECASE,
    ),
    "failure_or_retraction": re.compile(
        r"failed|failure|falsif|retracted|inconclusive|critical caveat|warning|pora[zk]|kryzys",
        re.IGNORECASE,
    ),
    "tuning_fit_risk": re.compile(r"fitting|fit|tuning|calibr|selection bias|retune", re.IGNORECASE),
    "derivation_claim": re.compile(r"derived|proof|q\.e\.d|wyprowadz|theorem", re.IGNORECASE),
}


@dataclass
class FileRow:
    path: str
    qids_pre1700: List[int]
    mtime_utc: str
    marker_counts: Dict[str, int]
    marker_presence: Dict[str, bool]


def read_text(path: Path) -> str:
    try:
        return path.read_text(encoding="utf-8", errors="ignore")
    except Exception:
        return ""


def qids_from_name(name: str) -> List[int]:
    return sorted({int(m.group(1)) for m in QW_RE.finditer(name)})


def qids_from_text_head(text: str, head_lines: int = 300) -> List[int]:
    lines = text.splitlines()[:head_lines]
    head = "\n".join(lines)
    return sorted({int(m.group(1)) for m in QW_RE.finditer(head)})


def file_qids_pre1700(path: Path, text: str) -> List[int]:
    q_name = qids_from_name(path.name)
    q_head = qids_from_text_head(text)
    qids = sorted({q for q in (q_name + q_head) if q < 1700})
    return qids


def bin_name_for_q(q: int) -> Optional[str]:
    for lo, hi, name in BINS:
        if lo <= q <= hi:
            return name
    return None


def classify_file(path: Path, text: str, qids_pre1700: List[int]) -> Optional[FileRow]:
    include_special = path.name.startswith("FULL_LOG") or path.name == "TOE_FINAL_DOCUMENTATION.tex"
    if not qids_pre1700 and not include_special:
        return None

    counts: Dict[str, int] = {}
    presence: Dict[str, bool] = {}
    for key, rgx in MARKERS.items():
        n = len(rgx.findall(text))
        counts[key] = n
        presence[key] = n > 0

    return FileRow(
        path=str(path.relative_to(ROOT)),
        qids_pre1700=qids_pre1700,
        mtime_utc=datetime.fromtimestamp(path.stat().st_mtime, tz=timezone.utc).isoformat(),
        marker_counts=counts,
        marker_presence=presence,
    )


def load_json(path: Path) -> Optional[dict]:
    if not path.exists():
        return None
    try:
        return json.loads(path.read_text(encoding="utf-8", errors="ignore"))
    except Exception:
        return None


def summarize_external_audits() -> Dict[str, dict]:
    refs = {}

    j1703 = load_json(ROOT / "report_qw1703_claims_vs_computation_audit.json")
    if j1703:
        refs["QW1703"] = {
            "verdict_hint": "CLAIMS_VS_COMPUTATION_GAP_PRESENT",
            "exact_mentions": j1703.get("declaration_density", {}).get("exact_or_0_00_mentions_in_core_docs"),
            "fit_mentions": j1703.get("declaration_density", {}).get("fitting_or_calibration_mentions_in_core_docs"),
            "weinberg_err_pct": j1703.get("recomputed_observables", {})
            .get("weinberg_from_alpha_over_12", {})
            .get("error_pct"),
            "alpha_err_pct": j1703.get("recomputed_observables", {})
            .get("fine_structure_inverse_alpha_geo_over_2beta", {})
            .get("error_pct"),
        }

    j1724 = load_json(ROOT / "report_qw1724_gw_method_audit.json")
    if j1724:
        refs["QW1724"] = {
            "verdict": j1724.get("verdict"),
            "risk_points": j1724.get("total_risk_points"),
            "issues_n": len(j1724.get("issues", [])),
        }

    j1858 = load_json(ROOT / "report_qw1858_full_log_extreme_qw420_1200_audit.json")
    if j1858:
        refs["QW1858"] = {
            "verdict": j1858.get("verdict"),
            "contradiction_count": j1858.get("contradiction_count"),
            "frozen_negative_hits": j1858.get("markers", {}).get("frozen_negative_hits"),
        }

    j1907 = load_json(ROOT / "report_qw1907_pre1700_tuning_boundary_audit.json")
    if j1907:
        refs["QW1907"] = {
            "overall": j1907.get("verdict", {}).get("overall"),
            "analysis_tuning_pre1700": j1907.get("verdict", {}).get("analysis_parameter_tuning_pre1700"),
            "kernel_external_retuning_signal": j1907.get("verdict", {}).get("kernel_core_retuning_pre1700_static_signal"),
            "rows_inference_simulation": j1907.get("counts", {}).get("rows_inference_simulation"),
            "rows_inference_external": j1907.get("counts", {}).get("rows_inference_external"),
        }

    j1960 = load_json(ROOT / "report_qw1960_mass_formula_derivation_audit.json")
    if j1960:
        refs["QW1960"] = {
            "verdict": j1960.get("verdict"),
            "gamma_arithmetic_inconsistency": j1960.get("flags", {}).get("arithmetic_inconsistency_gamma_origin"),
            "circularity_risk": j1960.get("flags", {}).get("q_mapping_circularity_risk"),
        }

    return refs


def aggregate_bin(rows: List[FileRow], q_lo: int, q_hi: int) -> dict:
    touched = []
    for r in rows:
        if any(q_lo <= q <= q_hi for q in r.qids_pre1700):
            touched.append(r)

    files_n = len(touched)
    if files_n == 0:
        return {
            "files": 0,
            "marker_file_presence": {},
            "marker_total_hits": {},
            "rigor_index": 0.0,
            "risk_level": "NO_DATA",
        }

    marker_presence_counts = defaultdict(int)
    marker_hit_totals = defaultdict(int)

    for r in touched:
        for k, v in r.marker_presence.items():
            if v:
                marker_presence_counts[k] += 1
        for k, n in r.marker_counts.items():
            marker_hit_totals[k] += int(n)

    pos = (
        marker_presence_counts["numeric_simulation"]
        + marker_presence_counts["inference_optimization"]
        + marker_presence_counts["rigor_controls"]
        + marker_presence_counts["external_data"]
        + marker_presence_counts["failure_or_retraction"]
    )
    neg = (
        marker_presence_counts["heuristic_or_placeholder"]
        + marker_presence_counts["tuning_fit_risk"]
        + marker_presence_counts["frozen_kernel_claim"]
    )

    rigor_index = (pos - neg) / float(pos + neg + 1)

    if rigor_index >= 0.35:
        risk = "LOW_TO_MODERATE"
    elif rigor_index >= 0.0:
        risk = "MODERATE"
    elif rigor_index >= -0.25:
        risk = "MODERATE_TO_HIGH"
    else:
        risk = "HIGH"

    return {
        "files": files_n,
        "marker_file_presence": dict(sorted(marker_presence_counts.items())),
        "marker_total_hits": dict(sorted(marker_hit_totals.items())),
        "rigor_index": rigor_index,
        "risk_level": risk,
    }


def aggregate_by_month(rows: List[FileRow]) -> Dict[str, dict]:
    month_rows: Dict[str, List[FileRow]] = defaultdict(list)
    for r in rows:
        if not r.qids_pre1700:
            continue
        month = r.mtime_utc[:7]
        month_rows[month].append(r)

    out: Dict[str, dict] = {}
    for month, rs in sorted(month_rows.items()):
        marker_presence_counts = defaultdict(int)
        marker_hit_totals = defaultdict(int)
        for r in rs:
            for k, v in r.marker_presence.items():
                if v:
                    marker_presence_counts[k] += 1
            for k, n in r.marker_counts.items():
                marker_hit_totals[k] += int(n)

        pos = (
            marker_presence_counts["numeric_simulation"]
            + marker_presence_counts["inference_optimization"]
            + marker_presence_counts["rigor_controls"]
            + marker_presence_counts["external_data"]
            + marker_presence_counts["failure_or_retraction"]
        )
        neg = (
            marker_presence_counts["heuristic_or_placeholder"]
            + marker_presence_counts["tuning_fit_risk"]
            + marker_presence_counts["frozen_kernel_claim"]
        )
        rigor_index = (pos - neg) / float(pos + neg + 1)

        out[month] = {
            "files": len(rs),
            "rigor_index": rigor_index,
            "marker_file_presence": dict(sorted(marker_presence_counts.items())),
            "marker_total_hits": dict(sorted(marker_hit_totals.items())),
        }
    return out


def evidence_snippets(path: Path, pattern: str, max_items: int = 2) -> List[str]:
    if not path.exists():
        return []
    out = []
    text = read_text(path)
    rgx = re.compile(pattern, re.IGNORECASE)
    for i, line in enumerate(text.splitlines(), start=1):
        if rgx.search(line):
            clipped = line.strip()
            if len(clipped) > 180:
                clipped = clipped[:180] + "..."
            out.append(f"{path.name}:{i}: {clipped}")
            if len(out) >= max_items:
                break
    return out


def main() -> None:
    rows: List[FileRow] = []

    for p in ROOT.rglob("*"):
        if not p.is_file():
            continue
        if p.suffix.lower() not in ALLOWED_SUFFIXES:
            continue
        text = read_text(p)
        if not text:
            continue

        qids_pre1700 = file_qids_pre1700(p, text)
        row = classify_file(p, text, qids_pre1700)
        if row is not None:
            rows.append(row)

    rows_sorted = sorted(rows, key=lambda r: (min(r.qids_pre1700) if r.qids_pre1700 else 9999, r.path))

    bin_summary = {}
    for lo, hi, name in BINS:
        bin_summary[name] = aggregate_bin(rows_sorted, lo, hi)
    month_summary = aggregate_by_month(rows_sorted)

    marker_global_presence = defaultdict(int)
    marker_global_hits = defaultdict(int)
    files_with_q = 0

    for r in rows_sorted:
        if r.qids_pre1700:
            files_with_q += 1
        for k, v in r.marker_presence.items():
            if v:
                marker_global_presence[k] += 1
        for k, n in r.marker_counts.items():
            marker_global_hits[k] += int(n)

    refs = summarize_external_audits()

    evidence = {
        "q1501_toy": evidence_snippets(ROOT / "QW_1501_Generative_Vacuum.py", r"hallucination|pareidolia|random", max_items=3),
        "q1521_rigorous": evidence_snippets(ROOT / "QW_1521_Rigorous_Comparison.py", r"no fitting|discrepancy|rigorous", max_items=3),
        "phase4_disclaimer": evidence_snippets(ROOT / "FULL_LOG_COMPRESSED_PHASE4_QW1611_1620.md", r"METHODOLOGICAL DISCLAIMER|tautological|heuristic", max_items=4),
        "phase6_retraction": evidence_snippets(ROOT / "FULL_LOG_PHASE6.md", r"RETRACTED|circular|inconclusive", max_items=3),
        "tex_claims": evidence_snippets(ROOT / "TOE_FINAL_DOCUMENTATION.tex", r"zero parameter fitting|Chronological Freezing|verified|EXACT", max_items=6),
    }

    # Conservative final verdict logic
    refs_1724_bad = refs.get("QW1724", {}).get("verdict") == "GW_PIPELINE_METHOD_HIGH_RISK_INCONCLUSIVE"
    refs_1960_bad = refs.get("QW1960", {}).get("verdict") == "DERIVATION_CONTAINS_MATERIAL_ERRORS_AND_CIRCULAR_STEPS"
    refs_1907_tuning = refs.get("QW1907", {}).get("analysis_tuning_pre1700") == "DETECTED"

    if refs_1724_bad and refs_1960_bad and refs_1907_tuning:
        verdict = "PRE1700_METHODOLOGY_MULTI_REGIME_PARTIALLY_RIGOROUS_NOT_FULLY_CLOSED"
    else:
        verdict = "PRE1700_METHODOLOGY_REQUIRES_MANUAL_REVIEW"

    output = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "scope": {
            "q_max": 1699,
            "file_suffixes": sorted(ALLOWED_SUFFIXES),
            "rows_total_included": len(rows_sorted),
            "rows_with_explicit_pre1700_qid": files_with_q,
        },
        "global_markers": {
            "file_presence": dict(sorted(marker_global_presence.items())),
            "total_hits": dict(sorted(marker_global_hits.items())),
        },
        "chronology_bins": bin_summary,
        "chronology_months_mtime": month_summary,
        "reference_audits": refs,
        "evidence_snippets": evidence,
        "verdict": verdict,
    }

    OUT_JSON.write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")

    lines: List[str] = []
    lines.append("# RAPORT QW-2006: PRE-QW1700 METHODOLOGY DEEP AUDIT (EN/PL)")
    lines.append("")
    lines.append(f"- Date UTC: {output['generated_utc']}")
    lines.append(f"- Verdict: **{verdict}**")
    lines.append(f"- Scope rows: {len(rows_sorted)} (explicit pre1700 qid: {files_with_q})")
    lines.append("")
    lines.append("## EN: Executive Methodology Verdict")
    lines.append("- Pre-1700 is not one methodologically uniform campaign; it is a sequence of regimes with different rigor levels.")
    lines.append("- Strong numerical and control motifs exist, but they coexist with heuristic/placeholder layers and tuning-heavy segments.")
    lines.append("- Existing audits (QW1703/QW1724/QW1858/QW1907/QW1960) jointly indicate unresolved methodological closure risk.")
    lines.append("")
    lines.append("## PL: Werdykt metodologiczny")
    lines.append("- Badania pre-1700 nie sa jedna jednorodna metodologia; to kilka etapow o roznej jakosci rygoru.")
    lines.append("- Wystepuja realne komponenty numeryczne i kontrolne, ale rownolegle sa warstwy heurystyczne/placeholder oraz odcinki tuningowe.")
    lines.append("- Audyty QW1703/QW1724/QW1858/QW1907/QW1960 razem wskazuja, ze domkniecie metodologiczne nie jest pelne.")
    lines.append("")
    lines.append("## Chronology Bins")
    for lo, hi, name in BINS:
        b = bin_summary[name]
        lines.append(
            f"- {name} ({lo}-{hi}): files={b['files']}, rigor_index={b['rigor_index']:.3f}, risk={b['risk_level']}"
        )
    lines.append("")
    lines.append("## Chronology By File Date (mtime, UTC month)")
    for month, row in month_summary.items():
        lines.append(f"- {month}: files={row['files']}, rigor_index={row['rigor_index']:.3f}")
    lines.append("")
    lines.append("## Cross-Checks To Prior Audits")
    for k in ["QW1703", "QW1724", "QW1858", "QW1907", "QW1960"]:
        if k in refs:
            lines.append(f"- {k}: `{json.dumps(refs[k], ensure_ascii=False)}`")
    lines.append("")
    lines.append("## Evidence Snippets (selected)")
    for block, snippets in evidence.items():
        lines.append(f"- {block}:")
        if not snippets:
            lines.append("  - none")
        else:
            for s in snippets:
                lines.append(f"  - {s}")
    lines.append("")
    lines.append("## Artifacts")
    lines.append(f"- JSON: `{OUT_JSON.name}`")

    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-2006] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-2006] Saved MD:   {OUT_MD.name}")
    print(f"[QW-2006] Verdict: {verdict}")


if __name__ == "__main__":
    main()
