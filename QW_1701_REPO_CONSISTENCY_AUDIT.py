#!/usr/bin/env python3
"""
QW-1701: Repository-wide scientific consistency audit.

Purpose:
1) Inventory all research scripts/reports.
2) Map QW IDs across .py/.md/.json artifacts.
3) Detect contradictory verdicts and high-risk claim patterns.
4) Flag methodological risk markers (optimization/fitting/calibration hints).
"""

from __future__ import annotations

import json
import re
from collections import Counter, defaultdict
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Iterable, List, Set


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1701_repo_consistency_audit.json"
OUT_MD = ROOT / "RAPORT_QW1701_REPO_CONSISTENCY_AUDIT.md"

EXCLUDE_DIRS = {".git", ".venv", "__pycache__"}

QW_RE = re.compile(r"\bQW[\s\-_]?V?(\d{1,4})\b", re.IGNORECASE)
REPORT_NUM_RE = re.compile(r"\breport[_\-](\d{1,4})\b", re.IGNORECASE)
RAPORT_QW_RE = re.compile(r"\braport[_\-]?qw[\s\-_]?(\d{1,4})\b", re.IGNORECASE)
PERCENT_RE = re.compile(r"(-?\d+(?:\.\d+)?)\s*%")

SUMMARY_NAME_RE = re.compile(
    r"(summary|podsumowanie|synteza|status|release|report|raport|audit|verify|walkthrough|index|final)",
    re.IGNORECASE,
)

POSITIVE_WORDS = (
    "success",
    "sukces",
    "confirmed",
    "potwierdz",
    "verified",
    "zweryfik",
    "pass",
    "robust",
    "breakthrough",
    "działa",
    "works",
    "udowodn",
    "exact",
)

NEGATIVE_WORDS = (
    "fail",
    "failed",
    "falsified",
    "porażk",
    "obalon",
    "nie działa",
    "inconsistent",
    "limitation",
    "brak",
    "tautolog",
    "numerolog",
    "artifact",
    "error",
    "critical flaw",
)

OPTIMIZATION_HINTS = (
    "optimiz",
    "minimize",
    "curve_fit",
    "least_squares",
    "differential_evolution",
    "basinhopping",
    "nelder",
    "l-bfgs",
    "calibration",
    "fitted",
)

CALIBRATION_HINTS = (
    "experimental",
    "codata",
    "sm_masses",
    "target",
    "reference",
    "ground truth",
    "from data",
)

SM_VALUE_PATTERNS = (
    r"0\.000511",
    r"0\.511e-3",
    r"137\.0\d+",
    r"0\.2312\d*",
    r"91\.188",
    r"80\.379",
    r"1776\.9",
    r"105\.7",
)


@dataclass
class MdAssessment:
    path: str
    mtime_utc: str
    qw_ids: List[int]
    positive_hits: int
    negative_hits: int
    exact_claims: int
    tautology_mentions: int
    fitting_mentions: int
    calibration_mentions: int
    percent_values: List[float]
    is_summary_like: bool
    polarity: str


def iter_files(ext: str) -> Iterable[Path]:
    for p in ROOT.rglob(f"*{ext}"):
        if any(part in EXCLUDE_DIRS for part in p.parts):
            continue
        yield p


def extract_qw_ids(text: str) -> Set[int]:
    out: Set[int] = set()
    for rx in (QW_RE, RAPORT_QW_RE, REPORT_NUM_RE):
        out.update(int(m.group(1)) for m in rx.finditer(text))
    return out


def safe_read_text(path: Path, max_chars: int = 200_000) -> str:
    data = path.read_text(encoding="utf-8", errors="ignore")
    if len(data) > max_chars:
        return data[:max_chars]
    return data


def count_hits(text_lower: str, words: Iterable[str]) -> int:
    return sum(text_lower.count(w) for w in words)


def md_polarity(pos_hits: int, neg_hits: int) -> str:
    if pos_hits > 0 and neg_hits > 0:
        return "mixed"
    if pos_hits > 0:
        return "positive"
    if neg_hits > 0:
        return "negative"
    return "neutral"


def assess_md(path: Path) -> MdAssessment:
    text = safe_read_text(path)
    text_lower = text.lower()
    mtime = datetime.fromtimestamp(path.stat().st_mtime, tz=timezone.utc).isoformat()
    qw_ids = sorted(extract_qw_ids(path.name + "\n" + text))
    pos_hits = count_hits(text_lower, POSITIVE_WORDS)
    neg_hits = count_hits(text_lower, NEGATIVE_WORDS)
    exact_claims = text_lower.count("exact") + text_lower.count("0.00%")
    tautology_mentions = text_lower.count("tautolog")
    fitting_mentions = text_lower.count("fitting") + text_lower.count("fit")
    calibration_mentions = text_lower.count("calibr")
    perc = [float(m.group(1)) for m in PERCENT_RE.finditer(text)]
    return MdAssessment(
        path=str(path.relative_to(ROOT)),
        mtime_utc=mtime,
        qw_ids=qw_ids,
        positive_hits=pos_hits,
        negative_hits=neg_hits,
        exact_claims=exact_claims,
        tautology_mentions=tautology_mentions,
        fitting_mentions=fitting_mentions,
        calibration_mentions=calibration_mentions,
        percent_values=perc[:200],
        is_summary_like=bool(SUMMARY_NAME_RE.search(path.name)),
        polarity=md_polarity(pos_hits, neg_hits),
    )


def assess_py(path: Path) -> Dict[str, object]:
    text = safe_read_text(path)
    text_lower = text.lower()
    mtime = datetime.fromtimestamp(path.stat().st_mtime, tz=timezone.utc).isoformat()
    ids = sorted(extract_qw_ids(path.name + "\n" + text))
    has_opt = any(h in text_lower for h in OPTIMIZATION_HINTS)
    has_calib = any(h in text_lower for h in CALIBRATION_HINTS)
    sm_values = [pat for pat in SM_VALUE_PATTERNS if re.search(pat, text)]
    report_writes = sorted(
        set(
            m.group(1)
            for m in re.finditer(
                r"""['"]([A-Za-z0-9_#\-\s:.()]+(?:report|raport)[A-Za-z0-9_#\-\s:.()]*\.(?:json|md))['"]""",
                text,
                flags=re.IGNORECASE,
            )
        )
    )
    return {
        "path": str(path.relative_to(ROOT)),
        "mtime_utc": mtime,
        "qw_ids": ids,
        "uses_optimization_or_fitting": has_opt,
        "mentions_calibration_or_experiment": has_calib,
        "hardcoded_sm_value_patterns": sm_values,
        "report_targets_in_code": report_writes,
    }


def flatten_json_metrics(obj, prefix="") -> Dict[str, float]:
    out: Dict[str, float] = {}
    if isinstance(obj, dict):
        for k, v in obj.items():
            key = f"{prefix}.{k}" if prefix else str(k)
            out.update(flatten_json_metrics(v, key))
    elif isinstance(obj, list):
        for i, v in enumerate(obj):
            key = f"{prefix}[{i}]"
            out.update(flatten_json_metrics(v, key))
    elif isinstance(obj, (int, float)):
        out[prefix] = float(obj)
    return out


def main() -> None:
    py_files = sorted(iter_files(".py"))
    md_files = sorted(iter_files(".md"))
    json_files = sorted(iter_files(".json"))

    py_assess = [assess_py(p) for p in py_files]
    md_assess = [assess_md(p) for p in md_files]

    qw_to_py: Dict[int, List[str]] = defaultdict(list)
    qw_to_md: Dict[int, List[MdAssessment]] = defaultdict(list)
    qw_to_json: Dict[int, List[str]] = defaultdict(list)

    for rec in py_assess:
        for q in rec["qw_ids"]:
            qw_to_py[q].append(rec["path"])

    for rec in md_assess:
        for q in rec.qw_ids:
            qw_to_md[q].append(rec)

    for p in json_files:
        ids = extract_qw_ids(p.name)
        if not ids:
            # Fallback: try content (small JSON files)
            try:
                ids = extract_qw_ids(safe_read_text(p, max_chars=80_000))
            except Exception:
                ids = set()
        for q in ids:
            qw_to_json[q].append(str(p.relative_to(ROOT)))

    all_qw_ids = sorted(set(qw_to_py) | set(qw_to_md) | set(qw_to_json))

    contradictions = []
    for q in all_qw_ids:
        md_list = qw_to_md.get(q, [])
        polarities = {m.polarity for m in md_list}
        exact_docs = [m.path for m in md_list if m.exact_claims > 0]
        critical_docs = [
            m.path
            for m in md_list
            if (m.negative_hits > 0 or m.tautology_mentions > 0 or m.fitting_mentions > 0)
        ]
        contradictory = (
            ("positive" in polarities and "negative" in polarities)
            or ("mixed" in polarities and ("positive" in polarities or "negative" in polarities))
            or (len(exact_docs) > 0 and len(critical_docs) > 0)
        )
        if contradictory:
            contradictions.append(
                {
                    "qw_id": q,
                    "n_py": len(qw_to_py.get(q, [])),
                    "n_md": len(md_list),
                    "n_json": len(qw_to_json.get(q, [])),
                    "polarities": sorted(polarities),
                    "exact_docs": exact_docs[:8],
                    "critical_docs": critical_docs[:8],
                    "timeline_top_20": sorted(
                        [
                            {
                                "path": m.path,
                                "mtime_utc": m.mtime_utc,
                                "polarity": m.polarity,
                                "exact_claims": m.exact_claims,
                                "negative_hits": m.negative_hits,
                                "tautology_mentions": m.tautology_mentions,
                                "fitting_mentions": m.fitting_mentions,
                            }
                            for m in md_list
                        ],
                        key=lambda x: x["mtime_utc"],
                    )[:20],
                }
            )

    # Repository-level markers
    optimization_scripts = [r for r in py_assess if r["uses_optimization_or_fitting"]]
    calibrated_scripts = [r for r in py_assess if r["mentions_calibration_or_experiment"]]
    sm_hardcoded_scripts = [r for r in py_assess if r["hardcoded_sm_value_patterns"]]

    md_summary_like = [m for m in md_assess if m.is_summary_like]
    md_with_exact = [m for m in md_assess if m.exact_claims > 0]
    md_with_tautology = [m for m in md_assess if m.tautology_mentions > 0]
    md_with_fitting = [m for m in md_assess if m.fitting_mentions > 0]

    # Chronological aggregation by file modification date (UTC month)
    def month_key(iso_ts: str) -> str:
        return iso_ts[:7]

    md_timeline = defaultdict(lambda: Counter())
    for m in md_assess:
        mk = month_key(m.mtime_utc)
        md_timeline[mk]["md_docs"] += 1
        md_timeline[mk]["positive_hits"] += m.positive_hits
        md_timeline[mk]["negative_hits"] += m.negative_hits
        md_timeline[mk]["exact_claims"] += m.exact_claims
        md_timeline[mk]["tautology_mentions"] += m.tautology_mentions
        md_timeline[mk]["fitting_mentions"] += m.fitting_mentions

    py_timeline = defaultdict(lambda: Counter())
    for r in py_assess:
        mk = month_key(r["mtime_utc"])
        py_timeline[mk]["py_scripts"] += 1
        py_timeline[mk]["with_optimization"] += int(r["uses_optimization_or_fitting"])
        py_timeline[mk]["with_calibration"] += int(r["mentions_calibration_or_experiment"])
        py_timeline[mk]["with_hardcoded_sm"] += int(bool(r["hardcoded_sm_value_patterns"]))

    all_times = [datetime.fromtimestamp(p.stat().st_mtime, tz=timezone.utc) for p in (py_files + md_files + json_files)]
    full_date_min = min(all_times).isoformat()
    full_date_max = max(all_times).isoformat()

    per_qw_status = []
    for q in all_qw_ids:
        md_list = qw_to_md.get(q, [])
        pos = sum(m.positive_hits for m in md_list)
        neg = sum(m.negative_hits for m in md_list)
        exact = sum(m.exact_claims for m in md_list)
        taut = sum(m.tautology_mentions for m in md_list)
        fit = sum(m.fitting_mentions for m in md_list)
        per_qw_status.append(
            {
                "qw_id": q,
                "n_py": len(qw_to_py.get(q, [])),
                "n_md": len(md_list),
                "n_json": len(qw_to_json.get(q, [])),
                "positive_hits": pos,
                "negative_hits": neg,
                "exact_claims": exact,
                "tautology_mentions": taut,
                "fitting_mentions": fit,
            }
        )

    per_qw_status.sort(key=lambda x: (-(x["negative_hits"] + x["tautology_mentions"]), -x["n_md"], x["qw_id"]))

    # Parse JSON metrics for quick sanity
    json_metric_samples = []
    for p in json_files[:120]:
        try:
            obj = json.loads(safe_read_text(p))
            flat = flatten_json_metrics(obj)
            metric_subset = {
                k: v
                for k, v in flat.items()
                if any(tok in k.lower() for tok in ("error", "closure", "unit", "corr", "r2", "p_value", "mass", "vev"))
            }
            if metric_subset:
                top_items = sorted(metric_subset.items())[:20]
                json_metric_samples.append(
                    {
                        "path": str(p.relative_to(ROOT)),
                        "n_metrics": len(metric_subset),
                        "sample": top_items,
                    }
                )
        except Exception:
            continue
        if len(json_metric_samples) >= 25:
            break

    coverage_scripts_with_qw = sum(1 for r in py_assess if r["qw_ids"])
    coverage_scripts_with_any_report_target = sum(1 for r in py_assess if r["report_targets_in_code"])
    coverage_scripts_qw_with_mapped_report = sum(
        1
        for r in py_assess
        if r["qw_ids"] and any((q in qw_to_md or q in qw_to_json) for q in r["qw_ids"])
    )

    output = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "root": str(ROOT),
        "date_range_utc": {
            "earliest_file_mtime": full_date_min,
            "latest_file_mtime": full_date_max,
        },
        "file_counts": {
            "py": len(py_files),
            "md": len(md_files),
            "json": len(json_files),
        },
        "coverage": {
            "scripts_with_qw_ids": coverage_scripts_with_qw,
            "scripts_with_report_targets_in_code": coverage_scripts_with_any_report_target,
            "scripts_qw_with_mapped_md_or_json": coverage_scripts_qw_with_mapped_report,
        },
        "methodology_markers": {
            "scripts_with_optimization_or_fitting_markers": len(optimization_scripts),
            "scripts_with_calibration_or_experimental_markers": len(calibrated_scripts),
            "scripts_with_hardcoded_sm_values": len(sm_hardcoded_scripts),
        },
        "md_markers": {
            "summary_like_docs": len(md_summary_like),
            "docs_with_exact_claims": len(md_with_exact),
            "docs_with_tautology_mentions": len(md_with_tautology),
            "docs_with_fitting_mentions": len(md_with_fitting),
        },
        "timeline_by_month": {
            "md": {k: dict(v) for k, v in sorted(md_timeline.items())},
            "py": {k: dict(v) for k, v in sorted(py_timeline.items())},
        },
        "contradictions": {
            "count": len(contradictions),
            "items_top_80": contradictions[:80],
        },
        "high_risk_examples": {
            "hardcoded_sm_scripts_top_30": sm_hardcoded_scripts[:30],
            "optimization_scripts_top_30": optimization_scripts[:30],
            "calibrated_scripts_top_30": calibrated_scripts[:30],
            "exact_and_tautology_docs_top_40": [
                {
                    "path": m.path,
                    "exact_claims": m.exact_claims,
                    "tautology_mentions": m.tautology_mentions,
                    "fitting_mentions": m.fitting_mentions,
                    "polarity": m.polarity,
                }
                for m in md_assess
                if m.exact_claims > 0 and (m.tautology_mentions > 0 or m.fitting_mentions > 0)
            ][:40],
        },
        "per_qw_status_top_200": per_qw_status[:200],
        "json_metric_samples": json_metric_samples,
    }

    OUT_JSON.write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")

    contradiction_ids = [x["qw_id"] for x in contradictions[:20]]
    top_hardcoded = [r["path"] for r in sm_hardcoded_scripts[:20]]
    top_exact_taut = [
        m.path
        for m in md_assess
        if m.exact_claims > 0 and (m.tautology_mentions > 0 or m.fitting_mentions > 0)
    ][:20]

    lines = [
        "# RAPORT QW-1701: REPO CONSISTENCY AUDIT",
        "",
        f"- Data UTC: {output['generated_utc']}",
        f"- Zakres: {len(py_files)} skryptów `.py`, {len(md_files)} dokumentów `.md`, {len(json_files)} raportów `.json`",
        f"- Zakres dat plików (mtime UTC): {full_date_min} → {full_date_max}",
        "",
        "## 1) Pokrycie i mapowanie",
        f"- Skrypty z identyfikatorem QW: {coverage_scripts_with_qw}",
        f"- Skrypty zapisujące raport w kodzie: {coverage_scripts_with_any_report_target}",
        f"- Skrypty QW z powiązanym raportem `.md/.json`: {coverage_scripts_qw_with_mapped_report}",
        "",
        "## 2) Markery metodologiczne",
        f"- Skrypty z markerami optymalizacji/fittingu: {len(optimization_scripts)}",
        f"- Skrypty z markerami kalibracji/odniesień eksperymentalnych: {len(calibrated_scripts)}",
        f"- Skrypty z hardcoded wartościami SM: {len(sm_hardcoded_scripts)}",
        "",
        "## 3) Markery raportowe",
        f"- Dokumenty typu synteza/podsumowanie/status: {len(md_summary_like)}",
        f"- Dokumenty z deklaracjami `EXACT` / `0.00%`: {len(md_with_exact)}",
        f"- Dokumenty z wzmianką o tautologii: {len(md_with_tautology)}",
        f"- Dokumenty z wzmianką o fittingu: {len(md_with_fitting)}",
        "",
        "## 4) Chronologia (wg dat plików)",
        "- Miesięczne trendy (MD):",
    ]
    for month, vals in sorted(md_timeline.items()):
        lines.append(
            f"  - {month}: docs={vals['md_docs']}, pos={vals['positive_hits']}, neg={vals['negative_hits']}, exact={vals['exact_claims']}, taut={vals['tautology_mentions']}, fit={vals['fitting_mentions']}"
        )
    lines.append("- Miesięczne trendy (PY):")
    for month, vals in sorted(py_timeline.items()):
        lines.append(
            f"  - {month}: py={vals['py_scripts']}, opt={vals['with_optimization']}, calib={vals['with_calibration']}, hardcoded_sm={vals['with_hardcoded_sm']}"
        )
    lines.extend(
        [
            "",
            "## 5) Niespójności między raportami",
        f"- Liczba QW z sygnałem sprzeczności (success vs fail, lub exact vs krytyka): {len(contradictions)}",
        f"- Przykładowe QW: {', '.join(f'QW-{q}' for q in contradiction_ids) if contradiction_ids else 'brak'}",
        "",
        "## 6) Wysokie ryzyko metodologiczne (próbka)",
        "- Skrypty z hardcoded SM (próbka):",
        ]
    )
    lines.extend([f"  - `{p}`" for p in top_hardcoded] or ["  - brak"])
    lines.append("- Dokumenty z jednoczesnym `EXACT` i wzmianką `tautolog/fitting` (próbka):")
    lines.extend([f"  - `{p}`" for p in top_exact_taut] or ["  - brak"])
    lines.extend(
        [
            "",
            "## 7) Artefakty",
            f"- JSON szczegółowy: `{OUT_JSON.name}`",
            "",
            "## Wniosek (techniczny)",
            "Repo zawiera równolegle silne deklaracje potwierdzeń i liczne raporty krytyczne/falsyfikacyjne. "
            "W praktyce status teorii jest heterogeniczny i wymaga wersjonowanego rejestru hipotez z rozróżnieniem: "
            "wynik predykcyjny vs wynik dopasowany/diagnostyczny.",
        ]
    )

    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1701] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1701] Saved MD:   {OUT_MD.name}")
    print(f"[QW-1701] Contradictions detected: {len(contradictions)}")


if __name__ == "__main__":
    main()
