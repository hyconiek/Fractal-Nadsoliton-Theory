#!/usr/bin/env python3
"""
QW-1728: Chronology and consistency audit of TOE_FINAL_DOCUMENTATION.tex.

Goal:
1) Extract chronological milestones from the TeX source.
2) Build a date-ordered ledger.
3) Detect internal narrative inconsistencies that impact scientific rigor.
"""

from __future__ import annotations

import json
import re
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parent
TEX = ROOT / "TOE_FINAL_DOCUMENTATION.tex"
OUT_JSON = ROOT / "report_qw1728_tex_chronology_audit.json"
OUT_MD = ROOT / "RAPORT_QW1728_TEX_CHRONOLOGY_AUDIT.md"


def find_release_history(text: str) -> list[dict]:
    # Matches lines like:
    # \item \textbf{v4.4} (2026-01-06): "...."
    pattern = re.compile(
        r"\\item\s+\\textbf\{(?P<ver>v[0-9][0-9A-Za-z\.\-]*)\}\s+\((?P<date>[^)]+)\):\s+(?P<title>.+)"
    )
    out = []
    for m in pattern.finditer(text):
        out.append(
            {
                "version": m.group("ver").strip(),
                "date_raw": m.group("date").strip(),
                "title_raw": m.group("title").strip(),
            }
        )
    return out


def find_phase_blocks(lines: list[str]) -> list[dict]:
    phase_pat = re.compile(r"\\section\{(Phase[^}]*)\}")
    out = []
    for i, ln in enumerate(lines, start=1):
        m = phase_pat.search(ln)
        if m:
            out.append({"line": i, "label": m.group(1).strip()})
    return out


def find_claim_lines(lines: list[str]) -> dict:
    claims = {
        "direct_evidence_claims": [],
        "no_direct_evidence_claims": [],
        "gr_consistent_claims": [],
        "fin_confirmed_claims": [],
    }
    for i, ln in enumerate(lines, start=1):
        low = ln.lower()
        if "first direct empirical evidence" in low or "providing the first direct empirical evidence" in low:
            claims["direct_evidence_claims"].append({"line": i, "text": ln.strip()})
        if "no direct evidence for fin has been recovered" in low:
            claims["no_direct_evidence_claims"].append({"line": i, "text": ln.strip()})
        if "fully consistent with general relativity" in low or "consistent with general relativity" in low:
            claims["gr_consistent_claims"].append({"line": i, "text": ln.strip()})
        if "confirms the presence" in low and "raw gravitational-wave strain" in low:
            claims["fin_confirmed_claims"].append({"line": i, "text": ln.strip()})
    return claims


def parse_date_key(date_raw: str) -> tuple[int, int, int]:
    # Prefer YYYY-MM-DD, fallback to year-only parse.
    m = re.match(r"(\d{4})-(\d{2})-(\d{2})", date_raw)
    if m:
        return int(m.group(1)), int(m.group(2)), int(m.group(3))
    m2 = re.search(r"(20\d{2})", date_raw)
    if m2:
        return int(m2.group(1)), 1, 1
    return (0, 0, 0)


def main() -> None:
    if not TEX.exists():
        raise FileNotFoundError(f"Missing file: {TEX}")

    text = TEX.read_text(encoding="utf-8", errors="ignore")
    lines = text.splitlines()

    releases = find_release_history(text)
    releases_sorted = sorted(releases, key=lambda r: parse_date_key(r["date_raw"]))
    phase_blocks = find_phase_blocks(lines)
    claims = find_claim_lines(lines)

    issues = []

    # Issue A: coexistence of "no direct evidence" and "first direct empirical evidence"
    if claims["no_direct_evidence_claims"] and claims["direct_evidence_claims"]:
        issues.append(
            {
                "severity": "critical",
                "title": "Sprzeczny status dowodu empirycznego",
                "detail": (
                    "Dokument zawiera jednoczesnie tezy: 'no direct evidence for FIN has been recovered' "
                    "oraz 'first direct empirical evidence'."
                ),
                "lines_no_direct": [x["line"] for x in claims["no_direct_evidence_claims"]],
                "lines_direct": [x["line"] for x in claims["direct_evidence_claims"]],
            }
        )

    # Issue B: mixed epistemic status in one document branch (UV-consistent vs confirmed detection)
    if claims["gr_consistent_claims"] and claims["fin_confirmed_claims"]:
        issues.append(
            {
                "severity": "major",
                "title": "Mieszanie poziomow pewnosci (IR-consistency vs confirmed detection)",
                "detail": (
                    "W jednej wersji dokumentu wystepuja jednoczesnie: "
                    "zgodnosc z GR bez detekcji oraz twierdzenie o potwierdzonej detekcji."
                ),
                "lines_gr_consistent": [x["line"] for x in claims["gr_consistent_claims"]],
                "lines_fin_confirmed": [x["line"] for x in claims["fin_confirmed_claims"]],
            }
        )

    # Issue C: release branch overlap (v4.4 claims can conflict with phase 6-15 audit summary)
    if releases_sorted:
        # Check if same version appears multiple times with potentially conflicting titles.
        by_ver = {}
        for r in releases_sorted:
            by_ver.setdefault(r["version"], []).append(r["title_raw"])
        for ver, titles in by_ver.items():
            unique_titles = sorted(set(titles))
            if len(unique_titles) > 1:
                issues.append(
                    {
                        "severity": "minor",
                        "title": f"Wersja {ver} ma wiele roznych opisow",
                        "detail": f"Unikalne opisy: {unique_titles}",
                    }
                )

    risk_points = 0
    sev_w = {"critical": 4, "major": 2, "minor": 1}
    for it in issues:
        risk_points += sev_w.get(it["severity"], 1)

    if risk_points <= 2:
        verdict = "TEX_CHRONOLOGY_COHERENT"
    elif risk_points <= 6:
        verdict = "TEX_CHRONOLOGY_PARTIALLY_INCONSISTENT"
    else:
        verdict = "TEX_CHRONOLOGY_INCONSISTENT_REQUIRES_REFACTOR"

    output = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "source_file": str(TEX),
        "line_count": len(lines),
        "release_history_entries": releases_sorted,
        "phase_sections": phase_blocks,
        "claim_markers": claims,
        "issues": issues,
        "risk_points": risk_points,
        "verdict": verdict,
    }
    OUT_JSON.write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")

    md = [
        "# RAPORT QW-1728: TEX CHRONOLOGY AUDIT",
        "",
        f"- Data UTC: {output['generated_utc']}",
        f"- Plik: `{TEX.name}`",
        f"- Liczba linii: {len(lines)}",
        f"- Risk points: {risk_points}",
        f"- Werdykt: **{verdict}**",
        "",
        "## Release history (ordered)",
    ]
    for r in releases_sorted:
        md.append(f"- {r['date_raw']} | {r['version']} | {r['title_raw']}")

    md.extend(["", "## Sekcje fazowe w dokumencie"])
    for p in phase_blocks:
        md.append(f"- linia {p['line']}: {p['label']}")

    md.extend(["", "## Wykryte niespojnosci"])
    if not issues:
        md.append("- Brak krytycznych niespojnosci w osi czasu.")
    else:
        for i, it in enumerate(issues, start=1):
            md.append(f"- [{i}] ({it['severity']}) {it['title']}: {it['detail']}")

    md.extend(["", "## Artefakty", f"- JSON szczegolowy: `{OUT_JSON.name}`"])
    OUT_MD.write_text("\n".join(md) + "\n", encoding="utf-8")
    print(f"[QW-1728] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1728] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
