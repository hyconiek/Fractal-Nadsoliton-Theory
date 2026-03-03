#!/usr/bin/env python3
"""
QW-1871: Primary node evidence corpus extraction (chronological, non-circular).

Builds a line-level corpus for node claims from legacy research files while
excluding auto-generated audit/report loops.
"""

from __future__ import annotations

import json
import math
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1871_primary_node_evidence_corpus.json"
OUT_MD = ROOT / "RAPORT_QW1871_PRIMARY_NODE_EVIDENCE_CORPUS.md"

TEXT_EXTS = {".md", ".tex", ".py", ".json", ".txt"}

QW_RE = re.compile(r"QW[-_ ]?(\d{1,4}|V\d{1,3})\b", re.IGNORECASE)

NODE_A_RE = re.compile(r"2\s*,\s*5\s*,\s*8\s*,\s*11|\{\s*2\s*,\s*5\s*,\s*8\s*,\s*11\s*\}|2\s+5\s+8\s+11", re.IGNORECASE)
NODE_B_RE = re.compile(r"2\s*,\s*8\s*,\s*14|\{\s*2\s*,\s*8\s*,\s*14\s*\}|2\s+8\s+14", re.IGNORECASE)
NODE_A_FORMULA_RE = re.compile(r"2\s*\+\s*3\s*[nN]\b")

PHI_PI6_RE = re.compile(r"phi\s*(=|≈|~)?\s*((np\.)?pi\s*/\s*6|0\.5236|π\s*/\s*6)", re.IGNORECASE)
OMEGA_PI4_RE = re.compile(r"omega\s*(=|≈|~)?\s*(np\.pi\s*/\s*4|pi\s*/\s*4|0\.7854)", re.IGNORECASE)
BETA_001_RE = re.compile(r"beta[_\s-]*tors[^\n\r]{0,40}\b0\.0?1\b", re.IGNORECASE)

PHI_ZERO_RE = re.compile(r"phi\s*(=|≈|~)?\s*0(\D|$)", re.IGNORECASE)
BETA_005_RE = re.compile(r"beta[_\s-]*tors[^\n\r]{0,40}\b0\.0?5\b", re.IGNORECASE)

GENERATED_RE = re.compile(
    r"^(report_|raport_|plan_naprawa_toe_roboczy|qw_17\d\d|qw_18\d\d|qw_186\d|qw_187\d)",
    re.IGNORECASE,
)


def is_text_file(path: Path) -> bool:
    parts = {p.lower() for p in path.parts}
    if ".venv" in parts or ".git" in parts or "__pycache__" in parts:
        return False
    return path.is_file() and path.suffix.lower() in TEXT_EXTS


def read_text_safe(path: Path) -> str:
    try:
        return path.read_text(encoding="utf-8", errors="ignore")
    except Exception:
        return ""


def softmax(scores: Dict[str, float]) -> Dict[str, float]:
    if not scores:
        return {}
    m = max(scores.values())
    e = {k: math.exp(v - m) for k, v in scores.items()}
    s = sum(e.values())
    return {k: (v / s if s > 0 else 0.0) for k, v in e.items()}


def file_class(path: Path) -> str:
    name = path.name.lower()
    if GENERATED_RE.match(name):
        return "generated_or_meta"
    if name.startswith("full_log"):
        return "compressed_log"
    if "summary" in name or "podsum" in name or "streszcz" in name:
        return "summary_context"
    return "legacy_primary_candidate"


def main() -> None:
    files = [p for p in ROOT.rglob("*") if is_text_file(p)]

    all_rows = []
    primary_rows = []
    file_rows = []

    primary_files = []

    for path in files:
        rel = str(path.relative_to(ROOT))
        fclass = file_class(path)

        text = read_text_safe(path)
        if not text:
            continue

        lines = text.splitlines()
        if not lines:
            continue

        cnt_a = 0
        cnt_b = 0
        cnt_a_formula = 0

        line_hits = []
        for i, ln in enumerate(lines, start=1):
            c_a = len(NODE_A_RE.findall(ln))
            c_b = len(NODE_B_RE.findall(ln))
            c_af = len(NODE_A_FORMULA_RE.findall(ln))
            if c_a == 0 and c_b == 0 and c_af == 0:
                continue

            cnt_a += c_a
            cnt_b += c_b
            cnt_a_formula += c_af

            row = {
                "file": rel,
                "line": i,
                "class": fclass,
                "text": ln[:240],
                "hits": {
                    "nodes_a": c_a,
                    "nodes_b": c_b,
                    "nodes_a_formula": c_af,
                },
            }
            all_rows.append(row)
            line_hits.append(row)

        if (cnt_a + cnt_b + cnt_a_formula) == 0:
            continue

        low = text.lower()
        has_phi = bool(PHI_PI6_RE.search(low))
        has_omega = bool(OMEGA_PI4_RE.search(low))
        has_beta = bool(BETA_001_RE.search(low))
        has_phi0 = bool(PHI_ZERO_RE.search(low))
        has_beta5 = bool(BETA_005_RE.search(low))

        q_refs = sorted(set(m.group(0) for m in QW_RE.finditer(text)))

        row_f = {
            "file": rel,
            "class": fclass,
            "mtime_utc": datetime.fromtimestamp(path.stat().st_mtime, tz=timezone.utc).isoformat(),
            "node_counts": {
                "a": cnt_a,
                "b": cnt_b,
                "a_formula": cnt_a_formula,
            },
            "canonical_markers": {
                "phi_pi6": has_phi,
                "omega_pi4": has_omega,
                "beta_001": has_beta,
            },
            "contradictory_markers": {
                "phi_zero": has_phi0,
                "beta_005": has_beta5,
            },
            "qw_refs_head": q_refs[:40],
            "line_hits_head": line_hits[:50],
        }

        file_rows.append(row_f)

        if fclass == "legacy_primary_candidate":
            primary_files.append(row_f)
            for h in line_hits:
                primary_rows.append(h)

    primary_files_sorted = sorted(primary_files, key=lambda x: x["mtime_utc"])

    if primary_files_sorted:
        t0 = datetime.fromisoformat(primary_files_sorted[0]["mtime_utc"]).timestamp()
        t1 = datetime.fromisoformat(primary_files_sorted[-1]["mtime_utc"]).timestamp()
    else:
        t0 = 0.0
        t1 = 1.0

    model_support = {
        "M_A_2_5_8_11_or_2plus3n": 0.0,
        "M_B_2_8_14": 0.0,
        "M_MIX_A_and_B": 0.0,
    }

    evidence_weighted = []
    for f in primary_files_sorted:
        ts = datetime.fromisoformat(f["mtime_utc"]).timestamp()
        if t1 <= t0:
            rec = 1.0
        else:
            rec = 0.70 + 0.30 * ((ts - t0) / (t1 - t0))

        pos = sum(1 for v in f["canonical_markers"].values() if v)
        neg = sum(1 for v in f["contradictory_markers"].values() if v)
        quality = max(0.15, 1.0 + 0.30 * pos - 0.45 * neg)

        a = f["node_counts"]["a"] + f["node_counts"]["a_formula"]
        b = f["node_counts"]["b"]
        mix = min(a, b)

        w = rec * quality

        model_support["M_A_2_5_8_11_or_2plus3n"] += w * a
        model_support["M_B_2_8_14"] += w * b
        model_support["M_MIX_A_and_B"] += w * mix

        evidence_weighted.append(
            {
                "file": f["file"],
                "weight": w,
                "a": a,
                "b": b,
                "mix": mix,
                "quality": quality,
                "recency": rec,
                "mtime_utc": f["mtime_utc"],
            }
        )

    scores = {
        "M_A_2_5_8_11_or_2plus3n": model_support["M_A_2_5_8_11_or_2plus3n"] - 0.55 * model_support["M_MIX_A_and_B"],
        "M_B_2_8_14": model_support["M_B_2_8_14"] - 0.55 * model_support["M_MIX_A_and_B"],
        "M_MIX_A_and_B": model_support["M_MIX_A_and_B"],
    }
    posterior = softmax(scores)

    winner = max(posterior, key=posterior.get) if posterior else "none"
    win_p = posterior.get(winner, 0.0)

    if win_p >= 0.75:
        verdict = "PRIMARY_NODE_MODEL_STRONGLY_SELECTED"
    elif win_p >= 0.55:
        verdict = "PRIMARY_NODE_MODEL_PARTIAL_SELECTION"
    else:
        verdict = "PRIMARY_NODE_MODEL_NOT_RESOLVED"

    evidence_weighted_sorted = sorted(
        evidence_weighted,
        key=lambda x: x["weight"] * (x["a"] + x["b"] + 0.5 * x["mix"]),
        reverse=True,
    )

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "summary": {
            "n_files_scanned": len(files),
            "n_files_with_any_node_hit": len(file_rows),
            "n_primary_files_with_node_hit": len(primary_files_sorted),
            "n_primary_line_hits": len(primary_rows),
        },
        "timeline": {
            "primary_earliest_utc": primary_files_sorted[0]["mtime_utc"] if primary_files_sorted else None,
            "primary_latest_utc": primary_files_sorted[-1]["mtime_utc"] if primary_files_sorted else None,
        },
        "model_support_raw": model_support,
        "model_scores": scores,
        "model_posterior": posterior,
        "winner": winner,
        "winner_posterior": win_p,
        "primary_files": primary_files_sorted,
        "top_weighted_evidence": evidence_weighted_sorted[:80],
        "all_file_rows": file_rows,
        "all_line_hits_head": all_rows[:1000],
        "verdict": verdict,
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1871: PRIMARY NODE EVIDENCE CORPUS",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- Winner: **{winner}** (posterior={win_p:.3f})",
        "",
        "## Summary",
        f"- files scanned: {out['summary']['n_files_scanned']}",
        f"- files with any node hit: {out['summary']['n_files_with_any_node_hit']}",
        f"- primary files with node hit: {out['summary']['n_primary_files_with_node_hit']}",
        f"- primary line hits: {out['summary']['n_primary_line_hits']}",
        "",
        "## Timeline",
        f"- earliest primary evidence: {out['timeline']['primary_earliest_utc']}",
        f"- latest primary evidence: {out['timeline']['primary_latest_utc']}",
        "",
        "## Model Scores",
    ]

    for k in ["M_A_2_5_8_11_or_2plus3n", "M_B_2_8_14", "M_MIX_A_and_B"]:
        lines.append(
            f"- {k}: score={scores[k]:.3f}, posterior={posterior.get(k, 0.0):.3f}, raw_support={model_support[k]:.3f}"
        )

    lines += ["", "## Top Weighted Primary Evidence"]
    for x in evidence_weighted_sorted[:20]:
        lines.append(
            f"- {x['file']}: weight={x['weight']:.3f}, a={x['a']}, b={x['b']}, mix={x['mix']}, mtime={x['mtime_utc']}"
        )

    lines += ["", "## Artifacts", f"- JSON: `{OUT_JSON.name}`"]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1871] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1871] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
