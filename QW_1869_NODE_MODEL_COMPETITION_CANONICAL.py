#!/usr/bin/env python3
"""
QW-1869: Node-model competition under canonical kernel constraints.

Compares candidate node narratives in QW-700..1600:
- M_A: 2,5,8,11 and/or 2+3n
- M_B: 2,8,14
- M_MIX: explicit coexistence of A and B markers
"""

from __future__ import annotations

import json
import math
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1869_node_model_competition_canonical.json"
OUT_MD = ROOT / "RAPORT_QW1869_NODE_MODEL_COMPETITION_CANONICAL.md"

TEXT_EXTS = {".md", ".tex", ".json", ".py"}
QW_RE = re.compile(r"QW[-_ ]?(\d{1,4})\b", re.IGNORECASE)

NODES_A_RE = re.compile(r"2\s*,\s*5\s*,\s*8\s*,\s*11|2\s+5\s+8\s+11", re.IGNORECASE)
NODES_A_FORMULA_RE = re.compile(r"2\s*\+\s*3\s*[nN]\b")
NODES_B_RE = re.compile(r"2\s*,\s*8\s*,\s*14|2\s+8\s+14", re.IGNORECASE)

PHI_PI6_RE = re.compile(r"phi\s*(=|≈|~)?\s*((np\.)?pi\s*/\s*6|0\.5236|π\s*/\s*6)", re.IGNORECASE)
PHI_ZERO_RE = re.compile(r"phi\s*(=|≈|~)?\s*0(\D|$)", re.IGNORECASE)
OMEGA_PI4_RE = re.compile(r"omega\s*(=|≈|~)?\s*(np\.pi\s*/\s*4|pi\s*/\s*4|0\.7854)", re.IGNORECASE)
BETA_001_RE = re.compile(r"beta[_\s-]*tors[^\n\r]{0,40}\b0\.0?1\b", re.IGNORECASE)
BETA_005_RE = re.compile(r"beta[_\s-]*tors[^\n\r]{0,40}\b0\.0?5\b", re.IGNORECASE)


def read_json(name: str) -> Dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def read_text_safe(path: Path) -> str:
    try:
        return path.read_text(encoding="utf-8", errors="ignore")
    except Exception:
        return ""


def is_text_file(path: Path) -> bool:
    return path.is_file() and path.suffix.lower() in TEXT_EXTS


def softmax(d: Dict[str, float]) -> Dict[str, float]:
    if not d:
        return {}
    m = max(d.values())
    ex = {k: math.exp(v - m) for k, v in d.items()}
    s = sum(ex.values())
    return {k: (v / s if s > 0 else 0.0) for k, v in ex.items()}


def main() -> None:
    d1860 = read_json("report_qw1860_full_log_trust_catalog.json")

    positive_files = {x["file"] for x in d1860.get("positive_context", [])}
    negative_files = {x["file"] for x in d1860.get("negative_evidence", [])}
    excluded_files = {x["file"] for x in d1860.get("excluded_contradictory", [])}

    files = [p for p in ROOT.rglob("*") if is_text_file(p)]
    mtimes = [p.stat().st_mtime for p in files]
    mt_min = min(mtimes) if mtimes else 0.0
    mt_max = max(mtimes) if mtimes else 1.0

    def chronology_weight(path: Path) -> float:
        if mt_max <= mt_min:
            return 1.0
        frac = (path.stat().st_mtime - mt_min) / (mt_max - mt_min)
        return 0.70 + 0.30 * max(0.0, min(1.0, frac))

    def trust_weight(path: Path) -> float:
        name = path.name
        if name in positive_files:
            return 1.00
        if name in negative_files:
            return 0.40
        if name in excluded_files:
            return 0.20
        if name.startswith("FULL_LOG"):
            return 0.65
        return 0.85

    model_stats = {
        "M_A_2_5_8_11_or_2plus3n": {"base": 0.0, "coherence": 0.0, "contradiction": 0.0},
        "M_B_2_8_14": {"base": 0.0, "coherence": 0.0, "contradiction": 0.0},
        "M_MIX_A_and_B": {"base": 0.0, "coherence": 0.0, "contradiction": 0.0},
    }

    top_files = {k: [] for k in model_stats.keys()}

    scanned = 0
    used = 0

    for path in files:
        scanned += 1
        name_low = path.name.lower()
        # Exclude auto-generated/meta artifacts to avoid circular evidence loops.
        if (
            name_low.startswith("report_qw")
            or name_low.startswith("raport_qw")
            or name_low.startswith("qw_")
            or name_low.startswith("plan_naprawa")
        ):
            continue

        text = read_text_safe(path)
        if not text:
            continue

        qids = {int(m.group(1)) for m in QW_RE.finditer(text) if 700 <= int(m.group(1)) <= 1600}
        if not qids:
            continue

        low = text.lower()

        cnt_a = len(NODES_A_RE.findall(text)) + len(NODES_A_FORMULA_RE.findall(text))
        cnt_b = len(NODES_B_RE.findall(text))
        has_a = cnt_a > 0
        has_b = cnt_b > 0
        has_mix = has_a and has_b

        if not (has_a or has_b):
            continue

        used += 1

        has_phi = bool(PHI_PI6_RE.search(low))
        has_omega = bool(OMEGA_PI4_RE.search(low))
        has_beta = bool(BETA_001_RE.search(low))
        has_phi_zero = bool(PHI_ZERO_RE.search(low))
        has_beta_005 = bool(BETA_005_RE.search(low))

        canonical_hits = int(has_phi) + int(has_omega) + int(has_beta)
        contradiction_hits = int(has_phi_zero) + int(has_beta_005)

        w = trust_weight(path) * chronology_weight(path)
        rel = str(path.relative_to(ROOT))

        def add(model: str, mention_count: int, extra_contra: int) -> None:
            base = w * (1.0 + math.log1p(max(0, mention_count)))
            coh = w * float(canonical_hits)
            contra = w * float(contradiction_hits + extra_contra)
            net = base + 1.30 * coh - 1.40 * contra

            model_stats[model]["base"] += base
            model_stats[model]["coherence"] += coh
            model_stats[model]["contradiction"] += contra

            top_files[model].append(
                {
                    "file": rel,
                    "qw_count_in_range": len(qids),
                    "weight": w,
                    "base": base,
                    "coherence": coh,
                    "contradiction": contra,
                    "net": net,
                    "mtime_utc": datetime.fromtimestamp(path.stat().st_mtime, tz=timezone.utc).isoformat(),
                }
            )

        if has_a:
            add("M_A_2_5_8_11_or_2plus3n", cnt_a, 1 if has_mix else 0)
        if has_b:
            add("M_B_2_8_14", cnt_b, 1 if has_mix else 0)
        if has_mix:
            add("M_MIX_A_and_B", min(cnt_a, cnt_b), 0)

    scores = {}
    for k, v in model_stats.items():
        scores[k] = v["base"] + 1.30 * v["coherence"] - 1.40 * v["contradiction"]

    post = softmax(scores)
    winner = max(post, key=post.get) if post else "none"
    winner_p = post.get(winner, 0.0)

    if winner_p >= 0.75 and scores.get(winner, 0.0) > 0:
        verdict = "NODE_MODEL_STRONGLY_SELECTED"
    elif winner_p >= 0.55:
        verdict = "NODE_MODEL_PARTIAL_SELECTION"
    else:
        verdict = "NODE_MODEL_NOT_RESOLVED"

    top_files_sorted = {
        k: sorted(v, key=lambda x: x["net"], reverse=True)[:30]
        for k, v in top_files.items()
    }

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "scan": {
            "files_scanned": scanned,
            "files_used_with_node_markers": used,
        },
        "model_stats": model_stats,
        "scores": scores,
        "posterior": post,
        "winner": winner,
        "winner_posterior": winner_p,
        "top_files": top_files_sorted,
        "verdict": verdict,
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1869: NODE MODEL COMPETITION (CANONICAL)",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- Winner: **{winner}** (posterior={winner_p:.3f})",
        "",
        "## Scan",
        f"- files scanned: {scanned}",
        f"- files used with node markers: {used}",
        "",
        "## Model Scores",
    ]

    for k in model_stats.keys():
        lines.append(
            f"- {k}: score={scores[k]:.3f}, posterior={post.get(k, 0.0):.3f}, "
            f"base={model_stats[k]['base']:.3f}, coherence={model_stats[k]['coherence']:.3f}, contradiction={model_stats[k]['contradiction']:.3f}"
        )

    lines += ["", "## Top Evidence Files (winner)"]
    for row in top_files_sorted.get(winner, [])[:15]:
        lines.append(
            f"- {row['file']}: net={row['net']:.3f}, coherence={row['coherence']:.3f}, contradiction={row['contradiction']:.3f}"
        )

    lines += ["", "## Artifacts", f"- JSON: `{OUT_JSON.name}`"]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1869] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1869] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
