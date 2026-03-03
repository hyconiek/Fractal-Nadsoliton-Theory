#!/usr/bin/env python3
"""
QW-1861: Canonical kernel reconstruction (QW-700..1600) under strict evidence filtering.

Combines:
- QW-1856 canonical rows,
- QW-1857 quarantine logic,
- QW-1860 full-log trust catalog,
- file-date weighting (mtime chronology).

Goal:
Recover the most supported kernel tuple from non-contradictory evidence traces.
"""

from __future__ import annotations

import json
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Set, Tuple


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1861_canonical_kernel_reconstruction_700_1600.json"
OUT_MD = ROOT / "RAPORT_QW1861_CANONICAL_KERNEL_RECONSTRUCTION_700_1600.md"

TEXT_EXTS = {".md", ".tex", ".json", ".py"}

QW_RE = re.compile(r"QW[-_ ]?(\d{1,4})\b", re.IGNORECASE)
PHI_PI6_RE = re.compile(r"phi\s*(=|≈|~)?\s*((np\.)?pi\s*/\s*6|0\.5236|π\s*/\s*6)", re.IGNORECASE)
PHI_ZERO_RE = re.compile(r"phi\s*(=|≈|~)?\s*0(\D|$)", re.IGNORECASE)
PHI_NUM_RE = re.compile(r"phi\s*(=|≈|~)?\s*(-?\d+(?:\.\d+)?)", re.IGNORECASE)

BETA_001_RE = re.compile(r"beta[_\s-]*tors[^\n\r]{0,40}\b0\.0?1\b", re.IGNORECASE)
BETA_005_RE = re.compile(r"beta[_\s-]*tors[^\n\r]{0,40}\b0\.0?5\b", re.IGNORECASE)
BETA_NUM_RE = re.compile(r"beta[_\s-]*tors[^\n\r]{0,40}?(-?\d+(?:\.\d+)?)", re.IGNORECASE)

OMEGA_PI4_RE = re.compile(r"omega\s*(=|≈|~)?\s*(np\.pi\s*/\s*4|pi\s*/\s*4|0\.7854)", re.IGNORECASE)
OMEGA_NUM_RE = re.compile(r"omega\s*(=|≈|~)?\s*(-?\d+(?:\.\d+)?)", re.IGNORECASE)

NODES_A_RE = re.compile(r"2\s*,\s*5\s*,\s*8\s*,\s*11|2\s+5\s+8\s+11", re.IGNORECASE)
NODES_B_RE = re.compile(r"2\s*,\s*8\s*,\s*14|2\s+8\s+14", re.IGNORECASE)
NODES_A_FORMULA_RE = re.compile(r"2\s*\+\s*3\s*[nN]\b")


def read_json(name: str) -> Dict:
    p = ROOT / name
    return json.loads(p.read_text(encoding="utf-8"))


def read_text_safe(path: Path) -> str:
    try:
        return path.read_text(encoding="utf-8", errors="ignore")
    except Exception:
        return ""


def is_text_file(path: Path) -> bool:
    return path.is_file() and path.suffix.lower() in TEXT_EXTS


def empty_axis_map() -> Dict[str, Dict[str, float]]:
    return {
        "phi": {"pi_over_6": 0.0, "zero": 0.0, "other": 0.0},
        "beta_tors": {"0.01": 0.0, "0.05": 0.0, "other": 0.0},
        "omega": {"pi_over_4": 0.0, "other_numeric": 0.0},
        "nodes": {"2,5,8,11_or_2+3n": 0.0, "2,8,14": 0.0},
    }


def detect_markers(window_text: str) -> Set[Tuple[str, str]]:
    markers: Set[Tuple[str, str]] = set()

    if PHI_PI6_RE.search(window_text):
        markers.add(("phi", "pi_over_6"))
    if PHI_ZERO_RE.search(window_text):
        markers.add(("phi", "zero"))

    for m in PHI_NUM_RE.finditer(window_text):
        try:
            v = float(m.group(2))
        except Exception:
            continue
        if abs(v - 0.523599) <= 0.05 or abs(v) <= 0.08:
            continue
        markers.add(("phi", "other"))

    if BETA_001_RE.search(window_text):
        markers.add(("beta_tors", "0.01"))
    if BETA_005_RE.search(window_text):
        markers.add(("beta_tors", "0.05"))

    for m in BETA_NUM_RE.finditer(window_text):
        try:
            v = float(m.group(1))
        except Exception:
            continue
        if abs(v - 0.01) <= 0.005 or abs(v - 0.05) <= 0.01:
            continue
        markers.add(("beta_tors", "other"))

    if OMEGA_PI4_RE.search(window_text):
        markers.add(("omega", "pi_over_4"))

    for m in OMEGA_NUM_RE.finditer(window_text):
        try:
            v = float(m.group(2))
        except Exception:
            continue
        if 0.74 <= abs(v) <= 0.83:
            continue
        markers.add(("omega", "other_numeric"))

    if NODES_A_RE.search(window_text) or NODES_A_FORMULA_RE.search(window_text):
        markers.add(("nodes", "2,5,8,11_or_2+3n"))
    if NODES_B_RE.search(window_text):
        markers.add(("nodes", "2,8,14"))

    return markers


def count_markers_in_lines(lines: List[str]) -> Dict[Tuple[str, str], int]:
    counts: Dict[Tuple[str, str], int] = {}
    for line in lines:
        markers = detect_markers(line.lower())
        for m in markers:
            counts[m] = counts.get(m, 0) + 1
    return counts


def axis_decision(axis_vals: Dict[str, float]) -> Dict[str, float]:
    total = float(sum(axis_vals.values()))
    if total <= 0:
        return {
            "winner": "none",
            "winner_weight": 0.0,
            "total_weight": 0.0,
            "consensus": 0.0,
            "runner_up_weight": 0.0,
            "conflict_ratio": 0.0,
        }

    ranked = sorted(axis_vals.items(), key=lambda x: x[1], reverse=True)
    winner, w = ranked[0]
    runner = ranked[1][1] if len(ranked) > 1 else 0.0
    return {
        "winner": winner,
        "winner_weight": w,
        "total_weight": total,
        "consensus": w / total,
        "runner_up_weight": runner,
        "conflict_ratio": runner / total if total > 0 else 0.0,
    }


def main() -> None:
    d1856 = read_json("report_qw1856_kernel_freeze_conflict_ledger_700_1600.json")
    d1857 = read_json("report_qw1857_kernel_freeze_chronology_stress_test_700_1600.json")
    d1860 = read_json("report_qw1860_full_log_trust_catalog.json")

    canonical_qws = sorted(
        {
            int(row["qw"])
            for row in d1856.get("canonical_rows", [])
            if 700 <= int(row["qw"]) <= 1600
        }
    )

    quarantine_qws = {
        int(row["qw"])
        for row in d1857.get("rows", [])
        if row.get("quarantine_from_canonical_freeze")
    }

    input_qw_set = {q for q in canonical_qws if q not in quarantine_qws}

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

    def base_weight(path: Path) -> float:
        name = path.name
        if name in excluded_files:
            return 0.0
        if name in positive_files:
            return 1.0
        if name in negative_files:
            return 0.40
        if name.startswith("FULL_LOG"):
            return 0.65
        return 0.85

    global_axis = empty_axis_map()
    per_q: Dict[int, Dict[str, Dict[str, float]]] = {q: empty_axis_map() for q in input_qw_set}
    file_scores: Dict[str, float] = {}
    evidence_hits: List[Dict] = []

    scanned_files = 0
    used_files = 0

    for path in files:
        scanned_files += 1
        b_w = base_weight(path)
        if b_w <= 0.0:
            continue

        text = read_text_safe(path)
        if not text:
            continue

        lines = text.splitlines()
        if not lines:
            continue

        file_qws = {
            int(m.group(1))
            for m in QW_RE.finditer(text)
            if int(m.group(1)) in input_qw_set
        }
        if not file_qws:
            continue

        c_w = chronology_weight(path)
        final_w = b_w * c_w

        local_hit = False

        for i, line in enumerate(lines):
            q_line = {
                int(m.group(1))
                for m in QW_RE.finditer(line)
                if int(m.group(1)) in input_qw_set
            }
            if not q_line:
                continue

            lo = max(0, i - 2)
            hi = min(len(lines), i + 3)
            window_text = " ".join(lines[lo:hi]).lower()
            markers = detect_markers(window_text)
            if not markers:
                continue

            for q in q_line:
                for axis, value in markers:
                    global_axis[axis][value] += final_w
                    per_q[q][axis][value] += final_w

                    if len(evidence_hits) < 1200:
                        evidence_hits.append(
                            {
                                "qw": q,
                                "file": str(path.relative_to(ROOT)),
                                "line": i + 1,
                                "axis": axis,
                                "value": value,
                                "weight": final_w,
                                "mtime_utc": datetime.fromtimestamp(path.stat().st_mtime, tz=timezone.utc).isoformat(),
                            }
                        )

                    file_scores[str(path.relative_to(ROOT))] = file_scores.get(str(path.relative_to(ROOT)), 0.0) + final_w
                    local_hit = True

        # Diffuse attribution fallback:
        # if local QW-line coupling is absent, still attribute weak evidence from
        # file-level markers across referenced canonical QWs.
        if not local_hit:
            marker_counts = count_markers_in_lines(lines)
            marker_counts = {k: min(v, 5) for k, v in marker_counts.items() if v > 0}
            if marker_counts:
                atten = 0.35
                q_share = atten * final_w / max(1, len(file_qws))
                rel = str(path.relative_to(ROOT))

                for q in file_qws:
                    for (axis, value), cnt in marker_counts.items():
                        w = q_share * float(cnt)
                        global_axis[axis][value] += w
                        per_q[q][axis][value] += w
                        file_scores[rel] = file_scores.get(rel, 0.0) + w

                        if len(evidence_hits) < 1200:
                            evidence_hits.append(
                                {
                                    "qw": q,
                                    "file": rel,
                                    "line": None,
                                    "axis": axis,
                                    "value": value,
                                    "weight": w,
                                    "attribution": "diffuse_file_level",
                                    "marker_count_capped": cnt,
                                    "mtime_utc": datetime.fromtimestamp(path.stat().st_mtime, tz=timezone.utc).isoformat(),
                                }
                            )
                local_hit = True

        if local_hit:
            used_files += 1

    axis_phi = axis_decision(global_axis["phi"])
    axis_beta = axis_decision(global_axis["beta_tors"])
    axis_omega = axis_decision(global_axis["omega"])
    axis_nodes = axis_decision(global_axis["nodes"])

    q_rows = []
    for q in sorted(input_qw_set):
        q_axis = per_q[q]
        q_phi = axis_decision(q_axis["phi"])
        q_beta = axis_decision(q_axis["beta_tors"])
        q_omega = axis_decision(q_axis["omega"])
        q_nodes = axis_decision(q_axis["nodes"])
        q_total = (
            q_phi["total_weight"]
            + q_beta["total_weight"]
            + q_omega["total_weight"]
            + q_nodes["total_weight"]
        )

        q_rows.append(
            {
                "qw": q,
                "total_weight": q_total,
                "phi": q_phi,
                "beta_tors": q_beta,
                "omega": q_omega,
                "nodes": q_nodes,
                "has_any_evidence": q_total > 0.0,
            }
        )

    q_rows_sorted = sorted(q_rows, key=lambda x: x["total_weight"], reverse=True)
    n_qw_with_evidence = sum(1 for x in q_rows if x["has_any_evidence"])

    min_consensus = min(
        axis_phi["consensus"],
        axis_beta["consensus"],
        axis_omega["consensus"],
        axis_nodes["consensus"],
    )

    evidence_cover = (n_qw_with_evidence / len(input_qw_set)) if input_qw_set else 0.0
    max_conflict = max(
        axis_phi["conflict_ratio"],
        axis_beta["conflict_ratio"],
        axis_omega["conflict_ratio"],
        axis_nodes["conflict_ratio"],
    )

    if min_consensus >= 0.80 and evidence_cover >= 0.60:
        if max_conflict <= 0.10:
            verdict = "STRICT_CANONICAL_KERNEL_RECONSTRUCTED"
        else:
            verdict = "CANONICAL_KERNEL_RECONSTRUCTED_WITH_RESIDUAL_CONFLICT"
    elif min_consensus >= 0.60 and evidence_cover >= 0.40:
        verdict = "CANONICAL_KERNEL_PARTIAL_RECONSTRUCTION"
    else:
        verdict = "CANONICAL_KERNEL_RECONSTRUCTION_NOT_CLOSED"

    top_files = [
        {"file": k, "score": v}
        for k, v in sorted(file_scores.items(), key=lambda x: x[1], reverse=True)[:40]
    ]

    evidence_times = sorted(x["mtime_utc"] for x in evidence_hits)
    earliest = evidence_times[0] if evidence_times else None
    latest = evidence_times[-1] if evidence_times else None

    canonical_tuple = {
        "phi": axis_phi["winner"],
        "beta_tors": axis_beta["winner"],
        "omega": axis_omega["winner"],
        "nodes": axis_nodes["winner"],
    }

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "inputs": {
            "canonical_qws_1856": len(canonical_qws),
            "quarantine_qws_1857": len(quarantine_qws),
            "usable_qws": len(input_qw_set),
            "trust_files": {
                "positive": sorted(positive_files),
                "negative": sorted(negative_files),
                "excluded": sorted(excluded_files),
            },
        },
        "scan": {
            "files_scanned": scanned_files,
            "files_used_with_markers": used_files,
            "evidence_time_earliest_utc": earliest,
            "evidence_time_latest_utc": latest,
        },
        "global_axis_weights": global_axis,
        "axis_decisions": {
            "phi": axis_phi,
            "beta_tors": axis_beta,
            "omega": axis_omega,
            "nodes": axis_nodes,
        },
        "coverage": {
            "n_qw_with_evidence": n_qw_with_evidence,
            "n_qw_usable": len(input_qw_set),
            "evidence_coverage_fraction": evidence_cover,
        },
        "consistency": {
            "min_consensus": min_consensus,
            "max_conflict_ratio": max_conflict,
        },
        "canonical_tuple": canonical_tuple,
        "top_files": top_files,
        "q_rows_sorted": q_rows_sorted,
        "evidence_hits_head": evidence_hits[:300],
        "verdict": verdict,
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1861: CANONICAL KERNEL RECONSTRUCTION (QW-700..1600)",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        "",
        "## Inputs",
        f"- canonical QWs from 1856: {out['inputs']['canonical_qws_1856']}",
        f"- quarantined QWs from 1857: {out['inputs']['quarantine_qws_1857']}",
        f"- usable QWs: {out['inputs']['usable_qws']}",
        "",
        "## Scan",
        f"- files scanned: {out['scan']['files_scanned']}",
        f"- files with markers used: {out['scan']['files_used_with_markers']}",
        f"- evidence date span (UTC): {earliest} .. {latest}",
        "",
        "## Axis Decisions",
        f"- phi: winner={axis_phi['winner']} | consensus={axis_phi['consensus']:.3f} | conflict={axis_phi['conflict_ratio']:.3f}",
        f"- beta_tors: winner={axis_beta['winner']} | consensus={axis_beta['consensus']:.3f} | conflict={axis_beta['conflict_ratio']:.3f}",
        f"- omega: winner={axis_omega['winner']} | consensus={axis_omega['consensus']:.3f} | conflict={axis_omega['conflict_ratio']:.3f}",
        f"- nodes: winner={axis_nodes['winner']} | consensus={axis_nodes['consensus']:.3f} | conflict={axis_nodes['conflict_ratio']:.3f}",
        "",
        "## Coverage/Consistency",
        f"- QW evidence coverage: {n_qw_with_evidence}/{len(input_qw_set)} ({evidence_cover:.3f})",
        f"- min consensus: {min_consensus:.3f}",
        f"- max conflict ratio: {max_conflict:.3f}",
        "",
        "## Reconstructed Tuple",
        f"- phi={canonical_tuple['phi']}",
        f"- beta_tors={canonical_tuple['beta_tors']}",
        f"- omega={canonical_tuple['omega']}",
        f"- nodes={canonical_tuple['nodes']}",
        "",
        "## Top Files",
    ]

    for x in top_files[:20]:
        lines.append(f"- {x['file']}: score={x['score']:.2f}")

    lines += ["", "## Top QW Rows"]
    for row in q_rows_sorted[:20]:
        lines.append(
            f"- QW-{row['qw']}: total={row['total_weight']:.2f}, "
            f"phi={row['phi']['winner']}({row['phi']['consensus']:.2f}), "
            f"beta={row['beta_tors']['winner']}({row['beta_tors']['consensus']:.2f}), "
            f"omega={row['omega']['winner']}({row['omega']['consensus']:.2f}), "
            f"nodes={row['nodes']['winner']}({row['nodes']['consensus']:.2f})"
        )

    lines += ["", "## Artifacts", f"- JSON: `{OUT_JSON.name}`"]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1861] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1861] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
