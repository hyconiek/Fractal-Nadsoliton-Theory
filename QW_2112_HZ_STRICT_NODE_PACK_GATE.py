#!/usr/bin/env python3
"""
QW-2112: strict H(z) candidate-pack gate for T3 closure.

Purpose:
- validate externally collected H(z) candidate nodes with per-node provenance,
- merge with baseline nodes in a non-destructive way,
- check strict-ready identifiability thresholds before downstream gates.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np


ROOT = Path(__file__).resolve().parent

BASE_CSV_DEFAULT = ROOT / "external_hz_nodes_qw2099.csv"
CANDIDATE_CSV_DEFAULT = ROOT / "external_hz_nodes_qw2112_candidates.csv"
MERGED_OUT_CSV = ROOT / "external_hz_nodes_qw2099.strict_candidate_merged.csv"
MERGED_OUT_META = ROOT / "external_hz_nodes_qw2099.strict_candidate_merged.metadata.json"
TEMPLATE_CSV = ROOT / "external_hz_nodes_qw2112_candidates.template.csv"

OUT_JSON = ROOT / "report_qw2112_hz_strict_node_pack_gate.json"
OUT_MD = ROOT / "RAPORT_QW2112_HZ_STRICT_NODE_PACK_GATE.md"


def _is_https_url(url: str) -> bool:
    return str(url).strip().lower().startswith("https://")


def _read_rows(path: Path, require_provenance: bool) -> Tuple[List[Dict], str]:
    needed = {"z", "h_km_s_mpc", "sigma_total"}
    if require_provenance:
        needed = needed | {"citation", "reference_url", "source_version"}
    rows: List[Dict] = []
    with path.open("r", encoding="utf-8", newline="") as f:
        r = csv.DictReader(f)
        header = set(r.fieldnames or [])
        missing = sorted(needed - header)
        if missing:
            return [], f"Missing columns: {missing}"
        for i, row in enumerate(r, start=2):
            try:
                z = float(row["z"])
                h = float(row["h_km_s_mpc"])
                s = float(row["sigma_total"])
            except Exception as exc:  # noqa: BLE001
                return [], f"Numeric parse error at line {i}: {exc}"
            if not (math.isfinite(z) and math.isfinite(h) and math.isfinite(s)):
                return [], f"Non-finite value at line {i}"
            if h <= 0.0 or s <= 0.0:
                return [], f"Non-positive h/sigma at line {i}"

            rec = {"z": z, "h_km_s_mpc": h, "sigma_total": s}
            if require_provenance:
                rec["citation"] = str(row["citation"]).strip()
                rec["reference_url"] = str(row["reference_url"]).strip()
                rec["source_version"] = str(row["source_version"]).strip()
                rec["source_dataset_name"] = str(row.get("source_dataset_name", "")).strip()
            rows.append(rec)
    return rows, ""


def _write_template_if_missing(path: Path) -> None:
    if path.exists():
        return
    with path.open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(
            [
                "z",
                "h_km_s_mpc",
                "sigma_total",
                "citation",
                "reference_url",
                "source_version",
                "source_dataset_name",
            ]
        )


def _metrics(rows: List[Dict], omega_m: float = 0.315, omega_r: float = 9e-5) -> Dict[str, float]:
    if not rows:
        return {
            "n_nodes": 0,
            "z_span": float("nan"),
            "e_span": float("nan"),
            "design_condition_number": float("nan"),
        }
    z = np.array([r["z"] for r in rows], dtype=float)
    e = np.array([omega_m * (1.0 + zz) ** 3 + omega_r * (1.0 + zz) ** 4 for zz in z], dtype=float)
    a = np.column_stack((e, np.ones_like(e)))
    return {
        "n_nodes": int(len(rows)),
        "z_span": float(np.max(z) - np.min(z)),
        "e_span": float(np.max(e) - np.min(e)),
        "design_condition_number": float(np.linalg.cond(a)),
    }


def _merge_unique(base: List[Dict], add: List[Dict], eps: float = 1e-9) -> Tuple[List[Dict], int]:
    out = [dict(r) for r in base]
    added_count = 0
    for cand in add:
        zc = float(cand["z"])
        if any(abs(float(r["z"]) - zc) <= eps for r in out):
            continue
        out.append(dict(cand))
        added_count += 1
    out.sort(key=lambda r: float(r["z"]))
    return out, added_count


def main() -> None:
    p = argparse.ArgumentParser(description="QW-2112 H(z) strict candidate-pack gate")
    p.add_argument("--base-csv", default=str(BASE_CSV_DEFAULT), help="Baseline H(z) CSV.")
    p.add_argument(
        "--candidate-csv",
        default=str(CANDIDATE_CSV_DEFAULT),
        help="Candidate external H(z) nodes with provenance columns.",
    )
    args = p.parse_args()

    base_csv = Path(args.base_csv).resolve()
    cand_csv = Path(args.candidate_csv).resolve()

    _write_template_if_missing(TEMPLATE_CSV)

    base_exists = base_csv.exists()
    cand_exists = cand_csv.exists()

    base_rows: List[Dict] = []
    cand_rows: List[Dict] = []
    base_err = ""
    cand_err = ""

    if base_exists:
        base_rows, base_err = _read_rows(base_csv, require_provenance=False)
    if cand_exists:
        cand_rows, cand_err = _read_rows(cand_csv, require_provenance=True)

    provenance_per_row_ok = bool(
        cand_rows
        and all(
            bool(r["citation"])
            and bool(r["source_version"])
            and _is_https_url(r["reference_url"])
            for r in cand_rows
        )
    )
    candidate_rows_ge_2 = len(cand_rows) >= 2

    merged_rows: List[Dict] = []
    added_unique = 0
    if base_rows and cand_rows:
        merged_rows, added_unique = _merge_unique(base_rows, cand_rows)
    elif base_rows:
        merged_rows = list(base_rows)

    m = _metrics(merged_rows)
    n_nodes_ge_5 = bool(m["n_nodes"] >= 5)
    z_span_ge_0p8 = bool(math.isfinite(m["z_span"]) and m["z_span"] >= 0.8)
    e_span_ge_1p0 = bool(math.isfinite(m["e_span"]) and m["e_span"] >= 1.0)
    cond_lt_8 = bool(math.isfinite(m["design_condition_number"]) and m["design_condition_number"] < 8.0)
    has_high_z_node = bool(any(float(r["z"]) >= 1.18 for r in merged_rows))

    flags = {
        "base_csv_exists": base_exists,
        "base_csv_parse_ok": bool(base_exists and not base_err),
        "candidate_csv_exists": cand_exists,
        "candidate_csv_parse_ok": bool(cand_exists and not cand_err),
        "candidate_rows_ge_2": candidate_rows_ge_2,
        "candidate_provenance_per_row_complete": provenance_per_row_ok,
        "merged_added_unique_ge_2": bool(added_unique >= 2),
        "merged_n_nodes_ge_5": n_nodes_ge_5,
        "merged_z_span_ge_0p8": z_span_ge_0p8,
        "merged_e_span_ge_1p0": e_span_ge_1p0,
        "merged_design_condition_lt_8": cond_lt_8,
        "merged_has_node_z_ge_1p18": has_high_z_node,
    }

    all_ready = all(bool(v) for v in flags.values())
    verdict = "HZ_STRICT_NODE_PACK_READY" if all_ready else "HZ_STRICT_NODE_PACK_PENDING"

    if merged_rows:
        with MERGED_OUT_CSV.open("w", encoding="utf-8", newline="") as f:
            w = csv.DictWriter(f, fieldnames=["z", "h_km_s_mpc", "sigma_total"])
            w.writeheader()
            for r in merged_rows:
                w.writerow(
                    {
                        "z": f"{float(r['z']):.8g}",
                        "h_km_s_mpc": f"{float(r['h_km_s_mpc']):.8g}",
                        "sigma_total": f"{float(r['sigma_total']):.8g}",
                    }
                )

        out_meta = {
            "source_dataset_name": "Merged strict candidate H(z) snapshot",
            "citation": "Multiple external sources; see candidate rows in report_qw2112_hz_strict_node_pack_gate.json",
            "reference_url": "https://arxiv.org/abs/1607.03155",
            "source_version": "QW2112_MERGED_CANDIDATE_V1",
            "provenance_anchor_free": True,
            "seeded_from_registry": False,
            "notes": "Generated by QW-2112 from baseline + externally provided candidate nodes.",
            "record_count": int(len(merged_rows)),
            "columns_schema": ["z", "h_km_s_mpc", "sigma_total"],
        }
        MERGED_OUT_META.write_text(json.dumps(out_meta, ensure_ascii=False, indent=2), encoding="utf-8")

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "inputs": {
            "base_csv": str(base_csv),
            "candidate_csv": str(cand_csv),
            "candidate_template_csv": str(TEMPLATE_CSV),
        },
        "errors": {
            "base_error": base_err,
            "candidate_error": cand_err,
        },
        "counts": {
            "n_base_rows": len(base_rows),
            "n_candidate_rows": len(cand_rows),
            "n_added_unique": int(added_unique),
            "n_merged_rows": len(merged_rows),
        },
        "merged_metrics": m,
        "candidate_rows_preview": cand_rows[:10],
        "flags": flags,
        "pass_count": sum(1 for v in flags.values() if bool(v)),
        "total_flags": len(flags),
        "output_artifacts": {
            "merged_csv": str(MERGED_OUT_CSV) if merged_rows else None,
            "merged_metadata": str(MERGED_OUT_META) if merged_rows else None,
        },
        "verdict": verdict,
        "required_next_step": (
            "USE_MERGED_HZ_FILE_IN_QW2099_AND_RERUN_QW2102_QW2090_QW2104_QW2105_QW2094"
            if all_ready
            else "COLLECT_REAL_EXTERNAL_HZ_NODES_WITH_PROVENANCE_AND_RERUN_QW2112"
        ),
    }
    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-2112: HZ STRICT NODE PACK GATE",
        "",
        f"- Date UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- pass_count: `{out['pass_count']}/{out['total_flags']}`",
        "",
        "## Counts",
        f"- n_base_rows: `{len(base_rows)}`",
        f"- n_candidate_rows: `{len(cand_rows)}`",
        f"- n_added_unique: `{added_unique}`",
        f"- n_merged_rows: `{len(merged_rows)}`",
        "",
        "## Merged Metrics",
        f"- n_nodes: `{m['n_nodes']}`",
        f"- z_span: `{m['z_span']}`",
        f"- e_span: `{m['e_span']}`",
        f"- design_condition_number: `{m['design_condition_number']}`",
        "",
        "## Artifact",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-2112] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-2112] Saved MD:   {OUT_MD.name}")
    print(f"[QW-2112] verdict={verdict} pass_count={out['pass_count']}/{out['total_flags']}")


if __name__ == "__main__":
    main()
