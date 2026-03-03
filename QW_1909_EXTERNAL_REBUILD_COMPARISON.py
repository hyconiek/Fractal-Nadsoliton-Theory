#!/usr/bin/env python3
"""
QW-1909: External-source rebuild comparison gate.

Compares multiple external-source candidate datasets under fixed QW-1853 criteria.
No threshold changes, no post-hoc model changes.
"""

from __future__ import annotations

import importlib.util
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List

import pandas as pd


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1909_external_rebuild_comparison.json"
OUT_MD = ROOT / "RAPORT_QW1909_EXTERNAL_REBUILD_COMPARISON.md"


def load_json(name: str) -> Dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def externality_ok(provider: str, ext_stmt: str) -> bool:
    text = ext_stmt.lower()
    return (
        provider not in {"INTERNAL_PROXY", "INTERNAL"}
        and ("not independent" not in text)
        and ("internal proxy" not in text)
        and ("independent" in text)
    )


def candidate_eval(candidate_dir: Path, pta_thr: Dict[str, float], gw_thr: Dict[str, float], m1853) -> Dict:
    out = {"candidate_dir": str(candidate_dir), "exists": candidate_dir.exists()}
    if not candidate_dir.exists():
        out["error"] = "candidate_missing"
        return out

    manifest_path = candidate_dir / "manifest.json"
    if not manifest_path.exists():
        out["error"] = "manifest_missing"
        return out

    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    role_map = {x.get("role"): x for x in manifest.get("files", []) if isinstance(x, dict)}
    if "pta_pairs" not in role_map or "gw_windows" not in role_map:
        out["error"] = "missing_roles"
        return out

    pta_path = candidate_dir / str(role_map["pta_pairs"]["path"])
    gw_path = candidate_dir / str(role_map["gw_windows"]["path"])
    if not pta_path.exists() or not gw_path.exists():
        out["error"] = "data_file_missing"
        return out

    df_pta = pd.read_csv(pta_path)
    df_gw = pd.read_csv(gw_path)

    pta_eval = m1853.eval_pta_v2(df_pta, pta_thr)
    gw_eval = m1853.eval_gw(df_gw, gw_thr)

    pta_all = all(bool(v) for v in pta_eval["pass_flags"].values())
    gw_all = all(bool(v) for v in gw_eval["pass_flags"].values())
    hard_gate = "PASS" if (pta_all and gw_all) else ("PARTIAL" if (pta_all or gw_all) else "FAIL")

    provider = str(manifest.get("dataset", {}).get("provider", ""))
    ext_stmt = str(manifest.get("dataset", {}).get("externality_statement", ""))
    ext_ok = externality_ok(provider, ext_stmt)

    out.update(
        {
            "dataset_id": manifest.get("dataset", {}).get("dataset_id"),
            "provider": provider,
            "externality_ok": ext_ok,
            "pta_summary": pta_eval["summary"],
            "pta_pass_flags": pta_eval["pass_flags"],
            "gw_summary": gw_eval["summary"],
            "gw_pass_flags": gw_eval["pass_flags"],
            "pta_all": pta_all,
            "gw_all": gw_all,
            "hard_gate": hard_gate,
            "row_counts": {"pta": int(len(df_pta)), "gw": int(len(df_gw))},
        }
    )
    return out


def main() -> None:
    d1839 = load_json("report_qw1839_joint_confirmatory_prereg_protocol.json")
    d1850 = load_json("report_qw1850_pta_v2_prereg_protocol.json")

    m1853 = load_module(ROOT / "QW_1853_JOINT_EXTERNAL_CONFIRMATORY_V2.py", "qw1853_cmp")
    pta_thr = d1850["protocol"]["pta_v2_protocol"]["thresholds"]
    gw_thr = d1839["protocol"]["gw_protocol"]["thresholds"]

    candidates = [
        ROOT / "external_confirmatory_v2" / "confirmatory_dataset_external_source_rebuild",
        ROOT / "external_confirmatory_v2" / "confirmatory_dataset_external_source_rebuild_v2_1831cfg",
    ]

    rows: List[Dict] = []
    for c in candidates:
        rows.append(candidate_eval(c, pta_thr=pta_thr, gw_thr=gw_thr, m1853=m1853))

    valid = [r for r in rows if "error" not in r]

    # Ranking: hard gate PASS>PARTIAL>FAIL, then GW AUC, then PTA prob
    rank_map = {"PASS": 2, "PARTIAL": 1, "FAIL": 0}
    valid_sorted = sorted(
        valid,
        key=lambda r: (
            rank_map.get(r["hard_gate"], -1),
            float(r["gw_summary"]["calibrated_mean_auc"]),
            float(r["pta_summary"]["prob_pair_mean_gain_positive"]),
        ),
        reverse=True,
    )

    best = valid_sorted[0] if valid_sorted else None
    best_hard = best["hard_gate"] if best else "FAIL"

    if best and best_hard == "PASS" and best["externality_ok"]:
        verdict = "BEST_CANDIDATE_PASSES_JOINT_AND_EXTERNALITY"
        readiness = "EMPIRICAL_CLOSURE_PASS_READY"
    elif best and best_hard in {"PASS", "PARTIAL"}:
        verdict = "BEST_CANDIDATE_PARTIAL_OR_EXTERNALITY_BLOCKED"
        readiness = "EMPIRICAL_CLOSURE_NOT_READY"
    else:
        verdict = "ALL_CANDIDATES_FAIL_JOINT_GATE"
        readiness = "EMPIRICAL_CLOSURE_FAIL"

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "threshold_sources": {
            "pta_protocol_sha256": d1850.get("protocol_sha256"),
            "gw_protocol_sha256": d1839.get("protocol_sha256"),
        },
        "candidates": rows,
        "best_candidate": best,
        "verdict": verdict,
        "readiness": readiness,
    }
    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1909: EXTERNAL REBUILD COMPARISON",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- readiness: **{readiness}**",
        "",
        "## Candidates",
    ]
    for r in rows:
        if "error" in r:
            lines.append(f"- `{r['candidate_dir']}`: error={r['error']}")
            continue
        lines.append(
            f"- `{r['candidate_dir']}`: hard={r['hard_gate']}, pta_all={r['pta_all']}, gw_all={r['gw_all']}, "
            f"externality_ok={r['externality_ok']}, pta_prob={r['pta_summary']['prob_pair_mean_gain_positive']:.3f}, "
            f"gw_auc={r['gw_summary']['calibrated_mean_auc']:.3f}"
        )

    if best:
        lines.extend(
            [
                "",
                "## Best Candidate",
                f"- dir: `{best['candidate_dir']}`",
                f"- hard_gate: {best['hard_gate']}",
                f"- externality_ok: {best['externality_ok']}",
                f"- PTA mean gain: {best['pta_summary']['mean_pair_mean_gain']:.6f}",
                f"- PTA prob positive: {best['pta_summary']['prob_pair_mean_gain_positive']:.3f}",
                f"- GW mean AUC: {best['gw_summary']['calibrated_mean_auc']:.3f}",
                f"- GW mean adv: {best['gw_summary']['calibrated_mean_adv']:.3f}",
            ]
        )

    lines.extend(["", "## Artifacts", f"- JSON: `{OUT_JSON.name}`"])
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1909] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1909] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()

