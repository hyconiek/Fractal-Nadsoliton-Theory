#!/usr/bin/env python3
"""
QW-1902: Empirical closure gate with externality check.

Separates metric performance from dataset independence.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict

ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1902_empirical_closure_gate.json"
OUT_MD = ROOT / "RAPORT_QW1902_EMPIRICAL_CLOSURE_GATE.md"


def read_json(name: str) -> Dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def to_score(v: float, lo: float, hi: float) -> float:
    if hi <= lo:
        return 0.0
    x = (float(v) - lo) / (hi - lo)
    return max(0.0, min(1.0, x))


def main() -> None:
    d1898 = read_json("report_qw1898_empirical_bridge_precondition_gate.json")
    d1899 = read_json("report_qw1899_external_detector_protocol_design.json")
    d1852 = read_json("report_qw1852_external_confirmatory_data_precheck.json")
    d1853 = read_json("report_qw1853_joint_external_confirmatory_v2.json")

    candidate_dir = Path(d1852.get("validation", {}).get("candidate_dir", ""))
    manifest_path = candidate_dir / "manifest.json"

    manifest = {}
    if manifest_path.exists():
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))

    provider = str(manifest.get("dataset", {}).get("provider", ""))
    ext_stmt = str(manifest.get("dataset", {}).get("externality_statement", ""))

    text = ext_stmt.lower()
    externality_ok = (
        provider not in {"INTERNAL_PROXY", "INTERNAL"}
        and ("not independent" not in text)
        and ("internal proxy" not in text)
        and ("independent" in text)
    )

    joint_flags = d1853.get("joint_pass_flags", {})
    pta_all = bool(joint_flags.get("pta_v2_all", False))
    gw_all = bool(joint_flags.get("gw_all", False))

    hard_joint = d1853.get("hard_gate", "PARTIAL")
    hard_ok = hard_joint == "PASS"

    # Score from achieved metrics.
    pta_sum = d1853.get("pta_v2", {}).get("summary", {})
    gw_sum = d1853.get("gw", {}).get("summary", {})

    metric_score = (
        0.35 * to_score(gw_sum.get("calibrated_mean_auc", 0.0), 0.60, 0.95)
        + 0.20 * to_score(gw_sum.get("calibrated_mean_adv", 0.0), 0.05, 0.80)
        + 0.25 * to_score(pta_sum.get("prob_pair_mean_gain_positive", 0.0), 0.50, 0.75)
        + 0.20 * to_score(pta_sum.get("one_sided_lower95_prob_pair_mean_gain_positive", 0.0), 0.45, 0.65)
    )

    precondition_ok = d1898.get("readiness") == "EMPIRICAL_BRIDGE_PRECONDITION_PASS"

    if precondition_ok and hard_ok and externality_ok:
        readiness = "EMPIRICAL_CLOSURE_PASS"
    elif precondition_ok and not externality_ok:
        readiness = "INTERNAL_REHEARSAL_ONLY_NOT_CONFIRMATORY"
    elif precondition_ok:
        readiness = "EMPIRICAL_CLOSURE_PARTIAL"
    else:
        readiness = "EMPIRICAL_CLOSURE_FAIL"

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "inputs": {
            "precondition_readiness": d1898.get("readiness"),
            "protocol_sha256": d1899.get("protocol_sha256"),
            "candidate_dir": str(candidate_dir),
            "provider": provider,
            "externality_statement": ext_stmt,
            "externality_ok": externality_ok,
            "joint_flags": {
                "pta_all": pta_all,
                "gw_all": gw_all,
                "hard_gate": hard_joint,
            },
        },
        "metric_score": metric_score,
        "readiness": readiness,
        "verdict": "EMPIRICAL_CLOSURE_GATE_COMPLETE",
        "required_next_step": (
            "DELIVER_TRULY_EXTERNAL_FROZEN_DATASET_AND_RERUN_1852_1853"
            if not externality_ok
            else "PREPARE_FINAL_EMPIRICAL_PUBLICATION_PACKAGE"
        ),
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1902: EMPIRICAL CLOSURE GATE",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{out['verdict']}**",
        f"- readiness: **{readiness}**",
        f"- metric_score: {metric_score:.3f}",
        "",
        "## Externality",
        f"- provider: `{provider}`",
        f"- externality_ok: {externality_ok}",
        "",
        "## Joint Execution Flags",
        f"- hard_gate: {hard_joint}",
        f"- pta_all: {pta_all}",
        f"- gw_all: {gw_all}",
        "",
        "## Required Next Step",
        f"- {out['required_next_step']}",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]

    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1902] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1902] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
