#!/usr/bin/env python3
"""
QW-2002: Single-kernel triple-sector closure gate v3.

Bridges stage-B closure status with QW-2001 lockable-triad result.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict


ROOT = Path(__file__).resolve().parent
IN_QW1934 = ROOT / "report_qw1934_strict_closure_gate_solo.json"
IN_QW2001 = ROOT / "report_qw2001_bounded_gw_triad_lockable_gate.json"

OUT_JSON = ROOT / "report_qw2002_single_kernel_triple_sector_closure_gate_v3.json"
OUT_MD = ROOT / "RAPORT_QW2002_SINGLE_KERNEL_TRIPLE_SECTOR_CLOSURE_GATE_V3.md"


def load(path: Path) -> Dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    d1934 = load(IN_QW1934)
    d2001 = load(IN_QW2001)

    stage_b_solo_closed = bool(d1934.get("readiness") == "TOE_STAGE_B_SOLO_CLOSED")
    lockable_pass = bool(d2001.get("verdict") == "BOUNDED_GW_TRIAD_LOCKABLE_PASS")

    if stage_b_solo_closed and lockable_pass:
        readiness = "TOE_STAGE_C_SINGLE_KERNEL_CLOSED_LOCKABLE_INTERNAL"
        verdict = "SINGLE_KERNEL_TRIPLE_SECTOR_CLOSURE_PASS_V3"
        required_next_step = "FREEZE_FULL_PARAMETER_VECTOR_AND_RUN_TRUE_EXTERNAL_CONFIRMATORY_PACKAGE"
    elif stage_b_solo_closed:
        readiness = "TOE_STAGE_C_BLOCKED_LOCKABILITY_NOT_REACHED"
        verdict = "SINGLE_KERNEL_TRIPLE_SECTOR_CLOSURE_FAIL_V3"
        required_next_step = "IMPROVE_BOOTSTRAP_LOCKABILITY_WITHOUT_SECTOR_RETUNE"
    else:
        readiness = "TOE_STAGE_B_NOT_CLOSED_CANNOT_CLOSE_STAGE_C"
        verdict = "SINGLE_KERNEL_TRIPLE_SECTOR_CLOSURE_BLOCKED_BY_STAGE_B_V3"
        required_next_step = "RECOVER_STAGE_B_SOLO_CLOSURE"

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "inputs": {
            "q1934_readiness": d1934.get("readiness"),
            "q1934_verdict": d1934.get("verdict"),
            "q2001_verdict": d2001.get("verdict"),
        },
        "flags": {
            "stage_b_solo_closed": stage_b_solo_closed,
            "q2001_lockable_pass": lockable_pass,
        },
        "evidence_snapshot": {
            "q2001_deterministic_all_pass": d2001.get("deterministic", {}).get("all_pass"),
            "q2001_boot5000_min": d2001.get("bootstrap", {}).get("triad_pass_rate_5000_min"),
            "q2001_local_r01": d2001.get("local_deterministic_pass_rates", {}).get("r01"),
            "q2001_local_r02": d2001.get("local_deterministic_pass_rates", {}).get("r02"),
            "q2001_local_r05": d2001.get("local_deterministic_pass_rates", {}).get("r05"),
        },
        "readiness": readiness,
        "verdict": verdict,
        "required_next_step": required_next_step,
    }
    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-2002: SINGLE KERNEL TRIPLE-SECTOR CLOSURE GATE V3",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- Readiness: **{readiness}**",
        "",
        "## Inputs",
        f"- QW-1934 readiness/verdict: {out['inputs']['q1934_readiness']} / {out['inputs']['q1934_verdict']}",
        f"- QW-2001 verdict: {out['inputs']['q2001_verdict']}",
        "",
        "## Evidence Snapshot",
        f"- deterministic_all_pass: {out['evidence_snapshot']['q2001_deterministic_all_pass']}",
        f"- boot5000_min: {out['evidence_snapshot']['q2001_boot5000_min']}",
        f"- local pass rates r01/r02/r05: {out['evidence_snapshot']['q2001_local_r01']} / {out['evidence_snapshot']['q2001_local_r02']} / {out['evidence_snapshot']['q2001_local_r05']}",
        "",
        "## Required Next Step",
        f"- {required_next_step}",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-2002] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-2002] Saved MD:   {OUT_MD.name}")
    print(f"[QW-2002] readiness={readiness} verdict={verdict}")


if __name__ == "__main__":
    main()
