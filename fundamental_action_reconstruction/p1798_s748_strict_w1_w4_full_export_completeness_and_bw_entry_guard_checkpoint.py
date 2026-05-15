#!/usr/bin/env python3
from __future__ import annotations
import json
from datetime import date
from pathlib import Path

ROOT = Path(__file__).resolve().parent
OUT_POLICY = ROOT / "generated" / "p1798_s748_strict_w1_w4_full_export_completeness_and_bw_entry_guard_checkpoint.json"
OUT_MATRIX = ROOT / "generated" / "p1798_s748_w1_w4_full_export_matrix_template.json"

policy = {
    "packet_id": "P1798_S748",
    "as_of": str(date(2026, 5, 15)),
    "bw_entry_requirements": [
        "W1_FULL_EXPORT_VERIFIED",
        "W2_FULL_EXPORT_VERIFIED",
        "W3_FULL_EXPORT_VERIFIED",
        "W4_FULL_EXPORT_VERIFIED",
        "common_freeze_id",
        "no_open_obstruction_in_w_matrix",
    ],
    "bw_entry_fail_effect": {
        "G_BW_ENTRY_ALLOWED": False,
        "G_BW": "OPEN_OBSTRUCTION_WITH_TRACE",
        "G_BRST": "LOCKED",
        "G_CUT": "LOCKED",
    },
}

row = {
    "artifact_id": "REQUIRED",
    "classification": {
        "representation": "NONPROXY",
        "artifact_level": "FULL_EXPORT"
    },
    "freeze_id": "FREEZE_ID_COMMON",
    "verification": {
        "check_command": "python3 verify_export.py --artifact REQUIRED",
        "check_digest": "HASH_REQUIRED"
    },
    "status": "OPEN_OBSTRUCTION_WITH_TRACE"
}

matrix = {
    "freeze_id_common": "FREEZE_ID_COMMON",
    "W1": row,
    "W2": row,
    "W3": row,
    "W4": row
}

for p, d in ((OUT_POLICY, policy), (OUT_MATRIX, matrix)):
    p.parent.mkdir(parents=True, exist_ok=True)
    with p.open("w", encoding="utf-8") as f:
        json.dump(d, f, indent=2, ensure_ascii=False)
        f.write("\n")
    print(p)
