#!/usr/bin/env python3
from __future__ import annotations
import json
from datetime import date
from pathlib import Path

ROOT = Path(__file__).resolve().parent
OUT_POLICY = ROOT / "generated" / "p1799_s749_strict_sv1_to_sv5_state_update_immutability_audit_checkpoint.json"
OUT_TEMPLATE = ROOT / "generated" / "p1799_s749_state_update_immutability_audit_input_template.json"

policy = {
    "packet_id": "P1799_S749",
    "as_of": str(date(2026, 5, 15)),
    "allowed_verdicts": [
        "AUDIT_PASS_IMMUTABILITY_CONFIRMED",
        "AUDIT_FAIL_OPEN_OBSTRUCTION_WITH_TRACE",
    ],
    "immutability_checks": [
        "SV6_unchanged",
        "SV7_unchanged",
        "SV8_unchanged",
    ],
    "gate_lock_checks": [
        "if_G_BW_not_PASS_ZERO_then_G_BRST_LOCKED",
        "if_G_BW_not_PASS_ZERO_then_G_CUT_LOCKED",
    ],
    "fail_effect": "BLOCK_GATE_PROMOTION_FOR_RUN",
}

template = {
    "state_before": {
        "SV6": "OPEN_LOCKED",
        "SV7": "OPEN_LOCKED",
        "SV8": "OPEN_LOCKED"
    },
    "state_after": {
        "SV6": "OPEN_LOCKED",
        "SV7": "OPEN_LOCKED",
        "SV8": "OPEN_LOCKED"
    },
    "gates": {
        "G_BW": "OPEN_OBSTRUCTION_WITH_TRACE",
        "G_BRST": "LOCKED",
        "G_CUT": "LOCKED"
    }
}

for p, d in ((OUT_POLICY, policy), (OUT_TEMPLATE, template)):
    p.parent.mkdir(parents=True, exist_ok=True)
    with p.open("w", encoding="utf-8") as f:
        json.dump(d, f, indent=2, ensure_ascii=False)
        f.write("\n")
    print(p)
