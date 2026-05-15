#!/usr/bin/env python3
from __future__ import annotations
import json
from datetime import date
from pathlib import Path

ROOT = Path(__file__).resolve().parent
OUT_POLICY = ROOT / "generated" / "p1804_s754_strict_one_shot_sv1_to_tg3_execution_bundle_verdict_checkpoint.json"
OUT_TEMPLATE = ROOT / "generated" / "p1804_s754_one_shot_sv1_to_tg3_bundle_input_template.json"

policy = {
    "packet_id": "P1804_S754",
    "as_of": str(date(2026, 5, 15)),
    "bundle_pass_requirements": [
        "SV1_to_SV5_all_pass",
        "TG1_PASS",
        "TG2_PASS",
        "TG3_PASS",
        "Helmholtz_PASS",
        "Renorm_PASS",
        "Background_consistency_PASS",
        "shared_freeze_id",
        "all_check_digests_present"
    ],
    "bundle_fail_verdict": "BUNDLE_OPEN_OBSTRUCTION_WITH_TRACE"
}

template = {
    "freeze_id_common": "FREEZE_ID_REQUIRED",
    "sv_section": {"SV1": "OPEN", "SV2": "OPEN", "SV3": "OPEN", "SV4": "OPEN", "SV5": "OPEN"},
    "tg_section": {"TG1": "OPEN", "TG2": "OPEN", "TG3": "OPEN"},
    "global_closures": {
        "Helmholtz": "OPEN",
        "renorm": "OPEN",
        "background_consistency": "OPEN"
    },
    "missing_witnesses": [],
    "first_fail_trace": ""
}

for p, d in ((OUT_POLICY, policy), (OUT_TEMPLATE, template)):
    p.parent.mkdir(parents=True, exist_ok=True)
    with p.open("w", encoding="utf-8") as f:
        json.dump(d, f, indent=2, ensure_ascii=False)
        f.write("\n")
    print(p)
