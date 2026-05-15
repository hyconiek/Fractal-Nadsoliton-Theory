#!/usr/bin/env python3
from __future__ import annotations
import json
from datetime import date
from pathlib import Path

ROOT = Path(__file__).resolve().parent
OUT_POLICY = ROOT / "generated" / "p1800_s750_strict_bw_componentwise_divergence_trace_evidence_pack_checkpoint.json"
OUT_TEMPLATE = ROOT / "generated" / "p1800_s750_bw_componentwise_divergence_trace_evidence_pack_template.json"

policy = {
    "packet_id": "P1800_S750",
    "as_of": str(date(2026, 5, 15)),
    "required_components": ["B1", "B2", "B3", "C1", "C2"],
    "required_divergence_trace": True,
    "allowed_verdicts": ["PASS_ZERO", "OPEN_OBSTRUCTION_WITH_TRACE"],
    "pass_zero_requirements": [
        "all_components_final_zero",
        "divergence_trace_final_zero",
        "no_freeze_or_index_mismatch",
        "all_check_digests_present"
    ]
}

template = {
    "freeze_id_common": "FREEZE_ID_REQUIRED",
    "background_family_id": "BG_FAMILY_REQUIRED",
    "index_convention_id": "INDEX_CONVENTION_REQUIRED",
    "classification": {
        "scope": "LOCAL",
        "representation": "NONPROXY",
        "artifact_level": "FULL_EXPORT"
    },
    "component_residuals": {
        k: {"symbolic": "EXPR", "final_simplified": "EXPR", "check_digest": "HASH"}
        for k in ["B1", "B2", "B3", "C1", "C2"]
    },
    "divergence_trace": {
        "symbolic": "EXPR",
        "final_simplified": "EXPR",
        "check_digest": "HASH"
    },
    "verdict_candidate": "OPEN_OBSTRUCTION_WITH_TRACE"
}

for p, d in ((OUT_POLICY, policy), (OUT_TEMPLATE, template)):
    p.parent.mkdir(parents=True, exist_ok=True)
    with p.open("w", encoding="utf-8") as f:
        json.dump(d, f, indent=2, ensure_ascii=False)
        f.write("\n")
    print(p)
