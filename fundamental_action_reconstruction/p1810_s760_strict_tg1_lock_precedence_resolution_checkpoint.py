#!/usr/bin/env python3
from __future__ import annotations
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
G = ROOT / "generated"
OUT = G / "p1810_s760_strict_tg1_lock_precedence_resolution_checkpoint.json"

sources = {
    "p1805": json.loads((G/"p1805_s755_strict_bw_entry_nonproxy_covariant_export_completeness_audit_checkpoint.json").read_text()).get("gate_vector",{}).get("TG1_BW","OPEN_LOCKED_BY_UNIFIED_NONPROXY_RESIDUAL"),
    "p1808": json.loads((G/"p1808_s758_strict_w1_full_export_gate_reconciliation_checkpoint.json").read_text()).get("gate_vector",{}).get("TG1_BW","OPEN_LOCKED_BY_UNIFIED_NONPROXY_RESIDUAL"),
    "p1809": json.loads((G/"p1809_s759_strict_nonproxy_export_semantic_level_reconciliation_checkpoint.json").read_text()).get("derived_gate_vector",{}).get("TG1_BW","OPEN_LOCKED_BY_UNIFIED_NONPROXY_RESIDUAL")
}

rank = {
    "OPEN_LOCKED_BY_W1_CONSISTENCY_CONFLICT": 3,
    "OPEN_LOCKED_BY_UNIFIED_NONPROXY_RESIDUAL": 2,
    "PASS_ZERO": 1,
}

resolved = sorted(sources.values(), key=lambda s: rank.get(s,0), reverse=True)[0]
out = {
    "packet_id":"P1810",
    "stage_id":"S760",
    "status":"PASS_ZERO_PRECEDENCE_APPLIED",
    "tg1_sources":sources,
    "tg1_resolved":resolved,
    "rule":"MAX_BLOCKING_PRECEDENCE"
}
OUT.write_text(json.dumps(out,indent=2)+"\n")
print(str(OUT))
