#!/usr/bin/env python3
from __future__ import annotations
import json
from pathlib import Path
from typing import Any

ROOT=Path(__file__).resolve().parent
GEN=ROOT/"generated"
INP=GEN/"p2147_s1097_strict_cmp2_real_data_required_rerun_checkpoint.json"
OUT=GEN/"p2148_s1098_strict_cmp2_external_data_delivery_blocker_report.json"
MD=GEN/"p2148_s1098_strict_cmp2_external_data_delivery_blocker_report.md"


def load(p:Path)->dict[str,Any]:
    if not p.exists(): return {"_missing":str(p)}
    return json.loads(p.read_text(encoding='utf-8'))


def main()->None:
    GEN.mkdir(exist_ok=True)
    p=load(INP)
    chk=(p.get('real_data_required_rerun_checkpoint',{}) or {})
    r=(chk.get('result_kinds',{}) or {})
    blocked=[k for k in ['p2132','p2133','p2134','p2135'] if not str(r.get(k,'OPEN')).startswith('PASS_')]
    ext_present=bool(chk.get('extension_present',False))

    payload={
      "schema_version":"p2148_s1098_v1","packet_id":"P2148","stage_id":"S1098","produced_by":Path(__file__).name,
      "timestamp_utc":"2026-05-26T00:00:00+00:00","status":"OPEN_PARTIAL_PROGRESS_WITH_TRACE",
      "result_kind":"OPEN_STRICT_CMP2_EXTERNAL_DATA_DELIVERY_BLOCKER_REPORT" if blocked else "PASS_STRICT_CMP2_EXTERNAL_DATA_DELIVERY_BLOCKER_REPORT",
      "external_data_delivery_blocker_report":{
        "source_checkpoint":str(INP.relative_to(ROOT)),
        "extension_present":ext_present,
        "blocked_stages":blocked,
        "required_external_artifact":"generated/cmp2_backend_rows_extension_v1.json",
        "scope_limit":"blocker report only; requires external data delivery",
      },
      "recommended_next_honest_step":{
        "id":"P2149_candidate",
        "goal":"deliver real extension rows then rerun P2147 and confirm p2132-p2135 PASS before any CI stability interpretation"
      },
      "gatekeeper_checks":{
        "external_data_required":True,
        "no_toe_closure_claimed":True,
        "full_d3_covariance_transport_proven":False,
        "full_cutkosky_closure_proven":False,
        "c3_theorem_proven":False
      },
      "global_status":"OPEN_OBSTRUCTION_WITH_TRACE"
    }
    OUT.write_text(json.dumps(payload,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    MD.write_text(f"# P2148 S1098\n\n- Result kind: `{payload['result_kind']}`\n- Blocked stages: `{blocked}`\n",encoding='utf-8')

if __name__=='__main__':
    main()
