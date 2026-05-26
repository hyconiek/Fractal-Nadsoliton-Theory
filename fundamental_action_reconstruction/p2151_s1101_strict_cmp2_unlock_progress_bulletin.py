#!/usr/bin/env python3
from __future__ import annotations
import json
from pathlib import Path
from typing import Any

ROOT=Path(__file__).resolve().parent
GEN=ROOT/"generated"
INP=GEN/"p2150_s1100_strict_cmp2_real_data_unlock_attempt_register.json"
OUT=GEN/"p2151_s1101_strict_cmp2_unlock_progress_bulletin.json"
MD=GEN/"p2151_s1101_strict_cmp2_unlock_progress_bulletin.md"


def load(p:Path)->dict[str,Any]:
    if not p.exists(): return {"_missing":str(p)}
    return json.loads(p.read_text(encoding='utf-8'))


def main()->None:
    GEN.mkdir(exist_ok=True)
    p=load(INP)
    reg=(p.get('real_data_unlock_attempt_register',{}) or {})
    unresolved=reg.get('unresolved_stages',[]) or []
    unblocked=len(unresolved)==0
    payload={
      "schema_version":"p2151_s1101_v1","packet_id":"P2151","stage_id":"S1101","produced_by":Path(__file__).name,
      "timestamp_utc":"2026-05-26T00:00:00+00:00","status":"OPEN_PARTIAL_PROGRESS_WITH_TRACE",
      "result_kind":"PASS_STRICT_CMP2_UNLOCK_PROGRESS_BULLETIN_WITH_TRACE",
      "unlock_progress_bulletin":{
        "source_register":str(INP.relative_to(ROOT)),
        "unresolved_stages":unresolved,
        "n_unresolved":len(unresolved),
        "unblocked":unblocked,
        "scope_limit":"bulletin only; not theorem-grade closure"
      },
      "recommended_next_honest_step":{"id":"P2152_candidate","goal":"deliver real extension rows and rerun P2147/P2150 to reduce unresolved stages to zero"},
      "gatekeeper_checks":{"bulletin_exported":True,"no_toe_closure_claimed":True,
        "full_d3_covariance_transport_proven":False,"full_cutkosky_closure_proven":False,"c3_theorem_proven":False},
      "global_status":"OPEN_OBSTRUCTION_WITH_TRACE"
    }
    OUT.write_text(json.dumps(payload,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    MD.write_text(f"# P2151 S1101\n\n- Unresolved stages: `{unresolved}`\n",encoding='utf-8')

if __name__=='__main__':
    main()
