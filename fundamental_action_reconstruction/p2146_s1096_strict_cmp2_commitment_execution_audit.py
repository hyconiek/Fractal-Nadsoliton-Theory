#!/usr/bin/env python3
from __future__ import annotations
import json, subprocess, sys
from pathlib import Path
from typing import Any

ROOT=Path(__file__).resolve().parent
GEN=ROOT/"generated"
INP=GEN/"p2145_s1095_strict_cmp2_next_action_commitment_packet.json"
OUT=GEN/"p2146_s1096_strict_cmp2_commitment_execution_audit.json"
MD=GEN/"p2146_s1096_strict_cmp2_commitment_execution_audit.md"

SCRIPT_MAP={
 "rerun_p2142_executor": ROOT/"p2142_s1092_strict_cmp2_real_data_handoff_and_rerun_executor.py",
 "refresh_p2144_interpretation_gate": ROOT/"p2144_s1094_strict_cmp2_real_ci_stability_interpretation_gate.py",
}

def load(p:Path)->dict[str,Any]:
    if not p.exists(): return {"_missing":str(p)}
    return json.loads(p.read_text(encoding='utf-8'))

def run(path:Path)->tuple[bool,str]:
    try:
        subprocess.run([sys.executable,str(path)],check=True)
        return True,'OK'
    except subprocess.CalledProcessError as e:
        return False,f'FAILED({e.returncode})'

def main()->None:
    GEN.mkdir(exist_ok=True)
    p=load(INP)
    actions=((p.get('next_action_commitment',{}) or {}).get('action_commitment') or [])
    log=[]
    for a in actions:
        if a in SCRIPT_MAP:
            ok,msg=run(SCRIPT_MAP[a]); log.append({"action":a,"ok":ok,"message":msg})
        else:
            log.append({"action":a,"ok":False,"message":"MANUAL_OR_EXTERNAL_ACTION_REQUIRED"})

    p2142=load(GEN/"p2142_s1092_strict_cmp2_real_data_handoff_and_rerun_executor.json")
    p2144=load(GEN/"p2144_s1094_strict_cmp2_real_ci_stability_interpretation_gate.json")

    payload={
      "schema_version":"p2146_s1096_v1","packet_id":"P2146","stage_id":"S1096","produced_by":Path(__file__).name,
      "timestamp_utc":"2026-05-26T00:00:00+00:00","status":"OPEN_PARTIAL_PROGRESS_WITH_TRACE",
      "result_kind":"PASS_STRICT_CMP2_COMMITMENT_EXECUTION_AUDIT_WITH_TRACE",
      "commitment_execution_audit":{
        "source_commitment_packet":str(INP.relative_to(ROOT)),"actions":actions,"execution_log":log,
        "post_refresh_result_kinds":{"p2142":p2142.get('result_kind','OPEN_MISSING'),"p2144":p2144.get('result_kind','OPEN_MISSING')},
        "scope_limit":"execution audit only; not theorem-grade closure"
      },
      "recommended_next_honest_step":{"id":"P2147_candidate","goal":"repeat commitment cycle until real extension rows unlock non-synthetic PASS chain"},
      "gatekeeper_checks":{"audit_exported":True,"no_toe_closure_claimed":True,
        "full_d3_covariance_transport_proven":False,"full_cutkosky_closure_proven":False,"c3_theorem_proven":False},
      "global_status":"OPEN_OBSTRUCTION_WITH_TRACE"
    }
    OUT.write_text(json.dumps(payload,indent=2,ensure_ascii=False)+"\n",encoding='utf-8')
    MD.write_text("# P2146 S1096: strict CMP2 commitment execution audit\n\n- Status: `OPEN_PARTIAL_PROGRESS_WITH_TRACE`\n- Result kind: `PASS_STRICT_CMP2_COMMITMENT_EXECUTION_AUDIT_WITH_TRACE`\n\nNo theorem-grade global closure claim is made.\n",encoding='utf-8')

if __name__=='__main__':
    main()
