#!/usr/bin/env python3
from __future__ import annotations
import json
from pathlib import Path
from typing import Any

ROOT=Path(__file__).resolve().parent
GEN=ROOT/"generated"
IN={
 'p2142':GEN/'p2142_s1092_strict_cmp2_real_data_handoff_and_rerun_executor.json',
 'p2143':GEN/'p2143_s1093_strict_cmp2_real_ci_stability_readiness_memo.json',
 'p2144':GEN/'p2144_s1094_strict_cmp2_real_ci_stability_interpretation_gate.json',
 'p2147':GEN/'p2147_s1097_strict_cmp2_real_data_required_rerun_checkpoint.json',
 'p2150':GEN/'p2150_s1100_strict_cmp2_real_data_unlock_attempt_register.json',
}
OUT=GEN/'p2154_s1104_strict_cmp2_internal_milestone_freeze_packet.json'
MD=GEN/'p2154_s1104_strict_cmp2_internal_milestone_freeze_packet.md'

def load(p:Path)->dict[str,Any]:
    if not p.exists(): return {'_missing':str(p),'result_kind':'OPEN_MISSING_ARTIFACT'}
    return json.loads(p.read_text(encoding='utf-8'))

def main()->None:
    GEN.mkdir(exist_ok=True)
    rk={k:load(v).get('result_kind','OPEN_MISSING_ARTIFACT') for k,v in IN.items()}
    ready=all(str(rk[k]).startswith('PASS_') for k in ['p2142','p2143','p2144','p2147','p2150'])
    payload={
      'schema_version':'p2154_s1104_v1','packet_id':'P2154','stage_id':'S1104','produced_by':Path(__file__).name,
      'timestamp_utc':'2026-05-26T00:00:00+00:00','status':'OPEN_PARTIAL_PROGRESS_WITH_TRACE',
      'result_kind':'PASS_STRICT_CMP2_INTERNAL_MILESTONE_FREEZE_PACKET_WITH_TRACE' if ready else 'OPEN_STRICT_CMP2_INTERNAL_MILESTONE_FREEZE_PACKET_BLOCKED',
      'internal_milestone_freeze_packet':{
        'milestone_name':'CMP2_NONSYNTHETIC_INTERNAL_MILESTONE_V1',
        'component_result_kinds':rk,
        'ready_for_internal_freeze':ready,
        'closure_label':'NOT_THEOREM_GRADE_TOE_CLOSURE',
        'scope_limit':'internal CMP2 process milestone only; not global theorem-grade closure'
      },
      'recommended_next_honest_step':{'id':'P2155_candidate','goal':'run D3/C3 transport theorem-gap formalization packet on strict lane and quantify unresolved theorem objects'},
      'gatekeeper_checks':{'milestone_packet_exported':True,'no_toe_closure_claimed':True,'full_d3_covariance_transport_proven':False,'full_cutkosky_closure_proven':False,'c3_theorem_proven':False},
      'global_status':'OPEN_OBSTRUCTION_WITH_TRACE'
    }
    OUT.write_text(json.dumps(payload,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    MD.write_text('# P2154 S1104: strict CMP2 internal milestone freeze packet\n\n'
                  f"- Result kind: `{payload['result_kind']}`\n"
                  f"- Ready for internal freeze: `{ready}`\n"
                  '- Closure label: `NOT_THEOREM_GRADE_TOE_CLOSURE`\n',encoding='utf-8')

if __name__=='__main__':
    main()
