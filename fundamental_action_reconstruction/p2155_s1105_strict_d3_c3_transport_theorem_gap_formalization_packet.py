#!/usr/bin/env python3
from __future__ import annotations
import json
from pathlib import Path
from typing import Any

ROOT=Path(__file__).resolve().parent
GEN=ROOT/'generated'
INS=[
 GEN/'p2144_s1094_strict_cmp2_real_ci_stability_interpretation_gate.json',
 GEN/'p2146_s1096_strict_cmp2_commitment_execution_audit.json',
 GEN/'p2147_s1097_strict_cmp2_real_data_required_rerun_checkpoint.json',
 GEN/'p2150_s1100_strict_cmp2_real_data_unlock_attempt_register.json',
]
OUT=GEN/'p2155_s1105_strict_d3_c3_transport_theorem_gap_formalization_packet.json'
MD=GEN/'p2155_s1105_strict_d3_c3_transport_theorem_gap_formalization_packet.md'

def load(p:Path)->dict[str,Any]:
    if not p.exists(): return {'_missing':str(p)}
    return json.loads(p.read_text(encoding='utf-8'))

def main()->None:
    GEN.mkdir(exist_ok=True)
    objs=[load(p) for p in INS]
    checks=[(o.get('gatekeeper_checks',{}) or {}) for o in objs]
    unresolved=[]
    for key in ['full_d3_covariance_transport_proven','c3_theorem_proven']:
        if not all(bool(c.get(key,False)) for c in checks):
            unresolved.append(key)
    payload={
      'schema_version':'p2155_s1105_v1','packet_id':'P2155','stage_id':'S1105','produced_by':Path(__file__).name,
      'timestamp_utc':'2026-05-26T00:00:00+00:00','status':'OPEN_PARTIAL_PROGRESS_WITH_TRACE',
      'result_kind':'OPEN_STRICT_D3_C3_TRANSPORT_THEOREM_GAP_FORMALIZATION_PACKET',
      'd3_c3_transport_theorem_gap_formalization_packet':{
        'source_packets':[str(p.relative_to(ROOT)) for p in INS],
        'unresolved_theorem_objects':unresolved,
        'formal_target':'construct strict-lane transport theorem object that flips both flags to true without legacy-role transfer',
        'scope_limit':'gap formalization only; no theorem-grade closure'
      },
      'recommended_next_honest_step':{'id':'P2156_candidate','goal':'implement symbolic transport witness constructor (FRW/Bianchi-I shared finite-part lock) with executable JSON witness export'},
      'gatekeeper_checks':{'gap_packet_exported':True,'no_toe_closure_claimed':True,'full_d3_covariance_transport_proven':False,'full_cutkosky_closure_proven':False,'c3_theorem_proven':False},
      'global_status':'OPEN_OBSTRUCTION_WITH_TRACE'
    }
    OUT.write_text(json.dumps(payload,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    MD.write_text('# P2155 S1105: strict D3/C3 transport theorem-gap formalization packet\n\n'
                  f"- Unresolved theorem objects: `{unresolved}`\n"
                  '- No theorem-grade closure claim.\n',encoding='utf-8')

if __name__=='__main__':
    main()
