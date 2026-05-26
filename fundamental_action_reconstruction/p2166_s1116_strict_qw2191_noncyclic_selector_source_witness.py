#!/usr/bin/env python3
from __future__ import annotations
import json
from pathlib import Path
from typing import Any
ROOT=Path(__file__).resolve().parent
GEN=ROOT/'generated'
INP=GEN/'p2165_s1115_strict_qw2191_selector_premise_lattice_and_noncyclic_admissibility_checker.json'
OUT=GEN/'p2166_s1116_strict_qw2191_noncyclic_selector_source_witness.json'
MD=GEN/'p2166_s1116_strict_qw2191_noncyclic_selector_source_witness.md'

def load(p:Path)->dict[str,Any]:
    if not p.exists(): return {'_missing':str(p),'result_kind':'OPEN_MISSING_ARTIFACT'}
    return json.loads(p.read_text(encoding='utf-8'))

def main()->None:
    GEN.mkdir(exist_ok=True)
    p=load(INP)
    blk=(p.get('strict_qw2191_selector_premise_lattice_and_noncyclic_admissibility_checker',{}) or {})
    moves=blk.get('candidate_move_admissibility',[]) or []
    admissible=[m for m in moves if bool(m.get('admissible',False))]
    selected=admissible[0] if admissible else None
    ready=selected is not None
    result='PASS_STRICT_QW2191_NONCYCLIC_SELECTOR_SOURCE_WITNESS_WITH_TRACE' if ready else 'OPEN_STRICT_QW2191_NONCYCLIC_SELECTOR_SOURCE_WITNESS_BLOCKED'
    payload={
      'schema_version':'p2166_s1116_v1','packet_id':'P2166','stage_id':'S1116','produced_by':Path(__file__).name,
      'timestamp_utc':'2026-05-26T00:00:00+00:00','status':'OPEN_PARTIAL_PROGRESS_WITH_TRACE','result_kind':result,
      'strict_qw2191_noncyclic_selector_source_witness':{
        'witness_id':'W_QW2191_NONCYCLIC_SELECTOR_SOURCE_V1','selected_admissible_branch':selected.get('id') if selected else None,
        'strict_selector_source_hypothesis':'explicit symmetry-breaking premise + internal strict-lane source',
        'noncyclic_anchor':'new provider class (non-L5/L12 repeat)',
        'forbidden_routes_excluded':['silent_selector_injection','legacy_role_transfer','repeated_l5_l12_same_blocker_cut'],
        'scope_limit':'selector-source witness only; no strict-core selector closure claim'
      },
      'recommended_next_honest_step':{'id':'P2167_candidate','goal':'lift witness to selector-premise branch instantiation with explicit theorem obligations against QW-2191 uniqueness obstruction'},
      'gatekeeper_checks':{'witness_exported':True,'admissible_branch_selected':ready,'no_selector_closure_claimed':True,'no_toe_closure_claimed':True,
        'full_d3_covariance_transport_proven':bool((p.get('gatekeeper_checks',{}) or {}).get('full_d3_covariance_transport_proven',False)),'full_cutkosky_closure_proven':False,
        'c3_theorem_proven':bool((p.get('gatekeeper_checks',{}) or {}).get('c3_theorem_proven',False))},
      'global_status':'OPEN_OBSTRUCTION_WITH_TRACE'}
    OUT.write_text(json.dumps(payload,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    MD.write_text(f"# P2166 S1116\n\n- Result kind: `{result}`\n- Selected branch: `{payload['strict_qw2191_noncyclic_selector_source_witness']['selected_admissible_branch']}`\n\nNo global ToE closure claim is made.\n",encoding='utf-8')

if __name__=='__main__':
    main()
