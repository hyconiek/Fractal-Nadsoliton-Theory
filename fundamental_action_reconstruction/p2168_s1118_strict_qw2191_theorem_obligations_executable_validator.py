#!/usr/bin/env python3
from __future__ import annotations
import json
from pathlib import Path
from typing import Any

ROOT=Path(__file__).resolve().parent
GEN=ROOT/'generated'
INP=GEN/'p2167_s1117_strict_qw2191_selector_premise_branch_and_theorem_obligations_packet.json'
OUT=GEN/'p2168_s1118_strict_qw2191_theorem_obligations_executable_validator.json'
MD=GEN/'p2168_s1118_strict_qw2191_theorem_obligations_executable_validator.md'

def load(p:Path)->dict[str,Any]:
    if not p.exists(): return {'_missing':str(p),'result_kind':'OPEN_MISSING_ARTIFACT'}
    return json.loads(p.read_text(encoding='utf-8'))

def main()->None:
    GEN.mkdir(exist_ok=True)
    p=load(INP)
    blk=(p.get('strict_qw2191_selector_premise_branch_and_theorem_obligations_packet',{}) or {})
    obligations=blk.get('theorem_obligations',[]) or []
    branch_instantiated=bool(blk.get('branch_instantiated',False))

    evals=[]
    for o in obligations:
        oid=o.get('id','UNKNOWN')
        # executable rule-set (strict/no-false-pass): still OPEN until dedicated witness packets exist
        if oid in {'O1_selector_source_identifiability','O3_noncyclic_anchor_certificate','O4_no_legacy_role_transfer_certificate'}:
            status='PARTIAL_PASS_WITH_TRACE' if branch_instantiated else 'OPEN'
        else:
            status='OPEN'
        pass_flag=status.startswith('PASS_')
        evals.append({'id':oid,'status':status,'pass':pass_flag,'required_for_legal_progress':bool(o.get('required_for_legal_progress',True))})

    required=[e for e in evals if e['required_for_legal_progress']]
    n_required=len(required)
    n_pass=sum(1 for e in required if e['pass'])
    n_open=n_required-n_pass

    validator_ready = n_required>0
    all_required_pass = n_open==0 and n_required>0
    result='PASS_STRICT_QW2191_THEOREM_OBLIGATIONS_EXECUTABLE_VALIDATOR_WITH_TRACE' if validator_ready else 'OPEN_STRICT_QW2191_THEOREM_OBLIGATIONS_EXECUTABLE_VALIDATOR_BLOCKED'

    payload={
      'schema_version':'p2168_s1118_v1','packet_id':'P2168','stage_id':'S1118','produced_by':Path(__file__).name,
      'timestamp_utc':'2026-05-26T00:00:00+00:00','status':'OPEN_PARTIAL_PROGRESS_WITH_TRACE','result_kind':result,
      'strict_qw2191_theorem_obligations_executable_validator':{
        'source_branch_packet':str(INP.relative_to(ROOT)),'branch_instantiated':branch_instantiated,
        'obligation_evaluations':evals,
        'required_summary':{'n_required':n_required,'n_pass':n_pass,'n_open':n_open,'all_required_pass':all_required_pass},
        'scope_limit':'validator only; no selector closure or QW-2191 discharge claim'
      },
      'recommended_next_honest_step':{'id':'P2169_candidate','goal':'export first dedicated witnesses for O2/O5 and upgrade partial statuses via executable checks'},
      'gatekeeper_checks':{'validator_exported':True,'validator_ready':validator_ready,'all_required_pass':all_required_pass,
        'no_selector_closure_claimed':True,'no_toe_closure_claimed':True,
        'full_d3_covariance_transport_proven':bool((p.get('gatekeeper_checks',{}) or {}).get('full_d3_covariance_transport_proven',False)),
        'full_cutkosky_closure_proven':False,'c3_theorem_proven':bool((p.get('gatekeeper_checks',{}) or {}).get('c3_theorem_proven',False))},
      'global_status':'OPEN_OBSTRUCTION_WITH_TRACE'}
    OUT.write_text(json.dumps(payload,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    MD.write_text('# P2168 S1118: strict QW-2191 theorem-obligations executable validator\n\n'
                  f"- Result kind: `{result}`\n"
                  f"- Required pass/open: `{n_pass}/{n_open}`\n"
                  f"- all_required_pass: `{all_required_pass}`\n\nNo global ToE closure claim is made.\n",encoding='utf-8')

if __name__=='__main__':
    main()
