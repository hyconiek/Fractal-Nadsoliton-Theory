#!/usr/bin/env python3
from __future__ import annotations
import json, subprocess, sys
from pathlib import Path
from typing import Any

ROOT=Path(__file__).resolve().parent
GEN=ROOT/"generated"
EXT=GEN/"cmp2_backend_rows_extension_v1.json"
OUT=GEN/"p2147_s1097_strict_cmp2_real_data_required_rerun_checkpoint.json"
MD=GEN/"p2147_s1097_strict_cmp2_real_data_required_rerun_checkpoint.md"

PIPE=[
 ROOT/"p2132_s1082_strict_cmp2_real_support_acquisition_gate.py",
 ROOT/"p2133_s1083_strict_cmp2_real_extension_merge_contract.py",
 ROOT/"p2134_s1084_strict_cmp2_nonsynthetic_rerun_comparison_audit.py",
 ROOT/"p2135_s1085_strict_cmp2_merged_real_block_variant_stability_audit.py",
 ROOT/"p2141_s1091_strict_cmp2_flag_outcome_matrix.py",
 ROOT/"p2146_s1096_strict_cmp2_commitment_execution_audit.py",
]

def load(p:Path)->dict[str,Any]:
    if not p.exists(): return {"_missing":str(p),"result_kind":"OPEN_MISSING_ARTIFACT"}
    return json.loads(p.read_text(encoding='utf-8'))

def run(p:Path):
    try:
        subprocess.run([sys.executable,str(p)],check=True)
        return True,'OK'
    except subprocess.CalledProcessError as e:
        return False,f'FAILED({e.returncode})'

def main()->None:
    GEN.mkdir(exist_ok=True)
    log=[]
    for s in PIPE:
        ok,msg=run(s); log.append({"script":s.name,"ok":ok,"message":msg})
    r={
      'p2132':load(GEN/'p2132_s1082_strict_cmp2_real_support_acquisition_gate.json').get('result_kind'),
      'p2133':load(GEN/'p2133_s1083_strict_cmp2_real_extension_merge_contract.json').get('result_kind'),
      'p2134':load(GEN/'p2134_s1084_strict_cmp2_nonsynthetic_rerun_comparison_audit.json').get('result_kind'),
      'p2135':load(GEN/'p2135_s1085_strict_cmp2_merged_real_block_variant_stability_audit.json').get('result_kind'),
      'p2141':load(GEN/'p2141_s1091_strict_cmp2_flag_outcome_matrix.json').get('result_kind'),
      'p2146':load(GEN/'p2146_s1096_strict_cmp2_commitment_execution_audit.json').get('result_kind'),
    }
    unblocked=all(str(r[k]).startswith('PASS_') for k in ['p2132','p2133','p2134','p2135'])
    payload={
      "schema_version":"p2147_s1097_v1","packet_id":"P2147","stage_id":"S1097","produced_by":Path(__file__).name,
      "timestamp_utc":"2026-05-26T00:00:00+00:00","status":"OPEN_PARTIAL_PROGRESS_WITH_TRACE",
      "result_kind":"PASS_STRICT_CMP2_REAL_DATA_REQUIRED_RERUN_CHECKPOINT_WITH_TRACE" if unblocked else "OPEN_STRICT_CMP2_REAL_DATA_REQUIRED_RERUN_CHECKPOINT_BLOCKED",
      "real_data_required_rerun_checkpoint":{
        "extension_file":str(EXT.relative_to(ROOT)),"extension_present":EXT.exists(),"pipeline_log":log,
        "result_kinds":r,"nonsynthetic_unblocked":unblocked,
        "scope_limit":"checkpoint only; requires real external data for unlock"
      },
      "recommended_next_honest_step":{"id":"P2148_candidate","goal":"provide real cmp2_backend_rows_extension_v1.json then rerun P2147 and interpret only after p2132-p2135 PASS"},
      "gatekeeper_checks":{"extension_present":EXT.exists(),"p2132_pass":str(r['p2132']).startswith('PASS_'),"p2133_pass":str(r['p2133']).startswith('PASS_'),"p2134_pass":str(r['p2134']).startswith('PASS_'),"p2135_pass":str(r['p2135']).startswith('PASS_'),"no_toe_closure_claimed":True,
      "full_d3_covariance_transport_proven":False,"full_cutkosky_closure_proven":False,"c3_theorem_proven":False},
      "global_status":"OPEN_OBSTRUCTION_WITH_TRACE"
    }
    OUT.write_text(json.dumps(payload,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    MD.write_text(f"# P2147 S1097\n\n- Result kind: `{payload['result_kind']}`\n- Extension present: `{EXT.exists()}`\n",encoding='utf-8')

if __name__=='__main__':
    main()
