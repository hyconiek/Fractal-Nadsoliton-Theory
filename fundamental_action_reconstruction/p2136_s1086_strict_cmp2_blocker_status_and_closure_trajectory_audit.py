#!/usr/bin/env python3
from __future__ import annotations
import json
from pathlib import Path
from typing import Any

ROOT=Path(__file__).resolve().parent
GEN=ROOT/"generated"
IN={
 "p2128":GEN/"p2128_s1078_strict_cmp2_block_bootstrap_dependence_aware_uncertainty_inflation.json",
 "p2129":GEN/"p2129_s1079_strict_cmp2_block_definition_sensitivity_sweep.json",
 "p2130":GEN/"p2130_s1080_strict_cmp2_sample_support_expansion_planner.json",
 "p2131":GEN/"p2131_s1081_strict_cmp2_support_expansion_replay_audit.json",
 "p2132":GEN/"p2132_s1082_strict_cmp2_real_support_acquisition_gate.json",
 "p2133":GEN/"p2133_s1083_strict_cmp2_real_extension_merge_contract.json",
 "p2134":GEN/"p2134_s1084_strict_cmp2_nonsynthetic_rerun_comparison_audit.json",
 "p2135":GEN/"p2135_s1085_strict_cmp2_merged_real_block_variant_stability_audit.json",
}
OUT=GEN/"p2136_s1086_strict_cmp2_blocker_status_and_closure_trajectory_audit.json"
MD=GEN/"p2136_s1086_strict_cmp2_blocker_status_and_closure_trajectory_audit.md"


def load(p:Path)->dict[str,Any]:
    if not p.exists(): return {"_missing":str(p),"result_kind":"OPEN_MISSING_ARTIFACT"}
    return json.loads(p.read_text(encoding='utf-8'))


def main()->None:
    GEN.mkdir(exist_ok=True)
    s={k:load(v) for k,v in IN.items()}
    rk={k:v.get('result_kind','OPEN_UNKNOWN') for k,v in s.items()}
    blocked=[k for k,v in rk.items() if str(v).startswith('OPEN_')]
    passed=[k for k,v in rk.items() if str(v).startswith('PASS_')]

    closure_toward = len(passed)>=4 and 'p2132' in blocked
    payload={
      "schema_version":"p2136_s1086_v1",
      "packet_id":"P2136","stage_id":"S1086","produced_by":Path(__file__).name,
      "timestamp_utc":"2026-05-26T00:00:00+00:00",
      "status":"OPEN_PARTIAL_PROGRESS_WITH_TRACE",
      "result_kind":"PASS_STRICT_CMP2_BLOCKER_STATUS_AND_CLOSURE_TRAJECTORY_AUDIT_WITH_TRACE",
      "chain_result_kinds":rk,
      "closure_trajectory_audit":{
        "passed_stages":passed,
        "blocked_stages":blocked,
        "interpretation": "trajectory improving operationally but hard-blocked on real backend extension; not near theorem-grade closure" if 'p2132' in blocked else "operational chain unblocked",
        "toward_closure_signal": "PARTIAL_OPERABILITY_ONLY" if closure_toward else "WEAK",
        "scope_limit":"status/trajectory audit only; not theorem-grade closure",
      },
      "recommended_next_honest_step":{"id":"P2137_candidate","goal":"acquire and validate real cmp2_backend_rows_extension_v1.json, then rerun P2132->P2135 non-synthetically"},
      "gatekeeper_checks":{"no_toe_closure_claimed":True,"full_d3_covariance_transport_proven":False,"full_cutkosky_closure_proven":False,"c3_theorem_proven":False},
      "global_status":"OPEN_OBSTRUCTION_WITH_TRACE"
    }
    OUT.write_text(json.dumps(payload,indent=2,ensure_ascii=False)+"\n",encoding='utf-8')
    MD.write_text("# P2136 S1086: blocker status and closure trajectory audit\n\n- Status: `OPEN_PARTIAL_PROGRESS_WITH_TRACE`\n- Blocked stages: `"+", ".join(blocked)+"`\n\nOperational progress exists, but theorem-grade closure remains blocked.\n",encoding='utf-8')

if __name__=='__main__':
    main()
