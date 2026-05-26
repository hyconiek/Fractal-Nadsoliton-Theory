#!/usr/bin/env python3
from __future__ import annotations
import json
from pathlib import Path
from typing import Any

ROOT=Path(__file__).resolve().parent
GEN=ROOT/"generated"
BASE=GEN/"p2126_s1076_strict_cmp2_backend_evidence_weighted_posterior_calibration_audit.json"
EXT=GEN/"cmp2_backend_rows_extension_v1.json"
OUT=GEN/"p2133_s1083_strict_cmp2_real_extension_merge_contract.json"
MD=GEN/"p2133_s1083_strict_cmp2_real_extension_merge_contract.md"


def load(p:Path)->dict[str,Any]:
    if not p.exists(): return {"_missing":str(p)}
    return json.loads(p.read_text(encoding='utf-8'))


def main()->None:
    GEN.mkdir(exist_ok=True)
    b=load(BASE); e=load(EXT)
    rows=((b.get("backend_evidence_weighted_posterior_predictive_calibration_audit",{}) or {}).get("rows") or [])
    ext_rows=e.get("rows") or []
    required={"cmp_bin_index","backend_s","posterior_weights_backend_evidence","posterior_predictive_covered"}
    valid=[r for r in ext_rows if isinstance(r,dict) and required.issubset(r.keys())]
    merged=rows+valid
    merge_ready=len(valid)==len(ext_rows) and len(valid)>0 and len(rows)>0
    payload={
      "schema_version":"p2133_s1083_v1",
      "packet_id":"P2133","stage_id":"S1083","produced_by":Path(__file__).name,
      "timestamp_utc":"2026-05-26T00:00:00+00:00",
      "status":"OPEN_PARTIAL_PROGRESS_WITH_TRACE",
      "result_kind":"PASS_STRICT_CMP2_REAL_EXTENSION_MERGE_CONTRACT_WITH_TRACE" if merge_ready else "OPEN_STRICT_CMP2_REAL_EXTENSION_MERGE_CONTRACT_BLOCKED",
      "real_extension_merge_contract":{
        "base_rows":len(rows),"extension_rows":len(ext_rows),"valid_extension_rows":len(valid),"merged_rows":len(merged),
        "required_row_keys":sorted(required),"merge_ready":merge_ready,
        "scope_limit":"merge-contract/readiness only; no theorem-grade closure"
      },
      "recommended_next_honest_step":{"id":"P2134_candidate","goal":"if merge-ready, execute non-synthetic rerun of P2127/P2128/P2129/P2131 on merged dataset"},
      "gatekeeper_checks":{"merge_ready":merge_ready,"no_toe_closure_claimed":True,
        "full_d3_covariance_transport_proven":False,"full_cutkosky_closure_proven":False,"c3_theorem_proven":False},
      "global_status":"OPEN_OBSTRUCTION_WITH_TRACE"
    }
    OUT.write_text(json.dumps(payload,indent=2,ensure_ascii=False)+"\n",encoding='utf-8')
    MD.write_text("# P2133 S1083: strict CMP2 real-extension merge contract\n\n- Status: `OPEN_PARTIAL_PROGRESS_WITH_TRACE`\n- Result kind: `"+payload['result_kind']+"`\n\nNo theorem-grade closure claim.\n",encoding='utf-8')

if __name__=='__main__':
    main()
