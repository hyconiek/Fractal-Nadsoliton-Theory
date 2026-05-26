#!/usr/bin/env python3
from __future__ import annotations
import json
from pathlib import Path
from typing import Any

ROOT=Path(__file__).resolve().parent
GEN=ROOT/"generated"
INP=GEN/"p2144_s1094_strict_cmp2_real_ci_stability_interpretation_gate.json"
OUT=GEN/"p2145_s1095_strict_cmp2_next_action_commitment_packet.json"
MD=GEN/"p2145_s1095_strict_cmp2_next_action_commitment_packet.md"


def load(p:Path)->dict[str,Any]:
    if not p.exists(): return {"_missing":str(p)}
    return json.loads(p.read_text(encoding='utf-8'))


def main()->None:
    GEN.mkdir(exist_ok=True)
    p=load(INP)
    gate=(p.get('real_ci_stability_interpretation_gate',{}) or {})
    interpretation=gate.get('interpretation','UNKNOWN')
    ready=bool(gate.get('ready',False))

    if ready:
        commitment=[
            "freeze_non_synthetic_ci_stability_interpretation",
            "archive_cmp2_operational_chain_state",
            "advance_to_next_formal_obstruction_family",
        ]
    else:
        commitment=[
            "deliver_real_cmp2_backend_rows_extension_v1_json",
            "rerun_p2142_executor",
            "refresh_p2144_interpretation_gate",
        ]

    payload={
      "schema_version":"p2145_s1095_v1",
      "packet_id":"P2145","stage_id":"S1095","produced_by":Path(__file__).name,
      "timestamp_utc":"2026-05-26T00:00:00+00:00",
      "status":"OPEN_PARTIAL_PROGRESS_WITH_TRACE",
      "result_kind":"PASS_STRICT_CMP2_NEXT_ACTION_COMMITMENT_PACKET_WITH_TRACE",
      "next_action_commitment":{
        "source_interpretation_gate":str(INP.relative_to(ROOT)),
        "ready":ready,
        "interpretation":interpretation,
        "action_commitment":commitment,
        "scope_limit":"action commitment packet only; not theorem-grade closure",
      },
      "recommended_next_honest_step":{"id":"P2146_candidate","goal":"execute action_commitment and re-run P2142->P2145 for updated flags"},
      "gatekeeper_checks":{"commitment_exported":True,"no_toe_closure_claimed":True,
        "full_d3_covariance_transport_proven":False,"full_cutkosky_closure_proven":False,"c3_theorem_proven":False},
      "global_status":"OPEN_OBSTRUCTION_WITH_TRACE"
    }

    OUT.write_text(json.dumps(payload,indent=2,ensure_ascii=False)+"\n",encoding='utf-8')
    MD.write_text(
      "# P2145 S1095: strict CMP2 next action commitment packet\n\n"
      f"- Status: `{payload['status']}`\n"
      f"- Result kind: `{payload['result_kind']}`\n"
      f"- Interpretation: `{interpretation}`\n\n"
      "No theorem-grade global closure claim is made.\n",encoding='utf-8')

if __name__=='__main__':
    main()
