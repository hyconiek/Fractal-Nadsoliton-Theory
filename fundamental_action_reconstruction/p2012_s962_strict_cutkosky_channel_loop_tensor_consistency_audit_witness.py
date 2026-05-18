#!/usr/bin/env python3
"""P2012 S962 strict Cutkosky channel loop-tensor consistency audit witness.

Next honest step after P2011: audit consistency between P2010 backend-channel-map
covariance and P2011 channel-loop-tensor placeholders, then build blended
covariance scan to localize map-vs-tensor tension.
"""
from __future__ import annotations
import json, platform
from pathlib import Path
from typing import Any
import numpy as np
import scipy.linalg as la

ROOT=Path(__file__).resolve().parent
GEN=ROOT/"generated"
OUT=GEN/"p2012_s962_strict_cutkosky_channel_loop_tensor_consistency_audit_witness.json"
TS="2026-05-18T00:00:00+00:00"

def load(name:str)->dict[str,Any]:
    p=GEN/name
    if not p.exists(): return {"_missing":name,"status":"OPEN_MISSING_ARTIFACT"}
    return json.loads(p.read_text(encoding='utf-8'))

def main()->None:
    GEN.mkdir(exist_ok=True)
    p2010=load("p2010_s960_strict_cutkosky_backend_channel_tensor_map_classifier_witness.json")
    p2011=load("p2011_s961_strict_cutkosky_loop_amplitude_tensor_placeholder_and_coupled_scan_witness.json")
    p2004=load("p2004_s954_strict_cutkosky_gx_loop_amplitude_and_robustness_refresh_witness.json")

    C4_map=np.array(p2010.get("backend_tensor_objects",{}).get("C4_channel_covariance",[]),dtype=float)
    C4_loop=np.array(p2011.get("channel_loop_tensor_placeholders",{}).get("covariance",[]),dtype=float)

    rows=p2004.get("grid_table",[])
    delta=np.array([float(r["Delta_opt"]) for r in rows],dtype=float)
    ch={k:np.array([float(r["Cut_channels"].get(k,0.0)) for r in rows],dtype=float) for k in ["gg","gh","hh","gx"]}
    base_l2=float(la.norm(delta,2))

    # consistency audit metrics
    fro_map=float(la.norm(C4_map,'fro')) if C4_map.size else float('nan')
    fro_loop=float(la.norm(C4_loop,'fro')) if C4_loop.size else float('nan')
    fro_diff=float(la.norm(C4_map-C4_loop,'fro')) if C4_map.shape==C4_loop.shape else float('nan')
    rel_diff=fro_diff/(fro_map+1e-18) if np.isfinite(fro_diff) else float('nan')

    # blended covariance family C(λ)=λ*C_map+(1-λ)*C_loop
    scans=[]; counts={"MISSING_CHANNEL_PRESSURE_SUPPORTED":0,"STRUCTURAL_OBSTRUCTION_PRESSURE_SUPPORTED":0,"MIXED_OR_INCONCLUSIVE":0}
    for lam in np.linspace(0,1,11):
        C=lam*C4_map+(1-lam)*C4_loop
        C+=1e-18*np.eye(4)
        sigma=np.sqrt(np.diag(C)); sigma=sigma/(np.linalg.norm(sigma)+1e-15)
        # evaluate at t=+0.5 and t=-0.5 to probe both pressure directions
        for t in (-0.5,0.5):
            scales=1+t*sigma
            d=delta.copy()
            for i,k in enumerate(["gg","gh","hh","gx"]): d-=(scales[i]-1.0)*ch[k]
            ratio=float(la.norm(d,2)/base_l2)
            if ratio<0.95: cls="MISSING_CHANNEL_PRESSURE_SUPPORTED"
            elif ratio>1.05: cls="STRUCTURAL_OBSTRUCTION_PRESSURE_SUPPORTED"
            else: cls="MIXED_OR_INCONCLUSIVE"
            counts[cls]+=1
            scans.append({"lambda":float(lam),"t":float(t),"l2_ratio_vs_p2004":ratio,"classifier":cls})

    total=len(scans); dominant=max(counts,key=lambda k:counts[k]); share=counts[dominant]/total if total else 0.0
    gate={
      "p2010_present":p2010.get("result_kind")=="PASS_BACKEND_CHANNEL_MAP_COUPLED_CLASSIFIER_WITNESS",
      "p2011_present":p2011.get("result_kind")=="PASS_LOOP_TENSOR_PLACEHOLDER_COUPLED_SCAN_WITNESS",
      "covariance_shapes_match":C4_map.shape==(4,4) and C4_loop.shape==(4,4),
      "scan_nonempty":total>0,
      "share_bounded":0<=share<=1,
    }

    out={"ledger_id":"P2012_S962_STRICT_CUTKOSKY_CHANNEL_LOOP_TENSOR_CONSISTENCY_AUDIT_WITNESS",
         "packet_id":"P2012","stage_id":"S962","produced_by":Path(__file__).name,
         "timestamp_utc":TS,"route":"strict_only","depends_on":{"p2010_present":gate["p2010_present"],"p2011_present":gate["p2011_present"]},
         "consistency_audit":{"fro_map":fro_map,"fro_loop":fro_loop,"fro_diff":fro_diff,"relative_diff_vs_map":rel_diff},
         "blended_scan":{"count":total,"table":scans,"counts":counts,"dominant_classifier":dominant,"dominant_share":share},
         "delta_stats":{"l2_base_p2004":base_l2},"gatekeeper_checks":gate,
         "result_kind":"PASS_CHANNEL_LOOP_TENSOR_CONSISTENCY_AUDIT_WITNESS" if all(gate.values()) else "OPEN_OBSTRUCTION_WITH_TRACE",
         "status":"OPEN_OBSTRUCTION_WITH_TRACE",
         "false_pass_guard":"Consistency audit is diagnostics-only and does not claim theorem-grade closure.",
         "next_honest_step":"Import explicit backend loop-amplitude tensors and replace placeholders before repeating blended scans.",
         "lay_explanation":"Sprawdziliśmy, jak bardzo różnią się dwa sposoby budowy kowariancji kanałów i jak to wpływa na wynik klasyfikacji.",
         "environment":{"python":platform.python_version(),"numpy":np.__version__,"scipy":__import__('scipy').__version__}}
    OUT.write_text(json.dumps(out,indent=2,ensure_ascii=False)+"\n",encoding='utf-8')
    print(f"[P2012] wrote witness: {OUT}")

if __name__=="__main__":
    main()
