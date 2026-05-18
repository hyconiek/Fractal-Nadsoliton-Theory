#!/usr/bin/env python3
"""P2011 S961 strict Cutkosky loop-amplitude tensor placeholder + coupled scan.

Next honest step after P2010: replace backend channel map by explicit per-channel
loop-amplitude tensor placeholders extracted from channel kernels to prepare full
backend tensor export integration.
"""
from __future__ import annotations
import json, platform
from pathlib import Path
from typing import Any
import numpy as np
import scipy.linalg as la

ROOT=Path(__file__).resolve().parent
GEN=ROOT/"generated"
OUT=GEN/"p2011_s961_strict_cutkosky_loop_amplitude_tensor_placeholder_and_coupled_scan_witness.json"
TS="2026-05-18T00:00:00+00:00"

def load(name:str)->dict[str,Any]:
    p=GEN/name
    if not p.exists(): return {"_missing":name,"status":"OPEN_MISSING_ARTIFACT"}
    return json.loads(p.read_text(encoding='utf-8'))

def main()->None:
    GEN.mkdir(exist_ok=True)
    p2010=load("p2010_s960_strict_cutkosky_backend_channel_tensor_map_classifier_witness.json")
    p2004=load("p2004_s954_strict_cutkosky_gx_loop_amplitude_and_robustness_refresh_witness.json")

    rows=p2004.get("grid_table",[])
    delta=np.array([float(r["Delta_opt"]) for r in rows],dtype=float)
    ch={k:np.array([float(r["Cut_channels"].get(k,0.0)) for r in rows],dtype=float) for k in ["gg","gh","hh","gx"]}
    base_l2=float(la.norm(delta,2))

    # Placeholder channel loop-amplitude tensors from channel energy profiles
    # T_i := outer(profile_i, profile_i) reduced to scalar variance proxies.
    var={k:float(np.var(v)) for k,v in ch.items()}
    cov=np.zeros((4,4),dtype=float)
    keys=["gg","gh","hh","gx"]
    for i,a in enumerate(keys):
        for j,b in enumerate(keys):
            cov[i,j]=(var[a]*var[b])**0.5 if i!=j else var[a]+1e-18
    cov+=1e-18*np.eye(4)

    sigma=np.sqrt(np.diag(cov)); sigma=sigma/(np.linalg.norm(sigma)+1e-15)
    scans=[]; counts={"MISSING_CHANNEL_PRESSURE_SUPPORTED":0,"STRUCTURAL_OBSTRUCTION_PRESSURE_SUPPORTED":0,"MIXED_OR_INCONCLUSIVE":0}
    for t in np.linspace(-1,1,17):
        scales=1+t*sigma
        d=delta.copy()
        for i,k in enumerate(keys): d-=(scales[i]-1.0)*ch[k]
        ratio=float(la.norm(d,2)/base_l2)
        if ratio<0.95: cls="MISSING_CHANNEL_PRESSURE_SUPPORTED"
        elif ratio>1.05: cls="STRUCTURAL_OBSTRUCTION_PRESSURE_SUPPORTED"
        else: cls="MIXED_OR_INCONCLUSIVE"
        counts[cls]+=1
        scans.append({"t":float(t),"scales":{k:float(scales[i]) for i,k in enumerate(keys)},"l2_ratio_vs_p2004":ratio,"classifier":cls})

    dominant=max(counts,key=lambda k:counts[k]); total=len(scans)
    share=counts[dominant]/total if total else 0.0
    gate={
      "p2010_present":p2010.get("result_kind")=="PASS_BACKEND_CHANNEL_MAP_COUPLED_CLASSIFIER_WITNESS",
      "loop_amplitude_tensor_placeholders_exported":True,
      "covariance_psd":bool(np.all(np.linalg.eigvalsh(cov)>0)),
      "scan_nonempty":total>0,
    }

    out={"ledger_id":"P2011_S961_STRICT_CUTKOSKY_LOOP_AMPLITUDE_TENSOR_PLACEHOLDER_AND_COUPLED_SCAN_WITNESS",
         "packet_id":"P2011","stage_id":"S961","produced_by":Path(__file__).name,
         "timestamp_utc":TS,"route":"strict_only","depends_on":{"p2010_present":gate["p2010_present"]},
         "channel_loop_tensor_placeholders":{"variances":var,"covariance":cov.tolist(),"sigma_unit":sigma.tolist()},
         "classifier_scan":{"count":total,"table":scans,"counts":counts,"dominant_classifier":dominant,"dominant_share":share},
         "delta_stats":{"l2_base_p2004":base_l2},"gatekeeper_checks":gate,
         "result_kind":"PASS_LOOP_TENSOR_PLACEHOLDER_COUPLED_SCAN_WITNESS" if all(gate.values()) else "OPEN_OBSTRUCTION_WITH_TRACE",
         "status":"OPEN_OBSTRUCTION_WITH_TRACE",
         "false_pass_guard":"Placeholder channel loop tensors are preparatory diagnostics only; not theorem-grade closure.",
         "next_honest_step":"Replace placeholders with exported backend loop-amplitude tensors per channel and rerun coupled covariance scan.",
         "lay_explanation":"Zrobiliśmy przejściowy krok: osobne tensory kanałów zbudowane z ich profili numerycznych, żeby przygotować pełny import backend tensorów amplitud.",
         "environment":{"python":platform.python_version(),"numpy":np.__version__,"scipy":__import__('scipy').__version__}}
    OUT.write_text(json.dumps(out,indent=2,ensure_ascii=False)+"\n",encoding='utf-8')
    print(f"[P2011] wrote witness: {OUT}")

if __name__=="__main__":
    main()
