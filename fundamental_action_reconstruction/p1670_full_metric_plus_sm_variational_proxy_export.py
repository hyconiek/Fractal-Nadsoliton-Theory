#!/usr/bin/env python3
from __future__ import annotations
import json, math
from pathlib import Path

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "generated" / "p1670_s620_full_metric_plus_sm_variational_proxy_export.json"

points=[(0.60,1.10,0.18,1.0,0.00),(0.72,1.31,0.18575,1.0,0.08),(0.85,1.45,0.19,1.0,0.12)]
rows=[]
for beta,eta,omega,A,R in points:
    Mpl2=A*(1+beta); cR2=A*beta/(1+eta); cRic2=A*beta*eta/(1+eta); cRiem2=A*beta*(1+eta)/4
    muH2=A*omega**2; lambdaH=(1+eta**2)/(1+beta); xiHR=beta/(1+beta); ZA=1+beta**2
    alphaQ=cR2+cRic2+cRiem2
    # sample field values
    h=1.0+0.1*beta; boxh=0.2
    scalar_res=boxh + muH2*h + lambdaH*h**3 - 2*xiHR*R*h
    # metric proxy residual with T/Mpl2 using simple Ttrace sample
    T_over_M = 0.3*(1+beta)/Mpl2
    metric_res=(1+alphaQ*R)-T_over_M
    gauge_res=ZA*0.25-0.25*ZA  # constructed identity proxy
    rows.append({
      "input":{"beta":beta,"eta":eta,"omega":omega,"A":A,"R":R},
      "coeff":{"Mpl2":Mpl2,"cR2":cR2,"cRic2":cRic2,"cRiem2":cRiem2,"muH2":muH2,"lambdaH":lambdaH,"xiHR":xiHR,"ZA":ZA},
      "eom_residuals":{"metric_proxy":metric_res,"scalar_proxy":scalar_res,"gauge_proxy":gauge_res}
    })

payload={
 "checkpoint":"P1670_S620_FULL_METRIC_PLUS_SM_VARIATIONAL_PROXY_EXPORT",
 "strict_only":True,
 "legacy_bridge_used":False,
 "num_points":len(rows),
 "rows":rows,
 "status":"OPEN_OBLIGATION",
 "open_obligations":["full_tensor_variation_nonproxy_export","spin2_spin0_unitarity_theorem","qg_renormalization_background_independence_proof"]
}
OUT.write_text(json.dumps(payload,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
print(f"Wrote {OUT}")
