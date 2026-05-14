#!/usr/bin/env python3
from __future__ import annotations
import json, math
from pathlib import Path

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "generated" / "p1669_s619_multipoint_full_l_total_variational_proxy.json"

points = [
    {"beta":0.60,"eta":1.10,"omega":0.18,"A":1.0},
    {"beta":0.72,"eta":1.31,"omega":0.18575,"A":1.0},
    {"beta":0.85,"eta":1.45,"omega":0.19,"A":1.0},
]

rows=[]
max_err=0.0
passes=0
for p in points:
    beta,eta,omega,A = p["beta"],p["eta"],p["omega"],p["A"]
    Mpl2=A*(1+beta); muH2=A*omega**2; lambdaH=(1+eta**2)/(1+beta); xiHR=beta/(1+beta)
    cRic2=A*beta*eta/(1+eta)
    # forward residuals
    C1=muH2-A*omega**2
    C2=lambdaH*(1+beta)-(1+eta**2)
    C3=xiHR*(1+beta)-beta
    C4=cRic2*(1+eta)-A*beta*eta
    # reverse
    beta_r=xiHR/(1-xiHR); A_r=Mpl2/(1+beta_r); omega_r=math.sqrt(muH2/A_r); eta_r=math.sqrt(max(lambdaH*(1+beta_r)-1,0.0))
    err=max(abs(beta_r-beta),abs(A_r-A),abs(omega_r-omega),abs(eta_r-eta),abs(C1),abs(C2),abs(C3),abs(C4))
    ok=err<1e-12
    if ok: passes+=1
    max_err=max(max_err,err)
    rows.append({"input":p,"reverse":{"beta":beta_r,"A":A_r,"omega":omega_r,"eta":eta_r},"max_err":err,"local_pass":ok})

payload={
  "checkpoint":"P1669_S619_MULTIPOINT_FULL_L_TOTAL_VARIATIONAL_PROXY",
  "strict_only":True,
  "legacy_bridge_used":False,
  "num_points":len(points),
  "num_local_pass":passes,
  "max_error":max_err,
  "rows":rows,
  "status":"OPEN_OBLIGATION",
  "open_obligations":["full_metric_tensor_variation_export","nonproxy_multifield_eom_theorem","qg_renormalization_unitarity_background_independence_closure"]
}
OUT.write_text(json.dumps(payload,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
print(f"Wrote {OUT}")
