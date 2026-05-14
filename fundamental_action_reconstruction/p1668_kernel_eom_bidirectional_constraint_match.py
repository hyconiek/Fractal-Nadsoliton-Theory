#!/usr/bin/env python3
from __future__ import annotations
import json, math
from pathlib import Path

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "generated" / "p1668_s618_kernel_eom_bidirectional_constraint_match.json"

omega, phi, beta, eta, A = 0.18575, 4.5231, 0.72, 1.31, 1.0
Mpl2 = A*(1+beta)
cR2 = A*beta/(1+eta)
cRic2 = A*beta*eta/(1+eta)
cRiem2 = A*beta*(1+eta)/4
muH2 = A*omega**2
lambdaH = (1+eta**2)/(1+beta)
xiHR = beta/(1+beta)

# forward EOM proxy constraints (symbolic relations as numeric checks)
C1 = muH2 - A*omega**2
C2 = lambdaH*(1+beta) - (1+eta**2)
C3 = xiHR*(1+beta) - beta
C4 = cRic2*(1+eta) - A*beta*eta

# reverse recovery
beta_r = xiHR/(1-xiHR)
A_r = Mpl2/(1+beta_r)
omega_r = math.sqrt(muH2/A_r)
eta_r = math.sqrt(max(lambdaH*(1+beta_r)-1,0.0))

errs = {
  "beta_abs": abs(beta_r-beta),
  "A_abs": abs(A_r-A),
  "omega_abs": abs(omega_r-omega),
  "eta_abs": abs(eta_r-eta),
}
res = {"C1": C1, "C2": C2, "C3": C3, "C4": C4}
local_pass = max(abs(v) for v in list(errs.values()) + [abs(x) for x in res.values()]) < 1e-12

payload = {
 "checkpoint":"P1668_S618_KERNEL_EOM_BIDIRECTIONAL_CONSTRAINT_MATCH",
 "strict_only": True,
 "legacy_bridge_used": False,
 "kernel_input": {"omega":omega,"phi":phi,"beta":beta,"eta":eta,"A":A},
 "coefficients": {"Mpl2":Mpl2,"cR2":cR2,"cRic2":cRic2,"cRiem2":cRiem2,"muH2":muH2,"lambdaH":lambdaH,"xiHR":xiHR},
 "forward_constraint_residuals": res,
 "reverse_recovery": {"beta":beta_r,"A":A_r,"omega":omega_r,"eta":eta_r,"abs_errors":errs},
 "local_pass": local_pass,
 "status":"OPEN_OBLIGATION",
 "open_obligations":["full_metric_and_sm_variational_eom_export","nonproxy_unitarity_theorem","nonproxy_renormalization_closure"]
}
OUT.write_text(json.dumps(payload,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
print(f"Wrote {OUT}")
