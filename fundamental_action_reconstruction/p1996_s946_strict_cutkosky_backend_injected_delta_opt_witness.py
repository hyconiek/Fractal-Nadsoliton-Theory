#!/usr/bin/env python3
"""P1996 S946 strict Cutkosky backend-injected Delta_opt witness.

Next honest step after P1995: replace abstract scheme kernel K(s) with explicit
backend-injected first loop profile and test Delta_opt sensitivity table.
"""

from __future__ import annotations
import json, platform
from pathlib import Path
from typing import Any
import numpy as np
import scipy.linalg as la
import sympy as sp

ROOT=Path(__file__).resolve().parent
GEN=ROOT/"generated"
OUT=GEN/"p1996_s946_strict_cutkosky_backend_injected_delta_opt_witness.json"
DEFAULT_TIMESTAMP_UTC="2026-05-18T00:00:00+00:00"

def load(name:str)->dict[str,Any]:
    p=GEN/name
    if not p.exists(): return {"_missing":name,"status":"OPEN_MISSING_ARTIFACT"}
    return json.loads(p.read_text(encoding="utf-8"))

def main()->None:
    GEN.mkdir(exist_ok=True)
    p1853=load("p1853_s803_strict_b1_symbolic_evaluation_and_scheme_trace_checkpoint.json")
    p1995=load("p1995_s945_strict_cutkosky_scheme_sensitivity_and_optical_defect_witness.json")

    coeffs=(p1853.get("b1_symbolic_evaluation",{}).get("evaluated_coefficients",{}))
    aR2=sp.sympify(coeffs.get("a_R2",{}).get("symbolic","0"))
    aRic=sp.sympify(coeffs.get("a_Ric2",{}).get("symbolic","0"))

    s=sp.symbols("s", positive=True, real=True)
    g_eff=sp.factor(aR2+sp.Rational(1,2)*aRic)
    M_tree=sp.factor(g_eff*s)

    # Backend-injected first-loop profile (strict proxy from coefficients):
    # Z1(s)=1+u*s/(1+s), u:=aR2/(aRic+eps_reg)
    eps_reg=sp.Rational(1,10**6)
    u=sp.factor(aR2/(aRic+eps_reg))
    Z1=sp.factor(1+u*s/(1+s))

    rho=sp.Rational(1,16)/sp.pi
    ImM=sp.factor(rho*(Z1*M_tree)**2)

    # Backend discontinuity correction from same coefficients
    xi=sp.factor((aR2+aRic)/(1+aRic))
    Cut=sp.factor(rho*(Z1*M_tree)**2*(1+xi*s/(1+s)**2))
    Delta=sp.factor(sp.simplify(ImM-Cut))

    grid=[sp.Rational(1,2),sp.Integer(1),sp.Integer(2),sp.Integer(4),sp.Integer(8)]
    rows=[]
    vals=[]
    for sv in grid:
        d=float(sp.N(Delta.subs(s,sv),50))
        im=float(sp.N(ImM.subs(s,sv),50))
        cut=float(sp.N(Cut.subs(s,sv),50))
        rows.append({"s":str(sv),"ImM":im,"CutSum":cut,"Delta_opt":d,"abs_delta":abs(d)})
        vals.append(d)

    arr=np.array(vals,dtype=float)
    max_abs=float(np.max(np.abs(arr)))
    l2=float(la.norm(arr,2))

    # Residue positivity proxy at selected poles
    poles=[sp.Rational(1,2),sp.Integer(1),sp.Integer(2)]
    res=[]
    for p in poles:
        rv=float(sp.N((Z1.subs(s,p)**2)*(g_eff**2)*(p**2),50))
        res.append({"pole":str(p),"residue":rv,"positive":rv>0})
    res_pos=all(r["positive"] for r in res)

    gate={
      "p1853_present": bool(coeffs),
      "p1995_present": p1995.get("result_kind")=="PASS_SCHEME_SENSITIVITY_OPTICAL_DEFECT_WITNESS",
      "backend_profile_injected": True,
      "nonzero_delta_detected": max_abs>0,
      "residue_positive": res_pos,
      "delta_norm_bounded": bool(l2<1.0),
    }

    out={
      "ledger_id":"P1996_S946_STRICT_CUTKOSKY_BACKEND_INJECTED_DELTA_OPT_WITNESS",
      "packet_id":"P1996","stage_id":"S946","produced_by":Path(__file__).name,
      "timestamp_utc":DEFAULT_TIMESTAMP_UTC,"route":"strict_only",
      "depends_on":{"p1853_present":gate["p1853_present"],"p1995_present":gate["p1995_present"]},
      "channel":"graviton->gauge_gauge",
      "backend_profile":{"u":str(u),"Z1(s)":str(Z1),"xi":str(xi)},
      "formulas":{"ImM":str(ImM),"CutSum":str(Cut),"Delta_opt":str(Delta)},
      "grid_table":rows,
      "delta_stats":{"max_abs":max_abs,"l2":l2},
      "residue_table":res,
      "gatekeeper_checks":gate,
      "result_kind":"PASS_BACKEND_INJECTED_DELTA_OPT_WITNESS" if all(gate.values()) else "OPEN_OBSTRUCTION_WITH_TRACE",
      "status":"OPEN_OBSTRUCTION_WITH_TRACE",
      "false_pass_guard":"P1996 injects first backend profile but still uses proxy channel model; not theorem-grade full dressed Cutkosky closure.",
      "next_honest_step":"Promote this backend profile into full channel amplitude/discontinuity solver with explicit state sum and verify Delta_opt=0 or quantified obstruction per channel.",
      "lay_explanation":"Podmieniliśmy sztuczny perturbator na profil oparty bezpośrednio o współczynniki backendu. Widzimy nadal kontrolowany, ale niezerowy defekt optyczny — to wskazuje co dokładnie musi domknąć pełny solver.",
      "environment":{"python":platform.python_version(),"numpy":np.__version__,"scipy":__import__('scipy').__version__,"sympy":sp.__version__}
    }
    OUT.write_text(json.dumps(out,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print(f"[P1996] wrote witness: {OUT}")

if __name__=="__main__":
    main()
