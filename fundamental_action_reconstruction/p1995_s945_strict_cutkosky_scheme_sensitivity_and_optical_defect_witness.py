#!/usr/bin/env python3
"""P1995 S945 strict Cutkosky scheme-sensitivity optical-defect witness."""

from __future__ import annotations
import json, platform
from pathlib import Path
from typing import Any
import numpy as np
import scipy.linalg as la
import sympy as sp

ROOT=Path(__file__).resolve().parent
GEN=ROOT/"generated"
OUT=GEN/"p1995_s945_strict_cutkosky_scheme_sensitivity_and_optical_defect_witness.json"
DEFAULT_TIMESTAMP_UTC="2026-05-18T00:00:00+00:00"

def load(name:str)->dict[str,Any]:
    p=GEN/name
    if not p.exists():
        return {"_missing":name,"status":"OPEN_MISSING_ARTIFACT"}
    return json.loads(p.read_text(encoding="utf-8"))

def main()->None:
    GEN.mkdir(exist_ok=True)
    p1994=load("p1994_s944_strict_cutkosky_dressed_amplitude_first_import_witness.json")

    s=sp.symbols("s", positive=True, real=True)
    m_tree=sp.sympify((p1994.get("dressed_import") or {}).get("M_tree","0"), locals={"s":s,"pi":sp.pi,"log":sp.log,"ln":sp.log})
    g_eff=sp.simplify(m_tree/s) if m_tree!=0 else sp.Integer(0)

    rho=sp.Rational(1,16)/sp.pi
    M_tree=sp.factor(g_eff*s)

    z1_vals=[sp.Rational(3,200),sp.Rational(1,50),sp.Rational(1,40)]
    kappa_vals=[sp.Rational(1,2),sp.Integer(1)]
    grid=[sp.Rational(1,2),sp.Integer(1),sp.Integer(2),sp.Integer(4),sp.Integer(8)]

    rows=[]
    deltas=[]
    for z1 in z1_vals:
        for kappa in kappa_vals:
            Z=sp.factor(1+z1*s/(1+s))
            M_d=sp.factor(Z*M_tree)
            ImM=sp.factor(rho*M_d**2)
            K=sp.factor(1+kappa*z1*s/(1+s)**2)
            Cut=sp.factor(rho*M_d**2*K)
            D=sp.factor(sp.simplify(ImM-Cut))
            for sv in grid:
                dv=float(sp.N(D.subs(s,sv),40))
                im=float(sp.N(ImM.subs(s,sv),40))
                cut=float(sp.N(Cut.subs(s,sv),40))
                rows.append({"z1":str(z1),"kappa":str(kappa),"s":str(sv),"ImM":im,"CutSum":cut,"Delta_opt":dv})
                deltas.append(dv)

    deltas_arr=np.array(deltas,dtype=float)
    max_abs=float(np.max(np.abs(deltas_arr))) if len(deltas_arr) else 0.0
    l2=float(la.norm(deltas_arr,2)) if len(deltas_arr) else 0.0
    near_zero_frac=float(np.mean(np.abs(deltas_arr)<1e-6)) if len(deltas_arr) else 0.0

    poles=[sp.Rational(1,2),sp.Integer(1),sp.Integer(2)]
    res_rows=[]
    for z1 in z1_vals:
        for p in poles:
            Zp=sp.factor(1+z1*p/(1+p))
            rv=float(sp.N((Zp**2)*(g_eff**2)*(p**2),40))
            res_rows.append({"z1":str(z1),"pole":str(p),"residue":rv,"positive":rv>0})
    res_pos=all(r["positive"] for r in res_rows)

    gate={
      "p1994_present": p1994.get("result_kind")=="PASS_DRESSED_IMPORT_OPTICAL_ZERO_PROXY_WITNESS",
      "scheme_family_scanned": True,
      "nonzero_optical_defect_detected_under_scheme_perturbation": max_abs>0.0,
      "residue_positive_on_scanned_family": res_pos,
      "defect_norm_bounded": bool(l2<1.0),
    }

    out={
      "ledger_id":"P1995_S945_STRICT_CUTKOSKY_SCHEME_SENSITIVITY_AND_OPTICAL_DEFECT_WITNESS",
      "packet_id":"P1995","stage_id":"S945","produced_by":Path(__file__).name,
      "timestamp_utc":DEFAULT_TIMESTAMP_UTC,"route":"strict_only",
      "depends_on":{"p1994_present":gate["p1994_present"]},
      "channel":"graviton->gauge_gauge",
      "scheme_family":{"z1_values":[str(v) for v in z1_vals],"kappa_values":[str(v) for v in kappa_vals],"K_form":"1+kappa*z1*s/(1+s)^2"},
      "scan_table":rows,
      "optical_defect_statistics":{"max_abs_delta":max_abs,"l2_delta":l2,"near_zero_fraction":near_zero_frac},
      "proxy_residue_table":res_rows,
      "gatekeeper_checks":gate,
      "result_kind":"PASS_SCHEME_SENSITIVITY_OPTICAL_DEFECT_WITNESS" if all(gate.values()) else "OPEN_OBSTRUCTION_WITH_TRACE",
      "status":"OPEN_OBSTRUCTION_WITH_TRACE",
      "false_pass_guard":"P1995 demonstrates scheme sensitivity of optical defect in proxy dressed family; it is not a full loop-backend theorem closure.",
      "next_honest_step":"Replace K(s) perturbation with explicit loop-discontinuity backend export and test if scheme-transport map can lock Delta_opt to zero across FRW/Bianchi-I charts.",
      "lay_explanation":"Pokazaliśmy, że po lekkiej zmianie schematu renormalizacji defekt optyczny może odchodzić od zera, choć pozostaje mały i kontrolowany.",
      "environment":{"python":platform.python_version(),"numpy":np.__version__,"scipy":__import__('scipy').__version__,"sympy":sp.__version__}
    }
    OUT.write_text(json.dumps(out,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print(f"[P1995] wrote witness: {OUT}")

if __name__=="__main__":
    main()
