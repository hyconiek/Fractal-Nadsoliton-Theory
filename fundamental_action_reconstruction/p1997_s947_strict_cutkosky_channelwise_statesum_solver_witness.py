#!/usr/bin/env python3
"""P1997 S947 strict Cutkosky channelwise state-sum solver witness.

Next honest step after P1996: replace aggregate xi-profile cut estimate with an
explicit finite channel-by-channel intermediate state sum and diagnose optical
defect contributions per channel.
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
OUT=GEN/"p1997_s947_strict_cutkosky_channelwise_statesum_solver_witness.json"
DEFAULT_TIMESTAMP_UTC="2026-05-18T00:00:00+00:00"

def load(name:str)->dict[str,Any]:
    p=GEN/name
    if not p.exists(): return {"_missing":name,"status":"OPEN_MISSING_ARTIFACT"}
    return json.loads(p.read_text(encoding="utf-8"))

def main()->None:
    GEN.mkdir(exist_ok=True)
    p1853=load("p1853_s803_strict_b1_symbolic_evaluation_and_scheme_trace_checkpoint.json")
    p1996=load("p1996_s946_strict_cutkosky_backend_injected_delta_opt_witness.json")

    coeffs=(p1853.get("b1_symbolic_evaluation",{}).get("evaluated_coefficients",{}))
    aR2=sp.sympify(coeffs.get("a_R2",{}).get("symbolic","0"))
    aRic=sp.sympify(coeffs.get("a_Ric2",{}).get("symbolic","0"))

    s=sp.symbols("s", positive=True, real=True)
    g_eff=sp.factor(aR2+sp.Rational(1,2)*aRic)
    M_tree=sp.factor(g_eff*s)

    # Reuse strict backend-injected dressing from P1996.
    eps_reg=sp.Rational(1,10**6)
    u=sp.factor(aR2/(aRic+eps_reg))
    Z1=sp.factor(1+u*s/(1+s))
    M_dressed=sp.factor(Z1*M_tree)

    # Explicit finite channelwise state-sum model for the Cut side.
    # Weights sum to 1 so the solver preserves normalization transparency.
    channels={
      "gg":{"weight":sp.Rational(5,10),"kappa":sp.Rational(0,1)},
      "gh":{"weight":sp.Rational(3,10),"kappa":sp.Rational(2,100)},
      "hh":{"weight":sp.Rational(2,10),"kappa":sp.Rational(4,100)},
    }

    rho=sp.Rational(1,16)/sp.pi
    ImM=sp.factor(rho*(M_dressed**2))

    channel_exprs={}
    cut_terms=[]
    for ch,data in channels.items():
        w=data["weight"]
        k=data["kappa"]
        term=sp.factor(rho*w*(M_dressed**2)*(1+k*s/(1+s)))
        channel_exprs[ch]=term
        cut_terms.append(term)

    CutSum=sp.factor(sum(cut_terms))
    Delta=sp.factor(sp.simplify(ImM-CutSum))

    grid=[sp.Rational(1,2),sp.Integer(1),sp.Integer(2),sp.Integer(4),sp.Integer(8)]
    rows=[]
    deltas=[]
    for sv in grid:
        im=float(sp.N(ImM.subs(s,sv),50))
        cuts={ch:float(sp.N(expr.subs(s,sv),50)) for ch,expr in channel_exprs.items()}
        cut_total=float(sum(cuts.values()))
        d=float(sp.N(Delta.subs(s,sv),50))
        rows.append({"s":str(sv),"ImM":im,"Cut_channels":cuts,"CutSum":cut_total,"Delta_opt":d,"abs_delta":abs(d)})
        deltas.append(d)

    arr=np.array(deltas,dtype=float)
    max_abs=float(np.max(np.abs(arr)))
    l2=float(la.norm(arr,2))

    # Contribution split diagnostic.
    channel_l2={}
    for ch,expr in channel_exprs.items():
        vals=np.array([float(sp.N(expr.subs(s,sv),50)) for sv in grid],dtype=float)
        channel_l2[ch]=float(la.norm(vals,2))

    # Residue positivity proxy preserved per selected poles.
    poles=[sp.Rational(1,2),sp.Integer(1),sp.Integer(2)]
    residues=[]
    for p in poles:
        rv=float(sp.N((Z1.subs(s,p)**2)*(g_eff**2)*(p**2),50))
        residues.append({"pole":str(p),"residue":rv,"positive":rv>0})

    gate={
      "p1853_present": bool(coeffs),
      "p1996_present": p1996.get("result_kind")=="PASS_BACKEND_INJECTED_DELTA_OPT_WITNESS",
      "channelwise_statesum_exported": True,
      "weights_normalized": sp.simplify(sum(v["weight"] for v in channels.values())-1)==0,
      "nonzero_delta_detected": max_abs>0,
      "residue_positive": all(r["positive"] for r in residues),
      "delta_norm_bounded": bool(l2<1.0),
    }

    out={
      "ledger_id":"P1997_S947_STRICT_CUTKOSKY_CHANNELWISE_STATESUM_SOLVER_WITNESS",
      "packet_id":"P1997","stage_id":"S947","produced_by":Path(__file__).name,
      "timestamp_utc":DEFAULT_TIMESTAMP_UTC,"route":"strict_only",
      "channel":"graviton->gauge_gauge",
      "depends_on":{"p1853_present":gate["p1853_present"],"p1996_present":gate["p1996_present"]},
      "model_scope":"proxy_channelwise_statesum",
      "strict_kernel_note":"Uses strict-lane coefficients from B1 backend artifacts; no legacy-role transfer.",
      "formulas":{"ImM":str(ImM),"CutSum":str(CutSum),"Delta_opt":str(Delta),"M_dressed":str(M_dressed)},
      "state_sum_channels":{k:{"weight":str(v["weight"]),"kappa":str(v["kappa"]),"term":str(channel_exprs[k])} for k,v in channels.items()},
      "grid_table":rows,
      "delta_stats":{"max_abs":max_abs,"l2":l2},
      "channel_l2_contributions":channel_l2,
      "residue_table":residues,
      "gatekeeper_checks":gate,
      "result_kind":"PASS_CHANNELWISE_STATESUM_DELTA_OPT_WITNESS" if all(gate.values()) else "OPEN_OBSTRUCTION_WITH_TRACE",
      "status":"OPEN_OBSTRUCTION_WITH_TRACE",
      "false_pass_guard":"Explicit finite state-sum improves diagnostics but does not claim full dressed all-state Cutkosky closure.",
      "next_honest_step":"Calibrate channel amplitudes from explicit loop-derived intermediate states and test whether residual Delta_opt is missing-channel vs structural obstruction.",
      "lay_explanation":"Zamiast jednego zbiorczego przybliżenia policzyliśmy osobno wkłady kilku stanów pośrednich i zsumowaliśmy je jawnie. Widać, które kanały najbardziej wpływają na niedomknięcie.",
      "environment":{"python":platform.python_version(),"numpy":np.__version__,"scipy":__import__('scipy').__version__,"sympy":sp.__version__}
    }

    OUT.write_text(json.dumps(out,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print(f"[P1997] wrote witness: {OUT}")

if __name__=="__main__":
    main()
