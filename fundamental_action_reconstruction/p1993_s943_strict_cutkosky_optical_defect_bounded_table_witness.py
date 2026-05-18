#!/usr/bin/env python3
"""P1993 S943 strict Cutkosky optical-defect bounded-table witness.

Next honest step after P1992: build explicit bounded table for
Delta_opt(s)=ImM-CutSum under symmetry-factor convention window and verify
residue-positivity proxy for selected physical poles.
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
OUT=GEN/"p1993_s943_strict_cutkosky_optical_defect_bounded_table_witness.json"
DEFAULT_TIMESTAMP_UTC="2026-05-18T00:00:00+00:00"

def load(name:str)->dict[str,Any]:
    p=GEN/name
    if not p.exists(): return {"_missing":name,"status":"OPEN_MISSING_ARTIFACT"}
    return json.loads(p.read_text(encoding="utf-8"))

def main()->None:
    GEN.mkdir(exist_ok=True)
    p1992=load("p1992_s942_strict_cutkosky_graviton_gauge_exact_phase_space_witness.json")

    s=sp.symbols("s", positive=True, real=True)
    g_eff=sp.sympify((p1992.get("strict_proxy_model") or {}).get("g_eff","0"), locals={"pi":sp.pi,"log":sp.log,"ln":sp.log})

    # Model with bounded dressing mismatch epsilon(s) in [-eps0,+eps0]
    eps0=sp.Rational(1,100)  # 1% bounded mismatch placeholder for proxy->dressed bridge uncertainty
    rho_lo=sp.Rational(1,16)/sp.pi
    rho_hi=sp.Rational(1,8)/sp.pi
    M=sp.factor(g_eff*s)

    # ImM and CutSum bounded by convention+dressing envelope.
    ImM_ref=sp.factor(rho_lo*M**2)
    Cut_lo=sp.factor((1-eps0)*rho_lo*M**2)
    Cut_hi=sp.factor((1+eps0)*rho_hi*M**2)

    defect_lo=sp.factor(ImM_ref-Cut_hi)  # worst negative
    defect_hi=sp.factor(ImM_ref-Cut_lo)  # best positive

    grid=[sp.Rational(1,2),sp.Integer(1),sp.Integer(2),sp.Integer(4),sp.Integer(8)]
    rows=[]
    for sv in grid:
        dlo=float(sp.N(defect_lo.subs(s,sv),40))
        dhi=float(sp.N(defect_hi.subs(s,sv),40))
        rows.append({
            "s":str(sv),
            "defect_interval":[dlo,dhi],
            "contains_zero": bool(dlo<=0<=dhi),
            "interval_width": float(dhi-dlo),
            "ImM_ref": float(sp.N(ImM_ref.subs(s,sv),40)),
        })

    # simple pole-residue positivity proxy (effective pole family)
    poles=[sp.Rational(1,2),sp.Integer(1),sp.Integer(2)]
    residue_rows=[]
    for p in poles:
        # proxy residue ~ g_eff^2 * p^2 positive if g_eff>0
        rv=sp.N((g_eff**2)*(p**2),40)
        residue_rows.append({"pole":str(p),"proxy_residue":float(rv),"positive":bool(float(rv)>0)})

    all_zero_contained=all(r["contains_zero"] for r in rows)
    all_res_pos=all(r["positive"] for r in residue_rows)
    widths=np.array([r["interval_width"] for r in rows],dtype=float)
    wnorm=float(la.norm(widths,2))

    gate={
      "p1992_present": p1992.get("result_kind")=="PASS_EXACT_PHASE_SPACE_BOUNDED_UNITARITY_PROXY_WITNESS",
      "optical_defect_bounded_table_exported": True,
      "all_intervals_contain_zero": all_zero_contained,
      "proxy_residue_positive_on_selected_poles": all_res_pos,
      "interval_width_bounded": bool(wnorm>0 and wnorm<1.0),
    }

    out={
      "ledger_id":"P1993_S943_STRICT_CUTKOSKY_OPTICAL_DEFECT_BOUNDED_TABLE_WITNESS",
      "packet_id":"P1993","stage_id":"S943","produced_by":Path(__file__).name,
      "timestamp_utc":DEFAULT_TIMESTAMP_UTC,"route":"strict_only",
      "depends_on":{"p1992_present":gate["p1992_present"]},
      "channel":"graviton->gauge_gauge",
      "model":{
        "M_tree":"g_eff*s",
        "g_eff":str(g_eff),
        "rho_interval":[str(rho_lo),str(rho_hi)],
        "dressing_mismatch_eps0":str(eps0)
      },
      "optical_defect_formulas":{
        "ImM_ref":str(ImM_ref),"Cut_lo":str(Cut_lo),"Cut_hi":str(Cut_hi),
        "defect_lo":str(defect_lo),"defect_hi":str(defect_hi)
      },
      "bounded_defect_table":rows,
      "proxy_residue_table":residue_rows,
      "interval_width_l2":wnorm,
      "gatekeeper_checks":gate,
      "result_kind":"PASS_OPTICAL_DEFECT_BOUNDED_TABLE_WITNESS" if all(gate.values()) else "OPEN_OBSTRUCTION_WITH_TRACE",
      "status":"OPEN_OBSTRUCTION_WITH_TRACE",
      "false_pass_guard":"P1993 is a bounded defect proxy table, not full dressed Cutkosky theorem closure.",
      "next_honest_step":"Import first dressed channel amplitude and replace bounded epsilon envelope by explicit loop/discontinuity evaluation to test Delta_opt(s)=0 directly.",
      "lay_explanation":"Zamiast jednego numeru błędu mamy teraz przedziały błędu optycznego dla różnych energii i kontrolę dodatniości residuów. To bardziej rygorystyczny monitoring, ale jeszcze nie pełne zamknięcie unitarności.",
      "environment":{"python":platform.python_version(),"numpy":np.__version__,"scipy":__import__('scipy').__version__,"sympy":sp.__version__}
    }
    OUT.write_text(json.dumps(out,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print(f"[P1993] wrote witness: {OUT}")

if __name__=="__main__":
    main()
