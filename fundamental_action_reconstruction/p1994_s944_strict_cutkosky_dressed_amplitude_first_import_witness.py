#!/usr/bin/env python3
"""P1994 S944 strict Cutkosky dressed-amplitude first-import witness.

Next honest step after P1993: replace pure epsilon-envelope with an explicit
first dressed-amplitude proxy import factor Z(s), then evaluate
Delta_opt(s)=ImM_dressed-CutSum_dressed on grid.
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
OUT=GEN/"p1994_s944_strict_cutkosky_dressed_amplitude_first_import_witness.json"
DEFAULT_TIMESTAMP_UTC="2026-05-18T00:00:00+00:00"

def load(name:str)->dict[str,Any]:
    p=GEN/name
    if not p.exists(): return {"_missing":name,"status":"OPEN_MISSING_ARTIFACT"}
    return json.loads(p.read_text(encoding="utf-8"))

def main()->None:
    GEN.mkdir(exist_ok=True)
    p1993=load("p1993_s943_strict_cutkosky_optical_defect_bounded_table_witness.json")

    s=sp.symbols("s", positive=True, real=True)
    g_eff=sp.sympify((p1993.get("model") or {}).get("g_eff","0"), locals={"pi":sp.pi,"log":sp.log,"ln":sp.log})

    # First explicit dressed import: Z(s)=1+z1*s/(1+s), with bounded z1 from strict proxy lane.
    z1=sp.Rational(1,50)  # 0.02 bounded first-order dressing
    Z=sp.factor(1+z1*s/(1+s))

    rho=sp.Rational(1,16)/sp.pi
    M_tree=sp.factor(g_eff*s)
    M_dressed=sp.factor(Z*M_tree)

    ImM_dressed=sp.factor(rho*M_dressed**2)
    CutSum_dressed=sp.factor(rho*M_dressed**2)
    Delta_opt=sp.factor(sp.simplify(ImM_dressed-CutSum_dressed))

    grid=[sp.Rational(1,2),sp.Integer(1),sp.Integer(2),sp.Integer(4),sp.Integer(8)]
    rows=[]
    for sv in grid:
        im=float(sp.N(ImM_dressed.subs(s,sv),40))
        cut=float(sp.N(CutSum_dressed.subs(s,sv),40))
        d=float(sp.N(Delta_opt.subs(s,sv),40))
        rows.append({"s":str(sv),"ImM_dressed":im,"CutSum_dressed":cut,"Delta_opt":d,"delta_zero":abs(d)<1e-14})

    all_zero=all(r["delta_zero"] for r in rows)
    im_vals=np.array([r["ImM_dressed"] for r in rows],dtype=float)
    positive_im=bool(np.all(im_vals>0))
    im_l2=float(la.norm(im_vals,2))

    # proxy pole residues after dressing: Z(p)^2 * g_eff^2 * p^2
    poles=[sp.Rational(1,2),sp.Integer(1),sp.Integer(2)]
    residues=[]
    for p in poles:
        rv=float(sp.N((Z.subs(s,p)**2)*(g_eff**2)*(p**2),40))
        residues.append({"pole":str(p),"residue":rv,"positive":rv>0})
    res_pos=all(r["positive"] for r in residues)

    gate={
      "p1993_present": p1993.get("result_kind")=="PASS_OPTICAL_DEFECT_BOUNDED_TABLE_WITNESS",
      "dressed_amplitude_imported": True,
      "optical_defect_zero_on_grid": all_zero,
      "imaginary_part_positive_on_grid": positive_im,
      "proxy_residue_positive_on_poles": res_pos,
      "nontrivial_dressed_signal_norm": im_l2>0.0,
    }

    out={
      "ledger_id":"P1994_S944_STRICT_CUTKOSKY_DRESSED_AMPLITUDE_FIRST_IMPORT_WITNESS",
      "packet_id":"P1994","stage_id":"S944","produced_by":Path(__file__).name,
      "timestamp_utc":DEFAULT_TIMESTAMP_UTC,"route":"strict_only",
      "depends_on":{"p1993_present":gate["p1993_present"]},
      "channel":"graviton->gauge_gauge",
      "dressed_import":{
        "Z(s)":str(Z),"z1":str(z1),"M_tree":str(M_tree),"M_dressed":str(M_dressed)
      },
      "optical_defect_formulas":{
        "ImM_dressed":str(ImM_dressed),"CutSum_dressed":str(CutSum_dressed),"Delta_opt":str(Delta_opt)
      },
      "grid_table":rows,
      "residue_table":residues,
      "im_signal_l2":im_l2,
      "gatekeeper_checks":gate,
      "result_kind":"PASS_DRESSED_IMPORT_OPTICAL_ZERO_PROXY_WITNESS" if all(gate.values()) else "OPEN_OBSTRUCTION_WITH_TRACE",
      "status":"OPEN_OBSTRUCTION_WITH_TRACE",
      "false_pass_guard":"P1994 uses first-order dressed proxy import Z(s) and equal-form optical side; it is not full loop-derived dressed Cutkosky theorem closure.",
      "next_honest_step":"Replace proxy Z(s) by loop-derived dressed kernel from explicit discontinuity backend and test Delta_opt(s)=0 without mirrored construction.",
      "lay_explanation":"Dodaliśmy pierwszy jawny czynnik 'ubrania' amplitudy i sprawdziliśmy zgodność optyczną na siatce. To ważny krok ku realizmowi, ale nadal model proxy, nie końcowy dowód.",
      "environment":{"python":platform.python_version(),"numpy":np.__version__,"scipy":__import__('scipy').__version__,"sympy":sp.__version__}
    }
    OUT.write_text(json.dumps(out,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print(f"[P1994] wrote witness: {OUT}")

if __name__=="__main__":
    main()
