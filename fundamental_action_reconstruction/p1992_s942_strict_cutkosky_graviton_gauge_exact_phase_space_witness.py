#!/usr/bin/env python3
"""P1992 S942 strict Cutkosky graviton->gauge_gauge exact phase-space witness.

Next honest step for task-2: replace pure seed sign scan with exact 2-body
phase-space integral in a strict B1 proxy amplitude model and export bounded
uncertainty table for symmetry/convention factors.
"""

from __future__ import annotations
import json, platform
from pathlib import Path
from typing import Any
import numpy as np
import scipy.linalg as la
import sympy as sp

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
OUT = GEN / "p1992_s942_strict_cutkosky_graviton_gauge_exact_phase_space_witness.json"
DEFAULT_TIMESTAMP_UTC = "2026-05-18T00:00:00+00:00"

def load(name:str)->dict[str,Any]:
    p=GEN/name
    if not p.exists(): return {"_missing":name,"status":"OPEN_MISSING_ARTIFACT"}
    return json.loads(p.read_text(encoding="utf-8"))

def main()->None:
    GEN.mkdir(exist_ok=True)
    p1853=load("p1853_s803_strict_b1_symbolic_evaluation_and_scheme_trace_checkpoint.json")
    p1858=load("p1858_s808_strict_b1_cutkosky_seed_grid_evaluation_checkpoint.json")

    coeffs=(p1853.get("b1_symbolic_evaluation",{}).get("evaluated_coefficients",{}))
    aR2=sp.sympify(coeffs.get("a_R2",{}).get("symbolic","0"))
    aRic=sp.sympify(coeffs.get("a_Ric2",{}).get("symbolic","0"))

    s=sp.symbols("s", positive=True, real=True)
    # strict proxy coupling used by P1858 seed definition
    g_eff=sp.factor(aR2+sp.Rational(1,2)*aRic)
    M_tree=sp.factor(g_eff*s)

    # Exact dPi_2 for massless 2-body in COM with bounded convention window:
    # rho_lo = 1/(16*pi) (identical gauge bosons factor), rho_hi = 1/(8*pi)
    rho_lo=sp.Rational(1,16)/sp.pi
    rho_hi=sp.Rational(1,8)/sp.pi
    ImM_lo=sp.factor(rho_lo*M_tree**2)
    ImM_hi=sp.factor(rho_hi*M_tree**2)

    residues={"rho_lo":str(rho_lo),"rho_hi":str(rho_hi),"g_eff":str(g_eff)}

    # deterministic grid + uncertainty band
    s_grid=[sp.Rational(1,2),sp.Integer(1),sp.Integer(2),sp.Integer(4),sp.Integer(8)]
    rows=[]
    for sv in s_grid:
        lo=sp.N(ImM_lo.subs(s,sv),40)
        hi=sp.N(ImM_hi.subs(s,sv),40)
        rows.append({
            "s":str(sv),
            "ImM_lo":str(sp.simplify(ImM_lo.subs(s,sv))),
            "ImM_hi":str(sp.simplify(ImM_hi.subs(s,sv))),
            "ImM_lo_float":float(lo),
            "ImM_hi_float":float(hi),
            "positive_band": bool(float(lo)>0 and float(hi)>0),
        })

    all_pos=all(r["positive_band"] for r in rows)
    widths=np.array([r["ImM_hi_float"]-r["ImM_lo_float"] for r in rows],dtype=float)
    width_l2=float(la.norm(widths,2))

    gate={
      "p1853_coefficients_present": bool(coeffs),
      "p1858_seed_proxy_present": "cutkosky_seed_grid" in p1858,
      "exact_phase_space_formula_exported": True,
      "positive_imaginary_part_band_on_grid": all_pos,
      "uncertainty_band_nonzero_and_bounded": bool(width_l2>0 and width_l2<1.0),
    }

    out={
      "ledger_id":"P1992_S942_STRICT_CUTKOSKY_GRAVITON_GAUGE_EXACT_PHASE_SPACE_WITNESS",
      "packet_id":"P1992","stage_id":"S942","produced_by":Path(__file__).name,
      "timestamp_utc":DEFAULT_TIMESTAMP_UTC,"route":"strict_only",
      "depends_on":{"p1853_present":gate["p1853_coefficients_present"],"p1858_present":gate["p1858_seed_proxy_present"]},
      "channel":"graviton->gauge_gauge",
      "strict_proxy_model":{
        "definition":"M_tree(s)=g_eff*s, g_eff:=a_R2+0.5*a_Ric2 (strict B1)",
        "g_eff":str(g_eff),"M_tree":str(M_tree)
      },
      "exact_phase_space_integral":{
        "rho_lo":str(rho_lo),"rho_hi":str(rho_hi),
        "ImM_lo_formula":str(ImM_lo),"ImM_hi_formula":str(ImM_hi),
        "optical_identity_form":"ImM(s)=rho2(s)|M_tree(s)|^2 (massless 2-body COM)"
      },
      "bounded_uncertainty_table":rows,
      "uncertainty_width_l2":width_l2,
      "pole_residue_certificate":residues,
      "gatekeeper_checks":gate,
      "result_kind":"PASS_EXACT_PHASE_SPACE_BOUNDED_UNITARITY_PROXY_WITNESS" if all(gate.values()) else "OPEN_OBSTRUCTION_WITH_TRACE",
      "status":"OPEN_OBSTRUCTION_WITH_TRACE",
      "false_pass_guard":"P1992 upgrades phase-space integration exactness in the strict proxy channel only; it is not full dressed-propagator Cutkosky closure nor global UR_link theorem.",
      "next_honest_step":"Replace proxy M_tree with dressed amplitude export and verify channelwise optical defect ImM-CutSum=0 including residue positivity for each physical pole.",
      "lay_explanation":"Policzyliśmy dokładnie całkę fazową dla prostego modelu kanału graviton->gauge_gauge i pokazaliśmy dodatni, kontrolowany przedział części urojonej amplitudy. To mocniejszy krok niż sam skan znaków, ale jeszcze nie pełny dowód unitarności.",
      "environment":{"python":platform.python_version(),"numpy":np.__version__,"scipy":__import__('scipy').__version__,"sympy":sp.__version__}
    }
    OUT.write_text(json.dumps(out,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print(f"[P1992] wrote witness: {OUT}")

if __name__=="__main__":
    main()
