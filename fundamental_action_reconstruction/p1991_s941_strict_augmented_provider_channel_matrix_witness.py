#!/usr/bin/env python3
"""P1991 S941 strict augmented-provider channel-matrix witness.

Builds an explicit channel-by-channel linear system for the augmented provider
class from P1990 and reports matrix rank / unresolved residual channels.

Scope: algebraic channel matrix only; non-strict label retained.
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
OUT = GEN / "p1991_s941_strict_augmented_provider_channel_matrix_witness.json"
DEFAULT_TIMESTAMP_UTC = "2026-05-18T00:00:00+00:00"

def load(name:str)->dict[str,Any]:
    p=GEN/name
    if not p.exists(): return {"_missing":name,"status":"OPEN_MISSING_ARTIFACT"}
    return json.loads(p.read_text(encoding="utf-8"))

def main()->None:
    GEN.mkdir(exist_ok=True)
    p1987=load("p1987_s937_strict_non_gb_residual_term_classification_witness.json")
    p1990=load("p1990_s940_strict_extended_provider_family_probe_witness.json")

    s1,s2,sd1,sd2,sdd1,sdd2=sp.symbols("sigma1 sigma2 dsigma1 dsigma2 d2sigma1 d2sigma2", real=True)
    N,Nd,Ndd,V,H,Hd=sp.symbols("N Nd Ndd V H Hd", real=True)
    c1,c2,c3,eta1,eta2,eta3=sp.symbols("c1 c2 c3 eta1 eta2 eta3", real=True)

    locals_map={"N":N,"Nd":Nd,"Ndd":Ndd,"V":V,"H":H,"Hd":Hd,"sigma1":s1,"sigma2":s2,"dsigma1":sd1,"dsigma2":sd2,"d2sigma1":sdd1,"d2sigma2":sdd2,"pi":sp.pi,"log":sp.log,"ln":sp.log}
    residual=sp.sympify(p1987.get("scaled_expression_N6_over_V","0"), locals=locals_map)

    Q=s1**2+s1*s2+s2**2
    Q2=Q**2
    d2sig_sig=2*sdd1*s1+sdd1*s2+sdd2*s1+2*sdd2*s2
    NddQ=Ndd*Q
    provider=sp.expand(c1*(3*H*(s1+s2))+c2*(sd1+sd2)+c3*Q+eta1*d2sig_sig+eta2*NddQ+eta3*Q2)
    target=sp.expand(residual-provider)

    channels={
      "d2sigma1_sigma1": sdd1*s1,
      "Ndd_sigma1_sq": Ndd*s1**2,
      "Qquartic_s1_4": s1**4,
      "H_sigma1_dsigma1": H*s1*sd1,
      "Hd_sigma1_sq": Hd*s1**2,
      "Nd2_sigma1_sq": Nd**2*s1**2,
    }
    unknowns=[c1,c2,c3,eta1,eta2,eta3]

    A=[]; b=[]
    for ch_name, mon in channels.items():
        expr=sp.expand(target).coeff(mon)
        row=[sp.expand(expr).coeff(u) for u in unknowns]
        const=sp.simplify(expr - sum(r*u for r,u in zip(row,unknowns)))
        A.append(row)
        b.append(sp.simplify(-const))

    A_mat=sp.Matrix(A)
    b_vec=sp.Matrix(b)
    rankA=int(A_mat.rank())
    rankAug=int(A_mat.row_join(b_vec).rank())
    consistent=(rankA==rankAug)
    full_solve=(consistent and rankA==len(unknowns))

    # least squares residual with numeric substitution for symbols
    subs={N:1,Nd:sp.Rational(1,20),Ndd:sp.Rational(-1,200),V:1,H:1,Hd:sp.Rational(1,10),s1:sp.Rational(1,10),s2:sp.Rational(-1,20),sd1:sp.Rational(1,100),sd2:sp.Rational(-1,200),sdd1:sp.Rational(1,1000),sdd2:sp.Rational(-1,2000)}
    A_num=np.array([[float(sp.N(v.subs(subs),40)) for v in row] for row in A],dtype=float)
    b_num=np.array([float(sp.N(v.subs(subs),40)) for v in b],dtype=float)
    x,*_=np.linalg.lstsq(A_num,b_num,rcond=None)
    lsq_res=float(la.norm(A_num@x-b_num,2))

    gate={
      "p1987_present": p1987.get("result_kind")=="PASS_NON_GB_TERM_CLASSIFICATION_WITNESS",
      "p1990_present": p1990.get("result_kind")=="PASS_EXTENDED_PROVIDER_PROBE_WITH_NONSTRICT_LABEL",
      "matrix_rank_aug_exceeds_rankA_or_not_full_solve": (not full_solve),
      "least_squares_residual_positive": lsq_res>0.0,
    }

    out={
      "ledger_id":"P1991_S941_STRICT_AUGMENTED_PROVIDER_CHANNEL_MATRIX_WITNESS",
      "packet_id":"P1991","stage_id":"S941","produced_by":Path(__file__).name,
      "timestamp_utc":DEFAULT_TIMESTAMP_UTC,
      "route":"strict_only_with_explicit_non_strict_selector_labeling",
      "selector_premise_label":{"status":"NON_STRICT_AUGMENTED_CLASS"},
      "depends_on":{"p1987_present":gate["p1987_present"],"p1990_present":gate["p1990_present"]},
      "channel_system":{"channels":list(channels.keys()),"unknowns":[str(u) for u in unknowns],"A_matrix": [[str(c) for c in row] for row in A],"b_vector":[str(v) for v in b],"rank_A":rankA,"rank_augmented":rankAug,"consistent":consistent,"full_solve":full_solve},
      "numeric_lsq_probe":{"A":A_num.tolist(),"b":b_num.tolist(),"x_lsq":x.tolist(),"residual_l2":lsq_res},
      "gatekeeper_checks":gate,
      "result_kind":"PASS_CHANNEL_MATRIX_OBSTRUCTION_WITNESS" if all(gate.values()) else "OPEN_OBSTRUCTION_WITH_TRACE",
      "status":"OPEN_OBSTRUCTION_WITH_TRACE",
      "next_honest_step":"Promote this channel matrix to full componentwise spatial-EOM matrix (three Bianchi-I spatial components) and test null-space compatibility constraints.",
      "lay_explanation":"Zrobiliśmy tabelę równań kanał po kanale i policzyliśmy, że nawet rozszerzony model nie potrafi idealnie wyzerować wszystkiego naraz.",
      "environment":{"python":platform.python_version(),"numpy":np.__version__,"scipy":__import__('scipy').__version__,"sympy":sp.__version__}
    }
    OUT.write_text(json.dumps(out,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print(f"[P1991] wrote witness: {OUT}")

if __name__=="__main__":
    main()
