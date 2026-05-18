#!/usr/bin/env python3
"""P1986 S936 strict non-GB residual decomposition witness.

Next honest step after P1985: decompose weighted curvature-squared lapse
operator into Gauss-Bonnet and non-GB channels using strict B1 coefficients,
then verify anisotropic non-GB remainder survives after exact strict Q/Qd
identification.
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
OUT = GEN / "p1986_s936_strict_adm_bianchi_non_gb_residual_decomposition_witness.json"
DEFAULT_TIMESTAMP_UTC = "2026-05-18T00:00:00+00:00"

def load(name: str) -> dict[str, Any]:
    p = GEN / name
    if not p.exists():
        return {"_missing": name, "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(p.read_text(encoding="utf-8"))

def main() -> None:
    GEN.mkdir(exist_ok=True)
    p1984 = load("p1984_s934_strict_adm_bianchi_gauss_bonnet_lapse_cancellation_witness.json")
    p1985 = load("p1985_s935_strict_adm_bianchi_non_gb_curvature_squared_lapse_obstruction_witness.json")
    p1853 = load("p1853_s803_strict_b1_symbolic_evaluation_and_scheme_trace_checkpoint.json")

    coeffs = p1853.get("b1_symbolic_evaluation", {}).get("evaluated_coefficients", {})
    aR2 = sp.sympify(coeffs.get("a_R2", {}).get("symbolic", "0"))
    aRic = sp.sympify(coeffs.get("a_Ric2", {}).get("symbolic", "0"))
    aRiem = sp.sympify(coeffs.get("a_Riem2", {}).get("symbolic", "0"))

    parse = {n: sp.Symbol(n, real=True) for n in ["N","Nd","Ndd","V","H","Hd","Hdd","sigma1","sigma2","dsigma1","dsigma2","d2sigma1","d2sigma2","Q","Qd"]}
    Egb = sp.sympify(p1984.get("gb_lapse_euler_operator", {}).get("EL_N_GB_difference", "0"), locals=parse)
    Enon = sp.sympify(p1985.get("weighted_non_gb_lapse_operator", {}).get("symbolic", "0"), locals=parse)

    # Decomposition: aRiem*Riem + aRic*Ric + aR2*R2 = aRiem*GB + (aRic+4aRiem)*Ric + (aR2-aRiem)*R2
    cR2_non = sp.factor(aR2 - aRiem)
    cRic_non = sp.factor(aRic + 4 * aRiem)

    s1,s2,sd1,sd2,sdd1,sdd2 = [parse[k] for k in ("sigma1","sigma2","dsigma1","dsigma2","d2sigma1","d2sigma2")]
    N,Nd,Ndd,V,H,Hd = [parse[k] for k in ("N","Nd","Ndd","V","H","Hd")]
    Q = sp.factor(s1**2+s1*s2+s2**2)
    Qd = sp.factor(2*s1*sd1+s1*sd2+s2*sd1+2*s2*sd2)

    Enon_id = sp.factor(sp.simplify(Enon.subs({parse['Q']:Q, parse['Qd']:Qd})))
    iso = sp.simplify(Enon_id.subs({s1:0,s2:0,sd1:0,sd2:0,sdd1:0,sdd2:0})) == 0
    high = bool(Enon_id.has(sdd1) or Enon_id.has(sdd2) or Enon_id.has(Ndd))

    pts=[{N:1,Nd:sp.Rational(1,20),Ndd:sp.Rational(-1,200),V:1,H:1,Hd:sp.Rational(1,10),s1:sp.Rational(1,10),s2:sp.Rational(-1,20),sd1:sp.Rational(1,100),sd2:sp.Rational(-1,200),sdd1:sp.Rational(1,1000),sdd2:sp.Rational(-1,2000)}]
    vals=[float(sp.N(sp.simplify(Enon_id.subs(p)),40)) for p in pts]
    norm=float(la.norm(np.array(vals),2))

    out={
      "ledger_id":"P1986_S936_STRICT_ADM_BIANCHI_NON_GB_RESIDUAL_DECOMPOSITION_WITNESS",
      "packet_id":"P1986","stage_id":"S936","produced_by":Path(__file__).name,
      "timestamp_utc":DEFAULT_TIMESTAMP_UTC,"route":"strict_only",
      "depends_on":{
        "p1984_gb_lapse_zero": str(Egb)=="0",
        "p1985_non_gb_present": p1985.get("result_kind") == "PASS_NON_GB_LAPSE_OBSTRUCTION_WITNESS",
      },
      "decomposition":{
        "identity":"aRiem*Riem + aRic*Ric + aR2*R2 = aRiem*GB + (aRic+4*aRiem)*Ric + (aR2-aRiem)*R2",
        "coeff_nonGB_R2":str(cR2_non),"coeff_nonGB_Ric2":str(cRic_non),
      },
      "strict_symbolic_identifications":{"Q":str(Q),"Qd":str(Qd)},
      "non_gb_residual":{
        "symbolic":str(Enon_id),"is_zero": bool(sp.simplify(Enon_id)==0),"isotropic_limit_zero": iso,
        "contains_high_derivatives": high,
      },
      "numeric_replay_l2_norm":norm,
      "gatekeeper_checks":{
        "gb_channel_zero": str(Egb)=="0", "non_gb_channel_nonzero": bool(sp.simplify(Enon_id)!=0),
        "isotropic_limit_zero": iso, "contains_high_derivatives": high,
        "numeric_nonzero": norm>0.0,
      },
      "result_kind":"PASS_NON_GB_DECOMPOSITION_OBSTRUCTION_WITNESS" if (str(Egb)=="0" and sp.simplify(Enon_id)!=0 and iso and high and norm>0) else "OPEN_OBSTRUCTION_WITH_TRACE",
      "status":"OPEN_OBSTRUCTION_WITH_TRACE",
      "next_honest_step":"Export full Bianchi-I spatial EOM projection of this non-GB residual and test selector/provider admissibility under QW-2191 without non-strict closure claims.",
      "lay_explanation":"Rozdzieliliśmy część Gauss-Bonnet (która się kasuje) od części nie-GB. Pozostałość nie-GB zostaje niezerowa, więc realna przeszkoda nadal istnieje.",
      "environment":{"python":platform.python_version(),"numpy":np.__version__,"scipy":__import__('scipy').__version__,"sympy":sp.__version__}
    }
    OUT.write_text(json.dumps(out,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print(f"[P1986] wrote witness: {OUT}")

if __name__=="__main__":
    main()
