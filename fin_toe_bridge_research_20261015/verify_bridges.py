#!/usr/bin/env python3
"""Small new bridge tests; exact rational checking after interval input binding.

Not a proof of physical locality, quantum mechanics, or a native operator limit.
"""
import argparse
from copy import deepcopy
from fractions import Fraction as F
import hashlib
import importlib.util
import json
from pathlib import Path
import math
import numpy as np
import mpmath as mp
from scipy.linalg import expm
from scipy.special import softmax
import sympy as sy

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
SPECTRA=ROOT/"fin_rank7_followup/inputs/fin_handoff_audit/results.json"
SEED=ROOT/"fin_physical_bridge_campaign_01/PHYS-002/safe_seed.py"

class I:
    def __init__(self,a,b=None):
        self.lo=F(a); self.hi=F(a if b is None else b)
        assert self.lo<=self.hi
    def __add__(self,other):
        other=other if isinstance(other,I) else I(other)
        return I(self.lo+other.lo,self.hi+other.hi)
    __radd__=__add__
    def __neg__(self):return I(-self.hi,-self.lo)
    def __sub__(self,other):return self+-other if isinstance(other,I) else self+I(-F(other))
    def __rsub__(self,other):return I(other)+-self
    def __mul__(self,other):
        other=other if isinstance(other,I) else I(other)
        vals=[a*b for a in (self.lo,self.hi) for b in (other.lo,other.hi)]
        return I(min(vals),max(vals))
    __rmul__=__mul__
    def __truediv__(self,other):
        other=other if isinstance(other,I) else I(other)
        assert not other.lo<=0<=other.hi
        return self*I(1/other.hi,1/other.lo)
    def public(self,digits=18):
        scale=10**digits
        a=self.lo.numerator*scale//self.lo.denominator
        b=-((-self.hi.numerator*scale)//self.hi.denominator)
        def dec(n):
            prefix="-" if n<0 else ""; n=abs(n)
            return prefix+str(n//scale)+"."+str(n%scale).zfill(digits)
        result=[dec(a),dec(b)]
        assert F(result[0])<=self.lo<=self.hi<=F(result[1])
        return result

def dyadic(t):
    sign,mantissa,exponent,_=t
    value=F(int(mantissa))*F(2)**int(exponent)
    return -value if sign else value

def from_iv(x):return I(*(dyadic(t) for t in x._mpi_))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def build_certificate():
    mp.iv.dps=80; iv=mp.iv
    raw=json.loads(SPECTRA.read_text())["exact"]["laplacian_intervals"]
    lam={k:I(*raw[k]) for k in (3,4,5,6)}
    # Rational enclosure, independently checked by squaring.
    sqrt3=I("1.7320508075688772935274463415058723669428052538103806280558069794",
            "1.7320508075688772935274463415058723669428052538103806280558069795")
    assert sqrt3.lo**2<3<sqrt3.hi**2
    slope_left=(lam[5]-lam[4])/(sqrt3-1)
    slope_right=(lam[6]-lam[5])/(2-sqrt3)
    defect=slope_right-slope_left
    assert defect.lo>0
    weights={d:iv.cos(iv.mpf(".18575")*d+iv.mpf(".1625"))/(1+iv.mpf(d)**(iv.mpf(9)/5)) for d in range(1,7)}
    binding={}
    for k in (3,4,5,6):
        val=2*sum(weights[d]*(1-iv.cos(iv.pi*k*d/6)) for d in range(1,6))+weights[6]*(1-iv.cos(iv.pi*k))
        exact=from_iv(val)
        assert lam[k].lo<=exact.lo<=exact.hi<=lam[k].hi
        binding[str(k)]=exact.public(40)
    public_w={str(d):from_iv(weights[d]).public(40) for d in range(1,7)}
    # Reconstruct the public-input certificate in exact rational arithmetic.
    w={d:I(*public_w[str(d)]) for d in range(1,7)}
    minor=w[6]*w[3]-w[4]*w[5]
    assert minor.hi<0
    result={
        "date":"2026-10-15", "scope":"two finite class-specific locality obstructions, not all-locality no-go",
        "spectral_path":str(SPECTRA.relative_to(ROOT)),"spectral_sha256":sha(SPECTRA),
        "arithmetic":"exact Fraction intervals after inherited rational spectra and fresh mpmath.iv exact-decimal kernel binding",
        "model":"q12 supplied positive full Laplacian L12, distinct from projected A7 and strict infinite absolute-kernel model",
        "spectral_inputs":{str(k):[str(lam[k].lo),str(lam[k].hi)] for k in lam},
        "fresh_kernel_spectral_binding":binding,
        "sqrt3":[str(sqrt3.lo),str(sqrt3.hi)],
        "bernstein_obstruction":{
            "base":"nearest-neighbour C12 Laplacian; mu4=3, mu5=2+sqrt3, mu6=4",
            "left_slope":slope_left.public(),"right_slope":slope_right.public(),
            "concavity_defect":defect.public(),
            "conclusion":"L12 is not f(L_cycle) for any Bernstein f; positive scalar clock rescaling and added constant mass do not repair concavity",
            "rank7_separate_exact_obstruction":"A7 has equal zero values at distinct positive base modes1,2 and positive values later; no nonconstant Bernstein f can have that plateau"
        },
        "kernel_weight_inputs":public_w,
        "circular_planar_obstruction":{
            "response_convention":"M=-L12, off-diagonal M_ij=w_d",
            "P":[0,1],"Q":[6,4],"circular_order":[0,1,4,6],
            "determinant":"w6*w3-w4*w5","enclosure":minor.public(),
            "conclusion":"No finite positive circular-planar resistor realization with all ports on one disk boundary in the natural cyclic order",
            "not_excluded":["other port orders", "nonplanar or multiple-boundary graphs", "matrix-valued or constrained bulk", "general emergent locality"]
        },
        "not_proved":["physical origin of A7 or L12", "unique bulk", "physical spatial dimension", "operator-limit theorem", "QM/GR/SM/ToE" ]
    }
    check_certificate(result)
    return result

def check_certificate(c):
    """No numpy, mpmath or producer transcendental reevaluation is used here."""
    l={int(k):I(*v) for k,v in c["spectral_inputs"].items()}
    r=I(*c["sqrt3"]); assert r.lo*r.lo<3<r.hi*r.hi
    d=(l[6]-l[5])/(2-r)-(l[5]-l[4])/(r-1)
    displayed=I(*c["bernstein_obstruction"]["concavity_defect"])
    assert 0<displayed.lo<=d.lo<=d.hi<=displayed.hi
    w={int(k):I(*v) for k,v in c["kernel_weight_inputs"].items()}
    m=w[6]*w[3]-w[4]*w[5]
    shown=I(*c["circular_planar_obstruction"]["enclosure"])
    assert shown.lo<=m.lo<=m.hi<=shown.hi<0
    assert c["circular_planar_obstruction"]["P"]==[0,1] and c["circular_planar_obstruction"]["Q"]==[6,4]
    assert c["circular_planar_obstruction"]["circular_order"]==[0,1,4,6]

def main():
    parser=argparse.ArgumentParser(); parser.add_argument("--output",type=Path); parser.add_argument("--check",type=Path)
    args=parser.parse_args(); checks=[]
    def ck(name,condition,**detail):
        checks.append({"name":name,"passed":bool(condition),**detail}); assert condition,name
    cert=build_certificate()
    ck("exact rational public certificate",True)
    bad=deepcopy(cert); bad["circular_planar_obstruction"]["enclosure"]=["0","1"]
    try:
        check_certificate(bad)
    except AssertionError:rejected=True
    else:rejected=False
    ck("positive instead of negative minor certificate rejected",rejected)
    # Necessary signs on legitimate analysis controls, not fitted FIN repairs.
    x=[I(3),2+I(*cert["sqrt3"]),I(4)]
    linear=[v for v in x]
    ck("linear Bernstein control has zero concavity defect",
       ((linear[2]-linear[1])/(x[2]-x[1])-(linear[1]-linear[0])/(x[1]-x[0])).lo<=0<=
       ((linear[2]-linear[1])/(x[2]-x[1])-(linear[1]-linear[0])/(x[1]-x[0])).hi)
    ck("uncrossed two-edge circular planar control has positive minor", F(1)*1-F(0)*0==1)
    # Exact critical ground-state transform, kappa=chi*lambda.
    z,kappa=sy.symbols("z kappa",real=True)
    V=z**4/12-kappa*z**2/2; psi=sy.exp(-V/2)
    potential=sy.expand(sy.diff(V,z)**2/4-sy.diff(V,z,2)/2)
    ck("critical transformed ground state exact",sy.simplify(-sy.diff(psi,z,2)+potential*psi)==0)
    ck("critical sextic potential exact at lambda0",sy.simplify(potential.subs(kappa,0)-(z**6/36-z**2/2))==0)
    # Independent labelled N2 parent; small algebra fixture, not MT asymptotics.
    spec=importlib.util.spec_from_file_location("seed",SEED); seed=importlib.util.module_from_spec(spec); spec.loader.exec_module(seed)
    _,_,lam,_,A=seed.build_rank7(); g=12/lam[6]
    states=[(a,b) for a in range(12) for b in range(12)]; L=np.zeros((144,144))
    for row,(a,b) in enumerate(states):
        for slot,other,old in ((0,b,a),(1,a,b)):
            probs=softmax(g/2*A[:,other])
            for new in range(12):
                if new==old:continue
                target=12*new+b if slot==0 else 12*a+new
                L[row,target]+=probs[new]; L[row,row]-=probs[new]
    logpi=np.array([g/2*A[a,b] for a,b in states]); pi=softmax(logpi)
    ck("native labelled generator and DB fixture",np.max(abs(L.sum(1)))<2e-14 and np.max(abs(pi[:,None]*L-(pi[:,None]*L).T))<3e-16)
    f1=np.array([math.cos(2*math.pi*a/12) for a,b in states]); f3=np.array([math.cos(6*math.pi*b/12) for a,b in states])
    P1=expm(.2*L); P2=expm(.5*L)
    histories=np.c_[P1@f1,P2@f3,P1@(f1*(expm(.3*L)@f3))]
    gram=histories.T@(pi[:,None]*histories)
    ck("future-cylinder reflection Gram fixture",np.linalg.eigvalsh(gram).min()>-2e-14)
    S=np.sqrt(pi)[:,None]*L/np.sqrt(pi)[None,:]
    ck("native transfer positivity and ground vector fixture",np.max(abs(S-S.T))<2e-14 and np.linalg.eigvalsh(-S).min()>-2e-14 and np.max(abs(S@np.sqrt(pi)))<2e-14)
    # A change of exponential changes predictions, even at the g0 reset control.
    H=np.eye(12)-np.ones((12,12))/12
    heat=expm(-2*math.pi*H)[0,0]; U=expm(-2j*math.pi*H)
    ck("heat is not physical unitary evolution",abs(abs(U[0,0])**2-1)<2e-13 and heat<.09)
    result={"date":"2026-10-15","certificate":cert,"checks":checks,
            "diagnostic_scope":"finite floating fixtures and symbolic identities; no mechanized asymptotic proof",
            "critical_potential":str(potential),"ground_state":str(psi),
            "input_hashes":{str(p.relative_to(ROOT)):sha(p) for p in (
                SPECTRA,SEED,Path(__file__),
                ROOT/"fin_critical_metastability_review_20261015/PROOF_COMPLETIONS.md",
                ROOT/"fin_intrinsic_observation_review_20261015/PROOF_COMPLETIONS_AND_CORRECTIONS.md",
                ROOT/"fin_uniform_critical_review_20261014/INTAKE_20261014.md",
                ROOT/"fin_kernel_observability_review_20261015/PROOF_COMPLETIONS_AND_ARITHMETIC.md",
                ROOT/"FIN_POST_PHY_TASK_CATALOG.md",
                ROOT/"fin_post_phy_review/RESEARCH_INTAKE_20260924.md")}}
    if args.check:
        saved=json.loads(args.check.read_text()); check_certificate(saved["certificate"])
        assert saved["certificate"]==cert
        assert all(saved["input_hashes"][p]==h for p,h in result["input_hashes"].items())
    if args.output:args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
    print(json.dumps({"checks":len(checks),"concavity_defect":cert["bernstein_obstruction"]["concavity_defect"],
                      "circular_minor":cert["circular_planar_obstruction"]["enclosure"],"status":"PASS_SCOPED_BRIDGE_CHECKS"},indent=2))

if __name__=="__main__":main()
