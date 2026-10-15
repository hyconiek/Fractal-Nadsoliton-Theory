#!/usr/bin/env python3
"""Independent stdlib-only rational checker of the two class obstructions.

The transcendental input enclosure still relies on the declared interval backend.
"""
from copy import deepcopy
from fractions import Fraction
import hashlib
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent

def plus(x,y):return (x[0]+y[0],x[1]+y[1])
def neg(x):return (-x[1],-x[0])
def minus(x,y):return plus(x,neg(y))
def times(x,y):
    products=[a*b for a in x for b in y]
    return (min(products),max(products))
def divide(x,y):
    assert y[0]>0
    return times(x,(1/y[1],1/y[0]))
def interval(values):return tuple(Fraction(v) for v in values)

def verify(c):
    raw=json.loads((ROOT/c["spectral_path"]).read_text())["exact"]["laplacian_intervals"]
    assert hashlib.sha256((ROOT/c["spectral_path"]).read_bytes()).hexdigest()==c["spectral_sha256"]
    l={int(k):interval(v) for k,v in c["spectral_inputs"].items()}
    assert all(l[k]==interval(raw[k]) for k in l)
    s=interval(c["sqrt3"])
    assert s[0]**2<3<s[1]**2
    one=(Fraction(1),Fraction(1)); two=(Fraction(2),Fraction(2))
    left=divide(minus(l[5],l[4]),minus(s,one))
    right=divide(minus(l[6],l[5]),minus(two,s))
    defect=minus(right,left)
    public=interval(c["bernstein_obstruction"]["concavity_defect"])
    assert 0<public[0]<=defect[0]<=defect[1]<=public[1]
    w={int(k):interval(v) for k,v in c["kernel_weight_inputs"].items()}
    determinant=minus(times(w[6],w[3]),times(w[4],w[5]))
    public=interval(c["circular_planar_obstruction"]["enclosure"])
    assert public[0]<=determinant[0]<=determinant[1]<=public[1]<0
    info=c["circular_planar_obstruction"]
    assert info["P"]+info["Q"][::-1]==info["circular_order"]==[0,1,4,6]

if __name__=="__main__":
    saved=json.loads((HERE/"verification.json").read_text())
    verify(saved["certificate"])
    bad=deepcopy(saved["certificate"])
    bad["bernstein_obstruction"]["concavity_defect"]=["-1","0"]
    try:verify(bad)
    except AssertionError:pass
    else:raise AssertionError("invalid public sign accepted")
    print("PASS: two exact-rational sign certificates, inputs/order and altered-sign rejection")
