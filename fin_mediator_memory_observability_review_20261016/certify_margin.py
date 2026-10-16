#!/usr/bin/env python3
"""Independent stable-factor interval certificate, not a finite-N theorem."""
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import mpmath as mp

ROOT=Path(__file__).resolve().parents[1]
SPECTRA=ROOT/'fin_rank7_followup/inputs/fin_handoff_audit/results.json'

def dyadic(t):
    sign,m,e,_=t
    v=F(int(m)*(1<<int(e))) if e>=0 else F(int(m),1<<int(-e))
    return -v if sign else v

def ends(x):return tuple(dyadic(t) for t in x._mpi_)

def public(x,digits=24):
    lo,hi=ends(x); scale=10**digits
    a=lo.numerator*scale//lo.denominator
    b=-((-hi.numerator*scale)//hi.denominator)
    def dec(n):
        s='-' if n<0 else ''; n=abs(n)
        return s+str(n//scale)+'.'+str(n%scale).zfill(digits)
    result=[dec(a),dec(b)]
    assert F(result[0])<=lo<=hi<=F(result[1])
    return result

def build():
    mp.iv.dps=80; iv=mp.iv
    def rat(s):
        x=F(s); return iv.mpf(x.numerator)/x.denominator
    boxes=json.loads(SPECTRA.read_text())['exact']['laplacian_intervals']
    lam={r:iv.mpf([rat(boxes[r][0]).a,rat(boxes[r][1]).b]) for r in (3,6)}
    # Exact-decimal target and mathematical pi: independent paired-distance DFT.
    for r in (3,6):
        terms=[]
        for d in range(1,7):
            w=iv.cos(rat('0.18575')*d+rat('0.1625'))/(1+iv.mpf(d)**rat('1.8'))
            terms.append((2 if d<6 else 1)*w*(1-iv.cos(2*iv.pi*r*d/12)))
        low,high=ends(sum(terms)); assert F(boxes[r][0])<=low<=high<=F(boxes[r][1])
    c=lam[3]/lam[6]; gamma=1-c; q=iv.exp(-1)
    assert F('0.16257283698729')<ends(gamma)[0]<ends(gamma)[1]<F('0.16257283698731')
    result={'date':'2026-10-16','status':'VALIDATED_IDEAL_FILTERED_MARGIN_UNDER_INPUT_AND_IV_ASSUMPTIONS',
        'spectral_path':str(SPECTRA.relative_to(ROOT)),
        'spectral_sha256':hashlib.sha256(SPECTRA.read_bytes()).hexdigest(),
        'method':'fresh interval paired-distance spectrum binding; stable factored determinant; exact raw-dyadic outward decimals',
        'gamma':public(gamma),'q':public(q),'alpha':'1/2','spacing':'1',
        'scope':'ideal normalized sign-law functional, not finite-N error, detector calibration or finite-data decision'}
    for k in (1,4):
        root=iv.sqrt((1-k)**2+4*k*c); slow=(1+k-root)/2; fast=(1+k+root)/2
        ws=(fast-gamma)/(fast-slow); wf=(gamma-slow)/(fast-slow)
        assert ends(slow)[0]>0 and ends(slow)[1]<ends(gamma)[0]
        assert ends(fast)[0]>1 and ends(ws)[0]>0 and ends(wf)[0]>0
        x,y=iv.exp(-slow),iv.exp(-fast)
        norm=1/gamma+1
        A,B=ws/gamma/norm,wf/gamma/norm
        h=A*B*x*y*(x-q)*(y-q)*(x-y)**2
        assert ends(h)[1]<0
        result[f'kappa_{k}']={'beta_s':public(slow),'beta_f':public(fast),'H':public(h)}
    return result

if __name__=='__main__':
    p=argparse.ArgumentParser(); p.add_argument('--output',type=Path); a=p.parse_args(); data=build()
    if a.output:a.output.write_text(json.dumps(data,indent=2)+'\n')
    print(json.dumps(data,indent=2))
