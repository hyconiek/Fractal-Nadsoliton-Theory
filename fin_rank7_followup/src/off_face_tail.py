from __future__ import annotations
import json
from pathlib import Path
import mpmath as mp
ROOT=Path(__file__).resolve().parents[1]
mp.iv.dps=60

def iv(a,b,c,d): return mp.iv.mpf([mp.mpf(a)/b,mp.mpf(c)/d])
L3=iv(196140686197643,10**14,39228137239529,2*10**13)
L4=iv(5498922123333,2500000000000,109978442466661,50000000000000)
L5=iv(57465156801977,25000000000000,22986062720791,10000000000000)
L6=iv(234218204114629,100000000000000,234218204114631,100000000000000)

def run():
    t=mp.iv.mpf(1)/2**11; pi=mp.iv.pi; tr=mp.iv.mpf(0)
    rows=[]
    for j in range(1,12):
        c3=mp.iv.cos(2*pi*3*j/12); c4=mp.iv.cos(2*pi*4*j/12); c5=mp.iv.cos(2*pi*5*j/12); c6=(-1)**j
        ratio=t**(2*(1-c5))
        d2=L3/6*(c3-1)**2+L4/6*(c4-1)**2+L5/6*(c5-1)**2+L6/12*(c6-1)**2
        tr += ratio*d2
        rows.append({'j':j,'ratio_upper':str(ratio.b),'distance2_upper':str(d2.b)})
    sigma=(2*L3*(L4+L5)-L4*L5)/(24*L3)
    ok=tr.b < 2*sigma.a
    return {'task':'R7P-069','tail_domain':'t=exp(-J5/2)<=2^-11, arbitrary J3,J4,J6>=0',
      'argument':'anchor j=0; every non-anchor weight ratio is <=t^(2(1-c5_j)); trace Cov <= sum ratio_j ||F_j-F_0||^2; PSD gives lambda2<=trace/2',
      'trace_upper':str(tr.b),'two_sigma_lower':str((2*sigma).a),'strict':bool(ok),'rows':rows,
      'other_certified_subdomains':['boundary-Ising closure R7P-055','intraparity W R7P-063','extreme face R7P-017--024','local off-face cone rho=1/8192 R7P-068'],
      'residual':'compactified core with t>2^-11 excluding the certified boundary/face/local sets remains unresolved',
      'global_4D_ceiling':'NOT_PROVED'}

def main():
    r=run(); (ROOT/'results/R7P-069_072_off_face_partial.json').write_text(json.dumps(r,indent=2)+'\n'); print(json.dumps({'strict':r['strict'],'trace_upper':r['trace_upper'],'two_sigma_lower':r['two_sigma_lower']},indent=2))
if __name__=='__main__': main()
