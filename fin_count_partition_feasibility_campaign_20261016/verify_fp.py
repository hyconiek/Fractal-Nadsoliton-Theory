#!/usr/bin/env python3
"""Bounded independent diagnostics: not an all-N q12/MT feasibility proof."""
import json
from fractions import Fraction as F
from pathlib import Path
import mpmath as mp

def execute():
 checks=[]
 def check(name, truth, detail=None):
  checks.append({'name':name,'pass':bool(truth),'detail':detail})
  if not truth: raise AssertionError(name)
 # Exact rational finite directed two-edge chain: parity 0 -> 1 -> 2.
 edges=[(0,1),(1,2)]
 A=[[F(0) for e in edges] for x in range(3)]
 B=[[F(0) for e in edges] for x in range(3)]
 for j,(u,v) in enumerate(edges):
  A[u][j]=F(1); A[v][j]=F(-1)
  B[u][j]=B[v][j]=F(1)
 c=[F(2),F(3)];s=[sum(A[i][j]*c[j] for j in range(2)) for i in range(3)]
 check('primal diverges correctly',s==[2,1,-3])
 check('total source zero',sum(s)==0)
 check('phase cut1 current',s[0]==c[0])
 check('phase cut2 current',sum(s[:2])==c[1])
 load=[sum(B[i][j]*c[j] for j in range(2)) for i in range(3)]
 check('rate load',load==[2,5,3])
 check('primal feasible at cap5',all(v<=5 for v in load))
 check('primal not capacity-compatible at cap4',any(v>4 for v in load))
 phi=[F(-1),F(0),F(1)];beta=[F(0),F(1),F(0)]
 check('dual edge inequality',all(phi[u]-phi[v]+beta[u]+beta[v]>=0 for u,v in edges))
 dual=sum(phi[i]*s[i] for i in range(3))+4*sum(beta)
 check('negative cap4 dual',dual==F(-1))
 check('cap5 dual zero',sum(phi[i]*s[i] for i in range(3))+5*sum(beta)==0)
 # Integer harmonic 3 and parity rotation/reflection fixture.
 def four(n):
  C=sum(n[j]*(1,0,-1,0)[j%4] for j in range(12))
  S=sum(n[j]*(0,1,0,-1)[j%4] for j in range(12))
  v=sum((1 if j%2==0 else -1)*n[j] for j in range(12))
  return C,S,v
 n=[19,14,13,20,8,17,24,21,16,9,15,18]
 nrot=n[-1:]+n[:-1];nref=[n[(-j)%12] for j in range(12)]
 C,S,v=four(n);Cr,Sr,vr=four(nrot);Cf,Sf,vf=four(nref)
 check('C,S exact integers',all(isinstance(i,int) for i in [C,S,Cr,Sr]))
 check('rotation rotates radial 3-mode',C*C+S*S==Cr*Cr+Sr*Sr)
 check('rotation reverses parity',vr==-v)
 check('reflection preserves radius',Cf*Cf+Sf*Sf==C*C+S*S)
 check('reflection preserves parity',vf==v)
 # Exact source cut and conditional bracket identity on toy chain.
 masses=[F(2),F(3),F(5)];Wtot=sum(masses)
 d=[s[i]/masses[i] for i in range(3)]
 cut0=sum(masses[i]/Wtot*d[i] for i in [0]);cut1=sum(masses[i]/Wtot*d[i] for i in [0,1]);
 check('Gibbs weighted cut source 0',cut0==c[0]/Wtot)
 check('Gibbs weighted cut source 1',cut1==c[1]/Wtot)
 # Ideal limiting scalar ODE, gamma does not enter this diagnostic.
 mp.mp.dps=55
 normalizer=2*mp.quad(lambda z:mp.exp(-z**4/12),[0,mp.inf])
 rho=lambda z:mp.exp(-z**4/12)/normalizer
 p=lambda z:mp.mpf('0.5')+mp.exp(-z*z)/8
 q=2*mp.quad(lambda z:mp.exp(-z*z)*rho(z),[0,mp.inf])
 M=2*mp.quad(lambda z:z*z*mp.exp(-z*z)*rho(z),[0,mp.inf])
 PA=mp.mpf('0.5')+q/8;PB=1-PA
 JAt0=mp.quad(lambda z:(z**3/3)*p(z)*rho(z),[0,mp.inf])
 JBt0=mp.quad(lambda z:(z**3/3)*(1-p(z))*rho(z),[0,mp.inf])
 gA0=2*JAt0/(p(0)*rho(0));gB0=2*JBt0/((1-p(0))*rho(0))
 cap=gB0/4
 tol=mp.mpf('1e-35')
 check('rho normalized',abs(2*mp.quad(rho,[0,mp.inf])-1)<tol)
 check('quartic fourth moment',abs(2*mp.quad(lambda z:z**4*rho(z),[0,mp.inf])-3)<tol)
 check('positive phase masses',0<PB<PA<1)
 check('positive ideal M',M>0)
 check('phase mixture at z0',abs(p(0)*gA0+(1-p(0))*gB0-2)<tol)
 check('nontrivial phase conditional noise',gA0<2<gB0)
 check('necessary cap strict >1/2',cap>mp.mpf('0.5'))
 check('finite cap number reproducible',abs(cap-mp.mpf('0.647021987471151126747518283966135604712868869'))<mp.mpf('1e-32'))
 meanA=2-M/(2*PA);meanB=2+M/(2*PB)
 check('weighted stationary mean conserved',abs(PA*meanA+PB*meanB-2)<tol)
 check('phase mean residual weighted positive = M',abs(PA*abs(meanA-2)+PB*abs(meanB-2)-M)<tol)
 vals={'mass_A':mp.nstr(PA,42),'mass_B':mp.nstr(PB,42),'M':mp.nstr(M,42),
       'gA_at_zero':mp.nstr(gA0,42),'gB_at_zero':mp.nstr(gB0,42),'necessary_C_lower':mp.nstr(cap,42),
       'phase_mean_G_A':mp.nstr(meanA,42),'phase_mean_G_B':mp.nstr(meanB,42)}
 out={'scope':'bounded exact toy algebra and high-precision noncertified ideal quadratures only, no actual q12 flow solution',
      'checks':checks,'passed':sum(x['pass'] for x in checks),'total':len(checks),'numerical_diagnostics':vals,
      'critical_unproved': ['all-N actual q12/MT nonnegative vertexwise conductance solution','all-N exact dual impossibility']}
 (Path(__file__).parent/'FINITE_DIAGNOSTICS.json').write_text(json.dumps(out,indent=2,ensure_ascii=False)+'\n')
 print(f"PASS {out['passed']}/{out['total']} bounded checks; required ideal C >= {vals['necessary_C_lower']}")
if __name__=='__main__':execute()
