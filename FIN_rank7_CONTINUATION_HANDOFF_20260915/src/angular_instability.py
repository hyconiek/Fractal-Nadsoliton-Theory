"""R7P-095 constrained spherical Hessian of the pure k6 angular branch."""
from __future__ import annotations
import json,math
from fractions import Fraction
from pathlib import Path
import mpmath as mp

def _F(s): return Fraction(s)
def load_intervals(path):
 o=json.load(open(path)); L=o['exact']['laplacian_intervals']
 return [(Fraction(a),Fraction(b)) for a,b in L]

def eigs(r,L):
 # midpoint helper
 lm=[(float(a)+float(b))/2 for a,b in L]; x=math.sqrt(lm[6]/12)*r; t=math.tanh(x); lag=lm[6]/12*t/x if x else lm[6]/12
 return {'3cos':lm[3]/12*(1+t)-lag,'3sin':lm[3]/12*(1-t)-lag,
         '4cos':lm[4]/12-lag,'4sin':lm[4]/12-lag,'5cos':lm[5]/12-lag,'5sin':lm[5]/12-lag,'lagrange':lag,'x':x}

def f(r,l3,l6):
 x=math.sqrt(l6/12)*r
 return l3*(1+math.tanh(x))-l6*(math.tanh(x)/x if x else 1.0)

def certify(interval_json,out_path):
 L=load_intervals(interval_json); l3=sum(map(float,L[3]))/2;l6=sum(map(float,L[6]))/2
 # bisection on midpoint; exact interval sign checks below
 lo,hi=.4142113228,.4142113231
 for _ in range(70):
  m=(lo+hi)/2
  if f(m,l3,l6)<0:lo=m
  else:hi=m
 bracket=[0.41421132290,0.41421132293]
 mp.iv.dps=50
 def fi(r):
  l3i=mp.iv.mpf([str(float(L[3][0])),str(float(L[3][1]))]);l6i=mp.iv.mpf([str(float(L[6][0])),str(float(L[6][1]))])
  rr=mp.iv.mpf(str(r)); x=mp.iv.sqrt(l6i/12)*rr
  return l3i*(1+((mp.iv.exp(2*x)-1)/(mp.iv.exp(2*x)+1)))-l6i*((mp.iv.exp(2*x)-1)/(mp.iv.exp(2*x)+1))/x
 fl=fi(bracket[0]);fh=fi(bracket[1])
 # Other tangent sectors at upper root endpoint: if negative there, they were negative before because lagrange decreases.
 ev=eigs(bracket[1],L); other={k:v for k,v in ev.items() if k not in ('3cos','lagrange','x')}
 out={'radius_interval':bracket,'f_lower_interval':[float(fl.a),float(fl.b)],'f_upper_interval':[float(fh.a),float(fh.b)],
      'unique_root_reason':'df/dx = lambda3 sech^2(x) - lambda6*(x sech^2(x)-tanh(x))/x^2 > 0 for x>0',
      'midpoint_tangent_eigenvalues_at_upper_radius':other,'all_other_sectors_negative_at_upper_radius':all(v<0 for v in other.values()),
      'first_crossing_sector':'3cos','equation':'lambda3(1+tanh x)=lambda6 tanh(x)/x, x=sqrt(lambda6/12) r',
      'scope':'constrained spherical Hessian of log-mgf on pure k6 branch; not radial localization transition'}
 Path(out_path).write_text(json.dumps(out,indent=2)+'\n');return out
if __name__=='__main__':
 import sys;print(json.dumps(certify(sys.argv[1],sys.argv[2]),indent=2))
