from pathlib import Path
from fractions import Fraction as F
import sys,itertools,numpy as np
ROOT=Path('/mnt/data/fin_rank7_next_campaign');sys.path[:0]=[str(ROOT/'src')]
import target_p_trace_tight_rounded as tr
import target_p_trace_cover as b
TAU=F(67,250)
# exact conservative quadratic Q: P=w^T Q w >= D^2(trace-2tau) using D2 upper.
Q=[[F(0) for _ in range(7)] for __ in range(7)]
for i in range(7): Q[i][i]=-2*TAU
for i in range(7):
 for j in range(i+1,7): Q[i][j]=Q[j][i]=(tr.D2[i][j]-4*TAU)/2
# Bernstein grid and integer feature arrays for degree 2 in each variable
K=np.array(list(itertools.product((0,1,2),repeat=7)),dtype=np.int8)
KF=K/2.0
SQ=(K==2).astype(float) # choose(k,2): 1 iff k=2
CROSS=np.stack([(K[:,i]*K[:,j])/4.0 for i in range(7) for j in range(i+1,7)],axis=1)

def coeffs_float(cell):
 wb=tr.weights(cell); l=np.array([float(x[0]) for x in wb]); d=np.array([float(x[1]-x[0]) for x in wb]); q=np.array([[float(Q[i][j]) for j in range(7)] for i in range(7)])
 c0=float(l@q@l); lin=2*d*(q@l); sq=np.diag(q)*d*d
 cross=np.array([2*q[i,j]*d[i]*d[j] for i in range(7) for j in range(i+1,7)])
 vals=c0+KF@lin+SQ@sq+CROSS@cross
 return float(vals.max()),int(vals.argmax())

def coeff_at_exact(cell,krow):
 wb=tr.weights(cell); l=[x[0] for x in wb];d=[x[1]-x[0] for x in wb]
 c0=sum(l[i]*Q[i][j]*l[j] for i in range(7) for j in range(7))
 lin=[2*d[i]*sum(Q[i][j]*l[j] for j in range(7)) for i in range(7)]
 sq=[Q[i][i]*d[i]*d[i] for i in range(7)]
 v=c0
 for i,k in enumerate(krow): v += lin[i]*F(int(k),2) + (sq[i] if k==2 else 0)
 for i in range(7):
  for j in range(i+1,7): v += 2*Q[i][j]*d[i]*d[j]*F(int(krow[i]*krow[j]),4)
 return v

def certify_float_candidate(cell,margin=1e-12):
 mx,idx=coeffs_float(cell)
 if mx>=-margin:return {'ok':False,'float_max':mx,'idx':idx}
 # exact check all coefficients only if float suggests pass; stop on first positive.
 kbest=None;best=None
 for row in K:
  v=coeff_at_exact(cell,row)
  if best is None or v>best:best=v;kbest=row.tolist()
  if v>0:return {'ok':False,'float_max':mx,'exact_positive':str(v),'idx':row.tolist()}
 return {'ok':True,'float_max':mx,'exact_max':str(best),'idx':kbest}
