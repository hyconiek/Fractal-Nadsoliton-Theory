#!/usr/bin/env python3
"""Finite 12-bit ternary exponential-family helpers for ST2442--ST2531."""
import itertools,math,numpy as np
from fin_st2172_st2261_common import ROOT,TRIANGLES,cycle_field,strict_A_W,write_packet,write_round

def data():
 _,W=strict_A_W();tau=cycle_field(W);X=np.array(list(itertools.product((-1,1),repeat=12)),dtype=np.int8)
 Y=np.array([[int(x[i]*x[j]*x[k]) for i,j,k in TRIANGLES] for x in X],dtype=float);S=Y@tau
 return tau,X,Y,S
def law(theta):
 tau,X,Y,S=data();z=theta*S;m=z.max();p=np.exp(z-m);p/=p.sum()
 tri=p@Y;mean=p@X;pairs=np.array([np.sum(p*X[:,i]*X[:,j]) for i in range(12) for j in range(i+1,12)])
 return {'theta':theta,'logZ':float(m+math.log(np.exp(z-m).sum())),'S_mean':float(p@S),'S_var':float(p@(S**2)-(p@S)**2),
  'tri_min':float(tri.min()),'tri_max':float(tri.max()),'mean_max_abs':float(abs(mean).max()),'pair_max_abs':float(abs(pairs).max()),
  'linear_response_residual':float(abs(tri/theta-tau).max()) if theta else 0.0}
