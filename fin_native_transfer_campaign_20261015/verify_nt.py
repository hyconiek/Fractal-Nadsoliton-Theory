#!/usr/bin/env python3
"""Bounded independent numerical/algebraic diagnostics for FIN NT-001..003.
No stochastic convergence theorem is mechanically proved by this file.
"""
import math,json,itertools,hashlib,os
import numpy as np
import sympy as sp
from scipy.linalg import expm

Q=12
omega,phi,eta=0.18575,0.16250,1.8
W=np.zeros((Q,Q))
for i in range(Q):
    for j in range(Q):
        if i!=j:
            d=min((i-j)%Q,(j-i)%Q)
            W[i,j]=math.cos(omega*d+phi)/(1+d**eta)
L=np.diag(W.sum(axis=1))-W
lam=np.fft.fft(L[0]).real
x=np.arange(Q)
cols=[]
for k in (3,4,5):
    cols.extend([np.sqrt(lam[k]/6)*np.cos(2*np.pi*k*x/Q),np.sqrt(lam[k]/6)*np.sin(2*np.pi*k*x/Q)])
cols.append(np.sqrt(lam[6]/12)*(-1.)**x)
X=np.array(cols).T
A=X@X.T
g_c=12/lam[6]

def checks():
  N=2 # deliberately tiny GLOBAL Gibbs convention fixture; NOT MT local window.
  states=[]
  for i in range(Q):
      for j in range(i,Q):
          n=np.zeros(Q,dtype=int);n[i]+=1;n[j]+=1;states.append(n)
  idx={tuple(n):i for i,n in enumerate(states)}
  assert len(states)==78
  def sf(y):
      e=np.exp(y-np.max(y));return e/e.sum()
  generator=np.zeros((len(states),len(states)))
  for ix,n in enumerate(states):
      for i in range(Q):
          if n[i]==0:continue
          field=(g_c/N)*(A@(n-np.eye(Q,dtype=int)[i]))
          probs=sf(field)
          for j in range(Q):
              if i==j:continue
              nxt=n.copy();nxt[i]-=1;nxt[j]+=1
              generator[ix,idx[tuple(nxt)]]+=n[i]*probs[j]
      generator[ix,ix]=-np.sum(generator[ix])
  weights=[]
  for n in states:
      multinomial=math.factorial(N)/math.prod(math.factorial(int(v)) for v in n)
      weights.append(multinomial*math.exp((g_c/(2*N))*float(n@A@n)))
  pi=np.array(weights);pi/=pi.sum()
  reversible_err=float(np.max(abs(pi[:,None]*generator -pi[None,:]*generator.T)))
  stationary_err=float(np.max(abs(pi@generator)))
  row_err=float(np.max(abs(generator.sum(1))))
  D=np.sqrt(pi); invD=1/D
  H=-(D[:,None]*generator)*invD[None,:]
  sym_err=float(np.max(abs(H-H.T)))
  ground_err=float(np.linalg.norm(H@D))
  eigen=np.linalg.eigvalsh(H)
  f=np.cos(np.array([float(n@((-1.)**x))/N for n in states])*0.7)
  t=0.45
  P=expm(t*generator);heat=expm(-t*H)
  corr=float(np.sum(pi*f*(P@f)))
  heat_corr=float((D*f)@(heat@(D*f)))
  corr_err=abs(corr-heat_corr)
  def E_form(f):return 0.5*sum(pi[i]*generator[i,j]*(f[j]-f[i])**2 for i in range(len(states)) for j in range(len(states)) if j!=i)
  E_direct=float(-f@(pi*(generator@f)))
  E_edge=float(E_form(f))
  # Two time-cylinder at s=.17 and t=.45 with complex functions; reflect positivity.
  s=0.17
  Ps=expm(s*generator); Pdt=expm((t-s)*generator)
  ai=np.array([np.cos(np.arange(78)*.12)+1j*np.sin(np.arange(78)*.24),np.sin(np.arange(78)*.22)-1j*np.cos(np.arange(78)*.21)]).T
  bi=np.array([np.sin(np.arange(78)*.07)+.3j*np.cos(np.arange(78)*.13),np.cos(np.arange(78)*.09)+.2j*np.sin(np.arange(78)*.20)]).T
  U=np.column_stack([Ps@(ai[:,k]*(Pdt@bi[:,k])) for k in range(2)])
  Gram=U.conj().T@(pi[:,None]*U)
  gram_min=float(np.linalg.eigvalsh(Gram).min())
  # Directly expand E|u(X0)|² through summation over all initial states.
  gram_diag_direct=float(sum(pi[q]*abs(sum(Ps[q,i]*ai[i,0]*sum(Pdt[i,j]*bi[j,0] for j in range(78)) for i in range(78)))**2 for q in range(78)))
  gram_diag_err=abs(gram_diag_direct-Gram[0,0].real)
  # Symbolic supersymmetric operator check, including lambda variable a.
  z,a=sp.symbols('z a', real=True)
  U_p=z**3/sp.Integer(3)-a*z
  V=sp.expand(U_p**2/sp.Integer(4)-sp.diff(U_p,z)/2)
  expected=z**6/36-a*z**4/6+(a*a/4-sp.Rational(1,2))*z*z+a/2
  psi0=sp.exp(-z**4/24+a*z**2/4)
  opground=sp.simplify(-sp.diff(psi0,z,2)+V*psi0)
  sigma=(-1.)**x
  vals={'state_count':78, 'rank_A7':int(np.linalg.matrix_rank(A,1e-10)), 'critical_gain_float':g_c,'reversibility_max':reversible_err,'stationary_max':stationary_err,'generator_rows_max':row_err,'H_symmetry_max':sym_err,'ground_null_norm':ground_err,'H_lowest_eigenvalue':float(eigen.min()),'finite_correlation_classical':corr,'finite_correlation_transfer':heat_corr,'correlation_max_error':corr_err,'edge_form_error':abs(E_direct-E_edge),'reflection_Gram_eig_min':gram_min,'direct_cylinder_reflection_error':gram_diag_err,'symbolic_potential_match':bool(sp.simplify(V-expected)==0),'symbolic_ground_zero':bool(opground==0),'lambda0_potential_match':bool(sp.simplify(V.subs(a,0)-(z**6/36-z*z/2))==0)}
  assert vals['rank_A7']==7
  for k in ['reversibility_max','stationary_max','generator_rows_max','H_symmetry_max','ground_null_norm','correlation_max_error','edge_form_error','direct_cylinder_reflection_error']:
      assert vals[k]<1e-9,(k,vals[k])
  assert vals['H_lowest_eigenvalue']> -1e-9
  assert vals['reflection_Gram_eig_min']>-1e-11
  for k in ['symbolic_potential_match','symbolic_ground_zero','lambda0_potential_match']:
      assert vals[k]
  vals['bounded_fixtures_passed']=14
  vals['scientific_scope']='N=2 full Gibbs q12 fixture + symbolic 1D operator; not stationary local critical process or stochastic path theorem.'
  return vals
if __name__=='__main__':
  result=checks(); path=os.path.join(os.path.dirname(__file__),'verification.json')
  with open(path,'w') as f:json.dump(result,f,indent=2,ensure_ascii=False)
  print(json.dumps(result,indent=2))
