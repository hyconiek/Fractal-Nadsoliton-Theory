import numpy as np
from scipy.optimize import minimize
Z=np.load('/mnt/data/FIN_NEXT_RESEARCH_20260925/quartic_tensor_certificate.npz')
alphas=np.array(Z['alphas'],int);g=3.7183448981203875;c=Z['pred'][0]+g*Z['pred'][1]
# restrict to aligned real coordinates 3c,4c,5c,6 = original indices 0,2,4,6
sel=[0,2,4,6]
terms=[]
for a,cc in zip(alphas,c):
 if abs(cc)>1e-14 and all(a[i]==0 for i in range(7) if i not in sel):
  terms.append((np.array([a[i] for i in sel],int),float(cc)))
def vgh(x):
 x=np.asarray(x,float);v=0.;gr=np.zeros(4);H=np.zeros((4,4))
 for a,cc in terms:
  # value
  mon=cc
  for i,p in enumerate(a):
   if p: mon*=x[i]**p
  v+=mon
  for i,p in enumerate(a):
   if not p: continue
   t=cc*p
   for k,q in enumerate(a):
    e=q-(1 if k==i else 0)
    if e:t*=x[k]**e
   gr[i]+=t
   for j,q in enumerate(a):
    if not q-(1 if j==i else 0):
     # diagonal requires p>=2; off diag q can be >=1
     pass
    coef=p*(q-(1 if j==i else 0))
    if coef==0: continue
    t2=cc*coef
    for k,r in enumerate(a):
     e=r-(1 if k==i else 0)-(1 if k==j else 0)
     if e:t2*=x[k]**e
    H[i,j]+=t2
 return v,gr,H
cons={'type':'eq','fun':lambda x:np.dot(x,x)-1,'jac':lambda x:2*x}
x0=np.array([.29322386,.48635378,.65210399,.50223516])
r=minimize(lambda x:-vgh(x)[0],x0,jac=lambda x:-vgh(x)[1],constraints=[cons],method='SLSQP',options={'ftol':1e-15,'maxiter':2000})
x=r.x/np.linalg.norm(r.x);v,gr,H=vgh(x)
lagres=gr-4*v*x
# tangent basis via QR of null(x^T)
U,S,Vh=np.linalg.svd(x.reshape(1,-1));T=Vh[1:].T
HL=T.T@(H-4*v*np.eye(4))@T
print('success',r.success,r.message)
print('x',*[f'{q:.16g}' for q in x])
print('C',f'{v:.17g}')
print('lagrange_residual',np.max(np.abs(lagres)))
print('tangent_lagrangian_eigs',np.linalg.eigvalsh(HL).tolist())
print('pure_k5_unit',float(sum(cc for a,cc in terms if tuple(a)==(0,0,4,0))))
print('gain_ratio',v/float(sum(cc for a,cc in terms if tuple(a)==(0,0,4,0))))
print('polynomial_terms')
for a,cc in terms: print(tuple(a),repr(cc))
