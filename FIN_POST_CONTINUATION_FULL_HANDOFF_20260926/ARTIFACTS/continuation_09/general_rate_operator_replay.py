import itertools, functools, numpy as np
Q=12; u=1/Q; j=np.arange(Q)
def ec(k): return (np.cos(2*np.pi*k*j/Q)/np.sqrt(Q) if k==6 else np.sqrt(2/Q)*np.cos(2*np.pi*k*j/Q))
def es(k): return np.sqrt(2/Q)*np.sin(2*np.pi*k*j/Q)
E=np.column_stack([ec(3),es(3),ec(4),es(4),ec(5),es(5),ec(6)])
PH=np.eye(Q)-np.ones((Q,Q))/Q-E@E.T
rng=np.random.default_rng(123456)
c=rng.normal(size=7); c/=np.linalg.norm(c); phi=E@c
h=PH@(phi*phi); m=np.mean(phi*phi)
K=rng.normal(size=(Q,Q)); np.fill_diagonal(K,0)
a=K.sum(axis=0)-K.sum(axis=1)
D={(i,k):tuple((np.eye(Q,dtype=int)[k]-np.eye(Q,dtype=int)[i]).tolist()) for i in range(Q) for k in range(Q) if i!=k}
def add(s,d): return tuple(x+y for x,y in zip(s,d))
@functools.lru_cache(None)
def f(s): return float(phi@np.asarray(s,float))**4
def seqval(seq):
 @functools.lru_cache(None)
 def rec(pos,s):
  if pos==3:return f(s)
  base=rec(pos+1,s); op=seq[pos]; d=np.asarray(s,float); z=PH@d; out=0.0
  for i in range(Q):
   for k in range(Q):
    if i==k: continue
    rate=u*u if op=='G' else (u*z[i] if op=='D' else K[i,k])
    out += rate*(rec(pos+1,add(s,D[(i,k)]))-base)
  return out
 return rec(0,(0,)*Q)
vals={''.join(p):seqval(p) for p in set(itertools.permutations('GDC'))}
pred=36*m*np.dot(a,h)
print('permutations',vals)
print('sum',sum(vals.values()))
print('prediction',pred)
print('residual',sum(vals.values())-pred)
assert abs(sum(vals.values())-pred)<2e-13
