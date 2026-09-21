"""R7P-013 small exact-rational transcendental enclosures.

Proof decisions use Fraction endpoints and explicit Taylor/geometric remainders.
The routines are intentionally conservative and scoped to moderate finite inputs.
"""
from __future__ import annotations
from fractions import Fraction as F
import math

class QI:
    def __init__(self,lo,hi=None):
        self.lo=F(lo); self.hi=F(lo if hi is None else hi)
        if self.lo>self.hi: raise ValueError('reversed interval')
    def __repr__(self): return f'QI({self.lo},{self.hi})'
    @staticmethod
    def cast(x): return x if isinstance(x,QI) else QI(x)
    def __add__(self,o): return add(self,QI.cast(o))
    __radd__=__add__
    def __neg__(self): return QI(-self.hi,-self.lo)
    def __sub__(self,o): return self+(-QI.cast(o))
    def __rsub__(self,o): return QI.cast(o)-self
    def __mul__(self,o): return mul(self,QI.cast(o))
    __rmul__=__mul__
    def __truediv__(self,o): return div(self,QI.cast(o))
    def __rtruediv__(self,o): return div(QI.cast(o),self)
    def __pow__(self,n):
        if n<0:return QI(1)/(self**(-n))
        if n==0:return QI(1)
        out=QI(1)
        for _ in range(n):out=out*self
        return out

def reciprocal(I:QI):
    if I.lo<=0<=I.hi: raise ValueError('zero denominator')
    return QI(1/I.hi,1/I.lo)

def add(A:QI,B:QI): return QI(A.lo+B.lo,A.hi+B.hi)
def mul(A:QI,B:QI):
    xs=[A.lo*B.lo,A.lo*B.hi,A.hi*B.lo,A.hi*B.hi]
    return QI(min(xs),max(xs))
def div(A:QI,B:QI): return mul(A,reciprocal(B))

def sqrt_point(x:F,digits=30):
    x=F(x)
    if x<0: raise ValueError('sqrt negative')
    if x==0:return QI(0)
    Q=10**digits
    A=x.numerator*Q*Q; B=x.denominator
    n=math.isqrt(A//B)
    while (n+1)*(n+1)*B<=A:n+=1
    while n*n*B>A:n-=1
    if n*n*B==A:return QI(F(n,Q))
    return QI(F(n,Q),F(n+1,Q))

def sqrt_interval(I:QI,digits=30):
    if I.lo<0: raise ValueError('sqrt interval negative')
    lo=sqrt_point(I.lo,digits).lo; hi=sqrt_point(I.hi,digits).hi
    return QI(lo,hi)

def _exp_nonneg_point(x:F):
    x=F(x)
    if x<0: raise ValueError
    N=max(80,int(math.ceil(float(x)))+60)
    term=F(1); s=term
    for k in range(1,N+1):
        term=term*x/k; s+=term
    next_term=term*x/(N+1)
    ratio=x/(N+2)
    if ratio>=1: raise ValueError('increase Taylor order')
    rem=next_term/(1-ratio)
    return QI(s,s+rem)

def exp_point(x:F):
    x=F(x)
    if x>=0:return _exp_nonneg_point(x)
    return reciprocal(_exp_nonneg_point(-x))

def exp_interval(I:QI):
    return QI(exp_point(I.lo).lo,exp_point(I.hi).hi)

def _log_ge1_point(x:F,N=180):
    if x<1: raise ValueError
    if x==1:return QI(0)
    z=(x-1)/(x+1)
    s=F(0)
    for k in range(N): s+=z**(2*k+1)/(2*k+1)
    lo=2*s
    tail=2*z**(2*N+1)/((2*N+1)*(1-z*z))
    return QI(lo,lo+tail)

def log_point(x:F):
    x=F(x)
    if x<=0: raise ValueError('log nonpositive')
    if x>=1:return _log_ge1_point(x)
    J=_log_ge1_point(1/x)
    return QI(-J.hi,-J.lo)

def log_interval(I:QI):
    if I.lo<=0: raise ValueError('log nonpositive interval')
    return QI(log_point(I.lo).lo,log_point(I.hi).hi)

def tanh_point(x:F):
    E=exp_point(2*F(x)); one=QI(1)
    return div(QI(E.lo-1,E.hi-1),QI(E.lo+1,E.hi+1))

def tanh_interval(I:QI):
    return QI(tanh_point(I.lo).lo,tanh_point(I.hi).hi)

def sech_point(x:F):
    x=abs(F(x)); ep=exp_point(x); em=reciprocal(ep)
    return div(QI(2),add(ep,em))

def sech_interval(I:QI):
    # sech is even and decreases with |x|
    mx=max(abs(I.lo),abs(I.hi)); mn=F(0) if I.lo<=0<=I.hi else min(abs(I.lo),abs(I.hi))
    return QI(sech_point(mx).lo,sech_point(mn).hi)
